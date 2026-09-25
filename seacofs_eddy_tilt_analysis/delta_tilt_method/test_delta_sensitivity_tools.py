"""Synthetic invariants for the isolated temporal-weight experiment."""
import sys
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'seacofs_eddy_dataset_modular' / 'src'))
from seacofs_eddy_dataset.core.tilt import compute_weighted_tilt
from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.stages import tilt as tilt_stage
from delta_sensitivity_tools import increments, fit_snapshot, lifetime, maximum_span, scan_projected_extent_differences


def fixture():
    rows = []
    for day in range(20):
        z = np.arange(0, (800 if day % 3 else 600) + 1, 20.)
        rows.append(pd.DataFrame({'Eddy': 1, 'Day': day, 'Depth': z,
            'xc': 100+day+z*(.003+.001*np.sin(day)),
            'yc': 400-day+z*.01+np.sin(z/90+day)*.2}))
    return pd.concat(rows, ignore_index=True)


class SensitivityTests(unittest.TestCase):
    def test_stages_summarise_reference_limited_estimates(self):
        config = PipelineConfig(raw={'paths': {'output_root': '/unused'},
            'parallel': {'skip_existing': False}}, config_path=Path('/unused/config.yaml'))
        for all_shallow in [False, True]:
            p = fixture().loc[lambda d: d.Depth.le(200) if all_shallow
                else (~d.Day.eq(10) | d.Depth.le(200))]
            with patch.object(Path, 'exists', return_value=True), \
                 patch.object(tilt_stage.pd, 'read_parquet', return_value=p), \
                 patch.object(tilt_stage, 'write_partition') as write, \
                 patch('builtins.print'):
                tilt_stage.run(config)
                result = write.call_args.args[0]
            self.assertTrue(result.loc[result.Day.eq(10), ['TiltDis', 'TiltDir']].isna().all().all())
            with patch.object(Path, 'exists', return_value=True), \
                 patch.object(tilt_stage.pd, 'read_parquet', side_effect=[result, p]), \
                 patch.object(tilt_stage, 'write_partition') as write, \
                 patch('builtins.print'), warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                tilt_stage.run_analysis(config)
                summary = write.call_args_list[0].args[0].iloc[0]
            self.assertEqual(summary.n_valid_tilts, result.TiltDis.notna().sum())
            if all_shallow:
                self.assertTrue(np.isnan(summary.mean_tilt_distance_km))
                self.assertTrue(np.isnan(summary.median_tilt_distance_km))
            else:
                self.assertAlmostEqual(summary.mean_tilt_distance_km, result.TiltDis.mean())

    def test_production_matches_notebook_with_reference_depth_limits(self):
        p = fixture()
        p = p.loc[~p.Day.eq(10) | p.Depth.le(200)]
        p = p.loc[~p.Day.eq(11) | p.Depth.between(100, 400)]
        p = p.loc[~p.Day.eq(8)]
        dx, dy = increments(p)
        for num, sigma in [(5, 1.0), (6, None), (5, 1.5)]:
            result = compute_weighted_tilt(p, 1, num=num, temporal_sigma_days=sigma)
            for row in result.itertuples():
                expected = fit_snapshot(dx, dy, row.Day, sigma=sigma, half_window=num // 2)
                if expected is None:
                    self.assertTrue(np.isnan([row.TiltDis, row.TiltDir]).all())
                else:
                    np.testing.assert_allclose([row.TiltDis, row.TiltDir],
                        [expected['TiltDis'], expected['TiltDir']], atol=1e-9)

    def test_production_ignores_neighbours_below_reference(self):
        p = fixture().loc[lambda d: ~d.Day.eq(10) | d.Depth.le(400)].copy()
        before = compute_weighted_tilt(p, 1).set_index('Day').loc[10]
        p.loc[p.Depth.gt(400), 'xc'] += 10000
        p.loc[p.Depth.gt(400), 'yc'] -= 9000
        after = compute_weighted_tilt(p, 1).set_index('Day').loc[10]
        np.testing.assert_allclose(before, after, atol=1e-9)

    def test_production_profile_outside_depth_cap_has_no_estimate(self):
        p = fixture()
        p.loc[p.Day.eq(10), 'Depth'] += 2000
        result = compute_weighted_tilt(p, 1).set_index('Day').loc[10]
        self.assertTrue(result[['TiltDis', 'TiltDir']].isna().all())

    def test_five_day_gaussian_matches_production_when_support_is_equal(self):
        p = fixture().loc[lambda d: d.Depth.le(600)]
        dx, dy = increments(p)
        actual = compute_weighted_tilt(p, 1)
        self.assertEqual(actual.Day.tolist(), list(range(2, 18)))
        for row in actual.itertuples():
            expected = fit_snapshot(dx, dy, row.Day, sigma=1.0, half_window=2)
            self.assertIsNotNone(expected)
            np.testing.assert_allclose(
                [row.TiltDis, row.TiltDir],
                [expected['TiltDis'], expected['TiltDir']], rtol=1e-9, atol=1e-9)

    def test_gaussian_central_day_contribution_and_window(self):
        z = np.arange(0, 301, 10.)
        p = pd.concat([pd.DataFrame({'Eddy': 1, 'Day': day, 'Depth': z,
            'xc': z * (0.01 if day == 2 else 0.), 'yc': z * 0.})
            for day in range(5)], ignore_index=True)
        result = compute_weighted_tilt(p, 1)
        a = np.exp(-0.5 * np.arange(-2, 3)**2)
        self.assertEqual(result.Day.tolist(), [2])
        # Upper-interval labels span 0..290 m, preserving existing convention.
        self.assertAlmostEqual(result.TiltDis.iloc[0], 290 * .01 * a[2] / a.sum())
        self.assertAlmostEqual(result.TiltDir.iloc[0], 290.)
        self.assertTrue(compute_weighted_tilt(p.loc[p.Day.lt(4)], 1).empty)

    def test_invalid_sigma(self):
        for sigma in [0, -1, np.nan, np.inf]:
            with self.assertRaises(ValueError):
                compute_weighted_tilt(fixture(), 1, temporal_sigma_days=sigma)

    def test_equal_matches_production_when_support_is_equal(self):
        p = fixture().loc[lambda d: d.Depth.le(600)]
        table, _, _ = lifetime(p)
        baseline = compute_weighted_tilt(p, 1, num=6, temporal_sigma_days=None)
        joined = table.merge(baseline, on='Day')
        np.testing.assert_allclose(joined['Equal'], joined.TiltDis, rtol=1e-9, atol=1e-9, equal_nan=True)

    def test_shallow_reference_cannot_borrow_deeper_neighbours(self):
        p = fixture()
        p = p.loc[~p.Day.eq(10) | p.Depth.le(200)]
        dx, dy = increments(p)
        for sigma in [None, 1.0, 1.5]:
            # Neighbours reach 800 m, but cannot rescue the 190 m fit-label span.
            self.assertIsNone(fit_snapshot(dx, dy, 10, sigma=sigma))
            fit = fit_snapshot(dx, dy, 10, sigma=sigma, min_depth_range=100)
            self.assertEqual(fit['trace'].Depth.max(), 190)
            self.assertEqual(fit['ends'][1, 2], 190)
        # Changing only deeper neighbour centres has no effect on the estimate.
        old = fit_snapshot(dx, dy, 10, sigma=1, min_depth_range=100)
        p.loc[p.Depth.gt(200), 'xc'] += 10000
        p.loc[p.Depth.gt(200), 'yc'] -= 9000
        x2, y2 = increments(p)
        new = fit_snapshot(x2, y2, 10, sigma=1, min_depth_range=100)
        np.testing.assert_allclose(old['ends'], new['ends'])
        np.testing.assert_allclose(old['trace'], new['trace'])

    def test_missing_reference_day_has_no_lifetime_estimate(self):
        p = fixture().loc[lambda d: ~d.Day.eq(10)]
        table, dx, dy = lifetime(p)
        self.assertIsNone(fit_snapshot(dx, dy, 10))
        row = table.loc[table.Day.eq(10)].drop(columns='Day')
        self.assertTrue(row.isna().to_numpy().all())

    def test_nonzero_top_is_clipped_before_cumulative_sum(self):
        # Reference support 100..405 m: full intervals have labels 100..390.
        # Every day follows the same straight line. No increments from 0..100
        # may be included in the cumulative sum at the reference's first row.
        rows = []
        for day in range(7):
            z = np.r_[np.arange(100, 401, 10), 405.] if day == 3 else np.arange(0, 801, 10.)
            rows.append(pd.DataFrame({'Eddy':1, 'Day':day, 'Depth':z,
                                     'xc':day+z*.01, 'yc':day*2+z*.02}))
        dx, dy = increments(pd.concat(rows, ignore_index=True))
        fit = fit_snapshot(dx, dy, 3, sigma=1)
        self.assertEqual(fit['trace'].Depth.min(), 100)
        self.assertEqual(fit['trace'].Depth.max(), 390)
        np.testing.assert_allclose(fit['trace'][['x','y']].iloc[0], [.1,.2])
        self.assertAlmostEqual(fit['TiltDis'], 290*np.hypot(.01,.02))

    def test_lifetime_limits_each_day_independently(self):
        p = fixture()
        p = p.loc[~p.Day.eq(10) | p.Depth.le(400)]
        dx, dy = increments(p)
        for day in range(3,17):
            fit = fit_snapshot(dx, dy, day, sigma=1)
            self.assertIsNotNone(fit)
            deepest = p.loc[p.Day.eq(day), 'Depth'].max()
            self.assertLessEqual(fit['trace'].Depth.max()+10, deepest)
            self.assertLessEqual(fit['ends'][1,2]+10, deepest)

    def test_depth_weights_unchanged_and_translation_invariant(self):
        p = fixture()
        dx, dy = increments(p)
        a = fit_snapshot(dx, dy, 10)
        b = fit_snapshot(dx, dy, 10, 1.5)
        np.testing.assert_allclose(a['trace'].weight, b['trace'].weight)
        p.xc += p.Day**2
        p.yc -= p.Day*123
        dx, dy = increments(p)
        c = fit_snapshot(dx, dy, 10, 1.5)
        np.testing.assert_allclose(b['TiltDis'], c['TiltDis'], atol=1e-8)

    def test_focus_day_contribution(self):
        p = fixture()
        dx, dy = increments(p)
        # Isolate one depth-increment pulse on the central day.
        dx.loc[:, :] = 0.
        dy.loc[:, :] = 0.
        dx.loc[:, 10] = 1.
        fit = fit_snapshot(dx, dy, 10, 1.5)
        a = np.exp(-.5*(np.arange(-3,4)/1.5)**2)
        expected = a[3]/a.sum()
        self.assertGreater(expected, .25)
        self.assertAlmostEqual(fit['trace'].x.iloc[0], expected)

    def test_ranked_projected_difference_known_straight_profiles(self):
        z = np.arange(0, 401, 10.)
        profiles = pd.concat([pd.DataFrame({'Eddy': eddy, 'Day': day,
            'Depth': z, 'xc': 0*z+day, 'yc': z*slope+day})
            for eddy, slope in [(1, .01), (2, .04)] for day in range(9)], ignore_index=True)
        progress = []
        ranked = scan_projected_extent_differences(profiles, progress=lambda *x:progress.append(x))
        self.assertEqual(len(ranked), 10)
        self.assertEqual(ranked.iloc[:5].Eddy.tolist(), [2]*5)
        np.testing.assert_allclose(ranked.iloc[:5].MaxProjectedExtent, 16)
        np.testing.assert_allclose(ranked.iloc[:5].TiltDis, 15.6)
        np.testing.assert_allclose(ranked.iloc[:5].AbsoluteDifference, .4)
        self.assertEqual(ranked.Rank.tolist(), list(range(1,11)))
        self.assertTrue(ranked.AbsoluteDifference.is_monotonic_decreasing)
        self.assertEqual(progress[-1], (2,2,10))
        # Input order cannot change the ranking.
        pd.testing.assert_frame_equal(ranked,
            scan_projected_extent_differences(profiles.sample(frac=1, random_state=3)))

    def test_ranker_excludes_missing_shallow_and_zero_tilt_days(self):
        p = fixture()
        p = p.loc[(~p.Day.eq(10) | p.Depth.le(200)) & ~p.Day.eq(11)]
        ranked = scan_projected_extent_differences(p)
        self.assertNotIn(10, ranked.Day.values)
        self.assertNotIn(11, ranked.Day.values)
        p[['xc','yc']] = 0.
        self.assertTrue(scan_projected_extent_differences(p).empty)

    def test_span_is_pairwise_not_surface_relative(self):
        p = pd.DataFrame({'Depth':[0,100,200], 'xc':[0,-3,4], 'yc':[0,0,0]})
        self.assertEqual(maximum_span(p), 7.)


if __name__ == '__main__':
    unittest.main()
