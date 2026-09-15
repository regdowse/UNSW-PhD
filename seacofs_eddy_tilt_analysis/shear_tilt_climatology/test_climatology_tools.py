"""Deterministic scientific regression tests; run with unittest discovery."""
import unittest
import numpy as np
import pandas as pd
from climatology_tools import classify_regimes, relative_geometry, lagged_pairs, cluster_summary, wrap


def sample(days=range(10), bearings=None):
    n = len(days)
    d = pd.DataFrame(dict(Cyc=['AE']*n, Eddy=[1]*n, Day=list(days),
                         TiltDis=[10.]*n, TiltDir=[90.]*n if bearings is None else bearings,
                         h=[4000.]*n, water_depth_m=[4000.]*n,
                         PV_grad_plan_mag=[1.]*n, PV_grad_topo_mag=[.1]*n))
    for f in ['clim', 'full']:
        for level in ['200', '500', 'surface']:
            d[f'{f}_{level}_east_ms'] = .1 if level != '500' else .04
            d[f'{f}_{level}_north_ms'] = 0.
    return d


def geom(d):
    return relative_geometry(classify_regimes(d, window=1, min_periods=1))


class ScientificTests(unittest.TestCase):
    def test_cardinal_and_hemisphere_independent_conventions(self):
        d = geom(sample(range(4), [90, 270, 0, 180]))
        np.testing.assert_allclose(d.parallel_km, [10, -10, 0, 0], atol=1e-12)
        np.testing.assert_allclose(d.left_km, [0, 0, 10, -10], atol=1e-12)
        np.testing.assert_allclose(d.angle_deg, [0, -180, 90, -90], atol=1e-12)
        d2 = sample(range(4), [90, 270, 0, 180]); d2.Cyc = 'CE'
        np.testing.assert_allclose(geom(d2).parallel_km, d.parallel_km)

    def test_rotating_shear_stationary_tilt(self):
        d = sample(range(2))
        for level in ['200', 'surface']:
            d.loc[1, f'clim_{level}_east_ms'] = 0.
            d.loc[1, f'clim_{level}_north_ms'] = .1
        d.loc[1, 'clim_500_east_ms'] = 0.
        d.loc[1, 'clim_500_north_ms'] = .04
        p = lagged_pairs(geom(d), 1).iloc[0]
        self.assertAlmostEqual(p.tilt_rotation_deg_day, 0)
        self.assertAlmostEqual(p.shear_rotation_deg_day, 90)
        self.assertAlmostEqual(p.relative_rotation_deg_day, -90)
        self.assertAlmostEqual(p.fixed_shear_alignment_change_day, 0)

    def test_no_bridging_gaps_regime_changes_or_invalid_days(self):
        for how in ['gap', 'regime', 'invalid']:
            d = sample(range(5))
            if how == 'gap': d = d.drop(index=2)
            if how == 'regime': d.loc[2, 'PV_grad_topo_mag'] = 10.
            if how == 'invalid': d.loc[2, 'TiltDis'] = 0.
            p = lagged_pairs(geom(d), 4)
            self.assertTrue(p.empty, how)

    def test_trailing_no_future_leakage_and_gap_reset(self):
        a = sample(range(12)); b = a.copy(); b.loc[b.Day > 6, 'PV_grad_topo_mag'] = 100.
        x, y = classify_regimes(a), classify_regimes(b)
        self.assertEqual(x.loc[6, 'regime'], y.loc[6, 'regime'])
        g = classify_regimes(sample([0,1,2,3,4,10,11,12,13,14]))
        self.assertEqual(g.loc[5, 'regime'], 'unknown')
        self.assertEqual(g.loc[9, 'regime'], 'planetary')

    def test_actual_future_change_and_wrap(self):
        d = sample(range(2), [359., 1.]); d.loc[1, 'TiltDis'] = 12.
        p = lagged_pairs(geom(d), 1).iloc[0]
        self.assertAlmostEqual(p.growth_km_day, 2.)
        self.assertAlmostEqual(p.tilt_rotation_deg_day, -2.)
        self.assertAlmostEqual(float(wrap(p.relative_rotation_deg_day -
                                         p.tilt_rotation_deg_day + p.shear_rotation_deg_day)), 0.)

    def test_track_keys_and_duplicates(self):
        a = sample(range(2)); b = a.copy(); b.Cyc = 'CE'
        self.assertEqual(len(lagged_pairs(geom(pd.concat([a,b])), 1)), 2)
        with self.assertRaises(ValueError): geom(pd.concat([a,a]))

    def test_weighting_does_not_replace_days_by_eddy_means(self):
        d = pd.DataFrame(dict(Cyc=['AE']*11, Eddy=[1]*10+[2], regime=['planetary']*11,
                             response=[0.]*10+[10.]))
        day = cluster_summary(d, ['response'], min_eddies=2, min_rows=1, n_boot=100)
        eq = cluster_summary(d, ['response'], weighting='eddy', min_eddies=2, min_rows=1, n_boot=100)
        self.assertAlmostEqual(day['mean'].iloc[0], 10/11)
        self.assertAlmostEqual(eq['mean'].iloc[0], 5.)
        self.assertEqual(day.eddies.iloc[0], 2)
        low = cluster_summary(d, ['response'])
        self.assertFalse(low.supported.iloc[0]); self.assertTrue(np.isnan(low.ci_low.iloc[0]))

    def test_zero_topo_valid_missing_ratio_invalid_shallow_column(self):
        d = sample(range(3)); d.loc[0,'PV_grad_topo_mag']=0.
        d.loc[1,'PV_grad_topo_mag']=np.nan; d.loc[2,'water_depth_m']=300.
        g = geom(d)
        self.assertEqual(g.regime.iloc[0], 'planetary')
        self.assertEqual(g.regime.iloc[1], 'unknown')
        self.assertEqual(g.eligible.tolist(), [True,False,False])


if __name__ == '__main__':
    unittest.main()
