"""Synthetic invariants for the isolated temporal-weight experiment."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'seacofs_eddy_dataset_modular' / 'src'))
from seacofs_eddy_dataset.core.tilt import compute_weighted_tilt
from delta_sensitivity_tools import increments, fit_snapshot, lifetime, maximum_span


def fixture():
    rows = []
    for day in range(20):
        z = np.arange(0, (800 if day % 3 else 600) + 1, 20.)
        rows.append(pd.DataFrame({'Eddy': 1, 'Day': day, 'Depth': z,
            'xc': 100+day+z*(.003+.001*np.sin(day)),
            'yc': 400-day+z*.01+np.sin(z/90+day)*.2}))
    return pd.concat(rows, ignore_index=True)


class SensitivityTests(unittest.TestCase):
    def test_equal_matches_production_with_missing_days_and_depths(self):
        p = fixture()
        p = p.loc[~p.Day.isin([5, 8])]
        table, _, _ = lifetime(p)
        baseline = compute_weighted_tilt(p, 1)
        joined = table.merge(baseline, on='Day')
        np.testing.assert_allclose(joined['Equal'], joined.TiltDis, rtol=1e-9, atol=1e-9, equal_nan=True)

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

    def test_span_is_pairwise_not_surface_relative(self):
        p = pd.DataFrame({'Depth':[0,100,200], 'xc':[0,-3,4], 'yc':[0,0,0]})
        self.assertEqual(maximum_span(p), 7.)


if __name__ == '__main__':
    unittest.main()
