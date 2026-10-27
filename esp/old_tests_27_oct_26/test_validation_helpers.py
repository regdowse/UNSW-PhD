"""Fast scientific/diagnostic checks; no Katana data or backend required."""
import types
import unittest

import numpy as np
import pandas as pd

import validation_helpers as h


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.cases = h.make_cases()

    def test_geometry_and_units(self):
        q = h.shape_matrix(1.7, 35)
        self.assertAlmostEqual(np.linalg.det(q), 1)
        a, theta = h.geometry(q)
        self.assertAlmostEqual(a, 1.7)
        self.assertAlmostEqual(theta, 35)
        case = self.cases['SOLO']
        rc, om = case['rc'], case['omega']
        uv = h.velocity(np.array([[rc / np.sqrt(2), 0]]), [0, 0], om, rc, np.eye(2))
        self.assertAlmostEqual(abs(uv[0, 1]), case['speed_ref'])
        self.assertLess(uv[0, 1], 0)  # Southern Hemisphere cyclone

    def test_core_vorticity_and_divergence(self):
        case = self.cases['LATTE']
        d = 1.
        points = np.array([[d, 0], [-d, 0], [0, d], [0, -d]])
        u = h.velocity(points, [0, 0], case['omega'], case['rc'], case['q'])
        curl = (u[0, 1] - u[1, 1] - u[2, 0] + u[3, 0]) / (2 * d)
        divergence = (u[0, 0] - u[1, 0] + u[2, 1] - u[3, 1]) / (2 * d)
        self.assertAlmostEqual(curl, case['vorticity'], places=12)
        self.assertAlmostEqual(divergence, 0, places=12)

    def test_axial_wrap(self):
        self.assertEqual(h.axial_error(179, 1), -2)
        self.assertEqual(h.axial_error(5, 185), 0)

    def test_spatial_covariance(self):
        xy = np.array([[0., 0.], [1000., 0.], [0., 3000.]])
        factor = h.noise_factor(xy, 2000.)
        covariance = factor @ factor.T
        np.testing.assert_allclose(np.diag(covariance), 1)
        self.assertAlmostEqual(covariance[0, 1], np.exp(-.125), places=9)

    def test_shared_intersection(self):
        case = self.cases['DOPPIO']
        i, j = h.inner_indices(case, 30000)
        self.assertEqual(len(np.intersect1d(i, j)), 1)
        np.testing.assert_equal(case['xy'][np.intersect1d(i, j)][0], [10000., 10000.])

    def test_invalid_geometry_is_not_success(self):
        case = self.cases['LATTE']
        backend = types.SimpleNamespace(latte=lambda *a: (0., 0., 1., np.diag([1., -1.]), 1., 1.))
        row = h.fit_case(backend, case, case['uv'], 30000)
        self.assertFalse(row['valid'])
        self.assertEqual(row['status'], 'invalid_ellipse')

    def test_fallback_is_not_success(self):
        case = self.cases['LATTE']
        backend = types.SimpleNamespace(
            latte=lambda *a: (0., 0., case['vorticity'], case['q'], case['omega'], 1.),
            out_core_param_fit=lambda *a, **k: (k['Rc0'], -.5 * k['Omega0'] * k['Rc0']**2, k['Omega0']))
        row = h.fit_case(backend, case, case['uv'], 30000)
        self.assertFalse(row['valid'])
        self.assertEqual(row['status'], 'outer_seed_return')

    def test_outer_refinement_used_for_vorticity(self):
        case = self.cases['LATTE']
        backend = types.SimpleNamespace(
            latte=lambda *a: (0., 0., 0., case['q'], .9 * case['omega'], 1.),
            out_core_param_fit=lambda *a, **k: (case['rc'], -.5 * case['omega'] * case['rc']**2, case['omega']))
        row = h.fit_case(backend, case, case['uv'], 30000)
        self.assertTrue(row['valid'])
        self.assertAlmostEqual(row['vorticity_error_pct'], 0)
        self.assertAlmostEqual(row['angle_error_deg'], 0)

    def test_failures_remain_in_denominator(self):
        rows = pd.DataFrame([
            dict(experiment='noise', method='LATTE', noise_kind='independent', level=.05,
                 valid=True, centre_error_pct=1.),
            dict(experiment='noise', method='LATTE', noise_kind='independent', level=.05,
                 valid=False, centre_error_pct=np.nan)])
        row = h.summarise(rows).iloc[0]
        self.assertEqual(row.n_attempted, 2)
        self.assertEqual(row.valid_fraction, .5)
        self.assertEqual(row.centre_error_pct_median, 1.)


if __name__ == '__main__':
    unittest.main()
