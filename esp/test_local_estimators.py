import unittest
from unittest.mock import patch
import numpy as np
import local_estimators as le
import validation_helpers as vh

class LocalEstimatorTests(unittest.TestCase):
    def setUp(self):
        self.case=vh.make_cases()['SOLO']

    def test_offset_initialisation(self):
        c=self.case;xy=c['xy'];uv=c['uv']
        x,y,_,q,om=le.solo_local_initialisation(xy[:,0],uv[:,0],uv[:,1],
            centre_guess=10000.,search_half_width=60000.)
        self.assertLess(abs(x),100.)
        self.assertLess(abs(y+20000.),100.)
        np.testing.assert_equal(q,np.eye(2))
        self.assertLess(om,0)

    def test_no_crossing_is_explicit(self):
        with self.assertRaisesRegex(le.FitFailure,'no_supported_crossing'):
            le.solo_local_initialisation(np.arange(10.)*1000,np.ones(10),np.ones(10),
                centre_guess=4000,search_half_width=5000)

    def test_outer_matches_exact_truth(self):
        c=self.case;xy=c['xy'];uv=c['uv']
        fit=le.outer_fit_diagnostic(*xy.T,*uv.T,0,0,c['q'],Omega0=.9*c['omega'])
        self.assertEqual(fit['status'],'ok')
        self.assertAlmostEqual(fit['Rc_m']/c['rc'],1,places=6)
        self.assertAlmostEqual(fit['Omega']/c['omega'],1,places=6)

    def test_radius_cap_retains_attempt(self):
        c=self.case;xy=c['xy'];uv=c['uv']
        fit=le.outer_fit_diagnostic(*xy.T,*uv.T,0,0,c['q'],Omega0=c['omega'],Rc_max=80000.)
        self.assertEqual(fit['status'],'outer_radius_limit')
        self.assertGreater(fit['Rc_m'],80000.)

    def test_optimiser_failure_is_explicit(self):
        c=self.case;xy=c['xy'];uv=c['uv']
        with patch.object(le,'curve_fit',side_effect=RuntimeError('test failure')):
            fit=le.outer_fit_diagnostic(*xy.T,*uv.T,0,0,c['q'],Omega0=c['omega'])
        self.assertEqual(fit['status'],'outer_optimisation_failed')
        self.assertTrue(np.isnan(fit['Rc_m']))

    def test_original_and_local_get_identical_noise(self):
        captured=[]
        def fit(backend,case,uv,*a,**k):
            captured.append(uv.copy())
            return dict(valid=True,status='ok')
        settings=dict(experiment='noise',levels=[0,.05],repeats=2,seed=17)
        with patch.object(vh,'fit_case',side_effect=fit):
            vh.run_experiment(None,{'SOLO':self.case},**settings)
            original=captured.copy();captured.clear()
            vh.run_experiment(None,{'SOLO':self.case},solo_local=True,**settings)
        self.assertEqual(len(original),len(captured))
        for a,b in zip(original,captured):np.testing.assert_array_equal(a,b)

if __name__=='__main__':unittest.main()
