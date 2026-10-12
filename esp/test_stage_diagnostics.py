import unittest
from unittest.mock import patch
import types
import numpy as np
import validation_helpers as vh
import departure_helpers as dh

class StageDiagnosticsTests(unittest.TestCase):
    def test_outer_failure_preserves_inner_scores(self):
        c=vh.make_cases()['LATTE']
        backend=types.SimpleNamespace(latte=lambda *a:(0.,0.,c['vorticity'],c['q'],c['omega'],1.))
        for status in ('outer_radius_limit','outer_optimisation_failed'):
            with patch('local_estimators.outer_fit_diagnostic', return_value=dict(status=status,Rc_m=1.1e6,Omega=c['omega'],psi0=1.)):
                r=vh.fit_case(backend,c,c['uv'],30000.,diagnostic_outer=True)
            self.assertTrue(r['inner_valid']);self.assertFalse(r['valid'])
            self.assertAlmostEqual(r['centre_error_km'],0)
            self.assertAlmostEqual(r['alpha_error_pct'],0)
            self.assertNotIn('vorticity_error_pct',r)
            self.assertEqual(r['status'],status)
    def test_inner_failure_has_no_inner_scores(self):
        c=vh.make_cases()['LATTE']
        backend=types.SimpleNamespace(latte=lambda *a:(0.,0.,1.,np.diag([1.,-1.]),1.,1.))
        r=vh.fit_case(backend,c,c['uv'],30000.,diagnostic_outer=True)
        self.assertFalse(r['inner_valid']);self.assertFalse(r['valid'])
        self.assertNotIn('centre_error_km',r)
    def test_plot_uses_stage_specific_validity(self):
        import pandas as pd
        import matplotlib.pyplot as plt
        r=pd.DataFrame([dict(test='profile',method='LATTE',departure=0.,inner_valid=True,valid=False,centre_error_km=2.,peak_radius_error_pct=20.,eddy_velocity_error_pct=30.)])
        f=dh.plot_test(r)
        self.assertEqual(f.axes[0].lines[0].get_ydata()[0],2.)
        self.assertTrue(np.isnan(f.axes[1].lines[0].get_ydata()[0]))
        self.assertTrue(np.isnan(f.axes[2].lines[0].get_ydata()[0]))
        plt.close(f)
if __name__=='__main__':unittest.main()
