import unittest
import numpy as np
import pandas as pd
import obs_clim_tools as t


def profiles(eddy=1, day=1, depth=(0, 50, 200, 400, 660), dx=.03, dy=.01):
    return pd.DataFrame([dict(Eddy=eddy, Day=day, Depth=z, xc=10+dx*z,
                              yc=20+dy*z, w=-5e-5) for z in depth])


class ComparisonTests(unittest.TestCase):
    def test_physical_ro_conversion(self):
        o = t.observations()
        np.testing.assert_allclose(o['zeta_s-1']/o['f_s-1'], [0.74, 1.30])
        self.assertTrue((o['zeta_s-1'] < 0).all())

    def test_direction_depth_and_rotation(self):
        r = t.profile_metrics(profiles(), rotation_rad=np.pi/2).iloc[0]
        self.assertTrue(r.profile_ok)
        self.assertEqual(r.fit_top_m, 50)
        self.assertEqual(r.fit_bottom_m, 660)
        self.assertAlmostEqual(r.rate_km_per_100m, np.sqrt(10))
        # Increasing depth (+x,+y) rotated 90 deg -> (-east,+north).
        # Deep-to-shallow therefore (+east,-north), SE.
        self.assertAlmostEqual(r.bearing_deg, np.degrees(np.arctan2(1,-3)))
        self.assertAlmostEqual(r.fit_r2, 1)

    def test_no_shallow_or_remote_endpoint(self):
        for depths in [(0,50,200,599), (0,50,200,750)]:
            self.assertFalse(t.profile_metrics(profiles(depth=depths)).iloc[0].profile_ok)

    def test_peak_cap_and_sign(self):
        p = profiles()
        p.loc[p.Depth.eq(660), 'w'] = -1
        p.loc[p.Depth.eq(400), 'w'] = -9e-5
        r = t.profile_metrics(p).iloc[0]
        self.assertEqual(r.w_peak, -9e-5)
        self.assertEqual(r.w_peak_depth_m, 400)

    def test_join_both_keys_and_no_direction_selection(self):
        o = t.observations().iloc[:1].copy()
        p = profiles(); p['w'] = o['zeta_s-1'].iloc[0]
        s = pd.DataFrame([dict(Eddy=e,Day=1,Cyc=c,lat=o.lat.iloc[0],lon=o.lon.iloc[0],
                               TiltDis=10,TiltDir=90) for e,c in [(1,'CE'),(2,'CE'),(3,'AE')]])
        g,a = t.match(s,t.profile_metrics(p),o,lambda *a: np.zeros(len(a[0])))
        self.assertEqual(g['Murphy'].Eddy.tolist(), [1])
        self.assertEqual(a.iloc[-1].days, 1)
        self.assertEqual(len(t.summary(g,o)), 3)

    def test_duplicates_fail(self):
        p=profiles()
        with self.assertRaises(ValueError):
            t.profile_metrics(pd.concat([p,p.iloc[:1]]))

    def test_zero_tilt_and_wrapping(self):
        r=t.profile_metrics(profiles(dx=0,dy=0)).iloc[0]
        self.assertEqual(r.rate_km_per_100m, 0)
        self.assertTrue(np.isnan(r.bearing_deg))
        self.assertEqual(t.angle_difference(1,359), 2)

    def test_empty_candidates(self):
        o=t.observations()
        p=profiles().iloc[:0]
        s=pd.DataFrame([dict(Eddy=1,Day=1,Cyc='CE',lat=-28,lon=154,TiltDis=2,TiltDir=20)])
        g,a=t.match(s,t.profile_metrics(p),o,lambda *a: np.zeros(len(a[0])))
        self.assertTrue(g['Murphy'].empty)
        self.assertTrue(t.summary(g,o).model_median.isna().all())
        fig=t.plot_comparison(g,o)
        import matplotlib.pyplot as plt
        plt.close(fig)


if __name__ == '__main__':
    unittest.main()
