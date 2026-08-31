"""Independent checks of probability weighting, pairing and sparse samples."""
import unittest
import numpy as np
import pandas as pd
import publication_tools as pub


class PublicationTests(unittest.TestCase):
    def test_weighting_changes_estimand_not_bootstrap_unit(self):
        sums=np.array([[1.,0.],[0.,9.]])
        counts=np.array([1.,9.])
        equal=pub.cluster_interval(sums, counts, min_eddies=2, n_boot=50)
        obs=pub.cluster_interval(sums, counts, equal_eddy=False, min_eddies=2, n_boot=50)
        np.testing.assert_allclose(equal[0],[.5,.5])
        np.testing.assert_allclose(obs[0],[.1,.9])
        for estimate,low,high in [equal,obs]:
            self.assertAlmostEqual(estimate.sum(),1.)
            self.assertTrue(np.all(low>=0) and np.all(high<=1))

    def test_histogram_endpoints_and_equal_total_mass(self):
        g=pd.DataFrame({'Eddy':[1]+[2]*9,'AlignmentDeg':[0]+[90]*9})
        h=pub.histogram_estimate(g,np.arange(0,91,5),min_eddies=2,n_boot=50)
        self.assertAlmostEqual(h.estimate.sum(),100.)
        self.assertAlmostEqual(h.estimate.iloc[0],50.)
        self.assertAlmostEqual(h.estimate.iloc[-1],50.)

    def test_day_support_applies_inside_class(self):
        frame=pd.DataFrame({'Eddy':[1]*10+[2]*5,'AxisRatio':[1.2]*9+[2.2]+[2.2]*5})
        part=pub.eligible(frame.loc[frame.AxisRatio>2],min_days=5)
        self.assertEqual(part.Eddy.unique().tolist(),[2])

    def test_paired_depth_differences_and_sparse_plot(self):
        records=[]
        for eddy in range(4):
            for day in range(5):
                for depth in [0,100]:
                    q11,q22=(.25,1) if depth==0 else (1,.25)
                    records.append(dict(Eddy=eddy,Day=day,Cyc='AE',ShapeDepth=depth,
                                        TiltDis=20,TiltDir=110,q11=q11,q12=0,q22=q22))
        import ellipse_tilt_tools as et
        g=et.ellipse_geometry(pd.DataFrame(records))
        g['AlignmentDeg']=np.abs(et.axial_difference(g.TiltDir,g.MajorBearing))
        g['AlignmentCos2']=np.cos(np.deg2rad(2*g.AlignmentDeg))
        r=pub.direction_data(g,depths=(0,100),n_boot=30,min_eddies=3)
        ae=r['depth_contrasts'].query("Cyc=='AE'").iloc[0]
        self.assertAlmostEqual(ae.estimate,-2.)
        self.assertAlmostEqual(ae.low,-2.)
        self.assertAlmostEqual(ae.high,-2.)
        self.assertEqual(ae.eddies,4)
        self.assertTrue(r['depth'].query("Cyc=='CE'").estimate.isna().all())
        fig=pub.plot_direction(r)
        self.assertEqual(len(fig.axes),6)
        import matplotlib.pyplot as plt
        plt.close(fig)


if __name__=='__main__':
    unittest.main()
