import unittest
import numpy as np
import departure_helpers as dh
import validation_helpers as vh

class DepartureTests(unittest.TestCase):
    def setUp(self):self.c=vh.make_cases()['LATTE']
    def test_gaussian_limit(self):
        np.testing.assert_allclose(dh.vortex_velocity(self.c['xy'],self.c),self.c['uv'])
        np.testing.assert_allclose(dh.vortex_velocity(self.c['xy'],self.c,1e6),self.c['uv'],atol=1e-5)
    def test_peak_location_and_core_gradient(self):
        c=self.c
        for p in [np.inf,8,4,2]:
            circular=dict(c,q=np.eye(2))
            r=c['rc']/np.sqrt(2)*np.array([.999,1,1.001])
            u=dh.vortex_velocity(np.c_[r,np.zeros(3)],circular,p)
            self.assertEqual(np.argmax(np.linalg.norm(u,axis=1)),1)
            near=dh.vortex_velocity(np.array([[.1,0.]]),circular,p)
            self.assertAlmostEqual(near[0,1]/.1,c['omega'],places=12)
    def test_divergence_and_shear_curl(self):
        c=self.c;centre=np.array([14000.,-7000.]);d=.1
        xy=centre+np.array([[d,0],[-d,0],[0,d],[0,-d]])
        for angle in [0,35,90]:
            b=dh.background_velocity(xy,c,'shear',.5,angle)
            divergence=(b[0,0]-b[1,0]+b[2,1]-b[3,1])/(2*d)
            curl=(b[0,1]-b[1,1]-b[2,0]+b[3,0])/(2*d)
            self.assertAlmostEqual(divergence,0,places=11)
            self.assertAlmostEqual(curl,-.5*c['speed_ref']/c['rc'],places=11)
        for p in [8,4,2]:
            u=dh.vortex_velocity(xy,c,p)
            self.assertAlmostEqual((u[0,0]-u[1,0]+u[2,1]-u[3,1])/(2*d),0,places=10)
    def test_constant_background_and_zero_control(self):
        c=self.c
        b=dh.background_velocity(c['xy'],c,'current',.1,90)
        np.testing.assert_allclose(b[:,0],0,atol=1e-14)
        np.testing.assert_allclose(b[:,1],.1*c['speed_ref'])
        for kind in ['current','shear']:
            np.testing.assert_array_equal(dh.background_velocity(c['xy'],c,kind,0,35),0)
    def test_evaluation_excludes_observations(self):
        for c in vh.make_cases().values():
            xy=dh.evaluation_points(c)
            distance2=np.sum((xy[:,None]-c['xy'][None,:])**2,axis=2)
            self.assertGreater(distance2.min(),1e-6)

if __name__=='__main__':unittest.main()
