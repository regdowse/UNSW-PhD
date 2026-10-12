"""Focused objective/geometry checks; no external backend or field files required."""
import unittest
import numpy as np
import pandas as pd
import platform_weighting as p
import validation_helpers as v


def observations(gaussian=False):
    rng=np.random.default_rng(8)
    xy=rng.uniform(-85000,85000,(150,2))
    q=v.shape_matrix(1.35,27)
    centre=np.array([2500.,-1800.]);om=-2e-5
    if gaussian:
        uv=v.velocity(xy,centre,om,75000.,q)
    else:
        qr=(xy-centre)@q
        uv=om*np.column_stack((-qr[:,1],qr[:,0]))
    d=pd.DataFrame(dict(x=xy[:,0],y=xy[:,1],u=uv[:,0],v=uv[:,1],
        platform=np.where(np.arange(150)<110,'ADCP','Satellite'),id=np.arange(150)))
    return d,p.checked_inner(*centre,q,om)


class WeightingTests(unittest.TestCase):
    def test_stage_mass(self):
        d,_=observations();w=p.platform_weights(d.platform,.3)
        self.assertAlmostEqual(w[d.platform=='ADCP'].sum(),1)
        self.assertAlmostEqual(w[d.platform=='Satellite'].sum(),.3)
        for bad in [-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):p.platform_weights(d.platform,bad)

    def test_exact_inner_recovery(self):
        d,truth=observations()
        for weight in [0,.1,1,10]:
            fit=p.weighted_inner(d,weight)
            np.testing.assert_allclose([fit['xc'],fit['yc']],[truth['xc'],truth['yc']],atol=1e-7)
            np.testing.assert_allclose(fit['Q'],truth['Q'],atol=1e-10)
            self.assertAlmostEqual(fit['Omega'],truth['Omega'],12)

    def test_zero_inner_weight_ignores_satellite(self):
        d,_=observations();expected=p.weighted_inner(d,0)
        d.loc[d.platform=='Satellite',['u','v']]+=20
        actual=p.weighted_inner(d,0)
        np.testing.assert_allclose(actual['Q'],expected['Q'])
        self.assertAlmostEqual(actual['xc'],expected['xc'])

    def test_replication_does_not_change_platform_objective(self):
        d,_=observations();d.loc[d.platform=='Satellite','u']+=.08
        duplicate=pd.concat([d,d[d.platform=='ADCP'],d[d.platform=='ADCP']],ignore_index=True)
        a,b=p.weighted_inner(d,.5),p.weighted_inner(duplicate,.5)
        np.testing.assert_allclose(a['Q'],b['Q'],atol=1e-10)
        np.testing.assert_allclose([a['xc'],a['yc']],[b['xc'],b['yc']],atol=1e-6)

    def test_outer_replication_and_zero_weight(self):
        d,inner=observations(True)
        d.loc[d.platform=='Satellite','u']*=.8
        duplicate=pd.concat([d,d[d.platform=='ADCP']],ignore_index=True)
        a=p.weighted_outer(d,inner,.5,1e6);b=p.weighted_outer(duplicate,inner,.5,1e6)
        self.assertEqual(a['status'],'ok')
        np.testing.assert_allclose([a['Rc'],a['Omega']],[b['Rc'],b['Omega']],rtol=1e-5)
        zero=p.weighted_outer(d,inner,0,1e6)
        only=p.weighted_outer(d[d.platform=='ADCP'],inner,0,1e6)
        np.testing.assert_allclose([zero['Rc'],zero['Omega']],[only['Rc'],only['Omega']],rtol=1e-10)
        self.assertTrue((zero['kept'].platform=='ADCP').all())

    def test_equal_point_outer_equivalence(self):
        d,inner=observations(True)
        d.loc[d.platform=='Satellite','u']*=.85
        kept,_,_=p.outer_samples(d,inner)
        ratio=(kept.platform=='Satellite').sum()/(kept.platform=='ADCP').sum()
        a=p.weighted_outer(d,inner,ratio,1e6)
        b=p.weighted_outer(d,inner,1,1e6,per_observation=True)
        np.testing.assert_allclose([a['Rc'],a['Omega']],[b['Rc'],b['Omega']],rtol=1e-10)

    def test_degenerate_data_rejected(self):
        d,_=observations();d[['x','y']]=0
        with self.assertRaises(p.FitError):p.weighted_inner(d,1)
        d,_=observations();d['u']=d.x*1e-5;d['v']=-d.y*1e-5
        with self.assertRaises(p.FitError):p.weighted_inner(d,1)

    def test_temporal_interpolation_no_extrapolation(self):
        t=pd.to_datetime(['2023-10-19','2023-10-21']);a=np.array([[0.,np.nan],[4.,5.]])
        np.testing.assert_allclose(p.interpolate_field(t,a,'2023-10-20'),[2,np.nan])
        np.testing.assert_allclose(p.interpolate_field(t,a,'2023-10-19'),a[0])
        with self.assertRaises(ValueError):p.interpolate_field(t,a,'2023-10-18')

    def test_outer_limit_is_explicit(self):
        d,inner=observations(True);r=p.weighted_outer(d,inner,1,10000)
        self.assertEqual(r['status'],'outer_radius_limit')
        self.assertGreater(r['Rc'],10000)

if __name__=='__main__':unittest.main()
