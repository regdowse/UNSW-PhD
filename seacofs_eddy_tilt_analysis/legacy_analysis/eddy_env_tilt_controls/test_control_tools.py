"""Scientific and cache-join regression checks using small independent fixtures."""
import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
import control_tools as c


def profile_fixture():
    rows=[]
    for eddy, depths in [(1,[0,100,300]),(2,[0,100])]:
        for z in depths:
            rows.append(dict(Eddy=eddy,Day=10,Depth=z,xc=z/100,yc=0,
                             w=-1e-5*(1+z/100),Omega=-1e-5,q11=1.,q12=0.,q22=4.,Rc=20.,
                             f=-1e-4,TiltDis=3.,Cyc='CE'))
    return c.features(pd.DataFrame(rows))


class ControlTests(unittest.TestCase):
    def test_magnitude_and_ellipse(self):
        d=profile_fixture()
        np.testing.assert_allclose(d.AR,2)
        self.assertAlmostEqual(d.Ro_abs.iloc[0],.1)
        self.assertAlmostEqual(d.rotation_speed_scale.iloc[0],.2)
        bad=d.iloc[:1].copy();bad.q22=-1
        self.assertTrue(c.features(bad).AR.isna().all())

    def test_profile_summaries(self):
        d=c.profile_summaries(profile_fixture()).set_index('Eddy')
        self.assertAlmostEqual(d.loc[1,'abs_zeta_max'],4e-5)
        self.assertEqual(d.loc[1,'abs_zeta_peak_depth_m'],300)
        self.assertAlmostEqual(d.loc[1,'abs_zeta_depth_mean'],2.5e-5)
        self.assertAlmostEqual(d.loc[1,'profile_straightness'],1)
        self.assertEqual(d.loc[1,'profile_bend_km'],0)
        self.assertEqual(d.loc[1,'profile_reversal_fraction'],0)
        self.assertTrue(np.isnan(d.loc[2,'profile_bend_km']))

    def test_missing_levels_not_bridged(self):
        d=profile_fixture();d.loc[1,'abs_zeta']=np.nan
        s=c.profile_summaries(d).iloc[0]
        self.assertTrue(np.isnan(s.abs_zeta_depth_mean))
        self.assertAlmostEqual(s.abs_zeta_max,4e-5)

    def test_matched_days_not_just_matching_day_numbers(self):
        d=profile_fixture()
        available=c.depth_sample(d,[0,100,300],'AR')
        matched=c.depth_sample(d,[0,100,300],'AR',True)
        self.assertEqual(len(available),5)
        self.assertEqual(set(matched.Eddy),{1})
        self.assertEqual(len(matched),3)

    def test_duplicate_merge_rejected(self):
        with self.assertRaises(ValueError):
            c.merge(pd.DataFrame({'Eddy':[1],'Day':[1]}),pd.DataFrame({'Eddy':[1,1],'Day':[1,1]}))

    def test_velocity_rms_and_signed_extrema(self):
        raw=pd.DataFrame(dict(Eddy=[1,1],Day=[1,1],Depth=[0,100],fraction=[1.5,1.5],
                              n_valid=[10,30],coverage=[1,1],w_mean=[1.,3.],w_rms=[2.,4.],
                              w_max=[3.,5.],w_min=[-9.,-1.]))
        _,s=c.vertical_velocity(raw)
        self.assertAlmostEqual(s.vv_mean.iloc[0],2.5)
        self.assertAlmostEqual(s.vv_rms.iloc[0],np.sqrt(13))
        self.assertEqual(s.vv_abs_max.iloc[0],9)
        self.assertTrue(np.isnan(s.vv_deep_minus_upper.iloc[0]))

    def test_bins_include_maximum_and_constant_predictor(self):
        d=pd.DataFrame({'Eddy':range(20),'Cyc':'AE','x':range(20),'TiltDis':range(20)})
        t=c.binned(d,'x',bins=2,min_eddies=1,n_boot=10)
        self.assertEqual(t.eddy_days.sum(),20)
        self.assertAlmostEqual(t['median'].iloc[-1],14.5)
        d.x=1
        self.assertTrue(c.binned(d,'x').empty)

    def test_n2_metric_specific_coverage(self):
        raw=pd.DataFrame({'Eddy':[1],'Day':[1]})
        cache=raw.copy()
        for metric in c.N2_METRICS:
            cache[metric]=2.
            coverage=metric.replace('_core_s2','_core_valid_fraction') if metric.endswith('_core_s2') else metric+'_valid_fraction'
            cache[coverage]=1.
        cache['N2_200m_core_valid_fraction']=.1
        got=c.add_n2(raw,cache)
        self.assertTrue(got.N2_200m_core_s2.isna().all())
        self.assertEqual(got.N2_500m_core_s2.iloc[0],2)

    def test_age_is_elapsed_not_duration(self):
        d=profile_fixture().iloc[:2].copy();d.Day=[10,14];d['Age']=21
        d['ic']=0;d['jc']=0;d['Date']=pd.to_datetime(['2000-01-01','2000-01-05'])
        d=d.drop(columns='Depth')
        out=c.prepare_surface(d,SimpleNamespace(f=np.array([[-1e-4]])))
        np.testing.assert_allclose(out.age_fraction,[0,1])
        np.testing.assert_allclose(out.age_days,[0,4])

    def test_adjusted_model_recovers_known_effect(self):
        rng=np.random.default_rng(72);n=240
        x=rng.normal(size=n)
        d=pd.DataFrame(dict(Eddy=np.repeat(np.arange(80),3),Day=np.tile(np.arange(3),80),
                            Cyc='AE',x=x,TiltDis=np.exp(2+.2*x+rng.normal(0,.03,n))-1))
        fit=c.adjusted(d,['x']).iloc[0]
        self.assertEqual(fit.status,'ok')
        self.assertAlmostEqual(fit.slope,.2*d.x.std(),delta=.015)
        self.assertLess(fit.ci_low,fit.slope)
        self.assertGreater(fit.ci_high,fit.slope)

if __name__=='__main__':unittest.main()
