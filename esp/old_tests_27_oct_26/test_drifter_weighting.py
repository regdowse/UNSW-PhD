import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import drifter_weighting as dw
import platform_weighting as pw
import validation_helpers as vh


def field():
    rng=np.random.default_rng(22)
    xy=rng.uniform(-80000,80000,(150,2))
    q=vh.shape_matrix(1.2,20)
    uv=vh.velocity(xy,[0,0],-2e-5,75000,q)
    data=pd.DataFrame(dict(x=xy[:,0],y=xy[:,1],u=uv[:,0],v=uv[:,1],id=np.arange(150),
        platform=np.where(np.arange(150)<75,'Drifter','Satellite'),
        drifter_id=np.where(np.arange(150)<75,'1',''),time=pd.Timestamp('2023-10-29 06:00'),depth_m=15.))
    return data,pw.checked_inner(0,0,q,-2e-5)


class DrifterTests(unittest.TestCase):
    def test_reference_label_preserves_numerics(self):
        d,_=field();adcp=d.assign(platform=d.platform.replace({'Drifter':'ADCP'}))
        a=pw.weighted_inner(adcp,.3);b=pw.weighted_inner(d,.3,reference_platform='Drifter')
        np.testing.assert_allclose(a['Q'],b['Q'])
        self.assertEqual(a['xc'],b['xc'])
        _,inner=field()
        oa=pw.weighted_outer(adcp,inner,.3,1e6)
        ob=pw.weighted_outer(d,inner,.3,1e6,reference_platform='Drifter')
        np.testing.assert_allclose([oa['Rc'],oa['Omega']],[ob['Rc'],ob['Omega']])
        np.testing.assert_allclose(pw.platform_weights(adcp.platform,.3),
            pw.platform_weights(d.platform,.3,reference_platform='Drifter'))

    def test_fixed_masks_and_independent_outer_weight(self):
        d,inner=field();core=d[pw.rho2(d,inner)<35000**2].id.tolist()
        ref=dw.finish_fit(d,inner,core,d.id.tolist(),'reference',original=True)
        a=dw.weighted_fit(d,ref,'a',0,.1);b=dw.weighted_fit(d,ref,'b',0,1)
        self.assertTrue(a['summary']['valid']);self.assertTrue(b['summary']['valid'])
        self.assertEqual(a['inner_ids'],core);self.assertEqual(b['outer_ids'],ref['outer_ids'])
        np.testing.assert_allclose(a['inner']['Q'],b['inner']['Q'])
        self.assertEqual(a['inner']['xc'],b['inner']['xc'])

    def test_inner_only_has_no_full_profile(self):
        d,inner=field();f=dw.finish_fit(d,inner,d.id.tolist(),[],'inner',original=True,inner_only=True)
        self.assertTrue(f['summary']['inner_valid']);self.assertFalse(f['summary']['valid'])
        self.assertTrue(np.isnan(f['summary']['Rc_km']))
        self.assertTrue(pw.residuals([f],d).empty)
        self.assertEqual(len(dw.inner_residuals([f],d)),2)

    def test_missing_drifter_cannot_be_drifter_led(self):
        d,inner=field();d=d[d.platform=='Satellite']
        ref=dw.finish_fit(d,inner,d.id.tolist(),d.id.tolist(),'reference',original=True)
        f=dw.weighted_fit(d,ref,'weighted',0,1)
        self.assertFalse(f['summary']['inner_valid'])

    def test_snapshot_window_is_inclusive_and_no_extrapolation(self):
        t=pd.Timestamp('2023-10-29 06:00');d,_=field();d=d.iloc[:4].copy()
        d['time']=[t-pd.Timedelta(hours=h) for h in [25,24,-24,-25]]
        sat=dict(daily=pd.DatetimeIndex([t]),time=pd.DatetimeIndex([t-pd.Timedelta(hours=12),t+pd.Timedelta(hours=12)]),
                 uv=np.array([[[[0.,0.]]],[[[2.,4.]]]]),x=np.array([0.]),y=np.array([0.]))
        s=dw.snapshot(d,sat,t)
        self.assertEqual((s.platform=='Drifter').sum(),2)
        np.testing.assert_allclose(s[s.platform=='Satellite'][['u','v']],[ [1.,2.] ])
        with self.assertRaises(ValueError):dw.snapshot(d,sat,t+pd.Timedelta(days=1))

    def test_trajectory_second_pass_uses_previous_q(self):
        t=pd.date_range('2023-10-19 06:00',periods=2)
        data=pd.DataFrame(dict(x=[0.,3000.,6000.],y=0.,u=0.,v=0.,id=[0,1,2],platform='Satellite'))
        calls=[]
        def inner(backend,d):
            calls.append(d.id.tolist())
            q=np.diag([4.,.25]) if len(calls)==3 else np.diag([.25,4.])
            return pw.checked_inner(0,0,q,-1e-5)
        with patch.object(dw,'snapshot',return_value=data),patch.object(pw,'backend_inner',side_effect=inner):
            dw.track_to_day(None,None,{'daily':t},t[-1],core_radius=5000)
        self.assertEqual(calls[0],[0,1])
        self.assertEqual(calls[-1],[0,1,2])

    def test_drifter_only_outer_fit_does_not_prove_peak_support(self):
        d,inner=field()
        d=d[(d.platform=='Drifter')&(pw.rho2(d,inner)<30000**2)]
        result=pw.weighted_outer(d,inner,0,1e6,reference_platform='Drifter')
        self.assertEqual(result['status'],'ok')
        self.assertFalse(result['peak_bracketed'])

    def test_reference_coverage_counts(self):
        d,inner=field();ref=dw.finish_fit(d,inner,d.id.tolist(),d.id.tolist(),'reference',original=True)
        c=dw.coverage(d,ref)
        np.testing.assert_array_equal(c.n_below_peak+c.n_above_peak,c.n)
        self.assertEqual(c.n.sum(),len(d))

if __name__=='__main__':unittest.main()
