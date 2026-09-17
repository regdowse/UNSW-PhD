"""Statistics, units, orientation and full-fit adapter checks with analytic fields."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'seacofs_eddy_dataset_modular'/'src'))
import composite_breakdown_tools as cbt


def model(x,y,xc,yc,Q,omega,rc):
    dx=x-xc;dy=y-yc
    f=omega*np.exp(-(Q[0,0]*dx**2+2*Q[0,1]*dx*dy+Q[1,1]*dy**2)/rc**2)
    return -f*(Q[0,1]*dx+Q[1,1]*dy),f*(Q[0,0]*dx+Q[0,1]*dy)


def analytic_backend(xc,yc,Q,omega_si,rc):
    """Known inner solution plus a numerical outer fit, to check adapter plumbing.

    The installed Katana ESP implementation is not available in local tests.
    """
    calls=[]
    def doppio(x1,y1,u1,v1,x2,y2,u2,v2):
        assert np.allclose(y1,y1[0]) and np.allclose(x2,x2[0])
        raw=omega_si*1000.
        expected=model(x1,y1,xc,yc,Q,raw,rc)
        np.testing.assert_allclose(u1,expected[0]);np.testing.assert_allclose(v1,expected[1])
        expected=model(x2,y2,xc,yc,Q,raw,rc)
        np.testing.assert_allclose(u2,expected[0]);np.testing.assert_allclose(v2,expected[1])
        calls.append((x2[0],y1[0]))
        return xc,yc,raw*np.trace(Q),Q.copy(),raw
    def outer(x,y,u,v,xc,yc,Q,w):
        def residual(p):
            uh,vh=model(x,y,xc,yc,Q,p[0],np.exp(p[1]))
            return np.r_[uh-u,vh-v]
        fit=least_squares(residual,[w/np.trace(Q),np.log(40.)],xtol=1e-12,ftol=1e-12,gtol=1e-12)
        om,r=fit.x[0],np.exp(fit.x[1])
        return r,-.5*om*r*r,om
    return SimpleNamespace(doppio=doppio,out_core_param_fit=outer,model_uv_at_xy=model,calls=calls)


class CentreTests(unittest.TestCase):
    def test_day_weighting_variance_covariance_and_counts(self):
        members=pd.DataFrame(dict(Eddy=[1,1,1,2],Day=[1,2,3,1],Depth=[200]*4,
                                  xc=[0.,2.,4.,10.],yc=[0.,4.,8.,20.]))
        result=cbt.centre_statistics(members,[200],n_boot=100).iloc[0]
        self.assertAlmostEqual(result.mean_xc,4.)
        self.assertAlmostEqual(result.var_xc,np.var([0,2,4,10],ddof=1))
        self.assertAlmostEqual(result.cov_xy,2*result.var_xc)
        self.assertEqual(result.n_eddies,2)
        self.assertAlmostEqual(result.mean_tilt_x,-4.)
        self.assertLessEqual(result.mean_xc_ci_low,result.mean_xc)
        self.assertGreaterEqual(result.mean_xc_ci_high,result.mean_xc)

    def test_single_member_and_zero_reference(self):
        members=pd.DataFrame(dict(Eddy=[1,2,1],Day=[1,1,1],Depth=[0,0,200],xc=[0.,0.,5.],yc=[0.,0.,2.]))
        result=cbt.centre_statistics(members,[0,200],n_boot=100)
        self.assertEqual(result.var_xc.iloc[0],0)
        self.assertEqual(result.mean_xc_ci_low.iloc[0],0)
        self.assertTrue(np.isnan(result.var_xc.iloc[1]))
        self.assertTrue(np.isnan(result.mean_xc_ci_low.iloc[1]))

    def test_duplicate_members_rejected(self):
        members=pd.DataFrame(dict(Eddy=[1,1],Day=[1,1],Depth=[0,0],xc=[0,0],yc=[0,0]))
        with self.assertRaises(ValueError):cbt.centre_statistics(members,[0],n_boot=100)

    def test_full_fit_units_and_non_square_xy_grid(self):
        x=np.arange(-180,181,5.);y=np.arange(-150,151,5.)
        X,Y=np.meshgrid(x,y);Q=np.array([[1.2,.2],[.2,1.04/1.2]])
        backend=analytic_backend(10.,-5.,Q,-1e-5,60.)
        u,v=model(X*1000,Y*1000,10000.,-5000.,Q,-1e-5,60000.)
        result,uf,vf=cbt.fit_composite(X,Y,u[...,None],v[...,None],[0.],np.full((*X.shape,1),5),backend,expected_sign=-1)
        row=result.iloc[0]
        self.assertTrue(row.fit_ok,msg=row.reason)
        self.assertAlmostEqual(row.Omega,-1e-5,places=10)
        self.assertAlmostEqual(row.Rc,60.,places=5)
        self.assertAlmostEqual(row.psi0,18000.,places=3)
        self.assertAlmostEqual(row.w,-1e-5*np.trace(Q),places=10)
        self.assertAlmostEqual(row.xc,10.)
        np.testing.assert_allclose(uf[...,0],u,atol=1e-8)
        np.testing.assert_allclose(vf[...,0],v,atol=1e-8)
        self.assertGreater(row.vector_r2,.999999)

    def test_tracking_uses_last_success_and_failure_rows(self):
        x=np.arange(-200,201,10.);X,Y=np.meshgrid(x,x)
        Q=np.eye(2);backend=analytic_backend(10.,-10.,Q,1e-5,60.)
        u,v=model(X*1000,Y*1000,10000.,-10000.,Q,1e-5,60000.)
        uu=np.stack([u,np.zeros_like(u),u],axis=-1);vv=np.stack([v,np.zeros_like(v),v],axis=-1)
        # Reject unsupported middle depth before calling the backend.
        counts=np.ones_like(uu,int);counts[...,1]=0
        result,uf,_=cbt.fit_composite(X,Y,uu,vv,[0,200,500],counts,backend,expected_sign=1)
        np.testing.assert_equal(result.fit_ok.to_numpy(),[True,False,True])
        self.assertTrue(np.isnan(uf[...,1]).all())
        self.assertEqual(backend.calls,[(0.,0.),(10.,-10.)])
        self.assertEqual(result.xc_from_fit_reference.iloc[2],0.)
        self.assertTrue(result.reason.iloc[1])

    def test_wrong_polarity_and_outer_failure(self):
        x=np.arange(-200,201,10.);X,Y=np.meshgrid(x,x)
        backend=analytic_backend(0.,0.,np.eye(2),1e-5,60.)
        u,v=model(X*1000,Y*1000,0,0,np.eye(2),1e-5,60000)
        args=(X,Y,u[...,None],v[...,None],[0],np.ones((*X.shape,1),int),backend)
        result,_,_=cbt.fit_composite(*args,expected_sign=-1)
        self.assertFalse(result.fit_ok.iloc[0]);self.assertIn('polarity',result.reason.iloc[0])
        backend.out_core_param_fit=lambda *args: (np.nan,np.nan,np.nan)
        result,uf,_=cbt.fit_composite(*args)
        self.assertTrue(result.inner_ok.iloc[0]);self.assertFalse(result.fit_ok.iloc[0])
        self.assertTrue(np.isnan(uf).all())

    def test_plot_spread_and_mean_uncertainty(self):
        members=pd.DataFrame(dict(Eddy=[1,2,1,2],Day=[1,1,1,1],Depth=[0,0,200,200],xc=[0,0,3,5],yc=[0,0,1,2]))
        stats=cbt.centre_statistics(members,[0,200],n_boot=100)
        fits=pd.DataFrame(dict(Depth=[0,200],xc=[1,4],yc=[2,3],fit_ok=[True,False]))
        fig=cbt.plot_centrelines(stats,fits)
        self.assertEqual(len(fig.axes),4)
        fig.canvas.draw();plt.close(fig)

if __name__=='__main__':unittest.main()
