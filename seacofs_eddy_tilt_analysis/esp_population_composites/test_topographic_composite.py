"""Check rotation geometry, signed frames, composite support, and notebook cells."""
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'seacofs_eddy_dataset_modular'/'src'))
import composite_comparison_tools as ccomp
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import topographic_composite_tools as tct
import planetary_composite_tools as pct


def model(x,y,xc,yc,Q,omega,rc):
    dx,dy=x-xc,y-yc
    qx=Q[0,0]*dx+Q[0,1]*dy
    qy=Q[1,0]*dx+Q[1,1]*dy
    factor=omega*np.exp(-(dx*qx+dy*qy)/(2*rc**2))
    return -qy*factor,qx*factor


def unavailable_inner(*args):
    raise ValueError('synthetic notebook smoke test: no inner-fit solution')


def fixture():
    rows=[];profiles=[]
    grid_angle=.23
    Q=np.array([[1.4,.25],[.25,.8]])
    for e in range(1,9):
        phi=(e-1)*np.pi/4
        B=np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
        q=B@Q@B.T
        theta=(90-np.degrees(phi+grid_angle))%360
        for day in range(1,3):
            rows.append(dict(Eddy=e,Day=day,Cyc='AE' if e<=4 else 'CE',Ro=.2,
                lon=154.5 if day==1 else 155.,PV_grad_topo_mag=1e-10,PV_grad_topo_theta=theta,topo_plan_ratio=2.,h=2000.,
                dhdx=np.sin(np.deg2rad(theta)),dhdy=np.cos(np.deg2rad(theta))))
            for z in ([0,200,500,1000] if e%2 else [0,200,500,1000,1500]):
                centre=B@np.array([z/100,z/250])+[100,200]
                profiles.append(dict(Eddy=e,Day=day,Depth=z,xc=centre[0],yc=centre[1],
                    q11=q[0,0],q12=q[0,1],q22=q[1,1],Rc=30.,Omega=.0001 if e<=4 else -.0001))
    return pd.DataFrame(rows),pd.DataFrame(profiles),grid_angle,Q


class TopographicTests(unittest.TestCase):
    def test_signed_axes_and_invalid_frames(self):
        s=pd.DataFrame(dict(Eddy=[1,2,3,4],Day=[1]*4,PV_grad_topo_mag=[1,1,0,1],
            PV_grad_topo_theta=[0,180,90,np.nan],dhdx=[0]*4,dhdy=[1]*4))
        frames=tct.reference_frames(s,0)
        np.testing.assert_array_equal(frames.frame_valid,[True,True,False,False])
        np.testing.assert_allclose(frames.gradient_dot_downslope_cos[:2],[1,-1])
        # Northward gradient: north -> +along; west -> +perp. Opposite gradient reverses both.
        x,y=tct.rotate_xy(np.array([0.,-2.]),np.array([3.,0.]),frames.frame_angle_rad.iloc[0])
        np.testing.assert_allclose(x,[3,0],atol=1e-12)
        np.testing.assert_allclose(y,[0,2],atol=1e-12)
        x,y=tct.rotate_xy(0.,3.,frames.frame_angle_rad.iloc[1])
        self.assertAlmostEqual(x,-3.)

    def test_inverse_grid_velocity_and_anisotropic_ellipse(self):
        s,v,angle,Q=fixture()
        # Different geographical centres/ellipses become the same aligned field.
        s=s.loc[s.Eddy.le(4)]
        frames=tct.reference_frames(s,angle)
        X,Y=np.meshgrid(np.linspace(-60,60,13),np.linspace(-50,50,11))
        u,w,z,counts,support,audit,m=tct.composite_days(frames,v,X,Y,esp=SimpleNamespace(model_uv_at_xy=model))
        for k,depth in enumerate(z):
            expected=model(X*1000,Y*1000,depth/100*1000,depth/250*1000,Q,.0001,30000)
            np.testing.assert_allclose(u[...,k],expected[0],atol=1e-14)
            np.testing.assert_allclose(w[...,k],expected[1],atol=1e-14)
        np.testing.assert_allclose(m.along_km,m.Depth/100,atol=1e-12)
        np.testing.assert_allclose(m.perp_km,m.Depth/250,atol=1e-12)
        np.testing.assert_array_equal(support.eddy_days,[8,8,8,8,4])
        np.testing.assert_array_equal(counts[0,0],support.eddy_days)
        self.assertTrue(audit.included.all())

    def test_centres_only_matches_finite_reconstruction(self):
        s,v,angle,_=fixture();frames=tct.reference_frames(s,angle)
        audit,centres=pct.prepare_profiles(frames,v)
        m=tct.rotate_members(centres,frames)
        X,Y=np.meshgrid([-10.,0.,10.],[-10.,0.,10.])
        *_,reconstructed=tct.composite_days(frames,v,X,Y,esp=SimpleNamespace(model_uv_at_xy=model))
        cols=['Eddy','Day','Depth','along_km','perp_km','distance_km']
        pd.testing.assert_frame_equal(m[cols].sort_values(cols[:3]).reset_index(drop=True),
            reconstructed[cols].sort_values(cols[:3]).reset_index(drop=True),check_dtype=False)
        stats,boot=tct.summarise(m,100)
        np.testing.assert_allclose(stats.mean_along,stats.Depth/100,atol=1e-12)
        self.assertIn('cov_along_perp',stats)
        self.assertNotIn('bearing_deg',stats)
        np.testing.assert_allclose(stats.relative_angle_deg.iloc[1:],np.degrees(np.arctan2(1/250,1/100)))
        self.assertTrue(np.isnan(stats.relative_angle_deg.iloc[0]))

    def test_finite_support_and_missing_reference(self):
        s,v,angle,_=fixture();frames=tct.reference_frames(s,angle)
        v=v.loc[~(v.Eddy.eq(1)&v.Depth.eq(0))]
        def missing(*args):
            u,w=model(*args)
            # All fields at this first grid location excluded, not counted as zeros.
            u[0,0]=np.nan
            return u,w
        X,Y=np.meshgrid([-10.,0.,10.],[-10.,0.,10.])
        u,w,z,count,support,audit,m=tct.composite_days(frames,v,X,Y,esp=SimpleNamespace(model_uv_at_xy=missing))
        self.assertTrue(np.isnan(u[0,0]).all())
        self.assertTrue((count[0,0]==0).all())
        self.assertFalse(audit.loc[audit.Eddy.eq(1)].included.any())
        np.testing.assert_array_equal(support.eddy_days,m.groupby('Depth').size().to_numpy())

    def test_notebook_end_to_end(self):
        s,v,angle,_=fixture()
        grid=SimpleNamespace(angle=angle,z_r=np.full((151,151,6),-500.))
        tilt=SimpleNamespace(Paths=lambda:SimpleNamespace(grid='test',z_r='test'),
            load_grid=lambda *a:grid,load_tilt_tables=lambda *a,**kw:(s.copy(),None),
            add_pv_gradient_terms=lambda data,*a,**kw:data,load_vert=lambda *a:v)
        ns=dict(np=np,pd=pd,plt=plt,pct=pct,tct=tct,tilt=tilt,ccomp=ccomp,
            esp=SimpleNamespace(model_uv_at_xy=model,doppio=unavailable_inner,out_core_param_fit=lambda *a:None),display=lambda *a:None)
        cells=json.loads(Path(__file__).with_name('topographic_composite_tilt.ipynb').read_text())['cells']
        for c in cells[2:]:
            if c['cell_type']=='code':
                exec(''.join(c['source']),ns)
                plt.close('all')
        self.assertEqual(len(ns['results']),6)
        self.assertTrue(ns['centre_stats'].Depth.gt(1000).any())
        self.assertEqual(len(ns['frame_audit']),len(s))
        # Exercise the explicitly supported centre-only mode as well.
        ns['RUN_VELOCITY_COMPOSITES']=False
        exec(''.join(cells[4]['source']),ns)
        self.assertNotIn('u',ns['results']['all','AE'])

if __name__=='__main__':unittest.main()
