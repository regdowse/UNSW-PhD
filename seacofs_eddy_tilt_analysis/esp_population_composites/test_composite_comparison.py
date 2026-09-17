"""Boundary, shared-depth, full-refit and geographic-section regressions."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'seacofs_eddy_dataset_modular'/'src'))
import composite_comparison_tools as cc
from test_breakdown_centrelines import analytic_backend,model
from test_topographic_composite import fixture
import topographic_composite_tools as tct
import planetary_composite_tools as pct
from types import SimpleNamespace


def result(depths=(0.,480.,1400.),n=(3,3,3)):
    depths=np.array(depths)
    return dict(stats=pd.DataFrame(dict(Depth=depths,n_eddies=n,n_eddy_days=n,
                distance_km=depths/100,distance_ci_low=depths/100-.1,distance_ci_high=depths/100+.1,
                mean_member_distance_km=depths/80)),
                members=pd.DataFrame(dict(Eddy=[1,2,3],Day=[1,1,1])))


class ComparisonTests(unittest.TestCase):
    def test_shelf_boundary_and_missing(self):
        df=pd.DataFrame({'lon':[154.74,154.75,154.76,np.nan]})
        np.testing.assert_equal(cc.shelf_classes(df).shelf_class.to_numpy(),['On-shelf','Off-shelf','Off-shelf','Unknown'])
        self.assertNotIn('shelf_class',df)

    def test_shared_depth_and_distinct_endpoints(self):
        collections={'P':{('shallow','AE'):result((0,480),(3,3)),('deep','CE'):result()}}
        table=cc.tilt_summary(collections)
        shared=table.loc[table.summary.eq('Shared reference depth')]
        np.testing.assert_allclose(shared.Depth,480)
        end=table.loc[table.summary.eq('Deepest supported depth')]
        np.testing.assert_allclose(end.Depth,[480,1400])
        self.assertTrue(table.status.eq('ok').all())
        # No individual nearest-depth substitution when populations have different levels.
        collections['P']['deep','CE']=result((0,490,1400))
        table=cc.tilt_summary(collections)
        self.assertTrue(table.loc[table.summary.eq('Shared reference depth')].distance_km.isna().all())
        # Even a shared level is unavailable if outside the explicit tolerance.
        collections['P']['shallow','AE']=result((0,490,1400))
        table=cc.tilt_summary(collections,target_depth=800,depth_tolerance=100)
        self.assertTrue(table.loc[table.summary.eq('Shared reference depth')].Depth.isna().all())

    def test_rotated_sections_and_support(self):
        x=np.linspace(-20,20,21);y=np.linspace(-30,30,31)
        X,Y=np.meshgrid(x,y);angle=.3;c,s=np.cos(angle),np.sin(angle)
        E,N=c*X-s*Y,s*X+c*Y
        ue,vn=3*E-2*N,2*E+3*N
        ug,vg=c*ue+s*vn,-s*ue+c*vn
        r=result((0,480,1400),(3,3,1))
        r.update(depths=np.array([0,480,1400]),u=np.repeat(ug[...,None],3,axis=2),
                 v=np.repeat(vg[...,None],3,axis=2),counts=np.ones((*X.shape,3),int))
        cuts=cc.section_data(r,X,Y,angle,min_eddies=2)
        np.testing.assert_allclose(cuts['xcut'][:2],np.tile(2*x,(2,1)),atol=1e-12)
        np.testing.assert_allclose(cuts['ycut'][:2],np.tile(-2*y,(2,1)),atol=1e-12)
        self.assertTrue(np.isnan(cuts['xcut'][2]).all())
        self.assertTrue(np.isnan(cuts['ycut'][2]).all())

    def test_full_composite_fit_and_failure_gap(self):
        x=np.arange(-200,201,10.);X,Y=np.meshgrid(x,x);Q=np.eye(2)
        backend=analytic_backend(0,0,Q,1e-5,60)
        u,v=model(X*1000,Y*1000,0,0,Q,1e-5,60000)
        r=result((0,480,1400),(3,1,3))
        r.update(depths=np.array([0,480,1400]),u=np.repeat(u[...,None],3,axis=2),
                 v=np.repeat(v[...,None],3,axis=2),counts=np.ones((*X.shape,3),int))
        fits,audit=cc.fit_collections({'P':{('deep','AE'):r}},X,Y,backend)
        np.testing.assert_equal(fits.fit_ok.to_numpy(),[True,False,True])
        np.testing.assert_allclose(fits.loc[fits.fit_ok,'Rc'],60,atol=1e-4)
        np.testing.assert_allclose(fits.loc[fits.fit_ok,'Omega'],1e-5,atol=1e-10)
        self.assertIn('minimum',fits.reason.iloc[1])
        self.assertEqual(audit.successful_depths.iloc[0],2)
        self.assertEqual(len(r['composite_fit']),3)
        cc.plot_fit_profiles(fits,'P').canvas.draw();plt.close('all')
        cc.plot_sections({('deep','AE'):r},X,Y,'P','deep').canvas.draw();plt.close('all')
        skipped,audit=cc.fit_collections({'P':{('deep','AE'):result()}},X,Y,backend)
        self.assertTrue(skipped.empty)
        self.assertIn('disabled',audit.status.iloc[0])

    def test_shelf_centres_and_field_members_agree(self):
        s,v,angle,_=fixture()
        s['lon']=np.where(s.Day.eq(1),154.,155.)
        s=cc.shelf_classes(tct.reference_frames(s,angle))
        audit,centres=pct.prepare_profiles(s,v)
        X,Y=np.meshgrid([-20.,0,20],[-20.,0,20])
        from test_topographic_composite import model as backend_model
        esp=SimpleNamespace(model_uv_at_xy=backend_model)
        fast=cc.shelf_composites(audit,centres,v,X,Y,esp,velocity=False,n_boot=100)
        full=cc.shelf_composites(audit,centres,v,X,Y,esp,velocity=True,n_boot=100)
        for shelf in fast:
            self.assertEqual(len(fast[shelf]),6)
            for key in fast[shelf]:
                pd.testing.assert_frame_equal(fast[shelf][key]['stats'],full[shelf][key]['stats'],check_exact=False,atol=1e-12)
        on=full['On-shelf']['all','AE']['members'];off=full['Off-shelf']['all','AE']['members']
        self.assertTrue(on.Day.eq(1).all());self.assertTrue(off.Day.eq(2).all())

if __name__=='__main__':unittest.main()
