"""Scientific invariants and a complete synthetic run of the regional notebook."""
import sys,json,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'seacofs_eddy_dataset_modular'/'src'))
import regional_composite_tools as rct
import planetary_composite_tools as pct
import composite_comparison_tools as ccomp


def model(x,y,xc,yc,Q,omega,rc):
    f=omega*np.exp(-((x-xc)**2+(y-yc)**2)/rc**2)
    return -(y-yc)*f,(x-xc)*f


def no_inner(*a):raise ValueError('deliberate synthetic refit failure')


def fixture():
    rows=[];profiles=[];e=0
    for region in rct.REGIONS:
        for cyc in ['AE','CE']:
            for cohort in ['shallow','deep']:
                for ro in [.3,.5]:
                    for member in range(2):
                        e+=1
                        for day in [1,2]:
                            rows.append(dict(Eddy=e,Day=day,Region=region,Cyc=cyc,Rc=40.+member*5,Ro=-ro if cyc=='AE' else ro,
                                topo_plan_ratio=np.log(3) if member else -np.log(3),lon=154. if day==1 else 155.,h=1800.))
                            for z in ([0.,200.,500.,850.] if cohort=='shallow' else [0.,200.,500.,850.,1500.]):
                                profiles.append(dict(Eddy=e,Day=day,Depth=z,xc=100+z/100+member,yc=100+z/200,
                                    Rc=25.+member*5,Omega=1e-5 if cyc=='AE' else -1e-5,q11=1.,q12=0.,q22=1.))
    return pd.DataFrame(rows),pd.DataFrame(profiles)


def grid_fixture():
    x,y=np.meshgrid(np.arange(0.,301.,10),np.arange(0.,301.,10),indexing='ij')
    return SimpleNamespace(angle=.2,X_grid=x,Y_grid=y,h=np.full_like(x,1000.),mask_rho=np.ones_like(x,bool))


class RegionalTests(unittest.TestCase):
    def test_signed_rossby_boundary_and_inventory(self):
        s,v=fixture();selection,audit,p=rct.prepare_population(s,v)
        inv=rct.inventory(audit)
        self.assertEqual(len(inv),60)
        self.assertTrue(inv.eddies.gt(0).all())
        self.assertTrue(audit.loc[audit.Ro.abs().eq(.5)].Ro_class.eq('high').all())
        self.assertTrue(audit.loc[audit.Ro.abs().eq(.3)].Ro_class.eq('low').all())
        s.loc[0,'Ro']=np.nan
        _,a,p=rct.prepare_population(s,v)
        self.assertEqual(a.Ro_class.eq('unknown').sum(),1)
        row=a.loc[a.Ro_class.eq('unknown')].iloc[0]
        self.assertEqual(len(rct.memberships(row)),1)
        s.loc[0,'Region']='outside'
        selection,a,p=rct.prepare_population(s,v)
        self.assertEqual(selection.selection_status.eq('outside named regions').sum(),1)

    def test_one_reconstruction_and_match_previous_composite(self):
        s,v=fixture();s=s.loc[s.Region.eq('U1')&s.Cyc.eq('AE')]
        _,audit,p=rct.prepare_population(s,v)
        X,Y=np.meshgrid(np.arange(-40,41,10.),np.arange(-40,41,10.))
        calls=[]
        def counted(*args):calls.append(1);return model(*args)
        results,days=rct.build_composites(audit,p,X,Y,SimpleNamespace(model_uv_at_xy=counted),.2,progress_every=0)
        self.assertEqual(len(calls),len(p))
        self.assertEqual(len(results),5)
        self.assertEqual(len(days),len(s))
        u,w,z,count,support,_,m=pct.composite_days(s,v,X,Y,return_centres=True,esp=SimpleNamespace(model_uv_at_xy=model))
        base=results['U1','all','all','AE']
        np.testing.assert_allclose(base['u'],u);np.testing.assert_allclose(base['v'],w)
        np.testing.assert_array_equal(base['counts'],count)
        for r in results.values():
            np.testing.assert_array_equal(r['support'].eddy_days,r['members'].groupby('Depth').size().to_numpy())

    def test_normalise_before_averaging(self):
        m=pd.DataFrame(dict(Eddy=[1,2],Day=[1,1],Depth=[500.,500.],east_km=[10.,30.],north_km=[0.,0.],
                            distance_km=[10.,30.],surface_Rc_km=[10.,30.]))
        s=rct.normalised_statistics(m,100)
        self.assertAlmostEqual(s.distance_Rc.iloc[0],1.)
        m.loc[1,'surface_Rc_km']=10.
        s=rct.normalised_statistics(m,100)
        self.assertAlmostEqual(s.distance_Rc.iloc[0],2.)
        m.loc[1,'surface_Rc_km']=60.
        s=rct.normalised_statistics(m,100)
        self.assertAlmostEqual(s.distance_Rc.iloc[0],.75)
        self.assertNotAlmostEqual(s.distance_Rc.iloc[0],20/35)
        self.assertIn('var_east_Rc2',s)
        self.assertNotIn('distance_km',s)

    def test_fixed_membership_exact_depths_and_missing(self):
        m=pd.DataFrame(dict(Eddy=[1,1,1,2,2],Day=[1]*5,Depth=[0.,200.,500.,0.,200.],
            xc=[0,1,2,0,3],yc=[0]*5,surface_Rc_km=[20.]*5))
        m=pct.geographic_members(m,0.)
        key=('U1','shallow','low','AE');results={key:dict(members=m)}
        mapping=rct.resolve_depths([0,200,500],[200,500],10)
        fixed,audit=rct.fixed_membership(results,mapping,[200,500],[200,500],100)
        self.assertEqual(audit.matched_eddies.iloc[0],1)
        np.testing.assert_array_equal(fixed[key]['stats'].n_eddy_days,[1,1,1])
        self.assertTrue(fixed[key]['stats'].mean_east_ci_low.isna().all())
        mapping=rct.resolve_depths([0,200,500],[200,700],10)
        fixed,audit=rct.fixed_membership(results,mapping,[200,700],[200,700],100)
        self.assertFalse(fixed);self.assertIn('unresolved',audit.status.iloc[0])

    def test_bottom_audit_does_not_filter(self):
        s,v=fixture();_,a,p=rct.prepare_population(s,v)
        checked=rct.bottom_audit(p,grid_fixture())
        self.assertEqual(len(checked),len(p))
        self.assertTrue(checked.loc[checked.Depth.gt(1000)].below_local_bottom.eq(1).all())
        self.assertTrue(checked.loc[checked.Depth.le(1000)].below_local_bottom.eq(0).all())

    def test_nonfinite_velocity_members_and_empty_slots(self):
        s,v=fixture();s=s.iloc[:2].copy();_,audit,p=rct.prepare_population(s,v)
        X,Y=np.meshgrid([-20.,0.,20.],[-20.,0.,20.])
        def masked(*args):
            u,w=model(*args)
            if args[2]>5000:u[:]=np.nan
            return u,w
        results,days=rct.build_composites(audit,p,X,Y,SimpleNamespace(model_uv_at_xy=masked),0.,progress_every=0)
        self.assertTrue(all(r['members'].Depth.max()==500 for r in results.values()))
        inv=rct.inventory(audit)
        self.assertEqual(inv.eddy_days.eq(0).sum(),58)
        rct.add_statistics(results,100)
        table=rct.summary_at_depth(results,850)
        self.assertTrue(table.status.eq('depth unavailable').all())

    def test_notebook_end_to_end(self):
        s,v=fixture();grid=grid_fixture()
        tilt=SimpleNamespace(Paths=lambda:SimpleNamespace(grid='test',z_r='test'),
            load_grid=lambda *a:grid,load_tilt_tables=lambda *a,**kw:(s.copy(),None),
            add_pv_gradient_terms=lambda data,*a,**kw:data,add_region_labels=lambda data,*a:data,
            load_vert=lambda *a:v)
        ns=dict(np=np,pd=pd,plt=plt,Path=Path,sys=sys,json=json,rct=rct,pct=pct,ccomp=ccomp,tilt=tilt,
            esp=SimpleNamespace(model_uv_at_xy=model,doppio=no_inner,out_core_param_fit=lambda *a:None),
            display=lambda *a:None)
        cells=json.loads(Path(__file__).with_name('regional_composite_tilt.ipynb').read_text())['cells']
        for c in cells[2:]:
            if c['cell_type']=='code':
                source=''.join(c['source']).replace('BOOTSTRAPS=500','BOOTSTRAPS=100').replace(
                    'LOAD_STRATIFICATION_DIAGNOSTICS=True','LOAD_STRATIFICATION_DIAGNOSTICS=False')
                exec(source,ns);plt.close('all')
        self.assertEqual(len(ns['regional_results']),42)
        self.assertEqual(len(ns['matched_results']),36)
        self.assertEqual(len(ns['composite_fit_audit']),24)
        self.assertEqual(ns['population_inventory'].shape[0],42)
        self.assertTrue(ns['centre_stats'].below_bottom_fraction.max()>0)
        self.assertFalse(ns['population_context'].empty)

if __name__=='__main__':unittest.main()
