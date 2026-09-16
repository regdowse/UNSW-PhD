import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import population_tools as p


def fixture():
    surface=pd.DataFrame([dict(Eddy=e,Day=d,Cyc='AE',xc=0.,yc=0.,Rc=1.,
                             q11=1.,q12=0.,q22=1.,axis_x=1.,axis_y=0.)
                          for e,ds in [(1,[1,2,3]),(2,[1])] for d in ds])
    vertical=pd.DataFrame([dict(Eddy=r.Eddy,Day=r.Day,Depth=z,xc=(0 if z==0 else r.Eddy),
                               yc=0.,Rc=1.,Omega=1.,q11=1.,q12=0.,q22=1.)
                          for r in surface.itertuples() for z in [0.,200.]])
    pv=surface[p.KEYS].copy();pv['topo_plan_ratio']=-2.
    return surface,vertical,pv


def velocity(x,y,xc,yc,q,omega,rc):
    # Analytic circular vortex backend for tests only, same signature as ESP.
    g=np.exp(-((x-xc)**2+(y-yc)**2)/rc**2)
    return -(y-yc)/rc*g*omega,(x-xc)/rc*g*omega


class Tests(unittest.TestCase):
    def test_matched_depths(self):
        s,v,pv=fixture();v=v.drop(v[(v.Eddy==2)&(v.Depth==200)].index)
        a,_=p.audit_population(s,v,pv,[0,200])
        self.assertFalse(a.loc[a.Eddy==2,'eligible'].iloc[0])
        self.assertTrue(a.loc[a.Eddy==1,'eligible'].all())

    def test_invalid_geometry_and_polarity(self):
        s,v,pv=fixture();v.loc[0,'q12']=2
        v.loc[v.Eddy==2,'Omega']=-1
        a,_=p.audit_population(s,v,pv,[0,200])
        self.assertFalse(a.iloc[0].eligible)
        self.assertFalse(a.iloc[-1].eligible)

    def test_duplicate_keys_fail(self):
        s,v,pv=fixture()
        with self.assertRaises(ValueError):p.audit_population(s,pd.concat([v,v.iloc[:1]]),pv,[0,200])

    def test_depths_never_interpolated(self):
        _,v,_=fixture()
        np.testing.assert_equal(p.select_depths(v,[0,190]),[0,200])
        with self.assertRaises(ValueError):p.select_depths(v,[0,500])

    def test_missing_pv_and_small_tilt(self):
        s,v,pv=fixture();s['TiltDis']=0;pv.loc[0,'topo_plan_ratio']=np.nan
        a,_=p.audit_population(s,v,pv,[0,200])
        self.assertFalse(a.iloc[0].eligible)
        self.assertTrue(a.iloc[1].eligible)

    def test_rotation_preserves_tilt_and_vectors(self):
        s,v,_=fixture();row=s.iloc[0]
        angle=.4;B,G=p.rotation_basis(row,angle,'geographic')
        np.testing.assert_allclose(G,np.eye(2),atol=1e-12)
        field,c=p.reconstruct_day(row,v[(v.Eddy==1)&(v.Day==1)],[0,200],np.array([-1.,0,1]),
                                   angle,velocity,frame='geographic')
        np.testing.assert_allclose(c[0],[0,0])
        np.testing.assert_allclose(c[1],[np.cos(angle),-np.sin(angle)])
        self.assertGreater(np.linalg.norm(c[1]),0)
        row=row.copy();row.axis_x=0;row.axis_y=1
        _,G=p.rotation_basis(row,0,'onshore')
        np.testing.assert_equal(G,[[0,1],[-1,0]])

    def test_ocean_mask(self):
        g=SimpleNamespace(x_grid=np.array([-1.,1.]),y_grid=np.array([-1.,1.]),
                          h=np.array([[100.,300.],[300.,300.]]),mask_rho=np.array([[1,1],[0,1]]))
        mask=p.ocean_mask(g)
        np.testing.assert_equal(mask(np.array([-1,-1,1,5]),np.array([-1,1,-1,0]),200),[False,True,False,False])

    def test_equal_eddy_not_days_and_colliding_days(self):
        s,v,pv=fixture();a,pr=p.audit_population(s,v,pv,[0,200]);a['group']='AE_Planetary'
        grid=SimpleNamespace(angle=0)
        with tempfile.TemporaryDirectory() as t:
            manifest=p.build_members(a,pr,[0,200],np.array([-1.,0,1]),grid,velocity,t,mask_ocean=False)
            self.assertEqual(len(manifest),2)
            result=p.summarise_members(manifest,t,n_boot=200,min_members=2)['AE_Planetary']
            self.assertAlmostEqual(result['centre']['mean'][1,0],1.5)
            self.assertNotAlmostEqual(result['centre']['mean'][1,0],1.25)
            self.assertEqual(result['days'],4)
            self.assertEqual(result['centre']['support'][1,0],2)
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            for name,fig in p.plot_summary({'AE_Planetary':result},[0,200],np.array([-1.,0,1]))+p.plot_sections_3d(result,[0,200],np.array([-1.,0,1])):
                fig.savefig(Path(t)/(name+'.png'));plt.close(fig)

    def test_intervals_and_missing_support(self):
        values=np.array([[1.,np.nan],[3.,8.]])
        draws=np.array([[2,0],[0,2],[1,1]])
        r=p.bootstrap_mean(values,draws,min_members=2)
        self.assertAlmostEqual(r['mean'][0],2)
        self.assertLess(r['low'][0],2);self.assertGreater(r['high'][0],2)
        self.assertTrue(np.isnan(r['mean'][1]))
        np.testing.assert_equal(r['support'],[2,1])

    def test_shared_ids_have_identical_intervals(self):
        with tempfile.TemporaryDirectory() as t:
            rows=[]
            for group in ['a','b']:
                for e in range(3):
                    f=f'{group}{e}.npz'
                    np.savez(Path(t)/f,mean=np.array([e*1.]),centres=np.array([[e*1.,0]]))
                    rows.append(dict(group=group,Eddy=e,file=f,days=1,frame='geographic'))
            r=p.summarise_members(pd.DataFrame(rows),t,n_boot=100,min_members=2)
            np.testing.assert_equal(r['a']['centre']['low'],r['b']['centre']['low'])
            np.testing.assert_equal(r['a']['centre']['high'],r['b']['centre']['high'])



class ReportingTests(unittest.TestCase):
    def result(self):
        centre={'mean':np.array([[0.,0.],[-.05,.02]]),
                'low':np.array([[0.,0.],[-.08,-.01]]),
                'high':np.array([[0.,0.],[-.02,.05]]),
                'support':np.full((2,2),30)}
        return {'centre':centre,'frame':'onshore','days':100,'eddies':30,
                'members':np.array([centre['mean']])}

    def test_interval_sign_and_zero(self):
        table=p.centre_interval_table({'AE_Topographic':self.result()},[0,500])
        r=table.loc[(table.depth_m==500)&(table.component=='onshore')].iloc[0]
        self.assertAlmostEqual(r.tilt_mean,.05)
        self.assertAlmostEqual(r.tilt_ci_low,.02)
        self.assertAlmostEqual(r.tilt_ci_high,.08)
        self.assertTrue(r.excludes_zero)
        self.assertFalse(table.loc[(table.depth_m==500)&(table.component=='alongshore'),'excludes_zero'].iloc[0])

    def test_all_levels_preset(self):
        v=pd.DataFrame({'Depth':[0,5,20,50,181.264,300,515.416,660]})
        np.testing.assert_equal(p.depth_preset(v,'full_500m'),[0,5,20,50,181.264,300,515.416])
        np.testing.assert_equal(p.depth_preset(v,'shallow_200m'),[0,5,20,50,181.264])

    def test_zoom_ignores_extreme_members(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        r=self.result();r['members']=np.array([[[0,0],[1000,-1000]]])
        fig=p.plot_zoomed_centres({'AE_Topographic':r},[0,500])
        self.assertLess(fig.axes[0].get_xlim()[1],.1)
        np.testing.assert_equal(fig.axes[0].get_ylim(),[500,0])
        plt.close(fig)

    def test_contrast_sign_and_support(self):
        with tempfile.TemporaryDirectory() as t:
            rows=[]
            for cyc,value in [('AE',1.),('CE',3.)]:
                for e in range(3):
                    file=f'{cyc}{e}.npz'
                    np.savez(Path(t)/file,centres=np.array([[0.,0.],[value,0.]]))
                    rows.append(dict(group=f'{cyc}_Planetary',Eddy=e+(0 if cyc=='AE' else 3),
                                     file=file,days=1,frame='geographic'))
            r=p.centre_contrasts(pd.DataFrame(rows),t,[0,500],n_boot=100,min_members=2)
            row=r.loc[(r.depth_m==500)&(r.component=='east')].iloc[0]
            self.assertAlmostEqual(row['mean'],2.)
            self.assertAlmostEqual(row.ci_low,2.)
            self.assertAlmostEqual(row.ci_high,2.)
            r=p.centre_contrasts(pd.DataFrame(rows),t,[0,500],n_boot=100,min_members=4)
            self.assertTrue(r['mean'].isna().all())

    def test_section_limits(self):
        import matplotlib.pyplot as plt
        r=self.result()
        r['field']={'mean':np.ones((2,2,3,3))}
        figs=p.plot_sections_3d(r,np.array([0.,515.4]),np.array([-1,0,1]))
        for ax in figs[0][1].axes[:2]:
            np.testing.assert_allclose(ax.get_ylim(),[515.4,0])
        for _,fig in figs:plt.close(fig)

if __name__=='__main__': unittest.main()
