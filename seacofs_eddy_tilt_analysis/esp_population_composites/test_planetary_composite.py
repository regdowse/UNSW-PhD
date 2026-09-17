"""Synthetic scientific invariants and executable notebook workflow."""
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import planetary_composite_tools as pct


def fixture():
    sample=[];vertical=[]
    for e in range(1,9):
        for day in range(1,1+e%3+1):
            sample.append(dict(Eddy=e,Day=day,Cyc='AE' if e<5 else 'CE',Ro=.06*e))
            zs=[0,200,500,1000] if e%2 else [0,200,500,1000,1500,2000]
            for z in zs:
                vertical.append(dict(Eddy=e,Day=day,Depth=z,xc=100+z/100*(1 if e<5 else -1),
                    yc=200+z/200+day*z/1000,Rc=30.,Omega=.0001,q11=1.,q12=0.,q22=1.))
    return pd.DataFrame(sample),pd.DataFrame(vertical)


def model(x,y,xc,yc,q,omega,rc):
    factor=omega*np.exp(-((x-xc)**2+(y-yc)**2)/(2*rc**2))
    return -(y-yc)*factor,(x-xc)*factor


class PlanetaryTests(unittest.TestCase):
    def test_extent_reference_and_missing_surface(self):
        sample,v=fixture()
        v=v.loc[~(v.Eddy.eq(1)&v.Depth.eq(0))]
        audit,centred=pct.prepare_profiles(sample,v)
        self.assertTrue(audit.loc[audit.Eddy.eq(1)].profile_status.eq('no valid surface fit').all())
        self.assertTrue(audit.loc[audit.Eddy.eq(3)].extent_group.eq('shallow').all())
        self.assertTrue(audit.loc[audit.Eddy.eq(2)].extent_group.eq('deep').all())
        np.testing.assert_allclose(centred.loc[centred.Depth.eq(0),['xc','yc']],0)

    def test_rotation_bearing_and_cancellation(self):
        m=pd.DataFrame(dict(Eddy=[1,2],Day=[1,1],Depth=[500,500],xc=[10.,-10.],yc=[0.,0.]))
        g=pct.geographic_members(m,np.pi/2)
        np.testing.assert_allclose(g.north_km,[10,-10])
        np.testing.assert_allclose(g.bearing_deg,[0,180],atol=1e-12)
        s,b=pct.summarise(g,100)
        self.assertAlmostEqual(s.distance_km.iloc[0],0)
        self.assertTrue(np.isnan(s.bearing_deg.iloc[0]))
        self.assertAlmostEqual(s.mean_member_distance_km.iloc[0],10)
        self.assertAlmostEqual(s.coherence.iloc[0],0)
        self.assertEqual(len(b),100)

    def test_exact_matched_membership(self):
        sample,v=fixture();_,m=pct.prepare_profiles(sample,v)
        matched=pct.matched_members(m,[0,500,1500])
        self.assertTrue((matched.Eddy%2==0).all())
        self.assertTrue(pct.matched_members(m,[0,501]).empty)

    def test_reconstruction_records_and_day_weighting(self):
        sample,v=fixture()
        x,y=np.meshgrid(np.arange(-40,41,10),np.arange(-40,41,10))
        def partial(*args):
            u,w=model(*args)
            if args[2]>15000:u[:]=np.nan
            return u,w
        u,w,z,counts,support,audit,m=pct.composite_days(sample,v,x,y,return_centres=True,esp=SimpleNamespace(model_uv_at_xy=partial))
        self.assertFalse(m.xc.gt(15).any())
        actual=m.groupby('Depth').size().reindex(z,fill_value=0)
        np.testing.assert_equal(actual.to_numpy(),support.eddy_days.to_numpy())
        g=pct.geographic_members(m,.2);stats,_=pct.summarise(g,100)
        expected=g.groupby('Depth').east_km.mean()
        np.testing.assert_allclose(stats.mean_east,expected)

    def test_notebook_end_to_end_synthetic(self):
        sample,v=fixture()
        ns=dict(np=np,pd=pd,plt=plt,pct=pct,esp=SimpleNamespace(model_uv_at_xy=model),
            grid=SimpleNamespace(angle=.2),planetary=sample,vertical=v,display=lambda *a:None,
            SPLIT_DEPTH_M=1000.,TARGET_DEPTHS_M=[200.,500.,1000.,1500.,2000.],
            BOOTSTRAPS=100,SEED=731,MIN_DIRECTION_DISTANCE_KM=0.,MIN_PLOT_EDDIES=2,
            WIDTH_KM=100.,RES_KM=10.,RUN_VELOCITY_COMPOSITES=True)
        ns['profile_audit'],ns['centred_profiles']=pct.prepare_profiles(sample,v)
        cells=json.loads(Path(__file__).with_name('planetary_composite_tilt.ipynb').read_text())['cells']
        for cell in cells[4:]:
            if cell['cell_type']=='code':
                exec(''.join(cell['source']),ns)
                plt.close('all')
        self.assertEqual(len(ns['results']),6)
        self.assertFalse(ns['centre_members'].empty)
        self.assertTrue({'mean_east','mean_north','distance_km'}.issubset(ns['centre_stats']))
        self.assertTrue(ns['centre_stats'].Depth.gt(1000).any())

if __name__=='__main__':unittest.main()
