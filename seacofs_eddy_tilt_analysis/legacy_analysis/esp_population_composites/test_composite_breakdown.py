"""Execute the notebook's own helpers with synthetic ESP output."""
import json
import unittest
from pathlib import Path
from types import SimpleNamespace
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NB=Path(__file__).with_name('composite_breakdown.ipynb')
CELLS=json.loads(NB.read_text())['cells']

def model(x,y,xc,yc,q,omega,rc):
    return np.full_like(x,omega),np.full_like(y,2*omega)

class BreakdownTests(unittest.TestCase):
    def setUp(self):
        self.ns=dict(np=np,pd=pd,plt=plt,esp=SimpleNamespace(model_uv_at_xy=model))
        helper=next(''.join(c['source']) for c in CELLS if ''.join(c['source']).startswith('def valid_profiles'))
        exec(helper,self.ns)
        # First selected row has no profile. Day 10 belongs to two distinct eddies.
        self.sample=pd.DataFrame({'Eddy':[9,1,1,2],'Day':[0,10,11,10]})
        self.vertical=pd.DataFrame([dict(Eddy=e,Day=d,Depth=z,xc=e+z/100,yc=d,
                                    Rc=20.,Omega=w,q11=1.,q12=0.,q22=1.)
                              for e,d,zs,w in [(1,10,[0,200],1.),(1,11,[0],3.),(2,10,[0,500],5.)]
                              for z in zs])
        self.X,self.Y=np.meshgrid([-1.,0,1],[-1.,0,1])

    def test_all_days_variable_depths_and_missing_first(self):
        u,v,depths,counts,support,audit=self.ns['composite_days'](self.sample,self.vertical,self.X,self.Y)
        np.testing.assert_equal(depths,[0,200,500])
        np.testing.assert_allclose(u[0,0],[3,1,5])
        np.testing.assert_allclose(v[0,0],[6,2,10])
        np.testing.assert_equal(counts[0,0],[3,1,1])
        self.assertEqual(audit.included.sum(),3)
        self.assertFalse(audit.iloc[0].included)
        np.testing.assert_equal(support.eddy_days.to_numpy(),[3,1,1])

    def test_selection_and_no_mutation(self):
        subset=self.sample.loc[self.sample.Eddy.eq(2)]
        before=self.vertical.copy(deep=True)
        row,profile=self.ns['choose_case'](subset,self.vertical)
        self.assertEqual(row.Eddy,2)
        shifted=self.ns['recenter'](profile)
        np.testing.assert_allclose(shifted.xc,[0,5])
        pd.testing.assert_frame_equal(self.vertical,before)

    def test_absent_reference_and_empty_sample(self):
        v=self.vertical.loc[~((self.vertical.Eddy==2)&(self.vertical.Depth==0))]
        *_,audit=self.ns['composite_days'](self.sample,v,self.X,self.Y)
        self.assertEqual(audit.included.sum(),2)
        with self.assertRaises(ValueError):
            self.ns['choose_case'](self.sample.iloc[:1],v)
        with self.assertRaises(ValueError):
            self.ns['composite_days'](self.sample.iloc[:0],v,self.X,self.Y)

    def test_duplicate_and_invalid_rows(self):
        with self.assertRaises(ValueError):
            self.ns['composite_days'](pd.concat([self.sample,self.sample.iloc[:1]]),self.vertical,self.X,self.Y)
        v=self.vertical.copy();v.loc[v.Depth==500,'Rc']=np.nan
        _,_,depths,*_=self.ns['composite_days'](self.sample,v,self.X,self.Y)
        np.testing.assert_equal(depths,[0,200])

    def test_finite_pixel_denominators(self):
        def masked(*args):
            u,v=model(*args)
            if args[-2]==3: u[0,0]=np.nan
            return u,v
        self.ns['esp']=SimpleNamespace(model_uv_at_xy=masked)
        u,v,_,counts,*_=self.ns['composite_days'](self.sample,self.vertical,self.X,self.Y)
        self.assertEqual(counts[0,0,0],2)
        self.assertEqual(counts[1,1,0],3)
        self.assertAlmostEqual(u[0,0,0],3.)
        self.assertAlmostEqual(v[0,0,0],6.)

    def test_notebook_composite_and_plot_cells(self):
        self.ns.update(sample=self.sample,vertical=self.vertical,X=self.X,Y=self.Y,
                       REFERENCE_DEPTH_M=0.,MAX_COMPOSITE_DEPTH_M=None,display=lambda *a:None)
        for c in CELLS:
            source=''.join(c['source'])
            if source.startswith('u_ave, v_ave, depths') or source.startswith('# Plot at one actual'):
                exec(source,self.ns)
        self.assertEqual(self.ns['cnt'],3)
        plt.close('all')

if __name__=='__main__':unittest.main()
