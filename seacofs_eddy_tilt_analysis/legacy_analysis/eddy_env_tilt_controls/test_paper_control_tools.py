import unittest
import numpy as np
import pandas as pd
import paper_control_tools as p


class PaperControlTests(unittest.TestCase):
    def test_fit_screen_and_axis_ratio(self):
        surface=pd.DataFrame(dict(Eddy=[1,2,3],Day=[1,1,1],Cyc=['AE','CE','AE'],TiltDis=[2.,3.,4.]))
        vertical=pd.DataFrame(dict(Eddy=[1,2,3],Day=[1,1,1],Depth=[0,0,0],Rc=[30,30000,30],
                                   Omega=[1e-5,-1e-5,1e-5],q11=[1,1,-1],q12=[0,0,0],q22=[4,4,4]))
        clean,audit=p.prepare_depths(vertical,surface)
        self.assertEqual(clean.Eddy.tolist(),[1]);self.assertEqual(clean.AR.iloc[0],2)
        self.assertEqual(audit.bad_radius.sum(),1);self.assertEqual(audit.bad_ellipse.sum(),1)

    def test_matching_uses_eddy_and_day(self):
        d=pd.DataFrame(dict(Eddy=[1,1,2],Day=[1,1,1],Depth=[0,100,0],
                            abs_Omega=[1,2,1],AR=[1,2,1],TiltDis=[3,3,4]))
        out=p.matched_depths(d,[0,100])
        self.assertEqual(out.Eddy.tolist(),[1,1])
        d.loc[1,'AR']=np.nan
        self.assertTrue(p.matched_depths(d,[0,100]).empty)

    def test_rank_is_within_polarity_not_pooled(self):
        rows=[]
        for cyc in ['AE','CE']:
            for i in range(12):
                for z in [0,100]:
                    rows.append(dict(Eddy=i+(100 if cyc=='CE' else 0),Day=1,Depth=z,Cyc=cyc,
                                     TiltDis=i,abs_Omega=12-i,AR=i+1))
        rank=p.rank_depths(pd.DataFrame(rows),min_eddies=10)
        self.assertTrue(np.allclose(rank.mean_abs_rho,1))
        self.assertTrue(np.allclose(rank.loc[rank.metric.eq('abs_Omega'),'CE_rho'],-1))
        selected=p.select_depths(rank,n=2,min_spacing=75)
        self.assertEqual(selected['AR'],[0.,100.])

    def test_upper_bin_edge_and_suppressed_gaps(self):
        d=pd.DataFrame(dict(Eddy=range(8),Cyc=['AE']*8,x=range(8),TiltDis=range(8)))
        t=p.median_table(d,'x',[0,3,7],min_eddies=1)
        self.assertEqual(t.eddy_days.sum(),8)
        self.assertAlmostEqual(t.loc[t.Cyc.eq('AE')].iloc[-1]['median'],5.5)
        sparse=p.median_table(d,'x',[0,1,7],min_eddies=3)
        self.assertTrue(np.isnan(sparse.iloc[0]['median']))
        self.assertTrue(sparse.iloc[1].displayed)

    def test_zero_pv_excluded_from_both_panels(self):
        surface=pd.DataFrame(dict(Eddy=[1,2],Day=[1,1],Cyc=['AE','CE'],TiltDis=[4,5]))
        pv=pd.DataFrame(dict(Eddy=[1,2],Day=[1,1],beta=[2e-11,2e-11],PV_grad_mag=[0,1e-13]))
        d,audit=p.environment_data(surface,pv)
        self.assertEqual(d.Eddy.tolist(),[2]);self.assertEqual(d.log10_PV_grad_mag.iloc[0],-13)
        self.assertEqual(audit.kept_days.sum(),1)

    def test_duplicates_rejected(self):
        d=pd.DataFrame(dict(Eddy=[1,1],Day=[1,1]))
        with self.assertRaises(ValueError):p.check_keys(d)

    def test_candidate_support_and_unique_resolution(self):
        rows=[]
        for z in [0,105,515]:
            for cyc in ['AE','CE']:
                for i in range(12 if z<500 else 3):
                    rows.append(dict(Eddy=i+(100 if cyc=='CE' else 0),Day=1,Cyc=cyc,Depth=z))
        a=p.candidate_depths(pd.DataFrame(rows),[0,100,110,500,1000],tolerance=65,min_eddies=5)
        self.assertEqual(a.eligible.tolist(),[True,True,False,False,False])

    def test_override_must_use_eligible_exact_depth(self):
        ranking=pd.DataFrame(dict(metric=['abs_Omega','AR'],depth_m=[105.6,105.6],mean_abs_rho=[.5,.5]))
        with self.assertRaises(ValueError):p.select_depths(ranking,override={'abs_Omega':[100],'AR':[100]})

if __name__=='__main__':unittest.main()
