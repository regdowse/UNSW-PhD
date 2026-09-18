"""Checks for contributing-member Ro splits and pooled S/U/D weighting."""
import unittest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import composite_comparison_tools as cc
import regional_composite_tools as rct
import planetary_composite_tools as pct
from types import SimpleNamespace


class RossbyRegionTests(unittest.TestCase):
    def test_split_members_preserves_signed_threshold_and_frame(self):
        meta=pd.DataFrame(dict(Eddy=[1,2,3,4],Day=[1]*4,Ro=[-.49,-.5,.8,np.nan]))
        members=pd.DataFrame(dict(Eddy=[1,2,3,4],Day=[1]*4,Depth=[500.]*4,
            xc=[10,20,30,40],yc=[0,0,0,0],along_km=[0,0,0,0],perp_km=[10,20,30,40]))
        members=pct.geographic_members(members,0.)
        original={('shallow','AE'):dict(members=members)}
        split,audit=cc.rossby_centre_results(original,meta,n_boot=100)
        self.assertEqual(len(split['shallow','low','AE']['members']),1)
        self.assertEqual(len(split['shallow','high','AE']['members']),2)
        self.assertAlmostEqual(split['shallow','high','AE']['stats'].mean_east.iloc[0],25.)
        self.assertTrue(audit.unknown_Ro_days.eq(1).all())
        topo,_=cc.rossby_centre_results(original,meta,frame='topographic',n_boot=100)
        self.assertAlmostEqual(topo['shallow','high','AE']['stats'].mean_perp.iloc[0],25.)
        self.assertEqual(topo['shallow','high','AE']['stats'].mean_along.iloc[0],0.)
        fig,axs=cc.plot_rossby_centres(split,min_eddies=1)
        labelled={line.get_label():line.get_linestyle() for line in axs[0,0].lines}
        self.assertEqual(labelled['AE low |Ro|'],'-');self.assertEqual(labelled['AE high |Ro|'],'--')
        fig.canvas.draw();plt.close('all')

    def test_pool_days_not_means_and_retain_unknown_in_pooled(self):
        s=pd.DataFrame(dict(Eddy=[1,1,1,2],Day=[1,2,3,1],Region=['S1']*3+['S2'],
            Cyc=['AE']*4,Ro=[.2,-.3,-.5,np.nan]))
        v=pd.DataFrame([dict(Eddy=e,Day=d,Depth=z,xc=100+(offset if z else 0),yc=100,
            Rc=20.,Omega=1e-5,q11=1.,q12=0.,q22=1.)
            for (e,d,offset) in [(1,1,0),(1,2,2),(1,3,4),(2,1,10)] for z in [0,500]])
        s=rct.combine_region_labels(s)
        _,audit,profiles=rct.prepare_population(s,v,regions=rct.REGION_GROUPS)
        inv=rct.inventory(audit,regions=rct.REGION_GROUPS,pool_rossby=True)
        self.assertEqual(len(inv),42)
        self.assertEqual(inv.loc[inv.Region.eq('S')&inv.cohort.eq('shallow')&inv.Ro_class.eq('all')&inv.Cyc.eq('AE'),'eddy_days'].iloc[0],4)
        X,Y=np.meshgrid([-1.,0,1],[-1.,0,1])
        results,_=rct.build_composites(audit,profiles,X,Y,None,0.,velocity=False,pool_rossby=True,progress_every=0)
        rct.add_statistics(results,100)
        pooled=results['S','shallow','all','AE']['stats'].set_index('Depth').loc[500]
        self.assertAlmostEqual(pooled.mean_east,4.)
        self.assertNotAlmostEqual(pooled.mean_east,6.)
        self.assertEqual(pooled.n_eddies,2)
        for norm in [False,True]:
            for split in [False,True]:
                fig,axs=rct.plot_combined_regions(results,split_rossby=split,normalised=norm,min_eddies=1)
                self.assertEqual(axs.shape,(2,3));self.assertEqual(axs[0,0].get_title(),'S — shallow')
                self.assertEqual(axs[1,2].get_title(),'D — deep')
                fig.canvas.draw();plt.close(fig)

if __name__=='__main__':unittest.main()
