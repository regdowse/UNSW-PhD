import numpy as np
import pandas as pd
import region_overlay_tools as rot


def test_caps_full_depth_classification_pooling_and_surface_radius():
    audit=pd.DataFrame(dict(Eddy=[1,2,3],Day=[1,1,1],Region=['S','S','S'],
                            extent_group=['shallow','deep','deep'],Cyc=['AE']*3))
    surface=audit[['Eddy','Day']].assign(Rc=[2.,4.,np.nan])
    m=pd.DataFrame([dict(Eddy=e,Day=1,Depth=z,east_km=float(e)*z/100,
                        north_km=0.,surface_Rc_km=1000.)
                    for e in [1,2,3] for z in ([0,500,1000] if e==1 else [0,500,1000,1500])])
    stats,inventory=rot.compare_regions(m,audit,surface,n_boot=100)
    assert stats.loc[stats.cohort.isin(['all','shallow']),'Depth'].max()==1000
    assert stats.loc[stats.cohort.eq('deep'),'Depth'].max()==1500
    all500=stats.loc[(stats.cohort=='all')&(stats.Depth==500)].set_index('unit')
    assert all500.loc['km','mean_east']==10. # equal-day pooled raw displacements
    assert all500.loc['Rc','mean_east']==2.5 # (5/2 + 10/4)/2, not profile Rc
    assert all500.loc['Rc','n_eddy_days']==2
    assert all500.loc['km','n_eddy_days']==3
    assert inventory.query("cohort == 'all' and region == 'S' and Cyc == 'AE' and unit == 'Rc'").excluded_radius_days.iloc[0]==1
    deep500=stats.query("cohort == 'deep' and Depth == 500 and unit == 'km'")
    assert deep500.n_eddy_days.iloc[0]==2 # retain original deep classification
    pv=m.rename(columns={'east_km':'along_km','north_km':'perp_km'})
    pstats,_=rot.compare_regions(pv,audit,surface,frame='pv',n_boot=100)
    pd.testing.assert_frame_equal(stats,pstats)
