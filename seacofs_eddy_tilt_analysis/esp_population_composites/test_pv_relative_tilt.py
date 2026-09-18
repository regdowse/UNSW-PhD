"""Scientific invariants for the centre-only PV-relative notebook."""
import numpy as np
import pandas as pd
import pytest
import pv_relative_tilt_tools as pvt


def population():
    surface=[];profiles=[]
    for e in range(24):
        for day in [1,2]:
            deep=e%2==0
            surface.append(dict(Eddy=e,Day=day,Cyc='AE' if e%4<2 else 'CE',
                Region=['S1','U2','D1'][e%3],Ro=.2 if day==1 else .8,
                PV_grad_mag=1.,PV_grad_theta=90.,
                PV_grad_plan_mag=[1.,1.,.1][e%3],PV_grad_topo_mag=[.1,1.,1.][e%3]))
            for z in ([0.,200.,500.,1000.,1500.] if deep else [0.,200.,500.]):
                profiles.append(dict(Eddy=e,Day=day,Depth=z,xc=10+(1+e)*z/1000,
                    yc=20+day*z/1000,Rc=10.+day,Omega=1e-5,q11=1.,q12=0.,q22=1.))
    return pd.DataFrame(surface),pd.DataFrame(profiles)


def test_bearings_rotation_and_surface_reference():
    s,v=population()
    s=s.iloc[:1].copy();v=v.loc[v.Eddy.eq(0)&v.Day.eq(1)].copy()
    # Rotate a geographic eastward displacement into a model grid rotated 30deg.
    angle=np.pi/6
    v['xc']=10+np.cos(angle)*v.Depth/100
    v['yc']=20-np.sin(angle)*v.Depth/100
    _,m,_=pvt.prepare_members(s,v,angle)
    np.testing.assert_allclose(m.along_km,m.Depth/100,atol=1e-12)
    np.testing.assert_allclose(m.perp_km,0,atol=1e-12)
    s['PV_grad_theta']=0. # north: east is negative CCW-perpendicular
    _,m,_=pvt.prepare_members(s,v,angle)
    np.testing.assert_allclose(m.along_km,0,atol=1e-12)
    np.testing.assert_allclose(m.perp_km,-m.Depth/100,atol=1e-12)


def test_pairwise_spread_depth_limit_and_pooled_medians():
    s,v=population()
    v.loc[(v.Eddy==0)&(v.Day==1)&(v.Depth==200),'xc']=0
    v.loc[(v.Eddy==0)&(v.Day==1)&(v.Depth==1000),'xc']=30
    v.loc[v.Depth>1000,'xc']=10000 # excluded from spread, retained in tilt profile
    a,m,t=pvt.prepare_members(s,v,0.)
    first=a.loc[(a.Eddy==0)&(a.Day==1)].iloc[0]
    assert first.centre_spread_km>30 # span (-10 to +20), plus y separation
    assert first.centre_spread_km<31
    for extent in ['shallow','deep']:
        d=a.loc[a.extent_group.eq(extent)]
        threshold=t.set_index('extent_group').loc[extent,'spread_threshold_km']
        assert threshold==d.centre_spread_km.median()
        assert (d.loc[d.spread_class.eq('low'),'centre_spread_km']<threshold).all()
        assert (d.loc[d.spread_class.eq('high'),'centre_spread_km']>=threshold).all()
    assert m.loc[m.Depth.eq(1500),'along_km'].max()>9000
    np.testing.assert_allclose(m.along_Rc,m.along_km/m.surface_Rc_km)


def test_exclusions_unknowns_and_dominance():
    s,v=population()
    s.loc[0,'PV_grad_mag']=0
    s.loc[1,'Ro']=np.nan;s.loc[1,'Region']='outside'
    a,m,t=pvt.prepare_members(s,v,0.)
    assert a.iloc[0].analysis_status=='invalid PV reference'
    assert a.iloc[1].analysis_status=='included'
    assert a.iloc[1].Ro_class=='unknown'
    assert set(a.PV_regime)=={'planetary','mixed','topographic'}
    assert len(m.loc[(m.Eddy==0)&(m.Day==2)])>0
    assert len(pvt.select_group(m,('regional','S','deep','high','AE')).loc[lambda x:x.Eddy.eq(0)])==0
    with pytest.raises(ValueError,match='Duplicate'):
        pvt.prepare_members(pd.concat([s,s.iloc[:1]]),v,0.)


def test_no_subsurface_spread_and_empty_population():
    s,v=population();v=v.loc[v.Depth.eq(0)]
    a,m,t=pvt.prepare_members(s,v,0.)
    assert a.spread_class.eq('unknown').all()
    s['PV_grad_mag']=np.nan
    a,m,t=pvt.prepare_members(s,v,0.)
    assert m.empty
    stats,inventory=pvt.analyse(m,n_boot=20)
    assert stats.empty and inventory.n_eddies.eq(0).all()


def test_track_bootstrap_day_weighting_neutral_and_normalisation():
    # Unequal day counts: mean=3, not the equal-track mean=2.
    m=pd.DataFrame(dict(Eddy=[1,2,2,2],Day=[1,1,2,3],Depth=[200.]*4,
        along_km=[0.,4.,4.,4.],perp_km=[2.,-2.,-2.,-2.],
        along_Rc=[0.,1.,1.,1.],perp_Rc=[1.,-.5,-.5,-.5],distance_km=[2.,5.,5.,5.]))
    s=pvt.summarise(m,1000,3).iloc[0]
    assert s.along_km==3 and s.along_Rc==.75
    assert s.along_km_ci_low==0 and s.along_km_ci_high==4
    assert s.neutral_days==1 and s.directional_days==3
    assert s.fraction_with==1 and np.isnan(s.fraction_with_ci_low) # only one directional track
    both=pd.concat([m,m.assign(Depth=500.)])
    out=pvt.summarise(both,100,3)
    assert out.along_km_ci_low.nunique()==1
    assert out.along_km_ci_high.nunique()==1
    pd.testing.assert_frame_equal(out,pvt.summarise(both,100,3))


def test_focused_groups_and_matched_depths():
    a,m,_=pvt.prepare_members(*population(),0.)
    stats,inventory=pvt.analyse(m,n_boot=20)
    assert len(inventory)==56
    assert not inventory.loc[inventory.family.eq('dominance'),'split'].ne('all').any()
    matched,mapping=pvt.matched_profiles(m,{'shallow':[200,500],'deep':[200,500,1000,1500]},0.)
    assert mapping.status.eq('ok').all()
    assert len(matched)==len(m)
    missing,mapping=pvt.matched_profiles(m,{'shallow':[9000]},0.)
    assert missing.empty and mapping.status.ne('ok').any()
    direction,depths=pvt.direction_at_depths(stats,[200,500,1000,1500],0.)
    assert set(direction.Depth)=={200,500,1000,1500}
    assert not len(direction.loc[direction.extent_group.eq('shallow')&direction.Depth.eq(1500)])


def test_notebook_end_to_end(monkeypatch, tmp_path):
    """Execute every code cell, mocking only file access and notebook display."""
    import json, sys
    from pathlib import Path
    from types import SimpleNamespace
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root))
    monkeypatch.chdir(root/'esp_population_composites')
    import seacofs_tilt_tools as tilt
    surface,vertical=population()
    monkeypatch.setattr(tilt,'load_grid',lambda *a:SimpleNamespace(angle=0.))
    monkeypatch.setattr(tilt,'load_tilt_tables',lambda *a,**k:(surface,None))
    monkeypatch.setattr(tilt,'load_vert',lambda *a:vertical)
    monkeypatch.setattr(tilt,'add_pv_gradient_terms',lambda df,*a,**k:df)
    monkeypatch.setattr(tilt,'add_region_labels',lambda df,*a:df)
    monkeypatch.setattr(pd.DataFrame,'to_parquet',lambda df,path,**kw:df.to_pickle(path))
    shown=[]
    def show():
        shown.extend(plt.get_fignums());plt.close('all')
    monkeypatch.setattr(plt,'show',show)
    nb=json.loads((root/'esp_population_composites/pv_relative_tilt.ipynb').read_text())
    ns={}
    for i,c in enumerate(nb['cells']):
        if c['cell_type']!='code':continue
        source=''.join(c['source']).replace('from IPython.display import display','display = lambda *args, **kw: None')
        source=source.replace('BOOTSTRAPS = 500','BOOTSTRAPS = 20').replace('SAVE_RESULTS = False','SAVE_RESULTS = True')
        exec(compile(source,f'cell_{i}','exec'),ns)
        ns['OUTPUT']=tmp_path
    assert len(shown)==16
    assert len(list(tmp_path.iterdir()))==12
    assert not ns['profile_stats'].empty
