import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import pv_dominance_tilt_tools as pdt


def fixture():
    # Full-rule complement, including factor-two equality and water-depth boundary.
    s=pd.DataFrame(dict(Eddy=range(1,9),Day=1,Cyc=['AE','CE']*4,
        TiltDis=[5.,6.,8.,9.,10.,4.99,8.,8.],h=[3500.,2500.,2500.,3500.,3000.,3500.,np.nan,3500.],
        PV_grad_plan_mag=[2.,1.,2.,1.,2.,2.,2.,np.nan],
        PV_grad_topo_mag=[1.,2.,1.,2.,1.,1.,1.,1.],PV_grad_mag=1.,PV_grad_theta=90.))
    v=pd.DataFrame([dict(Eddy=e,Day=1,Depth=z,xc=100+z/100,yc=50+z/200,
        Rc=20.,Omega=1e-5,q11=1.,q12=0.,q22=1.) for e in s.Eddy for z in [0.,200.,1000.,1500.]])
    return s,v


def test_selection_complement_boundaries_and_actual_cap():
    s,v=fixture();original=s.copy(deep=True)
    audit,m=pdt.prepare(s,v,0.)
    a=audit.set_index('Eddy')
    assert list(a.regime.loc[1:5])==['planetary','topographic','mixed','mixed','mixed']
    assert a.selection_status.loc[1]=='selected' # >=5 inclusive
    assert a.selection_status.loc[6]=='surface TiltDis below threshold'
    assert a.selection_status.loc[7]=='invalid dominance inputs or water depth'
    assert a.regime.loc[8]=='unclassified'
    assert set(m.Eddy)=={1,2,3,4,5}
    assert set(m.Depth)=={0.,200.,1000.} # actual data cut, including boundary
    assert (m.loc[m.Depth.eq(200),'distance_km']<5).all() # retain levels below 5 km
    np.testing.assert_allclose(m.along_km,m.Depth/100)
    np.testing.assert_allclose(m.perp_km,m.Depth/200)
    pd.testing.assert_frame_equal(s,original)


def test_reference_rotation_and_missing_surface_fit():
    s,v=fixture();s['PV_grad_theta']=0.
    # Grid +x points north when angle=pi/2: +x must project positively along.
    a,m=pdt.prepare(s,v,np.pi/2)
    np.testing.assert_allclose(m.along_km,m.Depth/100)
    np.testing.assert_allclose(m.perp_km,m.Depth/200)
    v=v.loc[~(v.Eddy.eq(1)&v.Depth.eq(0))]
    a,m=pdt.prepare(s,v,np.pi/2)
    assert a.set_index('Eddy').selection_status.loc[1]=='no valid surface fit in capped profile'
    assert not m.Eddy.eq(1).any()


def test_empty_selection_and_invalid_total_reference():
    s,v=fixture();s['PV_grad_mag']=0.
    a,m=pdt.prepare(s,v,0.)
    assert m.empty
    stats,inv=pdt.summarise(m)
    assert stats.empty and len(inv)==6 and inv.eddy_days.eq(0).all()


def test_duplicate_days_rejected():
    s,v=fixture()
    with pytest.raises(ValueError,match='unique'):
        pdt.prepare(pd.concat([s,s.iloc[:1]]),v,0.)


def test_notebook_end_to_end(monkeypatch):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root))
    monkeypatch.chdir(root/'esp_population_composites')
    import seacofs_tilt_tools as tilt
    s,v=fixture()
    # Two independent tracks of each polarity in every category.
    surfaces=[];profiles=[]
    for i in range(4):
        ss=s.copy();vv=v.copy();ss.Eddy+=i*10;vv.Eddy+=i*10
        if i%2: ss['Cyc']=ss.Cyc.map({'AE':'CE','CE':'AE'})
        surfaces.append(ss);profiles.append(vv)
    s=pd.concat(surfaces,ignore_index=True);v=pd.concat(profiles,ignore_index=True)
    monkeypatch.setattr(tilt,'load_grid',lambda *a:SimpleNamespace(angle=0.))
    monkeypatch.setattr(tilt,'load_tilt_tables',lambda *a,**k:(s,None))
    monkeypatch.setattr(tilt,'load_vert',lambda *a:v)
    monkeypatch.setattr(tilt,'add_pv_gradient_terms',lambda df,*a,**k:df)
    shown=[]
    def show():
        shown.extend(plt.get_fignums())
        plt.gcf().savefig('/tmp/pv_dominance_preview.png',dpi=100)
        plt.close('all')
    monkeypatch.setattr(plt,'show',show)
    n=json.loads(Path('pv_dominance_tilt.ipynb').read_text());ns={}
    for i,c in enumerate(n['cells']):
        if c['cell_type']!='code':continue
        source=''.join(c['source']).replace('from IPython.display import display','display=lambda *a,**k: None')
        source=source.replace('BOOTSTRAPS = 500','BOOTSTRAPS = 100')
        exec(compile(source,f'cell_{i}','exec'),ns)
    assert len(shown)==1
    assert ns['axes'].shape==(2,3)
    assert ns['centre_stats'].Depth.max()==1000
    assert ns['inventory'].eddies.ge(2).all()
