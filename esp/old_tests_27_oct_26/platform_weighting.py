"""Notebook 03: real-data loading and explicit local weighted ESP fits.

No writes to ESP_zonodo. Coordinates/radii in m, velocities in m/s.
Weights are platform influence choices, not inverse measurement-error variance.
"""
from pathlib import Path
import hashlib
import warnings
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.optimize import curve_fit
import validation_helpers as vh


class FitError(ValueError):
    pass


def interpolate_field(times, values, target):
    """Linear interpolation without extrapolation, preserving missing cells."""
    t=pd.DatetimeIndex(times)
    if len(t)==0 or len(t)!=len(values):
        raise ValueError('Interpolation needs nonempty, matching times and fields.')
    if not t.is_monotonic_increasing or t.has_duplicates:
        raise ValueError('Satellite timestamps must be unique and increasing.')
    target=pd.Timestamp(target)
    if target<t[0] or target>t[-1]:
        raise ValueError(f'Reference time {target} outside satellite coverage {t[0]} to {t[-1]}.')
    k=t.searchsorted(target)
    if k<len(t) and t[k]==target:
        return values[k].copy()
    f=(target-t[k-1])/(t[k]-t[k-1])
    return (1-f)*values[k-1]+f*values[k]


def load_observations(adcp_path, satellite_path, *, start, end, depth=37.,
                      latitude_limit=-37., lon0=154.3, lat0=-37.9):
    """Reproduce the original example including daily-06:00 satellite resampling.

    ADCP [start,end), latitude filter, nearest actual depth; mean time computed
    before velocity finite filtering, as in the original notebook. No translation,
    temporal correction, thinning, or inter-platform bias adjustment.
    """
    adcp_path,satellite_path=Path(adcp_path),Path(satellite_path)
    s=loadmat(adcp_path,simplify_cells=True)
    z=np.asarray(s['z_sadcp']).ravel(); iz=int(np.argmin(abs(z-depth)))
    time=pd.Timestamp('2000-01-01T12:00:00')+pd.to_timedelta(np.asarray(s['tim_sadcp']).ravel()-2451545.,unit='D')
    lat=np.asarray(s['lat_sadcp']).ravel();lon=np.asarray(s['lon_sadcp']).ravel()
    selected=(lat<=latitude_limit)&(time>=pd.Timestamp(start))&(time<pd.Timestamp(end))
    if not selected.any():raise ValueError('No ADCP observations in the selected window/region.')
    reference=time[selected].mean()
    def frame(platform,lon,lat,u,v,time,depth):
        return pd.DataFrame(dict(platform=platform,
            x=(lon-lon0)*111000*np.cos(np.deg2rad(lat0)),y=(lat-lat0)*111000,
            u=u,v=v,time=time,depth_m=depth))
    adcp=frame('ADCP',lon[selected],lat[selected],np.asarray(s['u_sadcp'])[iz,selected],
               np.asarray(s['v_sadcp'])[iz,selected],time[selected],float(z[iz]))
    # netCDF4 is used on Katana. h5py is a read-only fallback for local verification.
    try:
        import netCDF4
        with netCDF4.Dataset(satellite_path) as f:
            arrays={k:np.asarray(np.ma.filled(f[k][:],np.nan)) for k in ('TIME','LONGITUDE','LATITUDE','UCUR','VCUR')}
            units=f['TIME'].units
        reader='netCDF4'
    except ImportError:
        import h5py
        with h5py.File(satellite_path) as f:
            arrays={k:np.asarray(f[k][:]) for k in ('TIME','LONGITUDE','LATITUDE','UCUR','VCUR')}
            units=f['TIME'].attrs['units']
        reader='h5py'
    if isinstance(units,bytes):units=units.decode()
    if units.strip()!='days since 1985-01-01 00:00:00 UTC':
        raise ValueError(f'Unexpected satellite TIME units: {units}')
    st=pd.DatetimeIndex(pd.Timestamp('1985-01-01')+pd.to_timedelta(arrays['TIME'],'D'))
    uv=np.stack((arrays['UCUR'],arrays['VCUR']),axis=-1).astype(float)
    uv=np.where(abs(uv)<4.9,uv,np.nan)  # original sentinel/velocity filter
    # Same two-step time interpolation as original, not a new time-alignment method.
    begin=st[0].normalize()+pd.Timedelta(hours=6)
    if begin<st[0]:begin+=pd.Timedelta(days=1)
    finish=st[-1].normalize()+pd.Timedelta(hours=6)
    if finish>st[-1]:finish-=pd.Timedelta(days=1)
    daily=pd.date_range(begin,finish,freq='D')
    k=daily.searchsorted(reference)
    if len(daily)==0 or (k==0 and reference!=daily[0]) or k==len(daily):
        raise ValueError('ADCP mean time is not bracketed by daily satellite fields.')
    bracket=[daily[k]] if daily[k]==reference else [daily[k-1],daily[k]]
    native_used=set()
    for t in bracket:
        j=st.searchsorted(t)
        native_used.update([st[j]] if st[j]==t else [st[j-1],st[j]])
    fields=np.stack([interpolate_field(st,uv,t) for t in bracket])
    snapshot=interpolate_field(bracket,fields,reference)
    lon_grid,lat_grid=np.meshgrid(arrays['LONGITUDE'],arrays['LATITUDE'])
    if snapshot.shape!=lon_grid.shape+(2,):raise ValueError('Unexpected satellite variable dimensions.')
    sat=frame('Satellite',lon_grid.ravel(),lat_grid.ravel(),snapshot[...,0].ravel(),snapshot[...,1].ravel(),reference,0.)
    observations=pd.concat([adcp,sat],ignore_index=True)
    good=np.isfinite(observations[['x','y','u','v']]).all(axis=1)
    observations=observations.loc[good].reset_index(drop=True)
    observations['id']=np.arange(len(observations))
    provenance=dict(reference_time=str(reference),adcp_pre_qc_count=int(selected.sum()),
        adcp_selected_start=str(time[selected].min()),adcp_selected_end=str(time[selected].max()),
        actual_adcp_depth_m=float(z[iz]),satellite_depth_label='surface geostrophic (nominal 0 m)',
        satellite_native_start=str(st[0]),satellite_native_end=str(st[-1]),
        satellite_daily_bracket=[str(t) for t in bracket],satellite_reader=reader,
        satellite_native_used=[str(t) for t in sorted(native_used)],
        temporal_processing='native fields -> daily 06 UTC -> ADCP mean time; no translation correction',
        files={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (adcp_path,satellite_path)})
    return observations,provenance


def platform_weights(platforms, satellite_weight, *, reference_platform='ADCP'):
    """Per-point objective weights 1/N_reference, lambda/N_satellite.

    ADCP remains the default; notebook 04 explicitly passes Drifter.
    Counts are stage-specific, not effective independent sample sizes.
    """
    if not np.isfinite(satellite_weight) or satellite_weight<0:raise ValueError('Weights must be finite and nonnegative.')
    labels=np.asarray(platforms);w=np.zeros(len(labels))
    if not np.isin(labels,[reference_platform,'Satellite']).all():raise ValueError('Unknown platform label.')
    for name,mass in [(reference_platform,1.),('Satellite',satellite_weight)]:
        m=labels==name
        if m.any():w[m]=mass/m.sum()
    return w


def checked_inner(xc,yc,q,omega):
    q=np.asarray(q,float)
    if q.shape!=(2,2) or not np.isfinite(q).all() or not np.isfinite([xc,yc,omega]).all():
        raise FitError('nonfinite_inner')
    if not np.allclose(q,q.T) or np.linalg.eigvalsh(q).min()<=0 or not np.isclose(np.linalg.det(q),1,rtol=1e-5):
        raise FitError('invalid_inner_geometry')
    if omega==0:raise FitError('zero_rotation')
    return dict(xc=float(xc),yc=float(yc),Q=q,Omega=float(omega))


def weighted_inner(data, satellite_weight, *, reference_platform='ADCP'):
    """Exact weighted affine fit of the SAME divergence-free linear-core model.

    u=a+g11*x+g12*y; v=b+g21*x-g11*y. Solve five linear coefficients,
    then recover centre and Omega,Q. Coordinates are centred/scaled to km only
    internally. No free background is added: the affine intercept encodes centre.
    """
    w=platform_weights(data.platform,satellite_weight,reference_platform=reference_platform);active=w>0
    d=data.loc[active];w=w[active]
    if len(d)<3:raise FitError('insufficient_inner')
    origin=d[['x','y']].mean().to_numpy();xy=(d[['x','y']].to_numpy()-origin)/1000.
    x,y=xy.T;n=len(d)
    design=np.zeros((2*n,5))
    design[:n,0]=1;design[:n,2]=x;design[:n,3]=y
    design[n:,1]=1;design[n:,2]=-y;design[n:,4]=x
    weights=np.sqrt(np.r_[w,w])
    coef,_,rank,_=np.linalg.lstsq(design*weights[:,None],np.r_[d.u,d.v]*weights,rcond=None)
    if rank<5:raise FitError('rank_deficient_inner')
    a,b,g11,g12,g21=coef;A=np.array([[g11,g12],[g21,-g11]])/1000.
    det=np.linalg.det(A)
    if not np.isfinite(det) or det<=0:raise FitError('non_elliptic_inner_gradient')
    omega=np.sign(A[1,0])*np.sqrt(det)
    q=np.array([[A[1,0],-A[0,0]],[-A[0,0],-A[0,1]]])/omega
    centre=origin-np.linalg.solve(A,[a,b])
    return checked_inner(*centre,q,omega)


def backend_inner(backend,data):
    if len(data)<3:raise FitError('insufficient_inner')
    xc,yc,_,q,omega,_=backend.latte(*[data[k].to_numpy() for k in ('x','y','u','v')])
    return checked_inner(xc,yc,q,omega)


def rho2(data,inner):
    r=data[['x','y']].to_numpy()-[inner['xc'],inner['yc']]
    return np.sum((r@inner['Q'])*r,axis=1)


def outer_samples(data,inner):
    """Original tangential normalisation/sign filter; keep original global row IDs."""
    r=data[['x','y']].to_numpy()-[inner['xc'],inner['yc']];qr=r@inner['Q']
    r2=np.sum(qr*r,axis=1);norm2=np.sum(qr**2,axis=1)
    with np.errstate(divide='ignore',invalid='ignore'):
        vt=(-qr[:,1]*data.u.to_numpy()+qr[:,0]*data.v.to_numpy())*np.sqrt(r2)/norm2
    keep=np.isfinite(vt)&np.isfinite(r2)&(r2>0)&(norm2>0)&(vt*inner['Omega']>=0)
    return data.loc[keep].copy(),r2[keep],vt[keep]


def weighted_outer(data,inner,satellite_weight,rc_limit,per_observation=False, *, reference_platform='ADCP'):
    """Original Gaussian profile with explicit local platform weights/statuses."""
    kept,r2,vt=outer_samples(data,inner)
    w=np.ones(len(kept)) if per_observation else platform_weights(kept.platform,satellite_weight,reference_platform=reference_platform)
    active=w>0;kept=kept.loc[active];r2,vt,w=r2[active],vt[active],w[active]
    result=dict(status='insufficient_outer',kept=kept,Rc=np.nan,Omega=np.nan,psi0=np.nan)
    if len(kept)<3 or np.unique(r2).size<3:return result
    rc0=max(float(np.sqrt(r2[np.argmax(abs(vt))])*np.sqrt(2)),1e-6)
    def model(r2,omega,rc):return omega*np.sqrt(r2)*np.exp(-r2/rc**2)
    # Common rescaling leaves the objective minimiser unchanged and stabilises sigma.
    sigma=np.sqrt(w.max()/w)
    try:
        p,_=curve_fit(model,r2,vt,p0=[inner['Omega'],rc0],sigma=sigma,
            bounds=([-np.inf,1e-8],[np.inf,np.inf]),maxfev=10000)
    except (RuntimeError,ValueError,FloatingPointError) as exc:
        result.update(status='outer_optimisation_failed',message=str(exc));return result
    om,rc=p
    status='ok'
    if not np.isfinite(p).all() or rc<=0:status='nonfinite_outer'
    elif rc>rc_limit:status='outer_radius_limit'
    elif np.sign(om)!=np.sign(inner['Omega']):status='rotation_sign_change'
    result.update(status=status,Rc=rc,Omega=om,psi0=-.5*om*rc**2,Rc_seed=rc0,
        peak_bracketed=bool(np.min(r2)<=rc**2/2<=np.max(r2)))
    return result


def fit(backend,data,label,*,source='both',mode='original',inner_weight=0.,outer_weight=1.,
        fixed_masks=None,core_radius=35000.,outer_radius=100000.,rc_limit=1e6):
    """Original adaptive selection or weighted fits on fixed baseline masks.

    Original inner calls use the external backend. All outer fits use the local
    explicit-status function; equal per-observation objective for originals.
    """
    row=dict(label=label,source=source,mode=mode,inner_valid=False,valid=False,status='unattempted',
             satellite_weight_inner=inner_weight if mode=='weighted' else np.nan,
             satellite_weight_outer=outer_weight if mode=='weighted' else np.nan)
    out=dict(summary=row,inner=None,inner_ids=[],outer_ids=[],outer_kept_ids=[])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            chosen=data if source=='both' else data[data.platform==source]
            if mode=='original':
                first=chosen[(chosen.x**2+chosen.y**2)<=core_radius**2]
                provisional=backend_inner(backend,first)
                core=chosen[rho2(chosen,provisional)<=core_radius**2]
                inner=backend_inner(backend,core)
                outer=chosen[rho2(chosen,inner)<=outer_radius**2]
            elif mode=='weighted':
                if fixed_masks is None:raise ValueError('Weighted sensitivity needs fixed masks.')
                core=chosen[chosen.id.isin(fixed_masks['inner_ids'])]
                outer=chosen[chosen.id.isin(fixed_masks['outer_ids'])]
                inner=weighted_inner(core,inner_weight)
            else:raise ValueError('Unknown mode.')
            out.update(inner=inner,inner_ids=core.id.tolist(),outer_ids=outer.id.tolist())
            alpha,angle=vh.geometry(inner['Q'])
            row.update(inner_valid=True,xc_km=inner['xc']/1000,yc_km=inner['yc']/1000,
                alpha=alpha,angle_deg=angle,inner_Omega=inner['Omega'],
                inner_vorticity=inner['Omega']*np.trace(inner['Q']))
            for stage,d in [('inner',core),('outer_selected',outer)]:
                for platform in ('ADCP','Satellite'):row[f'n_{stage}_{platform}']=int((d.platform==platform).sum())
            if mode=='weighted' and inner_weight==0:row['n_inner_Satellite_active']=0
            else:row['n_inner_Satellite_active']=row['n_inner_Satellite']
            result=weighted_outer(outer,inner,outer_weight,rc_limit,per_observation=mode=='original')
            out['outer_kept_ids']=result['kept'].id.tolist()
            for platform in ('ADCP','Satellite'):
                row[f'n_outer_used_{platform}']=int((result['kept'].platform==platform).sum())
            row.update(status=result['status'],Rc_km=result['Rc']/1000,Omega=result['Omega'],
                       psi0=result['psi0'],peak_bracketed=result.get('peak_bracketed',False))
            for key in ('message','Rc_seed'):
                if key in result:row[key]=result[key]
            if result['status']!='ok':return out
            row.update(valid=True,vorticity=result['Omega']*np.trace(inner['Q']),
                       peak_radius_km=result['Rc']/np.sqrt(2)/1000)
            # Support ranges are diagnostic, not formal parameter uncertainty.
            for platform in ('ADCP','Satellite'):
                d=result['kept'][result['kept'].platform==platform]
                radii=np.sqrt(rho2(d,inner)) if len(d) else np.array([])
                row[f'outer_rho_min_km_{platform}']=radii.min()/1000 if len(radii) else np.nan
                row[f'outer_rho_max_km_{platform}']=radii.max()/1000 if len(radii) else np.nan
        except (FitError,ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            row.update(status=str(exc) if isinstance(exc,FitError) else 'exception',message=str(exc))
        finally:row['warnings']=' | '.join(sorted({str(w.message) for w in caught}))
    return out


def residuals(fits,evaluation):
    """Common fixed evaluation observations, unweighted vector errors per platform.

    These are descriptive residuals, not held-out validation or ground-truth errors.
    """
    rows=[]
    for f in fits:
        r=f['summary']
        if not r['valid']:continue
        inner=f['inner']
        for platform,d in evaluation.groupby('platform'):
            predicted=vh.velocity(d[['x','y']].to_numpy(),[inner['xc'],inner['yc']],r['Omega'],r['Rc_km']*1000,inner['Q'])
            error=predicted-d[['u','v']].to_numpy()
            rows.append(dict(label=r['label'],platform=platform,n=len(d),
                vector_rmse_ms=np.sqrt(np.mean(np.sum(error**2,axis=1))),
                bias_u_ms=error[:,0].mean(),bias_v_ms=error[:,1].mean()))
    return pd.DataFrame(rows)


def support_table(data,ids,stage):
    selected=data[data.id.isin(ids)]
    rows=[]
    for name,d in selected.groupby('platform'):
        rows.append(dict(stage=stage,platform=name,n=len(d),start=d.time.min(),end=d.time.max(),
            depth_m=d.depth_m.iloc[0],x_span_km=(d.x.max()-d.x.min())/1000,y_span_km=(d.y.max()-d.y.min())/1000))
    return pd.DataFrame(rows)


def plot_comparison(fits,data):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    colours={'ADCP':'tab:blue','Satellite':'tab:orange'}
    for name,d in data.groupby('platform'):
        axes[0].scatter(d.x/1000,d.y/1000,s=9,alpha=.4,color=colours[name],label=name)
    for f in fits:
        r=f['summary']
        if not r['inner_valid']:continue
        axes[0].scatter(r['xc_km'],r['yc_km'],marker='x',s=65,label=r['label'])
        if r['valid']:
            rho=np.linspace(0,150000,200)
            speed=abs(r['Omega'])*rho*np.exp(-rho**2/(r['Rc_km']*1000)**2)
            axes[1].plot(rho/1000,speed,label=r['label'])
    axes[0].set(xlabel='East (km)',ylabel='North (km)',aspect='equal',title='Common evaluation observations and fitted centres')
    axes[1].set(xlabel='Effective radius (km)',ylabel='Normalised tangential speed (m/s)',title='Fitted Gaussian profiles (curves may extrapolate)')
    for ax in axes:ax.legend(fontsize=7);ax.grid(alpha=.2)
    return fig


def plot_sweep(fits,stage,reference,evaluation):
    import matplotlib.pyplot as plt
    rows=pd.DataFrame([f['summary'] for f in fits])
    for key in ('xc_km','yc_km','Rc_km','alpha'):
        if key not in rows:rows[key]=np.nan
    x=rows[f'satellite_weight_{stage}'].to_numpy()
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    centre=np.hypot(rows.xc_km-reference['xc_km'],rows.yc_km-reference['yc_km']).where(rows.inner_valid)
    axes[0,0].plot(x,centre,'o-');axes[0,0].set_ylabel('Centre shift from original combined (km)')
    axes[0,1].plot(x,rows.Rc_km.where(rows.valid),'o-');axes[0,1].axhline(reference['Rc_km'],color='.5',ls='--');axes[0,1].set_ylabel('Gaussian Rc (km)')
    axes[1,0].plot(x,rows.alpha.where(rows.inner_valid),'o-');axes[1,0].axhline(reference['alpha'],color='.5',ls='--');axes[1,0].set_ylabel('Aspect ratio')
    rr=residuals(fits,evaluation)
    if len(rr):
        for platform in ('ADCP','Satellite'):
            y=rr[rr.platform==platform].set_index('label').vector_rmse_ms.reindex(rows.label)
            axes[1,1].plot(x,y,'o-',label=platform)
        axes[1,1].legend()
    axes[1,1].set_ylabel('Common-set vector RMSE (m/s)')
    for ax in axes.flat:
        ax.set_xscale('symlog',linthresh=.1)
        ax.set_xticks(x);ax.set_xticklabels([f'{v:g}' for v in x]);ax.grid(alpha=.2)
        ax.set_xlabel(f'Satellite weight, {stage} stage')
    fig.suptitle(f'{stage.capitalize()} weight sensitivity; other weight fixed; dashed lines = original combined')
    return fig


def fit_support(fits,data):
    """Per-fit/stage sample support, with active objective weights made explicit."""
    rows=[]
    for f in fits:
        r=f['summary']
        for stage,key in [('inner','inner_ids'),('outer','outer_kept_ids')]:
            selected=data[data.id.isin(f[key])]
            weight=r.get(f'satellite_weight_{stage}',np.nan)
            for name in ('ADCP','Satellite'):
                d=selected[selected.platform==name]
                active=len(d) if r['mode']=='original' or name=='ADCP' or weight>0 else 0
                mass=active if r['mode']=='original' else ((1. if name=='ADCP' else weight) if active else 0.)
                rows.append(dict(label=r['label'],stage=stage,platform=name,n_selected=len(d),n_active=active,
                    objective_mass=mass,weight_per_active_point=mass/active if active else 0.,
                    start=d.time.min(),end=d.time.max(),depth_m=d.depth_m.iloc[0] if len(d) else np.nan,
                    uncertainty='not estimated; influence weights are not error variances'))
    return pd.DataFrame(rows)
