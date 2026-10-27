"""Focused drifter/satellite example; original data conventions, local weighting.

Metres, m/s, naive UTC timestamps. External ESP files are read only.
"""
from pathlib import Path
import hashlib
import warnings
import numpy as np
import pandas as pd
import netCDF4
import platform_weighting as pw
import validation_helpers as vh


def load_inputs(drifter_file, satellite_file, cutoff='2024-01-27 06:00', lon0=154.3, lat0=-37.9):
    def numeric(ds, name):
        return np.asarray(np.ma.filled(ds[name][:].astype(float), np.nan))
    with netCDF4.Dataset(drifter_file) as ds:
        if ds['time'].units != 'seconds since 1970-01-01T00:00:00Z':
            raise ValueError('Unexpected drifter time units.')
        d = pd.DataFrame(dict(time=pd.to_datetime(numeric(ds,'time'),unit='s').floor('s'),
            longitude=numeric(ds,'longitude'), latitude=numeric(ds,'latitude'),
            u=numeric(ds,'ve'), v=numeric(ds,'vn'), drifter_id=np.asarray(ds['ID'][:]).astype(str)))
    d = d[(np.hypot(d.u,d.v)<100) & (d.time<=pd.Timestamp(cutoff))].copy()
    d['x']=(d.longitude-lon0)*111000*np.cos(np.deg2rad(lat0))
    d['y']=(d.latitude-lat0)*111000
    d=d[np.isfinite(d[['x','y','u','v']]).all(axis=1)].reset_index(drop=True)
    d['platform']='Drifter';d['depth_m']=15.  # nominal manuscript drogue depth, not measured depth
    with netCDF4.Dataset(satellite_file) as ds:
        if ds['TIME'].units != 'days since 1985-01-01 00:00:00 UTC':
            raise ValueError('Unexpected satellite time units.')
        times=pd.DatetimeIndex(pd.Timestamp('1985-01-01')+pd.to_timedelta(numeric(ds,'TIME'),'D'))
        lon,lat=np.meshgrid(numeric(ds,'LONGITUDE'),numeric(ds,'LATITUDE'))
        uv=np.stack([numeric(ds,'UCUR'),numeric(ds,'VCUR')],axis=-1)
    uv=np.where(abs(uv)<4.9,uv,np.nan)
    first=times[0].normalize()+pd.Timedelta(hours=6)
    if first<times[0]:first+=pd.Timedelta(days=1)
    last=min(pd.Timestamp(cutoff),times[-1])
    daily=pd.date_range(first,last,freq='D')
    satellite=dict(time=times,uv=uv,x=((lon-lon0)*111000*np.cos(np.deg2rad(lat0))).ravel(),
                   y=((lat-lat0)*111000).ravel(),daily=daily)
    info=dict(files={str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest()
                     for p in [drifter_file,satellite_file]},
        drifter_ids=d.drifter_id.unique().tolist(),drifter_start=str(d.time.min()),drifter_end=str(d.time.max()),
        cutoff=str(cutoff),depth_note='Drifter 15 m nominal manuscript depth; satellite surface geostrophic.',
        time_processing='Native satellite linearly interpolated to original daily 06 UTC; drifters within inclusive +/- window.')
    return d,satellite,info


def snapshot(drifters,satellite,time,half_window_hours=24):
    time=pd.Timestamp(time)
    if time not in satellite['daily']:raise ValueError('Choose an original daily 06 UTC timestamp within coverage.')
    dt=pd.Timedelta(hours=half_window_hours)
    if dt<=pd.Timedelta(0):raise ValueError('Half-window must be positive.')
    d=drifters[(drifters.time>=time-dt)&(drifters.time<=time+dt)].copy()
    uv=pw.interpolate_field(satellite['time'],satellite['uv'],time).reshape(-1,2)
    s=pd.DataFrame(dict(x=satellite['x'],y=satellite['y'],u=uv[:,0],v=uv[:,1],
                        time=time,platform='Satellite',depth_m=0.,drifter_id=''))
    s=s[np.isfinite(s[['x','y','u','v']]).all(axis=1)]
    data=pd.concat([d,s],ignore_index=True);data['id']=np.arange(len(data))
    return data


def track_to_day(backend,drifters,satellite,time,source='both',core_radius=35000.,outer_radius=100000.,half_window_hours=24):
    """Replay original cell 8 through selected day, including its lagged Q.

    Both passes use previous day's Q; only centre updates between passes.
    Satellite-only follows the same selection rule independently.
    """
    target=pd.Timestamp(time)
    if target not in satellite['daily']:raise ValueError('Target not on original daily time grid.')
    previous=dict(xc=0.,yc=0.,Q=np.eye(2),Omega=-1.)
    history=[]
    for t in satellite['daily'][satellite['daily']<=target]:
        data=snapshot(drifters,satellite,t,half_window_hours)
        chosen=data if source=='both' else data[data.platform=='Satellite']
        first=chosen[pw.rho2(chosen,previous)<=core_radius**2]
        provisional=pw.backend_inner(backend,first)
        selection=dict(provisional,Q=previous['Q'])
        core=chosen[pw.rho2(chosen,selection)<=core_radius**2]
        inner=pw.backend_inner(backend,core)
        jump=np.hypot(inner['xc']-provisional['xc'],inner['yc']-provisional['yc'])
        if jump>=100000 or inner['Omega']>=0:
            raise pw.FitError(f'Original trajectory acceptance failed at {t}; inspect this day.')
        history.append(dict(time=t,xc=inner['xc'],yc=inner['yc'],Omega=inner['Omega'],
                            alpha=vh.geometry(inner['Q'])[0],n_inner=len(core)))
        previous=inner
    outer=chosen[pw.rho2(chosen,inner)<=outer_radius**2]
    return data,inner,core.id.tolist(),outer.id.tolist(),pd.DataFrame(history)


def finish_fit(data,inner,core_ids,outer_ids,label,*,inner_weight=np.nan,outer_weight=np.nan,
               original=False,inner_only=False,rc_limit=1e6):
    """Keep valid inner estimates when outer fit fails; no outer fit for inner-only reference."""
    alpha,angle=vh.geometry(inner['Q'])
    r=dict(label=label,inner_valid=True,valid=False,status='inner_only',xc_km=inner['xc']/1000,
        yc_km=inner['yc']/1000,alpha=alpha,angle_deg=angle,inner_Omega=inner['Omega'],
        inner_vorticity=inner['Omega']*np.trace(inner['Q']),Rc_km=np.nan,Omega=np.nan,vorticity=np.nan,
        satellite_weight_inner=inner_weight,satellite_weight_outer=outer_weight,
        original=original,peak_bracketed=False,warnings='')
    f=dict(summary=r,inner=inner,inner_ids=core_ids,outer_ids=outer_ids,outer_kept_ids=[])
    for stage,ids in [('inner',core_ids),('outer_selected',outer_ids)]:
        for name in ['Drifter','Satellite']:
            n=int(((data.platform==name)&data.id.isin(ids)).sum())
            r[f'n_{stage}_{name}']=n
    if inner_only:return f
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        result=pw.weighted_outer(data[data.id.isin(outer_ids)],inner,outer_weight,rc_limit,
                                 per_observation=original,reference_platform='Drifter')
    r['warnings']=' | '.join(sorted({str(w.message) for w in caught}))
    f['outer_kept_ids']=result['kept'].id.tolist()
    r.update(status=result['status'],Rc_km=result['Rc']/1000,Omega=result['Omega'],
             psi0=result['psi0'],peak_bracketed=result.get('peak_bracketed',False),
             message=result.get('message',''))
    for name in ['Drifter','Satellite']:
        d=result['kept'][result['kept'].platform==name]
        radii=np.sqrt(pw.rho2(d,inner)) if len(d) else np.array([])
        r[f'n_outer_used_{name}']=len(d)
        r[f'rho_max_km_{name}']=radii.max()/1000 if len(d) else np.nan
    if result['status']=='ok':
        r.update(valid=True,vorticity=result['Omega']*np.trace(inner['Q']))
    return f


def weighted_fit(data,reference,label,inner_weight,outer_weight,rc_limit=1e6):
    """Fixed combined selection; zero outer weight is a diagnostic, not proof of support."""
    try:
        core=data[data.id.isin(reference['inner_ids'])]
        inner=pw.weighted_inner(core,inner_weight,reference_platform='Drifter')
        return finish_fit(data,inner,reference['inner_ids'],reference['outer_ids'],label,
            inner_weight=inner_weight,outer_weight=outer_weight,rc_limit=rc_limit)
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
        return dict(summary=dict(label=label,status=str(exc),inner_valid=False,valid=False,
                    satellite_weight_inner=inner_weight,satellite_weight_outer=outer_weight,warnings=''),
                    inner=None,inner_ids=reference['inner_ids'],outer_ids=reference['outer_ids'],outer_kept_ids=[])


def support(fits,data):
    rows=[]
    for f in fits:
        r=f['summary']
        for stage,key in [('inner','inner_ids'),('outer','outer_kept_ids')]:
            for name in ['Drifter','Satellite']:
                d=data[data.id.isin(f[key])&(data.platform==name)]
                active=len(d)
                weight=r.get(f'satellite_weight_{stage}',np.nan)
                if not r.get('original',False) and name=='Satellite' and weight==0:active=0
                mass=active if r.get('original',False) else ((1. if name=='Drifter' else weight) if active else 0.)
                rows.append(dict(label=r['label'],stage=stage,platform=name,n_selected=len(d),n_active=active,
                    objective_mass=mass,start=d.time.min(),end=d.time.max(),
                    depth_m=d.depth_m.iloc[0] if len(d) else np.nan,
                    n_drifters=d.drifter_id.nunique() if name=='Drifter' else 0))
    return pd.DataFrame(rows)


def coverage(data,reference):
    """One common reference ellipse/peak for both platforms; no sign or weight exclusion."""
    peak=reference['summary']['Rc_km']/np.sqrt(2)
    rows=[]
    for name,d in data.groupby('platform'):
        radii=np.sqrt(pw.rho2(d,reference['inner']))/1000
        rows.append(dict(platform=name,n=len(d),rho_min_km=radii.min(),rho_max_km=radii.max(),
                         reference_peak_km=peak,n_below_peak=int((radii<peak).sum()),
                         n_above_peak=int((radii>peak).sum())))
    return pd.DataFrame(rows)


def plot_coverage(data,reference,core_radius):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for name,d in data.groupby('platform'):
        if name=='Drifter':
            for ident,g in d.groupby('drifter_id'):
                g=g.sort_values('time');ax[0].plot(g.x/1000,g.y/1000,'o-',ms=3,label=f'Drifter {ident}')
        else:ax[0].scatter(d.x/1000,d.y/1000,s=12,label='Satellite',color='0.6')
        radii=np.sqrt(pw.rho2(d,reference['inner']))/1000
        ax[1].hist(radii,bins=np.linspace(0,max(110,radii.max()),23),density=True,
                   histtype='step',lw=2,label=name)
    inner=reference['inner'];e,vec=np.linalg.eigh(inner['Q'])
    theta=np.linspace(0,2*np.pi,200);unit=np.array([np.cos(theta),np.sin(theta)])
    for radius,label,color in [(core_radius/1000,'Core selection scale','tab:green'),(reference['summary']['Rc_km']/np.sqrt(2),'Reference speed peak','tab:red')]:
        ellipse=(vec@np.diag(1/np.sqrt(e))@unit)*radius+np.array([[inner['xc']/1000],[inner['yc']/1000]])
        ax[0].plot(*ellipse,ls='--',lw=1,label=label,color=color);ax[1].axvline(radius,ls='--',lw=1,label=label,color=color)
    ax[0].set(xlabel='East (km)',ylabel='North (km)',aspect='equal',title='Selected-day observations')
    ax[1].set(xlabel='Effective radius in original combined ellipse (km)',ylabel='Density (per km)',title='Radial support; each platform normalised separately')
    for a in ax:a.legend(fontsize=7);a.grid(alpha=.2)
    return fig


def plot_sweep(fits,stage,reference,evaluation):
    import matplotlib.pyplot as plt
    r=pd.DataFrame([f['summary'] for f in fits])
    for key in ['xc_km','yc_km','Rc_km','alpha']:
        if key not in r:r[key]=np.nan
    x=r[f'satellite_weight_{stage}']
    fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    ax[0,0].plot(x,np.hypot(r.xc_km-reference['xc_km'],r.yc_km-reference['yc_km']).where(r.inner_valid),'o-')
    ax[0,0].set_ylabel('Centre shift from original combined (km)')
    for a,col,ylabel in [(ax[0,1],'Rc_km','Gaussian Rc (km)'),(ax[1,0],'alpha','Aspect ratio')]:
        a.plot(x,r[col].where(r.valid if col=='Rc_km' else r.inner_valid),'o-')
        a.axhline(reference[col],ls='--',color='.5');a.set_ylabel(ylabel)
        if col=='Rc_km' and 'peak_bracketed' in r:
            unsupported=r.valid & ~r.peak_bracketed.fillna(False).astype(bool)
            if unsupported.any():
                a.plot(x[unsupported],r.loc[unsupported,col],'rx',ms=9,label='Peak not bracketed')
                a.legend(fontsize=8)
    residuals=pw.residuals(fits,evaluation)
    if len(residuals):
        for name,d in residuals.groupby('platform'):
            ax[1,1].plot(x,d.set_index('label').vector_rmse_ms.reindex(r.label),'o-',label=name)
        ax[1,1].legend()
    ax[1,1].set_ylabel('Common-set vector RMSE (m/s)')
    for a in ax.flat:
        a.set_xscale('symlog',linthresh=.1);a.set_xticks(x);a.set_xticklabels([f'{v:g}' for v in x])
        a.set_xlabel(f'Satellite {stage} weight');a.grid(alpha=.2)
    fig.suptitle(f'{stage.capitalize()} weight sensitivity; other weight fixed; dashed = original combined')
    return fig


def inner_residuals(fits,evaluation):
    """Linear-core residuals on a common core set, separate from full Gaussian errors."""
    rows=[]
    for f in fits:
        if not f['summary']['inner_valid']:continue
        inner=f['inner']
        for name,d in evaluation.groupby('platform'):
            qr=(d[['x','y']].to_numpy()-[inner['xc'],inner['yc']])@inner['Q']
            predicted=inner['Omega']*np.column_stack([-qr[:,1],qr[:,0]])
            error=predicted-d[['u','v']].to_numpy()
            rows.append(dict(label=f['summary']['label'],platform=name,n=len(d),
                vector_rmse_ms=np.sqrt(np.mean(np.sum(error**2,axis=1))),
                bias_u_ms=error[:,0].mean(),bias_v_ms=error[:,1].mean()))
    return pd.DataFrame(rows)
