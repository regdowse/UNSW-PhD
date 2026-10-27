"""Inner-core window sensitivity; stop at first failure, without restart logic."""
import warnings
import numpy as np
import pandas as pd
import drifter_weighting as dw
import platform_weighting as pw
import validation_helpers as vh

CONFIGURATIONS = ['Drifter only', 'Combined original', 'Combined weighted']


def common_times(drifters, satellite, windows):
    windows=np.asarray(windows,float)
    if not len(windows) or not np.isfinite(windows).all() or np.any(windows<=0):
        raise ValueError('Window durations must be finite and positive.')
    margin=pd.Timedelta(hours=float(windows.max())/2)
    times=satellite['daily']
    times=times[(times>=drifters.time.min()+margin)&(times<=drifters.time.max()-margin)]
    if not len(times):raise ValueError('No common dates with the largest full temporal window.')
    return times


def estimate(backend, data, weighted, satellite_weight):
    active=data if not weighted or satellite_weight>0 else data[data.platform=='Drifter']
    if len(active)<3:raise pw.FitError('insufficient_inner_samples')
    xy=active[['x','y']].to_numpy();xy=(xy-xy.mean(axis=0))/1000
    x,y=xy.T;n=len(x)
    design=np.zeros((2*n,5))
    design[:n,0]=1;design[:n,2]=x;design[:n,3]=y
    design[n:,1]=1;design[n:,2]=-y;design[n:,4]=x
    if np.linalg.matrix_rank(design)<5:raise pw.FitError('rank_deficient_inner')
    if weighted:return pw.weighted_inner(data,satellite_weight,reference_platform='Drifter')
    return pw.backend_inner(backend,data)


def run(backend,drifters,satellite,windows,times,*,satellite_weight=.1,core_radius=35000.,
        initial_centre=(0.,0.),max_refinement_shift=100000.):
    """Original drifter all-window fit; combined original lagged-Q two-pass selection.

    All configurations stop after the first failed estimate. Subsequent dates
    remain rows marked not_attempted_after_failure, not additional failed fits.
    Each combined run starts at the same approximate centre and identity Q.
    """
    if not np.isfinite(satellite_weight) or satellite_weight<0:
        raise ValueError('Satellite weight must be finite and nonnegative.')
    rows=[]
    for hours in windows:
        # Same time interpolation and temporal observations for every configuration.
        snapshots={t:dw.snapshot(drifters,satellite,t,float(hours)/2) for t in times}
        for configuration in CONFIGURATIONS:
            previous=dict(xc=initial_centre[0],yc=initial_centre[1],Q=np.eye(2),Omega=-1.)
            stopped=False
            for time in times:
                row=dict(configuration=configuration,window_hours=float(hours),time=time,
                    valid=False,status='not_attempted_after_failure',warnings='',
                    xc_km=np.nan,yc_km=np.nan,alpha=np.nan,angle_deg=np.nan,inner_vorticity=np.nan,
                    n_drifter_available=np.nan,n_drifter_used=np.nan,n_satellite_used=np.nan,n_drifters_used=np.nan)
                rows.append(row)
                if stopped:continue
                data=snapshots[time]
                row['n_drifter_available']=int((data.platform=='Drifter').sum())
                weighted=configuration=='Combined weighted'
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    try:
                        if configuration=='Drifter only':
                            core=data[data.platform=='Drifter']
                            inner=estimate(backend,core,False,satellite_weight)
                        else:
                            first=data[pw.rho2(data,previous)<=core_radius**2]
                            provisional=estimate(backend,first,weighted,satellite_weight)
                            selection=dict(provisional,Q=previous['Q'])
                            core=data[pw.rho2(data,selection)<=core_radius**2]
                            inner=estimate(backend,core,weighted,satellite_weight)
                            if np.hypot(inner['xc']-provisional['xc'],inner['yc']-provisional['yc'])>=max_refinement_shift:
                                raise pw.FitError('original_refinement_shift_limit')
                        if inner['Omega']>=0:raise pw.FitError('unexpected_rotation_sign')
                        alpha,angle=vh.geometry(inner['Q'])
                        active=core if not weighted or satellite_weight>0 else core[core.platform=='Drifter']
                        d=active[active.platform=='Drifter']
                        row.update(valid=True,status='ok',xc_km=inner['xc']/1000,yc_km=inner['yc']/1000,
                            alpha=alpha,angle_deg=angle,inner_vorticity=inner['Omega']*np.trace(inner['Q']),
                            n_drifter_used=len(d),n_satellite_used=int((active.platform=='Satellite').sum()),
                            n_drifters_used=d.drifter_id.nunique())
                        previous=inner
                    except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
                        row.update(status=str(exc));stopped=True
                    finally:row['warnings']=' | '.join(sorted({str(w.message) for w in caught}))
    return pd.DataFrame(rows)


def counts(results):
    return results.groupby(['configuration','window_hours'],sort=False).agg(
        days_with_result=('valid','sum'),comparison_days=('time','size')).reset_index()


def differences(results,reference_hours=48.):
    """Matched successful dates only; denominators expose missing comparisons."""
    rows=[]
    for name,group in results.groupby('configuration',sort=False):
        reference=group[(group.window_hours==reference_hours)&group.valid].set_index('time')
        for hours,d in group.groupby('window_hours',sort=False):
            pair=d[d.valid].set_index('time').join(reference,how='inner',rsuffix='_ref')
            r=dict(configuration=name,window_hours=hours,matched_days=len(pair),
                median_centre_difference_km=np.nan,max_centre_difference_km=np.nan,
                median_abs_alpha_difference=np.nan,median_abs_vorticity_difference_1e5=np.nan)
            if len(pair):
                shift=np.hypot(pair.xc_km-pair.xc_km_ref,pair.yc_km-pair.yc_km_ref)
                r.update(median_centre_difference_km=shift.median(),max_centre_difference_km=shift.max(),
                    median_abs_alpha_difference=(pair.alpha-pair.alpha_ref).abs().median(),
                    median_abs_vorticity_difference_1e5=(pair.inner_vorticity-pair.inner_vorticity_ref).abs().median()*1e5)
            rows.append(r)
    return pd.DataFrame(rows)


def plot_parameters(results):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig,axes=plt.subplots(4,3,figsize=(14,10),sharex=True,sharey='row',constrained_layout=True)
    for j,name in enumerate(CONFIGURATIONS):
        d=results[results.configuration==name]
        for hours,g in d.groupby('window_hours',sort=True):
            g=g.sort_values('time')
            for i,(column,scale) in enumerate([('xc_km',1),('yc_km',1),('alpha',1),('inner_vorticity',1e5)]):
                axes[i,j].plot(g.time,g[column].where(g.valid)*scale,label=f'{hours:g} h')
        axes[0,j].set_title(name)
        for ax in axes[:,j]:
            ax.grid(alpha=.2)
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
        axes[-1,j].set_xlabel('Analysis date (UTC)')
    for ax,label in zip(axes[:,0],['Centre east (km)','Centre north (km)','Aspect ratio',r'Inner vorticity ($10^{-5}$ s$^{-1}$)']):ax.set_ylabel(label)
    axes[0,0].legend(title='Total window',fontsize=8,ncol=2)
    return fig
