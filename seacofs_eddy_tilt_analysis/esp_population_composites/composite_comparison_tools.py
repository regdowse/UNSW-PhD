"""Shared shelf, summary, full-ESP and vertical-section comparisons."""
import numpy as np
import pandas as pd
import composite_breakdown_tools as cbt

COLORS = {'AE': 'firebrick', 'CE': 'royalblue'}


def shelf_classes(frame, shelf_lon=154.75):
    """Use the previous notebook's shelf_class rule; unknown is never off-shelf."""
    out = frame.copy()
    out['shelf_class'] = np.where(np.isfinite(out.lon),
        np.where(out.lon < shelf_lon, 'On-shelf', 'Off-shelf'), 'Unknown')
    return out


def shelf_composites(profile_audit, centred_profiles, vertical, X, Y, esp,
                     velocity=True, n_boot=500, seed=731):
    """Separate populations before reconstruction, retaining the original frame."""
    import planetary_composite_tools as pct
    import topographic_composite_tools as tct
    out = {}
    for shelf in ['On-shelf', 'Off-shelf']:
        out[shelf] = {}
        for cohort in ['all', 'shallow', 'deep']:
            for cyc in ['AE', 'CE']:
                selected = profile_audit.loc[profile_audit.shelf_class.eq(shelf) &
                    profile_audit.Cyc.eq(cyc) & profile_audit.profile_status.eq('usable')]
                if cohort != 'all': selected = selected.loc[selected.extent_group.eq(cohort)]
                if selected.empty:
                    print(f'{shelf}, {cohort}, {cyc}: no usable profiles')
                    continue
                if velocity:
                    u,v,depths,counts,support,audit,members = tct.composite_days(selected,vertical,X,Y,esp=esp)
                    fields = dict(u=u,v=v,depths=depths,counts=counts,support=support,audit=audit)
                else:
                    members = centred_profiles.merge(selected[['Eddy','Day']],on=['Eddy','Day'],validate='many_to_one')
                    members = tct.rotate_members(members[['Eddy','Day','Depth','xc','yc']],selected)
                    fields = {}
                stats,boot = tct.summarise(members,n_boot,seed)
                out[shelf][cohort,cyc] = dict(members=members,stats=stats,bootstrap=boot,**fields)
                print(f'{shelf}, {cohort}, {cyc}: {members.Eddy.nunique()} contributing eddies')
    return out


def plot_tilt_comparison(population, results, components=('along','perp'), min_eddies=2):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(11,8),constrained_layout=True)
    metrics=[f'mean_{c}' for c in components]+['distance_km']
    for i,cohort in enumerate(['shallow','deep']):
        for cyc in ['AE','CE']:
            r=results.get((cohort,cyc))
            if r is None: continue
            s=r['stats'];ok=s.n_eddies.ge(min_eddies)
            for ax,col in zip(axs[i],metrics):
                lo,hi=(('distance_ci_low','distance_ci_high') if col=='distance_km' else (col+'_ci_low',col+'_ci_high'))
                ax.plot(s[col].where(ok),s.Depth,color=COLORS[cyc],label=cyc)
                ax.fill_betweenx(s.Depth,s[lo].where(ok),s[hi].where(ok),color=COLORS[cyc],alpha=.18)
        for ax,col in zip(axs[i],metrics):
            ax.set(xlabel=col.replace('mean_','').replace('_km','')+' (km)',ylabel='Depth (m)',title=cohort)
            ax.invert_yaxis();ax.axvline(0,color='.5',lw=.5)
            if ax.lines:ax.legend()
    for j in range(3):
        limits=[ax.get_xlim() for ax in axs[:,j]]
        for ax in axs[:,j]:ax.set_xlim(min(x[0] for x in limits),max(x[1] for x in limits))
    fig.suptitle(f'{population}: mean constituent tilt and 95% eddy-bootstrap CI')
    return fig


def tilt_summary(collections, target_depth=500., depth_tolerance=100., min_eddies=2):
    """One shared exact depth plus each group's own deepest supported depth.

    The shared depth is selected once across all available shallow/deep groups,
    within the stated target tolerance. No per-group nearest-depth substitution.
    """
    items=[(p,g,c,r) for p,rs in collections.items() for (g,c),r in rs.items() if g in ('shallow','deep')]
    common=None
    for _,_,_,r in items:
        s=r['stats']
        z=s.loc[s.n_eddies.ge(min_eddies)&np.isfinite(s.distance_km)&s.Depth.gt(0),'Depth'].to_numpy()
        common=z if common is None else np.intersect1d(common,z)
    shared=np.nan
    if common is not None and len(common):
        candidate=float(common[np.argmin(abs(common-target_depth))])
        if abs(candidate-target_depth)<=depth_tolerance:shared=candidate
    rows=[]
    for p,g,c,r in items:
        s=r['stats'];m=r['members']
        valid=s.loc[s.n_eddies.ge(min_eddies)&np.isfinite(s.distance_km)&s.Depth.gt(0)].sort_values('Depth')
        for kind in ['Shared reference depth','Deepest supported depth']:
            take=(s.loc[s.Depth.eq(shared)] if kind=='Shared reference depth' else valid.tail(1))
            base=dict(population=p,cohort=g,Cyc=c,summary=kind,
                total_eddy_days=len(m[['Eddy','Day']].drop_duplicates()),total_eddies=m.Eddy.nunique(),
                Depth=np.nan,distance_km=np.nan,ci_low=np.nan,ci_high=np.nan,
                mean_member_distance_km=np.nan,n_eddy_days=0,n_eddies=0,
                status='no shared supported depth within tolerance' if kind=='Shared reference depth' else 'no supported subsurface depth')
            if len(take):
                row=take.iloc[0]
                base.update(Depth=row.Depth,distance_km=row.distance_km,ci_low=row.distance_ci_low,
                    ci_high=row.distance_ci_high,mean_member_distance_km=row.mean_member_distance_km,
                    n_eddy_days=int(row.n_eddy_days),n_eddies=int(row.n_eddies),status='ok')
            rows.append(base)
    return pd.DataFrame(rows)


def plot_tilt_summary(summary):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(13,max(4,len(summary)/5)),constrained_layout=True)
    for ax,kind in zip(axs,['Shared reference depth','Deepest supported depth']):
        part=summary.loc[summary.summary.eq(kind)].reset_index(drop=True)
        for i,row in part.iterrows():
            if row.status!='ok':continue
            ax.plot(row.distance_km,i,'o',color=COLORS[row.Cyc])
            if np.isfinite([row.ci_low,row.ci_high]).all():
                ax.hlines(i,row.ci_low,row.ci_high,color=COLORS[row.Cyc])
        labels=[f'{r.population} | {r.cohort} {r.Cyc} | '+(f'{r.Depth:g} m' if r.status=='ok' else 'unavailable') for r in part.itertuples()]
        ax.set(yticks=np.arange(len(part)),yticklabels=labels,xlabel='Mean-vector tilt distance (km)',title=kind)
        ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
    values=summary[['distance_km','ci_high']].to_numpy(float)
    finite=values[np.isfinite(values)]
    xmax=max(float(finite.max()) if len(finite) else 1.,1.)*1.1
    for ax in axs:
        ax.set_xlim(0,xmax);ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('Surface-to-depth composite tilt; bars = pointwise 95% bootstrap CI')
    return fig, axs


def fit_collections(collections,X,Y,esp,min_eddies=2,**fit_kwargs):
    """Fit the composite field itself, storing diagnostics in each result dict."""
    tables=[];audit=[]
    for population,rs in collections.items():
        for (cohort,cyc),r in rs.items():
            if cohort not in ('shallow','deep'):continue
            meta=dict(population=population,cohort=cohort,Cyc=cyc)
            if 'u' not in r:
                audit.append(dict(**meta,status='velocity reconstruction disabled',successful_depths=0,attempted_depths=0))
                continue
            support=r['stats'].set_index('Depth').n_eddies.reindex(r['depths']).fillna(0).to_numpy()
            counts=np.where(support[None,None,:]>=min_eddies,r['counts'],0)
            fits,_,_=cbt.fit_composite(X,Y,r['u'],r['v'],r['depths'],counts,esp,
                expected_sign=1 if cyc=='AE' else -1,**fit_kwargs)
            fits=fits.merge(r['stats'][['Depth','n_eddy_days','n_eddies']],on='Depth',how='left',validate='one_to_one')
            fits['fit_eligible']=support>=min_eddies
            fits.loc[~fits.fit_eligible,'reason']='fewer than minimum contributing eddies'
            r['composite_fit']=fits
            tables.append(fits.assign(**meta))
            audit.append(dict(**meta,status='attempted',successful_depths=int(fits.fit_ok.sum()),attempted_depths=int((support>=min_eddies).sum())))
    return (pd.concat(tables,ignore_index=True) if tables else pd.DataFrame(),pd.DataFrame(audit))


# def plot_fit_profiles(fits,population):
#     import matplotlib.pyplot as plt
#     fig,axs=plt.subplots(2,4,figsize=(14,8),constrained_layout=True)
#     for i,cohort in enumerate(['shallow','deep']):
#         for cyc in ['AE','CE']:
#             s=fits.loc[fits.population.eq(population)&fits.cohort.eq(cohort)&fits.Cyc.eq(cyc)].sort_values('Depth')
#             for ax,col,scale in zip(axs[i],['w','Omega','Rc','vector_r2'],[1e5,1e5,1,1]):
#                 ax.plot((s[col]*scale).where(s.fit_ok),s.Depth,color=COLORS[cyc],label=cyc)
#         for ax,label in zip(axs[i],[r'$\zeta$ ($10^{-5}$ s$^{-1}$)',r'$\Omega$ ($10^{-5}$ s$^{-1}$)',r'$R_c$ (km)','Outer-fit vector R²']):
#             ax.set(xlabel=label,ylabel='Depth (m)',title=cohort);ax.invert_yaxis();ax.axvline(0,color='.5',lw=.5)
#             if ax.lines:ax.legend()
#     for j in range(4):
#         limits=[ax.get_xlim() for ax in axs[:,j]]
#         for ax in axs[:,j]:ax.set_xlim(min(x[0] for x in limits),max(x[1] for x in limits))
#     fig.suptitle(f'{population}: full inner/outer fits to composite velocities (gaps = failed/unsupported fits)')
#     return fig
def plot_fit_profiles(fits,population):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(11,8),constrained_layout=True)
    for i,cohort in enumerate(['shallow','deep']):
        for cyc in ['AE','CE']:
            s=fits.loc[fits.population.eq(population)&fits.cohort.eq(cohort)&fits.Cyc.eq(cyc)].sort_values('Depth')
            for ax,col,scale in zip(axs[i],['w','Omega','Rc'],[1e5,1e5,1]):
                if col in ['w','Omega']:
                    ax.plot(np.abs((s[col]*scale).where(s.fit_ok)),s.Depth,color=COLORS[cyc],label=cyc)
                else:
                    ax.plot((s[col]*scale).where(s.fit_ok),s.Depth,color=COLORS[cyc],label=cyc)
        for ax,label in zip(axs[i],[r'$|\zeta|$ ($10^{-5}$ s$^{-1}$)',r'$|\Omega|$ ($10^{-5}$ s$^{-1}$)',r'$R_c$ (km)']):
            ax.set(xlabel=label,ylabel='Depth (m)',title=cohort)
            ax.invert_yaxis()
            ax.axvline(0,color='.5',lw=.5)
            if ax.lines: ax.legend()
    for j in range(3):
        limits=[ax.get_xlim() for ax in axs[:,j]]
        for ax in axs[:,j]:
            ax.set_xlim(min(x[0] for x in limits),max(x[1] for x in limits))
    fig.suptitle(f'{population}: full inner/outer fits to composite velocities (gaps = failed/unsupported fits)')
    return fig

def section_data(result,X,Y,rotation_rad=0.,min_eddies=2):
    """Two surface-centred vertical cuts in a properly rotated coordinate frame.

    Positive rotation maps input vector coordinates into the display frame.
    For planetary grid -> east/north, use grid.angle; topo already aligned -> 0.
    Returns transverse velocity on x cut and x velocity on y cut, [depth,coord].
    Bilinear interpolation is horizontal only; missing/unsupported stays NaN.
    """
    from scipy.interpolate import RegularGridInterpolator
    x,y=X[0,:],Y[:,0]
    if not np.allclose(X,x[None,:]) or not np.allclose(Y,y[:,None]):
        raise ValueError('Expected meshgrid [y,x] coordinates')
    c,s=np.cos(rotation_rad),np.sin(rotation_rad)
    # Target (x,0) -> input (c*x,-s*x); target (0,y) -> input (s*y,c*y).
    points_x=np.column_stack((-s*x,c*x))
    points_y=np.column_stack((c*y,s*y))
    z=result['depths'];support=result['stats'].set_index('Depth').n_eddies.reindex(z).fillna(0).to_numpy()
    mask=(result['counts']>0)&(support[None,None,:]>=min_eddies)
    u=np.where(mask,result['u'],np.nan);v=np.where(mask,result['v'],np.nan)
    ui=RegularGridInterpolator((y,x),u,bounds_error=False,fill_value=np.nan)
    vi=RegularGridInterpolator((y,x),v,bounds_error=False,fill_value=np.nan)
    # Rotate sampled vectors as well as the sampling locations.
    xcut=(s*ui(points_x)+c*vi(points_x)).T
    ycut=(c*ui(points_y)-s*vi(points_y)).T
    return dict(x=x,y=y,depths=z,xcut=xcut,ycut=ycut)


def plot_sections(results,X,Y,population,cohort,rotation_rad=0.,frame='geographic',min_eddies=2):
    import matplotlib.pyplot as plt
    cuts={cyc:section_data(results[cohort,cyc],X,Y,rotation_rad,min_eddies)
          for cyc in ['AE','CE'] if (cohort,cyc) in results and 'u' in results[cohort,cyc]}
    if not cuts:return None
    finite=[a[np.isfinite(a)] for d in cuts.values() for a in [d['xcut'],d['ycut']]]
    limit=max([float(np.max(abs(a))) for a in finite if len(a)]+[1e-12])
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    labels=(('Zonal (km)','Meridional velocity'),('Meridional (km)','Zonal velocity')) if frame=='geographic' else (
        ('Along-gradient (km)','Perpendicular velocity'),('Perpendicular (km)','Along-gradient velocity'))
    deepest=max(float(np.max(d['depths'])) for d in cuts.values())
    im=None
    for i,cyc in enumerate(['AE','CE']):
        for ax,axis,key,(xlabel,title) in zip(axs[i],['x','y'],['xcut','ycut'],labels):
            ax.set(xlabel=xlabel,ylabel='Depth (m)',title=f'{cyc}: {title}')
            if cyc not in cuts:
                ax.text(.5,.5,'No composite',ha='center',transform=ax.transAxes);continue
            d=cuts[cyc];section=d[key]
            im=ax.pcolormesh(d[axis],d['depths'],section,shading='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
            vals=section[np.isfinite(section)]
            if len(d['depths'])>1 and len(vals) and vals.min()<0<vals.max():
                ax.contour(d[axis],d['depths'],section,levels=[0],colors='k',linewidths=.6)
            ax.set_ylim(max(deepest,1.),0);ax.axvline(0,color='.5',ls=':',lw=.6)
    if im is not None:fig.colorbar(im,ax=list(axs.flat),label='Composite velocity (m/s)')
    fig.suptitle(f'{population}, {cohort}: vertical cuts through the surface reference (black = zero velocity)')
    return fig
