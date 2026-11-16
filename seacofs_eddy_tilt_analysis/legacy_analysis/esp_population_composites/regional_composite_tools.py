"""Six-region composites with day-level depth/Rossby groups and shared support."""
from itertools import product
import numpy as np
import pandas as pd
import planetary_composite_tools as pct
import composite_comparison_tools as ccomp

REGIONS = ('S1','S2','U1','U2','D1','D2')
KEYS = ('Region','cohort','Ro_class','Cyc')
COLORS = {'AE':'firebrick','CE':'royalblue'}


def prepare_population(surface, vertical, split_depth=1000., ro_split=.5, dominance_factor=2., shelf_lon=154.75,regions=REGIONS):
    """Classify days, retaining missing Ro in baselines and auditing exclusions."""
    if surface[['Eddy','Day']].duplicated().any():raise ValueError('Duplicate surface Eddy-Day keys')
    s=surface.copy()
    s['Ro_abs']=np.abs(s.Ro)
    s['Ro_class']=np.where(np.isfinite(s.Ro_abs),np.where(s.Ro_abs<ro_split,'low','high'),'unknown')
    ratio=s.topo_plan_ratio if 'topo_plan_ratio' in s else pd.Series(np.nan,index=s.index)
    s['PV_regime']=np.select([ratio.le(-np.log(dominance_factor)),ratio.ge(np.log(dominance_factor)),ratio.notna()],
                            ['planetary','topographic','mixed'],default='unknown')
    s['shelf_class']=np.where(np.isfinite(s.lon),np.where(s.lon<shelf_lon,'On-shelf','Off-shelf'),'Unknown') if 'lon' in s else 'Unknown'
    s['selection_status']=np.select([~s.Region.isin(regions),~s.Cyc.isin(['AE','CE'])],
        ['outside named regions','unrecognised polarity'],default='selected')
    audit,profiles=pct.prepare_profiles(s.loc[s.selection_status.eq('selected')],vertical,split_depth)
    radius=profiles.loc[profiles.Depth.eq(0),['Eddy','Day','Rc']].rename(columns={'Rc':'surface_Rc_km'})
    profiles=profiles.merge(radius,on=['Eddy','Day'],validate='many_to_one')
    audit=audit.merge(radius,on=['Eddy','Day'],how='left',validate='one_to_one')
    return s,audit,profiles


def bottom_audit(profiles,grid):
    """Nearest model-cell bathymetry at each original fitted centre; diagnostic only."""
    from scipy.spatial import cKDTree
    p=profiles.copy()
    x=(p.xc+p.surface_xc).to_numpy();y=(p.yc+p.surface_yc).to_numpy()
    inside=(x>=np.nanmin(grid.X_grid))&(x<=np.nanmax(grid.X_grid))&(y>=np.nanmin(grid.Y_grid))&(y<=np.nanmax(grid.Y_grid))
    tree=cKDTree(np.column_stack((grid.X_grid.ravel(),grid.Y_grid.ravel())))
    _,idx=tree.query(np.column_stack((x,y)))
    h=np.asarray(grid.h).ravel()[idx].astype(float)
    valid=inside&np.asarray(grid.mask_rho).ravel()[idx].astype(bool)&np.isfinite(h)&(h>0)
    p['centre_bottom_m']=np.where(valid,h,np.nan)
    p['below_local_bottom']=np.where(valid,p.Depth.to_numpy()>h,np.nan)
    return p


def memberships(row,pool_rossby=False):
    groups=[(row.Region,'all','all',row.Cyc)]
    if pool_rossby:groups.append((row.Region,row.extent_group,'all',row.Cyc))
    if row.Ro_class in ('low','high'):
        groups.append((row.Region,row.extent_group,row.Ro_class,row.Cyc))
    return groups


def inventory(audit,regions=REGIONS,pool_rossby=False):
    """Requested baseline/detailed slots plus optional pooled-depth slots, including empties."""
    rows=[];usable=audit.loc[audit.profile_status.eq('usable')]
    combos=[(r,'all','all',c) for r,c in product(regions,['AE','CE'])]
    combos += list(product(regions,['shallow','deep'],['low','high'],['AE','CE']))
    if pool_rossby:combos += list(product(regions,['shallow','deep'],['all'],['AE','CE']))
    for region,cohort,ro,cyc in combos:
        d=usable.loc[usable.Region.eq(region)&usable.Cyc.eq(cyc)]
        if cohort!='all':d=d.loc[d.extent_group.eq(cohort)]
        if ro!='all':d=d.loc[d.Ro_class.eq(ro)]
        rows.append(dict(Region=region,cohort=cohort,Ro_class=ro,Cyc=cyc,
            eddy_days=len(d),eddies=d.Eddy.nunique(),unknown_Ro_days=int(d.Ro_class.eq('unknown').sum()),
            median_max_depth_m=d.max_fit_depth_m.median(),median_abs_Ro=d.Ro_abs.median(),
            median_surface_Rc_km=d.surface_Rc_km.median()))
    return pd.DataFrame(rows)


def build_composites(audit,profiles,X,Y,esp,grid_angle,velocity=True,progress_every=1000,pool_rossby=False):
    """Reconstruct each day once, accumulating its baseline and detailed group.

    Centres are recorded only for levels contributing finite paired velocity
    somewhere, exactly as in composite_breakdown. Without velocity, valid fitted
    levels are used and explicitly labelled. Inputs profiles are surface-centred.
    """
    meta=audit.loc[audit.profile_status.eq('usable')].set_index(['Eddy','Day'])
    depths=np.sort(profiles.Depth.unique());shape=(*X.shape,len(depths))
    accum={};day_audit=[]
    for n,((eddy,day),p) in enumerate(profiles.groupby(['Eddy','Day'],sort=True),1):
        row=meta.loc[eddy,day];groups=memberships(row,pool_rossby);p=p.sort_values('Depth')
        for key in groups:
            if key not in accum:
                accum[key]=dict(records=[],depth_days=np.zeros(len(depths),dtype=np.int64))
                if velocity:accum[key].update(ut=np.zeros(shape),vt=np.zeros(shape),counts=np.zeros(shape,dtype=np.int64))
        if velocity:u,v=pct.reconstruct(p,X,Y,esp)
        levels=0
        for j,record in enumerate(p.itertuples(index=False)):
            k=np.searchsorted(depths,record.Depth)
            finite=(np.isfinite(u[...,j])&np.isfinite(v[...,j])) if velocity else None
            if velocity and not finite.any():continue
            item=dict(Eddy=eddy,Day=day,Depth=float(record.Depth),xc=float(record.xc),yc=float(record.yc),
                surface_Rc_km=float(record.surface_Rc_km),
                centre_bottom_m=getattr(record,'centre_bottom_m',np.nan),
                below_local_bottom=getattr(record,'below_local_bottom',np.nan))
            levels+=1
            for key in groups:
                a=accum[key];a['records'].append(item);a['depth_days'][k]+=1
                if velocity:
                    a['ut'][...,k]+=np.where(finite,u[...,j],0.)
                    a['vt'][...,k]+=np.where(finite,v[...,j],0.)
                    a['counts'][...,k]+=finite
        day_audit.append(dict(Eddy=eddy,Day=day,contributing_levels=levels,
                             reconstruction_status='included' if levels else 'no finite paired ESP velocity'))
        if progress_every and n%progress_every==0:print(f'Reconstructed {n} unique eddy-days',flush=True)
    results={}
    for key,a in accum.items():
        if not a['records']:continue
        m=pct.geographic_members(pd.DataFrame(a.pop('records')),grid_angle)
        active=a['depth_days']>0;z=depths[active]
        r=dict(members=m,depths=z,mode='finite ESP contributors' if velocity else 'valid fitted centres only')
        if velocity:
            counts=a['counts'][...,active]
            r.update(u=np.divide(a['ut'][...,active],counts,out=np.full(counts.shape,np.nan),where=counts>0),
                     v=np.divide(a['vt'][...,active],counts,out=np.full(counts.shape,np.nan),where=counts>0),counts=counts,
                     support=pd.DataFrame(dict(Depth=z,eddy_days=a['depth_days'][active],
                        min_cell_days=counts.min(axis=(0,1)),max_cell_days=counts.max(axis=(0,1)))))
        results[key]=r
    return results,pd.DataFrame(day_audit)


def normalised_statistics(members,n_boot=500,seed=731):
    """Normalise each day's vector by its own surface Rc BEFORE averaging."""
    m=members.copy()
    if not (np.isfinite(m.surface_Rc_km)&m.surface_Rc_km.gt(0)).all():raise ValueError('Need finite positive surface Rc')
    m['east_km']=m.east_km/m.surface_Rc_km;m['north_km']=m.north_km/m.surface_Rc_km
    m['distance_km']=np.hypot(m.east_km,m.north_km)
    s,_=pct.summarise(m,n_boot,seed)
    names={c:c.replace('mean_east','mean_east_Rc').replace('mean_north','mean_north_Rc') for c in s}
    names.update(distance_km='distance_Rc',distance_ci_low='distance_Rc_ci_low',distance_ci_high='distance_Rc_ci_high',
                 mean_member_distance_km='mean_member_distance_Rc',var_east='var_east_Rc2',var_north='var_north_Rc2',
                 cov_xy='cov_xy_Rc2',std_east='std_east_Rc',std_north='std_north_Rc')
    return s.rename(columns=names)


def add_statistics(results,n_boot=500,seed=731):
    for key,r in results.items():
        s,_=pct.summarise(r['members'],n_boot,seed)
        bottom=r['members'].groupby('Depth').agg(bottom_checked_days=('below_local_bottom','count'),
            below_bottom_fraction=('below_local_bottom','mean'))
        r['stats']=s.merge(bottom,on='Depth',how='left',validate='one_to_one')
        r['normalised_stats']=normalised_statistics(r['members'],n_boot,seed)
        print('Statistics:',key,flush=True)


def table(results,field='stats'):
    parts=[r[field].assign(**dict(zip(KEYS,key))) for key,r in results.items() if field in r]
    return pd.concat(parts,ignore_index=True) if parts else pd.DataFrame()


def population_context(audit,results):
    """Summaries use actual contributing days, not potentially failed selections."""
    rows=[]
    for key,r in results.items():
        days=r['members'][['Eddy','Day']].drop_duplicates().merge(audit,on=['Eddy','Day'],validate='one_to_one')
        row=dict(zip(KEYS,key));row.update(eddy_days=len(days),eddies=days.Eddy.nunique())
        for col in ['Ro_abs','surface_Rc_km','h','N2_200m_s2','N2_500m_s2']:
            if col in days:
                a=days[col].where(np.isfinite(days[col])).dropna()
                row[col+'_median']=a.median();row[col+'_available_days']=len(a)
        for name in ['planetary','topographic','mixed','unknown']:row['PV_'+name+'_fraction']=days.PV_regime.eq(name).mean()
        for name in ['On-shelf','Off-shelf','Unknown']:row[name+'_fraction']=days.shelf_class.eq(name).mean()
        rows.append(row)
    return pd.DataFrame(rows)


def resolve_depths(available,targets,tolerance=100.):
    """Resolve once on the common dataset depth grid; never per group."""
    available=np.unique(np.asarray(available,float));rows=[];used=set()
    for target in targets:
        nearest=float(available[np.argmin(abs(available-target))]) if len(available) else np.nan
        ok=np.isfinite(nearest) and abs(nearest-target)<=tolerance and nearest not in used
        if ok:used.add(nearest)
        rows.append(dict(target_m=target,Depth=nearest if ok else np.nan,
                         status='ok' if ok else 'no unique fitted level within tolerance'))
    return pd.DataFrame(rows)


def fixed_membership(results,depth_map,shallow_targets,deep_targets,n_boot=500,seed=731):
    """Require every specified exact depth for each day; report empty groups."""
    fixed={};audit=[]
    for key,r in results.items():
        if key[1]=='all':continue
        targets=shallow_targets if key[1]=='shallow' else deep_targets
        mapped=depth_map.set_index('target_m').reindex(targets)
        if mapped.Depth.isna().any():
            audit.append(dict(zip(KEYS,key),matched_days=0,matched_eddies=0,status='unresolved target depth'));continue
        z=np.unique(np.r_[0.,mapped.Depth.to_numpy()])
        m=pct.matched_members(r['members'],z)
        row=dict(zip(KEYS,key));row.update(matched_days=len(m[['Eddy','Day']].drop_duplicates()),matched_eddies=m.Eddy.nunique(),
            required_depths_m=', '.join(f'{x:g}' for x in z),status='ok' if len(m) else 'no days at every required depth')
        audit.append(row)
        if len(m):
            m=m.loc[m.Depth.isin(z)].copy();s,_=pct.summarise(m,n_boot,seed)
            fixed[key]=dict(members=m,stats=s,normalised_stats=normalised_statistics(m,n_boot,seed))
    return fixed,pd.DataFrame(audit)


def summary_at_depth(results,depth,min_eddies=2):
    """Report the same requested actual depth for every group, or explicit missing."""
    rows=[]
    for key,r in results.items():
        row=dict(zip(KEYS,key),Depth=depth,status='depth unavailable',distance_km=np.nan,
                 ci_low=np.nan,ci_high=np.nan,distance_Rc=np.nan,n_eddies=0,n_eddy_days=0)
        s=r['stats'].loc[r['stats'].Depth.eq(depth)]
        if len(s):
            a=s.iloc[0];b=r['normalised_stats'].loc[lambda x:x.Depth.eq(depth)].iloc[0]
            row.update(distance_km=a.distance_km,ci_low=a.distance_ci_low,ci_high=a.distance_ci_high,
                distance_Rc=b.distance_Rc,n_eddies=int(a.n_eddies),n_eddy_days=int(a.n_eddy_days),
                status='ok' if a.n_eddies>=min_eddies else 'sparse')
        rows.append(row)
    return pd.DataFrame(rows)


def plot_regional_profiles(results,cohort='all',ro_class='all',normalised=False,min_eddies=2,sparse_eddies=20):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(13,9),constrained_layout=True)
    field='normalised_stats' if normalised else 'stats'
    col,lo,hi=('distance_Rc','distance_Rc_ci_low','distance_Rc_ci_high') if normalised else ('distance_km','distance_ci_low','distance_ci_high')
    for ax,region in zip(axs.flat,REGIONS):
        found=False
        for cyc in ['AE','CE']:
            r=results.get((region,cohort,ro_class,cyc))
            if r is None:continue
            s=r[field];ok=s.n_eddies.ge(min_eddies);found=True
            ax.plot(s[col].where(ok),s.Depth,color=COLORS[cyc],label=cyc)
            ax.fill_betweenx(s.Depth,s[lo].where(ok),s[hi].where(ok),color=COLORS[cyc],alpha=.15)
            sparse=ok&s.n_eddies.lt(sparse_eddies)
            ax.scatter(s.loc[sparse,col],s.loc[sparse,'Depth'],facecolors='none',edgecolors=COLORS[cyc],s=28)
        ax.set(title=region,xlabel='Mean-vector tilt / surface Rc' if normalised else 'Mean-vector tilt (km)',ylabel='Depth (m)')
        ax.invert_yaxis();ax.grid(alpha=.15)
        if found:ax.legend()
        else:ax.text(.5,.5,'No contributors',ha='center',transform=ax.transAxes)
    limits=[ax.get_xlim() for ax in axs.flat];depths=[ax.get_ylim()[0] for ax in axs.flat]
    for ax in axs.flat:ax.set_xlim(0,max(1.,max(x[1] for x in limits)));ax.set_ylim(max(1.,max(depths)),0)
    fig.suptitle(f'{cohort}, {ro_class} Rossby: mean tilt ±95% CI; open circles <{sparse_eddies} eddies')
    return fig


def plot_group_vectors(results,region,cohort,normalised=False,min_eddies=2):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(12,5),constrained_layout=True)
    field='normalised_stats' if normalised else 'stats'
    cols=['mean_east_Rc','mean_north_Rc','distance_Rc'] if normalised else ['mean_east','mean_north','distance_km']
    for ro,style in [('low','-'),('high','--')]:
        for cyc in ['AE','CE']:
            r=results.get((region,cohort,ro,cyc))
            if r is None:continue
            s=r[field];ok=s.n_eddies.ge(min_eddies)
            for ax,col in zip(axs,cols):
                if col=='distance_km':lo,hi='distance_ci_low','distance_ci_high'
                else:lo,hi=col+'_ci_low',col+'_ci_high'
                ax.plot(s[col].where(ok),s.Depth,style,color=COLORS[cyc],label=f'{cyc} {ro}')
                ax.fill_betweenx(s.Depth,s[lo].where(ok),s[hi].where(ok),color=COLORS[cyc],alpha=.10)
    for ax,label in zip(axs,['East displacement','North displacement','Mean-vector distance']):
        ax.set(xlabel=label+(' / surface Rc' if normalised else ' (km)'),ylabel='Depth (m)');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.5)
    if axs[0].lines:axs[0].legend()
    fig.suptitle(f'{region}, {cohort}: low/high surface |Ro|')
    return fig


def comparison_sets(results,regions=REGIONS):
    """Adapter for the existing full-ESP fitting and section plotting routines."""
    return {f'{region} | {ro} Ro':{(cohort,cyc):r for (rg,cohort,rr,cyc),r in results.items()
            if rg==region and rr==ro and cohort!='all'} for region in regions for ro in ['low','high']}


def horizontal_maps(results,X,Y,region,cohort,ro_class,depth):
    import matplotlib.pyplot as plt
    fields=[]
    for cyc in ['AE','CE']:
        r=results.get((region,cohort,ro_class,cyc))
        if r is not None and 'u' in r and depth in r['depths']:
            k=int(np.flatnonzero(r['depths']==depth)[0]);f=np.hypot(r['u'][...,k],r['v'][...,k])
            if np.isfinite(f).any():fields.append((cyc,r,k,f))
    if not fields:return None
    vmax=max(np.nanmax(f) for _,_,_,f in fields)
    fig,axs=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
    for ax,cyc in zip(axs,['AE','CE']):
        item=next((x for x in fields if x[0]==cyc),None)
        ax.set(title=cyc,xlabel='Model-grid x (km)',ylabel='Model-grid y (km)',aspect='equal')
        if item is None:ax.text(.5,.5,'No field at this depth',ha='center',transform=ax.transAxes);continue
        _,r,k,f=item;im=ax.pcolormesh(X,Y,f,shading='auto',cmap='magma',vmin=0,vmax=vmax)
        step=max(1,X.shape[0]//11);sl=(slice(None,None,step),slice(None,None,step))
        ax.quiver(X[sl],Y[sl],r['u'][...,k][sl],r['v'][...,k][sl],color='white',scale=max(vmax,1e-12)*15,scale_units='width')
        m=r['members'].loc[lambda x:x.Depth.eq(depth)]
        ax.plot(m.xc.mean(),m.yc.mean(),'co',label='Mean centre');ax.plot(0,0,'w+');ax.legend()
    fig.colorbar(im,ax=list(axs),label='Speed of mean ESP velocity (m/s)')
    fig.suptitle(f'{region}, {cohort}, {ro_class} Ro, {depth:g} m')
    return fig


REGION_GROUPS = ('S','U','D')
REGION_GROUP_MAP = {r:r[0] for r in REGIONS}


def combine_region_labels(surface):
    """Pool raw days before statistics; keep original labels for diagnostics."""
    s=surface.copy()
    s['Subregion']=s.Region
    s['Region']=s.Subregion.map(REGION_GROUP_MAP)
    return s


def plot_combined_regions(results,split_rossby=False,normalised=False,min_eddies=2,sparse_eddies=20, figsize=(12,8), title=True, sharex=False):
    """2x3: shallow/deep rows and S/U/D columns, equal-day pooled estimates."""
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=figsize,constrained_layout=True, sharex=sharex)
    field='normalised_stats' if normalised else 'stats'
    col,lo,hi=('distance_Rc','distance_Rc_ci_low','distance_Rc_ci_high') if normalised else ('distance_km','distance_ci_low','distance_ci_high')
    styles=[('low','-'),('high','--')] if split_rossby else [('all','-')]
    for i,cohort in enumerate(['shallow','deep']):
        deepest=1.
        for j,region in enumerate(REGION_GROUPS):
            ax=axs[i,j];found=False
            for ro,style in styles:
                for cyc in ['AE','CE']:
                    r=results.get((region,cohort,ro,cyc))
                    if r is None:continue
                    s=r[field];ok=s.n_eddies.ge(min_eddies);found=True;deepest=max(deepest,float(s.Depth.max()))
                    ax.plot(s[col].where(ok),s.Depth,style,color=COLORS[cyc],label=f'{cyc} {ro}' if split_rossby else cyc)
                    ax.fill_betweenx(s.Depth,s[lo].where(ok),s[hi].where(ok),color=COLORS[cyc],alpha=.12)
                    sparse=ok&s.n_eddies.lt(sparse_eddies)
                    ax.scatter(s.loc[sparse,col],s.loc[sparse,'Depth'],s=20,facecolors='none',edgecolors=COLORS[cyc])
            if title: ax.set_title(f'{region} — {cohort}')
            if sharex and (i==1): 
                ax.set_xlabel='Tilt / surface Rc' if normalised else 'Tilt distance (km)'
            if not sharex:
                ax.set_xlabel='Tilt / surface Rc' if normalised else 'Tilt distance (km)'
            ax.set(ylabel='Depth (m)' if j==0 else '')
            if found:ax.legend(fontsize=8)
            else:ax.text(.5,.5,'No contributors',ha='center',transform=ax.transAxes)
        for ax in axs[i]:ax.set_ylim(deepest,0)
    xmax=max(max(ax.get_xlim()) for ax in axs.flat)
    for ax in axs.flat:ax.set_xlim(0,max(xmax,1e-6));ax.grid(alpha=.15)
    fig.suptitle(('Low |Ro| < 0.5 solid; high |Ro| ≥ 0.5 dashed' if split_rossby else 'All Rossby numbers pooled')+
                 f'; 95% CI; open circles <{sparse_eddies} eddies')
    return fig,axs


def plot_combined_region_components(results, ro_class='all', normalised=False,
                                    min_eddies=2, sparse_eddies=20):
    """Signed geographic components: shallow/deep rows, S/U/D columns.

    Each panel overlays AE/CE (colour) and east/north (solid/dashed).
    Uses existing constituent-centre means and whole-track bootstrap intervals;
    normalised statistics divide each constituent by its own surface Rc first.
    """
    import matplotlib.pyplot as plt
    if ro_class not in ('all', 'low', 'high'):
        raise ValueError("ro_class must be 'all', 'low', or 'high'")
    fig, axs = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    field = 'normalised_stats' if normalised else 'stats'
    suffix = '_Rc' if normalised else ''
    for i, cohort in enumerate(['shallow', 'deep']):
        deepest = 1.
        for j, region in enumerate(REGION_GROUPS):
            ax = axs[i, j]
            found = False
            for cyc in ['AE', 'CE']:
                result = results.get((region, cohort, ro_class, cyc))
                if result is None or result[field].empty:
                    continue
                s = result[field].sort_values('Depth')
                ok = s.n_eddies.ge(min_eddies)
                deepest = max(deepest, float(s.Depth.max()))
                found |= bool(ok.any())
                for component, label, style in [('east', 'zonal', '-'), ('north', 'meridional', '--')]:
                    col = f'mean_{component}{suffix}'
                    ax.plot(s[col].where(ok), s.Depth, style, color=COLORS[cyc], label=f'{cyc} {label}')
                    ax.fill_betweenx(s.Depth, s[col+'_ci_low'].where(ok),
                                     s[col+'_ci_high'].where(ok), color=COLORS[cyc], alpha=.10)
                    sparse = ok & s.n_eddies.lt(sparse_eddies)
                    ax.scatter(s.loc[sparse, col], s.loc[sparse, 'Depth'], s=18,
                               facecolors='none', edgecolors=COLORS[cyc])
            ax.axvline(0, color='.4', lw=.7)
            ax.grid(alpha=.15)
            ax.set(title=f'{region} — {cohort}',
                   xlabel='Signed displacement / surface Rc' if normalised else 'Signed displacement (km)',
                   ylabel='Depth (m)' if j == 0 else '')
            if found:
                ax.legend(fontsize=8)
            else:
                ax.text(.5, .5, 'Insufficient contributors', ha='center', transform=ax.transAxes)
        for ax in axs[i]:
            ax.set_ylim(deepest, 0)
    limit = max(1e-6, max(abs(v) for ax in axs.flat for v in ax.get_xlim()))
    for ax in axs.flat:
        ax.set_xlim(-limit, limit)
    title = 'All Rossby classes pooled' if ro_class == 'all' else f'{ro_class.capitalize()} Rossby class'
    fig.suptitle(title + '\nZonal solid (+east); meridional dashed (+north); 95% CI; '
                 + f'open circles <{sparse_eddies} eddies')
    return fig, axs
