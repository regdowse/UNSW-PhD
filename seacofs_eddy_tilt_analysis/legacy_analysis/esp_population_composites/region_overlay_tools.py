"""Region/polarity overlays with explicit data caps and surface-table radii."""
import numpy as np
import pandas as pd
import planetary_composite_tools as pct

COLORS = {'S': 'darkorange', 'U': 'seagreen', 'D': 'slateblue'}
KEYS = ['Eddy', 'Day']


def regional_members(results):
    """Read each contributing centre once from all-day, all-Rossby baselines."""
    parts = [r['members'] for (region, cohort, ro, cyc), r in results.items()
             if cohort == 'all' and ro == 'all']
    if not parts:
        raise ValueError('No all-day regional constituent centres')
    members = pd.concat(parts, ignore_index=True)
    if members.duplicated(KEYS + ['Depth']).any():
        raise ValueError('Duplicate constituent centres in regional baselines')
    return members


def compare_regions(members, audit, surface, frame='geographic', n_boot=500, seed=731):
    """Classify on full-profile audit, then cap all/shallow DATA at 1000 m.

    Deep keeps all contributing depths. All Rossby classes are pooled at the
    day level. Normalised means use surface.Rc (km), never profile-level Rc.
    Invalid/missing radii exclude only normalised estimates, with counts exposed.
    """
    if frame not in ('geographic', 'pv'):
        raise ValueError('frame must be geographic or pv')
    if members.duplicated(KEYS + ['Depth']).any():
        raise ValueError('Duplicate constituent centres')
    if audit[KEYS].duplicated().any() or surface[KEYS].duplicated().any():
        raise ValueError('Duplicate eddy-day metadata')
    region = 'Region_group' if 'Region_group' in audit else 'Region'
    metadata = audit[KEYS + [region, 'extent_group', 'Cyc']].rename(columns={region:'region'})
    metadata['region'] = metadata.region.map({r:r[0] for r in ['S','U','D','S1','S2','U1','U2','D1','D2']})
    x, y = ('east_km', 'north_km') if frame == 'geographic' else ('along_km','perp_km')
    m = members[KEYS + ['Depth', x, y]].rename(columns={x:'east_km', y:'north_km'})
    m = m.merge(metadata, on=KEYS, how='left', validate='many_to_one', indicator=True)
    if m['_merge'].ne('both').any():
        raise ValueError('Every constituent needs full-profile classification')
    m = m.drop(columns='_merge').merge(surface[KEYS+['Rc']].rename(columns={'Rc':'radius_km'}),
                                      on=KEYS, how='left', validate='many_to_one')
    m = m.loc[np.isfinite(m[['Depth','east_km','north_km']]).all(axis=1) & m.Depth.ge(0)].copy()
    if not m.extent_group.isin(['shallow','deep']).all():
        raise ValueError('Classify full profiles as shallow/deep before applying caps')
    rows, inventory = [], []
    for cohort in ['all','shallow','deep']:
        d = m if cohort == 'all' else m.loc[m.extent_group.eq(cohort)]
        if cohort != 'deep':
            d = d.loc[d.Depth.le(1000)].copy()  # filter before bootstrap/statistics
        for region in COLORS:
            for cyc in ['AE','CE']:
                group = d.loc[d.region.eq(region) & d.Cyc.eq(cyc)]
                for unit in ['km','Rc']:
                    selected = group.copy()
                    if unit == 'Rc':
                        selected = selected.loc[np.isfinite(selected.radius_km) & selected.radius_km.gt(0)].copy()
                        selected[['east_km','north_km']] = selected[['east_km','north_km']].div(selected.radius_km, axis=0)
                    inventory.append(dict(cohort=cohort,region=region,Cyc=cyc,unit=unit,
                        n_eddy_days=len(selected[KEYS].drop_duplicates()),n_eddies=selected.Eddy.nunique(),
                        max_depth_m=selected.Depth.max(),
                        excluded_radius_days=len(group[KEYS].drop_duplicates())-len(selected[KEYS].drop_duplicates())))
                    if selected.empty:
                        continue
                    selected['distance_km'] = np.hypot(selected.east_km, selected.north_km)
                    stats, _ = pct.summarise(selected, n_boot=n_boot, seed=seed)
                    rows.append(stats.assign(cohort=cohort,region=region,Cyc=cyc,unit=unit))
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(), pd.DataFrame(inventory))


def plot_comparison(stats, cohort='all', frame='geographic', min_eddies=2, sparse_eddies=20, figsize=None):
    """Rows: km / surface Rc; columns: geographic components+norm, or PV components."""
    import matplotlib.pyplot as plt
    if cohort not in ('all','shallow','deep') or frame not in ('geographic','pv'):
        raise ValueError('Invalid cohort or frame')
    labels = ['Zonal (+east)', 'Meridional (+north)', 'Tilt distance'] if frame == 'geographic' else ['Along PV gradient', 'Perpendicular (+CCW)']
    metrics = ['mean_east','mean_north','distance_km'][:len(labels)]
    if figsize is None:
        figsize=(4.5*len(labels),8)
    fig, axes = plt.subplots(2,len(labels),figsize=figsize,squeeze=False,constrained_layout=True)
    deepest = 1000. if cohort != 'deep' else 1.
    for i, unit in enumerate(['km','Rc']):
        for j, (metric,label) in enumerate(zip(metrics,labels)):
            ax=axes[i,j];found=False
            for region,color in COLORS.items():
                for cyc,style in [('AE','-'),('CE','--')]:
                    if stats.empty: continue
                    s=stats.loc[stats.cohort.eq(cohort)&stats.region.eq(region)&stats.Cyc.eq(cyc)&stats.unit.eq(unit)].sort_values('Depth')
                    if s.empty: continue
                    deepest=max(deepest,float(s.Depth.max()))
                    valid=s.n_eddies.ge(min_eddies);found |= valid.any()
                    lo,hi=('distance_ci_low','distance_ci_high') if metric=='distance_km' else (metric+'_ci_low',metric+'_ci_high')
                    ax.plot(s[metric].where(valid),s.Depth,style,color=color,label=f'{region} {cyc}')
                    ax.fill_betweenx(s.Depth,s[lo].where(valid),s[hi].where(valid),color=color,alpha=.09)
                    sparse=valid & s.n_eddies.lt(sparse_eddies)
                    ax.scatter(s.loc[sparse,metric],s.loc[sparse,'Depth'],s=14,facecolors='none',edgecolors=color)
            ax.set(xlabel=label+(' (km)' if unit=='km' else ' / surface Rc'),ylabel='Depth (m)' if j==0 else '')
            ax.axvline(0,color='.5',lw=.6);ax.grid(alpha=.15)
            if found:ax.legend(fontsize=8,ncol=2)
            else:ax.text(.5,.5,'Insufficient contributors',ha='center',transform=ax.transAxes)
        # Common signed scale for the two components within each units row.
        limit=max(1e-6,max(abs(v) for ax in axes[i,:2] for v in ax.get_xlim()))
        for ax in axes[i,:2]:ax.set_xlim(-limit,limit)
        if len(labels)==3:axes[i,2].set_xlim(left=0)
    for ax in axes.flat:ax.set_ylim(deepest,0)
    extent='data ≤1000 m' if cohort!='deep' else 'full contributing depth range'
    fig.suptitle(f'{cohort.capitalize()} eddy-days: {extent}; all Rossby classes\nAE solid; CE dashed; 95% track-bootstrap CI; open circles <{sparse_eddies} eddies')
    return fig,axes
