"""Constituent-centre displacement in each day's fixed environmental PV frame.

Equal-day estimates, whole-track bootstrap, exact fitted depths, no ESP velocity
reconstruction or interpolation. Bearings are clockwise from geographic north.
"""
from itertools import product
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import planetary_composite_tools as pct

KEYS = ['Eddy', 'Day']
GROUPS = ['family', 'panel', 'extent_group', 'split', 'Cyc']
COLORS = {'AE': 'firebrick', 'CE': 'royalblue'}
PANELS = {'overall': ['All'], 'regional': ['S', 'U', 'D'],
          'dominance': ['planetary', 'mixed', 'topographic'], 'spread': ['All'],
          'cancellation': ['All'], 'matched': ['All']}


def prepare_members(surface, vertical, grid_angle, split_depth=1000.,
                    spread_depth=1000., ro_split=.5, dominance_factor=2.,
                    cancellation_min=.1):
    """Return every input day's audit, eligible centres, and pooled median thresholds.

    Spread requires >=2 valid fitted levels in the common upper-column window;
    it is NOT extrapolated to 1000 m. Coverage is reported for comparison.
    Thresholds use eligible days, pool polarity/region/Ro/dominance, and are
    computed separately in shallow/deep cohorts. High includes median ties.
    """
    if not np.isscalar(grid_angle) or not np.isfinite(grid_angle):
        raise ValueError('Expected one finite model-grid angle in radians')
    if split_depth <= 0 or spread_depth <= 0 or dominance_factor <= 1:
        raise ValueError('Depth bounds must be positive; dominance_factor must exceed one')
    if not 0 <= cancellation_min <= 1:
        raise ValueError('cancellation_min must lie between zero and one')
    required = KEYS + ['Cyc', 'Region', 'Ro', 'PV_grad_mag', 'PV_grad_theta',
                       'PV_grad_plan_mag', 'PV_grad_topo_mag']
    missing = set(required) - set(surface)
    if missing:
        raise ValueError(f'Missing surface columns: {sorted(missing)}')
    if surface[KEYS].isna().any().any():
        raise ValueError('Missing surface Eddy/Day keys')
    audit, profiles = pct.prepare_profiles(surface, vertical, split_depth)
    audit['Region_group'] = audit.Region.map({r: r[0] for r in ['S','U','D','S1','S2','U1','U2','D1','D2']})
    audit['Ro_class'] = np.where(np.isfinite(audit.Ro), np.where(audit.Ro.abs() < ro_split, 'low', 'high'), 'unknown')
    plan, topo = audit.PV_grad_plan_mag, audit.PV_grad_topo_mag
    finite = np.isfinite(plan) & np.isfinite(topo) & plan.ge(0) & topo.ge(0) & (plan + topo).gt(0)
    audit['PV_regime'] = np.select([finite & (plan >= dominance_factor*topo),
        finite & (topo >= dominance_factor*plan), finite], ['planetary','topographic','mixed'], default='unknown')
    audit['cancellation_ratio'] = audit.PV_grad_mag.div((plan+topo).where(finite))
    audit['cancellation_ok'] = audit.cancellation_ratio.ge(cancellation_min) & np.isfinite(audit.cancellation_ratio)
    audit['frame_valid'] = np.isfinite(audit.PV_grad_mag) & audit.PV_grad_mag.gt(0) & np.isfinite(audit.PV_grad_theta)
    audit['analysis_status'] = np.select([
        ~audit.Cyc.isin(['AE','CE']), audit.profile_status.ne('usable'), ~audit.frame_valid],
        ['unrecognised polarity','no valid surface fit','invalid PV reference'], default='included')
    # Missing regions/Ro/dominance only exclude their particular subgroup plots.
    eligible = audit.loc[audit.analysis_status.eq('included')]
    profiles = profiles.merge(eligible[KEYS], on=KEYS, validate='many_to_one')
    radius = profiles.loc[profiles.Depth.eq(0), KEYS+['Rc']].rename(columns={'Rc':'surface_Rc_km'})
    spread_rows = []
    for (eddy, day), p in profiles.loc[profiles.Depth.le(spread_depth)].groupby(KEYS, sort=False):
        spread_rows.append(dict(Eddy=eddy, Day=day, spread_levels=len(p),
            spread_max_depth_m=p.Depth.max(),
            centre_spread_km=float(pdist(p[['xc','yc']]).max()) if len(p)>1 else np.nan))
    spreads = pd.DataFrame(spread_rows, columns=KEYS+['spread_levels','spread_max_depth_m','centre_spread_km'])
    audit = audit.merge(spreads, on=KEYS, how='left', validate='one_to_one').merge(radius, on=KEYS, how='left', validate='one_to_one')
    thresholds = (audit.loc[audit.analysis_status.eq('included')].groupby('extent_group')
        .agg(spread_threshold_km=('centre_spread_km','median'),
             classified_days=('centre_spread_km','count')).reindex(['shallow','deep']).reset_index())
    audit = audit.merge(thresholds[['extent_group','spread_threshold_km']], on='extent_group', how='left', validate='many_to_one')
    audit['spread_class'] = np.where(audit.centre_spread_km.notna(),
        np.where(audit.centre_spread_km < audit.spread_threshold_km, 'low','high'), 'unknown')
    fields = KEYS+['Cyc','Region_group','extent_group','Ro_class','PV_regime','spread_class',
                   'PV_grad_mag','PV_grad_theta','cancellation_ratio','cancellation_ok','surface_Rc_km']
    members = profiles.merge(audit[fields], on=KEYS, validate='many_to_one')
    # Gradient's mathematical angle measured from model-grid +x.
    angle = np.pi/2 - np.deg2rad(members.PV_grad_theta) - grid_angle
    members['along_km'] = np.cos(angle)*members.xc + np.sin(angle)*members.yc
    members['perp_km'] = -np.sin(angle)*members.xc + np.cos(angle)*members.yc
    members['along_Rc'] = members.along_km / members.surface_Rc_km
    members['perp_Rc'] = members.perp_km / members.surface_Rc_km
    members['distance_km'] = np.hypot(members.along_km, members.perp_km)
    return audit, members, thresholds


def group_specs(families=('overall','regional','dominance','spread','cancellation')):
    """A small set of comparisons, never a full Cartesian product of categories."""
    for family in families:
        splits = ['low','high'] if family in ('regional','spread') else (['all','retained'] if family=='cancellation' else ['all'])
        for panel, extent, split, cyc in product(PANELS[family], ['shallow','deep'], splits, ['AE','CE']):
            yield family, panel, extent, split, cyc


def select_group(members, key):
    family, panel, extent, split, cyc = key
    keep = members.extent_group.eq(extent) & members.Cyc.eq(cyc)
    if family == 'regional': keep &= members.Region_group.eq(panel) & members.Ro_class.eq(split)
    if family == 'dominance': keep &= members.PV_regime.eq(panel)
    if family == 'spread': keep &= members.spread_class.eq(split)
    if family == 'cancellation' and split == 'retained': keep &= members.cancellation_ok
    return members.loc[keep]


def summarise(members, n_boot=500, seed=731, direction_epsilon_km=1e-8):
    """Pointwise 95% CIs using the SAME whole-track draws across depths/metrics.

    Day-weighted means; normalisation occurs before averaging. With/against
    fractions exclude effectively zero parallel components; their denominator
    and the number of neutral vectors are reported separately. CIs require
    >=2 contributing eddies (directional CIs >=2 directional eddies).
    """
    if n_boot < 2: raise ValueError('n_boot must be at least two')
    if members.duplicated(KEYS+['Depth']).any(): raise ValueError('Duplicate centre samples')
    if members.empty: return pd.DataFrame()
    ids = np.sort(members.Eddy.unique())
    rng = np.random.default_rng(seed)
    # Process bootstrap batches: memory independent of total track count x draws.
    levels = []
    metrics = ['along_km','perp_km','along_Rc','perp_Rc']
    for depth, d in members.groupby('Depth', sort=True):
        signed = d.along_km.abs().gt(direction_epsilon_km)
        work = d[metrics+['Eddy']].copy()
        work['count'] = 1.
        work['directional'] = signed.astype(float)
        work['with'] = (d.along_km > direction_epsilon_km).astype(float)
        sums = work.groupby('Eddy').sum().reindex(ids, fill_value=0).to_numpy(float)
        levels.append((depth, d, signed, sums))
    estimates = [[] for _ in levels]
    for start in range(0, n_boot, 32):
        draws = rng.multinomial(len(ids), np.full(len(ids), 1/len(ids)), size=min(32,n_boot-start)).astype(float)
        for i, (_, _, _, sums) in enumerate(levels):
            totals = draws @ sums
            denom = totals[:,4:5]
            means = np.divide(totals[:,:4], denom, out=np.full((len(draws),4),np.nan), where=denom>0)
            fraction = np.divide(totals[:,6], totals[:,5], out=np.full(len(draws),np.nan), where=totals[:,5]>0)
            estimates[i].append(np.column_stack([means,fraction]))
    rows = []
    for (depth,d,signed,_), batches in zip(levels, estimates):
        boot = np.concatenate(batches)
        row = dict(Depth=depth,n_eddy_days=len(d),n_eddies=d.Eddy.nunique(),
                   directional_days=int(signed.sum()),directional_eddies=d.loc[signed,'Eddy'].nunique(),
                   neutral_days=int((~signed).sum()),mean_member_distance_km=d.distance_km.mean())
        for j, metric in enumerate(metrics+['fraction_with']):
            value = d[metric].mean() if j<4 else ((d.loc[signed,'along_km']>0).mean() if signed.any() else np.nan)
            count = row['n_eddies'] if j<4 else row['directional_eddies']
            finite = np.isfinite(boot[:,j])
            lo,hi = np.quantile(boot[finite,j],[.025,.975]) if count>=2 and finite.any() else (np.nan,np.nan)
            row.update({metric:value, metric+'_ci_low':lo, metric+'_ci_high':hi})
            if j<4: row[metric+'_variance'] = d[metric].var(ddof=1)
        row['fraction_against'] = 1-row['fraction_with']
        row['fraction_against_ci_low'] = 1-row['fraction_with_ci_high']
        row['fraction_against_ci_high'] = 1-row['fraction_with_ci_low']
        rows.append(row)
    return pd.DataFrame(rows)


def analyse(members, n_boot=500, seed=731, families=('overall','regional','dominance','spread','cancellation')):
    stats, inventory = [], []
    for key in group_specs(families):
        d = select_group(members, key)
        labels = dict(zip(GROUPS,key))
        inventory.append(dict(**labels,n_eddy_days=len(d[KEYS].drop_duplicates()),n_eddies=d.Eddy.nunique()))
        if len(d): stats.append(summarise(d,n_boot,seed).assign(**labels))
    return (pd.concat(stats,ignore_index=True) if stats else pd.DataFrame(), pd.DataFrame(inventory))


def matched_profiles(members, targets, tolerance=100.):
    """Match targets once on the global exact-depth grid; same days across depths.

    targets maps shallow/deep to target depths. Unresolved or duplicate mappings
    produce an explicit empty cohort, rather than silently dropping a target.
    """
    available = np.sort(members.Depth.unique())
    maps, parts = [], []
    for extent, requested in targets.items():
        resolved=[];ok=True
        for target in sorted(set([0.,*requested])):
            z = available[np.argmin(abs(available-target))] if len(available) else np.nan
            valid = np.isfinite(z) and abs(z-target)<=tolerance and z not in resolved
            maps.append(dict(extent_group=extent,target_m=target,Depth=z if valid else np.nan,status='ok' if valid else 'unresolved/duplicate fitted depth'))
            ok &= valid
            resolved.append(z)
        if ok:
            d=members.loc[members.extent_group.eq(extent)]
            matched=pct.matched_members(d,resolved)
            # Retain full profiles of these days for the sensitivity curves.
            parts.append(matched)
    return (pd.concat(parts,ignore_index=True) if parts else members.iloc[:0].copy(), pd.DataFrame(maps))


def plot_profiles(stats, family='overall', unit='km', component=None, min_eddies=2, sparse_eddies=20):
    """Two rows (shallow/deep); component columns or S/U/D and dominance columns."""
    import matplotlib.pyplot as plt
    if family in ('regional','dominance'):
        if component not in ('along','perp'): raise ValueError('Specify along or perp')
        columns=[(p,component) for p in PANELS[family]]
    else: columns=[('All','along'),('All','perp')]
    fig,axes=plt.subplots(2,len(columns),figsize=(4.2*len(columns),8),squeeze=False,constrained_layout=True)
    splits=['low','high'] if family in ('regional','spread') else (['all','retained'] if family=='cancellation' else ['all'])
    for i,extent in enumerate(['shallow','deep']):
        for j,(panel,comp) in enumerate(columns):
            ax=axes[i,j];metric=f'{comp}_{unit}';found=False
            for split,cyc in product(splits,['AE','CE']):
                if stats.empty: continue
                d=stats.loc[stats.family.eq(family)&stats.panel.eq(panel)&stats.extent_group.eq(extent)&stats.split.eq(split)&stats.Cyc.eq(cyc)].sort_values('Depth')
                if d.empty: continue
                valid=d.n_eddies.ge(min_eddies);found |= valid.any()
                label=cyc if split=='all' and len(splits)==1 else f'{cyc} {split}'
                style='--' if split in ('high','retained') else '-'
                ax.plot(d[metric].where(valid),d.Depth,style,color=COLORS[cyc],label=label)
                ax.fill_betweenx(d.Depth,d[metric+'_ci_low'].where(valid),d[metric+'_ci_high'].where(valid),color=COLORS[cyc],alpha=.12)
                sparse=valid & d.n_eddies.lt(sparse_eddies)
                ax.scatter(d.loc[sparse,metric],d.loc[sparse,'Depth'],s=14,facecolors='none',edgecolors=COLORS[cyc])
            ax.axvline(0,color='.4',lw=.7);ax.grid(alpha=.15)
            ax.set(title=f'{panel} | {extent}',xlabel=f'{comp.capitalize()} displacement '+('(km)' if unit=='km' else '/ surface Rc'),ylabel='Depth (m)')
            if found: ax.legend(fontsize=8)
            else: ax.text(.5,.5,'Insufficient contributors',transform=ax.transAxes,ha='center')
        deepest=max([max(ax.get_ylim()) for ax in axes[i]]+[1.])
        for ax in axes[i]: ax.set_ylim(deepest,0)
    limit=max([max(abs(x) for x in ax.get_xlim()) for ax in axes.flat]+[1e-3])
    for ax in axes.flat: ax.set_xlim(-limit,limit)
    fig.suptitle(f'{family.capitalize()}: signed displacement; pointwise 95% track-bootstrap CI\nOpen circles: fewer than {sparse_eddies} eddies')
    return fig,axes


def direction_at_depths(stats, targets=(200.,500.,1000.,1500.), tolerance=100.):
    """One common exact fitted depth per target for all baseline groups."""
    if stats.empty: return pd.DataFrame(),pd.DataFrame()
    s=stats.loc[stats.family.eq('overall')]
    available=np.sort(s.Depth.unique());rows=[];chosen=[]
    for target in targets:
        z=available[np.argmin(abs(available-target))] if len(available) else np.nan
        valid=np.isfinite(z) and z>0 and abs(z-target)<=tolerance and z not in chosen
        rows.append(dict(target_m=target,Depth=z if valid else np.nan,status='ok' if valid else 'unresolved/duplicate fitted depth'))
        if valid: chosen.append(z)
    result=s.loc[s.Depth.isin(chosen)].copy()
    return result,pd.DataFrame(rows)


def plot_direction(table, min_eddies=2):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    depths=np.sort(table.Depth.unique()) if len(table) else []
    for ax,extent in zip(axes,['shallow','deep']):
        for cyc,offset in [('AE',-.10),('CE',.10)]:
            d=table.loc[table.extent_group.eq(extent)&table.Cyc.eq(cyc)].set_index('Depth').reindex(depths) if len(table) else pd.DataFrame()
            if d.empty: continue
            ok=d.directional_eddies.ge(min_eddies)
            x=np.arange(len(depths))+offset
            ax.plot(x,d.fraction_with.where(ok),'o',color=COLORS[cyc],label=cyc)
            ax.vlines(x[ok],d.loc[ok,'fraction_with_ci_low'],d.loc[ok,'fraction_with_ci_high'],color=COLORS[cyc])
            for k,(_,r) in enumerate(d.iterrows()):
                if ok.iloc[k]: ax.annotate(f'{int(r.directional_eddies)}e/{int(r.directional_days)}d',(x[k],r.fraction_with),xytext=(-4 if cyc=='AE' else 4,-14 if r.fraction_with>.8 else 8),textcoords='offset points',ha='right' if cyc=='AE' else 'left',fontsize=7)
        ax.axhline(.5,color='.4',ls=':');ax.set(xlim=(-.5,max(len(depths)-.5,.5)),ylim=(0,1),title=extent,xticks=np.arange(len(depths)),xticklabels=[f'{z:g}' for z in depths],xlabel='Actual fitted depth (m)',ylabel='Fraction displaced with PV gradient')
        if ax.lines and len(table): ax.legend()
    fig.suptitle('With versus against: neutral parallel components excluded; 95% track-bootstrap CI')
    return fig,axes
