"""One upper-ocean PV-relative centreline comparison with explicit day selection."""
import numpy as np
import pandas as pd
import planetary_composite_tools as pct
import topographic_composite_tools as tct

KEYS = ['Eddy', 'Day']
REGIMES = ['planetary', 'mixed', 'topographic']
COLORS = {'AE': 'firebrick', 'CE': 'royalblue'}


def prepare(surface, vertical, grid_angle, min_tilt_km=5., dominance_factor=2.,
            water_depth_m=3000., max_profile_depth_m=1000.):
    """Filter whole days by surface.TiltDis, then composite only capped profiles.

    Mixed is the complement of BOTH complete magnitude+bathymetry selections
    among classifiable days. Missing/invalid gradients or water depths are
    excluded, not silently labelled mixed. Exactly boundary water depth is mixed.
    """
    required = KEYS + ['Cyc', 'TiltDis', 'h', 'PV_grad_plan_mag',
                       'PV_grad_topo_mag', 'PV_grad_mag', 'PV_grad_theta']
    missing = set(required) - set(surface)
    if missing:
        raise ValueError(f'Missing surface columns: {sorted(missing)}')
    if surface[KEYS].isna().any().any() or surface[KEYS].duplicated().any():
        raise ValueError('Surface requires unique nonmissing Eddy-Day keys')
    if not np.isscalar(grid_angle) or not np.isfinite(grid_angle):
        raise ValueError('Need a finite scalar model-grid angle in radians')
    if (not np.isfinite([min_tilt_km, dominance_factor, water_depth_m, max_profile_depth_m]).all()
            or min_tilt_km < 0 or dominance_factor <= 1 or min(water_depth_m, max_profile_depth_m) <= 0):
        raise ValueError('Invalid selection thresholds')
    audit = surface.copy()
    plan, topo = audit.PV_grad_plan_mag, audit.PV_grad_topo_mag
    valid = (np.isfinite(plan) & np.isfinite(topo) & plan.ge(0) & topo.ge(0)
             & (plan + topo).gt(0) & np.isfinite(audit.h) & audit.h.gt(0))
    planetary = valid & plan.ge(dominance_factor * topo) & audit.h.gt(water_depth_m)
    topographic = valid & topo.ge(dominance_factor * plan) & audit.h.lt(water_depth_m)
    audit['regime'] = np.select([planetary, topographic, valid],
                                ['planetary', 'topographic', 'mixed'], default='unclassified')
    frame_valid = np.isfinite(audit.PV_grad_mag) & audit.PV_grad_mag.gt(0) & np.isfinite(audit.PV_grad_theta)
    audit['selection_status'] = np.select([
        ~audit.Cyc.isin(['AE','CE']), ~np.isfinite(audit.TiltDis), audit.TiltDis.lt(min_tilt_km),
        ~valid, ~frame_valid], ['unrecognised polarity', 'missing/invalid surface TiltDis',
        'surface TiltDis below threshold', 'invalid dominance inputs or water depth',
        'invalid total PV reference'], default='selected')
    audit['frame_angle_rad'] = np.pi/2 - np.deg2rad(audit.PV_grad_theta) - grid_angle
    sample = audit.loc[audit.selection_status.eq('selected')]
    # Actual data cap: deep levels do not enter profile validity, means or CIs.
    capped = vertical.loc[vertical.Depth.between(0, max_profile_depth_m)].copy()
    profile_audit, centred = pct.prepare_profiles(sample, capped)
    audit = audit.merge(profile_audit[KEYS + ['profile_status']], on=KEYS,
                        how='left', validate='one_to_one')
    audit.loc[audit.selection_status.eq('selected') & audit.profile_status.ne('usable'),
              'selection_status'] = 'no valid surface fit in capped profile'
    members = tct.rotate_members(centred, sample)
    members = members.merge(sample[KEYS + ['Cyc', 'regime']], on=KEYS, validate='many_to_one')
    return audit, members


def summarise(members, n_boot=500, seed=731):
    """Equal-day centre means and pointwise whole-track bootstrap confidence bands."""
    parts, counts = [], []
    for regime in REGIMES:
        for cyc in ['AE', 'CE']:
            m = members.loc[members.regime.eq(regime) & members.Cyc.eq(cyc)]
            counts.append(dict(regime=regime, Cyc=cyc,
                eddy_days=len(m[KEYS].drop_duplicates()), eddies=m.Eddy.nunique(),
                deepest_contributing_m=m.Depth.max()))
            if not m.empty:
                s, _ = tct.summarise(m, n_boot=n_boot, seed=seed)
                parts.append(s.assign(regime=regime, Cyc=cyc))
    return (pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(), pd.DataFrame(counts))


def plot(stats, max_depth_m=1000., min_eddies=2, sparse_eddies=20):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharey=True, constrained_layout=True)
    for i, (metric, label) in enumerate([('mean_along', 'Along PV gradient'),
                                       ('mean_perp', 'Perpendicular (+CCW)')]):
        for j, regime in enumerate(REGIMES):
            ax = axes[i, j]; found = False
            for cyc, color in COLORS.items():
                if stats.empty: continue
                s = stats.loc[stats.regime.eq(regime) & stats.Cyc.eq(cyc)].sort_values('Depth')
                if s.empty: continue
                ok = s.n_eddies.ge(min_eddies); found |= bool(ok.any())
                ax.plot(s[metric].where(ok), s.Depth, color=color, label=cyc)
                ax.fill_betweenx(s.Depth, s[metric+'_ci_low'].where(ok),
                                 s[metric+'_ci_high'].where(ok), color=color, alpha=.15)
                sparse = ok & s.n_eddies.lt(sparse_eddies)
                ax.scatter(s.loc[sparse, metric], s.loc[sparse, 'Depth'],
                           facecolors='none', edgecolors=color, s=18)
            ax.axvline(0, color='.4', lw=.7); ax.grid(alpha=.15)
            ax.set(xlabel=label+' displacement (km)', ylabel='Depth (m)' if j==0 else '',
                   title=regime.capitalize() if i==0 else '', ylim=(max_depth_m, 0))
            if found: ax.legend()
            else: ax.text(.5, .5, 'Insufficient contributors', ha='center', transform=ax.transAxes)
        limit = max(1e-6, max(abs(v) for ax in axes[i] for v in ax.get_xlim()))
        for ax in axes[i]: ax.set_xlim(-limit, limit)
    fig.suptitle('Mean constituent displacement relative to environmental PV gradient\n'
                 + f'95% track-bootstrap CI; open circles <{sparse_eddies} eddies')
    return fig, axes
