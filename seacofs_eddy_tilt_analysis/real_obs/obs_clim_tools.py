"""Small, explicit observation/climatology comparison helpers.

Depth is positive down. Bearings are geographic compass, deep-to-shallow.
No production tilt or shared cache is modified.
"""
import numpy as np
import pandas as pd

KEYS = ['Eddy', 'Day']


def coriolis(lat):
    return 2 * 7.292115e-5 * np.sin(np.deg2rad(lat))


def angle_difference(a, b):
    return (np.asarray(a) - np.asarray(b) + 180) % 360 - 180


def unique(df, keys):
    if df.duplicated(keys).any():
        raise ValueError(f'Duplicate keys: {keys}')


def observations():
    """Figure 11 peaks are approximate visual readings, not digitised data."""
    out = pd.DataFrame([
        dict(name='Murphy', lat=-28.6168, lon=154.8351, radius_km=80.,
             paper_Ro_peak=0.37, rate_km_per_100m=3.3, bearing_deg=248.8,
             published_distance_km=23., date='2015-06-05 03:30'),
        dict(name='Freddy', lat=-32.5763, lon=153.2057, radius_km=17.5,
             paper_Ro_peak=0.65, rate_km_per_100m=1.5, bearing_deg=250.3,
             published_distance_km=10., date='2015-06-09 00:15'),
    ])
    out['f_s-1'] = coriolis(out.lat)
    out['paper_Omega_s-1'] = out.paper_Ro_peak * out['f_s-1']
    # V_theta=Omega*r -> curl(u)=2*Omega in a circular solid-body core.
    out['zeta_s-1'] = 2 * out['paper_Omega_s-1']
    out['Ro_zeta_target'] = 2 * out.paper_Ro_peak
    out['equivalent_600m_km'] = 6 * out.rate_km_per_100m
    return out


def profile_metrics(vertical, *, target_depth=600., depth_tolerance=100.,
                    fit_min_depth=50., rotation_rad=0.):
    """OLS x(z),y(z) at actual levels; no temporal smoothing/interpolation.

    Use the shallowest valid level >= target, within tolerance; never a
    shallower endpoint. Peak signed w comes strictly from 0..target_depth.
    OLS is not the paper's inverse-error-variance fit (errors unavailable).
    """
    unique(vertical, KEYS + ['Depth'])
    records = []
    for (eddy, day), p in vertical.groupby(KEYS, sort=False):
        p = p.sort_values('Depth')
        p = p.loc[np.isfinite(p[['Depth', 'xc', 'yc', 'w']]).all(axis=1)
                  & p.Depth.ge(0)].copy()
        row = dict(Eddy=eddy, Day=day, profile_ok=False, profile_status='no valid levels')
        if p.empty:
            records.append(row)
            continue
        row['max_valid_depth_m'] = p.Depth.max()
        peak = p.loc[p.Depth.le(target_depth)]
        if not peak.empty:
            peak = peak.iloc[np.argmax(np.abs(peak.w.to_numpy()))]
            row.update(w_peak=float(peak.w), w_peak_depth_m=float(peak.Depth))
        deep = p.loc[p.Depth.ge(target_depth) & p.Depth.le(target_depth + depth_tolerance)]
        if deep.empty:
            row['profile_status'] = 'no endpoint at/just below target'
        else:
            endpoint = deep.Depth.iloc[0]
            fit = p.loc[p.Depth.between(fit_min_depth, endpoint)]
            row.update(fit_top_m=fit.Depth.min(), fit_bottom_m=endpoint, fit_n=len(fit))
            if len(fit) < 3 or fit.Depth.min() > fit_min_depth + 100:
                row['profile_status'] = 'insufficient upper-column fit support'
            else:
                z = fit.Depth.to_numpy(float)
                xy = fit[['xc', 'yc']].to_numpy(float)
                zc = z - z.mean()
                slope = zc @ (xy - xy.mean(axis=0)) / (zc @ zc)
                residual = xy - (xy.mean(axis=0) + zc[:, None] * slope)
                sst = np.sum((xy - xy.mean(axis=0))**2)
                r2 = 1 - np.sum(residual**2) / sst if sst > 0 else np.nan
                # Reverse increasing-depth slope, then rotate grid axes to E/N.
                dx, dy = -slope * 100
                east = dx*np.cos(rotation_rad) - dy*np.sin(rotation_rad)
                north = dx*np.sin(rotation_rad) + dy*np.cos(rotation_rad)
                rate = float(np.hypot(east, north))
                row.update(profile_ok=True, profile_status='usable',
                           rate_km_per_100m=rate, fit_r2=r2,
                           bearing_deg=float(np.degrees(np.arctan2(east, north)) % 360)
                           if rate > 1e-10 else np.nan,
                           equivalent_600m_km=rate*6,
                           fit_span_m=float(z.max()-z.min()))
        records.append(row)
    columns = KEYS + ['profile_ok', 'profile_status', 'max_valid_depth_m', 'w_peak',
                     'w_peak_depth_m', 'fit_top_m', 'fit_bottom_m', 'fit_n', 'fit_r2',
                     'rate_km_per_100m', 'bearing_deg', 'equivalent_600m_km', 'fit_span_m']
    return pd.DataFrame(records).reindex(columns=columns)


def match(surface, metrics, obs, distance_function, *, radius_factor=1.5,
          ro_fraction=0.25):
    unique(surface, KEYS)
    unique(metrics, KEYS)
    s = surface.merge(metrics, on=KEYS, how='left', validate='one_to_one')
    s['Ro_peak'] = s.w_peak / coriolis(s.lat)
    groups, audits = {}, []
    for o in obs.itertuples(index=False):
        d = s.copy()
        d['distance_to_obs_km'] = distance_function(d.lat, d.lon, o.lat, o.lon)
        stages = [
            ('CE', d.Cyc.eq('CE')),
            ('within search circle', d.distance_to_obs_km.le(radius_factor*o.radius_km)),
            ('usable >=600m profile', d.profile_ok.eq(True)),
            ('peak Ro within tolerance', d.Ro_peak.between(
                o.Ro_zeta_target*(1-ro_fraction), o.Ro_zeta_target*(1+ro_fraction))),
        ]
        keep = pd.Series(True, index=d.index)
        for label, condition in stages:
            keep &= condition
            audits.append(dict(observation=o.name, stage=label,
                               days=int(keep.sum()), eddies=d.loc[keep, 'Eddy'].nunique()))
        groups[o.name] = d.loc[keep].copy()
    return groups, pd.DataFrame(audits)


def summary(groups, obs):
    rows = []
    for o in obs.itertuples(index=False):
        d = groups[o.name]
        for col, observed in [('rate_km_per_100m', o.rate_km_per_100m),
                              ('equivalent_600m_km', o.equivalent_600m_km),
                              ('TiltDis', o.published_distance_km)]:
            v = d[col].dropna()
            q = v.quantile([.25, .5, .75])
            rows.append(dict(observation=o.name, metric=col, observed=observed,
                             days=len(v), eddies=d.loc[v.index, 'Eddy'].nunique(),
                             model_q25=q.get(.25), model_median=q.get(.5), model_q75=q.get(.75),
                             observed_percentile=100*((v<observed).sum()+.5*(v==observed).sum())/len(v)
                             if len(v) else np.nan))
    return pd.DataFrame(rows)


def plot_comparison(groups, obs):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(obs), 2, figsize=(11, 3.5*len(obs)), squeeze=False)
    for ax, o in zip(axes, obs.itertuples(index=False)):
        d = groups[o.name]
        angles = angle_difference(d.bearing_deg.dropna(), o.bearing_deg)
        ax[0].hist(angles, bins=np.linspace(-180,180,25), color='steelblue', alpha=.8)
        ax[0].axvline(0, color='black', linestyle='--', label='Observed direction')
        ax[0].set(xlim=(-180,180), xlabel='Model minus observed bearing (degrees)', ylabel='Eddy-days')
        ax[1].hist(d.rate_km_per_100m.dropna(), bins=20, color='steelblue', alpha=.8)
        ax[1].axvline(o.rate_km_per_100m, color='black', linestyle='--', label='Observation')
        ax[1].set(xlabel='Upper-column tilt (km per 100 m)', ylabel='Eddy-days')
        for a in ax:
            a.set_title(f'{o.name}: {len(d):,} days, {d.Eddy.nunique():,} eddies')
            a.legend()
        if d.empty:
            for a in ax:
                a.text(.5,.5,'No matches; do not infer agreement', transform=a.transAxes, ha='center')
    fig.tight_layout()
    return fig
