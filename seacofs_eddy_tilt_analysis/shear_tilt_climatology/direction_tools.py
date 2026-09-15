"""Shear compass-direction tests; all angles Cartesian unless marked bearing."""
import numpy as np
import pandas as pd
import climatology_tools as ct

POPULATIONS = ('planetary', 'topographic', 'mixed', 'all')
SECTORS = ('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')


def add_directions(frame):
    out = frame.copy()
    out['shear_bearing'] = (90 - out.shear_angle) % 360
    out['sector'] = np.floor(((out.shear_bearing + 22.5) % 360) / 45).astype(int)
    out['left_probability'] = (out.left_fraction > 0).astype(float)
    # This explicitly signed metric is secondary; raw AE/CE signs remain primary.
    out['expected_side'] = np.where(out.Cyc.eq('AE'), out.left_fraction, -out.left_fraction)
    out['expected_probability'] = np.where(out.Cyc.eq('AE'), out.left_probability,
                                          (out.left_fraction < 0).astype(float))
    return out


def populations(data):
    # Geometry validity is independent of PV label. Includes unknown PV in all.
    good = data.loc[np.isfinite(data.left_fraction) & np.isfinite(data.parallel_fraction)].copy()
    good = add_directions(good)
    pieces = [good.assign(population='all')]
    for regime in POPULATIONS[:-1]:
        pieces.append(good.loc[good.regime.eq(regime)].assign(population=regime))
    return pd.concat(pieces, ignore_index=True)


def rebuild_episodes(data, population):
    """Rebuild on all rows: invalid days and calendar gaps cannot be bridged.

    All-population episodes may cross PV regimes intentionally; regime-specific
    episodes cannot. The original regime labels are retained.
    """
    q = data.sort_values(ct.KEYS + ['Day']).copy()
    valid = np.isfinite(q.left_fraction) & np.isfinite(q.parallel_fraction)
    q['eligible'] = valid if population == 'all' else valid & q.regime.eq(population)
    same = q[ct.KEYS].eq(q[ct.KEYS].shift()).all(axis=1)
    continuous = same & q.Day.diff().eq(1) & q.eligible & q.eligible.shift(fill_value=False)
    q['episode'] = (~continuous).cumsum()
    q['population'] = population
    return q


def shear_windows(data, length=14):
    """Nonoverlapping fixed-duration windows avoid track-length range bias."""
    rows = []
    for population in POPULATIONS:
        q = rebuild_episodes(data, population)
        for _, episode in q.loc[q.eligible].groupby('episode'):
            for start in range(0, len(episode) - length + 1, length):
                w = episode.iloc[start:start + length]
                angles = np.deg2rad(w.shear_angle.to_numpy())
                resultant = abs(np.exp(1j * angles).mean())
                turns = np.abs(ct.wrap(np.diff(w.shear_angle.to_numpy())))
                rows.append(dict(population=population, Cyc=w.Cyc.iloc[0], Eddy=w.Eddy.iloc[0],
                    Day=w.Day.iloc[0], circular_variance=1-resultant,
                    mean_abs_daily_turn=turns.mean(),
                    net_turn=abs(ct.wrap(w.shear_angle.iloc[-1]-w.shear_angle.iloc[0])),
                    median_lat=w.lat.median()))
    return pd.DataFrame(rows, columns=['population','Cyc','Eddy','Day','circular_variance',
                                     'mean_abs_daily_turn','net_turn','median_lat'])


def turning_pairs(data, lag=7):
    pieces = []
    for population in POPULATIONS:
        p = ct.lagged_pairs(rebuild_episodes(data, population), lag)
        p['abs_shear_turn'] = abs(p.shear_rotation_deg_day * lag)
        p['turn_bin'] = pd.cut(p.abs_shear_turn, [-.001,15,45,90,180],
                               labels=['0–15','15–45','45–90','90–180']).astype(object)
        p['start_left'] = p.left_fraction
        p['end_left'] = np.sin(np.deg2rad(p.future_angle_deg))
        # Compare new shear with the frozen start direction, on identical endpoints.
        p['end_left_old_shear'] = np.sin(np.deg2rad(ct.wrap(p.future_tilt_angle-p.shear_angle)))
        sign = np.where(p.Cyc.eq('AE'), 1., -1.)
        p['new_minus_old_expected'] = sign * (p.end_left-p.end_left_old_shear)
        p['both_expected_side'] = ((sign*p.start_left > 0) & (sign*p.end_left > 0)).astype(float)
        pieces.append(p)
    return pd.concat(pieces, ignore_index=True)


def sector_balanced(frame, n_boot=500, min_eddies=20, min_rows=50, seed=20260916):
    """Equal compass-sector mean with a JOINT whole-eddy bootstrap.

    Requires all eight sectors supported; cannot quietly drop missing sectors.
    A replicate lacking any sector is invalid and reported, not reweighted.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for (population, cyc), q in frame.groupby(['population','Cyc']):
        q = q.copy(); q['id'] = q.groupby(ct.KEYS).ngroup()
        sums = q.pivot_table(index='id', columns='sector', values='left_fraction', aggfunc='sum').reindex(columns=range(8)).fillna(0).to_numpy()
        counts = q.pivot_table(index='id', columns='sector', values='left_fraction', aggfunc='count').reindex(columns=range(8)).fillna(0).to_numpy()
        support = (counts.sum(axis=0)>=min_rows) & ((counts>0).sum(axis=0)>=min_eddies)
        mean = np.mean(sums.sum(axis=0)/counts.sum(axis=0)) if (counts.sum(axis=0)>0).all() else np.nan
        draws = []
        if support.all():
            for _ in range(n_boot):
                ix = rng.integers(0,len(sums),len(sums)); den = counts[ix].sum(axis=0)
                if (den>0).all(): draws.append(np.mean(sums[ix].sum(axis=0)/den))
        lo, hi = np.quantile(draws,[.025,.975]) if len(draws)>=.95*n_boot else (np.nan,np.nan)
        rows.append(dict(population=population,Cyc=cyc,mean=mean,ci_low=lo,ci_high=hi,
            supported_sectors=int(support.sum()),valid_bootstraps=len(draws),eddies=len(sums)))
    return pd.DataFrame(rows)
