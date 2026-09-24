"""Small, independent sensitivity experiment; production tilt code is unchanged."""
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

# METHODS = {'Equal': None, 'Gaussian 2 d': 2.0, 'Gaussian 1.5 d': 1.5}
METHODS = {
    'Equal': None,
    'Gaussian 1.5 d': 1.5,
    'Gaussian 1.0 d': 1.0,
}

def increments(track, depth_int=10, max_depth=1000):
    """Calendar-day columns and upper-interval depth labels match production."""
    days = np.arange(int(track.Day.min()), int(track.Day.max()) + 1)
    target = np.arange(0, max_depth + depth_int, depth_int)
    x, y = {}, {}
    for day in days:
        g = track.loc[track.Day.eq(day)].copy()
        g['Depth'] = g.Depth.abs()
        g = g.loc[g.Depth.le(max_depth)].sort_values('Depth')
        if g.Depth.duplicated().any():
            raise ValueError(f'Duplicate depths on Day {day}')
        if not np.isfinite(g[['Depth', 'xc', 'yc']].to_numpy()).all():
            raise ValueError(f'Nonfinite profile on Day {day}')
        levels = target[(target >= g.Depth.min()) & (target <= g.Depth.max())]
        if len(g) < 2 or len(levels) < 2:
            x[day] = pd.Series(np.nan, index=target[:-1])
            y[day] = pd.Series(np.nan, index=target[:-1])
        else:
            x[day] = pd.Series(np.diff(np.interp(levels, g.Depth, g.xc)), index=levels[:-1])
            y[day] = pd.Series(np.diff(np.interp(levels, g.Depth, g.yc)), index=levels[:-1])
    return pd.DataFrame(x), pd.DataFrame(y)


def fit_snapshot(dx, dy, day, sigma=None, half_window=3, min_depth_range=200,
                 min_points=5, bearing_offset=20.0, eps=1e-10):
    days = np.arange(day - half_window, day + half_window + 1)
    if days[0] < dx.columns.min() or days[-1] > dx.columns.max():
        return None
    x, y = dx.reindex(columns=days), dy.reindex(columns=days)
    a = pd.Series(1.0 if sigma is None else np.exp(-0.5*((days-day)/sigma)**2), index=days)
    def mean(d):
        return d.mul(a, axis=1).sum(axis=1, min_count=1) / d.notna().mul(a, axis=1).sum(axis=1).replace(0, np.nan)
    trace = pd.DataFrame({'x': mean(x).cumsum(), 'y': mean(y).cumsum(),
                          'weight': 1/(x.var(axis=1) + y.var(axis=1) + eps)})
    trace['Depth'] = trace.index
    trace = trace.replace([np.inf, -np.inf], np.nan).dropna()
    if len(trace) < min_points or trace.Depth.max()-trace.Depth.min() < min_depth_range:
        return None
    points = trace[['x', 'y', 'Depth']].to_numpy()
    w = trace.weight.to_numpy()/trace.weight.max()
    centre = np.average(points, axis=0, weights=w)
    _, _, vt = np.linalg.svd((points-centre)*np.sqrt(w[:, None]), full_matrices=False)
    v = vt[0]
    if not np.isfinite(v).all() or abs(v[2]) < eps:
        return None
    ends = centre + ((np.array([trace.Depth.min(), trace.Depth.max()])-centre[2])/v[2])[:, None]*v
    delta = ends[0, :2]-ends[1, :2]
    return dict(trace=trace, ends=ends, TiltDis=float(np.linalg.norm(delta)),
                TiltDir=float((np.degrees(np.arctan2(delta[0], delta[1]))+bearing_offset)%360))


def maximum_span(profile, max_depth=1000):
    """Maximum pairwise horizontal distance, not surface-to-deep displacement."""
    g = profile.loc[profile.Depth.abs().le(max_depth)]
    xy = g[['xc', 'yc']].to_numpy()
    return float(pdist(xy).max()) if len(xy) >= 2 and np.isfinite(xy).all() else np.nan


def lifetime(track, depth_int=10, max_depth=1000):
    dx, dy = increments(track, depth_int, max_depth)
    spans = {d: maximum_span(g, max_depth) for d, g in track.groupby('Day')}
    rows = []
    for day in dx.columns:
        row = {'Day': day, 'Maximum span': spans.get(day, np.nan)}
        for name, sigma in METHODS.items():
            fit = fit_snapshot(dx, dy, day, sigma)
            row[name] = np.nan if fit is None else fit['TiltDis']
        rows.append(row)
    return pd.DataFrame(rows), dx, dy
