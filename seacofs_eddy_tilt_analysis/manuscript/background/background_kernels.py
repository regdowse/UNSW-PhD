"""Local monthly climatology, physical-layer coverage and track bootstrap kernels."""
from pathlib import Path
import os

import netCDF4 as nc
import numpy as np
import pandas as pd

ORIGIN = pd.Timestamp('1990-01-01')
FIELDS = ('u_surface', 'v_surface', 'u_200', 'v_200', 'u_500', 'v_500')


def clean(values):
    values = np.asarray(np.ma.filled(values, np.nan), dtype=float)
    return np.where(abs(values) > 1e30, np.nan, values)


def sigma_layer_weights(z_r, h, depth_max_m):
    """Original midpoint-interface approximation, clipped at the requested depth.

    Physical centres are static source-dataset depths, not level numbers. Surface
    and bottom interfaces are 0 and h. We do not infer extra coverage from missing
    velocity values; weighted_mean requires the complete requested thickness.
    """
    depth = abs(np.asarray(z_r, float))
    if np.nanmedian(depth[..., 0]) > np.nanmedian(depth[..., -1]):
        depth = depth[..., ::-1]
    if np.nanmedian(depth[..., 0]) >= np.nanmedian(depth[..., -1]):
        raise ValueError('Cannot establish surface-first sigma ordering')
    finite = np.isfinite(depth).all(axis=-1) & np.isfinite(h) & (h > 0)
    if (np.diff(depth[finite], axis=-1) <= 0).any():
        raise ValueError('Nonmonotonic physical sigma depths')
    edges = np.empty(depth.shape[:-1] + (depth.shape[-1]+1,))
    edges[..., 0] = 0
    edges[..., 1:-1] = (depth[..., :-1]+depth[..., 1:])/2
    edges[..., -1] = h
    weights = np.maximum(np.diff(np.minimum(edges, depth_max_m), axis=-1), 0)
    weights[~finite] = 0
    return weights


def weighted_mean(values, weights, thickness):
    """Return NaN unless every metre in the requested layer has finite velocity."""
    valid = np.isfinite(values) & (weights > 0)
    coverage = np.sum(np.where(valid, weights, 0), axis=-1)
    numerator = np.sum(np.where(valid, values, 0)*weights, axis=-1)
    return np.divide(numerator, thickness, out=np.full(coverage.shape, np.nan),
                     where=np.isclose(coverage, thickness, rtol=0, atol=1e-7))


def window_weights(day, window):
    if window < 1 or window % 2 != 1:
        raise ValueError('Window must be a positive odd integer')
    date = ORIGIN + pd.Timedelta(days=int(day))
    # Original circular 365-day climatological calendar: Feb 29 maps to Feb 28.
    ref = pd.Timestamp(2001, date.month, min(date.day, 28) if date.month == 2 else date.day)
    dates = ref + pd.to_timedelta(np.arange(-(window//2), window//2+1), unit='D')
    return np.bincount(dates.month, minlength=13)[1:].astype(float)


def attach(table, climatology, window):
    out = table[['Eddy', 'Day', 'ic', 'jc']].copy()
    i, j = out.ic.to_numpy(int), out.jc.to_numpy(int)
    unique_days = np.sort(out.Day.unique())
    by_day = {day: window_weights(day, window) for day in unique_days}
    weights = np.vstack([by_day[day] for day in out.Day])
    for name in FIELDS:
        values = np.column_stack([climatology[f'm{m:02d}_{name}'][i,j] for m in range(1,13)])
        covered = np.isfinite(values) | (weights == 0)
        mean = np.sum(np.where(np.isfinite(values),values,0)*weights, axis=1)/window
        out[name] = np.where(covered.all(axis=1), mean, np.nan)
    out['u_200_500'] = (500*out.u_500-200*out.u_200)/300
    out['v_200_500'] = (500*out.v_500-200*out.v_200)/300
    out['surface_minus_deep_east_ms'] = out.u_surface-out.u_200_500
    return out


def atomic_npz(path, **arrays):
    path = Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix('.tmp')
    with temporary.open('wb') as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary,path)


def process_file(path, weights, partition_root):
    """Read one entire model file; fresh restartable monthly sums and counts."""
    target = Path(partition_root)/(Path(path).stem+'.npz')
    if target.exists(): return str(target), 'resumed'
    sums, counts = {}, {}
    with nc.Dataset(path) as ds:
        days = np.rint(clean(ds['ocean_time'][:])/86400).astype(int)
        for t, day in enumerate(days):
            # Geographic velocities are already on rho cells; storage is bottom-first.
            u = clean(ds['u_eastward'][t]).T[..., ::-1]
            v = clean(ds['v_northward'][t]).T[..., ::-1]
            fields = [u[...,0],v[...,0]]
            for depth in (200,500):
                fields.extend([weighted_mean(u,weights[depth],depth),
                               weighted_mean(v,weights[depth],depth)])
            month = (ORIGIN+pd.Timedelta(days=int(day))).month
            for name, field in zip(FIELDS, fields):
                key = f'm{month:02d}_{name}'
                if key not in sums:
                    sums[key]=np.zeros(field.shape);counts[key]=np.zeros(field.shape,dtype=np.int32)
                valid = np.isfinite(field)
                sums[key][valid] += field[valid];counts[key][valid] += 1
    atomic_npz(target, **{**{'sum_'+k:v for k,v in sums.items()},
                          **{'count_'+k:v for k,v in counts.items()}})
    return str(target), 'computed'


def reduce_partitions(paths, output):
    sums, counts = {}, {}
    for path in paths:
        with np.load(path) as part:
            for key in part.files:
                if not key.startswith('sum_'): continue
                suffix=key[4:]
                if suffix not in sums:
                    sums[suffix]=np.zeros_like(part[key]);counts[suffix]=np.zeros_like(part['count_'+suffix],dtype=np.int64)
                sums[suffix] += part[key];counts[suffix] += part['count_'+suffix]
    means = {k:np.divide(v,counts[k],out=np.full_like(v,np.nan),where=counts[k]>0) for k,v in sums.items()}
    expected = {f'm{m:02d}_{name}' for m in range(1,13) for name in FIELDS}
    if set(means) != expected: raise ValueError('Missing climatology months/fields')
    atomic_npz(output, **means, **{'count_'+k:v for k,v in counts.items()})
    return means


def cluster_median(frame, column, n_boot=5000, seed=42):
    """Pooled eddy-day median; whole-eddy bootstrap preserves track dependence."""
    q=frame[['Eddy',column]].replace([np.inf,-np.inf],np.nan).dropna()
    if q.empty: return np.nan,np.nan,np.nan
    clusters=[p[column].to_numpy(float) for _,p in q.groupby('Eddy',sort=True)]
    estimate=float(np.median(q[column]))
    if len(clusters)<2:return estimate,np.nan,np.nan
    rng=np.random.default_rng(seed);boot=np.empty(n_boot)
    for b in range(n_boot):
        chosen=rng.choice(len(clusters),len(clusters),replace=True)
        boot[b]=np.median(np.concatenate([clusters[k] for k in chosen]))
    return estimate,*np.quantile(boot,[.025,.975])
