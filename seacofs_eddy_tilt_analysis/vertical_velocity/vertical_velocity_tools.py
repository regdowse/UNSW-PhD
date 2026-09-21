"""Exploratory native-grid vertical-velocity sampling around fitted eddy cores."""

from pathlib import Path
from types import SimpleNamespace

import netCDF4 as nc
import numpy as np
import pandas as pd


def model_time_index(ds, day):
    days = np.rint(np.asarray(ds["ocean_time"][:], float) / 86400).astype(int)
    hits = np.flatnonzero(days == int(day))
    if len(hits) != 1:
        raise ValueError(f"Expected one sample for day {day}; found {len(hits)}")
    return int(hits[0])


def interface_depths(ds, time_index, ii, jj):
    """ROMS physical depth (negative metres) at local s_w interfaces."""
    i0, i1, j0, j1 = ii.min(), ii.max() + 1, jj.min(), jj.max() + 1
    h = np.asarray(ds["h"][i0:i1, j0:j1], float)[ii-i0, jj-j0]
    zeta = np.asarray(ds["zeta"][time_index, i0:i1, j0:j1], float)[ii-i0, jj-j0]
    s = np.asarray(ds["s_w"][:], float)[:, None]
    cs = np.asarray(ds["Cs_w"][:], float)[:, None]
    hc = float(np.asarray(ds["hc"][:]))
    transform = int(np.asarray(ds["Vtransform"][:]))
    if transform == 2:
        z0 = (hc * s + h[None, :] * cs) / (hc + h[None, :])
        return zeta[None, :] + (zeta[None, :] + h[None, :]) * z0
    if transform == 1:
        z0 = hc * (s - cs) + h[None, :] * cs
        return z0 + zeta[None, :] * (1 + z0 / h[None, :])
    raise ValueError(f"Unsupported Vtransform={transform}")


def core_indices(row, grid, fraction):
    from seacofs_tilt_tools import core_grid_indices

    local_grid = SimpleNamespace(x_grid=grid.x_grid, y_grid=grid.y_grid,
                                 mask_rho=grid.mask_rho)
    return core_grid_indices(row, local_grid, frac=float(fraction))


def sample_snapshot(ds, row, profile_depths, grid, fractions=(0.5, 1, 1.5),
                    max_depth_m=1000, max_mismatch_m=75):
    """Return per-depth extrema and selected maps for one fitted eddy-day.

    Samples the nearest *local* s_w interface at each rho cell. Full-column
    extrema are taken over the per-depth samples, not over interpolated data.
    """
    t = model_time_index(ds, row.Day)
    wvar = ds["w"]
    if tuple(wvar.dimensions[1:]) != ("s_w", "eta_rho", "xi_rho"):
        raise ValueError(f"Unexpected w dimensions: {wvar.dimensions}")
    if wvar.shape[2:] != grid.mask_rho.shape:
        raise ValueError("Model w and analysis grid have different horizontal shapes")
    results, maps = [], {}
    for fraction in fractions:
        ii, jj = core_indices(row, grid, fraction)
        if len(ii) == 0:
            continue
        i0, i1, j0, j1 = ii.min(), ii.max() + 1, jj.min(), jj.max() + 1
        local_i, local_j = ii - i0, jj - j0
        raw = wvar[t, :, i0:i1, j0:j1]
        vel = np.asarray(np.ma.filled(raw, np.nan), float)
        vel[np.abs(vel) >= 1e30] = np.nan
        z = interface_depths(ds, t, ii, jj)
        for depth in profile_depths:
            if depth < 0 or depth > max_depth_m:
                continue
            nearest = np.argmin(np.abs(z + depth), axis=0)
            actual_z = z[nearest, np.arange(len(ii))]
            mismatch = np.abs(actual_z + depth)
            values = vel[nearest, local_i, local_j]
            valid = np.isfinite(values) & np.isfinite(mismatch) & (mismatch <= max_mismatch_m)
            record = dict(Eddy=int(row.Eddy), Day=int(row.Day), fraction=float(fraction),
                          Depth=float(depth), n_core=len(ii), n_valid=int(valid.sum()),
                          coverage=float(valid.mean()), max_mismatch_m=float(np.nanmax(mismatch)))
            field = np.full((i1-i0, j1-j0), np.nan)
            field[local_i[valid], local_j[valid]] = values[valid]
            if valid.any():
                v = values[valid]
                pos = np.flatnonzero(valid)
                hi, lo = pos[np.argmax(v)], pos[np.argmin(v)]
                record.update(w_max=float(values[hi]), w_min=float(values[lo]),
                              w_mean=float(np.mean(v)), w_abs_max=float(values[pos[np.argmax(np.abs(v))]]),
                              x_max=float(grid.x_grid[ii[hi]]), y_max=float(grid.y_grid[jj[hi]]),
                              x_min=float(grid.x_grid[ii[lo]]), y_min=float(grid.y_grid[jj[lo]]),
                              z_max_m=float(-actual_z[hi]), z_min_m=float(-actual_z[lo]))
            results.append(record)
            maps[(float(fraction), float(depth))] = (field, (i0, i1, j0, j1))
    return pd.DataFrame(results), maps


def column_extrema(per_depth):
    """One upward and one downward extreme per eddy-day and core fraction."""
    if per_depth.empty:
        return pd.DataFrame()
    valid = per_depth.dropna(subset=["w_max", "w_min"])
    rows = []
    for (eddy, day, fraction), group in valid.groupby(["Eddy", "Day", "fraction"]):
        up = group.loc[group.w_max.idxmax()]
        down = group.loc[group.w_min.idxmin()]
        rows.append(dict(Eddy=eddy, Day=day, fraction=fraction,
                         w_max=up.w_max, max_depth_m=up.z_max_m,
                         w_min=down.w_min, min_depth_m=down.z_min_m,
                         w_abs_max=up.w_max if up.w_max >= abs(down.w_min) else down.w_min,
                         n_depths=len(group), min_coverage=group.coverage.min()))
    return pd.DataFrame(rows)
