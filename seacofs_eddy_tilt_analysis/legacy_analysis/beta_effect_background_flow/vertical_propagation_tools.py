"""Depth-specific eddy-centre propagation diagnostics from vertical profiles."""

from __future__ import annotations

import numpy as np
import pandas as pd


KM_DAY_TO_M_S = 1000.0 / 86400.0


def interpolate_centres(vertical, target_depths=(0.0, 200.0, 500.0)):
    """Interpolate each eddy-day centre horizontally to fixed physical depths."""

    rows = []
    for (eddy, day), profile in vertical.groupby(["Eddy", "Day"], sort=False):
        p = profile[["Depth", "xc", "yc"]].dropna().copy()
        p["depth_m"] = p["Depth"].abs()
        p = p.sort_values("depth_m").drop_duplicates("depth_m")
        if p.empty:
            continue
        z = p.depth_m.to_numpy(float)
        for target in target_depths:
            if target < z.min() or target > z.max():
                continue
            rows.append({
                "Eddy": eddy, "Day": day, "centre_depth_m": float(target),
                "centre_xc_km": float(np.interp(target, z, p.xc)),
                "centre_yc_km": float(np.interp(target, z, p.yc)),
            })
    return pd.DataFrame(rows)


def _rate(part, column, window=5):
    part = part.sort_values("Day")
    value = part[column].to_numpy(float)
    day = part["Day"].to_numpy(float)
    rate = np.full(len(part), np.nan)
    valid = np.isfinite(value) & np.isfinite(day)
    if valid.sum() >= 3:
        rate[valid] = np.gradient(value[valid], day[valid])
    return pd.Series(rate, index=part.index).rolling(window, center=True, min_periods=3).median()


def add_depth_propagation(centres, grid_angle, window=5):
    """Differentiate fixed-depth centre tracks and rotate them to east/north."""

    out = centres.sort_values(["Eddy", "centre_depth_m", "Day"]).copy()
    out["centre_dxc_km_day"] = np.nan
    out["centre_dyc_km_day"] = np.nan
    for indices in out.groupby(["Eddy", "centre_depth_m"], sort=False).groups.values():
        part = out.loc[indices]
        out.loc[indices, "centre_dxc_km_day"] = _rate(part, "centre_xc_km", window)
        out.loc[indices, "centre_dyc_km_day"] = _rate(part, "centre_yc_km", window)
    alpha = float(grid_angle)
    out["centre_east_ms"] = (
        out.centre_dxc_km_day * np.cos(alpha) - out.centre_dyc_km_day * np.sin(alpha)
    ) * KM_DAY_TO_M_S
    out["centre_north_ms"] = (
        out.centre_dxc_km_day * np.sin(alpha) + out.centre_dyc_km_day * np.cos(alpha)
    ) * KM_DAY_TO_M_S
    return out


def wide_depth_propagation(centres):
    """Return one row per eddy-day with fixed-depth velocity columns."""

    wide = centres.pivot(index=["Eddy", "Day"], columns="centre_depth_m",
                         values=["centre_east_ms", "centre_north_ms"])
    wide.columns = [f"centre_{component.removeprefix('centre_')}_{int(depth)}m"
                    for component, depth in wide.columns]
    return wide.reset_index()
