"""Statistical helpers for background-relative eddy propagation."""

from __future__ import annotations

import numpy as np
import pandas as pd


KM_DAY_TO_M_S = 1000.0 / 86400.0
METHODS = {
    "ann_surface": ("ann_surface_east_ms", "ann_surface_north_ms"),
    "ann_200": ("ann_200_east_ms", "ann_200_north_ms"),
    "ann_500": ("ann_500_east_ms", "ann_500_north_ms"),
    "clim_surface": ("clim_surface_east_ms", "clim_surface_north_ms"),
    "clim_200": ("clim_200_east_ms", "clim_200_north_ms"),
    "clim_500": ("clim_500_east_ms", "clim_500_north_ms"),
    "full_surface": ("full_surface_east_ms", "full_surface_north_ms"),
    "full_200": ("full_200_east_ms", "full_200_north_ms"),
    "full_500": ("full_500_east_ms", "full_500_north_ms"),
}


def _smoothed_rate(part, column, window=5):
    """Time-aware position derivative followed by centred robust smoothing."""

    part = part.sort_values("Day")
    value = part[column].to_numpy(dtype=float)
    time = part["Day"].to_numpy(dtype=float)
    rate = np.full(len(part), np.nan)
    valid = np.isfinite(value) & np.isfinite(time)
    if valid.sum() >= 3:
        rate[valid] = np.gradient(value[valid], time[valid])
    rate = pd.Series(rate, index=part.index)
    return rate.rolling(window, center=True, min_periods=3).median()


def add_track_velocity(df, grid_angle, window=5):
    """Add smoothed geographic track velocities from rotated-grid xc/yc."""

    out = df.sort_values(["Eddy", "Day"]).copy()
    out["dxc_km_day"] = np.nan
    out["dyc_km_day"] = np.nan
    for indices in out.groupby("Eddy", sort=False).groups.values():
        part = out.loc[indices]
        out.loc[indices, "dxc_km_day"] = _smoothed_rate(part, "xc", window)
        out.loc[indices, "dyc_km_day"] = _smoothed_rate(part, "yc", window)
    alpha = float(grid_angle)
    out["track_east_ms"] = (
        out["dxc_km_day"] * np.cos(alpha) - out["dyc_km_day"] * np.sin(alpha)
    ) * KM_DAY_TO_M_S
    out["track_north_ms"] = (
        out["dxc_km_day"] * np.sin(alpha) + out["dyc_km_day"] * np.cos(alpha)
    ) * KM_DAY_TO_M_S
    return out


def add_residual_velocities(df):
    """Subtract each background estimate from geographic track velocity."""

    out = df.copy()
    out["expected_sign"] = np.where(out["Cyc"].eq("AE"), 1.0, -1.0)
    for method, (east, north) in METHODS.items():
        out[f"{method}_residual_east_ms"] = out["track_east_ms"] - out[east]
        out[f"{method}_residual_north_ms"] = out["track_north_ms"] - out[north]
        out[f"{method}_aligned_residual_north_ms"] = (
            out["expected_sign"] * out[f"{method}_residual_north_ms"]
        )
    return out


def eddy_bootstrap_ci(df, column, group="Cyc", n_boot=5000, seed=42):
    """Median and whole-eddy bootstrap CI, so days are not pseudoreplicates."""

    rng = np.random.default_rng(seed)
    records = []
    for label, part in df.groupby(group):
        values = part.groupby("Eddy")[column].median().dropna().to_numpy()
        if not len(values):
            continue
        boot = np.median(rng.choice(values, (n_boot, len(values)), replace=True), axis=1)
        records.append({
            group: label,
            "column": column,
            "eddies": len(values),
            "median": float(np.median(values)),
            "ci_low": float(np.percentile(boot, 2.5)),
            "ci_high": float(np.percentile(boot, 97.5)),
        })
    return pd.DataFrame(records)


def propagation_summary(df, n_boot=5000, vec='north'):
    """Whole-eddy estimates for raw and all background-relative velocities."""

    columns = [f"track_{vec}_ms"] + [
        f"{method}_residual_{vec}_ms" for method in METHODS
    ]
    return pd.concat(
        [eddy_bootstrap_ci(df, column, n_boot=n_boot) for column in columns],
        ignore_index=True,
    )


def vorticity_budget_summary(df, method="ann_surface", n_boot=5000):
    """Compare observed dw/dt with the planetary-advection term -beta*v_res."""

    out = df.sort_values(["Eddy", "Day"]).copy()
    out["dw_dt_s2"] = out.groupby("Eddy")["w"].diff() / (
        out.groupby("Eddy")["Day"].diff() * 86400.0
    )
    out["beta_advection_s2"] = -out["beta"] * out[f"{method}_residual_north_ms"]
    out["budget_residual_s2"] = out["dw_dt_s2"] - out["beta_advection_s2"]
    return out, eddy_bootstrap_ci(out, "budget_residual_s2", n_boot=n_boot)
