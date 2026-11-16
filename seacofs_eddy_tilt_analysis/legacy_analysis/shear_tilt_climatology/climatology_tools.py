"""Population shear-relative geometry and episode-safe tilt evolution.

Only numpy/pandas/matplotlib are required. TiltDir is a true-north bearing of
the measured lower-to-upper coherent tilt; angles here are Cartesian, positive
counterclockwise (left) from upper-minus-lower background velocity.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEYS = ["Cyc", "Eddy"]
GROUPS = ["regime", "Cyc"]
REGIMES = ("planetary", "topographic")
COLORS = {"AE": "firebrick", "CE": "royalblue"}
ANGLE_EDGES = np.arange(-180., 181., 45.)
LABELS = {
    "parallel_fraction": "Along-shear tilt / tilt magnitude",
    "left_fraction": "Left-of-shear tilt / tilt magnitude",
    "growth_km_day": "Tilt growth (km/day)",
    "dparallel_km_day": "Along-shear tilt change (km/day)",
    "dleft_km_day": "Left-of-shear tilt change (km/day)",
    "tilt_rotation_deg_day": "Tilt rotation (degrees/day; positive = left)",
    "shear_rotation_deg_day": "Shear rotation (degrees/day; positive = left)",
    "relative_rotation_deg_day": "Shear-relative rotation (degrees/day)",
    "fixed_shear_alignment_change_day": "Alignment change / day (fixed initial shear)",
    "alignment_change_day": "Alignment change / day (moving shear)",
    "dpv_parallel_km_day": "Along-PV-gradient tilt change (km/day)",
    "dpv_left_km_day": "Left-of-PV-gradient tilt change (km/day)",
}


def wrap(angle):
    return (np.asarray(angle) + 180.) % 360. - 180.


def validate_days(df):
    if df[KEYS + ["Day"]].isna().any().any():
        raise ValueError("Cyc, Eddy and Day must be non-missing")
    if df.duplicated(KEYS + ["Day"]).any():
        raise ValueError("Duplicate Cyc-Eddy-Day observations")
    days = pd.to_numeric(df.Day, errors="raise").to_numpy(float)
    if not np.isfinite(days).all() or not np.allclose(days, np.round(days)):
        raise ValueError("Day must contain finite integer model days")


def classify_regimes(df, factor=2., window=7, min_periods=5,
                     planetary_depth=3000., smoothing="trailing"):
    """Smooth log dominance only inside uninterrupted daily tracks.

    Trailing is primary: future observations cannot affect today's regime.
    Centred is a prespecified sensitivity matching older analyses. Invalid
    current-day ratios remain unknown even if neighbours are available.
    """
    validate_days(df)
    if factor <= 1 or window < 1 or not 1 <= min_periods <= window:
        raise ValueError("Invalid dominance/smoothing settings")
    if smoothing not in ("trailing", "centered"):
        raise ValueError("smoothing must be trailing or centered")
    out = df.sort_values(KEYS + ["Day"]).reset_index(drop=True).copy()
    out["Day"] = out.Day.astype(int)
    gap = out.groupby(KEYS, sort=False).Day.diff().ne(1)
    out["track_segment"] = gap.cumsum()
    # Recompute the ratio from magnitudes, avoiding stale cached sign conventions.
    p = pd.to_numeric(out.PV_grad_plan_mag, errors="coerce")
    t = pd.to_numeric(out.PV_grad_topo_mag, errors="coerce")
    good = np.isfinite(p) & np.isfinite(t) & (p > 0) & (t >= 0)
    # Zero topography is valid planetary dominance; cap only for rolling arithmetic.
    ratio = np.log(np.maximum(t, np.finfo(float).tiny)) - np.log(p.where(p > 0))
    out["log_ratio"] = ratio.where(good)
    out["log_ratio_smooth"] = out.groupby("track_segment").log_ratio.transform(
        lambda s: s.rolling(window, min_periods=min_periods,
                            center=smoothing == "centered").median())
    valid = good & out.log_ratio_smooth.notna() & np.isfinite(out.h)
    out["regime"] = "unknown"
    out.loc[valid, "regime"] = "mixed"
    out.loc[valid & out.log_ratio_smooth.le(-np.log(factor)) &
            out.h.ge(planetary_depth), "regime"] = "planetary"
    out.loc[valid & out.log_ratio_smooth.ge(np.log(factor)), "regime"] = "topographic"
    return out


def relative_geometry(df, family="clim", definition="upper_deep",
                      min_shear=.005, min_tilt=5., min_water_depth=500.):
    """Retain ALL rows for episode construction, marking unusable ones.

    Recover the 200-500 m mean only where the cached column covers 500 m.
    water_depth_m is centre bathymetry, not core-mean h. Full layer coverage
    and finite inputs are required; missing data never become zero shear.
    """
    if definition not in ("upper_deep", "surface_deep"):
        raise ValueError("Unknown shear definition")
    validate_days(df)
    out = df.sort_values(KEYS + ["Day"]).reset_index(drop=True).copy()
    for axis in ("east", "north"):
        a = out[f"{family}_200_{axis}_ms"]
        b = out[f"{family}_500_{axis}_ms"]
        deep = (500 * b - 200 * a) / 300
        upper = a if definition == "upper_deep" else out[f"{family}_surface_{axis}_ms"]
        out[f"shear_{axis}"] = upper - deep
    out["shear_ms"] = np.hypot(out.shear_east, out.shear_north)
    # Explicit bearing conversion; measured TiltDis/TiltDir are not modified.
    bearing = np.deg2rad(out.TiltDir)
    out["tilt_east"] = out.TiltDis * np.sin(bearing)
    out["tilt_north"] = out.TiltDis * np.cos(bearing)
    out["tilt_angle"] = np.degrees(np.arctan2(out.tilt_north, out.tilt_east))
    out["shear_angle"] = np.degrees(np.arctan2(out.shear_north, out.shear_east))
    valid = (out.TiltDis.ge(min_tilt) & out.TiltDis.gt(0) &
             out.TiltDir.between(0, 360, inclusive="both") &
             out.shear_ms.gt(min_shear) & out.water_depth_m.ge(min_water_depth))
    valid &= np.isfinite(out[["TiltDis", "TiltDir", "shear_ms"]]).all(axis=1)
    out["eligible"] = valid & out.regime.isin(REGIMES)
    den = out.shear_ms.where(valid)
    out["unit_east"] = out.shear_east / den
    out["unit_north"] = out.shear_north / den
    out["parallel_km"] = out.tilt_east * out.unit_east + out.tilt_north * out.unit_north
    out["left_km"] = -out.tilt_east * out.unit_north + out.tilt_north * out.unit_east
    out["parallel_fraction"] = out.parallel_km / out.TiltDis.where(out.TiltDis > 0)
    out["left_fraction"] = out.left_km / out.TiltDis.where(out.TiltDis > 0)
    out["angle_deg"] = wrap(out.tilt_angle - out.shear_angle)
    out.loc[~valid, "angle_deg"] = np.nan
    out["angle_bin"] = pd.cut(out.angle_deg, ANGLE_EDGES, right=False,
                              labels=(ANGLE_EDGES[:-1] + ANGLE_EDGES[1:]) / 2).astype(float)
    out["upshear"] = (out.parallel_fraction < 0).astype(float).where(valid)
    out["track_id"] = out.groupby(KEYS, sort=False).ngroup()
    same_track = out.track_id.eq(out.track_id.shift())
    continuous = (same_track & out.Day.diff().eq(1) & out.regime.eq(out.regime.shift()) &
                  out.eligible & out.eligible.shift(fill_value=False))
    out["episode"] = (~continuous).cumsum()
    return out


def lagged_pairs(df, lag=3):
    """Exact endpoints with every intervening day eligible in the same regime.

    Project displacement change on START shear axes. Decompose angular changes
    with wrapped endpoint rotations: relative = tilt rotation - shear rotation
    modulo 360. Angles cannot resolve complete turns between observations.
    """
    if not isinstance(lag, (int, np.integer)) or lag < 1:
        raise ValueError("lag must be a positive integer")
    validate_days(df)
    use = df.loc[df.eligible].copy()
    cols = KEYS + ["Day", "episode", "tilt_east", "tilt_north", "TiltDis",
                   "tilt_angle", "shear_angle", "angle_deg", "parallel_fraction"]
    future = use[cols].copy()
    future.Day -= lag
    future = future.rename(columns={c: "future_" + c for c in cols if c not in KEYS + ["Day"]})
    out = use.merge(future, on=KEYS + ["Day"], how="inner", validate="one_to_one")
    out = out.loc[out.episode.eq(out.future_episode)].copy()
    out["lag_days"] = lag
    de, dn = out.future_tilt_east - out.tilt_east, out.future_tilt_north - out.tilt_north
    out["dparallel_km_day"] = (de * out.unit_east + dn * out.unit_north) / lag
    out["dleft_km_day"] = (-de * out.unit_north + dn * out.unit_east) / lag
    out["growth_km_day"] = (out.future_TiltDis - out.TiltDis) / lag
    out["tilt_rotation_deg_day"] = wrap(out.future_tilt_angle - out.tilt_angle) / lag
    out["shear_rotation_deg_day"] = wrap(out.future_shear_angle - out.shear_angle) / lag
    out["relative_rotation_deg_day"] = wrap(out.future_angle_deg - out.angle_deg) / lag
    out["toward_downshear_deg_day"] = (np.abs(out.angle_deg) - np.abs(out.future_angle_deg)) / lag
    # Fixed start basis isolates tilt reorientation from changing shear direction.
    future_fixed = wrap(out.future_tilt_angle - out.shear_angle)
    out["fixed_shear_alignment_change_day"] = (
        np.cos(np.deg2rad(future_fixed)) - out.parallel_fraction) / lag
    out["alignment_change_day"] = (out.future_parallel_fraction - out.parallel_fraction) / lag
    if "PV_grad_theta" in out:
        pv = np.deg2rad(out.PV_grad_theta)
        out["dpv_parallel_km_day"] = (de * np.sin(pv) + dn * np.cos(pv)) / lag
        out["dpv_left_km_day"] = (-de * np.cos(pv) + dn * np.sin(pv)) / lag
    return out


def cluster_summary(df, metrics, groups=GROUPS, weighting="day",
                    n_boot=500, min_eddies=20, min_rows=50, seed=20260915):
    """Means and whole-eddy bootstrap CIs; preserve all qualifying observations.

    For day weighting, resample cluster sums/counts, NOT cluster means. Eddy
    weighting is an explicit alternative estimand. Counts are metric-specific.
    Cluster aggregation is sufficient statistics for inference only; it does
    not collapse episodes or change observations used in the point estimate.
    """
    if weighting not in ("day", "eddy") or n_boot < 1:
        raise ValueError("Invalid weighting or bootstrap count")
    rows = []
    rng = np.random.default_rng(seed)
    for key, part in df.groupby(list(groups), observed=True, dropna=True):
        key = key if isinstance(key, tuple) else (key,)
        for metric in metrics:
            q = part.loc[np.isfinite(part[metric]), KEYS + [metric]]
            agg = q.groupby(KEYS)[metric].agg(["sum", "count"])
            n = len(agg)
            supported = n >= min_eddies and len(q) >= min_rows
            point = lo = hi = np.nan
            if n:
                sums, counts = agg["sum"].to_numpy(), agg["count"].to_numpy()
                if weighting == "eddy":
                    sums, counts = sums / counts, np.ones(n)
                point = sums.sum() / counts.sum()
                if supported:
                    boot = np.empty(n_boot)
                    for b in range(n_boot):
                        ix = rng.integers(0, n, n)
                        boot[b] = sums[ix].sum() / counts[ix].sum()
                    lo, hi = np.quantile(boot, [.025, .975])
            rows.append(dict(zip(groups, key), metric=metric, mean=point,
                             ci_low=lo, ci_high=hi, eddies=n, observations=len(q),
                             supported=supported, weighting=weighting))
    return pd.DataFrame(rows, columns=list(groups) + ["metric", "mean", "ci_low", "ci_high",
                         "eddies", "observations", "supported", "weighting"])


def plot_angles(df, weighting="day"):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True,
                             constrained_layout=True)
    for i, regime in enumerate(REGIMES):
        for j, cyc in enumerate(("AE", "CE")):
            ax = axes[i, j]
            q = df.loc[df.eligible & df.regime.eq(regime) & df.Cyc.eq(cyc)].dropna(subset=["angle_deg"])
            if len(q):
                w = np.ones(len(q)) if weighting == "day" else 1 / q.groupby(KEYS).Day.transform("size")
                ax.hist(q.angle_deg, bins=np.arange(-180, 181, 15), weights=100 * w / np.sum(w),
                        color=COLORS[cyc], alpha=.75)
            ax.axvspan(-90, 90, color="grey", alpha=.08)
            ax.set(title=f'{regime.capitalize()} {cyc}: {len(q):,} days, {q[KEYS].drop_duplicates().shape[0]} eddies',
                   xlabel='Tilt relative to shear (degrees; positive = left)', ylabel='Probability per 15 degrees (%)',
                   xticks=[-180, -90, 0, 90, 180], xlim=(-180, 180))
    fig.suptitle(f'{weighting}-weighted orientation; shaded centre = downshear')
    return fig


def plot_response(summary, metric, ylabel=None):
    """Conditional response by initial orientation; unsupported bins are gaps."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True, constrained_layout=True)
    for ax, regime in zip(axes, REGIMES):
        for cyc in ("AE", "CE"):
            q = summary.loc[summary.regime.eq(regime) & summary.Cyc.eq(cyc) & summary.metric.eq(metric)].sort_values("angle_bin")
            if q.empty:
                continue
            q = q.set_index("angle_bin").reindex((ANGLE_EDGES[:-1] + ANGLE_EDGES[1:]) / 2)
            y = q["mean"].where(q.supported.eq(True))
            ax.plot(q.index, y, 'o-', color=COLORS[cyc], label=cyc)
            ax.fill_between(q.index.to_numpy(float), q.ci_low.to_numpy(float), q.ci_high.to_numpy(float),
                            color=COLORS[cyc], alpha=.15)
        ax.axhline(0, color='black', lw=.8)
        ax.set(title=regime.capitalize(), xlabel='Initial tilt relative to shear (degrees)',
               ylabel=ylabel or LABELS.get(metric, metric), xticks=[-180, -90, 0, 90, 180], xlim=(-180, 180))
        if ax.lines:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend()
    return fig
