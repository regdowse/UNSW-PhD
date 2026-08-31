"""Surface/depth ellipse geometry and eddy-clustered descriptive tilt analysis.

Q defines rho**2 = [dx, dy] Q [dx, dy]. Its smallest eigenvalue
therefore identifies the major axis. Bearings follow core/tilt.py, not
the mathematical-angle doubled encodings used by the ML feature table.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

KEY = ["Eddy", "Day"]
QCOLS = ["q11", "q12", "q22"]
COLOURS = {"AE": "#b74646", "CE": "#326d9b"}


def axial_difference(a, b):
    """Signed axis difference in [-90, 90); invariant to 180-degree flips."""
    return (np.asarray(a) - np.asarray(b) + 90.0) % 180.0 - 90.0


def ellipse_geometry(frame, bearing_offset_deg=20.0):
    """Derive geometry only for finite, positive-definite Q matrices."""
    out = frame.copy()
    q = out[QCOLS].to_numpy(float)
    good = np.isfinite(q).all(axis=1)
    matrices = np.zeros((len(out), 2, 2))
    matrices[:, 0, 0], matrices[:, 0, 1] = q[:, 0], q[:, 1]
    matrices[:, 1, 0], matrices[:, 1, 1] = q[:, 1], q[:, 2]
    # Replace missing matrices only for eigensolver safety; keep them invalid.
    matrices[~good] = np.eye(2)
    vals, vecs = np.linalg.eigh(matrices)
    good &= (vals[:, 0] > 0) & np.isfinite(vals).all(axis=1)
    ratio = np.full(len(out), np.nan)
    ratio[good] = np.sqrt(vals[good, 1] / vals[good, 0])
    major = vecs[:, :, 0]
    bearing = (np.degrees(np.arctan2(major[:, 0], major[:, 1]))
               + bearing_offset_deg) % 180.0
    bearing[~good | np.isclose(ratio, 1.0, rtol=1e-8)] = np.nan
    out["AxisRatio"] = ratio
    out["MajorBearing"] = bearing
    out["Q_valid"] = good
    return out


def sample_depth_geometry(profiles, depths=(50, 100, 200, 300), *,
                          method="interpolate", max_gap_m=100.0,
                          nearest_tolerance_m=25.0):
    """Sample Q, not wrapped angles; never extrapolate or bridge invalid Q.

    Exact levels are preferred. Linear interpolation uses adjacent valid SPD
    matrices only, with a bounded bracket width. Nearest sampling is a
    sensitivity alternative, also restricted to the observed depth range.
    Duplicate eddy/day/depth rows fail explicitly rather than being averaged.
    """
    if method not in {"interpolate", "nearest"}:
        raise ValueError("method must be interpolate or nearest")
    depths = tuple(float(z) for z in depths)
    if len(set(depths)) != len(depths) or any(z <= 0 for z in depths):
        raise ValueError("Subsurface target depths must be unique and positive")
    data = profiles[KEY + ["Depth"] + QCOLS].copy()
    if data[KEY].isna().any().any():
        raise ValueError("Missing profile keys")
    if (data.Depth.dropna() < 0).any():
        raise ValueError("Expected positive-down Depth, as in the modular pipeline")
    if data.duplicated(KEY + ["Depth"]).any():
        raise ValueError("Duplicate Eddy/Day/Depth rows")
    data = ellipse_geometry(data)
    records = []
    for (eddy, day), g in data.groupby(KEY, sort=False):
        g = g.loc[np.isfinite(g.Depth)].sort_values("Depth")
        z = g.Depth.to_numpy(float)
        q = g[QCOLS].to_numpy(float)
        valid = g.Q_valid.to_numpy()
        for target in depths:
            if not len(z) or target < z[0] or target > z[-1]:
                continue
            j = int(np.searchsorted(z, target))
            if j < len(z) and np.isclose(z[j], target, atol=1e-6, rtol=0):
                lo = hi = j
                used_method = "exact"
            elif method == "nearest":
                lo = hi = int(np.argmin(np.abs(z - target)))
                if abs(z[lo] - target) > nearest_tolerance_m:
                    continue
                used_method = "nearest"
            else:
                lo, hi = j - 1, j
                if lo < 0 or hi >= len(z) or z[hi] - z[lo] > max_gap_m:
                    continue
                used_method = "interpolate"
            if not valid[lo] or not valid[hi]:
                continue
            weight = 0.0 if lo == hi else (target - z[lo]) / (z[hi] - z[lo])
            qq = (1 - weight) * q[lo] + weight * q[hi]
            records.append([eddy, day, target, *qq, z[lo], z[hi], used_method])
    return pd.DataFrame(records, columns=KEY + ["ShapeDepth"] + QCOLS
                        + ["DepthLower", "DepthUpper", "DepthMethod"])


def build_analysis_table(surface, sampled, bearing_offset_deg=20.0):
    if surface[KEY].isna().any().any() or surface.duplicated(KEY).any():
        raise ValueError("Surface Eddy/Day keys must be complete and unique")
    if sampled.duplicated(KEY + ["ShapeDepth"]).any():
        raise ValueError("Duplicate sampled geometry keys")
    base = surface.copy()
    if "Region" in base:
        base["Sector"] = base.Region.map({"S1": "Upstream", "U1": "Upstream",
            "U2": "Upstream", "S2": "Downstream", "D1": "Downstream", "D2": "Downstream"})
    if "topo_plan_ratio" in base:
        r = base.topo_plan_ratio
        base["PVRegime"] = np.select([r < -np.log(2), r > np.log(2), r.notna()],
            ["Planetary", "Topographic", "Mixed"], default="Missing")
    surf = base.assign(ShapeDepth=0.0, DepthLower=0.0, DepthUpper=0.0, DepthMethod="surface")
    deep = sampled.merge(base.drop(columns=QCOLS), on=KEY, how="inner", validate="many_to_one")
    combined = pd.concat([surf, deep], ignore_index=True) if len(deep) else surf.reset_index(drop=True)
    out = ellipse_geometry(combined, bearing_offset_deg)
    out["SignedAlignment"] = axial_difference(out.TiltDir, out.MajorBearing)
    out["AlignmentDeg"] = out.SignedAlignment.abs()
    out["AlignmentCos2"] = np.cos(np.deg2rad(2 * out.SignedAlignment))
    out["AlignmentSin2"] = np.sin(np.deg2rad(2 * out.SignedAlignment))
    return out.sort_values(["ShapeDepth", "Eddy", "Day"]).reset_index(drop=True)


def select_rows(table, *, directional=False, min_ar=1.1, min_tilt=5.0, max_ar=5.0):
    keep = (table.Q_valid & np.isfinite(table.AxisRatio) & table.AxisRatio.between(1, max_ar)
            & np.isfinite(table.TiltDis) & table.TiltDis.ge(0))
    if directional:
        keep &= (table.AxisRatio.ge(min_ar) & table.TiltDis.ge(min_tilt)
                 & np.isfinite(table.AlignmentDeg))
    return table.loc[keep].copy()


def matched_depth_sample(table, depths):
    """Identical eddy-days at every requested depth, after outcome-specific QC."""
    selected = table.loc[table.ShapeDepth.isin(depths)]
    counts = selected.groupby(KEY).ShapeDepth.nunique()
    keys = counts[counts == len(depths)].reset_index()[KEY]
    return selected.merge(keys, on=KEY, how="inner", validate="many_to_one")


def equal_eddy_weights(frame):
    return 1.0 / frame.groupby("Eddy").Eddy.transform("size")


def bootstrap_stat(values, statistic, n_boot=1000, seed=42, min_eddies=10):
    """Rows must be independent eddy-level units (or sufficient statistics)."""
    values = np.asarray(values, float)
    if len(values) < min_eddies:
        return np.nan, np.nan, np.nan
    estimate = float(statistic(values))
    rng = np.random.default_rng(seed)
    draws = np.array([statistic(values[rng.integers(len(values), size=len(values))])
                      for _ in range(n_boot)])
    finite = draws[np.isfinite(draws)]
    if len(finite) < 0.8 * n_boot:
        return estimate, np.nan, np.nan
    low, high = np.percentile(finite, [2.5, 97.5])
    return estimate, float(low), float(high)


def _rho(a):
    if np.ptp(a[:, 0]) == 0 or np.ptp(a[:, 1]) == 0:
        return np.nan
    return float(spearmanr(a[:, 0], a[:, 1]).statistic)


def magnitude_summary(table, groups=("ShapeDepth", "Cyc"), n_boot=1000,
                      seed=42, min_days=5, min_eddies=10):
    """Between: Spearman of eddy medians. Within: equal-eddy demeaned slope.

    Whole-eddy bootstrap preserves serial dependence. Within slope is km per
    unit AR; it is descriptive, not adjusted for environmental confounders.
    """
    records = []
    for key, g in table.groupby(list(groups), observed=True):
        key = key if isinstance(key, tuple) else (key,)
        g = g.replace([np.inf, -np.inf], np.nan).dropna(subset=["AxisRatio", "TiltDis"])
        g = g.loc[g.groupby("Eddy").Eddy.transform("size") >= min_days]
        med = g.groupby("Eddy")[["AxisRatio", "TiltDis"]].median().to_numpy()
        between = bootstrap_stat(med, _rho, n_boot, seed, min_eddies)
        dx = g.AxisRatio - g.groupby("Eddy").AxisRatio.transform("mean")
        dy = g.TiltDis - g.groupby("Eddy").TiltDis.transform("mean")
        moments = pd.DataFrame({"Eddy": g.Eddy, "xy": dx * dy, "xx": dx**2})
        moments = moments.groupby("Eddy")[["xy", "xx"]].mean()
        moments = moments.loc[moments.xx > 1e-12]
        within = bootstrap_stat(moments.to_numpy(), lambda a: a[:, 0].mean() / a[:, 1].mean(),
                                n_boot, seed, min_eddies)
        for metric, result, n in [("Between-eddy Spearman", between, len(med)),
                                   ("Within-eddy slope (km / AR)", within, len(moments))]:
            records.append(dict(zip(groups, key), metric=metric, estimate=result[0],
                                low=result[1], high=result[2], eddies=n, observations=len(g)))
    return pd.DataFrame(records, columns=[*groups, "metric", "estimate", "low", "high", "eddies", "observations"])


def alignment_summary(table, groups=("ShapeDepth", "Cyc"), n_boot=1000,
                      seed=42, min_days=5, min_eddies=10):
    """CI on mean of eddy means, including signed intermediate orientations."""
    records = []
    for key, g in table.groupby(list(groups), observed=True):
        key = key if isinstance(key, tuple) else (key,)
        g = g.loc[g.groupby("Eddy").Eddy.transform("size") >= min_days]
        values = g.groupby("Eddy")[["AlignmentCos2", "AlignmentSin2", "AlignmentDeg"]].mean()
        for col in values:
            a = values[[col]].dropna().to_numpy()
            est, low, high = bootstrap_stat(a, np.mean, n_boot, seed, min_eddies)
            records.append(dict(zip(groups, key), metric=col, estimate=est, low=low,
                                high=high, eddies=len(a), observations=len(g)))
    return pd.DataFrame(records, columns=[*groups, "metric", "estimate", "low", "high", "eddies", "observations"])


def rotation_pairs(table, max_gap_days=1.0):
    """Consecutive qualified observations only; axial rotation of tilt and Q."""
    out = table.sort_values(["ShapeDepth", "Eddy", "Day"]).copy()
    groups = out.groupby(["ShapeDepth", "Eddy"], sort=False)
    dt = groups.Day.diff()
    out["AxisTurn"] = axial_difference(out.MajorBearing, groups.MajorBearing.shift())
    out["TiltTurn"] = axial_difference(out.TiltDir, groups.TiltDir.shift())
    out["TurnAgreement"] = np.cos(np.deg2rad(2 * (out.TiltTurn - out.AxisTurn)))
    out["DeltaAR"] = groups.AxisRatio.diff()
    out["DeltaTilt"] = groups.TiltDis.diff()
    return out.loc[dt.gt(0) & dt.le(max_gap_days)].copy()


def plot_alignment(table, depth=0):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True, sharey=True)
    bins = np.linspace(0, 90, 19)
    for ax, cyc in zip(axes, COLOURS):
        g = table.loc[table.ShapeDepth.eq(depth) & table.Cyc.eq(cyc)]
        if len(g):
            hist, _ = np.histogram(g.AlignmentDeg, bins=bins, weights=equal_eddy_weights(g))
            ax.stairs(100 * hist / hist.sum(), bins, fill=True, alpha=.55, color=COLOURS[cyc])
        ax.axhline(100 / (len(bins)-1), ls="--", color="0.4", lw=1)
        ax.set(title=f"{cyc}: {g.Eddy.nunique():,} eddies", xlabel="Tilt–major-axis angle (°)", xlim=(0, 90))
    axes[0].set_ylabel("Eddy-equal probability per bin (%)")
    fig.suptitle(f"Shape at {depth:g} m: 0° parallel · 90° perpendicular")
    fig.tight_layout()
    return fig


def plot_magnitude(table, depth=0):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    for ax, cyc in zip(axes, COLOURS):
        g = table.loc[table.ShapeDepth.eq(depth) & table.Cyc.eq(cyc)]
        med = g.groupby("Eddy")[["AxisRatio", "TiltDis"]].median()
        ax.scatter(med.AxisRatio, med.TiltDis, s=10, alpha=.25, color=COLOURS[cyc], rasterized=True)
        if len(med) >= 10 and med.AxisRatio.nunique() > 1:
            med = med.assign(bin=pd.qcut(med.AxisRatio, 6, duplicates="drop"))
            b = med.groupby("bin", observed=True).agg(x=("AxisRatio", "median"),
                y=("TiltDis", "median"), lo=("TiltDis", lambda x: x.quantile(.25)),
                hi=("TiltDis", lambda x: x.quantile(.75)))
            ax.plot(b.x, b.y, "o-", color=COLOURS[cyc])
            ax.fill_between(b.x, b.lo, b.hi, color=COLOURS[cyc], alpha=.15)
        ax.set(title=f"{cyc}: one point per eddy", xlabel="Median major/minor axis ratio")
    axes[0].set_ylabel("Median tilt distance (km)")
    fig.suptitle(f"Shape at {depth:g} m · shading = descriptive IQR, not confidence interval")
    fig.tight_layout()
    return fig


def plot_interaction(table, ar_edges=(1.1, 1.3, 1.6, 2, 3, 5.000001)):
    """Equal eddy weighting separately within each AR class; no binwise CIs."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    classes = pd.IntervalIndex.from_breaks(ar_edges, closed="left")
    colours = dict(zip(classes.left, plt.get_cmap("viridis")(np.linspace(.1, .85, len(classes)))))
    for ax, cyc in zip(axes, COLOURS):
        g = table.loc[table.ShapeDepth.eq(0) & table.Cyc.eq(cyc)].copy()
        g["ARClass"] = pd.cut(g.AxisRatio, ar_edges, right=False)
        curves, labels, line_colours = [], [], []
        for label, part in g.groupby("ARClass", observed=True):
            if part.Eddy.nunique() < 10:
                continue
            hist, edges = np.histogram(part.AlignmentDeg, bins=np.linspace(0, 90, 10),
                                       weights=equal_eddy_weights(part))
            curves.append(100 * hist / hist.sum())
            labels.append(f"{label.left:g}–{label.right:g} (n={part.Eddy.nunique()})")
            line_colours.append(colours.get(label.left, "0.4"))
        for curve, label, colour in zip(curves, labels, line_colours):
            ax.plot((edges[:-1]+edges[1:])/2, curve, "o-", ms=3, label=label, color=colour)
        ax.set(title=cyc, xlabel="Tilt–major-axis angle (°)", xlim=(0,90))
        if curves:
            ax.legend(title="Axis ratio; eddies", fontsize=7)
    axes[0].set_ylabel("Eddy-equal probability per 10° (%)")
    fig.tight_layout()
    return fig


def plot_depth_summary(summary, metric):
    fig, ax = plt.subplots(figsize=(7, 4))
    for offset, (cyc, colour) in zip([-.8, .8], COLOURS.items()):
        g = summary.loc[summary.metric.eq(metric) & summary.Cyc.eq(cyc)].sort_values("ShapeDepth")
        valid = g.dropna(subset=["estimate", "low", "high"])
        ax.plot(valid.ShapeDepth + offset, valid.estimate, "o-", color=colour, label=cyc)
        ax.vlines(valid.ShapeDepth + offset, valid.low, valid.high, color=colour)
    ax.axhline(0, color="0.5", ls="--", lw=1)
    ylabel = {"AlignmentCos2": "Mean cos(2 × tilt–axis offset)",
              "AlignmentSin2": "Mean sin(2 × tilt–axis offset)"}.get(metric, metric)
    ax.set(xlabel="Ellipse sampling depth (m)", ylabel=ylabel,
           title="Eddy-bootstrap 95% intervals (pointwise)")
    ax.legend()
    fig.tight_layout()
    return fig
