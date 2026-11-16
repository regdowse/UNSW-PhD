"""Utilities for comparing PV-gradient proxies along depth-resolved eddy spines.

The cached quantity is the shallow-water/topographic proxy grad[(f + zeta) / h]
sampled with the eddy ellipse at each depth.  It is not local three-dimensional
Ertel PV.  TiltDis and TiltDir remain the authoritative whole-column tilt.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


KEYS = ["Eddy", "Day"]
DEPTH_KEYS = ["Eddy", "Day", "Depth"]
PV_COLUMNS = [
    "PV_grad_plan_mag", "PV_grad_topo_mag", "PV_grad_mag",
    "PV_grad_plan_theta", "PV_grad_topo_theta", "PV_grad_theta",
    "topo_plan_ratio", "h", "dhdx", "dhdy", "Ro", "w", "PV",
]


def angle_difference(a, b):
    """Signed a-minus-b angular difference in [-180, 180) degrees."""
    return (np.asarray(a, dtype=float) - np.asarray(b, dtype=float) + 180.0) % 360.0 - 180.0


def validate_depth_tables(depth: pd.DataFrame, snapshot: pd.DataFrame | None = None):
    """Validate cache keys and the columns used throughout the notebooks."""
    required = set(DEPTH_KEYS + ["Cyc", "TiltDis", "TiltDir", "xc", "yc"] + PV_COLUMNS)
    missing = required - set(depth.columns)
    if missing:
        raise KeyError(f"Depth cache is missing {sorted(missing)}")
    if depth.duplicated(DEPTH_KEYS).any():
        raise ValueError("Depth cache contains duplicate Eddy-Day-Depth rows")
    if snapshot is not None:
        missing = set(KEYS + ["Cyc", "TiltDis", "TiltDir"] + PV_COLUMNS) - set(snapshot.columns)
        if missing:
            raise KeyError(f"Snapshot cache is missing {sorted(missing)}")
        if snapshot.duplicated(KEYS).any():
            raise ValueError("Snapshot cache contains duplicate Eddy-Day rows")
    return True


def nearest_cached_depths(depth: pd.DataFrame, targets=(0, 200, 500, 700, 1000)) -> np.ndarray:
    """Return unique cached levels nearest requested target depths."""
    available = np.sort(pd.to_numeric(depth["Depth"], errors="coerce").dropna().unique())
    if not len(available):
        raise ValueError("No finite depth levels are available")
    chosen = [available[np.abs(available - float(target)).argmin()] for target in targets]
    return np.asarray(list(dict.fromkeys(chosen)), dtype=float)


def surface_rows(depth: pd.DataFrame) -> pd.DataFrame:
    """Select the shallowest available cached level independently per snapshot."""
    ordered = depth.sort_values(DEPTH_KEYS)
    return ordered.groupby(KEYS, sort=False, as_index=False).first()


def add_regime(df: pd.DataFrame, dominance_factor=2.0) -> pd.DataFrame:
    """Classify planetary, mixed, and topographic regimes from the log ratio."""
    if dominance_factor <= 1:
        raise ValueError("dominance_factor must exceed one")
    out = df.copy()
    threshold = np.log(float(dominance_factor))
    ratio = pd.to_numeric(out["topo_plan_ratio"], errors="coerce")
    out["regime"] = pd.Categorical(
        np.select(
            [ratio <= -threshold, ratio >= threshold],
            ["Planetary", "Topographic"], default="Mixed",
        ), categories=["Planetary", "Mixed", "Topographic"], ordered=True,
    )
    out.loc[~np.isfinite(ratio), "regime"] = np.nan
    return out


def add_expected_alignment(df: pd.DataFrame, direction_col="PV_grad_theta") -> pd.DataFrame:
    """Add polarity-aware directional agreement with the signed PV gradient.

    AEs are expected along signed grad(PV), while CEs are expected opposite it.
    ``preference_error_deg`` is zero for the expected response, and
    ``preference_alignment`` is +1 for perfect agreement.
    """
    out = df.copy()
    raw = np.abs(angle_difference(out["TiltDir"], out[direction_col]))
    out["preference_error_deg"] = np.where(out["Cyc"].eq("CE"), np.abs(180.0 - raw), raw)
    out["preference_alignment"] = np.cos(np.deg2rad(out["preference_error_deg"]))
    return out


def add_surface_differences(depth: pd.DataFrame, dominance_factor=2.0) -> pd.DataFrame:
    """Add paired changes from each snapshot's shallowest cached level."""
    validate_depth_tables(depth)
    out = add_regime(depth, dominance_factor)
    surface = add_regime(surface_rows(depth), dominance_factor)
    reference = ["Depth", "xc", "yc", "regime"] + PV_COLUMNS
    surface = surface[KEYS + reference].rename(columns={c: f"surface_{c}" for c in reference})
    out = out.merge(surface, on=KEYS, how="left", validate="many_to_one")
    out["spine_displacement_km"] = np.hypot(
        out["xc"] - out["surface_xc"], out["yc"] - out["surface_yc"]
    )
    if "Rc" in out:
        out["spine_displacement_Rc"] = out["spine_displacement_km"] / out["Rc"]
    out["delta_h_m"] = out["h"] - out["surface_h"]
    out["abs_delta_h_m"] = out["delta_h_m"].abs()
    out["delta_topo_plan_ratio"] = out["topo_plan_ratio"] - out["surface_topo_plan_ratio"]
    for direction in ("PV_grad_theta", "PV_grad_topo_theta", "PV_grad_plan_theta"):
        out[f"surface_rotation_{direction}_deg"] = np.abs(
            angle_difference(out[direction], out[f"surface_{direction}"])
        )
    out["regime_changed"] = out["regime"].notna() & out["surface_regime"].notna() & out["regime"].ne(out["surface_regime"])
    return add_expected_alignment(out)


def matched_depth_rows(df: pd.DataFrame, depths) -> pd.DataFrame:
    """Restrict to Eddy-Day snapshots represented at every requested exact level."""
    depths = np.asarray(depths, dtype=float)
    use = df[df["Depth"].isin(depths)].copy()
    counts = use.groupby(KEYS)["Depth"].nunique()
    complete = counts[counts.eq(len(depths))].index
    return use.set_index(KEYS).loc[complete].reset_index()


def eddy_equal_depth_summary(df: pd.DataFrame, value: str) -> pd.DataFrame:
    """Depth summaries in which each eddy contributes one median value."""
    eddy = df.groupby(["Cyc", "Depth", "Eddy"], observed=True)[value].median().reset_index()
    return (eddy.groupby(["Cyc", "Depth"], observed=True)[value]
            .agg(median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75),
                 eddies="count").reset_index())


def directional_scorecard(df: pd.DataFrame, group=("Cyc", "Depth")) -> pd.DataFrame:
    """Eddy-equal polarity-aware directional metrics by group."""
    required = set(group) | {"Eddy", "preference_error_deg", "preference_alignment"}
    if missing := required - set(df.columns):
        raise KeyError(f"Missing scorecard columns: {sorted(missing)}")
    eddy = (df.groupby([*group, "Eddy"], observed=True)
            .agg(preference_error_deg=("preference_error_deg", "median"),
                 preference_alignment=("preference_alignment", "mean"))
            .reset_index())
    return (eddy.groupby(list(group), observed=True)
            .agg(eddies=("Eddy", "nunique"),
                 median_error_deg=("preference_error_deg", "median"),
                 q25_error_deg=("preference_error_deg", lambda x: x.quantile(.25)),
                 q75_error_deg=("preference_error_deg", lambda x: x.quantile(.75)),
                 mean_alignment=("preference_alignment", "mean"),
                 fraction_within_30=("preference_error_deg", lambda x: np.mean(x <= 30)),
                 fraction_within_45=("preference_error_deg", lambda x: np.mean(x <= 45)))
            .reset_index())


def regime_transition_matrix(df: pd.DataFrame, depth: float, cyc: str | None = None) -> pd.DataFrame:
    """Row-normalised surface-to-depth regime transition matrix."""
    part = df[np.isclose(df["Depth"], depth)].copy()
    if cyc is not None:
        part = part[part["Cyc"].eq(cyc)]
    labels = ["Planetary", "Mixed", "Topographic"]
    matrix = pd.crosstab(part["surface_regime"], part["regime"], normalize="index")
    return matrix.reindex(index=labels, columns=labels, fill_value=0.0)


def representation_scorecard(depth: pd.DataFrame, snapshot: pd.DataFrame, depths) -> pd.DataFrame:
    """Compare direction metrics for fixed levels and the column-vector average."""
    rows = []
    selected = matched_depth_rows(add_expected_alignment(depth), depths)
    for z, part in selected.groupby("Depth"):
        score = directional_scorecard(part).copy()
        score["representation"] = f"{z:g} m"
        rows.append(score)
    snap = add_expected_alignment(snapshot).copy()
    score = directional_scorecard(snap.assign(Depth=-1)).drop(columns="Depth")
    score["representation"] = "0–1000 m vector mean"
    rows.append(score)
    return pd.concat(rows, ignore_index=True)


def rank_mismatch_cases(df: pd.DataFrame, min_tilt_km=5.0) -> pd.DataFrame:
    """Rank snapshots where depth reveals topography missed at the surface."""
    use = df[
        df["surface_regime"].eq("Planetary")
        & df["regime"].eq("Topographic")
        & df["TiltDis"].ge(min_tilt_km)
    ].copy()
    use["mismatch_score"] = (
        use["delta_topo_plan_ratio"].clip(lower=0)
        * np.log1p(use["spine_displacement_km"].clip(lower=0))
        * np.log1p(use["TiltDis"].clip(lower=0))
    )
    return use.sort_values("mismatch_score", ascending=False)
