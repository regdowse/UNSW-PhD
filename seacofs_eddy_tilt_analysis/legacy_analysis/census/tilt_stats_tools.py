"""Descriptive statistics used by the tilt-census notebooks and manuscript.

The functions in this module deliberately keep the daily observations as the
population being described. Per-eddy statistics are calculated separately so
that the two sampling units are never mixed silently.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd


UPSTREAM_REGIONS = ("S1", "U1", "U2")
DOWNSTREAM_REGIONS = ("S2", "D1", "D2")
OFFSHORE_REGIONS = ("U1", "U2", "D1", "D2")
SHELF_REGIONS = ("S1", "S2")
OFFSHORE_UPSTREAM_REGIONS = ("U1", "U2")
OFFSHORE_DOWNSTREAM_REGIONS = ("D1", "D2")


def _require_columns(data: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = sorted(set(columns).difference(data.columns))
    if missing:
        raise KeyError(f"Missing required columns: {', '.join(missing)}")


def _safe_pct(numerator: float, denominator: float) -> float:
    return 100.0 * numerator / denominator if denominator else np.nan


def _fmt(value: float, digits: int = 1) -> str:
    return "NA" if pd.isna(value) else f"{value:.{digits}f}"


def prepare_tilt_data(
    data: pd.DataFrame,
    *,
    first_last_days: int = 3,
    upstream_regions: Sequence[str] = UPSTREAM_REGIONS,
    downstream_regions: Sequence[str] = DOWNSTREAM_REGIONS,
) -> pd.DataFrame:
    """Return a census-ready copy with eligibility and broad-region labels."""
    _require_columns(data, ["Cyc", "Eddy", "Day", "TiltDis", "Rc"])
    out = data.copy()
    keys = ["Cyc", "Eddy"]
    track_index = out.groupby(keys).cumcount()
    track_length = out.groupby(keys)["Day"].transform("size")
    out["tilt_possible"] = (
        track_index.ge(first_last_days)
        & track_index.lt(track_length - first_last_days)
    )
    out["valid_tilt"] = out["TiltDis"].notna()
    valid_rc = np.isfinite(out["Rc"]) & out["Rc"].gt(0)
    out["tilt_over_Rc"] = np.where(valid_rc, out["TiltDis"] / out["Rc"], np.nan)
    if "Region" in out:
        out["Sector"] = np.select(
            [out["Region"].isin(upstream_regions), out["Region"].isin(downstream_regions)],
            ["Upstream", "Downstream"],
            default=None,
        )
    return out


def _basic_stats(values: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "n": len(x),
        "mean_km": x.mean(),
        "median_km": x.median(),
        "q25_km": x.quantile(0.25),
        "q75_km": x.quantile(0.75),
        "p95_km": x.quantile(0.95),
        "maximum_km": x.max(),
        "over_40_pct": 100 * x.gt(40).mean() if len(x) else np.nan,
        "over_50_pct": 100 * x.gt(50).mean() if len(x) else np.nan,
        "over_100_pct": 100 * x.gt(100).mean() if len(x) else np.nan,
    }


def _ratio_stats(values: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "n": len(x),
        "mean_tilt_over_Rc": x.mean(),
        "median_tilt_over_Rc": x.median(),
        "q25_tilt_over_Rc": x.quantile(0.25),
        "q75_tilt_over_Rc": x.quantile(0.75),
        "p95_tilt_over_Rc": x.quantile(0.95),
        "maximum_tilt_over_Rc": x.max(),
        "over_0.5_Rc_pct": 100 * x.gt(0.5).mean() if len(x) else np.nan,
        "over_1_Rc_pct": 100 * x.gt(1).mean() if len(x) else np.nan,
    }


def _modal_interval(values: pd.Series, width: float = 10.0) -> tuple[float, float]:
    x = pd.to_numeric(values, errors="coerce").dropna()
    if x.empty:
        return np.nan, np.nan
    edges = np.arange(0, max(width * 2, np.ceil(x.max() / width) * width + width), width)
    counts, edges = np.histogram(x, bins=edges)
    peak = int(np.argmax(counts))
    return float(edges[peak]), float(edges[peak + 1])


def _grouped_basic_stats(data: pd.DataFrame, groups: Sequence[str], value: str) -> pd.DataFrame:
    rows = []
    for names, group in data.groupby(list(groups), observed=True, dropna=False):
        names = names if isinstance(names, tuple) else (names,)
        rows.append({**dict(zip(groups, names)), **_basic_stats(group[value])})
    return pd.DataFrame(rows).set_index(list(groups)).sort_index()


def _grouped_ratio_stats(data: pd.DataFrame, groups: Sequence[str], value: str) -> pd.DataFrame:
    rows = []
    for names, group in data.groupby(list(groups), observed=True, dropna=False):
        names = names if isinstance(names, tuple) else (names,)
        rows.append({**dict(zip(groups, names)), **_ratio_stats(group[value])})
    return pd.DataFrame(rows).set_index(list(groups)).sort_index()


def tilt_census_report(
    data: pd.DataFrame,
    *,
    first_last_days: int = 2,
) -> dict[str, object]:
    """Calculate daily and per-eddy tilt-distance census statistics."""
    prepared = prepare_tilt_data(data, first_last_days=first_last_days)
    valid = prepared.loc[prepared["valid_tilt"]].copy()
    keys = ["Cyc", "Eddy"]

    census_rows = []
    for cyc in ("All", "AE", "CE"):
        group = prepared if cyc == "All" else prepared.loc[prepared["Cyc"].eq(cyc)]
        possible = group["tilt_possible"].sum()
        obtained = (group["tilt_possible"] & group["valid_tilt"]).sum()
        census_rows.append({
            "Cyc": cyc,
            "unique_eddies": group[keys].drop_duplicates().shape[0],
            "eddy_days": len(group),
            "possible_tilt_days": possible,
            "valid_tilt_days": group["valid_tilt"].sum(),
            "raw_coverage_pct": _safe_pct(group["valid_tilt"].sum(), len(group)),
            "valid_possible_tilt_days": obtained,
            "possible_coverage_pct": _safe_pct(obtained, possible),
        })
    census = pd.DataFrame(census_rows).set_index("Cyc").reindex(["All", "AE", "CE"])

    daily = _grouped_basic_stats(valid, ["Cyc"], "TiltDis")
    daily = pd.concat(
        [pd.DataFrame([_basic_stats(valid["TiltDis"])], index=pd.Index(["All"], name="Cyc")), daily]
    )
    daily_ratio = _grouped_ratio_stats(valid, ["Cyc"], "tilt_over_Rc")
    daily_ratio = pd.concat([
        pd.DataFrame([_ratio_stats(valid["tilt_over_Rc"])], index=pd.Index(["All"], name="Cyc")),
        daily_ratio,
    ])
    modal_rows = []
    for cyc in ("All", "AE", "CE"):
        group = valid if cyc == "All" else valid.loc[valid["Cyc"].eq(cyc)]
        lower, upper = _modal_interval(group["TiltDis"], 10.0)
        modal_rows.append({"Cyc": cyc, "modal_bin_lower_km": lower, "modal_bin_upper_km": upper})
    modal_bins = pd.DataFrame(modal_rows).set_index("Cyc")
    eddy_summary = (
        valid.groupby(keys, observed=True)
        .agg(
            mean_tilt=("TiltDis", "mean"), median_tilt=("TiltDis", "median"),
            min_tilt=("TiltDis", "min"), max_tilt=("TiltDis", "max"), n_tilt=("TiltDis", "size"),
            mean_tilt_over_Rc=("tilt_over_Rc", "mean"),
            median_tilt_over_Rc=("tilt_over_Rc", "median"),
            max_tilt_over_Rc=("tilt_over_Rc", "max"),
        )
        .reset_index()
    )
    per_eddy = pd.concat(
        {
            metric: _grouped_basic_stats(eddy_summary, ["Cyc"], metric)
            for metric in ("mean_tilt", "median_tilt", "min_tilt", "max_tilt")
        },
        names=["metric"],
    )
    per_eddy_ratio = pd.concat(
        {
            metric: _grouped_ratio_stats(eddy_summary, ["Cyc"], metric)
            for metric in ("mean_tilt_over_Rc", "median_tilt_over_Rc", "max_tilt_over_Rc")
        },
        names=["metric"],
    )
    lifetime_events = (
        eddy_summary.assign(
            ever_over_50_km=eddy_summary["max_tilt"].gt(50),
            ever_over_100_km=eddy_summary["max_tilt"].gt(100),
            ever_over_1_Rc=eddy_summary["max_tilt_over_Rc"].gt(1),
        )
        .groupby("Cyc", observed=True)
        .agg(
            n_eddies=("Eddy", "size"),
            median_mean_tilt_km=("mean_tilt", "median"),
            median_mean_tilt_over_Rc=("mean_tilt_over_Rc", "median"),
            ever_over_50_km_pct=("ever_over_50_km", lambda x: 100 * x.mean()),
            ever_over_100_km_pct=("ever_over_100_km", lambda x: 100 * x.mean()),
            ever_over_1_Rc_pct=("ever_over_1_Rc", lambda x: 100 * x.mean()),
        )
    )

    sector = (
        _grouped_basic_stats(valid.dropna(subset=["Sector"]), ["Cyc", "Sector"], "TiltDis")
        if "Sector" in valid else pd.DataFrame()
    )
    sector_ratio = (
        _grouped_ratio_stats(valid.dropna(subset=["Sector"]), ["Cyc", "Sector"], "tilt_over_Rc")
        if "Sector" in valid else pd.DataFrame()
    )
    regional = (
        _grouped_basic_stats(valid.dropna(subset=["Region"]), ["Cyc", "Region"], "TiltDis")
        if "Region" in valid else pd.DataFrame()
    )
    regional_ratio = (
        _grouped_ratio_stats(valid.dropna(subset=["Region"]), ["Cyc", "Region"], "tilt_over_Rc")
        if "Region" in valid else pd.DataFrame()
    )
    latitude = pd.DataFrame()
    if "lat" in valid:
        lat_data = valid.assign(latitude_band=pd.cut(valid["lat"], np.arange(-40, -25, 1)))
        latitude = _grouped_basic_stats(
            lat_data.dropna(subset=["latitude_band"]), ["Cyc", "latitude_band"], "TiltDis"
        )

    return {
        "data": prepared,
        "valid": valid,
        "census": census,
        "daily": daily,
        "daily_ratio": daily_ratio,
        "modal_bins": modal_bins,
        "eddy_summary": eddy_summary,
        "per_eddy": per_eddy,
        "per_eddy_ratio": per_eddy_ratio,
        "lifetime_events": lifetime_events,
        "sector": sector,
        "sector_ratio": sector_ratio,
        "regional": regional,
        "regional_ratio": regional_ratio,
        "latitude": latitude,
    }


def tilt_census_sentences(report: dict[str, object]) -> list[str]:
    """Return copyable statements covering the manuscript's census claims."""
    census = report["census"]
    daily = report["daily"]
    ratio = report["daily_ratio"]
    modal = report["modal_bins"]
    lifetime = report["lifetime_events"]
    all_row, ae_count, ce_count = census.loc["All"], census.loc["AE"], census.loc["CE"]
    lines = [
        f"The dataset contained {int(all_row.unique_eddies):,} unique eddies "
        f"({int(ae_count.unique_eddies):,} AEs and {int(ce_count.unique_eddies):,} CEs) and "
        f"{int(all_row.eddy_days):,} eddy-days ({int(ae_count.eddy_days):,} AE-days and "
        f"{int(ce_count.eddy_days):,} CE-days).",
        f"The Delta method produced {int(all_row.valid_tilt_days):,} valid tilt estimates, "
        f"corresponding to {_fmt(all_row.raw_coverage_pct)}% of all eddy-days.",
        f"Valid tilt estimates were obtained for {int(ae_count.valid_tilt_days):,} AE-days and "
        f"{int(ce_count.valid_tilt_days):,} CE-days, covering {_fmt(ae_count.possible_coverage_pct)}% and "
        f"{_fmt(ce_count.possible_coverage_pct)}% of eligible days, respectively, after excluding the first and last two days of each track.",
    ]
    for cyc in ("AE", "CE"):
        s, r, m, life = daily.loc[cyc], ratio.loc[cyc], modal.loc[cyc], lifetime.loc[cyc]
        lines.extend([
            f"The modal 10-km {cyc} tilt-distance interval was {m.modal_bin_lower_km:g}–{m.modal_bin_upper_km:g} km.",
            f"The median {cyc} tilt distance was {_fmt(s.median_km)} km "
            f"(IQR {_fmt(s.q25_km)}–{_fmt(s.q75_km)} km; mean {_fmt(s.mean_km)} km).",
            f"{cyc} tilt exceeded 50 km on {_fmt(s.over_50_pct)}% of eddy-days and 100 km on {_fmt(s.over_100_pct)}% of eddy-days.",
            f"The median normalised {cyc} tilt was {_fmt(r.median_tilt_over_Rc, 2)}Rc "
            f"(IQR {_fmt(r.q25_tilt_over_Rc, 2)}–{_fmt(r.q75_tilt_over_Rc, 2)}Rc), and "
            f"tilt exceeded one core radius on {_fmt(r['over_1_Rc_pct'])}% of eddy-days.",
            f"The median of per-eddy mean tilt was {_fmt(life.median_mean_tilt_km)} km "
            f"({_fmt(life.median_mean_tilt_over_Rc, 2)}Rc).",
            f"Across individual {cyc}s, {_fmt(life.ever_over_50_km_pct)}% experienced at least one tilt above 50 km, "
            f"{_fmt(life.ever_over_100_km_pct)}% exceeded 100 km, and {_fmt(life.ever_over_1_Rc_pct)}% exceeded one core radius.",
        ])
    return lines


def circular_stats_deg(direction: pd.Series) -> tuple[float, float]:
    """Circular mean bearing and mean-resultant length for degrees in [0, 360)."""
    theta = np.deg2rad(pd.to_numeric(direction, errors="coerce").dropna().to_numpy())
    if not len(theta):
        return np.nan, np.nan
    x, y = np.cos(theta).mean(), np.sin(theta).mean()
    return float(np.rad2deg(np.arctan2(y, x)) % 360), float(np.hypot(x, y))


def prepare_direction_data(
    data: pd.DataFrame,
    *,
    minimum_tilt_km: float = 5.0,
    offshore_regions: Sequence[str] = OFFSHORE_REGIONS,
) -> pd.DataFrame:
    """Filter valid directions and add compass, polarity and offshore labels."""
    _require_columns(data, ["Cyc", "Eddy", "TiltDis", "TiltDir", "Region"])
    out = data.loc[
        data["TiltDis"].ge(minimum_tilt_km)
        & data["TiltDir"].notna()
        & data["Cyc"].isin(["AE", "CE"])
    ].copy()
    out["TiltDir"] %= 360
    if "Rc" in out:
        valid_rc = np.isfinite(out["Rc"]) & out["Rc"].gt(0)
        out["tilt_over_Rc"] = np.where(valid_rc, out["TiltDis"] / out["Rc"], np.nan)
    else:
        out["tilt_over_Rc"] = np.nan
    angle = out["TiltDir"]
    out["northward"] = angle.lt(90) | angle.ge(270)
    out["southward"] = angle.ge(90) & angle.lt(270)
    out["eastward"] = angle.lt(180)
    out["westward"] = angle.ge(180)
    out["NE"] = angle.lt(90)
    out["SE"] = angle.ge(90) & angle.lt(180)
    out["SW"] = angle.ge(180) & angle.lt(270)
    out["NW"] = angle.ge(270)
    out["is_offshore"] = out["Region"].isin(offshore_regions)
    out["is_shelf"] = out["Region"].isin(SHELF_REGIONS)
    out["OffshoreSector"] = np.select(
        [
            out["Region"].isin(OFFSHORE_UPSTREAM_REGIONS),
            out["Region"].isin(OFFSHORE_DOWNSTREAM_REGIONS),
        ],
        ["Upstream offshore", "Downstream offshore"],
        default=None,
    )
    return out


def _direction_stats(group: pd.DataFrame, cyc: str) -> dict[str, float]:
    mean_direction, concentration = circular_stats_deg(group["TiltDir"])
    equatorward = group["northward"] if cyc == "AE" else group["southward"]
    preferred_quadrant = group["NE"] if cyc == "AE" else group["SE"]
    return {
        "n": len(group),
        "n_eddies": group["Eddy"].nunique(),
        "mean_direction_deg": mean_direction,
        "directional_concentration_R": concentration,
        "equatorward_AE_or_poleward_CE_pct": 100 * equatorward.mean(),
        "eastward_pct": 100 * group["eastward"].mean(),
        "westward_pct": 100 * group["westward"].mean(),
        "preferred_quadrant_NE_AE_or_SE_CE_pct": 100 * preferred_quadrant.mean(),
        "NE_pct": 100 * group["NE"].mean(),
        "SE_pct": 100 * group["SE"].mean(),
        "SW_pct": 100 * group["SW"].mean(),
        "NW_pct": 100 * group["NW"].mean(),
        "median_tilt_km": group["TiltDis"].median(),
        "median_tilt_over_Rc": group["tilt_over_Rc"].median(),
        "over_40_km_pct": 100 * group["TiltDis"].gt(40).mean(),
        "over_1_Rc_pct": 100 * group["tilt_over_Rc"].gt(1).mean(),
    }


def _direction_table(data: pd.DataFrame, groups: Sequence[str]) -> pd.DataFrame:
    rows = []
    for names, group in data.groupby(list(groups), observed=True, dropna=False):
        names = names if isinstance(names, tuple) else (names,)
        labels = dict(zip(groups, names))
        rows.append({**labels, **_direction_stats(group, str(labels["Cyc"]))})
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).set_index(list(groups)).sort_index()


def direction_census_report(
    data: pd.DataFrame,
    *,
    minimum_tilt_km: float = 5.0,
    offshore_regions: Sequence[str] = OFFSHORE_REGIONS,
) -> dict[str, object]:
    """Summarise direction for all, shelf and offshore eddies.

    ``offshore`` pools U1, U2, D1 and D2. ``offshore_sector`` contrasts
    U1/U2 with D1/D2 without allowing S1 or S2 into either group.
    """
    prepared = prepare_direction_data(
        data, minimum_tilt_km=minimum_tilt_km, offshore_regions=offshore_regions
    )
    offshore = prepared.loc[prepared["is_offshore"]].copy()
    shelf = prepared.loc[prepared["is_shelf"]].copy()
    return {
        "data": prepared,
        "offshore_data": offshore,
        "shelf_data": shelf,
        "overall": _direction_table(prepared, ["Cyc"]),
        "shelf": _direction_table(shelf, ["Cyc"]),
        "offshore": _direction_table(offshore, ["Cyc"]),
        "regional": _direction_table(prepared.dropna(subset=["Region"]), ["Cyc", "Region"]),
        "shelf_regional": _direction_table(shelf, ["Cyc", "Region"]),
        "offshore_regional": _direction_table(offshore, ["Cyc", "Region"]),
        "offshore_sector": _direction_table(
            offshore.dropna(subset=["OffshoreSector"]), ["Cyc", "OffshoreSector"]
        ),
        "minimum_tilt_km": minimum_tilt_km,
        "offshore_regions": tuple(offshore_regions),
    }


def direction_census_sentences(report: dict[str, object]) -> list[str]:
    """Return copyable total, shelf and offshore direction statements."""
    lines = [
        f"Direction statistics exclude tilts below {report['minimum_tilt_km']:g} km, for which direction is poorly defined."
    ]
    populations = (
        ("All regions", report["overall"]),
        ("Shelf regions (S1 and S2)", report["shelf"]),
        ("Offshore regions (U1, U2, D1 and D2)", report["offshore"]),
    )
    for label, table in populations:
        for cyc in ("AE", "CE"):
            if cyc not in table.index:
                continue
            row = table.loc[cyc]
            preferred = "equatorward" if cyc == "AE" else "poleward"
            lines.append(
                f"{label} — {cyc}: {int(row.n):,} eddy-days from {int(row.n_eddies):,} eddies; "
                f"{_fmt(row.equatorward_AE_or_poleward_CE_pct)}% tilted {preferred}, "
                f"{_fmt(row.westward_pct)}% had a westward component; "
                f"circular mean {_fmt(row.mean_direction_deg)}° (R = {_fmt(row.directional_concentration_R, 2)})."
            )
    for scope, table in (("Shelf", report["shelf_regional"]), ("Offshore", report["offshore_regional"])):
        for (cyc, region), row in table.iterrows():
            preferred = "equatorward" if cyc == "AE" else "poleward"
            lines.append(
                f"{scope} region {region} — {cyc}: {_fmt(row.equatorward_AE_or_poleward_CE_pct)}% tilted {preferred}, "
                f"{_fmt(row.westward_pct)}% had a westward component, median tilt was "
                f"{_fmt(row.median_tilt_km)} km ({_fmt(row.median_tilt_over_Rc, 2)}Rc), and the circular mean was "
                f"{_fmt(row.mean_direction_deg)}° (R = {_fmt(row.directional_concentration_R, 2)})."
            )
    return lines


def angle_difference_deg(a: pd.Series, b: pd.Series) -> np.ndarray:
    """Signed smallest angular difference a - b in [-180, 180)."""
    return (pd.to_numeric(a, errors="coerce") - pd.to_numeric(b, errors="coerce") + 180) % 360 - 180


def _pv_table(data: pd.DataFrame, groups: Sequence[str]) -> pd.DataFrame:
    rows = []
    for names, group in data.groupby(list(groups), observed=True, dropna=False):
        names = names if isinstance(names, tuple) else (names,)
        diff = group["tilt_minus_expected_deg"].abs()
        rows.append({
            **dict(zip(groups, names)),
            "n": len(group),
            "median_absolute_misalignment_deg": diff.median(),
            "within_30_deg_pct": 100 * diff.le(30).mean(),
            "within_45_deg_pct": 100 * diff.le(45).mean(),
            "within_90_deg_pct": 100 * diff.le(90).mean(),
        })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).set_index(list(groups)).sort_index()


def pv_alignment_report(
    direction_data: pd.DataFrame,
    *,
    pv_direction_col: str = "PV_grad_theta",
) -> dict[str, pd.DataFrame]:
    """Test AE anti-alignment and CE alignment with the PV-gradient direction."""
    _require_columns(direction_data, ["Cyc", "TiltDir", "Region", pv_direction_col])
    data = direction_data.dropna(subset=[pv_direction_col]).copy()
    expected = np.where(data["Cyc"].eq("AE"), data[pv_direction_col] + 180, data[pv_direction_col]) % 360
    data["tilt_minus_expected_deg"] = angle_difference_deg(data["TiltDir"], expected)
    offshore = data.loc[data["Region"].isin(OFFSHORE_REGIONS)].copy()
    shelf = data.loc[data["Region"].isin(SHELF_REGIONS)].copy()
    return {
        "data": data,
        "overall": _pv_table(data, ["Cyc"]),
        "shelf": _pv_table(shelf, ["Cyc"]),
        "offshore": _pv_table(offshore, ["Cyc"]),
        "regional": _pv_table(data, ["Cyc", "Region"]),
        "shelf_regional": _pv_table(shelf, ["Cyc", "Region"]),
        "offshore_regional": _pv_table(offshore, ["Cyc", "Region"]),
    }
