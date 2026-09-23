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
    _require_columns(data, ["Cyc", "Eddy", "Day", "TiltDis"])
    out = data.copy()
    keys = ["Cyc", "Eddy"]
    track_index = out.groupby(keys).cumcount()
    track_length = out.groupby(keys)["Day"].transform("size")
    out["tilt_possible"] = (
        track_index.ge(first_last_days)
        & track_index.lt(track_length - first_last_days)
    )
    out["valid_tilt"] = out["TiltDis"].notna()
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
        "over_100_pct": 100 * x.gt(100).mean() if len(x) else np.nan,
    }


def _grouped_basic_stats(data: pd.DataFrame, groups: Sequence[str], value: str) -> pd.DataFrame:
    rows = []
    for names, group in data.groupby(list(groups), observed=True, dropna=False):
        names = names if isinstance(names, tuple) else (names,)
        rows.append({**dict(zip(groups, names)), **_basic_stats(group[value])})
    return pd.DataFrame(rows).set_index(list(groups)).sort_index()


def tilt_census_report(
    data: pd.DataFrame,
    *,
    first_last_days: int = 3,
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
            "valid_possible_tilt_days": obtained,
            "possible_coverage_pct": _safe_pct(obtained, possible),
        })
    census = pd.DataFrame(census_rows).set_index("Cyc").reindex(["All", "AE", "CE"])

    daily = _grouped_basic_stats(valid, ["Cyc"], "TiltDis")
    daily = pd.concat(
        [pd.DataFrame([_basic_stats(valid["TiltDis"])], index=pd.Index(["All"], name="Cyc")), daily]
    )
    eddy_summary = (
        valid.groupby(keys, observed=True)["TiltDis"]
        .agg(mean_tilt="mean", median_tilt="median", min_tilt="min", max_tilt="max", n_tilt="size")
        .reset_index()
    )
    per_eddy = pd.concat(
        {
            metric: _grouped_basic_stats(eddy_summary, ["Cyc"], metric)
            for metric in ("mean_tilt", "median_tilt", "min_tilt", "max_tilt")
        },
        names=["metric"],
    )

    sector = (
        _grouped_basic_stats(valid.dropna(subset=["Sector"]), ["Cyc", "Sector"], "TiltDis")
        if "Sector" in valid else pd.DataFrame()
    )
    regional = (
        _grouped_basic_stats(valid.dropna(subset=["Region"]), ["Cyc", "Region"], "TiltDis")
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
        "eddy_summary": eddy_summary,
        "per_eddy": per_eddy,
        "sector": sector,
        "regional": regional,
        "latitude": latitude,
    }


def tilt_census_sentences(report: dict[str, object]) -> list[str]:
    """Return short, copyable statements for the manuscript."""
    census = report["census"]
    daily = report["daily"]
    eddies = report["eddy_summary"]
    lines = []
    for cyc in census.index:
        c = census.loc[cyc]
        s = daily.loc[cyc]
        e = eddies if cyc == "All" else eddies.loc[eddies["Cyc"].eq(cyc)]
        lines.extend([
            f"{cyc} eddies: {int(c.unique_eddies):,}.",
            f"Valid {cyc} tilt estimates: {int(c.valid_tilt_days):,}; coverage of eligible eddy-days: {_fmt(c.possible_coverage_pct)}%.",
            f"Median {cyc} daily tilt distance: {_fmt(s.median_km)} km (IQR {_fmt(s.q25_km)}–{_fmt(s.q75_km)} km; maximum {_fmt(s.maximum_km)} km).",
            f"Mean lifetime-average {cyc} tilt distance: {_fmt(e.mean_tilt.mean())} km; maximum observed per-eddy tilt: {_fmt(e.max_tilt.max())} km.",
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
        "preferred_quadrant_NE_AE_or_SE_CE_pct": 100 * preferred_quadrant.mean(),
        "NE_pct": 100 * group["NE"].mean(),
        "SE_pct": 100 * group["SE"].mean(),
        "SW_pct": 100 * group["SW"].mean(),
        "NW_pct": 100 * group["NW"].mean(),
        "median_tilt_km": group["TiltDis"].median(),
        "over_40_km_pct": 100 * group["TiltDis"].gt(40).mean(),
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
    """Summarise direction for all eddies and explicitly for non-shelf eddies.

    ``offshore`` pools U1, U2, D1 and D2. ``offshore_sector`` contrasts
    U1/U2 with D1/D2 without allowing S1 or S2 into either group.
    """
    prepared = prepare_direction_data(
        data, minimum_tilt_km=minimum_tilt_km, offshore_regions=offshore_regions
    )
    offshore = prepared.loc[prepared["is_offshore"]].copy()
    return {
        "data": prepared,
        "offshore_data": offshore,
        "overall": _direction_table(prepared, ["Cyc"]),
        "offshore": _direction_table(offshore, ["Cyc"]),
        "regional": _direction_table(prepared.dropna(subset=["Region"]), ["Cyc", "Region"]),
        "offshore_regional": _direction_table(offshore, ["Cyc", "Region"]),
        "offshore_sector": _direction_table(
            offshore.dropna(subset=["OffshoreSector"]), ["Cyc", "OffshoreSector"]
        ),
        "minimum_tilt_km": minimum_tilt_km,
        "offshore_regions": tuple(offshore_regions),
    }


def direction_census_sentences(report: dict[str, object]) -> list[str]:
    """Return copyable total-population and non-shelf direction statements."""
    lines = [
        f"Direction statistics exclude tilts below {report['minimum_tilt_km']:g} km, for which direction is poorly defined."
    ]
    for label, table in (("All regions", report["overall"]), ("Non-shelf regions (U1, U2, D1 and D2)", report["offshore"])):
        for cyc in ("AE", "CE"):
            if cyc not in table.index:
                continue
            row = table.loc[cyc]
            preferred = "equatorward" if cyc == "AE" else "poleward"
            quadrant = "northeastward" if cyc == "AE" else "southeastward"
            lines.append(
                f"{label} — {cyc}: {int(row.n):,} eddy-days from {int(row.n_eddies):,} eddies; "
                f"{_fmt(row.equatorward_AE_or_poleward_CE_pct)}% tilted {preferred}, "
                f"{_fmt(row.eastward_pct)}% had an eastward component, and "
                f"{_fmt(row.preferred_quadrant_NE_AE_or_SE_CE_pct)}% tilted {quadrant}; "
                f"circular mean {_fmt(row.mean_direction_deg)}° (R = {_fmt(row.directional_concentration_R, 2)})."
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
    return {
        "data": data,
        "overall": _pv_table(data, ["Cyc"]),
        "offshore": _pv_table(offshore, ["Cyc"]),
        "regional": _pv_table(data, ["Cyc", "Region"]),
        "offshore_regional": _pv_table(offshore, ["Cyc", "Region"]),
    }
