"""Geometry helpers for testing rotation-dependent shear/lift eddy tilt.

The measured ``TiltDir`` and ``TiltDis`` fields are authoritative.  This
module only expresses that measured vector in coordinates defined by an
independently estimated environmental vertical-shear vector.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mechanism_tools import add_tilt_components, bearing_from_east_north


def add_200_500_layer_mean(df: pd.DataFrame, families=("clim", "full")) -> pd.DataFrame:
    """Recover the 200--500 m mean from cached 0--200 and 0--500 m means."""

    out = df.copy()
    for family in families:
        for component in ("east", "north"):
            upper = f"{family}_200_{component}_ms"
            full = f"{family}_500_{component}_ms"
            missing = {upper, full} - set(out.columns)
            if missing:
                raise KeyError(f"Cannot calculate {family} 200--500 m mean; missing {sorted(missing)}")
            out[f"{family}_200_500_{component}_ms"] = (
                500.0 * out[full] - 200.0 * out[upper]
            ) / 300.0
    return out


def add_shear_vectors(df: pd.DataFrame, families=("clim", "full")) -> pd.DataFrame:
    """Add the two predefined upper-minus-deep environmental shear vectors."""

    out = add_200_500_layer_mean(df, families=families)
    definitions = {
        "surface_deep": ("surface", "200_500"),
        "upper_deep": ("200", "200_500"),
    }
    for family in families:
        for name, (upper, deep) in definitions.items():
            for component in ("east", "north"):
                out[f"{family}_{name}_{component}_ms"] = (
                    out[f"{family}_{upper}_{component}_ms"]
                    - out[f"{family}_{deep}_{component}_ms"]
                )
    return out


def add_shear_relative_tilt(
    df: pd.DataFrame,
    shear_east: str,
    shear_north: str,
    *,
    prefix: str,
    minimum_shear_ms: float = 0.0,
) -> pd.DataFrame:
    """Project measured tilt parallel and leftward of upper-minus-deep shear.

    ``left`` is positive along k x shear = (-north, east).  No AE/CE sign is
    embedded in the geometry; polarity reversal is tested statistically.
    """

    out = add_tilt_components(df)
    missing = {shear_east, shear_north} - set(out.columns)
    if missing:
        raise KeyError(f"Missing shear columns: {sorted(missing)}")
    east = out[shear_east].astype(float)
    north = out[shear_north].astype(float)
    magnitude = np.hypot(east, north)
    valid = np.isfinite(magnitude) & (magnitude > float(minimum_shear_ms))
    unit_east = np.divide(east, magnitude, out=np.full(len(out), np.nan), where=valid)
    unit_north = np.divide(north, magnitude, out=np.full(len(out), np.nan), where=valid)

    out[f"{prefix}_shear_mag_ms"] = magnitude
    out[f"{prefix}_shear_dir"] = bearing_from_east_north(east, north)
    out[f"{prefix}_unit_east"] = unit_east
    out[f"{prefix}_unit_north"] = unit_north
    out[f"{prefix}_tilt_parallel_km"] = (
        out["tilt_east_km"] * unit_east + out["tilt_north_km"] * unit_north
    )
    out[f"{prefix}_tilt_left_km"] = (
        -out["tilt_east_km"] * unit_north + out["tilt_north_km"] * unit_east
    )
    out[f"{prefix}_tilt_parallel_fraction"] = np.divide(
        out[f"{prefix}_tilt_parallel_km"], out["TiltDis"],
        out=np.full(len(out), np.nan), where=out["TiltDis"].to_numpy(float) > 0,
    )
    out[f"{prefix}_tilt_left_fraction"] = np.divide(
        out[f"{prefix}_tilt_left_km"], out["TiltDis"],
        out=np.full(len(out), np.nan), where=out["TiltDis"].to_numpy(float) > 0,
    )
    out[f"{prefix}_offset_deg"] = np.degrees(np.arctan2(
        out[f"{prefix}_tilt_left_km"], out[f"{prefix}_tilt_parallel_km"]
    ))
    return out


def add_lagged_tilt_change(
    df: pd.DataFrame,
    *,
    lag_days: int,
    prefix: str,
) -> pd.DataFrame:
    """Project future-minus-current measured tilt onto the current shear axes.

    Only exact day separations are matched. This avoids treating gaps as a
    continuous daily tendency and does not recompute either tilt endpoint.
    """

    if lag_days <= 0:
        raise ValueError("lag_days must be positive")
    required = {
        "Eddy", "Day", "tilt_east_km", "tilt_north_km",
        f"{prefix}_unit_east", f"{prefix}_unit_north",
    }
    if missing := required - set(df.columns):
        raise KeyError(f"Missing lagged-change columns: {sorted(missing)}")
    if df.duplicated(["Eddy", "Day"]).any():
        raise ValueError("Eddy-Day rows must be unique")

    future = df[["Eddy", "Day", "tilt_east_km", "tilt_north_km"]].copy()
    future["Day"] = future["Day"] - lag_days
    future = future.rename(columns={
        "tilt_east_km": "future_tilt_east_km",
        "tilt_north_km": "future_tilt_north_km",
    })
    out = df.merge(future, on=["Eddy", "Day"], how="left", validate="one_to_one")
    de = out["future_tilt_east_km"] - out["tilt_east_km"]
    dn = out["future_tilt_north_km"] - out["tilt_north_km"]
    ue = out[f"{prefix}_unit_east"]
    un = out[f"{prefix}_unit_north"]
    stem = f"{prefix}_{lag_days}d_dtilt"
    out[f"{stem}_east_km"] = de
    out[f"{stem}_north_km"] = dn
    out[f"{stem}_parallel_km"] = de * ue + dn * un
    out[f"{stem}_left_km"] = -de * un + dn * ue
    out[f"{stem}_parallel_km_day"] = out[f"{stem}_parallel_km"] / lag_days
    out[f"{stem}_left_km_day"] = out[f"{stem}_left_km"] / lag_days
    return out


def eddy_equal_binned_summary(
    df: pd.DataFrame,
    predictor: str,
    response: str,
    *,
    bins: int = 6,
) -> pd.DataFrame:
    """Equal-count predictor bins with one median response per eddy per bin."""

    q = df[["Eddy", predictor, response]].dropna().copy()
    if q.empty:
        return pd.DataFrame(columns=["bin", "x", "y", "eddies", "observations"])
    q["bin"] = pd.qcut(q[predictor].rank(method="first"), bins, duplicates="drop")
    per_eddy = q.groupby(["bin", "Eddy"], observed=True).agg(
        x=(predictor, "median"), y=(response, "median"), observations=(response, "size")
    ).reset_index()
    return per_eddy.groupby("bin", observed=True).agg(
        x=("x", "median"), y=("y", "median"), eddies=("Eddy", "nunique"),
        observations=("observations", "sum"),
    ).reset_index()
