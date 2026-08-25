from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.doppio import nearest_ij
from seacofs_eddy_dataset.core.grid import read_reference_grid
from seacofs_eddy_dataset.io import write_partition
from seacofs_eddy_dataset.stages.detection import reference_grid_path


FINAL_COLUMNS = [
    "Eddy",
    "Day",
    "Cyc",
    "lon",
    "lat",
    "ic",
    "jc",
    "xc",
    "yc",
    "w",
    "Omega",
    "q11",
    "q12",
    "q22",
    "Rc",
    "psi0",
    "AR",
    "R",
    "Age",
    "Date",
    "fname",
]


def compute_ar_from_q_columns(df: pd.DataFrame) -> np.ndarray:
    q11 = df["q11"].to_numpy(dtype=float)
    q12 = df["q12"].to_numpy(dtype=float)
    q22 = df["q22"].to_numpy(dtype=float)
    trace = q11 + q22
    radius = np.sqrt((q11 - q22) ** 2 + 4 * q12**2)
    lam_min = 0.5 * (trace - radius)
    lam_max = 0.5 * (trace + radius)
    ar = np.sqrt(lam_max / np.maximum(lam_min, 1e-12))
    ar[lam_min <= 0] = np.nan
    return ar


def fill_missing_eddy_days(df: pd.DataFrame, max_outer_gap: float = np.inf) -> pd.DataFrame:
    out = []
    interp_cols = ["xc", "yc", "w", "Omega", "q11", "q12", "q22", "R"]
    outer_cols = ["Rc", "psi0"]

    for eddy, group in df.groupby("Eddy"):
        group = group.sort_values("Day").copy()
        full_days = pd.DataFrame(
            {
                "Eddy": eddy,
                "Day": np.arange(int(group["Day"].min()), int(group["Day"].max()) + 1),
            }
        )
        merged = full_days.merge(group, on=["Eddy", "Day"], how="left")

        existing_interp_cols = [col for col in interp_cols if col in merged.columns]
        merged[existing_interp_cols] = merged[existing_interp_cols].interpolate(
            method="linear", limit_direction="both"
        )

        if "Cyc" in merged.columns:
            merged["Cyc"] = merged["Cyc"].ffill().bfill()

        for col in outer_cols:
            if col not in merged.columns:
                continue
            interpolated = merged[col].interpolate(method="linear", limit_direction="both")
            is_nan = merged[col].isna().to_numpy()
            start = None
            for i, bad in enumerate(is_nan):
                if bad and start is None:
                    start = i
                elif not bad and start is not None:
                    if i - start <= max_outer_gap:
                        merged.loc[merged.index[start:i], col] = interpolated.iloc[start:i]
                    start = None
            if start is not None and len(merged) - start <= max_outer_gap:
                merged.loc[merged.index[start:], col] = interpolated.iloc[start:]

        out.append(merged)

    if not out:
        return df.iloc[0:0].copy()
    return pd.concat(out, ignore_index=True).sort_values(["Eddy", "Day"]).reset_index(drop=True)


def _add_lon_lat_ic_jc(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    grid = read_reference_grid(reference_grid_path(config))
    out = df.copy()
    points = np.column_stack((out["xc"].to_numpy(float), out["yc"].to_numpy(float)))
    lon_interp = RegularGridInterpolator(
        (grid.x_grid, grid.y_grid), grid.lon_rho, bounds_error=False, fill_value=np.nan
    )
    lat_interp = RegularGridInterpolator(
        (grid.x_grid, grid.y_grid), grid.lat_rho, bounds_error=False, fill_value=np.nan
    )
    out["lon"] = lon_interp(points)
    out["lat"] = lat_interp(points)
    ij = [nearest_ij(xc, yc, grid.X_grid, grid.Y_grid) for xc, yc in zip(out["xc"], out["yc"], strict=False)]
    out["ic"] = [i for i, _ in ij]
    out["jc"] = [j for _, j in ij]
    return out


def _fname_from_day(day: int, config: PipelineConfig) -> str:
    settings = config.raw.get("processing", {})
    file_start = int(settings.get("file_start_fnumber", 1461))
    file_stride = int(settings.get("file_stride_days", 30))
    day_offset = int(settings.get("file_day_offset", file_start + 1))
    fnumber = file_start + ((int(day) - day_offset) // file_stride) * file_stride
    return str(config.model_root / f"outer_avg_{fnumber:05}.nc")


def process_tracked_dataset(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    settings = config.raw.get("processing", {})
    min_duration_days = int(settings.get("min_duration_days", 21))
    ar_max = float(settings.get("ar_max", 5.0))
    rc_max = float(settings.get("rc_max_km", 300.0))
    rc_r_ratio_max = float(settings.get("rc_r_ratio_max", 1.75))
    omega_abs_max = float(settings.get("omega_abs_max", 5e-5))
    psi0_abs_max = float(settings.get("psi0_abs_max", 300.0))
    mean_r_min = float(settings.get("mean_radius_min_km", 15.0))
    max_outer_gap = settings.get("max_outer_gap_days", None)
    max_outer_gap = np.inf if max_outer_gap is None else float(max_outer_gap)

    out = df.copy().sort_values(["Eddy", "Day"]).reset_index(drop=True)
    smooth_cols = ["q11", "q12", "q22"]
    out[smooth_cols] = (
        out.groupby("Eddy", group_keys=False)[smooth_cols]
        .transform(lambda values: values.rolling(window=3, center=True, min_periods=1).mean())
    )
    out["AR"] = compute_ar_from_q_columns(out)

    keep = ["Eddy", "Day", "Cyc", "xc", "yc", "w", "Omega0", "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR", "R"]
    out = out[keep].sort_values(["Eddy", "Day"]).copy()
    out = out.groupby("Eddy").filter(lambda group: group["Day"].max() - group["Day"].min() >= min_duration_days)

    ar_mask = out["AR"] > ar_max
    out.loc[ar_mask, ["w", "Omega0", "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR"]] = np.nan

    poor_outer_fits = (
        (out["Rc"] > rc_max)
        | (out["Rc"] > rc_r_ratio_max * out["R"])
        | (out["Omega"].abs() > omega_abs_max)
        | (out["psi0"].abs() > psi0_abs_max)
    )
    out.loc[poor_outer_fits, "Omega"] = out.loc[poor_outer_fits, "Omega0"]
    out.loc[poor_outer_fits, ["Rc", "psi0"]] = np.nan
    out = out.drop(columns=["Omega0"])

    out = fill_missing_eddy_days(out, max_outer_gap=max_outer_gap)
    out["Eddy"] = out["Eddy"].rank(method="dense").astype(int)
    out = out.reset_index(drop=True)
    out = _add_lon_lat_ic_jc(out, config)

    check_cols = ["w", "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR", "R"]
    valid_eddies = out.groupby("Eddy")[check_cols].apply(lambda group: ~group.isna().all().any())
    out = out[out["Eddy"].isin(valid_eddies[valid_eddies].index)].copy().reset_index(drop=True)
    out = out.groupby("Eddy").filter(lambda group: group["R"].mean() > mean_r_min).reset_index(drop=True)

    out["Eddy"] = out["Eddy"].rank(method="dense").astype(int)
    out["Age"] = out.groupby("Eddy")["Eddy"].transform("count")
    out["Date"] = pd.Timestamp(settings.get("date_origin", "1990-01-01")) + pd.to_timedelta(out["Day"], unit="D")
    out["fname"] = [_fname_from_day(day, config) for day in out["Day"]]
    out["AR"] = compute_ar_from_q_columns(out)
    out["w"] = out["Omega"] * (out["q11"] + out["q22"])
    return out[FINAL_COLUMNS].sort_values(["Eddy", "Day"]).reset_index(drop=True)


def run(config: PipelineConfig) -> None:
    out_path = config.output_root / "processed" / "eddy_dataset_processed.parquet"
    if config.skip_existing and out_path.exists():
        print(out_path)
        return

    source = config.output_root / "tracked" / "eddy_tracks.parquet"
    if not source.exists():
        raise FileNotFoundError(f"Missing tracked dataset: {source}")

    df = pd.read_parquet(source)
    processed = process_tracked_dataset(df, config) if not df.empty else pd.DataFrame(columns=FINAL_COLUMNS)
    write_partition(processed, out_path)
    print(out_path)
