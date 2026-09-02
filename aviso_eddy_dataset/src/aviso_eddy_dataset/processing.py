from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from seacofs_eddy_dataset.core.doppio import nearest_ij
from seacofs_eddy_dataset.io import write_partition
from seacofs_eddy_dataset.stages.processing import compute_ar_from_q_columns, fill_missing_eddy_days

from .common import find_data_files, open_aviso, source_file_by_day
from .config import PipelineConfig
from .grid import build_grid


FINAL_COLUMNS = [
    "Eddy", "Day", "Date", "Cyc", "lon", "lat", "ic", "jc", "xc", "yc", "w",
    "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR", "R", "Age", "source_file",
]


def _reference_grid(config: PipelineConfig):
    path = find_data_files(config)[0]
    variables = config.raw.get("variables", {})
    with open_aviso(path, config) as dataset:
        return build_grid(
            dataset[variables.get("longitude", "longitude")].values,
            dataset[variables.get("latitude", "latitude")].values,
        )


def _add_geographic_coordinates(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    grid = _reference_grid(config)
    out = df.copy()
    points = np.column_stack((out["xc"].to_numpy(float), out["yc"].to_numpy(float)))
    out["lon"] = RegularGridInterpolator(
        (grid.x_grid, grid.y_grid),
        np.broadcast_to(grid.longitude[:, None], grid.X_grid.shape),
        bounds_error=False,
        fill_value=np.nan,
    )(points)
    out["lat"] = RegularGridInterpolator(
        (grid.x_grid, grid.y_grid),
        np.broadcast_to(grid.latitude[None, :], grid.Y_grid.shape),
        bounds_error=False,
        fill_value=np.nan,
    )(points)
    ij = [nearest_ij(x, y, grid.X_grid, grid.Y_grid) for x, y in zip(out.xc, out.yc)]
    out["ic"] = [item[0] for item in ij]
    out["jc"] = [item[1] for item in ij]
    return out


def process(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=FINAL_COLUMNS)
    settings = config.raw.get("processing", {})
    out = df.sort_values(["Eddy", "Day"]).copy()
    smooth_cols = ["q11", "q12", "q22"]
    out[smooth_cols] = out.groupby("Eddy", group_keys=False)[smooth_cols].transform(
        lambda values: values.rolling(window=3, center=True, min_periods=1).mean()
    )
    out["AR"] = compute_ar_from_q_columns(out)
    keep = [
        "Eddy", "Day", "Date", "source_file", "Cyc", "xc", "yc", "w", "Omega0",
        "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR", "R",
    ]
    out = out[keep]
    minimum_duration = int(settings.get("min_duration_days", 21))
    out = out.groupby("Eddy").filter(
        lambda group: group.Day.max() - group.Day.min() >= minimum_duration
    )
    if out.empty:
        return pd.DataFrame(columns=FINAL_COLUMNS)

    ar_bad = out.AR > float(settings.get("ar_max", 5.0))
    invalid_shape = ["w", "Omega0", "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR"]
    out.loc[ar_bad, invalid_shape] = np.nan
    poor_outer = (
        (out.Rc > float(settings.get("rc_max_km", 300.0)))
        | (out.Rc > float(settings.get("rc_r_ratio_max", 1.75)) * out.R)
        | (out.Omega.abs() > float(settings.get("omega_abs_max", 5e-5)))
        | (out.psi0.abs() > float(settings.get("psi0_abs_max", 300.0)))
    )
    out.loc[poor_outer, "Omega"] = out.loc[poor_outer, "Omega0"]
    out.loc[poor_outer, ["Rc", "psi0"]] = np.nan
    out = out.drop(columns="Omega0")

    source_by_day = source_file_by_day(config)
    max_gap = settings.get("max_outer_gap_days")
    out = fill_missing_eddy_days(out, np.inf if max_gap is None else float(max_gap))
    origin = pd.Timestamp(settings.get("date_origin", "1990-01-01"))
    out["Date"] = origin + pd.to_timedelta(out.Day, unit="D")
    out["source_file"] = out.Day.map(source_by_day)
    out = _add_geographic_coordinates(out, config)

    required = ["w", "Omega", "q11", "q12", "q22", "Rc", "psi0", "AR", "R"]
    valid = out.groupby("Eddy")[required].apply(lambda group: ~group.isna().all().any())
    out = out[out.Eddy.isin(valid[valid].index)]
    out = out.groupby("Eddy").filter(
        lambda group: group.R.mean() > float(settings.get("mean_radius_min_km", 15.0))
    )
    if out.empty:
        return pd.DataFrame(columns=FINAL_COLUMNS)

    out["Eddy"] = out.Eddy.rank(method="dense").astype(int)
    out["Age"] = out.groupby("Eddy").Eddy.transform("count")
    out["AR"] = compute_ar_from_q_columns(out)
    out["w"] = out.Omega * (out.q11 + out.q22)
    return out[FINAL_COLUMNS].sort_values(["Eddy", "Day"]).reset_index(drop=True)


def run(config: PipelineConfig) -> None:
    source = config.output_root / "tracked" / "eddy_tracks.parquet"
    if not source.exists():
        raise FileNotFoundError(f"Missing tracked dataset: {source}")
    output = config.output_root / "processed" / "eddy_dataset_processed.parquet"
    if config.skip_existing and output.exists():
        print(output)
        return
    write_partition(process(pd.read_parquet(source), config), output)
    print(output)
