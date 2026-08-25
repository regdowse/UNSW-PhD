from __future__ import annotations

import numpy as np

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.grid import read_reference_grid
from seacofs_eddy_dataset.core.tracking import clean_surface_eddies as clean_surface_eddies_core
from seacofs_eddy_dataset.core.tracking import tracking_kdtree_with_omega0
from seacofs_eddy_dataset.io import read_partitions, write_partition
from seacofs_eddy_dataset.stages.detection import reference_grid_path


def _surface_paths(config: PipelineConfig):
    return sorted((config.output_root / "surface_eddies").glob("fnumber=*.parquet"))


def clean_surface_eddies(config: PipelineConfig):
    df = read_partitions(_surface_paths(config))
    if df.empty:
        return df
    grid = read_reference_grid(reference_grid_path(config))
    settings = config.raw.get("surface_fit", {})
    return clean_surface_eddies_core(
        df,
        grid.X_grid,
        grid.Y_grid,
        threshold=float(settings.get("center_error_mad_threshold", 4.0)),
        omega0_thresh=float(settings.get("omega0_abs_threshold", 5e-5)),
    )


def track_eddies(config: PipelineConfig):
    df = clean_surface_eddies(config)
    if df.empty:
        return df

    scale = float(config.raw.get("surface_fit", {}).get("omega_units_scale", 1e-3))
    for col in ["w", "Omega0", "Omega"]:
        if col in df.columns:
            df[col] = df[col] * scale
    df["Cyc"] = df["nCyc"]
    df = df.sort_values(["Day", "yc", "xc"]).reset_index(drop=True)
    df["eddy_idx"] = df.groupby("Day").cumcount()
    first_count = int(df["Day"].eq(df["Day"].min()).sum())
    start_ids = np.arange(1, first_count + 1)

    settings = config.raw.get("tracking", {})
    return tracking_kdtree_with_omega0(
        df,
        start_ids=start_ids,
        next_num=first_count + 1,
        length_scale=float(settings.get("length_scale_km", 50.0)),
        omega0_scale=float(settings.get("omega0_scale", 1e-5)),
        radius_threshold=float(settings.get("radius_threshold", 1.0)),
        lookback=int(settings.get("lookback_days", 4)),
    )


def run(config: PipelineConfig) -> None:
    out_path = config.output_root / "tracked" / "eddy_tracks.parquet"
    if config.skip_existing and out_path.exists():
        print(out_path)
        return
    tracked = track_eddies(config)
    write_partition(tracked, out_path)
    print(out_path)
