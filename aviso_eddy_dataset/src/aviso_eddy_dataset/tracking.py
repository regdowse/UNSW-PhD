from __future__ import annotations

import numpy as np
import pandas as pd

from seacofs_eddy_dataset.core.tracking import clean_surface_eddies, tracking_kdtree_with_omega0
from seacofs_eddy_dataset.io import read_partitions, write_partition

from .common import find_data_files, open_aviso, partition_path
from .config import PipelineConfig
from .grid import build_grid


def _surface_paths(config: PipelineConfig):
    return [partition_path(config, "surface_eddies", path) for path in find_data_files(config)]


def clean(config: PipelineConfig) -> pd.DataFrame:
    paths = [path for path in _surface_paths(config) if path.exists()]
    surface = read_partitions(paths)
    if surface.empty:
        return surface

    first_file = find_data_files(config)[0]
    variables = config.raw.get("variables", {})
    with open_aviso(first_file, config) as dataset:
        grid = build_grid(
            dataset[variables.get("longitude", "longitude")].values,
            dataset[variables.get("latitude", "latitude")].values,
        )
    settings = config.raw.get("surface_fit", {})
    return clean_surface_eddies(
        surface,
        grid.X_grid,
        grid.Y_grid,
        threshold=float(settings.get("center_error_mad_threshold", 4.0)),
        omega0_thresh=float(settings.get("omega0_abs_threshold", 5e-5)),
    )


def track(config: PipelineConfig) -> pd.DataFrame:
    surface = clean(config)
    if surface.empty:
        return surface
    surface["Cyc"] = surface["nCyc"]
    surface = surface.sort_values(["Day", "yc", "xc"]).reset_index(drop=True)
    surface["eddy_idx"] = surface.groupby("Day").cumcount()
    first_count = int(surface["Day"].eq(surface["Day"].min()).sum())
    settings = config.raw.get("tracking", {})
    return tracking_kdtree_with_omega0(
        surface,
        start_ids=np.arange(1, first_count + 1),
        next_num=first_count + 1,
        length_scale=float(settings.get("length_scale_km", 50.0)),
        omega0_scale=float(settings.get("omega0_scale", 1e-5)),
        radius_threshold=float(settings.get("radius_threshold", 1.0)),
        lookback=int(settings.get("lookback_days", 4)),
    )


def run(config: PipelineConfig) -> None:
    output = config.output_root / "tracked" / "eddy_tracks.parquet"
    if config.skip_existing and output.exists():
        print(output)
        return
    write_partition(track(config), output)
    print(output)
