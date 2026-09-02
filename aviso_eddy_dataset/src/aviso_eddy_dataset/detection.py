from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from seacofs_eddy_dataset.core.nencioli import nencioli
from seacofs_eddy_dataset.io import write_partition

from .common import day_number, find_data_files, open_aviso, partition_path
from .config import PipelineConfig
from .grid import build_grid, native_velocity, regular_xy_grid


DETECTION_COLUMNS = ["Day", "Date", "source_file", "nxc", "nyc", "Cyc", "nic", "njc"]


def interpolate_for_nencioli(u, v, grid, X_new, Y_new):
    """Interpolate native unrotated AVISO velocities to the detection grid."""

    points = np.column_stack((X_new.ravel(), Y_new.ravel()))
    kwargs = {"method": "linear", "bounds_error": False, "fill_value": np.nan}
    u_interp = RegularGridInterpolator((grid.x_grid, grid.y_grid), u, **kwargs)
    v_interp = RegularGridInterpolator((grid.x_grid, grid.y_grid), v, **kwargs)
    return u_interp(points).reshape(X_new.shape), v_interp(points).reshape(X_new.shape)


def detect_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    variables = config.raw.get("variables", {})
    lon_name = variables.get("longitude", "longitude")
    lat_name = variables.get("latitude", "latitude")
    time_name = variables.get("time", "time")
    u_name = variables.get("u", "ugos")
    v_name = variables.get("v", "vgos")
    settings = config.raw.get("nencioli", {})
    a = int(settings.get("a", 4))
    b = int(settings.get("b", 3))
    resolution = float(settings.get("interpolation_resolution_km", 1.0))
    frames: list[pd.DataFrame] = []

    with open_aviso(path, config) as dataset:
        grid = build_grid(dataset[lon_name].values, dataset[lat_name].values)
        X_new, Y_new = regular_xy_grid(grid, resolution)
        native_tree = cKDTree(np.column_stack((grid.X_grid.ravel(), grid.Y_grid.ravel())))

        for time_index, time_value in enumerate(dataset[time_name].values):
            u = native_velocity(dataset[u_name].isel({time_name: time_index}), grid)
            v = native_velocity(dataset[v_name].isel({time_name: time_index}), grid)
            u_interp, v_interp = interpolate_for_nencioli(u, v, grid, X_new, Y_new)
            detected = nencioli(u_interp.T, v_interp.T, X_new.T, Y_new.T, a, b)[2]
            if len(detected) == 0:
                continue

            detected = detected[detected[:, 1].argsort()[::-1]]
            _, native_flat = native_tree.query(detected[:, :2])
            nic, njc = np.unravel_index(native_flat, grid.X_grid.shape)
            date = pd.Timestamp(time_value).normalize()
            frames.append(
                pd.DataFrame(
                    {
                        "Day": day_number(time_value, config),
                        "Date": date,
                        "source_file": str(path),
                        "nxc": detected[:, 0],
                        "nyc": detected[:, 1],
                        "Cyc": np.where(detected[:, 2] == -1, "CE", "AE"),
                        "nic": nic,
                        "njc": njc,
                    }
                )
            )

    if not frames:
        return pd.DataFrame(columns=DETECTION_COLUMNS)
    return pd.concat(frames, ignore_index=True)[DETECTION_COLUMNS]


def run_file(path: Path, config: PipelineConfig) -> Path:
    output = partition_path(config, "detections", path)
    if config.skip_existing and output.exists():
        return output
    return write_partition(detect_file(path, config), output)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(run_file)(path, config) for path in find_data_files(config)
    )
    for path in written:
        print(path)
