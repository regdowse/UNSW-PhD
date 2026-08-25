from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.grid import (
    fnumber_from_outer_avg,
    read_reference_grid,
    regular_xy_grid,
)
from seacofs_eddy_dataset.core.nencioli import nencioli
from seacofs_eddy_dataset.core.velocity import interpolate_uv_to_grid
from seacofs_eddy_dataset.io import partition_path, write_partition


DETECTION_COLUMNS = ["Day", "nxc", "nyc", "Cyc", "nic", "njc", "fnumber"]


def find_model_files(config: PipelineConfig) -> list[Path]:
    pattern = config.raw.get("files", {}).get("model_pattern", "outer_avg_*.nc")
    return sorted(config.model_root.glob(pattern))


def reference_grid_path(config: PipelineConfig) -> Path:
    reference = config.raw.get("files", {}).get("grid_reference", "outer_avg_01461.nc")
    return config.model_root / reference


def detection_output_path(path: Path, config: PipelineConfig) -> Path:
    fnumber = fnumber_from_outer_avg(path)
    return partition_path(config.output_root, "detections", f"fnumber={fnumber:05}")


def _empty_detection_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=DETECTION_COLUMNS)


def detect_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    """Detect Nencioli eddy candidates for one model file.

    This stage should load one `outer_avg_*.nc`, use `core.velocity` to rotate
    and interpolate surface velocities, then call `core.nencioli.nencioli`.

    Expected output columns include Day, nxc, nyc, Cyc, nic, njc, and fnumber.
    """
    fnumber = fnumber_from_outer_avg(path)
    grid = read_reference_grid(reference_grid_path(config))

    nencioli_config = config.raw.get("nencioli", {})
    a = int(nencioli_config.get("a", 4))
    b = int(nencioli_config.get("b", 3))
    resolution_km = float(nencioli_config.get("interpolation_resolution_km", 1.0))
    X_new, Y_new = regular_xy_grid(grid, resolution_km)

    source_points = np.column_stack((grid.X_grid.ravel(), grid.Y_grid.ravel()))
    source_tree = cKDTree(source_points)
    frames: list[pd.DataFrame] = []

    with nc.Dataset(path) as dataset:
        ocean_time = np.asarray(dataset.variables["ocean_time"][:].data, dtype=float) / 86400

        for t, day_value in enumerate(ocean_time):
            day = int(round(day_value))
            u = dataset["u_eastward"][t, -1, :, :].T
            v = dataset["v_northward"][t, -1, :, :].T
            u_interp, v_interp = interpolate_uv_to_grid(
                u, v, grid.x_grid, grid.y_grid, X_new, Y_new, grid.angle
            )

            detected = nencioli(u_interp.T, v_interp.T, X_new.T, Y_new.T, a, b)[2]
            if len(detected) == 0:
                continue

            detected = detected[detected[:, 1].argsort()[::-1]]
            _, idx = source_tree.query(detected[:, :2])
            nic, njc = np.unravel_index(idx, grid.X_grid.shape)

            frames.append(
                pd.DataFrame(
                    {
                        "Day": day,
                        "nxc": detected[:, 0],
                        "nyc": detected[:, 1],
                        "Cyc": np.where(detected[:, 2] == -1, "CE", "AE"),
                        "nic": nic,
                        "njc": njc,
                        "fnumber": fnumber,
                    }
                )
            )

    if not frames:
        return _empty_detection_frame()
    return pd.concat(frames, ignore_index=True)[DETECTION_COLUMNS]


def detect_and_write_file(path: Path, config: PipelineConfig) -> Path:
    out_path = detection_output_path(path, config)
    if config.skip_existing and out_path.exists():
        return out_path
    detections = detect_file(path, config)
    return write_partition(detections, out_path)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    files = find_model_files(config)
    if not files:
        raise FileNotFoundError(f"No model files found in {config.model_root}")

    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(detect_and_write_file)(path, config) for path in files
    )
    for path in written:
        print(path)
