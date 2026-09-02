from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from seacofs_eddy_dataset.core.doppio import find_directional_radii, transect_indexer
from seacofs_eddy_dataset.core.esp import load_doppio_functions
from seacofs_eddy_dataset.io import write_partition

from .common import day_number, find_data_files, open_aviso, partition_path
from .config import PipelineConfig
from .grid import build_grid, native_velocity


SURFACE_COLUMNS = [
    "Day", "Date", "source_file", "nxc", "nyc", "nCyc", "nic", "njc",
    "xc", "yc", "w", "q11", "q12", "q22", "Omega0", "Omega", "Rc", "psi0", "R",
]


def _bad_row(row) -> dict:
    return {
        "Day": row.Day, "Date": row.Date, "source_file": row.source_file,
        "nxc": row.nxc, "nyc": row.nyc, "nCyc": row.Cyc, "nic": row.nic, "njc": row.njc,
        "xc": np.nan, "yc": np.nan, "w": np.nan, "q11": np.nan, "q12": np.nan,
        "q22": np.nan, "Omega0": np.nan, "Omega": np.nan, "Rc": np.nan,
        "psi0": np.nan, "R": np.nan,
    }


def fit_detection(row, u, v, grid, doppio, out_core_param_fit, settings: dict) -> dict:
    """Fit one Nencioli candidate with DOPPIO on the native AVISO grid."""

    radius_km = float(settings.get("transect_radius_km", 30.0))
    rho_max = float(settings.get("rho_max_km", 200.0))
    rho_min = float(settings.get("rho_min_km", 30.0))
    out_core_fac = float(settings.get("out_core_fac", 1.75))
    local_limit_factor = float(settings.get("local_limit_factor", 3.0))
    omega_scale = float(settings.get("omega_units_scale", 1e-3))

    try:
        ic, jc = int(row.nic), int(row.njc)
        x1, y1, x2, y2, ii, jj = transect_indexer(ic, jc, grid.X_grid, grid.Y_grid, r=radius_km)
        if len(ii) < 3 or len(jj) < 3:
            return _bad_row(row)
        transects = (u[ii, jc], v[ii, jc], u[ic, jj], v[ic, jj])
        if any(np.count_nonzero(np.isfinite(values)) < 3 for values in transects):
            return _bad_row(row)
        xc, yc, w, q, omega0 = doppio(x1, y1, *transects[:2], x2, y2, *transects[2:])
        q = np.asarray(q, dtype=float)
        if q.shape != (2, 2) or not np.all(np.isfinite(q)):
            return _bad_row(row)

        radii = find_directional_radii(u, v, grid.X_grid, grid.Y_grid, xc, yc, q)
        radius_values = np.asarray(list(radii.values()), dtype=float)
        radius_values = radius_values[np.isfinite(radius_values)]
        if radius_values.size == 0:
            return _bad_row(row)
        radius = float(radius_values.mean())

        rho_limit = max(min(radius * out_core_fac, rho_max), rho_min)
        local_limit = rho_limit * local_limit_factor
        local = (
            (grid.X_grid >= xc - local_limit) & (grid.X_grid <= xc + local_limit)
            & (grid.Y_grid >= yc - local_limit) & (grid.Y_grid <= yc + local_limit)
        )
        dx = grid.X_grid[local] - xc
        dy = grid.Y_grid[local] - yc
        rho2 = q[0, 0] * dx**2 + 2 * q[0, 1] * dx * dy + q[1, 1] * dy**2
        rho = np.sqrt(np.where(rho2 >= 0, rho2, np.nan))
        mask = np.isfinite(rho) & (rho <= rho_limit)
        xloc, yloc = grid.X_grid[local][mask], grid.Y_grid[local][mask]
        uloc, vloc = u[local][mask], v[local][mask]
        finite = np.isfinite(xloc) & np.isfinite(yloc) & np.isfinite(uloc) & np.isfinite(vloc)
        if int(finite.sum()) < 10:
            return _bad_row(row)
        rc, psi0, omega = out_core_param_fit(
            xloc[finite], yloc[finite], uloc[finite], vloc[finite], xc, yc, q, w
        )
    except Exception:
        return _bad_row(row)

    return {
        "Day": int(row.Day), "Date": row.Date, "source_file": row.source_file,
        "nxc": float(row.nxc), "nyc": float(row.nyc), "nCyc": row.Cyc,
        "nic": int(row.nic), "njc": int(row.njc), "xc": float(xc), "yc": float(yc),
        "w": float(w) * omega_scale, "q11": float(q[0, 0]), "q12": float(q[0, 1]),
        "q22": float(q[1, 1]), "Omega0": float(omega0) * omega_scale,
        "Omega": float(omega) * omega_scale, "Rc": float(rc), "psi0": float(psi0),
        "R": radius,
    }


def fit_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    detection_path = partition_path(config, "detections", path)
    if not detection_path.exists():
        raise FileNotFoundError(f"Missing detection partition: {detection_path}")
    detections = pd.read_parquet(detection_path)
    if detections.empty:
        return pd.DataFrame(columns=SURFACE_COLUMNS)

    variables = config.raw.get("variables", {})
    lon_name = variables.get("longitude", "longitude")
    lat_name = variables.get("latitude", "latitude")
    time_name = variables.get("time", "time")
    u_name = variables.get("u", "ugos")
    v_name = variables.get("v", "vgos")
    doppio, out_core_param_fit = load_doppio_functions(config)
    settings = config.raw.get("surface_fit", {})
    rows: list[dict] = []

    with open_aviso(path, config) as dataset:
        grid = build_grid(dataset[lon_name].values, dataset[lat_name].values)
        for time_index, time_value in enumerate(dataset[time_name].values):
            day = day_number(time_value, config)
            candidates = detections.loc[detections["Day"].eq(day)]
            if candidates.empty:
                continue
            u = native_velocity(dataset[u_name].isel({time_name: time_index}), grid)
            v = native_velocity(dataset[v_name].isel({time_name: time_index}), grid)
            rows.extend(
                fit_detection(row, u, v, grid, doppio, out_core_param_fit, settings)
                for row in candidates.itertuples(index=False)
            )
    return pd.DataFrame(rows, columns=SURFACE_COLUMNS)


def run_file(path: Path, config: PipelineConfig) -> Path:
    output = partition_path(config, "surface_eddies", path)
    if config.skip_existing and output.exists():
        return output
    return write_partition(fit_file(path, config), output)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(run_file)(path, config) for path in find_data_files(config)
    )
    for path in written:
        print(path)
