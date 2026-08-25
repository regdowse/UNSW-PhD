from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

from seacofs_eddy_dataset.config import PipelineConfig
from seacofs_eddy_dataset.core.doppio import bad_doppio_row, find_directional_radii, transect_indexer
from seacofs_eddy_dataset.core.esp import load_doppio_functions
from seacofs_eddy_dataset.core.grid import fnumber_from_outer_avg, read_reference_grid
from seacofs_eddy_dataset.core.velocity import rotate_uv
from seacofs_eddy_dataset.core.vertical import interp_3d_to_reference_depths
from seacofs_eddy_dataset.io import partition_path, write_partition
from seacofs_eddy_dataset.stages.detection import find_model_files, reference_grid_path


SURFACE_COLUMNS = [
    "Day",
    "fnumber",
    "nxc",
    "nyc",
    "nCyc",
    "nic",
    "njc",
    "xc",
    "yc",
    "w",
    "q11",
    "q12",
    "q22",
    "Omega0",
    "Omega",
    "Rc",
    "psi0",
    "R",
]


def detection_path_for_file(path: Path, config: PipelineConfig) -> Path:
    return partition_path(config.output_root, "detections", f"fnumber={fnumber_from_outer_avg(path):05}")


def surface_output_path(path: Path, config: PipelineConfig) -> Path:
    return partition_path(config.output_root, "surface_eddies", f"fnumber={fnumber_from_outer_avg(path):05}")


def _row_with_q(row: dict, q) -> dict:
    out = dict(row)
    if np.asarray(q).shape == (2, 2):
        out["q11"] = float(q[0, 0])
        out["q12"] = float(q[0, 1])
        out["q22"] = float(q[1, 1])
    else:
        out["q11"] = np.nan
        out["q12"] = np.nan
        out["q22"] = np.nan
    out.pop("Q", None)
    return out


def _bad_row(row, fnumber: int) -> dict:
    base = bad_doppio_row(row, fnumber)
    base["nic"] = getattr(row, "nic", np.nan)
    base["njc"] = getattr(row, "njc", np.nan)
    return _row_with_q(base, np.nan)


def _target_depths_for_surface_check(z_r, settings: dict) -> np.ndarray: # check this
    configured_depths = settings.get("vertical_check_target_depths")
    if configured_depths is None:
        target_depths = np.abs(z_r[150,150,1:])
    else:
        target_depths = np.asarray(configured_depths, dtype=float)

    target_depths = np.sort(target_depths[np.isfinite(target_depths) & (target_depths > 0)])
    depth_threshold = float(settings.get("vertical_check_depth_m", 50.0))
    below = target_depths[target_depths <= depth_threshold]
    above = target_depths[target_depths >= depth_threshold]
    if above.size:
        return np.unique(np.concatenate([below, above[:1]]))
    return below


def _load_z_r(config: PipelineConfig):
    z_r_path = config.raw.get("paths", {}).get("z_r")
    if not z_r_path:
        raise KeyError("paths.z_r is required when surface_fit.require_vertical_profile is true")
    return np.load(Path(z_r_path).expanduser())


def _has_vertical_profile_to_depth(
    *,
    xc_surface: float,
    yc_surface: float,
    w_surface: float,
    row,
    u_depth,
    v_depth,
    target_depths,
    grid,
    doppio,
    radius_km: float,
    max_jump_km: float,
    depth_threshold_m: float,
) -> bool:
    xc_prev = xc_surface
    yc_prev = yc_surface
    ic = int(row.nic)
    jc = int(row.njc)

    for depth_index, target_depth in enumerate(target_depths):
        if (
            (xc_prev < radius_km)
            or (xc_prev > grid.X_grid.max() - radius_km)
            or (yc_prev < radius_km)
            or (yc_prev > grid.Y_grid.max() - radius_km)
        ):
            return False

        u2d, v2d = rotate_uv(u_depth[:, :, depth_index], v_depth[:, :, depth_index], grid.angle)
        x1, y1, x2, y2, ii, jj = transect_indexer(ic, jc, grid.X_grid, grid.Y_grid, r=radius_km)
        u1 = u2d[ii, jc]
        v1 = v2d[ii, jc]
        u2 = u2d[ic, jj]
        v2 = v2d[ic, jj]
        if any(np.all(np.isnan(values)) for values in [u1, v1, u2, v2]):
            return False

        try:
            xc, yc, w, _, _ = doppio(x1, y1, u1, v1, x2, y2, u2, v2)
        except Exception:
            return False

        if np.sign(w) != np.sign(w_surface) or np.hypot(xc - xc_prev, yc - yc_prev) > max_jump_km:
            return False

        if target_depth >= depth_threshold_m:
            return True

        xc_prev = xc
        yc_prev = yc

    return False


def _fit_detection_row(
    row,
    fnumber: int,
    u_rot,
    v_rot,
    grid,
    doppio,
    out_core_param_fit,
    radius_km: float,
    rho_max: float,
    rho_min: float,
    out_core_fac: float,
    omega_scale: float,
    local_limit_factor: float,
    vertical_check: dict | None = None,
) -> dict:
    try:
        x1, y1, x2, y2, ii, jj = transect_indexer(int(row.nic), int(row.njc), grid.X_grid, grid.Y_grid, r=radius_km)
        if len(ii) < 3 or len(jj) < 3:
            return _bad_row(row, fnumber)
        xc, yc, w, q, omega0 = doppio(
            x1,
            y1,
            u_rot[ii, int(row.njc)],
            v_rot[ii, int(row.njc)],
            x2,
            y2,
            u_rot[int(row.nic), jj],
            v_rot[int(row.nic), jj],
        )
        radii = find_directional_radii(u_rot, v_rot, grid.X_grid, grid.Y_grid, xc, yc, q)
        finite_radii = np.asarray([radii["up"], radii["right"], radii["down"], radii["left"]], dtype=float)
        finite_radii = finite_radii[np.isfinite(finite_radii)]
        if finite_radii.size == 0:
            return _bad_row(row, fnumber)
        radius = float(finite_radii.mean())

        rho_limit = max(min(radius * out_core_fac, rho_max), rho_min) 
        local_limit = rho_limit * local_limit_factor
        local = (grid.X_grid >= xc - local_limit) & (grid.X_grid <= xc + local_limit)
        local &= (grid.Y_grid >= yc - local_limit) & (grid.Y_grid <= yc + local_limit)
        if int(local.sum()) < 10:
            return _bad_row(row, fnumber)

        xloc = grid.X_grid[local]
        yloc = grid.Y_grid[local]
        uloc = u_rot[local]
        vloc = v_rot[local]
        rho = np.sqrt(np.maximum(q[0, 0] * (xloc - xc) ** 2 + 2 * q[0, 1] * (xloc - xc) * (yloc - yc) + q[1, 1] * (yloc - yc) ** 2, 0))
        mask = np.isfinite(rho) & (rho <= rho_limit)
        if int(mask.sum()) < 10:
            return _bad_row(row, fnumber)
        rc, psi0, omega = out_core_param_fit(xloc[mask], yloc[mask], uloc[mask], vloc[mask], xc, yc, q, w)
        if vertical_check is not None and not _has_vertical_profile_to_depth(
            xc_surface=xc,
            yc_surface=yc,
            w_surface=w,
            row=row,
            u_depth=vertical_check["u_depth"],
            v_depth=vertical_check["v_depth"],
            target_depths=vertical_check["target_depths"],
            grid=grid,
            doppio=doppio,
            radius_km=radius_km,
            max_jump_km=vertical_check["max_jump_km"],
            depth_threshold_m=vertical_check["depth_threshold_m"],
        ):
            return _bad_row(row, fnumber)
    except Exception:
        return _bad_row(row, fnumber)

    return _row_with_q(
        {
            "Day": int(row.Day),
            "fnumber": fnumber,
            "nxc": float(row.nxc),
            "nyc": float(row.nyc),
            "nCyc": row.Cyc,
            "nic": int(row.nic),
            "njc": int(row.njc),
            "xc": float(xc),
            "yc": float(yc),
            "w": float(w) * omega_scale,
            "Q": q,
            "Omega0": float(omega0) * omega_scale,
            "Omega": float(omega) * omega_scale,
            "Rc": float(rc),
            "psi0": float(psi0),
            "R": float(radius),
        },
        q,
    )


def fit_surface_file(path: Path, config: PipelineConfig) -> pd.DataFrame:
    fnumber = fnumber_from_outer_avg(path)
    detection_path = detection_path_for_file(path, config)
    if not detection_path.exists():
        raise FileNotFoundError(f"Missing detection partition for {path.name}: {detection_path}")

    detections = pd.read_parquet(detection_path)
    if detections.empty:
        return pd.DataFrame(columns=SURFACE_COLUMNS)

    grid = read_reference_grid(reference_grid_path(config))
    doppio, out_core_param_fit = load_doppio_functions(config)
    settings = config.raw.get("surface_fit", {})
    rho_max = float(settings.get("rho_max_km", 200.0))
    rho_min = float(settings.get("rho_min_km", 30.0))
    local_limit_factor = float(settings.get("local_limit_factor", 3.0))
    out_core_fac = float(settings.get("out_core_fac", 1.75))
    radius_km = float(settings.get("transect_radius_km", 30.0))
    omega_scale = float(settings.get("omega_units_scale", 1e-3))
    require_vertical_profile = bool(settings.get("require_vertical_profile", True))
    z_r = None
    target_depths = None
    if require_vertical_profile:
        z_r = _load_z_r(config)
        target_depths = _target_depths_for_surface_check(z_r, settings)
        if target_depths.size == 0:
            raise ValueError("No valid surface_fit vertical-check target depths found")
    rows = []

    with nc.Dataset(path) as dataset:
        for t, day_value in enumerate(np.asarray(dataset.variables["ocean_time"][:].data, dtype=float) / 86400):
            day = int(round(day_value))
            day_detections = detections.loc[detections["Day"].eq(day)]
            if day_detections.empty:
                continue
            u = dataset["u_eastward"][t, -1, :, :].T
            v = dataset["v_northward"][t, -1, :, :].T
            u_rot, v_rot = rotate_uv(u, v, grid.angle)
            vertical_check = None
            if require_vertical_profile:

                u3d = np.flip(dataset["u_eastward"][t].T.astype(float), axis=2)
                v3d = np.flip(dataset["v_northward"][t].T.astype(float), axis=2)

                vertical_check = {
                    "u_depth": interp_3d_to_reference_depths(u3d, z_r, target_depths),
                    "v_depth": interp_3d_to_reference_depths(v3d, z_r, target_depths),
                    "target_depths": target_depths,
                    "max_jump_km": float(settings.get("vertical_check_max_jump_km", 100.0)),
                    "depth_threshold_m": float(settings.get("vertical_check_depth_m", 50.0)),
                }
            rows.extend(
                _fit_detection_row(
                    row,
                    fnumber,
                    u_rot,
                    v_rot,
                    grid,
                    doppio,
                    out_core_param_fit,
                    radius_km,
                    rho_max,
                    rho_min,
                    out_core_fac,
                    omega_scale,
                    local_limit_factor,
                    vertical_check,
                )
                for row in day_detections.itertuples(index=False)
            )

    return pd.DataFrame(rows, columns=SURFACE_COLUMNS)


def fit_and_write_file(path: Path, config: PipelineConfig) -> Path:
    out_path = surface_output_path(path, config)
    if config.skip_existing and out_path.exists():
        return out_path
    return write_partition(fit_surface_file(path, config), out_path)


def run(config: PipelineConfig) -> None:
    from joblib import Parallel, delayed

    files = find_model_files(config)
    backend = config.raw.get("parallel", {}).get("backend", "process")
    prefer = "threads" if backend == "thread" else "processes"
    written = Parallel(n_jobs=config.workers, prefer=prefer)(
        delayed(fit_and_write_file)(path, config) for path in files
    )
    for path in written:
        print(path)
