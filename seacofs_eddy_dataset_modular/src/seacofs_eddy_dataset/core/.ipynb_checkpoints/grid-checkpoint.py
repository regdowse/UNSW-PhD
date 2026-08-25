from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np


@dataclass(frozen=True)
class Grid:
    lon_rho: np.ndarray
    lat_rho: np.ndarray
    x_grid: np.ndarray
    y_grid: np.ndarray
    X_grid: np.ndarray
    Y_grid: np.ndarray
    angle: float
    mask_rho: np.ndarray | None = None
    h: np.ndarray | None = None


def distance_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    earth_radius_km = 6357.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return earth_radius_km * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def read_reference_grid(path: str | Path) -> Grid:
    with nc.Dataset(path) as dataset:
        lon_rho = np.transpose(dataset.variables["lon_rho"], axes=(1, 0))
        lat_rho = np.transpose(dataset.variables["lat_rho"], axes=(1, 0))
        angle = float(dataset.variables["angle"][0, 0])

        mask_rho = None
        h = None
        if "mask_rho" in dataset.variables:
            mask_rho = np.transpose(dataset.variables["mask_rho"], axes=(1, 0))
        if "h" in dataset.variables:
            h = np.transpose(dataset.variables["h"], axes=(1, 0))

    i_mid = lon_rho.shape[0] // 2
    j_mid = lon_rho.shape[1] // 2
    dx = distance_km(
        lat_rho[:-1, j_mid],
        lon_rho[:-1, j_mid],
        lat_rho[1:, j_mid],
        lon_rho[1:, j_mid],
    )
    dy = distance_km(
        lat_rho[i_mid, :-1],
        lon_rho[i_mid, :-1],
        lat_rho[i_mid, 1:],
        lon_rho[i_mid, 1:],
    )
    x_grid = np.insert(np.cumsum(dx), 0, 0)
    y_grid = np.insert(np.cumsum(dy), 0, 0)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")
    return Grid(
        lon_rho=lon_rho,
        lat_rho=lat_rho,
        x_grid=x_grid,
        y_grid=y_grid,
        X_grid=X_grid,
        Y_grid=Y_grid,
        angle=angle,
        mask_rho=mask_rho,
        h=h,
    )


def regular_xy_grid(grid: Grid, resolution_km: float) -> tuple[np.ndarray, np.ndarray]:
    x_new = np.arange(0, grid.x_grid[-1], resolution_km)
    y_new = np.arange(0, grid.y_grid[-1], resolution_km)
    return np.meshgrid(x_new, y_new, indexing="ij")


def fnumber_from_outer_avg(path: str | Path) -> int:
    return int(Path(path).stem.split("_")[-1])
