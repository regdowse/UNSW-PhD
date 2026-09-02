from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from seacofs_eddy_dataset.core.grid import distance_km


@dataclass(frozen=True)
class AvisoGrid:
    longitude: np.ndarray
    latitude: np.ndarray
    x_grid: np.ndarray
    y_grid: np.ndarray
    X_grid: np.ndarray
    Y_grid: np.ndarray
    lon_order: np.ndarray
    lat_order: np.ndarray


def build_grid(longitude, latitude) -> AvisoGrid:
    """Build increasing native Cartesian axes from 1-D AVISO coordinates."""

    longitude = np.asarray(longitude, dtype=float)
    latitude = np.asarray(latitude, dtype=float)
    if longitude.ndim != 1 or latitude.ndim != 1:
        raise ValueError("AVISO longitude and latitude coordinates must be one-dimensional")
    if np.unique(longitude).size != longitude.size or np.unique(latitude).size != latitude.size:
        raise ValueError("AVISO coordinates must not contain duplicates")

    lon_order = np.argsort(longitude)
    lat_order = np.argsort(latitude)
    lon = longitude[lon_order]
    lat = latitude[lat_order]
    if lon.size < 2 or lat.size < 2:
        raise ValueError("AVISO grid needs at least two points on each axis")

    reference_lat = float(np.mean(lat))
    dx = distance_km(reference_lat, lon[:-1], reference_lat, lon[1:])
    reference_lon = float(np.mean(lon))
    dy = distance_km(lat[:-1], reference_lon, lat[1:], reference_lon)
    x_grid = np.insert(np.cumsum(dx), 0, 0.0)
    y_grid = np.insert(np.cumsum(dy), 0, 0.0)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")
    return AvisoGrid(lon, lat, x_grid, y_grid, X_grid, Y_grid, lon_order, lat_order)


def native_velocity(data_array, grid: AvisoGrid) -> np.ndarray:
    """Return an AVISO velocity slice in internal (longitude, latitude) order."""

    values = np.asarray(data_array.transpose("longitude", "latitude").values, dtype=float)
    values = values[np.ix_(grid.lon_order, grid.lat_order)]
    return np.where(np.abs(values) > 1e30, np.nan, values)


def regular_xy_grid(grid: AvisoGrid, resolution_km: float) -> tuple[np.ndarray, np.ndarray]:
    if resolution_km <= 0:
        raise ValueError("interpolation_resolution_km must be positive")
    x_new = np.arange(np.floor(grid.x_grid[-1] / resolution_km) + 1) * resolution_km
    y_new = np.arange(np.floor(grid.y_grid[-1] / resolution_km) + 1) * resolution_km
    return np.meshgrid(x_new, y_new, indexing="ij")
