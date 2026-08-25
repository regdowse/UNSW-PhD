from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def clean_fill_values(array, fill_threshold: float = 1e30) -> np.ndarray:
    values = np.asarray(array, dtype=float)
    return np.where(np.abs(values) > fill_threshold, np.nan, values)


def rotate_uv(u, v, angle: float) -> tuple[np.ndarray, np.ndarray]:
    u = clean_fill_values(u)
    v = clean_fill_values(v)
    u_rot = v * np.sin(angle) + u * np.cos(angle)
    v_rot = v * np.cos(angle) - u * np.sin(angle)
    return u_rot, v_rot


def interpolate_uv_to_grid(u, v, x_grid, y_grid, X_new, Y_new, angle: float):
    u_rot, v_rot = rotate_uv(u, v, angle)
    shape_new = X_new.shape
    new_points = np.column_stack((X_new.ravel(), Y_new.ravel()))

    interp_u = RegularGridInterpolator(
        (x_grid, y_grid), u_rot, method="linear", bounds_error=False, fill_value=np.nan
    )
    interp_v = RegularGridInterpolator(
        (x_grid, y_grid), v_rot, method="linear", bounds_error=False, fill_value=np.nan
    )
    return interp_u(new_points).reshape(shape_new), interp_v(new_points).reshape(shape_new)

