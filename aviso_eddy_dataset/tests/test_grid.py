import numpy as np
import xarray as xr

from aviso_eddy_dataset.grid import build_grid, native_velocity, regular_xy_grid


def test_grid_sorts_coordinates_and_velocity_consistently():
    longitude = np.array([151.0, 150.0])
    latitude = np.array([-30.0, -31.0, -32.0])
    grid = build_grid(longitude, latitude)
    velocity = xr.DataArray(
        np.array([[1, 2], [3, 4], [5, 6]]),
        dims=("latitude", "longitude"),
        coords={"latitude": latitude, "longitude": longitude},
    )
    result = native_velocity(velocity, grid)
    assert np.all(np.diff(grid.longitude) > 0)
    assert np.all(np.diff(grid.latitude) > 0)
    np.testing.assert_array_equal(result, [[6, 4, 2], [5, 3, 1]])


def test_regular_grid_uses_requested_resolution():
    grid = build_grid([150.0, 150.2], [-32.0, -31.8])
    X, Y = regular_xy_grid(grid, 1.0)
    assert X.shape == Y.shape
    np.testing.assert_allclose(np.diff(X[:, 0]), 1.0)
    np.testing.assert_allclose(np.diff(Y[0, :]), 1.0)
