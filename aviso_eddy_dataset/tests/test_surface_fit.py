from types import SimpleNamespace

import numpy as np

from aviso_eddy_dataset.grid import build_grid
from aviso_eddy_dataset.surface_fit import fit_detection


def test_doppio_receives_native_grid_transects():
    longitude = np.linspace(150.0, 151.0, 11)
    latitude = np.linspace(-32.0, -31.0, 11)
    grid = build_grid(longitude, latitude)
    xc = grid.x_grid[5]
    yc = grid.y_grid[5]
    dx = grid.X_grid - xc
    dy = grid.Y_grid - yc
    radius = np.hypot(dx, dy)
    speed = radius * np.exp(-(radius / 30.0) ** 2)
    u = np.divide(-dy * speed, radius, out=np.zeros_like(radius), where=radius > 0)
    v = np.divide(dx * speed, radius, out=np.zeros_like(radius), where=radius > 0)
    observed = {}

    def doppio(x1, y1, u1, v1, x2, y2, u2, v2):
        observed["x1"] = x1.copy()
        observed["u1"] = u1.copy()
        return xc, yc, 1.0, np.eye(2), 0.01

    def outer(*args):
        return 20.0, 1.0, 0.02

    row = SimpleNamespace(
        Day=1,
        Date=np.datetime64("1990-01-02"),
        source_file="source.nc",
        nxc=xc,
        nyc=yc,
        Cyc="AE",
        nic=5,
        njc=5,
    )
    result = fit_detection(
        row,
        u,
        v,
        grid,
        doppio,
        outer,
        {
            "transect_radius_km": 40.0,
            "rho_min_km": 30.0,
            "rho_max_km": 100.0,
            "local_limit_factor": 3.0,
            "out_core_fac": 1.75,
            "omega_units_scale": 1.0,
        },
    )
    assert len(observed["x1"]) < len(grid.x_grid)
    native_indices = np.flatnonzero((grid.x_grid >= xc - 40.0) & (grid.x_grid <= xc + 40.0))
    np.testing.assert_array_equal(observed["x1"], grid.x_grid[native_indices])
    np.testing.assert_array_equal(observed["u1"], u[native_indices, 5])
    assert np.isfinite(result["xc"])
