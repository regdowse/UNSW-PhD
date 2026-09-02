import numpy as np

from aviso_eddy_dataset.detection import interpolate_for_nencioli
from aviso_eddy_dataset.grid import build_grid, regular_xy_grid


def test_detection_interpolation_preserves_constant_unrotated_velocity():
    grid = build_grid([150.0, 150.1, 150.2], [-32.0, -31.9, -31.8])
    X, Y = regular_xy_grid(grid, 2.0)
    u = np.full(grid.X_grid.shape, 0.25)
    v = np.full(grid.X_grid.shape, -0.10)
    ui, vi = interpolate_for_nencioli(u, v, grid, X, Y)
    np.testing.assert_allclose(ui[np.isfinite(ui)], 0.25)
    np.testing.assert_allclose(vi[np.isfinite(vi)], -0.10)
