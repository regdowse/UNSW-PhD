import numpy as np
import pandas as pd

import seacofs_tilt_tools as tilt


def _grid(seamount=False):
    x = np.arange(-10, 11, dtype=float)
    y = np.arange(-10, 11, dtype=float)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    if seamount:
        h = 3000.0 - 1800.0 * np.exp(-(xx**2 + yy**2) / (2 * 2.5**2))
    else:
        h = 2500.0 + 20.0 * xx + 5.0 * yy
    return tilt.Grid(
        lon_rho=xx, lat_rho=yy, mask_rho=np.ones_like(xx, dtype=bool),
        h=h, f=-8e-5 + 2e-9 * yy, angle=0.0,
        z_r=np.empty((21, 21, 0)), x_grid=x, y_grid=y,
        X_grid=xx, Y_grid=yy,
    )


def _snapshot():
    return pd.DataFrame({
        "Eddy": [1], "Day": [1], "Cyc": ["AE"], "fname": ["x_01461"],
        "ic": [10], "jc": [10], "xc": [0.0], "yc": [0.0], "Rc": [10.0],
        "q11": [1.0], "q12": [0.0], "q22": [1.0], "w": [-2e-5],
        "TiltDis": [10.0], "TiltDir": [20.0],
    })


def test_frac_scales_axes_and_annulus_excludes_inner_core():
    row = next(_snapshot().itertuples(index=False))
    grid = _grid()
    full = set(zip(*tilt.core_grid_indices(row, grid, frac=1.0)))
    half = set(zip(*tilt.core_grid_indices(row, grid, frac=0.5)))
    annulus = set(zip(*tilt.core_grid_indices(row, grid, frac=1.0, inner_frac=0.5)))
    assert half < full
    assert annulus < full
    assert half | annulus == full
    assert (10, 10) not in annulus


def test_seamount_exposure_survives_vector_cancellation():
    result = tilt.add_pv_gradient_terms(
        _snapshot(), _grid(seamount=True), core_mean=True,
        frac=1.0, averaging="nonlinear",
    ).iloc[0]
    assert result.PV_grad_topo_mean_local_mag > result.PV_grad_topo_mag
    assert result.PV_grad_topo_coherence < 0.25
    assert result.PV_grad_topo_p90_local_mag > 0


def test_legacy_standard_columns_match_historical_formula():
    source = _snapshot()
    grid = _grid()
    result = tilt.add_pv_gradient_terms(
        source, grid, core_mean=True, frac=0.5, averaging="legacy"
    ).iloc[0]
    h = tilt.compute_core_mean(source, grid, fixed_field=grid.h, colname="h", frac=0.5)
    gx, gy = tilt.phys_grad(grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho)
    de = np.cos(grid.angle) * gx - np.sin(grid.angle) * gy
    dn = np.sin(grid.angle) * gx + np.cos(grid.angle) * gy
    mx = tilt.compute_core_mean(source, grid, fixed_field=de, colname="dhdx", frac=0.5)
    my = tilt.compute_core_mean(source, grid, fixed_field=dn, colname="dhdy", frac=0.5)
    omega_f = source.w.iloc[0] + grid.f[10, 10]
    np.testing.assert_allclose(result.h, h.h.iloc[0])
    np.testing.assert_allclose(result.PV_grad_topo_x, -omega_f * mx.dhdx.iloc[0] / h.h.iloc[0]**2)
    np.testing.assert_allclose(result.PV_grad_topo_y, -omega_f * my.dhdy.iloc[0] / h.h.iloc[0]**2)


def test_surface_only_options_are_rejected_for_cached_sources():
    try:
        tilt.add_pv_gradient_terms(source="depth", frac=0.5)
    except ValueError as exc:
        assert "source='original'" in str(exc)
    else:
        raise AssertionError("Expected surface-only option validation")


def test_invalid_annulus_is_rejected():
    for outer, inner in [(0, 0), (1, -0.1), (0.5, 0.5), (0.5, 0.8)]:
        try:
            tilt.add_pv_gradient_terms(
                _snapshot(), _grid(), core_mean=True,
                frac=outer, inner_frac=inner,
            )
        except ValueError:
            pass
        else:
            raise AssertionError((outer, inner))


def test_esp_gaussian_uses_saved_w_as_central_vorticity():
    result = tilt.add_pv_gradient_terms(
        _snapshot(), _grid(), core_mean=True, frac=0.1,
        surface_method="esp_gaussian",
    ).iloc[0]
    np.testing.assert_allclose(result.zeta_mean, _snapshot().w.iloc[0])
    assert result.PV_footprint_n == 1
    assert result.PV_weight_effective_n == 1


def test_esp_gaussian_internal_gradient_has_zero_net_in_symmetric_footprint():
    grid = _grid()
    grid.h[:] = 2500.0
    result = tilt.add_pv_gradient_terms(
        _snapshot(), grid, core_mean=True, frac=1.5,
        surface_method="esp_gaussian",
    ).iloc[0]
    assert result.PV_grad_eddy_mean_local_mag > 0
    assert result.PV_grad_eddy_mag < result.PV_grad_eddy_mean_local_mag * 1e-12
    np.testing.assert_allclose(result.PV_grad_full_x, result.PV_grad_x, atol=1e-25)
    np.testing.assert_allclose(result.PV_grad_full_y, result.PV_grad_y, rtol=1e-10)


def test_gaussian_weighting_has_smaller_effective_than_raw_sample_size():
    result = tilt.add_pv_gradient_terms(
        _snapshot(), _grid(seamount=True), core_mean=True, frac=2.0,
        surface_method="esp_gaussian",
    ).iloc[0]
    assert 0 < result.PV_weight_effective_n < result.PV_footprint_n
    assert result.pv_surface_method == "esp_gaussian"


def test_legacy_rejects_esp_gaussian_surface_method():
    try:
        tilt.add_pv_gradient_terms(
            _snapshot(), _grid(), core_mean=True,
            averaging="legacy", surface_method="esp_gaussian",
        )
    except ValueError as exc:
        assert "legacy" in str(exc)
    else:
        raise AssertionError("Expected incompatible surface-method validation")
