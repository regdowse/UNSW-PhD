import numpy as np
import pandas as pd

import seacofs_tilt_tools as tilt


def _grid():
    x = np.arange(9, dtype=float)
    y = np.arange(9, dtype=float)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    return tilt.Grid(
        lon_rho=xx,
        lat_rho=yy,
        mask_rho=np.ones_like(xx, dtype=bool),
        h=2000.0 + 20.0 * xx + 10.0 * yy,
        f=-8e-5 + 2e-9 * yy,
        angle=0.0,
        z_r=np.empty((9, 9, 0)),
        x_grid=x,
        y_grid=y,
        X_grid=xx,
        Y_grid=yy,
    )


def test_vertical_cell_weights_are_thickness_based():
    weights = tilt._vertical_cell_weights([100.0, 500.0, 900.0], 1000.0)
    np.testing.assert_allclose(weights, [300.0, 400.0, 300.0])


def test_depth_following_returns_snapshot_and_depth_tables():
    snapshots = pd.DataFrame({
        "Eddy": [7], "Day": [12], "Cyc": ["CE"],
        "TiltDis": [20.0], "TiltDir": [180.0],
    })
    vertical = pd.DataFrame({
        "Eddy": [7, 7, 7], "Day": [12, 12, 12],
        "Depth": [100.0, 500.0, 900.0],
        "xc": [3.0, 4.0, 5.0], "yc": [4.0, 4.0, 4.0],
        "Rc": [2.0, 2.0, 2.0],
        "q11": [1.0, 1.0, 1.0], "q12": [0.0, 0.0, 0.0],
        "q22": [1.0, 1.0, 1.0],
        "w": [-2e-5, -3e-5, -4e-5],
    })

    snapshot, depth = tilt.add_pv_gradient_terms(
        snapshots, _grid(), core_mean=True, depth_following=True,
        vertical=vertical, max_depth_m=1000.0,
    )

    assert len(snapshot) == 1
    assert len(depth) == 3
    np.testing.assert_allclose(depth["vertical_weight_m"], [300.0, 400.0, 300.0])
    assert snapshot.loc[0, "PV_depth_n"] == 3
    assert snapshot.loc[0, "PV_vertical_coverage_fraction"] == 1.0
    assert np.isfinite(depth["PV_grad_mag"]).all()
    assert np.isfinite(snapshot.loc[0, "PV_grad_theta"])
    np.testing.assert_allclose(
        snapshot.loc[0, "PV_grad_x"],
        snapshot.loc[0, "PV_grad_plan_x"] + snapshot.loc[0, "PV_grad_topo_x"],
    )
    np.testing.assert_allclose(
        snapshot.loc[0, "PV_grad_y"],
        snapshot.loc[0, "PV_grad_plan_y"] + snapshot.loc[0, "PV_grad_topo_y"],
    )


def test_depth_following_requires_core_mean_and_vertical_table():
    snapshots = pd.DataFrame({
        "Eddy": [1], "Day": [1], "TiltDis": [1.0], "TiltDir": [0.0]
    })
    try:
        tilt.add_pv_gradient_terms(snapshots, _grid(), depth_following=True)
    except ValueError as exc:
        assert "core_mean=True" in str(exc)
    else:
        raise AssertionError("Expected depth-following core-mean validation")


def test_default_surface_api_still_returns_one_dataframe():
    snapshots = pd.DataFrame({
        "Eddy": [1], "Day": [1], "Cyc": ["AE"], "fname": ["x_01461"],
        "ic": [4], "jc": [4], "w": [-2e-5],
        "TiltDis": [10.0], "TiltDir": [20.0],
    })
    result = tilt.add_pv_gradient_terms(snapshots, _grid(), core_mean=False)
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1
    assert np.isfinite(result.loc[0, "PV_grad_mag"])
