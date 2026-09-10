import numpy as np
import pandas as pd

from depth_resolved_pv_tilt import depth_pv_tools as dpt


def sample_depth():
    rows = []
    for eddy, cyc in [(1, "AE"), (2, "CE")]:
        for day in (1, 2):
            for depth, x, ratio, theta in [(0, 0, -1.0, 0), (500, 3, 0.0, 45), (1000, 6, 1.0, 180)]:
                rows.append({
                    "Eddy": eddy, "Day": day, "Depth": depth, "Cyc": cyc,
                    "TiltDis": 20.0, "TiltDir": 0.0, "xc": x, "yc": 0.0,
                    "Rc": 10.0, "PV_grad_plan_mag": 1.0,
                    "PV_grad_topo_mag": np.exp(ratio), "PV_grad_mag": 2.0,
                    "PV_grad_plan_theta": 0.0, "PV_grad_topo_theta": theta,
                    "PV_grad_theta": theta, "topo_plan_ratio": ratio,
                    "h": 2000 + depth, "dhdx": 0.01, "dhdy": 0.0,
                    "Ro": 0.2, "w": 1e-5, "PV": 1e-8,
                })
    return pd.DataFrame(rows)


def test_surface_differences_are_paired_by_snapshot():
    out = dpt.add_surface_differences(sample_depth(), dominance_factor=2)
    bottom = out[out.Depth.eq(1000)]
    np.testing.assert_allclose(bottom.spine_displacement_km, 6)
    np.testing.assert_allclose(bottom.delta_h_m, 1000)
    assert bottom.regime_changed.all()
    assert bottom.surface_regime.astype(str).eq("Planetary").all()
    assert bottom.regime.astype(str).eq("Topographic").all()


def test_expected_alignment_is_polarity_aware():
    df = pd.DataFrame({
        "Cyc": ["AE", "CE"], "TiltDir": [10.0, 190.0],
        "PV_grad_theta": [10.0, 10.0],
    })
    out = dpt.add_expected_alignment(df)
    np.testing.assert_allclose(out.preference_error_deg, 0)
    np.testing.assert_allclose(out.preference_alignment, 1)


def test_matched_depth_rows_drops_incomplete_snapshots():
    df = sample_depth()
    df = df[~((df.Eddy == 2) & (df.Day == 2) & (df.Depth == 500))]
    out = dpt.matched_depth_rows(df, [0, 500, 1000])
    assert len(out) == 9
    assert not ((out.Eddy == 2) & (out.Day == 2)).any()


def test_regime_matrix_is_row_normalised():
    df = dpt.add_surface_differences(sample_depth())
    matrix = dpt.regime_transition_matrix(df, 1000)
    np.testing.assert_allclose(matrix.loc["Planetary"].sum(), 1)


def test_nearest_depths_are_unique():
    result = dpt.nearest_cached_depths(sample_depth(), [0, 10, 500, 900, 1000])
    np.testing.assert_allclose(result, [0, 500, 1000])
