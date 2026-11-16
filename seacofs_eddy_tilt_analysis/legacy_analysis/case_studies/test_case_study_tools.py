import numpy as np
import pandas as pd

from case_study_tools import PVAlignmentConfig, add_pv_alignment_diagnostics


def _direction_rows():
    return pd.DataFrame(
        {
            "Eddy": [1, 2],
            "Day": [0, 0],
            "Cyc": ["AE", "CE"],
            "h": [4000.0, 4000.0],
            "TiltDis": [10.0, 10.0],
            "TiltDir": [30.0, 210.0],
            "PV_grad_theta": [30.0, 30.0],
            "topo_plan_ratio": [-2.0, -2.0],
        }
    )


def test_signed_pv_targets_are_ae_along_and_ce_opposite():
    config = PVAlignmentConfig(smooth_window=1, min_periods=1)
    result = add_pv_alignment_diagnostics(_direction_rows(), config)

    np.testing.assert_allclose(result["expected_PV_theta"], [30.0, 210.0])
    np.testing.assert_allclose(result["expected_direction_error"], [0.0, 0.0])
    assert result["expected_match"].tolist() == [True, True]
