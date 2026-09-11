import numpy as np
import pandas as pd

from surface_pv_footprint_sensitivity import footprint_tools as ft


def _data():
    rows = []
    for cyc, eddy, tilt_dir, pv_dir in [("AE", 1, 0.0, 0.0), ("CE", 2, 180.0, 0.0)]:
        for day in (1, 2):
            for footprint, kind, outer, inner in [
                ("filled_1", "Filled", 1.0, 0.0),
                ("annulus_0.5_1", "Annulus", 1.0, 0.5),
            ]:
                rows.append({
                    "Eddy": eddy, "Day": day, "Cyc": cyc,
                    "TiltDis": 20.0, "TiltDir": tilt_dir,
                    "PV_grad_theta": pv_dir,
                    "PV_grad_topo_mag": 1.0,
                    "PV_grad_plan_mag": 1.0,
                    "PV_grad_topo_mean_local_mag": 2.0,
                    "PV_grad_plan_mean_local_mag": 1.0,
                    "PV_grad_topo_p90_local_mag": 3.0,
                    "PV_grad_topo_coherence": 0.5,
                    "topo_plan_ratio": 0.0,
                    "footprint": footprint, "footprint_kind": kind,
                    "ellipse_frac": outer, "ellipse_inner_frac": inner,
                })
    return pd.DataFrame(rows)


def test_default_experiment_contains_filled_annulus_and_legacy():
    specs = ft.default_footprints()
    assert set(specs.kind) == {"Filled", "Annulus", "Legacy"}
    assert not specs.footprint.duplicated().any()


def test_analysis_terms_preserve_cancellation_and_polarity_convention():
    data = ft.add_analysis_terms(_data())
    np.testing.assert_allclose(data.topo_cancellation_gap, np.log(2))
    np.testing.assert_allclose(data.preference_error_deg, 0)
    np.testing.assert_allclose(data.preference_alignment, 1)


def test_scorecard_is_eddy_equal():
    data = ft.add_analysis_terms(_data())
    score = ft.footprint_scorecard(data)
    assert score.eddies.eq(1).all()
    np.testing.assert_allclose(score.median_error_deg, 0)


def test_reference_pairing_is_many_to_one():
    data = ft.add_analysis_terms(_data())
    paired = ft.paired_with_reference(data)
    assert len(paired) == len(data)
    assert "reference_PV_grad_topo_coherence" in paired
