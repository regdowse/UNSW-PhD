import numpy as np
import pandas as pd

import selected_case_tools as cases


def _rows():
    return pd.DataFrame({
        "Cyc":["AE"]*5+["CE"], "Eddy":[356]*5+[356],
        "Day":[1,2,3,4,5,1], "topo_plan_ratio":[-2,-1,0,1,2,3],
    })


def test_selection_uses_polarity_and_eddy_id():
    out=cases.select_tracks(_rows(),{"AE":[356]})
    assert len(out)==5 and out.Cyc.eq("AE").all()


def test_regime_threshold_is_factor_two():
    out=cases.add_regimes(_rows().iloc[:5],dominance_factor=2,smooth_window=1)
    assert out.regime.tolist()==["planetary","planetary","mixed","topographic","topographic"]


def test_suggested_days_cover_extremes_and_transition():
    days=cases.suggested_days(_rows().iloc[:5])
    assert days[0]==1 and days[-1]==5 and len(days)==3


def test_compass_endpoint_and_reverse_tilt():
    east=cases._vector_endpoint(0,0,2,90)
    south=cases._vector_endpoint(0,0,2,0,reverse=True)
    np.testing.assert_allclose(east,[2,0],atol=1e-12)
    np.testing.assert_allclose(south,[0,-2],atol=1e-12)
