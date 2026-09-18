import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch
import seacofs_tilt_tools as tilt
from test_surface_pv_footprints import _grid, _snapshot


def profiles():
    return pd.DataFrame(dict(Eddy=[1]*6, Day=[1]*6,
                            Depth=[0, 200, 1000, 1100, -10, 300],
                            w=[-2e-5, -6e-5, -6e-5, -9e-5, 1., np.nan]))


@pytest.mark.parametrize('options', [{}, dict(core_mean=True),
    dict(core_mean=True, averaging='legacy'),
    dict(core_mean=True, surface_method='esp_gaussian')])
def test_peak_propagates_to_all_terms_and_restores_surface(options):
    source = _snapshot()
    expected_source = source.copy(); expected_source['w'] = -6e-5
    expected = tilt.add_pv_gradient_terms(expected_source, _grid(), **options)
    actual = tilt.add_pv_gradient_terms(source, _grid(), vertical=profiles(),
                                       use_max_abs_w=True, **options)
    cols = ['w', 'Ro', 'PV', 'abs_vort'] + [c for c in expected if c.startswith('PV_grad')]
    pd.testing.assert_frame_equal(actual[cols], expected[cols])
    assert actual.w_surface.iloc[0] == source.w.iloc[0]
    assert actual.w_selected_depth_m.iloc[0] == 200
    assert actual.w_profile_n.iloc[0] == 3
    restored = tilt.add_pv_gradient_terms(actual, _grid(), **options)
    original = tilt.add_pv_gradient_terms(source, _grid(), **options)
    pd.testing.assert_frame_equal(restored[cols], original[cols])


def test_grouping_boundaries_missing_and_index():
    surface = pd.concat([_snapshot()]*3); surface['Day'] = [2, 1, 3]
    vertical = pd.concat([profiles(), pd.DataFrame(dict(Eddy=[1], Day=[2], Depth=[1000], w=[8e-5]))])
    result = tilt._select_profile_peak_w(surface, vertical, 1000)
    assert result.index.equals(surface.index)
    np.testing.assert_allclose(result.w, [8e-5, -6e-5, np.nan], equal_nan=True)
    assert result.w_selected_depth_m.iloc[0] == 1000
    assert result.w_profile_n.iloc[2] == 0
    empty = tilt._select_profile_peak_w(surface, vertical.iloc[:0], 1000)
    assert empty.w.isna().all()


def test_default_profile_loading():
    with patch.object(tilt, 'load_vert', return_value=profiles()) as loader:
        tilt.add_pv_gradient_terms(_snapshot(), _grid(), use_max_abs_w=True)
        loader.assert_called_once_with()


def test_peak_cache_is_one_file(tmp_path):
    path = tmp_path / 'surface.parquet'
    with patch.object(pd.DataFrame, 'to_parquet', lambda df, path, **kw: df.to_pickle(path)), patch.object(pd, 'read_parquet', pd.read_pickle):
        tilt.add_pv_gradient_terms(_snapshot(), _grid(), use_cache=False, cache_path=path)
        with pytest.raises(ValueError, match='overwrite the same'):
            tilt.add_pv_gradient_terms(use_cache=True, cache_path=path, use_max_abs_w=True)
        result = tilt.add_pv_gradient_terms(_snapshot(), _grid(), use_cache=False,
            cache_path=path, use_max_abs_w=True, vertical=profiles())
        with patch.object(tilt, 'load_vert', side_effect=AssertionError('cache must not load profiles')):
            cached = tilt.add_pv_gradient_terms(use_cache=True, cache_path=path, use_max_abs_w=True)
        pd.testing.assert_frame_equal(result, cached)
        assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize('options', [dict(source='depth'), dict(depth_following=True), dict(max_depth_m=-1)])
def test_invalid_modes(options):
    with pytest.raises(ValueError):
        tilt.add_pv_gradient_terms(use_max_abs_w=True, **options)
