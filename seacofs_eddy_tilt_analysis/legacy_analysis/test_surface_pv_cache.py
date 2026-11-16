from unittest.mock import patch
import pandas as pd
import unittest
import tempfile
from pathlib import Path
import seacofs_tilt_tools as tilt
from test_surface_pv_footprints import _grid, _snapshot



class CacheTests(unittest.TestCase):
    def setUp(self):
        self.folder=tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.tmp_path=Path(self.folder.name)
        # File lifecycle test with pickle transport: no local Parquet engine.
        for patcher in [patch.object(pd.DataFrame,'to_parquet',lambda df,path,**kw:df.to_pickle(path)),
                        patch.object(pd,'read_parquet',pd.read_pickle)]:
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_save_load_refresh(self):
        tmp_path=self.tmp_path
        path=tmp_path/'surface.parquet'
        options=dict(core_mean=True,surface_method='esp_gaussian',cache_path=path)
        original=tilt.add_pv_gradient_terms(_snapshot(),_grid(),use_cache=False,**options)
        with patch.object(tilt,'phys_grad',side_effect=AssertionError('must not calculate')):
            loaded=tilt.add_pv_gradient_terms(use_cache=True,**options)
        pd.testing.assert_frame_equal(original,loaded)
        changed=_snapshot();changed['w'] *= 2
        updated=tilt.add_pv_gradient_terms(changed,_grid(),use_cache=False,**options)
        assert updated.Ro.iloc[0] != original.Ro.iloc[0]
        pd.testing.assert_frame_equal(updated,tilt.add_pv_gradient_terms(use_cache=True,**options))


    def test_missing_and_mismatched_settings(self):
        tmp_path=self.tmp_path
        path=tmp_path/'surface.parquet'
        with self.assertRaisesRegex(FileNotFoundError,'use_cache=False'):
            tilt.add_pv_gradient_terms(use_cache=True,cache_path=path)
        tilt.add_pv_gradient_terms(_snapshot(),_grid(),use_cache=False,cache_path=path)
        with self.assertRaisesRegex(ValueError,'settings differ'):
            tilt.add_pv_gradient_terms(use_cache=True,cache_path=path,core_mean=True)
        with self.assertRaisesRegex(ValueError,'surface calculation'):
            tilt.add_pv_gradient_terms(use_cache=True,source='depth')


    def test_failed_write_preserves_cache(self):
        tmp_path=self.tmp_path
        path=tmp_path/'surface.parquet'
        tilt.add_pv_gradient_terms(_snapshot(),_grid(),use_cache=False,cache_path=path)
        old=path.read_bytes()
        with patch.object(pd.DataFrame,'to_parquet',side_effect=OSError('write failure')):
            with self.assertRaises(OSError):
                tilt.add_pv_gradient_terms(_snapshot(),_grid(),use_cache=False,cache_path=path)
        assert path.read_bytes()==old
        assert list(tmp_path.iterdir())==[path]


    def test_omitted_flag_does_not_save(self):
        tmp_path=self.tmp_path
        path=tmp_path/'not_created.parquet'
        tilt.add_pv_gradient_terms(_snapshot(),_grid(),cache_path=path)
        assert not path.exists()
