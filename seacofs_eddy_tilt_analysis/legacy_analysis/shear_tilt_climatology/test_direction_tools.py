"""Scientific invariants and end-to-end execution on synthetic data only."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import climatology_tools as ct
import direction_tools as dt
from test_notebook_smoke import synthetic_inputs

ROOT=Path(__file__).resolve().parent

class DirectionTests(unittest.TestCase):
    def geometry(self):
        raw=synthetic_inputs()
        return ct.relative_geometry(ct.classify_regimes(raw))

    def test_left_right_independent_of_geographic_bearing(self):
        q=self.geometry(); valid=np.isfinite(q.left_fraction)
        q=q.loc[valid].copy()
        q['shear_angle']=np.resize(np.arange(0,360,45),len(q))
        q['tilt_angle']=q.shear_angle+np.where(q.Cyc.eq('AE'),90,-90)
        q['left_fraction']=np.sin(np.deg2rad(q.tilt_angle-q.shear_angle))
        p=dt.populations(q)
        self.assertTrue(np.allclose(p.expected_side,1))
        self.assertEqual(set(p.sector),set(range(8)))
        # A fixed eastward tilt cannot mimic the same cross-shear sign in all directions.
        fixed=np.sin(np.deg2rad(0-np.arange(0,360,45)))
        self.assertTrue((fixed>.1).any() and (fixed<-.1).any())

    def test_all_includes_unknown_but_never_invalid_geometry(self):
        q=self.geometry(); q['regime']='unknown'
        q.loc[q.index[0],'left_fraction']=np.nan
        p=dt.populations(q)
        self.assertEqual(set(p.population),{'all'})
        self.assertEqual(len(p),np.isfinite(q.left_fraction).sum())
        w=dt.shear_windows(q,14)
        self.assertEqual(set(w.population),{'all'})

    def test_fixed_duration_and_gap_safety(self):
        q=self.geometry(); q=q.loc[(q.Cyc=='AE')&(q.Eddy=='AE_planetary_0')].copy()
        w=dt.shear_windows(q,14)
        self.assertEqual(len(w.query("population == 'planetary'")),2)
        broken=q.loc[q.Day.ne(20)]
        pairs=dt.turning_pairs(broken,7)
        self.assertFalse(((pairs.Day<20)&(pairs.Day+7>20)).any())
        w=dt.shear_windows(broken,14)
        self.assertFalse(((w.Day<20)&(w.Day+13>20)).any())
        q['shear_angle']=45.
        self.assertTrue(np.allclose(dt.shear_windows(q,14).circular_variance,0))

    def test_balancing_rejects_missing_directions(self):
        p=dt.populations(self.geometry())
        p=p.loc[p.sector.ne(7)]
        result=dt.sector_balanced(p,n_boot=10,min_eddies=1,min_rows=1)
        self.assertTrue(result['mean'].isna().all())
        self.assertTrue(result.ci_low.isna().all())

    def test_turning_new_vs_frozen_shear(self):
        q=self.geometry()
        q['tilt_angle']=q.shear_angle+np.where(q.Cyc.eq('AE'),90,-90)
        q['angle_deg']=np.where(q.Cyc.eq('AE'),90.,-90.)
        q['left_fraction']=np.where(q.Cyc.eq('AE'),1.,-1.)
        p=dt.turning_pairs(q,7)
        self.assertTrue((p.new_minus_old_expected>0).all())
        self.assertTrue(np.allclose(p.both_expected_side,1))

    def test_notebook_end_to_end(self):
        notebook=json.loads((ROOT/'02_shear_direction_preference.ipynb').read_text())
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(os.environ.get('DIRECTION_QA_DIR',tmp)); target.mkdir(parents=True,exist_ok=True)
            (target/'climatology_inputs.parquet').touch()
            (target/'input_metadata.json').write_text('{}')
            ns={'PosixPath':Path}; old=os.getcwd()
            try:
                os.chdir(ROOT)
                with patch.object(pd,'read_parquet',return_value=synthetic_inputs()),contextlib.redirect_stdout(io.StringIO()):
                    for i,c in enumerate(notebook['cells']):
                        if c['cell_type']!='code': continue
                        s=''.join(c['source']).replace("Path('/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/shear_tilt_climatology')",repr(target))
                        s=s.replace('BOOTSTRAPS = 500','BOOTSTRAPS = 10').replace('MIN_EDDIES = 20','MIN_EDDIES = 3').replace('MIN_ROWS = 50','MIN_ROWS = 5')
                        exec(compile(s,f'direction-cell-{i}','exec'),ns)
                        ns['display']=lambda *a,**k:None
            finally:
                os.chdir(old); plt.close('all')
            self.assertGreater(len(ns['turn_results']),0)
            self.assertEqual(set(ns['sector_results'].population),{'all','planetary','topographic'})
            self.assertTrue((ns['OUTPUT']/'primary_sector_scorecard.csv').exists())
            self.assertGreater(len(list(ns['OUTPUT'].glob('*.png'))),15)

if __name__=='__main__': unittest.main()
