"""Execute analysis notebook cells against synthetic in-memory data.

No HPC access or Parquet engine required; I/O is isolated in a temporary folder.
The production notebooks are not edited. Synthetic results are never published.
"""
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

ROOT = Path(__file__).resolve().parent


def synthetic_inputs():
    rng = np.random.default_rng(42)
    pieces = []
    for cyc in ('AE', 'CE'):
        for regime in ('planetary', 'topographic'):
            for eddy in range(12):
                n = 40
                days = np.arange(n)
                angle = rng.uniform(-180, 180) + np.cumsum(rng.normal(0, 8, n))
                shear = rng.uniform(-180, 180) + days * .6
                q = pd.DataFrame(dict(Cyc=cyc, Eddy=f'{cyc}_{regime}_{eddy}', Day=days,
                    TiltDis=np.maximum(6, 15 + np.cumsum(rng.normal(0, .3, n))),
                    TiltDir=(90-angle-shear)%360, h=4000., water_depth_m=4500.,
                    lat=rng.uniform(-38,-26), PV_grad_plan_mag=1.,
                    PV_grad_topo_mag=.1 if regime=='planetary' else 10.,
                    PV_grad_mag=1., PV_grad_theta=20., PV_grad_coherence=.8))
                for family in ('clim','full'):
                    for axis, func in [('east',np.cos),('north',np.sin)]:
                        q[f'{family}_200_{axis}_ms']=.08*func(np.deg2rad(shear))
                        q[f'{family}_500_{axis}_ms']=.032*func(np.deg2rad(shear))
                        q[f'{family}_surface_{axis}_ms']=.1*func(np.deg2rad(shear+10))
                pieces.append(q)
    return pd.concat(pieces,ignore_index=True)


class NotebookSmoke(unittest.TestCase):
    def test_notebook_syntax(self):
        for path in ROOT.glob('*.ipynb'):
            for i,cell in enumerate(json.loads(path.read_text())['cells']):
                if cell['cell_type']=='code':
                    compile(''.join(cell['source']),f'{path.name}:{i}','exec')

    def run_notebook(self, raw, destination):
        notebook=json.loads((ROOT/'01_shear_tilt_climatology.ipynb').read_text())
        tmp=Path(destination)
        (tmp/'climatology_inputs.parquet').touch()
        (tmp/'input_metadata.json').write_text('{}')
        namespace={}
        old=os.getcwd()
        try:
            os.chdir(ROOT)
            with patch.object(pd,'read_parquet',return_value=raw), contextlib.redirect_stdout(io.StringIO()):
                for i,cell in enumerate(notebook['cells']):
                    if cell['cell_type']!='code': continue
                    source=''.join(cell['source'])
                    source=source.replace("Path('/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/shear_tilt_climatology')",repr(tmp))
                    # Path repr is PosixPath(...); make portable constructor available.
                    namespace['PosixPath']=Path
                    source=source.replace('BOOTSTRAPS = 500','BOOTSTRAPS = 15')
                    source=source.replace('MIN_EDDIES = 20','MIN_EDDIES = 3')
                    source=source.replace('MIN_ROWS = 50','MIN_ROWS = 5')
                    exec(compile(source,f'notebook-cell-{i}','exec'),namespace)
                    namespace['display']=lambda *args,**kwargs:None
        finally:
            os.chdir(old)
        return namespace

    def test_full_analysis_and_sparse_population(self):
        # Optional QA destination lets a reviewer inspect synthetic figures.
        destination=os.environ.get('SHEAR_CLIMATOLOGY_QA_DIR')
        with tempfile.TemporaryDirectory() as scratch:
            target=Path(destination or scratch); target.mkdir(parents=True,exist_ok=True)
            ns=self.run_notebook(synthetic_inputs(),target)
            self.assertGreater(len(ns['pairs']),0)
            self.assertTrue((ns['OUTPUT']/'settings.json').exists())
            self.assertGreater(len(list(ns['OUTPUT'].glob('*.png'))),15)
            # Missing CE and no eligible lag pairs must give empty plots/tables, not crash.
            q=synthetic_inputs().query("Cyc == 'AE' and Day % 2 == 0").copy()
            # Test sparse outputs with minimal smoothing sufficient for isolated days.
            # Direct geometry/pair plotting covers no-pair condition because production
            # seven-day smoothing correctly rejects isolated input days entirely.
            import climatology_tools as ct
            g=ct.relative_geometry(ct.classify_regimes(q,window=1,min_periods=1))
            p=ct.lagged_pairs(g,3)
            self.assertTrue(p.empty)
            s=ct.cluster_summary(p,['growth_km_day'],groups=ct.GROUPS+['angle_bin'])
            fig=ct.plot_response(s,'growth_km_day')
            import matplotlib.pyplot as plt
            plt.close(fig)


if __name__=='__main__':
    unittest.main()
