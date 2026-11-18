"""Surface-fit environmental PV. Build explicitly; subsequent runs validate and load."""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd

import pv_kernels as kernels
import pv_plots

KEYS = ['Eddy', 'Day']
PREFIXES = ['PV_grad_plan', 'PV_grad_topo', 'PV_grad', 'PV_grad_eddy', 'PV_grad_full']


def identity(path):
    p = Path(path)
    digest = hashlib.sha256()
    with p.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return dict(path=str(p.resolve()), bytes=p.stat().st_size, sha256=digest.hexdigest())


def check_keys(table, keys=KEYS):
    if table[keys].isna().any().any() or table.duplicated(keys).any():
        raise ValueError(f'Missing/duplicate keys: {keys}')


def load_config(path):
    path = Path(path).resolve()
    c = json.loads(path.read_text())
    for section in ['inputs']:
        for k, value in c[section].items():
            p = Path(value).expanduser()
            c[section][k] = str(p if p.is_absolute() else path.parent/p)
    for k in ['product_dir', 'output_dir','figure_dir']:
        p = Path(c[k]).expanduser()
        c[k] = str(p if p.is_absolute() else path.parent/p)
    a = c['calculation']
    if a['vorticity_amplitude'] != 'surface' or a['footprint_fraction'] != 1:
        raise ValueError('This workflow requires surface amplitude and FRAC=1')
    if a['chunk_rows'] < 1 or a['peak_vorticity_max_depth_m'] <= 0:
        raise ValueError('Invalid chunk size or peak-depth range')
    return c


def load_grid(path):
    with nc.Dataset(path) as ds:
        g = SimpleNamespace(**{k: np.asarray(ds[k][:]).T for k in
                               ['lon_rho', 'lat_rho', 'mask_rho', 'h', 'f']})
        g.angle = float(ds['angle'][0, 0])
    i, j = np.array(g.h.shape)//2
    g.x_grid = np.r_[0, np.cumsum(kernels.distance_km(g.lat_rho[:-1,j],g.lon_rho[:-1,j],
                                                     g.lat_rho[1:,j],g.lon_rho[1:,j]))]
    g.y_grid = np.r_[0, np.cumsum(kernels.distance_km(g.lat_rho[i,:-1],g.lon_rho[i,:-1],
                                                     g.lat_rho[i,1:],g.lon_rho[i,1:]))]
    g.X_grid, g.Y_grid = np.meshgrid(g.x_grid, g.y_grid, indexing='ij')
    return g


def calculate(surface, vertical, grid, settings):
    """One derived row per source day; no tilt outcomes or smoothed variables stored."""
    check_keys(surface)
    check_keys(vertical, KEYS+['Depth'])
    required = KEYS+['xc', 'yc', 'Rc', 'q11', 'q12', 'q22', 'w']
    source = surface[required].copy().reset_index(drop=True)
    chunks = []
    size = settings['chunk_rows']
    # Grid derivatives are invariant across days: calculate once.
    gradients = pv_plots.environmental_grid_gradients(grid)
    for start in range(0, len(source), size):
        chunk = source.iloc[start:start+size]
        result = kernels._surface_pv_footprint_statistics(
            chunk, grid, frac=1., inner_frac=0., surface_method='esp_gaussian',
            dh_dE=gradients['dh_e'], dh_dN=gradients['dh_n'], beta=gradients['beta'])
        chunks.append(result)
        print(f'PV footprints: {min(start+size,len(source)):,}/{len(source):,}', flush=True)
    result = pd.concat(chunks).reindex(source.index)
    for prefix in PREFIXES:
        result[prefix+'_mag'] = np.hypot(result[prefix+'_x'], result[prefix+'_y'])
        result[prefix+'_theta'] = kernels.bearing_from_xy(result[prefix+'_x'], result[prefix+'_y'])
    result['w_surface'] = source.w
    result['Ro'] = source.w/result.f
    with np.errstate(divide='ignore', invalid='ignore'):
        result['topo_plan_ratio'] = np.log(result.PV_grad_topo_mag/result.PV_grad_plan_mag)
    peaks = kernels._select_profile_peak_w(source, vertical, settings['peak_vorticity_max_depth_m'])
    result['w_peak_signed_s1'] = peaks.w
    result['w_peak_abs_s1'] = peaks.w.abs()
    result['w_peak_depth_m'] = peaks.w_selected_depth_m
    result['w_peak_n'] = peaks.w_profile_n
    return pd.concat([source[KEYS], result], axis=1)


def signature(config):
    # Plot/report edits must not invalidate unchanged numerical products.
    functions = [load_grid, check_keys, calculate, pv_plots.environmental_grid_gradients]
    return dict(schema_version=2, calculation=config['calculation'],
                inputs={k: identity(p) for k,p in config['inputs'].items()},
                kernels=identity(kernels.__file__), keys=KEYS, vector_prefixes=PREFIXES,
                calculation_functions={fn.__name__: hashlib.sha256(
                    inspect.getsource(fn).encode()).hexdigest() for fn in functions})


def save_product(data, signature_, root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    # Publish manifest last, so incomplete writes cannot be mistaken for a valid run.
    with tempfile.NamedTemporaryFile(dir=root, suffix='.parquet', delete=False) as f:
        temporary = Path(f.name)
    try:
        data.to_parquet(temporary, index=False)
        temporary.replace(root/'surface_pv.parquet')
    finally:
        temporary.unlink(missing_ok=True)
    manifest = dict(signature=signature_, rows=len(data), product=identity(root/'surface_pv.parquet'),
                    created_utc=datetime.now(timezone.utc).isoformat())
    temporary_manifest = root/'manifest.json.tmp'
    temporary_manifest.write_text(json.dumps(manifest, indent=2)+'\n')
    temporary_manifest.replace(root/'manifest.json')


def load_product(config):
    root = Path(config['product_dir'])
    manifest = json.loads((root/'manifest.json').read_text())
    if manifest['signature'] != signature(config):
        raise ValueError('PV product inputs/settings/code changed; rebuild explicitly')
    if manifest['product'] != identity(root/'surface_pv.parquet'):
        raise ValueError('PV product file changed; rebuild explicitly')
    d = pd.read_parquet(root/'surface_pv.parquet')
    check_keys(d)
    if len(d) != manifest['rows']:
        raise ValueError('PV product row count mismatch')
    return d


def validate(data):
    check_keys(data)
    valid = data.PV_grad_mag.notna()
    for axis in ['x','y']:
        if not np.allclose(data.loc[valid,'PV_grad_'+axis],
                           data.loc[valid,'PV_grad_plan_'+axis]+data.loc[valid,'PV_grad_topo_'+axis],
                           rtol=1e-12,atol=1e-28):
            raise ValueError('PV component closure failed')
    for prefix in PREFIXES:
        ok = np.isfinite(data[prefix+'_mag'])
        if not np.allclose(data.loc[ok,prefix+'_mag'], np.hypot(
                data.loc[ok,prefix+'_x'],data.loc[ok,prefix+'_y']), rtol=1e-12, atol=0):
            raise ValueError('Vector magnitude failed')
        coherence = data[prefix+'_coherence'].dropna()
        if not coherence.between(-1e-12,1+1e-12).all():
            raise ValueError('Coherence outside [0,1]')


def run(config_path=None, build=False, rebuild=False):
    path = Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    c = load_config(path)
    grid = load_grid(c['inputs']['grid'])
    surface = pd.read_parquet(c['inputs']['surface'])
    check_keys(surface)
    if build:
        root = Path(c['product_dir'])
        if any((root/name).exists() for name in ['surface_pv.parquet','manifest.json']) and not rebuild:
            raise FileExistsError('Product exists; load normally or use --rebuild explicitly')
        vertical = pd.read_parquet(c['inputs']['vertical'])
        product = calculate(surface, vertical, grid, c['calculation'])
        validate(product)
        save_product(product, signature(c), root)
    else:
        product = load_product(c)
        validate(product)
        vertical = pd.read_parquet(c['inputs']['vertical'])
    if not surface[KEYS].equals(product[KEYS]):
        # Preserve exact source identity/order instead of silently accepting partial caches.
        raise ValueError('PV keys/order differ from source days')
    requested = pd.DataFrame(c['figure']['examples'])
    check_keys(requested)
    selected = requested.merge(surface, on=KEYS, how='left', validate='one_to_one', indicator=True)
    if not selected._merge.eq('both').all():
        raise ValueError('Requested example missing from source')
    selected = selected.drop(columns='_merge').merge(product, on=KEYS, validate='one_to_one')
    f = c['figure']
    # Keep the tagged Figure 2 helper's defaults and explicit final-cell presentation.
    fig, axes, references = pv_plots.plot_gradient_examples(
        selected, grid, vertical=vertical, figsize=f['size_inches'],
        vert_lim=f['profile_depth_limit_m'], min_pad_km=f['min_pad_km'],pad_rc=f['pad_rc'],
        max_local_arrows=f['max_local_arrows'],reference_percentile=f['reference_percentile'],
        local_arrow_length_inches=f['local_arrow_length_inches'],net_arrow_fraction=f['net_arrow_fraction'],
        xq=f['quiver_key_position'][0],yq=f['quiver_key_position'][1],title=False,quiv_clr='magenta')
    for ax, label in zip(axes, ['a)', 'b)']):
        ax.text(-.06,1.05,label,transform=ax.transAxes,ha='left',va='top',fontsize=12,fontweight='bold')
    out = Path(c['output_dir']);out.mkdir(parents=True,exist_ok=True)
    figure_dir=Path(c['figure_dir']);figure_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure_dir/'pv_measurement_demo.png',dpi=f['dpi'],bbox_inches='tight')
    cols = KEYS+['Cyc','PV_grad_mag','PV_grad_theta','PV_grad_mean_local_mag','PV_grad_coherence',
                'PV_footprint_n','w_surface','w_peak_signed_s1','w_peak_depth_m']
    examples = selected[cols].copy();examples['local_arrow_reference'] = references
    examples.to_csv(out/'examples.csv',index=False)
    status = product[KEYS].copy()
    status['status'] = np.where(np.isfinite(product.PV_grad_mag),'valid','no finite footprint gradient')
    status.to_csv(out/'sample_audit.csv',index=False)
    summary = surface[KEYS+['Cyc']].merge(product,on=KEYS,validate='one_to_one').groupby('Cyc').agg(
        eddy_days=('Day','size'),eddies=('Eddy','nunique'),finite_pv_days=('PV_grad_mag','count'),
        median_pv=('PV_grad_mag','median'),median_coherence=('PV_grad_coherence','median'))
    summary.to_csv(out/'summary.csv')
    lines = ['# Surface PV results','',
             'Surface ESP geometry and surface relative vorticity; Gaussian-envelope weights inside '
             'the maximum-velocity contour (FRAC=1). Local planetary and topographic vectors are '
             'evaluated before spatial averaging. No temporal smoothing or water-depth selection.','',
             'Primary PV_grad = grad(f)/h − (f+zeta) grad(h)/h². The internal grad(zeta)/h term '
             'and full gradient are separate diagnostics. Units: gradient m⁻² s⁻¹; w s⁻¹; depth m.','',
             f'{len(product):,} source eddy-days; {np.isfinite(product.PV_grad_mag).sum():,} finite '
             'environmental gradients. Invalid footprints remain in the audit.','']
    for row in examples.itertuples():
        lines.append(f'- {row.Cyc}, Eddy {row.Eddy}, Day {row.Day}: gradient {row.PV_grad_mag:.6e} '
                     f'm⁻² s⁻¹, bearing {row.PV_grad_theta:.2f}°, coherence {row.PV_grad_coherence:.3f}; '
                     f'{row.PV_footprint_n:g} contributing cells.')
    lines += ['', 'Figure 2 preserves the tagged original plotting style and fixed examples. Local '
              'magenta arrows have panel-specific scales; the green mean arrow has a fixed display '
              'length and does not encode magnitude. Both arrow types are rotated from geographic '
              'east/north into model-grid x/y using the grid angle (radians), correcting the original '
              'display. Stored components and bearings remain geographic; magnitudes are unchanged. '
              'The fitted-centre profile is illustrative.', '',
              'Peak absolute vorticity over 0–1000 m is saved with signed value, depth and support '
              'for future use; it does not enter PV. Primary values differ intentionally from the '
              'old shared peak-vorticity PV cache. This calculation is an environmental shallow-water '
              'proxy, not full three-dimensional Ertel PV.','',
              'Product: '+c['product_dir']+'/surface_pv.parquet',
              'Writing tables: examples.csv, summary.csv; exclusions: sample_audit.csv. '
              'No inferential test or uncertainty interval is implied by this method illustration.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try: revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError): revision=None
    manifest = dict(config=c, config_file=identity(path), code=identity(__file__),
                    plotting_code=identity(pv_plots.__file__), kernels=identity(kernels.__file__),
                    product_manifest=json.loads((Path(c['product_dir'])/'manifest.json').read_text()),
                    git_revision=revision, python=sys.version, numpy=np.__version__,pandas=pd.__version__,
                    platform=platform.platform(),created_utc=datetime.now(timezone.utc).isoformat())
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(output_dir=out,figure=fig,examples=examples,summary=summary)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    parser.add_argument('--build',action='store_true',help='Build new scratch product')
    parser.add_argument('--rebuild',action='store_true',help='Explicitly replace this new product')
    args=parser.parse_args()
    if args.rebuild and not args.build: parser.error('--rebuild requires --build')
    result=run(args.config,build=args.build,rebuild=args.rebuild)
    print(result['output_dir'])
