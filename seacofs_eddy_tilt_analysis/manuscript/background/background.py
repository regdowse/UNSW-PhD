"""Figure 11 and a compact matched-sample 31/61/91-day seasonal-window check."""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
location=str(ROOT/'pv_dominance')
if location not in sys.path:sys.path.insert(0,location)
import dominance
import background_kernels as kernels
import background_plots

pv=dominance.pv
KEYS=['Eddy','Day']
LEVELS={'surface':'Surface','200':'0–200 m','200_500':'200–500 m'}


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text())
    for name in ['dominance_config','model_root','sigma_depths','product_dir','output_dir','figure_dir']:
        p=Path(c[name]).expanduser();c[name]=str(p if p.is_absolute() else path.parent/p)
    a=c['calculation']
    if (a['upper_depth_m'],a['deep_depth_m'])!=(200,500):
        raise ValueError('This schema uses fixed 0–200 and 200–500 m layers')
    if a['primary_window_days']!=91 or sorted(a['windows_days'])!=[31,61,91]:
        raise ValueError('Agreed windows are 31/61/91, primary 91 days')
    if a['workers']<1 or c['bootstrap_replicates']<100:
        raise ValueError('Invalid worker/bootstrap count')
    return c


def contexts(c):
    dc=dominance.load_config(c['dominance_config'])
    pc=pv.load_config(dc['pv_config']);cc=dominance.census.load_config(dc['census_config'])
    if pc['inputs']['surface']!=cc['inputs']['surface'] or pc['inputs']['grid']!=cc['inputs']['grid']:
        raise ValueError('Source configurations disagree')
    return dc,pc,cc


def archive_files(root):
    paths=sorted(Path(root).glob('outer_avg_*.nc'))
    if not paths:raise FileNotFoundError('No model archive files')
    # The archive is ~0.8 TB. Record exact file paths, sizes, mtimes and ocean-time
    # values, rather than claiming a content hash of every raw model file.
    return [dict(path=str(p.resolve()),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths]


def inventory(records,grid):
    rows=[];all_days=[]
    for record in records:
        with nc.Dataset(record['path']) as ds:
            if ds['ocean_time'].units!='seconds since 1990-01-01 00:00:00':
                raise ValueError('Unexpected model time origin/units')
            times=kernels.clean(ds['ocean_time'][:])/86400
            if not np.isfinite(times).all() or not np.allclose(times,np.rint(times),rtol=0,atol=1e-8):
                raise ValueError('Model output is not on integer dataset days')
            days=np.rint(times).astype(int)
            if not len(days):raise ValueError('Empty model file')
            if ds['u_eastward'].dimensions!=('ocean_time','s_rho','eta_rho','xi_rho') or ds['v_northward'].shape!=ds['u_eastward'].shape:
                raise ValueError('Unexpected geographic velocity dimensions')
            if ds['u_eastward'].shape[2:]!=grid.h.T.shape:
                raise ValueError('Model grid shape changed')
            if not np.allclose(ds['h'][:].T,grid.h,rtol=0,atol=1e-7):
                raise ValueError('Model bathymetry changed across archive')
            all_days.extend(days.tolist());rows.append(dict(**record,days=days.tolist()))
    ordered=np.sort(all_days)
    if len(np.unique(ordered))!=len(ordered) or not np.all(np.diff(ordered)==1):
        raise ValueError('Duplicate model times or missing archive days; discuss before building')
    return rows


def signature(c):
    _,pc,_=contexts(c)
    return dict(schema_version=1,upper_depth_m=200,deep_depth_m=500,
        layer_method='static physical centres; midpoint interfaces; complete-layer finite coverage',
        source_surface=pv.identity(pc['inputs']['surface']),grid=pv.identity(pc['inputs']['grid']),
        sigma_depths=pv.identity(c['sigma_depths']),kernels=pv.identity(kernels.__file__),
        archive=archive_files(c['model_root']))


def write_json(path,value):
    path=Path(path);temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def build(c):
    from joblib import Parallel,delayed
    root=Path(c['product_dir']);manifest_path=root/'manifest.json'
    if manifest_path.exists():
        load_product(c);print('Validated existing independent product',flush=True);return
    sig=signature(c);state=root/'building.json'
    if state.exists():
        if json.loads(state.read_text())!=sig:
            raise ValueError('Incomplete build has different inputs/code; choose a new product_dir')
    else:
        if root.exists() and any(root.iterdir()):
            raise ValueError('Product directory is not empty; choose a fresh product_dir')
        root.mkdir(parents=True,exist_ok=True);write_json(state,sig)
    _,pc,_=contexts(c);grid=pv.load_grid(pc['inputs']['grid'])
    source=pd.read_parquet(pc['inputs']['surface']);pv.check_keys(source)
    inv_path=root/'archive_inventory.json'
    inv=json.loads(inv_path.read_text()) if inv_path.exists() else inventory(sig['archive'],grid)
    days=np.array([day for row in inv for day in row['days']])
    if not source.Day.isin(days).all():raise ValueError('Source days outside archive')
    write_json(inv_path,inv)
    z=np.load(c['sigma_depths'])
    if z.shape[:2]!=grid.h.shape:raise ValueError('Physical-depth/grid shape mismatch')
    weights={d:kernels.sigma_layer_weights(z,grid.h,d) for d in (200,500)}
    print(f'Archive: {len(inv)} files, {len(days)} continuous days, {days.min()}–{days.max()}; '
          f'reading all {z.shape[-1]} physical sigma levels.',flush=True)
    partitions=Parallel(n_jobs=c['calculation']['workers'],prefer='processes',verbose=10)(
        delayed(kernels.process_file)(row['path'],weights,root/'file_partitions') for row in inv)
    monthly_path=root/'monthly_climatology.npz'
    climatology=kernels.reduce_partitions([p for p,_ in partitions],monthly_path)
    print('Monthly climatology reduced; attaching all source eddy-days.',flush=True)
    products=[pv.identity(monthly_path),pv.identity(inv_path)]
    for window in c['calculation']['windows_days']:
        data=kernels.attach(source,climatology,window);pv.check_keys(data)
        path=root/f'background_{window}.parquet';temporary=path.with_suffix('.tmp')
        data.to_parquet(temporary,index=False);temporary.replace(path);products.append(pv.identity(path))
    if signature(c)!=sig:raise ValueError('Inputs changed during background build')
    write_json(manifest_path,dict(created_utc=datetime.now(timezone.utc).isoformat(),signature=sig,
        products=products,archive_days=len(days),first_day=int(days.min()),last_day=int(days.max()),
        computed_partitions=sum(status=='computed' for _,status in partitions),
        resumed_partitions=sum(status=='resumed' for _,status in partitions)))
    print(f'Fresh product complete: {root}',flush=True)


def load_product(c):
    root=Path(c['product_dir']);path=root/'manifest.json'
    if not path.exists():raise FileNotFoundError('Run background.py --build first')
    manifest=json.loads(path.read_text())
    if manifest['signature']!=signature(c):
        raise ValueError('Background product inputs/code changed; choose a fresh product_dir')
    for item in manifest['products']:
        if pv.identity(item['path'])!=item:raise ValueError('Background product content changed')
    return {w:pd.read_parquet(root/f'background_{w}.parquet') for w in c['calculation']['windows_days']},manifest


def select(surface,backgrounds,grid):
    """One identical complete-layer sample for all levels/windows/polarities."""
    audit=surface[KEYS+['Cyc','regime','target_population','direction_valid','h','ic','jc']].copy()
    target=surface.target_population.eq('deep_planetary')&surface.direction_valid
    i=surface.ic.to_numpy(float);j=surface.jc.to_numpy(float)
    cell_ok=np.isfinite(i)&np.isfinite(j)&(i==np.floor(i))&(j==np.floor(j))&(i>=0)&(i<grid.h.shape[0])&(j>=0)&(j<grid.h.shape[1])
    centre_h=np.full(len(surface),np.nan);wet=np.zeros(len(surface),bool)
    centre_h[cell_ok]=grid.h[i[cell_ok].astype(int),j[cell_ok].astype(int)]
    wet[cell_ok]=grid.mask_rho[i[cell_ok].astype(int),j[cell_ok].astype(int)]==1
    audit['centre_water_depth_m']=centre_h
    eligible=cell_ok&wet&np.isfinite(centre_h)&(centre_h>=500)
    all_complete=np.ones(len(surface),bool);joined={}
    for window,bg in backgrounds.items():
        pv.check_keys(bg)
        if set(map(tuple,bg[KEYS].to_numpy()))!=set(map(tuple,surface[KEYS].to_numpy())):
            raise ValueError('Background/source key coverage differs')
        # Retain source grid indices as authority; reject mismatched product sampling.
        d=surface.merge(bg,on=KEYS,validate='one_to_one',suffixes=('','_background'))
        if not (d.ic.eq(d.ic_background)&d.jc.eq(d.jc_background)).all():
            raise ValueError('Background centre grid indices disagree')
        finite=np.isfinite(d[list(kernels.FIELDS)].to_numpy()).all(axis=1)
        audit[f'complete_{window}']=finite;all_complete &= finite
        joined[window]=d
    audit['selection_status']=np.select([~target,~eligible,~all_complete],
        ['outside deep-planetary directional population','invalid centre cell or <500 m coverage',
         'incomplete background in one or more windows'],default='selected')
    kept=audit.selection_status.eq('selected').to_numpy()
    return audit,{w:d.loc[kept].copy() for w,d in joined.items()}


def summaries(samples,c):
    rows=[];contrasts=[]
    for window,d in samples.items():
        for cyc in ['AE','CE']:
            q=d[d.Cyc.eq(cyc)]
            for level in LEVELS:
                for component,prefix in [('east','u'),('north','v')]:
                    name=f'{prefix}_{level}'
                    seed=c['seed']+len(rows)%12
                    median,lo,hi=kernels.cluster_median(q,name,c['bootstrap_replicates'],seed)
                    rows.append(dict(window_days=window,Cyc=cyc,level=level,component=component,
                        eddy_days=len(q),eddies=q.Eddy.nunique(),median=median,ci_low=lo,ci_high=hi,bootstrap_seed=seed))
            name='surface_minus_deep_east_ms'
            seed=c['seed']+1000+len(contrasts)%2
            median,lo,hi=kernels.cluster_median(q,name,c['bootstrap_replicates'],seed)
            contrasts.append(dict(window_days=window,Cyc=cyc,eddy_days=len(q),eddies=q.Eddy.nunique(),
                median=median,ci_low=lo,ci_high=hi,q25=q[name].quantile(.25),q75=q[name].quantile(.75),
                positive_fraction=q[name].gt(0).mean(),bootstrap_seed=seed))
    return pd.DataFrame(rows),pd.DataFrame(contrasts)


def monthly_coverage(c,audit,manifest):
    """Expose actual daily support behind each monthly/layer mean at target cells."""
    days=np.arange(manifest['first_day'],manifest['last_day']+1)
    months=(kernels.ORIGIN+pd.to_timedelta(days,unit='D')).month
    expected=np.bincount(months,minlength=13)
    target=audit[audit.target_population.eq('deep_planetary')&audit.direction_valid]
    rows=[]
    with np.load(Path(c['product_dir'])/'monthly_climatology.npz') as monthly:
        for cyc in ['AE','CE']:
            q=target[target.Cyc.eq(cyc)]
            i=q.ic.to_numpy(int);j=q.jc.to_numpy(int)
            for month in range(1,13):
                for field in kernels.FIELDS:
                    counts=monthly[f'count_m{month:02d}_{field}'][i,j]
                    rows.append(dict(Cyc=cyc,month=month,field=field,eddy_days=len(q),
                        expected_archive_days=int(expected[month]),minimum_days=int(counts.min()),
                        maximum_days=int(counts.max()),incomplete_eddy_days=int((counts!=expected[month]).sum())))
    return pd.DataFrame(rows)


def run(config_path=None):
    path=Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    c=load_config(path);dc,pc,cc=contexts(c)
    backgrounds,product_manifest=load_product(c)
    product=pv.load_product(pc);pv.validate(product)
    source=pd.read_parquet(pc['inputs']['surface']);tilt=pd.read_parquet(cc['inputs']['tilt'])
    surface=dominance.prepare(source,tilt,product,dc);grid=pv.load_grid(pc['inputs']['grid'])
    audit,samples=select(surface,backgrounds,grid)
    if any(d.empty for d in samples.values()):raise ValueError('No complete background samples; inspect coverage')
    coverage=monthly_coverage(c,audit,product_manifest)
    stats,contrasts=summaries(samples,c)
    primary=c['calculation']['primary_window_days']
    with plt.rc_context(matplotlib.rcParamsDefault), plt.style.context(c['plot']['style']):
        fig,colour_limit=background_plots.plot(samples[primary],stats[stats.window_days.eq(primary)],
            grid,dc['water_depth_threshold_m'],c['plot'])
    out=Path(c['output_dir']);out.mkdir(parents=True,exist_ok=True)
    figure_dir=Path(c['figure_dir']);figure_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure_dir/f'plan_back_flow.png',dpi=c['plot']['dpi'],bbox_inches='tight')
    stats.to_csv(out/'flow_summary.csv',index=False);contrasts.to_csv(out/'window_comparison.csv',index=False)
    coverage.to_csv(out/'monthly_coverage.csv',index=False)
    audit.to_csv(out/'sample_audit.csv',index=False)
    exclusions=audit.groupby(['Cyc','selection_status'],observed=True,as_index=False).agg(eddy_days=('Day','size'),eddies=('Eddy','nunique'))
    exclusions.to_csv(out/'exclusions.csv',index=False)
    members=pd.concat([d.assign(window_days=w) for w,d in samples.items()],ignore_index=True)
    cols=KEYS+['Cyc','window_days','xc','yc','h','ic','jc','TiltDis','TiltDir']+list(kernels.FIELDS)+['u_200_500','v_200_500','surface_minus_deep_east_ms']
    members[cols].to_parquet(out/'members.parquet',index=False)
    # Paired differences are checks of window sensitivity, not extra paper panels.
    differences=[];reference=samples[primary].set_index(KEYS)
    for w,d in samples.items():
        aligned=d.set_index(KEYS).loc[reference.index]
        delta=aligned.surface_minus_deep_east_ms-reference.surface_minus_deep_east_ms
        for cyc in ['AE','CE']:
            diff=delta[reference.Cyc.eq(cyc)]
            differences.append(dict(window_days=w,Cyc=cyc,eddy_days=len(diff),median_paired_difference_ms=diff.median(),
                q25=diff.quantile(.25),q75=diff.quantile(.75),max_abs_difference_ms=diff.abs().max()))
    paired=pd.DataFrame(differences);paired.to_csv(out/'window_paired_differences.csv',index=False)
    first=kernels.ORIGIN+pd.Timedelta(days=product_manifest['first_day'])
    last=kernels.ORIGIN+pd.Timedelta(days=product_manifest['last_day'])
    lines=['# Planetary eddies: background flow','',
        f'Figure 11. Fresh archive build: {product_manifest["archive_days"]:,} continuous daily fields '
        f'({first.date()}–{last.date()}); monthly climatological means use all archive days, including '
        'days without eddies. Primary: centred 91-day seasonal climatology, weighted by calendar-month '
        'overlap. This is not an event-year running mean; 31/61-day checks use the same monthly fields.','',
        f'Population: unsmoothed factor-{dc["dominance_factor"]:g} planetary dominance; footprint-mean '
        f'h≥{dc["water_depth_threshold_m"]:g} m; saved tilt≥{dc["direction_min_tilt_km"]:g} km and shared '
        'directional QC. Centre-column water depth is checked separately for full 500 m coverage. '
        'All three windows and all displayed layers use identical eddy-days.','',
        'Velocities are geographic east/north at the source eddy centre grid cell. Layer means use '
        'the original static physical sigma-centre depths and midpoint-interface thickness '
        'approximation, bounded by surface and bathymetry. Surface denotes the top model level. '
        'Require finite full-layer coverage; do not renormalise partially missing depths. '
        '200–500 m mean = (500 U(0–500) − 200 U(0–200))/300. The primary contrast is '
        'surface minus 200–500 m zonal velocity, in m/s, not a vertical derivative.','',
        f'Pooled eddy-day medians; pointwise 95% CIs resample whole eddy tracks '
        f'({c["bootstrap_replicates"]} draws, seed {c["seed"]} with recorded statistic offsets). '
        'The median of a daily contrast differs from the contrast between layer medians.','',
        'Surface-minus-deep zonal contrast (median [95% CI], m/s):','']
    for row in contrasts.itertuples():
        lines.append(f'- {row.window_days}-day {row.Cyc}: {row.median:.4f} [{row.ci_low:.4f}, {row.ci_high:.4f}]; '
                     f'{row.eddy_days:,} days / {row.eddies:,} eddies; positive fraction {row.positive_fraction:.3f}.')
    lines+=['','Internal window checks: window_comparison.csv and window_paired_differences.csv. '
        'Do not select a window for stronger results. These are background-flow descriptions, '
        'not independent observations of shallow/deep eddy propagation or evidence of causation.','',
        f'Original three-panel style retained. Map colours saturate at the shared '
        f'{c["plot"]["colour_percentile"]:g}th absolute percentile (±{colour_limit:.4f} m/s); '
        'all mapped days remain included. The 3000 m contour shows local bathymetry, while '
        'population selection uses footprint-mean h.','',
        'Tables: flow_summary.csv (layer medians/CIs), window_comparison.csv (daily contrast and '
        'support), window_paired_differences.csv, exclusions.csv and sample_audit.csv. '
        'members.parquet retains all matched window values; manifest.json records config/code '
        'hashes and the independent scratch-product provenance. Outputs are descriptive.']
    lines+=['', 'Monthly coverage: monthly_coverage.csv reports original daily support for '
        'every layer/component/month at the targeted cells before window exclusions. '
        + ('All target cells have every archive day for every required layer/component.'
           if coverage.incomplete_eddy_days.eq(0).all() else
           'Some target cells have incomplete original daily support; review this before interpreting layer contrasts.')]
    lines+=['','Selection/exclusions:','']
    for row in exclusions.itertuples():
        lines.append(f'- {row.Cyc}, {row.selection_status}: {row.eddy_days:,} days / {row.eddies:,} eddies.')
    lines+=['','Window sensitivity (population medians; paired daily differences remain in the internal table):','']
    for cyc in ['AE','CE']:
        q=contrasts[contrasts.Cyc.eq(cyc)]
        lines.append(f'- {cyc}: median contrast range {q["median"].min():.4f}–{q["median"].max():.4f} m/s; '
                     f'all three lower 95% CI bounds above zero: {bool(q.ci_low.gt(0).all())}.')
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    write_json(out/'manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),config=c,
        config_file=pv.identity(path),dominance_config=dc,dominance_config_file=pv.identity(c['dominance_config']),
        pv_config=pv.identity(dc['pv_config']),census_config=pv.identity(dc['census_config']),
        tilt=pv.identity(cc['inputs']['tilt']),background_product=product_manifest,
        pv_product=json.loads((Path(pc['product_dir'])/'manifest.json').read_text()),
        code=[pv.identity(f) for f in [__file__,kernels.__file__,background_plots.__file__,dominance.__file__,pv.__file__,dominance.census.__file__,dominance.census_plots.__file__]],
        git_revision=revision,python=sys.version,numpy=np.__version__,pandas=pd.__version__,
        netCDF4=nc.__version__,matplotlib=matplotlib.__version__))
    return dict(output_dir=out,figure=fig,statistics=stats,window_comparison=contrasts)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    parser.add_argument('--build',action='store_true',help='Build/resume the fresh scratch product explicitly')
    args=parser.parse_args()
    if args.build:build(load_config(args.config))
    else:print(run(args.config)['output_dir'])
