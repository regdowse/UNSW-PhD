"""Figure 17: actual constituent centres in each day's environmental-PV frame."""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
location=str(ROOT/'pv_dominance')
if location not in sys.path:sys.path.insert(0,location)
import dominance
import composite_tools as tools
import composite_plots

pv=dominance.pv
KEYS=['Eddy','Day']
REGIMES=['planetary','mixed','topographic']


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text())
    for name in ['dominance_config','output_dir','figure_dir']:
        p=Path(c[name]).expanduser();c[name]=str(p if p.is_absolute() else path.parent/p)
    if not np.isfinite(c['max_profile_depth_m']) or c['max_profile_depth_m']<=0:
        raise ValueError('Need positive finite depth cap')
    if c['bootstrap_replicates']<100:raise ValueError('Need at least 100 bootstrap draws')
    if c['plot']['min_eddies']<2:raise ValueError('Need at least two eddies for displayed CIs')
    return c


def rotate_centres(x,y,bearing_deg,grid_angle):
    """Model x/y to along total PV and left (+CCW), using a proper rotation."""
    angle=np.pi/2-np.deg2rad(bearing_deg)-grid_angle
    c,s=np.cos(angle),np.sin(angle)
    return c*x+s*y,-s*x+c*y


def prepare(surface,vertical,grid_angle,depth_cap):
    """Take agreed day selections, then retain actual valid capped profile levels."""
    pv.check_keys(surface)
    if not np.isscalar(grid_angle) or not np.isfinite(grid_angle):
        raise ValueError('Need a finite scalar grid angle in radians')
    audit=surface.copy()
    target=(audit.regime.eq('mixed')|audit.target_population.isin(['deep_planetary','shallow_topographic']))
    audit['selection_status']=np.select([
        audit.regime.eq('unclassified'),~target,~audit.direction_valid],
        ['unclassified PV magnitudes','outside targeted depth range or invalid water depth',
         'fails tilt/PV directional QC'],default='selected')
    sample=audit[audit.selection_status.eq('selected')]
    capped=vertical[vertical.Depth.between(0,depth_cap)].copy()
    pv.check_keys(capped,KEYS+['Depth'])
    valid=tools.valid_profiles(capped).merge(sample[KEYS+['Cyc','regime','PV_grad_theta']],
                                            on=KEYS,validate='many_to_one')
    refs=valid.loc[valid.Depth.eq(0),KEYS+['xc','yc']].rename(columns={'xc':'surface_xc','yc':'surface_yc'})
    valid=valid.merge(refs,on=KEYS,validate='many_to_one')
    valid['dx_km']=valid.xc-valid.surface_xc;valid['dy_km']=valid.yc-valid.surface_yc
    valid['along_km'],valid['perp_km']=rotate_centres(valid.dx_km,valid.dy_km,valid.PV_grad_theta,grid_angle)
    valid['distance_km']=np.hypot(valid.along_km,valid.perp_km)
    support=valid.groupby(KEYS).Depth.agg(['min','max','size']).reset_index().rename(
        columns={'min':'shallowest_fit_m','max':'deepest_fit_m','size':'valid_levels'})
    audit=audit.merge(support,on=KEYS,how='left',validate='one_to_one')
    audit.loc[audit.selection_status.eq('selected') & audit.valid_levels.isna(),
              'selection_status']='no valid surface-referenced profile within depth cap'
    members=valid[KEYS+['Cyc','regime','Depth','dx_km','dy_km','along_km','perp_km','distance_km']].copy()
    return audit,members


def summarise(members,n_boot,seed):
    parts=[];inventory=[]
    for regime in REGIMES:
        for cyc in ['AE','CE']:
            d=members[members.regime.eq(regime)&members.Cyc.eq(cyc)]
            inventory.append(dict(regime=regime,Cyc=cyc,eddy_days=len(d[KEYS].drop_duplicates()),
                                  eddies=d.Eddy.nunique(),deepest_contributing_m=d.Depth.max()))
            if d.empty:continue
            adapter=d[KEYS+['Depth','along_km','perp_km']].rename(columns={'along_km':'xc','perp_km':'yc'})
            s=tools.centre_statistics(adapter,np.sort(d.Depth.unique()),n_boot=n_boot,seed=seed)
            s=s.drop(columns=['mean_tilt_x','mean_tilt_y'])
            rename={c:c.replace('xc','along').replace('yc','perp') for c in s if 'xc' in c or 'yc' in c}
            s=s.rename(columns=rename).rename(columns={'cov_xy':'cov_along_perp'})
            parts.append(s.assign(regime=regime,Cyc=cyc))
    if not parts:raise ValueError('No profiles passed selection; inspect day/profile QC')
    return pd.concat(parts,ignore_index=True),pd.DataFrame(inventory)


def run(config_path=None):
    path=Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    c=load_config(path);dc=dominance.load_config(c['dominance_config']);pvc=pv.load_config(dc['pv_config'])
    cc=dominance.census.load_config(dc['census_config'])
    if pvc['inputs']['surface']!=cc['inputs']['surface'] or pvc['inputs']['grid']!=cc['inputs']['grid']:
        raise ValueError('Source configurations disagree')
    source=pd.read_parquet(pvc['inputs']['surface']);tilt=pd.read_parquet(cc['inputs']['tilt'])
    product=pv.load_product(pvc);pv.validate(product)
    surface=dominance.prepare(source,tilt,product,dc)
    vertical=pd.read_parquet(pvc['inputs']['vertical'],columns=KEYS+['Depth','xc','yc','Rc','Omega','q11','q12','q22'])
    grid=pv.load_grid(pvc['inputs']['grid'])
    audit,members=prepare(surface,vertical,grid.angle,c['max_profile_depth_m'])
    stats,inventory=summarise(members,c['bootstrap_replicates'],c['seed'])
    with plt.rc_context(matplotlib.rcParamsDefault):
        fig,_=composite_plots.plot(stats,max_depth_m=c['max_profile_depth_m'],
            min_eddies=c['plot']['min_eddies'],sparse_eddies=c['plot']['sparse_eddies'],
            figsize=c['plot']['size_inches'])
    out=Path(c['output_dir']);out.mkdir(parents=True,exist_ok=True)
    figure_dir=Path(c['figure_dir']);figure_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(figure_dir/f'composites.png',dpi=c['plot']['dpi'],bbox_inches='tight')
    stats.to_csv(out/'centre_statistics.csv',index=False);inventory.to_csv(out/'inventory.csv',index=False)
    members.to_parquet(out/'members.parquet',index=False)
    audit_cols=KEYS+['Cyc','regime','target_population','h','direction_valid','depth_selection_status',
                    'selection_status','shallowest_fit_m','deepest_fit_m','valid_levels']
    audit[audit_cols].to_csv(out/'day_audit.csv',index=False)
    exclusions=audit.groupby(['Cyc','regime','selection_status'],as_index=False).agg(eddy_days=('Day','size'),eddies=('Eddy','nunique'))
    exclusions.to_csv(out/'exclusions.csv',index=False)
    cap=c['max_profile_depth_m'];threshold=dc['water_depth_threshold_m']
    lines=['# PV-relative centreline composites','',
        'Figure 17: original six-panel layout and AE-red/CE-blue curves. Columns are deep '
        f'planetary (h≥{threshold:g} m), magnitude-defined mixed (all depths), and shallow '
        f'topographic (h≤{threshold:g} m). Dominance uses unsmoothed component magnitudes '
        f'with factor {dc["dominance_factor"]:g}; excluded other-depth dominated days are not mixed. '
        'h is Gaussian-footprint-mean water depth, separate from fitted profile depth.','',
        f'Day selection uses saved tilt ≥{dc["direction_min_tilt_km"]:g} km and the shared valid '
        f'tilt/PV direction requirements. Use actual fitted levels within 0–{cap:g} m, including '
        'surface, with a valid surface fit. No interpolation, extrapolation or velocity '
        'reconstruction. Existing profile QC requires finite centres/ESP parameters, positive '
        'radius and positive-definite ellipse matrices.','',
        'Each day is referenced to its own fitted surface centre and rotated into the local '
        'total environmental-PV frame before averaging: positive along is with the signed '
        'gradient; positive perpendicular is 90° counterclockwise. These are surface-to-depth '
        'constituent-centre displacements, not the saved deep-to-shallow fitted tilt vector.','',
        'Each available eddy-day has equal weight at each depth. Pointwise 95% confidence '
        f'intervals resample whole eddy tracks ({c["bootstrap_replicates"]} draws, seed {c["seed"]}); '
        'the same track draw is reused across depths within a regime/polarity. Bands quantify '
        'uncertainty in the mean, not population spread or fit-parameter uncertainty. '
        'Support can change with depth. These are descriptive composites, not a causal test.','',
        'Included profiles:','']
    for row in inventory.itertuples():
        lines.append(f'- {row.Cyc} {row.regime}: {row.eddy_days:,} days / {row.eddies:,} eddies; '
                     f'deepest contributing fit {row.deepest_contributing_m:.2f} m.')
    lines+=['','Deepest supported displayed points (mean along, perpendicular; pointwise 95% CI, km):','']
    for (regime,cyc),d in stats.groupby(['regime','Cyc']):
        supported=d[d.n_eddies.ge(c['plot']['min_eddies'])].sort_values('Depth')
        if supported.empty:continue
        row=supported.iloc[-1]
        lines.append(f'- {cyc} {regime}, {row.Depth:.2f} m: along {row.mean_along:.2f} '
            f'[{row.mean_along_ci_low:.2f}, {row.mean_along_ci_high:.2f}], perpendicular '
            f'{row.mean_perp:.2f} [{row.mean_perp_ci_low:.2f}, {row.mean_perp_ci_high:.2f}]; '
            f'{int(row.n_eddy_days):,} days / {int(row.n_eddies):,} eddies.')
    lines+=['',f'Curves require ≥{c["plot"]["min_eddies"]} eddies; open symbols flag fewer than '
        f'{c["plot"]["sparse_eddies"]}. Track counts across regimes overlap and must not be '
        'summed as distinct eddies. Missing levels contribute no artificial zero or extra depth '
        'coverage. The plot axis extends to the configured cap; curves stop at actual data.', '',
        'Writing tables: centre_statistics.csv (means, CI, covariance/spread and depth support), '
        'inventory.csv, exclusions.csv; day_audit.csv records both dominance and depth selection. '
        'members.parquet preserves projected constituent centres. manifest.json records source, '
        'config/product/code identities and grid rotation. No original source or PV cache was rebuilt.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    dependencies=[__file__,tools.__file__,composite_plots.__file__,dominance.__file__,pv.__file__,pv.kernels.__file__,pv.pv_plots.__file__,dominance.census.__file__]
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),config=c,
        config_file=pv.identity(path),dominance_config=dc,dominance_config_file=pv.identity(c['dominance_config']),
        pv_config=pv.identity(dc['pv_config']),census_config=pv.identity(dc['census_config']),
        tilt=pv.identity(cc['inputs']['tilt']),grid_angle_rad=grid.angle,
        product_manifest=json.loads((Path(pvc['product_dir'])/'manifest.json').read_text()),
        code=[pv.identity(f) for f in dependencies],git_revision=revision,python=sys.version,
        numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(output_dir=out,figure=fig,inventory=inventory,statistics=stats)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();result=run(args.config);print(result['output_dir'])
