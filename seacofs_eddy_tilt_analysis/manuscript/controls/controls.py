"""Figures 13–16 and concise manuscript handoff; fixed fitted depths, no legacy imports."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
for folder in ['pv','census']:
    location=str(ROOT/folder)
    if location not in sys.path:sys.path.insert(0,location)
import pv
import census
import control_tools as ct
import shape_tools as st
import alignment_tools as at

KEYS=['Eddy','Day']


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text())
    for name in ['pv_config','census_config','output_dir','figure_dir']:
        p=Path(c[name]).expanduser();c[name]=str(p if p.is_absolute() else path.parent/p)
    for name in ['structure_depths_m','ellipse_depths_m']:
        z=np.asarray(c[name],float)
        if len(z)<2 or len(set(z))!=len(z) or not np.isfinite(z).all() or z.min()<0 or 0 not in z:
            raise ValueError('Depths must be unique, finite, nonnegative and include zero')
    if c['alignment']['bootstrap_replicates']<2:raise ValueError('Need at least two bootstrap resamples')
    return c


def require_levels(vertical,depths):
    """Fail if exact configured levels disappear; never rank/interpolate/substitute."""
    if (vertical.Depth.dropna()<0).any():raise ValueError('Expected positive-down fitted depths')
    missing=set(depths)-set(vertical.Depth.dropna().unique())
    if missing:raise ValueError(f'Configured fitted depths absent: {sorted(missing)}')


def environment_data(surface,product):
    ct.check_keys(surface);ct.check_keys(product)
    source=ct.identity(surface).merge(surface[KEYS+['lat']],on=KEYS,validate='one_to_one')
    d=source.merge(product[KEYS+['beta','PV_grad_mag']],on=KEYS,how='left',validate='one_to_one')
    d['log10_PV_grad_mag']=np.log10(d.PV_grad_mag.where(d.PV_grad_mag>0))
    d['keep']=np.isfinite(d[['lat','beta','log10_PV_grad_mag']]).all(axis=1)
    return d[d.keep].copy(),d[KEYS+['Cyc','keep']].copy()


def exact_shapes(surface,vertical,depths,bearing_offset):
    """Surface processed ellipse plus exact subsurface fits at configured depths."""
    require_levels(vertical,depths)
    pv.check_keys(vertical,KEYS+['Depth'])
    sampled=vertical.loc[vertical.Depth.isin([z for z in depths if z>0]),KEYS+['Depth']+st.QCOLS].copy()
    sampled['ShapeDepth']=sampled.Depth
    sampled['DepthLower']=sampled.Depth;sampled['DepthUpper']=sampled.Depth
    sampled['DepthMethod']='exact'
    base=surface[KEYS+['Cyc','TiltDis','TiltDir']+st.QCOLS]
    return st.build_analysis_table(base,sampled.drop(columns='Depth'),bearing_offset)


def fit_table(bins):
    """Fits describe displayed bin medians, not independent-day regressions."""
    rows=[]
    for (outcome,metric,cyc),d in bins.groupby(['outcome','metric','Cyc']):
        q=d.dropna(subset=['x','median'])
        if len(q)<2 or q.x.nunique()<2:continue
        coefficients=np.polyfit(q.x,q['median'],1);pred=np.polyval(coefficients,q.x)
        ss=np.sum((q['median']-q['median'].mean())**2)
        rows.append(dict(outcome=outcome,metric=metric,Cyc=cyc,bins=len(q),slope=coefficients[0],
                         intercept=coefficients[1],r2=1-np.sum((q['median']-pred)**2)/ss if ss>0 else np.nan))
    return pd.DataFrame(rows)


def run(config_path=None):
    path=Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    c=load_config(path);pvc=pv.load_config(c['pv_config']);cc=census.load_config(c['census_config'])
    if pvc['inputs']['surface']!=cc['inputs']['surface'] or pvc['inputs']['grid']!=cc['inputs']['grid']:
        raise ValueError('Source configurations disagree')
    source=pd.read_parquet(pvc['inputs']['surface']);tilt=pd.read_parquet(cc['inputs']['tilt'])
    pv.check_keys(source);pv.check_keys(tilt)
    if not pd.MultiIndex.from_frame(tilt[KEYS]).isin(pd.MultiIndex.from_frame(source[KEYS])).all():
        raise ValueError('Tilt keys absent from surface')
    surface=source.drop(columns=['TiltDis','TiltDir'],errors='ignore').merge(
        tilt[KEYS+['TiltDis','TiltDir']],on=KEYS,how='left',validate='one_to_one')
    product=pv.load_product(pvc);pv.validate(product)
    if not product[KEYS].equals(source[KEYS]):raise ValueError('PV key coverage/order changed')
    vertical=pd.read_parquet(pvc['inputs']['vertical'],columns=KEYS+['Depth','Rc','Omega']+st.QCOLS)
    require_levels(vertical,c['structure_depths_m'])
    qc=c['fit_qc']
    clean,qc_counts=ct.prepare_depths(vertical,surface,max_depth=qc['max_depth_m'],rc_max=qc['max_radius_km'],
                                    omega_max=qc['max_abs_omega_s1'],ar_max=qc['max_axis_ratio'])
    matched=ct.matched_depths(clean,c['structure_depths_m'])
    selected={metric:c['structure_depths_m'] for metric in ct.METRICS}
    matched_norm=ct.with_normalised_tilt(matched,surface)
    environment,environment_audit=environment_data(surface,product)
    environment_norm=ct.with_normalised_tilt(environment,surface)
    shapes=exact_shapes(surface,vertical,c['ellipse_depths_m'],c['alignment']['bearing_offset_deg'])
    a=c['alignment']
    alignment=at.direction_data(shapes,depths=c['ellipse_depths_m'],min_ar=a['min_axis_ratio'],
        min_tilt=a['min_tilt_km'],max_ar=a['max_axis_ratio'],min_days=a['min_days_per_eddy'],
        ar_edges=[np.inf if x is None else x for x in a['axis_ratio_edges']],equal_eddy=False,
        n_boot=a['bootstrap_replicates'],seed=a['seed'],min_eddies=a['min_eddies'])
    bins=c['binning'];plot=c['plot']
    with plt.rc_context(matplotlib.rcParamsDefault):
        raw,structure_bins=ct.selected_depth_figure(matched,selected,bins=bins['quantile_bins'],
                                                   min_eddies=bins['min_eddies_per_bin'])
        if plot['figure13_expand_to_iqr']:
            # The original first-panel set_ylim froze the shared upper limit.
            upper=structure_bins[['median','q75']].max().max()
            if np.isfinite(upper):
                for ax in raw.axes:ax.set_ylim(0,max(ax.get_ylim()[1],upper*1.05))
        normalised,structure_norm_bins=ct.selected_depth_figure(matched_norm,selected,bins=bins['quantile_bins'],
            min_eddies=bins['min_eddies_per_bin'],outcome='tilt_over_surface_Rc')
        fig,axs=plt.subplots(2,3,figsize=plot['figure16_size_inches'],layout='constrained')
        env_bins=[]
        for row,(table,outcome,limits) in enumerate([(environment,'TiltDis',plot['figure16_raw_ylim']),
                (environment_norm,'tilt_over_surface_Rc',plot['figure16_normalised_ylim'])]):
            _,_,b=ct.environment_figure2(table,axs=axs[row],bins=bins['quantile_bins'],
                min_eddies=bins['min_eddies_per_bin'],linear_fit=True,outcome=outcome)
            env_bins.append(b.assign(outcome=outcome))
            for ax in axs[row]:ax.set_ylim(limits)
        for ax,label in zip(axs.flat,['a)','b)','c)','d)','e)','f)']):
            ax.text(-.1,1.09,label,transform=ax.transAxes,ha='left',va='top',fontsize=11,fontweight='bold')
        ellipse=at.plot_direction(alignment,plot['ellipse_orientation'])
    figures=dict(eddy_tilt_control=raw,eddy_tilt_control_Rc=normalised,ellipse_tilt2=ellipse,env_tilt_control=fig)
    out=Path(c['output_dir']);out.mkdir(parents=True,exist_ok=True)
    figure_dir=Path(c['figure_dir']);figure_dir.mkdir(parents=True,exist_ok=True)
    for name,figure in figures.items():
        figure.savefig(figure_dir/f'{name}.png',dpi=plot['dpi'],bbox_inches='tight')
    environment_bins=pd.concat(env_bins,ignore_index=True)
    fits=fit_table(environment_bins)
    tables=dict(structure_bins=structure_bins,structure_normalised_bins=structure_norm_bins,
        structure_members=matched_norm,fit_qc=qc_counts,environment_bins=environment_bins,
        environment_audit=environment_audit,environment_fits=fits)
    for name in ['ar','depth','depth_contrasts','coverage']:tables['alignment_'+name]=alignment[name]
    for cyc,h in alignment['histograms'].items():tables['alignment_histogram_'+cyc]=h
    shape_audit=shapes[KEYS+['Cyc','ShapeDepth','Q_valid','AxisRatio','AlignmentDeg','TiltDis']].copy()
    shape_audit['direction_eligible']=(shape_audit.Q_valid & np.isfinite(shape_audit.AlignmentDeg)
        &np.isfinite(shape_audit.TiltDis)&shape_audit.TiltDis.ge(a['min_tilt_km'])
        &shape_audit.AxisRatio.between(a['min_axis_ratio'],a['max_axis_ratio']))
    tables['alignment_audit']=shape_audit
    population=[]
    for label,frame in [('structure_matched',matched),('structure_normalised',matched_norm.dropna(subset=['tilt_over_surface_Rc'])),
                        ('environment',environment),('environment_normalised',environment_norm.dropna(subset=['tilt_over_surface_Rc']))]:
        for cyc,d in frame.groupby('Cyc'):
            population.append(dict(analysis=label,Cyc=cyc,eddy_days=len(d[KEYS].drop_duplicates()),eddies=d.Eddy.nunique()))
    tables['sample_counts']=pd.DataFrame(population)
    for name,table in tables.items():table.to_csv(out/(name+'.csv'),index=False)
    lines=['# Controls results','',
        'Figures 13/14: fixed actual fitted depths '+', '.join(f'{z:.8f}' for z in c['structure_depths_m'])+
        ' m. The same eddy-days pass both property screens at all three depths. No depth ranking '
        'or interpolation. |Omega| is the fitted angular rotation parameter (s⁻¹), not speed or '
        'relative vorticity. Normalisation uses each day’s surface Rc before summary.','',
        'Figures 13/14/16 show eddy-day medians and IQR, with shared quantile bins across polarities; '
        f'bins require {bins["min_eddies_per_bin"]} distinct eddies. IQR is spread, not uncertainty.','',
        'Figure 13 keeps the original shared-axis styling; its upper limit expands to '
        'include every displayed IQR band (approved correction to inherited clipping).','',
        'Figure 15: actual ellipse depths '+', '.join(f'{z:.8f}' for z in c['ellipse_depths_m'])+
        ' m, without interpolation. Depth curves use identical eddy-days after per-depth directional '
        f'QC: tilt ≥{a["min_tilt_km"]:g} km and {a["min_axis_ratio"]:g}≤AR≤{a["max_axis_ratio"]:g}. '
        f'Each sample/class requires ≥{a["min_days_per_eddy"]} days per eddy and ≥{a["min_eddies"]} eddies '
        'per estimate. Surface histograms/AR classes and depth-matched curves have different '
        'documented populations. Means weight days equally; pointwise 95% intervals resample '
        f'whole eddies ({a["bootstrap_replicates"]} replicates, seed {a["seed"]}). Score cos(2φ) '
        'is +1 parallel and −1 perpendicular. Its zero reference is not an environment-matched '
        'independence test. Surface geometry is processed upstream; deeper geometry uses the '
        'existing fitted levels. Their smoothing/uncertainty can differ. '
        'The three panels overlay AE/CE: surface-angle probability, surface-axis-ratio '
        'alignment, and matched-depth alignment. Panel a shows curves without shading; '
        'b/c retain confidence bands/bars without markers and start at zero.','',
        'Figure 16 displays latitude in positive degrees south, decreasing towards the right; '
        'displayed latitude-fit slopes use this coordinate. Numerical bin/fit tables retain '
        'the original signed degrees north. It uses current source latitude and fresh '
        'surface-fit PV, without dominance '
        'or water-depth selection. Fits and R² describe displayed binned medians, not adjusted '
        'individual-day regressions or causal effects. Beta is df/dy, not beta/h; PV uses log10 '
        'of magnitude in m⁻² s⁻¹. '
        f'Figure 16 display limits: raw {plot["figure16_raw_ylim"]} km (approved upper limit 50 km '
        f'to avoid clipping IQR); normalised {plot["figure16_normalised_ylim"]}.','', 'Sample counts:','']
    for row in tables['sample_counts'].itertuples():
        lines.append(f'- {row.analysis}, {row.Cyc}: {row.eddy_days:,} days / {row.eddies:,} eddies.')
    lines+=['','Ellipse coverage:','']
    for row in alignment['coverage'].itertuples():
        lines.append(f'- {row.Cyc} {row.sample}: {row.observations:,} days / {row.eddies:,} eddies.')
    lines+=['','Surface major-axis alignment (eddy-day mean cos(2φ), pointwise 95% interval):','']
    for row in alignment['depth'][alignment['depth'].depth.eq(0)].itertuples():
        lines.append(f'- {row.Cyc}, depth-matched: {row.estimate:.4f} [{row.low:.4f}, {row.high:.4f}].')
    lines+=['','Writing tables: structure_bins/structure_normalised_bins (medians/IQR/counts), '
        'environment_bins/environment_fits (fits in original x units), alignment_histogram_*, '
        'alignment_ar, alignment_depth, alignment_depth_contrasts and sample_counts. '
        'fit_qc reports overlapping exclusion reasons after valid source-tilt selection; '
        'environment_audit and alignment_audit preserve day/depth eligibility. '
        'See manifest.json for source identities, exact settings and code hashes. No per-day '
        'significance test or mechanism attribution is implied by these associations.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    files=[__file__,ct.__file__,st.__file__,at.__file__,pv.__file__,pv.kernels.__file__,pv.pv_plots.__file__,census.__file__]
    manifest=dict(config=c,config_file=pv.identity(path),pv_config=pv.identity(c['pv_config']),
        census_config=pv.identity(c['census_config']),tilt=pv.identity(cc['inputs']['tilt']),
        pv_product_manifest=json.loads((Path(pvc['product_dir'])/'manifest.json').read_text()),
        code=[pv.identity(f) for f in files],git_revision=revision,python=sys.version,
        numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__,
        created_utc=datetime.now(timezone.utc).isoformat())
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(output_dir=out,figures=figures,tables=tables)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();result=run(args.config);print(result['output_dir'])
