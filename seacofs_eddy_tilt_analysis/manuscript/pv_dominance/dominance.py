"""PV populations and Figures 7–10/12; uses only manuscript code and fresh PV."""
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

ROOT = Path(__file__).resolve().parents[1]
for folder in ['pv','census']:
    location = str(ROOT/folder)
    if location not in sys.path: sys.path.insert(0,location)
import pv
import census
import census_plots
import dominance_plots

KEYS = ['Eddy','Day']
REGIMES = ['planetary','mixed','topographic','unclassified']


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text())
    for key in ['pv_config','census_config','output_dir','figure_dir']:
        p=Path(c[key]).expanduser();c[key]=str(p if p.is_absolute() else path.parent/p)
    if not np.isfinite(c['dominance_factor']) or c['dominance_factor']<=1:
        raise ValueError('Dominance factor must be finite and >1')
    if c['direction_min_tilt_km']<0 or not 0<c['alignment_tolerance_deg']<=90:
        raise ValueError('Invalid direction/alignment threshold')
    if c['direction_min_tilt_km']>=10:
        raise ValueError('The original rose magnitude classes require a minimum tilt below 10 km')
    if not np.isfinite(c['water_depth_threshold_m']) or c['water_depth_threshold_m']<=0:
        raise ValueError('Water-depth threshold must be finite and positive')
    return c


def classify(plan,topo,factor=2.):
    """No smoothing/depth conditions; zero/zero and invalid inputs are unclassified."""
    if not np.isfinite(factor) or factor<=1: raise ValueError('Factor must exceed one')
    plan,topo=np.broadcast_arrays(np.asarray(plan,float),np.asarray(topo,float))
    valid=np.isfinite(plan)&np.isfinite(topo)&(plan>=0)&(topo>=0)&((plan>0)|(topo>0))
    return np.select([valid&(plan>=factor*topo),valid&(topo>=factor*plan),valid],
                     ['planetary','topographic','mixed'],default='unclassified')


def prepare(surface,tilt,product,settings):
    for table in [surface,tilt,product]:pv.check_keys(table)
    source_index=pd.MultiIndex.from_frame(surface[KEYS])
    if not pd.MultiIndex.from_frame(product[KEYS]).sort_values().equals(source_index.sort_values()):
        raise ValueError('PV product does not cover exactly the source Eddy-Day keys')
    if not pd.MultiIndex.from_frame(tilt[KEYS]).isin(source_index).all():
        raise ValueError('Tilt has keys absent from surface')
    if not surface.Cyc.isin(['AE','CE']).all():raise ValueError('Unrecognised polarity')
    data=surface.drop(columns=['TiltDis','TiltDir'],errors='ignore').merge(
        tilt[KEYS+['TiltDis','TiltDir']],on=KEYS,how='left',validate='one_to_one')
    overlap=(set(data)&set(product))-set(KEYS)
    if overlap:raise ValueError(f'Unexpected source/PV column overlap: {sorted(overlap)}')
    data=data.merge(product,on=KEYS,validate='one_to_one')
    data['regime']=classify(data.PV_grad_plan_mag,data.PV_grad_topo_mag,settings['dominance_factor'])
    data['dom']=data.regime.map({'planetary':'plan','mixed':'mixed','topographic':'topo'})
    # Depth selects the target environments; it never redefines PV dominance.
    bound=settings['water_depth_threshold_m']
    valid_depth=np.isfinite(data.h)&data.h.gt(0)
    plan=data.regime.eq('planetary');topo=data.regime.eq('topographic')
    deep_plan=plan&valid_depth&data.h.ge(bound)
    shallow_topo=topo&valid_depth&data.h.le(bound)
    data['target_population']=np.select([deep_plan,shallow_topo],
        ['deep_planetary','shallow_topographic'],default='not_targeted')
    data['depth_selection_status']=np.select([
        (plan|topo)&~valid_depth,plan&valid_depth&data.h.lt(bound),
        topo&valid_depth&data.h.gt(bound),deep_plan|shallow_topo],
        ['invalid water depth','planetary but shallower than threshold',
         'topographic but deeper than threshold','depth eligible'],default='not targeted')
    with np.errstate(divide='ignore',invalid='ignore'):
        data['topo_plan_ratio']=np.log(data.PV_grad_topo_mag/data.PV_grad_plan_mag)
        data['log_PV_grad_mag']=np.log(data.PV_grad_mag.where(data.PV_grad_mag>0))
    data['topo_plan_ratio']=data.topo_plan_ratio.replace([np.inf,-np.inf],np.nan)
    data['TiltDir']=data.TiltDir%360
    data['valid_pv_direction']=np.isfinite(data.PV_grad_mag)&data.PV_grad_mag.gt(0)&np.isfinite(data.PV_grad_theta)
    data['direction_valid']=(np.isfinite(data.TiltDis)&data.TiltDis.ge(settings['direction_min_tilt_km'])
                             &np.isfinite(data.TiltDir)&data.valid_pv_direction)
    data['dtheta_PV_grad']=abs((data.TiltDir-data.PV_grad_theta+180)%360-180).where(data.direction_valid)
    data['shelf_class']=np.where(data.lon<settings['shelf_longitude_deg'],'On-shelf','Off-shelf')
    data.loc[~np.isfinite(data.lon),'shelf_class']=None
    return data


def tables(data,settings):
    rows=[];directions=[]
    for cyc in ['AE','CE']:
        group=data[data.Cyc.eq(cyc)];denominator=group.regime.ne('unclassified').sum()
        populations=[(regime,'all_depths',group[group.regime.eq(regime)]) for regime in REGIMES]
        populations += [(regime,label,group[group.target_population.eq(label)]) for regime,label in
                        [('planetary','deep_planetary'),('topographic','shallow_topographic')]]
        for regime,population,d in populations:
            qualified=d[d.direction_valid]
            if population=='all_depths':
                rows.append(dict(Cyc=cyc,regime=regime,eddy_days=len(d),eddies=d.Eddy.nunique(),
                    classified_days=denominator,pct_classified=100*len(d)/denominator if denominator and regime!='unclassified' else np.nan,
                    directional_days=len(qualified),directional_eddies=qualified.Eddy.nunique(),
                    finite_log_ratio_days=np.isfinite(d.topo_plan_ratio).sum()))
            for scope,selected in [('All',qualified)]+[(n,qualified[qualified.shelf_class.eq(n)]) for n in ['On-shelf','Off-shelf']]:
                tilt_mean,tr=census.circular_mean(selected.TiltDir)
                pv_mean,pr=census.circular_mean(selected.PV_grad_theta)
                directions.append(dict(Cyc=cyc,regime=regime,population=population,scope=scope,eddy_days=len(selected),
                    eddies=selected.Eddy.nunique(),tilt_circular_mean_deg=tilt_mean,tilt_R=tr,
                    pv_circular_mean_deg=pv_mean,pv_R=pr,median_tilt_km=selected.TiltDis.median(),
                    median_pv=selected.PV_grad_mag.median(),median_angle_difference_deg=selected.dtheta_PV_grad.median(),
                    aligned_pct=100*selected.dtheta_PV_grad.le(settings['alignment_tolerance_deg']).mean(),
                    opposed_pct=100*selected.dtheta_PV_grad.ge(180-settings['alignment_tolerance_deg']).mean()))
    return pd.DataFrame(rows),pd.DataFrame(directions)


def plot_pv_roses(data,grid,config,census_config):
    selection=data[data.valid_pv_direction & np.isfinite(data.log_PV_grad_mag)].copy()
    # The displayed roses live in model-grid x/y, not geographic east/north.
    offset=np.degrees(grid.angle)
    edges=config['plot']['pv_log_magnitude_edges'];edges=[-np.inf]+edges[1:-1]+[np.inf]
    r=census_config['regions'];f=config['plot']
    fig,_=census_plots.rose_plot(selection,grid,mag='log_PV_grad_mag',theta='PV_grad_theta',
        mag_bins=edges,frac=f['figure7_rose_size_factor'],direction_offset=offset,
        shelf_offset=r['shelf_extension_km'],lon_split=r['longitude_split_deg'],
        lat_split=r['latitude_split_deg'],legend_title=r'log$(|\nabla PV|)$',rtick_flag=False,
        cmaps=['Oranges','Purples'],map_norm_factor=f['figure7_map_normalisation'],show=False,
        figsize=f['figure7_size_inches'],region_bins=census.region_grid(grid,r))
    return fig


def run(config_path=None):
    path=Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    config=load_config(path);pvc=pv.load_config(config['pv_config']);cc=census.load_config(config['census_config'])
    if cc['inputs']['surface']!=pvc['inputs']['surface'] or cc['inputs']['grid']!=pvc['inputs']['grid']:
        raise ValueError('Census and PV configurations refer to different sources/grids')
    product=pv.load_product(pvc);pv.validate(product)
    surface=pd.read_parquet(pvc['inputs']['surface']);tilt=pd.read_parquet(cc['inputs']['tilt'])
    data=prepare(surface,tilt,product,config);grid=pv.load_grid(pvc['inputs']['grid'])
    # Same geographic regions used by the census; no second hidden region boundary.
    regions=census.prepare(surface,tilt,grid,cc)[KEYS+['Region']]
    data=data.drop(columns='Region',errors='ignore').merge(regions,on=KEYS,validate='one_to_one')
    inventory,direction=tables(data,config)
    qualified=data[data.direction_valid]
    planetary=qualified[qualified.target_population.eq('deep_planetary')]
    topographic=qualified[qualified.target_population.eq('shallow_topographic')]
    masks={
        'figure7': data.valid_pv_direction & np.isfinite(data.log_PV_grad_mag) & data.Region.notna(),
        'figure8': np.isfinite(data[['topo_plan_ratio','xc','yc']]).all(axis=1),
        'figure9_category': data.regime.ne('unclassified'),
        'figure9_histogram': data.regime.ne('unclassified') & np.isfinite(data.topo_plan_ratio),
        'figure10': data.direction_valid & data.target_population.eq('deep_planetary'),
        'figure12': data.direction_valid & data.target_population.eq('shallow_topographic')}
    selection_counts=[]
    for figure,mask in masks.items():
        for cyc in ['AE','CE']:
            s=data[mask & data.Cyc.eq(cyc)]
            selection_counts.append(dict(figure=figure,Cyc=cyc,eddy_days=len(s),eddies=s.Eddy.nunique(),
                                         excluded_days=int(data.Cyc.eq(cyc).sum()-len(s))))
    settings=dict(config['plot'],context_isobath_m=config['water_depth_threshold_m'],
        **{k:config[k] for k in ['dominance_factor','direction_min_tilt_km',
        'shelf_longitude_deg','alignment_tolerance_deg']})
    with plt.rc_context(matplotlib.rcParamsDefault):
        figures={'pv_gradient_rose':plot_pv_roses(data,grid,config,cc),
                 'beta_dom_map':dominance_plots.plot_map(data,grid,settings)}
    with plt.style.context('seaborn-v0_8-whitegrid'):
        figures.update(topo_plan_ratio=dominance_plots.plot_distribution(data,grid,settings),
                       plan_dom_rose=dominance_plots.plot_planetary(planetary,grid,settings),
                       topo_dom_dtheta=dominance_plots.plot_topographic(topographic,grid,settings))
    out=Path(config['output_dir']);out.mkdir(parents=True,exist_ok=True)
    figure_dir=Path(config['figure_dir']);figure_dir.mkdir(parents=True,exist_ok=True)
    for name,fig in figures.items():
        fig.savefig(figure_dir/f'{name}.png',dpi=config['plot']['dpi'],bbox_inches='tight')
    audit=data[KEYS+['Cyc','Region','regime','h','target_population','depth_selection_status',
                   'shelf_class','direction_valid','valid_pv_direction','topo_plan_ratio']]
    audit.to_csv(out/'sample_audit.csv',index=False)
    inventory.to_csv(out/'regime_counts.csv',index=False);direction.to_csv(out/'direction_summary.csv',index=False)
    pd.DataFrame(selection_counts).to_csv(out/'figure_samples.csv',index=False)
    depth_rows=[]
    for (cyc,regime,status),group in data.groupby(['Cyc','regime','depth_selection_status']):
        depth_rows.append(dict(Cyc=cyc,regime=regime,depth_selection_status=status,
            eddy_days=len(group),eddies=group.Eddy.nunique(),directional_days=int(group.direction_valid.sum()),
            directional_eddies=group.loc[group.direction_valid,'Eddy'].nunique()))
    depth_counts=pd.DataFrame(depth_rows)
    depth_counts.to_csv(out/'depth_selection_counts.csv',index=False)
    map_data=data[masks['figure8']]
    xb=dominance_plots.bin_edges_fd(map_data.xc.values,grid.X_grid)
    yb=dominance_plots.bin_edges_fd(map_data.yc.values,grid.Y_grid)
    map_rows=[]
    for cyc in ['AE','CE']:
        group=map_data[map_data.Cyc.eq(cyc)]
        mean=dominance_plots.binned_statistic_2d(group.xc,group.yc,group.topo_plan_ratio,
                          statistic='mean',bins=[xb,yb]).statistic
        count=np.histogram2d(group.xc,group.yc,bins=[xb,yb])[0]
        for i,j in zip(*np.where(count>0)):
            map_rows.append(dict(Cyc=cyc,x_low_km=xb[i],x_high_km=xb[i+1],y_low_km=yb[j],
                y_high_km=yb[j+1],eddy_days=int(count[i,j]),mean_log_ratio=mean[i,j]))
    pd.DataFrame(map_rows).to_csv(out/'map_bins.csv',index=False)
    rows=[]
    for cyc,d in data.groupby('Cyc'):
        for name in ['PV_grad_plan_mag','PV_grad_topo_mag','PV_grad_mag']:
            v=d[name].where(np.isfinite(d[name])).dropna()
            rows.append(dict(Cyc=cyc,metric=name,n=len(v),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75)))
    pd.DataFrame(rows).to_csv(out/'gradient_summary.csv',index=False)
    lines=['# PV dominance results','',
        'Surface-fit environmental PV. Categories use raw component magnitudes: planetary ≥ '
        f'{config["dominance_factor"]:g}× topographic; topographic ≥ {config["dominance_factor"]:g}× planetary; '
        'otherwise mixed. Invalid, negative or zero/zero magnitudes are unclassified. No smoothing '
        'in the classification. Percentages use classified days within each polarity.','',
        f'Target populations: deep planetary = planetary dominated and h ≥{config["water_depth_threshold_m"]:g} m; '
        f'shallow topographic = topographic dominated and h ≤{config["water_depth_threshold_m"]:g} m. '
        'Here h is the Gaussian-footprint-mean water depth used in the PV calculation, not profile depth. '
        'Other depth ranges are excluded from these targeted comparisons, never relabelled mixed. '
        'Missing/nonpositive water depths cannot enter a targeted population. Figures 7–9 retain all depths.','',
        f'Direction analyses require tilt ≥{config["direction_min_tilt_km"]:g} km, finite tilt/PV bearings '
        'and positive finite total PV magnitude. Figures 7–9 do not require a valid tilt. '
        'Figure 7 uses a positive finite total PV vector; Figure 8 requires a finite log ratio/position. '
        'Figure 9 histograms require finite log ratios, while category bars use all classified days.',
        '', 'All-depth dominance inventory (before depth selection):','']
    for cyc in ['AE','CE']:
        for row in inventory[(inventory.Cyc==cyc)&inventory.regime.ne('unclassified')].itertuples():
            lines.append(f'- {cyc} {row.regime}: {row.eddy_days:,} days / {row.eddies:,} eddies '
                         f'({row.pct_classified:.2f}%); directional sample {row.directional_days:,} days.')
    lines+=['','Targeted directional results (equal day weight; geographic bearings):','']
    for row in direction[(direction.scope=='All')&direction.population.ne('all_depths')].itertuples():
        lines.append(f'- {row.Cyc} {row.population}: {row.eddy_days:,} days / {row.eddies:,} eddies; '
                     f'mean tilt bearing {row.tilt_circular_mean_deg:.1f}° '
                     f'(R={row.tilt_R:.3f}); mean PV bearing {row.pv_circular_mean_deg:.1f}°; '
                     f'aligned within {config["alignment_tolerance_deg"]:g}°: {row.aligned_pct:.2f}%; '
                     f'opposed within {config["alignment_tolerance_deg"]:g}°: {row.opposed_pct:.2f}%.')
    lines+=['','Depth selection audit (counts overlap in eddies across groups):','']
    for row in depth_counts[depth_counts.depth_selection_status.ne('not targeted')].itertuples():
        lines.append(f'- {row.Cyc} {row.regime}, {row.depth_selection_status}: '
                     f'{row.eddy_days:,} days / {row.eddies:,} eddies; '
                     f'{row.directional_days:,} days also satisfy directional QC.')
    lines+=['','Figure notes:',
        '- 7: original orange/purple regional rose layout. Geographic bearings are rotated into '
        'model-grid coordinates for display; shared census region boundaries. Offshore roses retain '
        'the original peak-normalisation, so their lengths do not compare regional sample sizes.',
        '- 8: spatial mean of each day’s natural-log topographic/planetary magnitude ratio. '
        'Contours describe spatial bin means, not classifications of every day in the bin. '
        'Colours saturate at the original limits; counts and numerical tables retain all values.',
        '- 9: original three-panel distribution/category/magnitude layout; axes expand to include '
        'the unsmoothed distribution instead of hiding tails. Zero-component ratios are not logged.',
        '- 10: deep planetary-dominated days only, with the water-depth and directional criteria '
        'above. Original map/rose layout. Radial limits expand only where necessary to '
        'avoid clipping counts. Green/magenta arrows are circular mean PV/tilt directions.',
        '- 12: shallow topographic-dominated days only, with the water-depth and directional '
        'criteria above. Original map/mismatch/boxplot layout. On/off-shelf remains a '
        f'geographic longitude split at {config["shelf_longitude_deg"]:g}°E; equality is off-shelf. '
        f'The {config["water_depth_threshold_m"]:g} m isobath marks the common depth threshold; '
        'selection uses footprint-mean h rather than the contour at the centre position.','',
        'These are descriptive associations. Repeated days are not independent eddies; no '
        'significance test or causal claim is implied. Per-regime eddy counts overlap across days '
        'and must not be summed as distinct eddies. The tables report counts and estimands explicitly.','',
        'CSV tables: regime_counts, figure_samples (exact populations/exclusions), map_bins '
        '(spatial means/support), depth_selection_counts (retained/excluded days and eddies), '
        'direction_summary (explicit all-depth or target population, and longitude groups), gradient_summary, '
        'sample_audit. PV gradients: m⁻² s⁻¹; distances: km; direction differences: degrees. '
        'See manifest.json for configurations, source/product identity and code hashes.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    dependencies=[__file__,dominance_plots.__file__,pv.__file__,pv.kernels.__file__,pv.pv_plots.__file__,
                  census.__file__,census_plots.__file__]
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),config=config,
        config_file=pv.identity(path),pv_config=pv.identity(config['pv_config']),
        census_config=pv.identity(config['census_config']),tilt_input=pv.identity(cc['inputs']['tilt']),
        product_manifest=json.loads((Path(pvc['product_dir'])/'manifest.json').read_text()),
        code=[pv.identity(f) for f in dependencies],git_revision=revision,
        python=sys.version,numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(output_dir=out,figures=figures,inventory=inventory,direction_summary=direction)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    args=p.parse_args();result=run(args.config);print(result['output_dir'])
