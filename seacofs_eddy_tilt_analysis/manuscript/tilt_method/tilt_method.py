"""Pinned Figure 3: production-equivalent delta fit and transparent step tables."""
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
location=str(ROOT/'pv')
if location not in sys.path:sys.path.insert(0,location)
import pv
import tilt_kernels as kernels
import tilt_plots

KEYS=['Eddy','Day']


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text())
    for key,value in c['inputs'].items():
        p=Path(value).expanduser();c['inputs'][key]=str(p if p.is_absolute() else path.parent/p)
    p=Path(c['output_dir']).expanduser();c['output_dir']=str(p if p.is_absolute() else path.parent/p)
    a=c['calculation'];num=a['window_days']
    if not isinstance(num,int) or isinstance(num,bool) or num<3 or num%2!=1:
        raise ValueError('Use an odd window of at least three calendar days')
    for name in ['depth_interval_m','max_depth_m','min_depth_range_m','variance_epsilon']:
        if not np.isfinite(a[name]) or a[name]<=0:raise ValueError(f'Invalid {name}')
    sigma=a['temporal_sigma_days']
    if sigma is not None and (not np.isfinite(sigma) or sigma<=0):raise ValueError('Invalid temporal sigma')
    if a['min_points']<2 or not np.isfinite(a['bearing_offset_deg']):raise ValueError('Invalid fit controls')
    for name in KEYS:
        value=c['example'][name]
        if not isinstance(value,int) or isinstance(value,bool):raise ValueError('Example keys must be integers')
    return c


def calculate(profiles,day,a):
    """Recompute the illustration without altering the dataset's saved outcomes."""
    pv.check_keys(profiles,KEYS+['Depth'])
    profiles=profiles.assign(Depth=profiles.Depth.abs()).sort_values(['Day','Depth'])
    pv.check_keys(profiles,KEYS+['Depth'])
    dx,dy=kernels.increments(profiles,a['depth_interval_m'],a['max_depth_m'])
    fit=kernels.fit_snapshot(dx,dy,day,sigma=a['temporal_sigma_days'],half_window=a['window_days']//2,
        min_depth_range=a['min_depth_range_m'],min_points=a['min_points'],
        bearing_offset=a['bearing_offset_deg'],eps=a['variance_epsilon'])
    if fit is None or not np.isfinite([fit['TiltDis'],fit['TiltDir']]).all() or fit['TiltDis']<=a['variance_epsilon']:
        raise ValueError('Example has insufficient reference/window support or zero fitted tilt')
    window=profiles[profiles.Depth.le(a['max_depth_m'])].copy()
    ref=window[window.Day.eq(day)].sort_values('Depth')
    top,bottom=fit['ends'];axis=(bottom[:2]-top[:2])/fit['TiltDis']
    raw=ref[['xc','yc']].to_numpy();projection=(raw-raw[0])@axis
    extent=float(np.ptp(projection))
    intervals=[]
    days=np.arange(day-a['window_days']//2,day+a['window_days']//2+1)
    temporal=np.ones(len(days)) if a['temporal_sigma_days'] is None else np.exp(-.5*((days-day)/a['temporal_sigma_days'])**2)
    support=np.isfinite(dx[day])&np.isfinite(dy[day])
    for d,w in zip(days,temporal):
        x=dx.reindex(columns=days).loc[support,d];y=dy.reindex(columns=days).loc[support,d]
        intervals.append(pd.DataFrame(dict(Day=d,upper_depth_m=x.index,
            lower_depth_m=x.index+a['depth_interval_m'],dx_km=x.to_numpy(),dy_km=y.to_numpy(),temporal_weight=w)))
    daily=pd.concat(intervals,ignore_index=True)
    trace=fit['trace'].copy();origin=trace[['x','y']].iloc[0].to_numpy()@axis
    trace['projected_position_km']=trace[['x','y']].to_numpy()@axis-origin
    trace['relative_fit_weight']=trace.weight/trace.weight.max()
    endpoints=pd.DataFrame(fit['ends'],columns=['x_km','y_km','Depth'])
    endpoints.insert(0,'endpoint',['shallow','deep'])
    endpoints['projected_position_km']=fit['ends'][:,:2]@axis-origin
    projected=window.copy()
    first=projected.groupby('Day')[['xc','yc']].transform('first')
    projected['projected_position_km']=np.sum((projected[['xc','yc']]-first).to_numpy()*axis,axis=1)
    summary=dict(Eddy=int(profiles.Eddy.iloc[0]),Day=int(day),TiltDis_km=fit['TiltDis'],TiltDir_deg=fit['TiltDir'],
        MaxProjectedExtent_km=extent,extent_minus_tilt_km=extent-fit['TiltDis'],
        ReferenceDepthMin_m=ref.Depth.min(),ReferenceDepthMax_m=ref.Depth.max(),
        FitDepthMin_m=top[2],FitDepthMax_m=bottom[2],fit_points=len(trace),
        contributing_days_min=int(trace.contributing_days.min()),contributing_days_max=int(trace.contributing_days.max()))
    return dict(fit=fit,summary=summary,trace=trace,endpoints=endpoints,intervals=daily,projected=projected)


def run(config_path=None):
    path=Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    c=load_config(path);a=c['calculation'];eddy,day=[c['example'][k] for k in KEYS]
    half=a['window_days']//2
    profiles=pd.read_parquet(c['inputs']['profiles'],columns=KEYS+['Depth','xc','yc'],
        filters=[('Eddy','=',eddy),('Day','>=',day-half),('Day','<=',day+half)])
    result=calculate(profiles,day,a)
    saved=pd.read_parquet(c['inputs']['tilt']);pv.check_keys(saved)
    saved=saved[saved.Eddy.eq(eddy)&saved.Day.eq(day)]
    if len(saved)!=1:raise ValueError('Example has no unique saved tilt')
    s=result['summary']
    s['SavedTiltDis_km']=float(saved.TiltDis.iloc[0]);s['SavedTiltDir_deg']=float(saved.TiltDir.iloc[0])
    s['distance_difference_km']=s['TiltDis_km']-s['SavedTiltDis_km']
    s['bearing_difference_deg']=(s['TiltDir_deg']-s['SavedTiltDir_deg']+180)%360-180
    with plt.rc_context(matplotlib.rcParamsDefault):
        fig=tilt_plots.plot(eddy,day,profiles,result['fit'],dict(**a,**c['plot']))
    out=Path(c['output_dir']);out.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(out/f'delta_schematic.{ext}',dpi=c['plot']['dpi'],bbox_inches='tight')
    pd.DataFrame([s]).to_csv(out/'example.csv',index=False)
    profiles.to_csv(out/'source_profiles.csv',index=False)
    for name in ['trace','endpoints','intervals','projected']:result[name].to_csv(out/f'{name}.csv',index=False)
    coverage=profiles.assign(within_cap=profiles.Depth.abs().le(a['max_depth_m'])).groupby('Day').agg(
        source_levels=('Depth','size'),within_cap_levels=('within_cap','sum'),source_max_depth_m=('Depth','max'))
    coverage.to_csv(out/'profile_support.csv')
    lines=['# Delta tilt method illustration','',
        f'Figure 3: fixed Eddy {eddy}, reference Day {day}; calendar window {day-half}–{day+half}. '
        'This is an illustration, not a population sample; the directional-analysis ≥5 km '
        'threshold does not apply. Original three-panel layout, colours, curves, brackets and '
        'fit-weight ribbon retained.','',
        f'Recomputed tilt: {s["TiltDis_km"]:.6f} km, bearing {s["TiltDir_deg"]:.6f}°. '
        f'Saved tilt: {s["SavedTiltDis_km"]:.6f} km, {s["SavedTiltDir_deg"]:.6f}°. '
        f'Differences: {s["distance_difference_km"]:.3g} km and {s["bearing_difference_deg"]:.3g}°. '
        f'Reference-day maximum projected extent: {s["MaxProjectedExtent_km"]:.6f} km.','',
        f'Use confirmed centre fits within 0–{a["max_depth_m"]:g} m. Interpolate to '
        f'{a["depth_interval_m"]:g} m levels inside each day’s observed range; no extrapolation. '
        'Keep only intervals supported by the reference day before temporal averaging, '
        'variance weighting and cumulative reconstruction. Neighbouring days cannot extend '
        'the reference depth range.','',
        f'Gaussian temporal sigma {a["temporal_sigma_days"]} days; renormalise available-day '
        'weights at each depth. Depth weights = 1/(sample variance of dx + sample variance '
        f'of dy + {a["variance_epsilon"]:g}), using unweighted sample variance across days. '
        f'Contributing days per retained interval: {s["contributing_days_min"]}–{s["contributing_days_max"]}. '
        'Cumulatively sum mean increments and fit the leading singular vector of the '
        'weighted, centred (x, y, depth) coordinates, retaining production units: horizontal '
        'km and depth m.','',
        f'Reference fitted levels span {s["ReferenceDepthMin_m"]:.2f}–{s["ReferenceDepthMax_m"]:.2f} m '
        f'within the cap. Upper-interval fit labels span {s["FitDepthMin_m"]:g}–{s["FitDepthMax_m"]:g} m '
        f'({s["fit_points"]} points); each label denotes the top of a '
        f'{a["depth_interval_m"]:g} m increment. Minimum label range '
        f'{a["min_depth_range_m"]:g} m and minimum {a["min_points"]} fit points.','',
        'The fitted line’s horizontal endpoint separation gives tilt distance. The arrow '
        'points deep-to-shallow; geographic bearing includes the configured '
        f'{a["bearing_offset_deg"]:g}° offset. Display projection is positive shallow-to-deep '
        'along that line’s horizontal axis. Each daily profile is recentered independently. '
        'The red extent is max minus min reference-day position on this axis; it is distinct '
        'from maximum pairwise 2-D separation or maximum surface-relative distance.','',
        'The blue ribbon encodes relative depth fit weight, not uncertainty. No CI or '
        'significance test applies to this selected illustration. Production dataset outcomes '
        'and scratch products are unchanged.','',
        'Writing tables: example.csv, profile_support.csv. Detailed steps: source_profiles.csv, '
        'intervals.csv (including temporal weights), trace.csv (mean increments, variance, '
        'fit weights and support), endpoints.csv and projected.csv. manifest.json records '
        'inputs, settings and local code hashes.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=path.parent,text=True).strip()
    except (OSError,subprocess.CalledProcessError):revision=None
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),config=c,config_file=pv.identity(path),
        inputs={k:pv.identity(p) for k,p in c['inputs'].items()},
        code=[pv.identity(p) for p in [__file__,kernels.__file__,tilt_plots.__file__,pv.__file__]],
        git_revision=revision,python=sys.version,numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(output_dir=out,figure=fig,example=pd.DataFrame([s]))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();print(run(args.config)['output_dir'])
