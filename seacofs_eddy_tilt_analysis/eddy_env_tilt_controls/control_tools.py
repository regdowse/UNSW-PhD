"""Cache-only exploratory controls of measured tilt; no tilt re-estimation."""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.integrate import trapezoid

KEYS = ['Eddy', 'Day']
COLORS = {'AE': 'darkred', 'CE': 'navy'}
PROFILE_METRICS = ['abs_zeta', 'abs_Omega', 'AR', 'Rc', 'rotation_speed_scale', 'Ro_abs']
N2_METRICS = ['N2_200m_core_s2', 'N2_500m_core_s2',
              'N2_200m_integral_m_s2_core', 'N2_500m_integral_m_s2_core',
              'N2_200m_max_s2_core', 'N2_500m_max_s2_core',
              'N2_pycnocline_mean_s2_core', 'N2_pycnocline_max_s2_core',
              'pycnocline_depth_m_core', 'MLD_density_m_core']
LABELS = {'abs_zeta': '|relative vorticity| (s⁻¹)', 'abs_Omega': '|Omega| (s⁻¹)',
          'AR': 'Axis ratio', 'Rc': 'Core radius (km)', 'Ro_abs': '|zeta/f|',
          'rotation_speed_scale': '|Omega| Rc (m/s; speed scale)',
          'TiltDis': 'Tilt distance (km)', 'tilt_over_radius': 'Tilt distance / surface Rc',
          'vertical_extent_m': 'Deepest fitted level (m)', 'age_fraction': 'Fractional observed lifetime',
          'age_days': 'Days since first track observation', 'lat': 'Latitude (°)',
          'beta': 'Planetary beta (m⁻¹ s⁻¹)', 'h': 'Water depth (m)',
          'slope': 'Magnitude of mean bathymetric slope (m/m)',
          'PV_grad_coherence': 'Environmental PV-gradient coherence',
          'profile_straightness': 'Centreline chord / path length',
          'profile_bend_km': 'Maximum departure from endpoint line (km)',
          'profile_reversal_fraction': 'Fraction of centre steps opposing endpoint displacement',
          'shear_mag_ms': 'Upper-minus-deep background velocity difference (m/s)',
          'shear_angle_deg': 'Tilt–background-shear angle (°)',
          'vv_mean': 'Sample-weighted upward velocity (m/s)',
          'vv_rms': 'Sample-weighted vertical-velocity RMS (m/s)',
          'vv_max': 'Largest signed upward velocity (m/s)',
          'vv_min': 'Smallest signed upward velocity (m/s)',
          'vv_abs_max': 'Largest |vertical velocity| (m/s)'}


LABELS.update({
    'log10_PV_grad_mag': 'log10 environmental |∇PV| (m⁻² s⁻¹)',
    'log10_PV_grad_plan_mag': 'log10 planetary |∇PV| (m⁻² s⁻¹)',
    'log10_PV_grad_topo_mag': 'log10 topographic |∇PV| (m⁻² s⁻¹)',
    'log10_topo_plan_ratio': 'log10 (topographic / planetary magnitude)',
    'N2_200m_core_s2': 'Core mean N², 0–200 m (s⁻²)',
    'N2_500m_core_s2': 'Core mean N², 0–500 m (s⁻²)',
    'N2_200m_integral_m_s2_core': 'Core mean integrated N², 0–200 m (m s⁻²)',
    'N2_500m_integral_m_s2_core': 'Core mean integrated N², 0–500 m (m s⁻²)',
    'N2_200m_max_s2_core': 'Core mean maximum N², 0–200 m (s⁻²)',
    'N2_500m_max_s2_core': 'Core mean maximum N², 0–500 m (s⁻²)',
    'N2_pycnocline_mean_s2_core': 'Pycnocline-window mean N² (s⁻²)',
    'N2_pycnocline_max_s2_core': 'Pycnocline maximum N² (s⁻²)',
    'pycnocline_depth_m_core': 'Core mean pycnocline depth (m)',
    'MLD_density_m_core': 'Core mean mixed-layer depth (m)',
    'vv_deep_minus_upper': 'Deep minus upper core-mean vertical velocity (m/s)',
})
for _metric in PROFILE_METRICS:
    for _suffix, _description in [('_depth_mean','Depth mean'),('_max','Profile maximum'),
                                  ('_deep_minus_shallow','Deep minus shallow')]:
        LABELS[_metric+_suffix]=_description+'\n'+LABELS[_metric]
    LABELS[_metric+'_peak_depth_m']='Depth of maximum '+_metric+' (m)'


def unique(frame, keys=KEYS):
    if frame.duplicated(keys).any():
        raise ValueError(f'Duplicate {keys}: {frame.loc[frame.duplicated(keys, False), keys].head()}')


def merge(frame, cache):
    unique(frame); unique(cache)
    overlap = (set(frame) & set(cache)) - set(KEYS)
    if overlap:
        raise ValueError(f'Ambiguous cache columns: {sorted(overlap)}')
    return frame.merge(cache, on=KEYS, how='left', validate='one_to_one')


def finite(frame, columns):
    return frame.replace([np.inf, -np.inf], np.nan).dropna(subset=columns)


def features(frame):
    out = frame.copy()
    for col in ['w', 'Omega', 'Rc', 'q11', 'q12', 'q22']:
        if col not in out:
            raise ValueError(f'Missing fitted property {col}')
    a, b, c = (out[k].to_numpy(float) for k in ('q11', 'q12', 'q22'))
    rad = np.sqrt((a-c)**2 + 4*b*b)
    lo, hi = (a+c-rad)/2, (a+c+rad)/2
    valid = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (out.Rc > 0)
    out['AR'] = np.sqrt(np.divide(hi, lo, out=np.full(len(out), np.nan), where=valid))
    out['abs_zeta'] = out.w.abs()
    out['abs_Omega'] = out.Omega.abs()
    # A dimensional rotation-speed scale, not a reconstructed maximum velocity.
    out['rotation_speed_scale'] = out.abs_Omega * out.Rc * 1000
    out['Ro_abs'] = out.abs_zeta / out.f.abs().where(out.f.abs() > 1e-8)
    out.loc[~valid, PROFILE_METRICS] = np.nan
    return out


def prepare_surface(surface, grid):
    unique(surface)
    out = surface.copy()
    # Age in the processed table is duration, not elapsed age.
    days = out.groupby('Eddy').Day
    out['age_days'] = out.Day - days.transform('min')
    duration = days.transform('max') - days.transform('min')
    out['age_fraction'] = out.age_days / duration.where(duration > 0)
    out['track_duration_days'] = duration
    out['f'] = grid.f[out.ic.to_numpy(int), out.jc.to_numpy(int)]
    out = features(out)
    out['tilt_over_radius'] = out.TiltDis / out.Rc.where(out.Rc > 0)
    if 'Date' in out:
        month = pd.to_datetime(out.Date, errors='coerce').dt.month
        out['Season'] = month.map({12:'DJF',1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',
                                  6:'JJA',7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON'})
    out.loc[~np.isfinite(out.TiltDis) | (out.TiltDis < 0), ['TiltDis','tilt_over_radius']] = np.nan
    return out.loc[out.Cyc.isin(COLORS)].copy()


def prepare_profiles(vertical, surface, max_depth=1000):
    out = vertical.copy()
    out['Depth'] = out.Depth.abs()
    unique(out, KEYS + ['Depth'])
    out = out.loc[out.Depth.between(0, max_depth)].copy()
    meta = [c for c in ['Cyc','TiltDis','TiltDir','tilt_over_radius','f','Region','Season','lat'] if c in surface]
    out = out.drop(columns=meta, errors='ignore').merge(surface[KEYS+meta], on=KEYS,
                                                        how='inner', validate='many_to_one')
    return features(out).sort_values(KEYS+['Depth']).reset_index(drop=True)


def profile_summaries(profiles):
    rows = []
    for key, g in profiles.groupby(KEYS, sort=False):
        g = g.sort_values('Depth')
        z = g.Depth.to_numpy(float)
        row = dict(zip(KEYS, key))
        row.update(vertical_extent_m=z[-1], profile_min_depth_m=z[0], profile_levels=len(g))
        for col in PROFILE_METRICS:
            vals = g[col].to_numpy(float)
            ok = np.isfinite(vals)
            row[col+'_max'] = np.max(vals[ok]) if ok.any() else np.nan
            row[col+'_peak_depth_m'] = z[np.flatnonzero(ok)[np.argmax(vals[ok])]] if ok.any() else np.nan
            # Require the complete sampled profile for means/endpoint changes.
            complete = ok.all() and len(z) >= 2 and z[-1] > z[0]
            row[col+'_depth_mean'] = trapezoid(vals,z)/(z[-1]-z[0]) if complete else np.nan
            row[col+'_deep_minus_shallow'] = vals[-1]-vals[0] if complete else np.nan
        xy = g[['xc','yc']].to_numpy(float)
        for c in ['profile_straightness','profile_bend_km','profile_reversal_fraction']:
            row[c] = np.nan
        if len(xy) >= 3 and np.isfinite(xy).all():
            step = np.diff(xy,axis=0); chord = xy[-1]-xy[0]
            length = np.linalg.norm(chord); path = np.linalg.norm(step,axis=1).sum()
            if length > 1e-8 and path > 1e-8:
                direction = chord / length
                delta = xy-xy[0]
                normal = delta - np.outer(delta @ direction, direction)
                row['profile_straightness'] = length/path
                row['profile_bend_km'] = np.linalg.norm(normal,axis=1).max()
                nonzero = np.linalg.norm(step,axis=1) > 1e-8
                row['profile_reversal_fraction'] = np.mean((step[nonzero] @ direction) < 0)
        rows.append(row)
    return pd.DataFrame(rows, columns=KEYS if not rows else None)


def resolve_depths(profiles, targets=(0,100,200,500,1000), tolerance=65):
    levels = np.sort(profiles.Depth.dropna().unique())
    if not len(levels): raise ValueError('No profile depths')
    rows=[]
    for target in targets:
        actual=levels[np.argmin(abs(levels-target))]
        rows.append(dict(target_m=target, actual_m=actual, mismatch_m=abs(actual-target),
                         accepted=bool(abs(actual-target)<=tolerance)))
    audit=pd.DataFrame(rows)
    accepted=audit.loc[audit.accepted, 'actual_m']
    if accepted.duplicated().any(): raise ValueError('Targets resolve to duplicate levels; change targets')
    return audit


def depth_sample(profiles, levels, metric, matched=False):
    out=finite(profiles.loc[profiles.Depth.isin(levels)], [metric,'TiltDis'])
    if matched:
        counts=out.groupby(KEYS).Depth.nunique()
        keys=counts[counts.eq(len(levels))].reset_index()[KEYS]
        out=out.merge(keys,on=KEYS,validate='many_to_one')
    return out


def add_n2(surface, cache, min_coverage=.8):
    unique(cache)
    table=cache[KEYS+N2_METRICS].copy()
    for col in N2_METRICS:
        coverage = col.replace('_core_s2','_core_valid_fraction') if col.endswith('_core_s2') else col+'_valid_fraction'
        if coverage not in cache: raise ValueError(f'Missing stratification coverage: {coverage}')
        table.loc[~cache[coverage].ge(min_coverage),col]=np.nan
    return merge(surface,table)


def vertical_velocity(raw, fraction=1.5, max_depth=1000, min_cells=8, min_coverage=.7):
    required=set(KEYS+['Depth','fraction','n_valid','coverage','w_mean','w_rms','w_max','w_min'])
    if required-set(raw): raise ValueError(f'Missing vertical-velocity columns: {required-set(raw)}')
    d=raw.loc[np.isclose(raw.fraction,fraction)&raw.Depth.between(0,max_depth)&
              raw.n_valid.ge(min_cells)&raw.coverage.ge(min_coverage)].copy()
    unique(d,KEYS+['Depth'])
    d=finite(d,['w_mean','w_rms','w_max','w_min'])
    d=d.rename(columns={c:'vv_'+c[2:] for c in ['w_mean','w_rms','w_max','w_min']})
    d['vv_abs_max']=d[['vv_max','vv_min']].abs().max(axis=1)
    rows=[]
    for key,g in d.groupby(KEYS):
        row=dict(zip(KEYS,key))
        row.update(vv_mean=np.average(g.vv_mean,weights=g.n_valid),
                   vv_rms=np.sqrt(np.average(g.vv_rms**2,weights=g.n_valid)),
                   vv_max=g.vv_max.max(),vv_min=g.vv_min.min(),vv_abs_max=g.vv_abs_max.max(),
                   vv_levels=len(g),vv_max_depth_m=g.Depth.max(),vv_min_coverage=g.coverage.min())
        upper=g.loc[g.Depth.between(100,300)]; deep=g.loc[g.Depth.between(500,900)]
        row['vv_deep_minus_upper']= (np.average(deep.vv_mean,weights=deep.n_valid)-
                                    np.average(upper.vv_mean,weights=upper.n_valid)) if len(upper) and len(deep) else np.nan
        rows.append(row)
    return d,pd.DataFrame(rows,columns=KEYS if not rows else None)


def add_environment(surface, pv):
    cols=['h','beta','dhdx','dhdy','PV_grad_mag','PV_grad_plan_mag','PV_grad_topo_mag','PV_grad_coherence']
    missing=set(cols)-set(pv)
    if missing: raise ValueError(f'Missing environmental PV fields: {missing}')
    out=merge(surface.drop(columns=[c for c in cols if c in surface]),pv[KEYS+cols])
    out['slope']=np.hypot(out.dhdx,out.dhdy)
    for col in ['PV_grad_mag','PV_grad_plan_mag','PV_grad_topo_mag']:
        out['log10_'+col]=np.log10(out[col].where(out[col]>0))
    out['log10_topo_plan_ratio']=out.log10_PV_grad_topo_mag-out.log10_PV_grad_plan_mag
    return out


def add_shear(surface, background):
    cols=['clim_surface_east_ms','clim_surface_north_ms','clim_500_east_ms','clim_500_north_ms']
    out=merge(surface,background[KEYS+cols])
    east=out[cols[0]]-out[cols[2]]; north=out[cols[1]]-out[cols[3]]
    out['shear_mag_ms']=np.hypot(east,north)
    direction=np.degrees(np.arctan2(east,north))%360
    out['shear_angle_deg']=abs((out.TiltDir-direction+180)%360-180)
    out.loc[(out.shear_mag_ms<=1e-12)|(out.TiltDis<=0),'shear_angle_deg']=np.nan
    return out.drop(columns=cols)


def audit(frame, metrics):
    rows=[]
    for metric in metrics:
        if metric not in frame:
            rows.append(dict(metric=metric,status='unavailable')); continue
        for cyc,g in frame.groupby('Cyc'):
            d=finite(g,[metric,'TiltDis'])
            rows.append(dict(metric=metric,Cyc=cyc,eddy_days=len(d),eddies=d.Eddy.nunique(),
                             missing=len(g)-len(d),status='available'))
    return pd.DataFrame(rows)


def binned(frame,x,y='TiltDis',bins=8,min_eddies=10,n_boot=300,seed=731):
    """Daily median/IQR; pointwise CI resamples whole contributing eddies per bin."""
    d=finite(frame,[x,y]); rows=[]; rng=np.random.default_rng(seed)
    if d.empty:return pd.DataFrame()
    edges=np.unique(np.quantile(d[x],np.linspace(0,1,bins+1)))
    if len(edges)<2:return pd.DataFrame()
    d=d.copy(); d['_bin']=pd.cut(d[x],edges,include_lowest=True,labels=False)
    for (cyc,bin_id),g in d.groupby(['Cyc','_bin'],observed=True):
        groups=[a[y].to_numpy() for _,a in g.groupby('Eddy')]
        if len(groups)<min_eddies:continue
        values=[]
        for _ in range(n_boot):
            values.append(np.median(np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))])))
        lo,hi=np.quantile(values,[.025,.975]) if values else (np.nan,np.nan)
        rows.append(dict(Cyc=cyc,bin=int(bin_id),x=g[x].median(),median=g[y].median(),
                         q25=g[y].quantile(.25),q75=g[y].quantile(.75),ci_low=lo,ci_high=hi,
                         eddy_days=len(g),eddies=len(groups)))
    return pd.DataFrame(rows)


def panels(frame,metrics,y='TiltDis',n_boot=300,min_eddies=10,title=''):
    metrics=[m for m in metrics if m in frame and finite(frame,[m,y]).shape[0]]
    tables=[]
    # Small pages avoid unreadable figures for the full metric collection.
    for start in range(0,len(metrics),6):
        chunk=metrics[start:start+6]
        fig,axes=plt.subplots(int(np.ceil(len(chunk)/3)),min(3,len(chunk)),
                              figsize=(5*min(3,len(chunk)),3.8*int(np.ceil(len(chunk)/3))),squeeze=False)
        for ax,x in zip(axes.flat,chunk):
            t=binned(frame,x,y,n_boot=n_boot,min_eddies=min_eddies)
            if len(t):
                for cyc,g in t.groupby('Cyc'):
                    ax.fill_between(g.x,g.q25,g.q75,color=COLORS[cyc],alpha=.12)
                    ax.plot(g.x,g['median'],'o-',color=COLORS[cyc],label=cyc,ms=3)
                    ax.vlines(g.x,g.ci_low,g.ci_high,color=COLORS[cyc],lw=1)
                tables.append(t.assign(metric=x,outcome=y))
                ax.legend(frameon=False)
            else:ax.text(.5,.5,'Insufficient distinct eddies per bin',ha='center',transform=ax.transAxes)
            ax.set(xlabel=LABELS.get(x,x),ylabel=LABELS.get(y,y));ax.grid(alpha=.15)
        for ax in list(axes.flat)[len(chunk):]:ax.set_visible(False)
        fig.suptitle(title+'\nDaily median; shade = IQR; bars = pointwise eddy-bootstrap 95% CI')
        fig.tight_layout();plt.show()
    return pd.concat(tables,ignore_index=True) if tables else pd.DataFrame()


def correlations(frame,metrics,y='TiltDis',n_boot=300,seed=731):
    """Equal-eddy sensitivity: correlations of per-eddy medians, bootstrap eddies."""
    rng=np.random.default_rng(seed);rows=[]
    for metric in metrics:
        if metric not in frame:continue
        for cyc,g in frame.groupby('Cyc'):
            e=finite(g,[metric,y]).groupby('Eddy')[[metric,y]].median()
            if len(e)<10 or e[metric].nunique()<2 or e[y].nunique()<2:continue
            draws=[]
            for _ in range(n_boot):
                b=e.iloc[rng.integers(len(e),size=len(e))]
                if b[metric].nunique()>1 and b[y].nunique()>1:draws.append(spearmanr(b[metric],b[y]).statistic)
            lo,hi=np.quantile(draws,[.025,.975]) if draws else (np.nan,np.nan)
            rows.append(dict(metric=metric,Cyc=cyc,outcome=y,eddies=len(e),
                             rho=spearmanr(e[metric],e[y]).statistic,ci_low=lo,ci_high=hi))
    return pd.DataFrame(rows)


def adjusted(frame,metrics,y='TiltDis',min_eddies=30):
    """One focal predictor per polarity; equal-eddy WLS and clustered errors."""
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    rows=[]
    for metric in metrics:
        if metric not in frame:continue
        for cyc,g in frame.groupby('Cyc'):
            controls=[c for c in ['Rc','vertical_extent_m','lat'] if c in g and c!=metric]
            # beta and latitude describe the same geographic gradient.
            if metric=='beta':controls=[c for c in controls if c!='lat']
            cats=[c for c in ['Region','Season'] if c in g]
            needed=list(dict.fromkeys([metric,y,'Eddy']+controls+cats))
            d=finite(g[needed],[metric,y,'Eddy']+controls+cats).copy()
            d=d.loc[d[y]>=0]
            record=dict(metric=metric,Cyc=cyc,outcome=y,eddy_days=len(d),eddies=d.Eddy.nunique())
            if d.Eddy.nunique()<min_eddies or d[metric].nunique()<3:
                rows.append(dict(record,status='insufficient data'));continue
            X=pd.DataFrame(index=d.index)
            for c in [metric]+controls:
                sd=d[c].std()
                if sd>0:X['focal' if c==metric else c]=(d[c]-d[c].mean())/sd
            X=pd.concat([X,pd.get_dummies(d[cats],drop_first=True,dtype=float)],axis=1) if cats else X
            X=sm.add_constant(X,has_constant='add').astype(float)
            if np.linalg.matrix_rank(X.to_numpy())<X.shape[1]:
                rows.append(dict(record,status='rank deficient; inspect correlated controls'));continue
            weights=1/d.groupby('Eddy').Eddy.transform('size')
            fit=sm.WLS(np.log1p(d[y]),X,weights=weights).fit(cov_type='cluster',cov_kwds={'groups':d.Eddy})
            ci=fit.conf_int().loc['focal']
            rows.append(dict(record,status='ok',slope=fit.params.focal,ci_low=ci.iloc[0],ci_high=ci.iloc[1],
                             p=fit.pvalues.focal,controls=', '.join(controls+cats)))
    out=pd.DataFrame(rows)
    if len(out) and 'p' in out:
        ok=out.p.notna();out.loc[ok,'q_BH']=multipletests(out.loc[ok,'p'],method='fdr_bh')[1]
    return out


def save_results(directory,tables,settings,inputs):
    directory=Path(directory)/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    directory.mkdir(parents=True,exist_ok=False)
    for name,table in tables.items():table.to_csv(directory/(name+'.csv'),index=False)
    provenance=[]
    for p in inputs:
        p=Path(p);provenance.append(dict(path=str(p),exists=p.exists(),
            size=p.stat().st_size if p.exists() else None,mtime=p.stat().st_mtime if p.exists() else None))
    (directory/'settings.json').write_text(json.dumps(dict(settings=settings,inputs=provenance,tables=list(tables)),indent=2,default=str)+'\n')
    return directory


def depth_comparison(profiles,levels,metrics=PROFILE_METRICS,n_boot=300,min_eddies=10):
    tables=[];counts=[]
    for metric in metrics:
        if metric not in profiles:continue
        fig,axes=plt.subplots(2,len(levels),figsize=(4*len(levels),7),squeeze=False,sharey=True)
        for row,matched in enumerate([False,True]):
            d=depth_sample(profiles,levels,metric,matched)
            for ax,level in zip(axes[row],levels):
                g=d.loc[d.Depth.eq(level)]
                for cyc,part in g.groupby('Cyc'):
                    counts.append(dict(metric=metric,matched=matched,depth_m=level,Cyc=cyc,
                                       eddy_days=len(part),eddies=part.Eddy.nunique()))
                t=binned(g,metric,n_boot=n_boot,min_eddies=min_eddies)
                if len(t):
                    for cyc,p in t.groupby('Cyc'):
                        ax.fill_between(p.x,p.q25,p.q75,color=COLORS[cyc],alpha=.12)
                        ax.plot(p.x,p['median'],'o-',color=COLORS[cyc],ms=3,label=cyc)
                        ax.vlines(p.x,p.ci_low,p.ci_high,color=COLORS[cyc],lw=1)
                    tables.append(t.assign(metric=metric,matched=matched,depth_m=level))
                ax.set(title=f'{level:g} m · {g.Eddy.nunique()} eddies',xlabel=LABELS.get(metric,metric))
                if len(t):ax.legend(frameon=False)
            axes[row,0].set_ylabel(('Matched days' if matched else 'Available days')+'\nTilt distance (km)')
        fig.suptitle('Depth comparison: '+LABELS.get(metric,metric)+'\nShade = IQR; bars = eddy-bootstrap 95% CI')
        fig.tight_layout();plt.show()
    return (pd.concat(tables,ignore_index=True) if tables else pd.DataFrame(),pd.DataFrame(counts))
