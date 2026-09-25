"""Small plotting/data helpers for the two paper-control notebooks.

Curves are unadjusted eddy-day medians. Depth ranking uses correlations of
per-eddy medians within each polarity on a common Eddy-Day sample.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

KEYS = ['Eddy', 'Day']
COLORS = {'AE': '#b2182b', 'CE': '#2166ac'}
METRICS = ('abs_Omega', 'AR')
SCALE = {'abs_Omega': 1e5, 'AR': 1., 'beta': 1e11, 'log10_PV_grad_mag': 1.}
LABEL = {'abs_Omega': r'$|\Omega|$ ($10^{-5}$ s$^{-1}$)',
         'AR': r'Axis ratio, $\alpha$',
         'beta': r'$\beta$ ($10^{-11}$ m$^{-1}$ s$^{-1}$)',
         'log10_PV_grad_mag': r'$\log_{10}(|\nabla q_{\mathrm{env}}|\,/\,1\,\mathrm{m}^{-2}\mathrm{s}^{-1})$'}
LINESTYLES = ['-', '--', ':', '-.']


def check_keys(df, keys=KEYS):
    if df[keys].isna().any().any() or df.duplicated(keys).any():
        raise ValueError(f'Missing or duplicate keys: {keys}')


def identity(surface):
    check_keys(surface)
    d = surface[KEYS + ['Cyc', 'TiltDis']].copy()
    return d.loc[d.Cyc.isin(COLORS) & np.isfinite(d.TiltDis) & d.TiltDis.ge(0)]


def prepare_depths(vertical, surface, max_depth=1000., rc_max=300.,
                   omega_max=5e-5, ar_max=5.):
    """Screen fits before selecting depths; all flags remain in an audit table."""
    d = vertical.copy()
    d['Depth'] = d.Depth.abs()
    check_keys(d, KEYS + ['Depth'])
    d = d.loc[np.isfinite(d.Depth) & d.Depth.between(0, max_depth)].copy()
    d = d.drop(columns=['Cyc', 'TiltDis'], errors='ignore').merge(
        identity(surface), on=KEYS, how='inner', validate='many_to_one')
    q11, q12, q22 = (d[c].to_numpy(float) for c in ('q11', 'q12', 'q22'))
    disc = np.sqrt((q11-q22)**2 + 4*q12**2)
    lo, hi = (q11+q22-disc)/2, (q11+q22+disc)/2
    valid_q = np.isfinite(lo) & np.isfinite(hi) & (lo > 0)
    d['AR'] = np.sqrt(np.divide(hi, lo, out=np.full(len(d), np.nan), where=valid_q))
    d['abs_Omega'] = d.Omega.abs()
    d['bad_ellipse'] = ~valid_q | ~d.AR.between(1, ar_max)
    d['bad_radius'] = ~np.isfinite(d.Rc) | ~d.Rc.between(0, rc_max, inclusive='right')
    d['bad_rotation'] = ~np.isfinite(d.abs_Omega) | ~d.abs_Omega.between(0, omega_max, inclusive='right')
    d['keep'] = ~d[['bad_ellipse', 'bad_radius', 'bad_rotation']].any(axis=1)
    audit = d.groupby(['Depth', 'Cyc'], as_index=False).agg(
        rows=('Day', 'size'), eddies=('Eddy', 'nunique'), kept=('keep', 'sum'),
        bad_ellipse=('bad_ellipse', 'sum'), bad_radius=('bad_radius', 'sum'),
        bad_rotation=('bad_rotation', 'sum'))
    return d.loc[d.keep, KEYS+['Depth','Cyc','TiltDis',*METRICS]].copy(), audit


def candidate_depths(d, targets, tolerance=65., min_eddies=100, min_day_fraction=.5):
    """Resolve targets globally, then require support in both polarities."""
    levels = np.sort(d.Depth.unique())
    if not len(levels): raise ValueError('No valid fitted depths')
    reference = d.loc[d.Depth.eq(levels[0])].groupby('Cyc').size()
    seen = set(); rows = []
    for target in targets:
        z = float(levels[np.argmin(abs(levels-target))])
        part = d.loc[d.Depth.eq(z)]
        counts = part.groupby('Cyc').agg(eddy_days=('Day','size'), eddies=('Eddy','nunique'))
        good = abs(z-target) <= tolerance and z not in seen
        record = dict(target_m=target, depth_m=z, mismatch_m=abs(z-target))
        for cyc in COLORS:
            n = int(counts.loc[cyc,'eddy_days']) if cyc in counts.index else 0
            e = int(counts.loc[cyc,'eddies']) if cyc in counts.index else 0
            fraction = n/reference.get(cyc, np.inf)
            record.update({cyc+'_eddy_days':n, cyc+'_eddies':e, cyc+'_day_fraction':fraction})
            good = good and e >= min_eddies and fraction >= min_day_fraction
        record['eligible'] = bool(good)
        rows.append(record)
        if good: seen.add(z)
    return pd.DataFrame(rows)


def matched_depths(d, levels):
    """Hold all Eddy-Day observations fixed across the candidate depths."""
    levels = list(levels)
    if len(set(levels)) != len(levels) or len(levels) < 2:
        raise ValueError('Need at least two distinct candidate depths')
    check_keys(d, KEYS+['Depth'])
    use = d.loc[d.Depth.isin(levels)].replace([np.inf,-np.inf],np.nan).dropna(subset=[*METRICS,'TiltDis'])
    counts = use.groupby(KEYS).Depth.nunique()
    keys = counts[counts.eq(len(levels))].reset_index()[KEYS]
    return use.merge(keys, on=KEYS, validate='many_to_one')


def rank_depths(d, min_eddies=100):
    rows = []
    for z, part in d.groupby('Depth'):
        for metric in METRICS:
            row = dict(metric=metric, depth_m=z)
            good = True
            for cyc in COLORS:
                e = part.loc[part.Cyc.eq(cyc)].groupby('Eddy')[[metric,'TiltDis']].median().dropna()
                r = (spearmanr(e[metric], e.TiltDis).statistic
                     if len(e) >= min_eddies and e[metric].nunique()>1 and e.TiltDis.nunique()>1 else np.nan)
                row.update({cyc+'_rho':r, cyc+'_eddies':len(e)})
                good = good and np.isfinite(r)
            row['mean_abs_rho'] = np.mean([abs(row[c+'_rho']) for c in COLORS]) if good else np.nan
            rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty: raise ValueError('No matched eddy-days; reduce candidate depths or revise support settings')
    return table.sort_values(['metric','mean_abs_rho','depth_m'], ascending=[True,False,True])


def select_depths(ranking, n=3, min_spacing=75., override=None):
    if not 1 <= n <= len(LINESTYLES): raise ValueError('Select between one and four depths')
    selected = {}
    for metric in METRICS:
        eligible = ranking.loc[ranking.metric.eq(metric) & ranking.mean_abs_rho.notna()]
        if override is not None:
            chosen = [float(z) for z in override[metric]]
            if len(chosen)!=len(set(chosen)) or not set(chosen).issubset(set(eligible.depth_m)):
                raise ValueError('Overrides must be distinct, eligible actual depths printed in the ranking')
            if not 1 <= len(chosen) <= len(LINESTYLES): raise ValueError('Choose one to four depths')
        else:
            chosen=[]
            for z in eligible.depth_m:
                if all(abs(z-other) >= min_spacing for other in chosen): chosen.append(float(z))
                if len(chosen)==n: break
        if not chosen: raise ValueError(f'No eligible depth for {metric}')
        selected[metric] = sorted(chosen)
    return selected


def bin_edges(d, metric, bins=8):
    x = d[metric].to_numpy(float); x=x[np.isfinite(x)]
    if not len(x): raise ValueError(f'No finite values for {metric}')
    edges = np.unique(np.quantile(x, np.linspace(0,1,bins+1)))
    if len(edges)<2: raise ValueError(f'Constant predictor: {metric}')
    return edges


def median_table(d, metric, edges, min_eddies=20):
    """Shared quantile edges, daily median/IQR; keep gaps for suppressed bins."""
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=[metric,'TiltDis']).copy()
    d['_bin']=pd.cut(d[metric],edges,include_lowest=True,labels=False)
    rows=[]
    for cyc in COLORS:
        for i in range(len(edges)-1):
            g=d.loc[d.Cyc.eq(cyc)&d._bin.eq(i)]
            enough = g.Eddy.nunique() >= min_eddies
            rows.append(dict(Cyc=cyc,bin=i,eddy_days=len(g),eddies=g.Eddy.nunique(),
                             x=g[metric].median(),median=g.TiltDis.median() if enough else np.nan,
                             q25=g.TiltDis.quantile(.25) if enough else np.nan,
                             q75=g.TiltDis.quantile(.75) if enough else np.nan,
                             displayed=enough))
    return pd.DataFrame(rows)


def style():
    plt.rcParams.update({'font.size':9,'axes.labelsize':10,'axes.titlesize':10,
                         'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})


def curve(ax, table, metric, shade=True, linestyle='-'):
    for cyc in COLORS:
        g=table.loc[table.Cyc.eq(cyc)]
        x=g.x.to_numpy()*SCALE[metric]
        ax.plot(x,g['median'],color=COLORS[cyc],ls=linestyle,lw=1.7,marker='o',ms=2.7)
        if shade: ax.fill_between(x,g.q25,g.q75,color=COLORS[cyc],alpha=.10,lw=0)
    ax.set_xlabel(LABEL[metric]); ax.set_ylim(bottom=0)
    ax.grid(axis='y',alpha=.12,lw=.5)


def polarity_handles():
    return [Line2D([0],[0],color=color,lw=1.8,label=cyc) for cyc,color in COLORS.items()]


def selected_depth_figure(d,selected,bins=8,min_eddies=20):
    style(); ncols=max(map(len,selected.values()))
    fig,axs=plt.subplots(2,ncols,figsize=(3.15*ncols,5.5),squeeze=False,sharey=True,sharex='row',layout='constrained')
    tables=[]
    for row,metric in enumerate(METRICS):
        # Identical bins across depths for a given property, and across polarity.
        edges=bin_edges(d.loc[d.Depth.isin(selected[metric])],metric,bins)
        for col,z in enumerate(selected[metric]):
            ax=axs[row,col]; part=d.loc[d.Depth.eq(z)]
            t=median_table(part,metric,edges,min_eddies)
            tables.append(t.assign(metric=metric,depth_m=z))
            curve(ax,t,metric,shade=True)
            ax.set_title(f'({chr(97+row*ncols+col)}) {z:.0f} m',loc='left')
        for ax in axs[row,len(selected[metric]):]:ax.set_visible(False)
        axs[row,0].set_ylabel('Tilt distance (km)')
    axs[0,0].legend(handles=polarity_handles(),frameon=False)
    return fig,pd.concat(tables,ignore_index=True)


def overlay_figure(d,selected,bins=8,min_eddies=20):
    style();fig,axs=plt.subplots(1,2,figsize=(8,3.45),sharey=True,layout='constrained')
    tables=[]
    for i,(ax,metric) in enumerate(zip(axs,METRICS)):
        edges=bin_edges(d.loc[d.Depth.isin(selected[metric])],metric,bins)
        handles=[]
        for j,z in enumerate(selected[metric]):
            t=median_table(d.loc[d.Depth.eq(z)],metric,edges,min_eddies)
            tables.append(t.assign(metric=metric,depth_m=z))
            curve(ax,t,metric,shade=False,linestyle=LINESTYLES[j])
            handles.append(Line2D([0],[0],color='.3',ls=LINESTYLES[j],lw=1.5,label=f'{z:.0f} m'))
        ax.legend(handles=polarity_handles()+handles,frameon=False,ncol=2)
        ax.set_title(f'({chr(97+i)})',loc='left')
    axs[0].set_ylabel('Tilt distance (km)')
    return fig,pd.concat(tables,ignore_index=True)


def environment_data(surface,pv):
    check_keys(pv)
    d=identity(surface).merge(pv[KEYS+['beta','PV_grad_mag']],on=KEYS,how='left',validate='one_to_one')
    d['log10_PV_grad_mag']=np.log10(d.PV_grad_mag.where(d.PV_grad_mag>0))
    d=d.replace([np.inf,-np.inf],np.nan)
    # Both panels use the same observations.
    d['keep']=d[['beta','log10_PV_grad_mag']].notna().all(axis=1)
    audit=d.groupby('Cyc',as_index=False).agg(tilt_days=('Day','size'),kept_days=('keep','sum'))
    return d.loc[d.keep].copy(),audit


def environment_figure(d,bins=8,min_eddies=20):
    style();fig,axs=plt.subplots(1,2,figsize=(8,3.35),sharey=True,layout='constrained')
    tables=[]
    for i,(ax,metric) in enumerate(zip(axs,['beta','log10_PV_grad_mag'])):
        t=median_table(d,metric,bin_edges(d,metric,bins),min_eddies)
        curve(ax,t,metric,shade=True);tables.append(t.assign(metric=metric))
        ax.set_title(f'({chr(97+i)})',loc='left')
    axs[0].set_ylabel('Tilt distance (km)');axs[0].legend(handles=polarity_handles(),frameon=False)
    return fig,pd.concat(tables,ignore_index=True)


def save_figure(fig,output,name):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','png'):fig.savefig(output/f'{name}.{ext}',bbox_inches='tight')


def save_metadata(output,settings,inputs):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    sources=[]
    for path in inputs:
        p=Path(path)
        sources.append(dict(path=str(p),exists=p.exists(),size=p.stat().st_size if p.exists() else None,
                            mtime=p.stat().st_mtime if p.exists() else None))
    (output/'settings.json').write_text(json.dumps(dict(settings=settings,sources=sources),indent=2,default=str)+'\n')
