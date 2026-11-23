"""Small plotting/data helpers for the two paper-control notebooks.

Curves are unadjusted eddy-day medians. Fixed fitted depths are configured explicitly; no depth ranking.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from scipy.stats import spearmanr

KEYS = ['Eddy', 'Day']
COLORS = {'AE': '#b2182b', 'CE': '#2166ac'}
METRICS = ('abs_Omega', 'AR')
SCALE = {'lat': 1., 'abs_Omega': 1e5, 'AR': 1., 'beta': 1., 'log10_PV_grad_mag': 1.}
LABEL = {'abs_Omega': r'$|\Omega|$ ($10^{-5}$ s$^{-1}$)',
         'AR': r'Axis ratio', 'lat': r'Latitude (°N)',
         'beta': r'$\beta$ (m$^{-1}$ s$^{-1}$)',
         'log10_PV_grad_mag': r'$\log_{10}(|\nabla PV|) \,(\mathrm{m}^{-2}\ \mathrm{s}^{-1})$'}
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

def bin_edges(d, metric, bins=8):
    x = d[metric].to_numpy(float); x=x[np.isfinite(x)]
    if not len(x): raise ValueError(f'No finite values for {metric}')
    edges = np.unique(np.quantile(x, np.linspace(0,1,bins+1)))
    if len(edges)<2: raise ValueError(f'Constant predictor: {metric}')
    return edges

def median_table(d, metric, edges, min_eddies=20, outcome='TiltDis'):
    """Shared quantile edges, daily median/IQR; keep gaps for suppressed bins."""
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=[metric,outcome]).copy()
    d['_bin']=pd.cut(d[metric],edges,include_lowest=True,labels=False)
    rows=[]
    for cyc in COLORS:
        for i in range(len(edges)-1):
            g=d.loc[d.Cyc.eq(cyc)&d._bin.eq(i)]
            enough = g.Eddy.nunique() >= min_eddies
            rows.append(dict(Cyc=cyc,bin=i,eddy_days=len(g),eddies=g.Eddy.nunique(),
                             x=g[metric].median(),median=g[outcome].median() if enough else np.nan,
                             q25=g[outcome].quantile(.25) if enough else np.nan,
                             q75=g[outcome].quantile(.75) if enough else np.nan,
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
        ax.plot(x,g['median'],color=COLORS[cyc],ls=linestyle,lw=1.7,ms=5,marker='.')
        if shade: ax.fill_between(x,g.q25,g.q75,color=COLORS[cyc],alpha=.10,lw=0)
    ax.set_xlabel(LABEL[metric]); ax.set_ylim(bottom=0)
    ax.grid(axis='y',alpha=.12,lw=.5)

def polarity_handles():
    return [Line2D([0],[0],color=color,lw=1.8,label=cyc) for cyc,color in COLORS.items()]

def selected_depth_figure(d,selected,bins=8,min_eddies=20,outcome='TiltDis'):
    style(); ncols=max(map(len,selected.values()))
    fig,axs=plt.subplots(2,ncols,figsize=(3.15*ncols,5.5),squeeze=False,sharey=True,sharex='row',layout='constrained')
    tables=[]
    for row,metric in enumerate(METRICS):
        # Identical bins across depths for a given property, and across polarity.
        edges=bin_edges(d.loc[d.Depth.isin(selected[metric])],metric,bins)
        for col,z in enumerate(selected[metric]):
            ax=axs[row,col]; part=d.loc[d.Depth.eq(z)]
            t=median_table(part,metric,edges,min_eddies,outcome=outcome)
            tables.append(t.assign(metric=metric,depth_m=z))
            curve(ax,t,metric,shade=True)
            ax.set_title(f'({chr(97+row*ncols+col)}) {z:.0f} m',loc='left')
        for ax in axs[row,len(selected[metric]):]:ax.set_visible(False)
        axs[row,0].set_ylabel('Tilt distance / surface $R_c$' if outcome == 'tilt_over_surface_Rc' else 'Tilt distance (km)')
    axs[0,0].legend(handles=polarity_handles(),frameon=False)
    if outcome == 'tilt_over_surface_Rc':
        for ax in fig.axes:
            ax.autoscale(enable=True, axis='y')
            ax.set_ylim(bottom=0)
    return fig,pd.concat(tables,ignore_index=True)

def environment_panel(ax, table, metric, linear_fit=True):
    """Draw saved environmental summaries; latitude ticks/equations use degrees south.

    Input bin coordinates and numerical fit tables retain signed degrees north.
    """
    curve(ax,table,metric,shade=True)
    if metric == 'lat':
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f'{-value:g}'))
        ax.set_xlabel('Latitude (°S)')
    if linear_fit:
        for cyc,color in [('AE','darkred'),('CE','navy')]:
            q=table[table.Cyc.eq(cyc)].dropna(subset=['x','median'])
            p=np.polyfit(q.x,q['median'],1);yhat=np.polyval(p,q.x)
            r2=1-np.sum((q['median']-yhat)**2)/np.sum((q['median']-q['median'].mean())**2)
            def fmt(v):
                if v==0:return '0'
                e=int(np.floor(np.log10(abs(v))));a=v/10**e
                return fr'{a:.2f}\times10^{{{e}}}' if abs(e)>=3 else f'{v:.2f}'
            ax.plot(q.x,yhat,'--',color=color,lw=1.5,
                    label=fr'$y={fmt(-p[0] if metric == "lat" else p[0])}x{p[1]:+.2f}$, $R^2={r2:.2f}$')
        ax.legend(frameon=False,loc='upper left')

def environment_figure2(d,axs=None,bins=8,min_eddies=20,linear_fit=True,outcome='TiltDis'):
    style();tables=[]
    if axs is None:
        fig,axs=plt.subplots(1,3,figsize=(9,3.35),sharey=True,layout='constrained')
    else:
        fig=axs[0].figure

    for i,(ax,metric) in enumerate(zip(axs,['lat','beta','log10_PV_grad_mag'])):
        t=median_table(d,metric,bin_edges(d,metric,bins),min_eddies,outcome=outcome)
        environment_panel(ax,t,metric,linear_fit=linear_fit)
        tables.append(t.assign(metric=metric))

    axs[0].set_ylabel('Tilt distance / surface $R_c$' if outcome=='tilt_over_surface_Rc' else 'Tilt distance (km)')

    if outcome=='tilt_over_surface_Rc':
        for ax in axs:
            ax.autoscale(enable=True,axis='y');ax.set_ylim(bottom=0)

    return fig,axs,pd.concat(tables,ignore_index=True)

def with_normalised_tilt(d, surface):
    """Normalise each Eddy-Day by its surface Rc (both distances in km)."""
    check_keys(surface)
    radius = surface[KEYS+['Rc']].rename(columns={'Rc':'surface_Rc'})
    out = d.drop(columns=['surface_Rc','tilt_over_surface_Rc'], errors='ignore').merge(
        radius, on=KEYS, how='left', validate='many_to_one')
    valid = np.isfinite(out.surface_Rc) & out.surface_Rc.gt(0)
    out['tilt_over_surface_Rc'] = out.TiltDis / out.surface_Rc.where(valid)
    return out
