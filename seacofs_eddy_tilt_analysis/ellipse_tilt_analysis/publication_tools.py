"""Publication panels; all estimates come from supplied data, never PDF results."""
from __future__ import annotations

import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import ellipse_tilt_tools as et

STYLE = {'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 10,
         'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
         'axes.spines.top': False, 'axes.spines.right': False,
         'axes.linewidth': .7, 'lines.linewidth': 1.5,
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white'}


def eligible(frame, min_days=5):
    """Require repeated valid observations in the actual plotted sample/class."""
    return frame.loc[frame.groupby('Eddy').Eddy.transform('size').ge(min_days)].copy()


def cluster_interval(sums, counts, *, equal_eddy=True, n_boot=1000, seed=42,
                     min_eddies=20):
    """Ratio estimator and whole-eddy percentile CI, scalar or vector valued.

    sums[i] is the scalar/vector sum for eddy i; counts[i] its observation
    count. Equal weighting averages eddy means. Observation weighting pools
    eddy sums, but still resamples entire eddies for the confidence interval.
    """
    sums = np.asarray(sums, float)
    if sums.ndim == 1:
        sums = sums[:, None]
    counts = np.asarray(counts, float)
    if len(sums) != len(counts) or (counts <= 0).any():
        raise ValueError('One positive count is required per eddy')
    if n_boot < 2:
        raise ValueError('n_boot must be at least two')
    if len(sums) < min_eddies:
        blank = np.full(sums.shape[1], np.nan)
        return blank, blank.copy(), blank.copy()
    if equal_eddy:
        sums = sums / counts[:, None]
        counts = np.ones(len(counts))
    estimate = sums.sum(axis=0) / counts.sum()
    rng = np.random.default_rng(seed)
    boot = np.empty((n_boot, sums.shape[1]))
    for i in range(n_boot):
        indices = rng.integers(len(counts), size=len(counts))
        boot[i] = sums[indices].sum(axis=0) / counts[indices].sum()
    low, high = np.percentile(boot, [2.5, 97.5], axis=0)
    return estimate, low, high


def alignment_estimate(frame, **kwargs):
    grouped = frame.groupby('Eddy').AlignmentCos2.agg(['sum', 'count'])
    result = cluster_interval(grouped['sum'], grouped['count'], **kwargs)
    return dict(estimate=result[0][0], low=result[1][0], high=result[2][0],
                eddies=len(grouped), observations=len(frame))


def histogram_estimate(frame, edges, **kwargs):
    if frame.empty:
        hist = np.empty((0, len(edges)-1))
    else:
        eddy_codes, eddies = pd.factorize(frame.Eddy)
        bins = np.searchsorted(edges, frame.AlignmentDeg, side='right') - 1
        bins = np.minimum(bins, len(edges)-2)  # include exactly 90 degrees
        if (bins < 0).any() or (frame.AlignmentDeg > edges[-1]).any():
            raise ValueError('Alignment values outside histogram bounds')
        hist = np.zeros((len(eddies), len(edges)-1))
        np.add.at(hist, (eddy_codes, bins), 1)
    est, low, high = cluster_interval(hist, hist.sum(axis=1), **kwargs)
    return pd.DataFrame({'angle': (edges[:-1]+edges[1:])/2, 'estimate': est*100,
                         'low': low*100, 'high': high*100})


def direction_data(all_shapes, *, depths=(0, 50, 100, 200, 300), min_ar=1.1,
                   min_tilt=5, max_ar=5, min_days=5,
                   ar_edges=(1.1, 1.3, 1.6, 2, 3, 5), **bootstrap):
    """Separate surface, AR-class and complete-depth samples, with audit tables."""
    depths = tuple(depths)
    if 0 not in depths or len(depths) != len(set(depths)):
        raise ValueError('depths must be unique and include the surface (0)')
    edges = np.asarray(ar_edges, float)
    if np.any(np.diff(edges) <= 0) or edges[0] > min_ar or edges[-1] < max_ar:
        raise ValueError('AR edges must increase and cover the admitted AR range')
    selected = et.select_rows(all_shapes, directional=True, min_ar=min_ar,
                              min_tilt=min_tilt, max_ar=max_ar)
    selected = selected.loc[selected.TiltDis.gt(0)]
    surface = eligible(selected.loc[selected.ShapeDepth.eq(0)], min_days)
    matched = et.matched_depth_sample(selected, depths)
    # Count days once; matched sample has exactly one row per depth on each day.
    day_counts = matched.loc[matched.ShapeDepth.eq(0)].groupby('Eddy').size()
    matched = matched.loc[matched.Eddy.isin(day_counts[day_counts >= min_days].index)]
    histograms, ar_rows, depth_rows, delta_rows, coverage = {}, [], [], [], []
    for cyc in et.COLOURS:
        g = surface.loc[surface.Cyc.eq(cyc)]
        histograms[cyc] = histogram_estimate(g, np.arange(0, 91, 5), **bootstrap)
        coverage.append(dict(Cyc=cyc, sample='surface', eddies=g.Eddy.nunique(), observations=len(g)))
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            mask = g.AxisRatio.ge(lo) & (g.AxisRatio.le(hi) if i == len(edges)-2 else g.AxisRatio.lt(hi))
            part = eligible(g.loc[mask], min_days)
            ar_rows.append(dict(Cyc=cyc, lower=lo, upper=hi, x=(lo+hi)/2,
                                **alignment_estimate(part, **bootstrap)))
        gm = matched.loc[matched.Cyc.eq(cyc)]
        for depth in depths:
            part = gm.loc[gm.ShapeDepth.eq(depth)]
            depth_rows.append(dict(Cyc=cyc, depth=depth, **alignment_estimate(part, **bootstrap)))
        coverage.append(dict(Cyc=cyc, sample='all-depth matched', eddies=gm.Eddy.nunique(),
                             observations=len(gm)//len(depths)))
        # Paired depth-minus-surface contrasts on exactly the same eddy-days.
        baseline = gm.loc[gm.ShapeDepth.eq(0), et.KEY+['AlignmentCos2']]
        for depth in depths:
            if depth == 0:
                continue
            paired = gm.loc[gm.ShapeDepth.eq(depth), et.KEY+['AlignmentCos2']].merge(
                baseline, on=et.KEY, suffixes=('', '_surface'), validate='one_to_one')
            paired['AlignmentCos2'] -= paired.AlignmentCos2_surface
            delta_rows.append(dict(Cyc=cyc, depth=depth, **alignment_estimate(paired, **bootstrap)))
    return dict(histograms=histograms, ar=pd.DataFrame(ar_rows), depth=pd.DataFrame(depth_rows),
                depth_contrasts=pd.DataFrame(delta_rows), coverage=pd.DataFrame(coverage),
                settings=dict(depths=depths, min_ar=min_ar, min_tilt=min_tilt, max_ar=max_ar,
                              min_days=min_days, ar_edges=list(ar_edges), **bootstrap))


def _panel(ax, letter):
    ax.text(-.17, 1.04, letter, transform=ax.transAxes, weight='bold', fontsize=11)
    ax.tick_params(length=3, width=.6)
    ax.grid(axis='y', color='.92', linewidth=.5)
    ax.set_axisbelow(True)


def _line_interval(ax, x, g, colour):
    ax.fill_between(x, g.low.to_numpy(float), g.high.to_numpy(float), color=colour, alpha=.16, lw=0)
    ax.plot(x, g.estimate, '-o', color=colour, ms=4)
    ax.vlines(x, g.low, g.high, color=colour, lw=.8)


def plot_direction(data):
    """Six panels: rows=polarity, columns=surface distribution, AR, depth."""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), layout='constrained')
        score_rows = pd.concat([data['ar'], data['depth']])
        finite = score_rows[['low','high']].to_numpy(float)
        finite = finite[np.isfinite(finite)]
        ymin, ymax = (min(-.03, finite.min()-.04), max(.4, finite.max()+.05)) if len(finite) else (-.05, .6)
        histmax = max([12] + [h.high.max() for h in data['histograms'].values() if h.high.notna().any()]) * 1.12
        for row, (cyc, colour) in enumerate(et.COLOURS.items()):
            ax = axes[row, 0]
            h = data['histograms'][cyc]
            edges = np.arange(0, 91, 5)
            ax.stairs(h.estimate, edges, fill=True, color=colour, alpha=.24, linewidth=0)
            ax.plot(h.angle, h.estimate, color=colour)
            ax.fill_between(h.angle, h.low, h.high, color=colour, alpha=.2, lw=0)
            ax.axhline(100/18, ls='--', color='.4', lw=.8)
            n = data['coverage'].query("Cyc == @cyc and sample == 'surface'").eddies.iloc[0]
            ax.text(.96, .94, f'{cyc} · n = {n:,}', ha='right', va='top', transform=ax.transAxes, color=colour)
            ax.set(xlim=(0,90), ylim=(0,histmax), xticks=[0,30,60,90],
                   xlabel='Tilt–major-axis angle (°)', ylabel='Probability per 5° bin (%)')
            # Column titles intentionally hidden.
            g = data['ar'].loc[data['ar'].Cyc.eq(cyc)]
            ax = axes[row, 1]
            # Equally spaced categorical classes avoid implying a fitted continuous curve.
            x = np.arange(len(g))
            _line_interval(ax, x, g, colour)
            ax.set(xticks=x, xticklabels=[f'{r.lower:g}–{r.upper:g}' for r in g.itertuples()],
                   ylim=(ymin,ymax), xlabel='Surface axis-ratio class', ylabel='Major-axis alignment score')
            ax.tick_params(axis='x', labelsize=6.5)
            ax.axhline(0, color='.4', ls='--', lw=.8)
            # Column titles intentionally hidden.
            for xi, r in zip(x, g.itertuples()):
                if np.isfinite(r.estimate):
                    ax.annotate(f'{r.eddies:,}', (xi,r.high), xytext=(0,5), textcoords='offset points',
                                ha='center', fontsize=6, color='.35')
            d = data['depth'].loc[data['depth'].Cyc.eq(cyc)]
            ax = axes[row, 2]
            _line_interval(ax, d.depth.to_numpy(), d, colour)
            ax.set(xticks=data['settings']['depths'], ylim=(ymin,ymax),
                   xlabel='Ellipse depth (m)', ylabel='Major-axis alignment score')
            ax.axhline(0, color='.4', ls='--', lw=.8)
            n = data['coverage'].query("Cyc == @cyc and sample == 'all-depth matched'").eddies.iloc[0]
            label_y = .12 if d.estimate.mean() > (ymin+ymax)/2 else .94
            ax.text(.96,label_y,f'Matched n = {n:,}', transform=ax.transAxes,
                    va='top', ha='right', fontsize=7,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=1))
            # Column titles intentionally hidden.
        for letter, ax in zip(string.ascii_lowercase, axes.flat):
            _panel(ax, letter)
        return fig


def magnitude_data(all_shapes, *, max_ar=5, min_days=5, n_boot=1000, seed=42, min_eddies=20):
    g = et.select_rows(all_shapes.loc[all_shapes.ShapeDepth.eq(0)], max_ar=max_ar)
    g = eligible(g, min_days)
    med = g.groupby(['Cyc','Eddy'])[['AxisRatio','TiltDis']].median().reset_index()
    stats = et.magnitude_summary(g, n_boot=n_boot, seed=seed, min_days=min_days, min_eddies=min_eddies)
    curves = []
    for cyc in et.COLOURS:
        part = med.loc[med.Cyc.eq(cyc)].copy()
        if len(part) < min_eddies or part.AxisRatio.nunique() < 2:
            continue
        part['bin'] = pd.qcut(part.AxisRatio, 6, duplicates='drop')
        for label, b in part.groupby('bin', observed=True):
            est, low, high = et.bootstrap_stat(b[['TiltDis']].to_numpy(), np.median,
                                               n_boot=n_boot, seed=seed, min_eddies=min_eddies)
            curves.append(dict(Cyc=cyc, x=b.AxisRatio.median(), estimate=est, low=low,
                               high=high, eddies=len(b), lower=label.left, upper=label.right))
    return dict(medians=med, statistics=stats,
                curves=pd.DataFrame(curves, columns=['Cyc','x','estimate','low','high','eddies','lower','upper']))


def plot_magnitude(data):
    """No directional AR/tilt thresholds; full extent retained in scatter panels."""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7), layout='constrained')
        med = data['medians']
        xmax = max(2.5, med.AxisRatio.max()) if len(med) else 2.5
        ymax = max(10, med.TiltDis.max()*1.05) if len(med) else 10
        for ax, (cyc, colour) in zip(axes[:2], et.COLOURS.items()):
            part = med.loc[med.Cyc.eq(cyc)]
            ax.scatter(part.AxisRatio, part.TiltDis, s=5, alpha=.13, color=colour,
                       rasterized=True, linewidths=0)
            b = data['curves'].loc[data['curves'].Cyc.eq(cyc)]
            _line_interval(ax, b.x.to_numpy(), b, colour)
            stat = data['statistics'].loc[data['statistics'].Cyc.eq(cyc) &
                                           data['statistics'].metric.eq('Between-eddy Spearman')]
            rho = stat.estimate.iloc[0] if len(stat) else np.nan
            ax.set(xlim=(1,xmax*1.03), ylim=(0,ymax), xlabel='Median surface axis ratio',
                   ylabel='Median tilt distance (km)', title=f'{cyc} · n = {len(part):,}')
            ax.text(.96,.94, f'Spearman ρ = {rho:.2f}', ha='right', va='top', transform=ax.transAxes, fontsize=7,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=.85, pad=1))
            ax.xaxis.set_major_locator(MaxNLocator(4))
        ax=axes[2]
        for i, (cyc, colour) in enumerate(et.COLOURS.items()):
            stat = data['statistics'].loc[data['statistics'].Cyc.eq(cyc) &
                                           data['statistics'].metric.str.startswith('Within')]
            if len(stat):
                r=stat.iloc[0]
                ax.hlines(i,r.low,r.high,color=colour,lw=2)
                ax.plot(r.estimate,i,'o',color=colour,ms=5)
        ax.axvline(0,color='.4',ls='--',lw=.8)
        ax.set(yticks=[0,1],yticklabels=['AE','CE'],ylim=(1.6,-.6),
               xlabel='Within-eddy slope\n(km per unit axis ratio)', title='Within-eddy association')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        for letter, ax in zip(string.ascii_lowercase,axes):
            _panel(ax,letter)
        return fig


def ar_histogram_data(all_shapes, *, regions=None, min_ar=1.1, min_tilt=5,
                      max_ar=5, ar_edges=(1.1, 1.3, 1.6, 2, np.inf),
                      min_days=5, **bootstrap):
    """Surface-only 5-degree distributions, equal weight within each AR class.

    Region, direction QC and AR-class selection precede minimum-day filtering.
    All classes use the same angle bins. An infinite final AR edge still
    respects max_ar. Sparse classes remain in the audit with NaN estimates.
    """
    edges = np.asarray(ar_edges, float)
    if (len(edges) < 2 or not np.isfinite(edges[:-1]).all()
            or np.isnan(edges[-1]) or np.any(np.diff(edges) <= 0)
            or edges[0] > min_ar or edges[-1] < max_ar):
        raise ValueError('Increasing AR edges must cover the admitted AR range')
    surface = all_shapes.loc[all_shapes.ShapeDepth.eq(0)]
    if regions is not None:
        surface = surface.loc[surface.Region.isin(regions)]
    selected = et.select_rows(surface, directional=True, min_ar=min_ar,
                               min_tilt=min_tilt, max_ar=max_ar)
    selected = selected.loc[selected.TiltDis.gt(0)]
    tables, coverage = [], []
    for cyc in et.COLOURS:
        group = selected.loc[selected.Cyc.eq(cyc)]
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            use = group.AxisRatio.ge(lo) & (group.AxisRatio.le(hi) if i == len(edges)-2 else group.AxisRatio.lt(hi))
            part = eligible(group.loc[use], min_days)
            label = f'≥{lo:g}' if np.isinf(hi) else f'{lo:g}–{hi:g}'
            h = histogram_estimate(part, np.arange(0, 91, 5), **bootstrap)
            h = h.assign(Cyc=cyc, ar_class=i, label=label, lower=lo, upper=hi,
                         eddies=part.Eddy.nunique(), observations=len(part))
            tables.append(h)
            coverage.append(dict(Cyc=cyc, ar_class=i, label=label,
                                 eddies=part.Eddy.nunique(), observations=len(part),
                                 plotted=h.estimate.notna().any()))
    return dict(histograms=pd.concat(tables, ignore_index=True), coverage=pd.DataFrame(coverage),
                settings=dict(min_ar=min_ar, min_tilt=min_tilt, max_ar=max_ar,
                              ar_edges=list(ar_edges), min_days=min_days, regions=regions,
                              **bootstrap))


def plot_ar_histograms(data, *, title='All regions', show_ci=True, ymax=None):
    """One row, AE red shades / CE blue shades; darker means greater AR."""
    table = data['histograms']
    n_classes = len(data['settings']['ar_edges'])-1
    if ymax is None:
        peak = table['high' if show_ci else 'estimate'].max()
        ymax = max(10., float(peak)*1.18) if np.isfinite(peak) else 10.
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharex=True, sharey=True,
                                 layout='constrained')
        for ax, cyc, cmap in zip(axes, ['AE', 'CE'], ['Reds', 'Blues']):
            colours = plt.get_cmap(cmap)(np.linspace(.48, .9, n_classes))
            for i in range(n_classes):
                h = table.loc[table.Cyc.eq(cyc) & table.ar_class.eq(i)]
                if h.empty or not h.estimate.notna().any():
                    continue
                label = f'{h.label.iloc[0]}  (n={int(h.eddies.iloc[0]):,})'
                ax.plot(h.angle, h.estimate, color=colours[i], lw=1.6, label=label)
                if show_ci:
                    ax.fill_between(h.angle, h.low, h.high, color=colours[i], alpha=.1, lw=0)
            ax.axhline(100/18, color='.45', ls='--', lw=.8, zorder=0)
            ax.set(xlim=(0,90), ylim=(0,ymax), xticks=[0,15,30,45,60,75,90],
                   xlabel='Tilt–major-axis angle (°)', title=cyc)
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(handles, labels, title='Axis ratio (eddy count)', frameon=False,
                          loc='upper right', fontsize=7, title_fontsize=7)
            else:
                ax.text(.5,.5,'Insufficient sample',ha='center',transform=ax.transAxes)
        axes[0].set_ylabel('Probability per 5° bin (%)')
        for letter, ax in zip('ab', axes):
            _panel(ax, letter)
        fig.suptitle(title, fontsize=10)
        return fig
