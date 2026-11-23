"""Tagged manuscript figure styles, adapted locally for agreed PV selections."""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter
from scipy.stats import binned_statistic_2d
import census_plots as tilt

def bin_edges_fd(x, xgrid, rule='fd'): # Freedman-Diaconis (fg) rationale
    n = len(x)
    if n < 2: return np.array([np.min(x), np.max(x)])
    rng = np.ptp(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    std = np.std(x, ddof=1)

    # raw width (km)
    if rule.lower() == 'fd':
        h = 2 * (iqr if iqr > 0 else 1.349*std) / (n ** (1/3))
    else:  # 'scott'
        h = 3.5 * std / (n ** (1/3))

    # fallback if degenerate
    if not np.isfinite(h) or h <= 0:
        h = rng / max(10, np.sqrt(n))

    # snap to grid spacing
    base = _grid_step(xgrid)
    h = _nice_step(h, base)

    lo = np.floor(np.min(x) / h) * h
    hi = np.ceil(np.max(x) / h) * h
    return np.arange(lo, hi + h, h)

def _grid_step(G):
    gx = np.diff(np.sort(np.unique(G.ravel())))
    return np.nanmedian(gx[gx > 0])

def _nice_step(h, base):
    s = h / base
    for k in [1, 2, 2.5, 5, 10]:
        if s <= k: return k * base
    return np.ceil(s) * base

def plot_map(data, grid, settings):
    df_eddies = data

    cmap=LinearSegmentedColormap.from_list("orange_green",["forestgreen","white","darkorange"])
    d=df_eddies.dropna(subset=['topo_plan_ratio','xc','yc'])
    xbins=bin_edges_fd(d.xc.values,grid.X_grid,rule='fd')
    ybins=bin_edges_fd(d.yc.values,grid.Y_grid,rule='fd')
    xcen=(xbins[:-1]+xbins[1:])/2; ycen=(ybins[:-1]+ybins[1:])/2
    lim=settings['map_log_ratio_limit']; t=np.log(settings['dominance_factor'])

    fig=plt.figure(figsize=settings['figure8_size_inches'],constrained_layout=True)
    gs=fig.add_gridspec(1,3)
    ax0=fig.add_subplot(gs[0,0]); ax1=fig.add_subplot(gs[0,1],sharey=ax0); ax2=fig.add_subplot(gs[0,2])
    axs=[ax0,ax1,ax2]

    for ax,cyc in zip(axs[:2],['AE','CE']):
        dd=d[d.Cyc==cyc]
        H=binned_statistic_2d(dd.xc,dd.yc,dd.topo_plan_ratio,
                                   statistic='mean',bins=[xbins,ybins]).statistic.T
        ax.contourf(grid.X_grid,grid.Y_grid,np.where(grid.mask_rho==0,1,np.nan),
                    levels=[.5,1.5],colors=['r' if cyc=='AE' else 'b'],alpha=.25)
        ax.contour(grid.X_grid,grid.Y_grid,grid.mask_rho,levels=[.5],
                   colors='r' if cyc=='AE' else 'b',linewidths=2,zorder=15)
        ax.contour(grid.X_grid,grid.Y_grid,grid.h,levels=[4000],colors='k')
        m=ax.pcolormesh(xbins,ybins,H,cmap=cmap,vmin=-lim,vmax=lim,shading='auto')
        ax.contour(xcen,ycen,H,levels=[-t],linewidths=3,
                    linestyles='-', zorder=15, colors='darkolivegreen')
        ax.contour(xcen,ycen,H,levels=[t],linewidths=3,
                    linestyles='-', zorder=15, colors='darkgoldenrod')
        ax.text(100,900,cyc,fontweight='bold',color='r' if cyc=='AE' else 'b',fontsize=14)
        ax.set(xlabel='x (km)',facecolor='lightgrey',
               xlim=(grid.X_grid.min(),grid.X_grid.max()),
               ylim=(grid.Y_grid.min(),grid.Y_grid.max()),aspect='equal')

    axs[1].tick_params(axis='y',labelleft=False)
    ax0.set_ylabel('y (km)')
    cb=fig.colorbar(m,ax=axs[:2],location='right',shrink=.8)
    cb.set_label(r'Mean $\log(|\beta_t|/|\beta_p|)$')

    h=np.where(grid.mask_rho,grid.h/1e3,np.nan)
    m2=ax2.contourf(grid.X_grid,grid.Y_grid,h,levels=10,cmap='Greys_r',zorder=2)
    ax2.contour(grid.X_grid,grid.Y_grid,grid.h,levels=[4000],colors='k',zorder=3)
    ax2.set(xlabel='x (km)',ylabel='y (km)',
            xlim=(grid.X_grid.min(),grid.X_grid.max()),
            ylim=(grid.Y_grid.min(),grid.Y_grid.max()),aspect='equal')
    ax2.grid(True); ax2=tilt.lat_lon_contours(ax2,grid)
    cb2=fig.colorbar(m2,ax=ax2,location='right',shrink=.8)
    cb2.set_label('Bathymetry (km)')

    for ax,l in zip(axs,['a)','b)','c)']):
        ax.text(-.1,1.01,l,transform=ax.transAxes,ha='left',va='top',
                fontsize=12,fontweight='bold')

    return fig


def plot_distribution(data, grid, settings):
    surface = data[data.regime.ne('unclassified')].copy()
    DOMINANCE_FACTOR=settings['dominance_factor']
    fig,(ax,ax2,ax3)=plt.subplots(1,3,figsize=settings['figure9_size_inches'],gridspec_kw={'width_ratios':[2,.5,.7]})

    cols={'AE':'red','CE':'blue'}
    t=np.log(DOMINANCE_FACTOR); v=surface.topo_plan_ratio
    bw=t/3
    bins=np.arange(np.floor(v.min()/bw)*bw,np.ceil(v.max()/bw)*bw+bw,bw)

    for cyc in ['AE','CE']:
        x=surface.loc[surface.Cyc.eq(cyc),'topo_plan_ratio'].dropna()
        ax.hist(x,bins=bins,weights=np.ones(len(x))*100/len(x),color=cols[cyc],alpha=.5,label=cyc)

    ax.axvspan(v.min(),-t,color='lightgreen',alpha=.25,zorder=0)
    ax.axvspan(-t,t,color='lightgrey',alpha=.35,zorder=0)
    ax.axvspan(t,v.max(),color='moccasin',alpha=.4,zorder=0)
    ax.axvline(-t,color='k',ls='--',lw=1); ax.axvline(t,color='k',ls='--',lw=1)
    ax.set_xlim(v.min(),v.max())

    ticks=list(ax.get_xticks())+[-t,t]
    labels=[f'{x:g}' if abs(x)>t else '' for x in ax.get_xticks()]+[
        rf'$-\log({DOMINANCE_FACTOR:g})$',rf'$\log({DOMINANCE_FACTOR:g})$']
    ax.set_xticks(ticks,labels)

    xmin,xmax=ax.get_xlim(); trans=ax.get_xaxis_transform()
    ax.text((xmin-t)/2,.96,'Planetary',transform=trans,ha='center',va='top')
    ax.text(0,.96,'Mixed',transform=trans,ha='center',va='top')
    ax.text((t+xmax)/2,.96,'Topographic',transform=trans,ha='center',va='top')
    ax.set(xlabel=r'Log ratio between topographic and planetary beta magnitude' '\n'
                  r'$\log(|\beta_t|/|\beta_p|)$',ylabel='Population (%)',
           xlim=(v.min(),v.max()))
    ax.legend(frameon=False,loc='upper right')

    ax1 = ax

    doms=['plan','mixed','topo']; xx=np.arange(3); w=.35
    for i,cyc in enumerate(['AE','CE']):
        d=surface[surface.Cyc.eq(cyc)]
        p=100*d.dom.value_counts(normalize=True).reindex(doms,fill_value=0)
        ax2.bar(xx+(i-.5)*w,p,width=w,color=cols[cyc],alpha=.7,label=cyc)
    ax2.set(xticks=xx,xticklabels=['Plan.','Mixed','Topo.'],ylabel='Population (%)')
    ax2.legend(frameon=False)

    vars=['PV_grad_plan_mag','PV_grad_topo_mag','PV_grad_mag']
    centres=np.arange(3); off=.2
    pos=np.ravel(np.column_stack((centres-off,centres+off)))
    data=[surface.loc[surface.Cyc.eq(c),var].dropna() for var in vars for c in ['AE','CE']]

    bp=ax3.boxplot(data,positions=pos,widths=.32,patch_artist=True,showfliers=False)
    for patch,c in zip(bp['boxes'],['AE','CE']*3):
        patch.set_facecolor(cols[c]); patch.set_alpha(.5)
    for med in bp['medians']: med.set_color('k')

    ax3.set_xticks(centres,[r'$|\beta_p|$',r'$|\beta_t|$',r'$|\nabla PV|$'])
    ax3.set_ylabel('PV gradient magnitude')
    ax3.set_yscale('log')

    for ax,l in zip((ax1,ax2,ax3),['a)','b)','c)']):
        ax.set_title(l,loc='left',fontsize=12,fontweight='bold')

    plt.tight_layout()

    return fig


def plot_planetary(data, grid, settings):
    planetary = data
    MIN_DEPTH_m=settings['context_isobath_m']
    fig = plt.figure(figsize=settings['figure10_size_inches'], constrained_layout=True)

    axs = [
        fig.add_subplot(1, 4, 1),
        fig.add_subplot(1, 4, 2, projection='polar'),
        fig.add_subplot(1, 4, 3),
        fig.add_subplot(1, 4, 4, projection='polar')
    ]

    mag_bins = [settings['direction_min_tilt_km'], 10, 20, 30, 40, np.inf]

    for i, cyc in enumerate(['AE', 'CE']):
        ax, axw = axs[2*i], axs[2*i+1]
        df = planetary[planetary.Cyc == cyc].copy()
        cmap = 'Reds' if cyc == 'AE' else 'Blues'
        col = 'r' if cyc == 'AE' else 'dodgerblue'
        colors = getattr(plt.cm, cmap)(np.linspace(.15, 1, len(mag_bins)-1))

        bath = ax.contourf(
            grid.lon_rho, grid.lat_rho,
            np.where(grid.mask_rho, grid.h/1e3, np.nan),
            cmap='Greys_r'
        )

        ax.hist2d(df.lon, df.lat, bins=50, cmap=cmap, alpha=.6, cmin=2)
        ax.scatter(df.lon, df.lat, s=3, c=col, alpha=.5, edgecolors='none')
        ax.contour(grid.lon_rho, grid.lat_rho, grid.h,
                   levels=[MIN_DEPTH_m], colors='yellow')

        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
        gap = .3
        ax.set_xlim(planetary.lon.min()-gap, planetary.lon.max()+gap)
        ax.set_ylim(planetary.lat.min()-gap, planetary.lat.max()+gap)

        ax.set_aspect('equal', adjustable='box')

        df = df.dropna(subset=['TiltDir', 'TiltDis'])
        tilt.plot_windrose(axw, df, title='', mag_bins=mag_bins, colors=colors, step=settings['planetary_radial_tick_step'], rlim=None)#, rotate=20)

        counts = tilt.windrose_counts(df.TiltDir, df.TiltDis, mag_bins=mag_bins)[0]
        peak = counts.sum(axis=0).max()
        minimum = settings['planetary_min_radial_limits'][cyc]
        radial_step = settings['planetary_radial_tick_step']
        limit = minimum if peak <= minimum else np.ceil(peak/radial_step)*radial_step
        axw.set_yticks(np.arange(radial_step,limit+1e-9,radial_step))
        axw.set_ylim(0,limit)

        theta = np.deg2rad(
            planetary.loc[planetary.Cyc == cyc, 'PV_grad_theta'].dropna()
        )
        theta_mean = np.arctan2(np.mean(np.sin(theta)), np.mean(np.cos(theta)))
        rmax = axw.get_ylim()[1]
        axw.annotate(
            '',
            xy=(theta_mean, 0.9*rmax),
            xytext=(theta_mean, 0),
            arrowprops=dict(arrowstyle='->', lw=2.5, color='limegreen')
        )

        theta = np.deg2rad(
            planetary.loc[planetary.Cyc == cyc, 'TiltDir'].dropna()
        )
        theta_mean = np.arctan2(np.mean(np.sin(theta)), np.mean(np.cos(theta)))
        axw.annotate(
            '',
            xy=(theta_mean, 0.9*rmax),
            xytext=(theta_mean, 0),
            arrowprops=dict(arrowstyle='->', lw=2.5, color='magenta')
        )
    
        axw.legend(
            title=f'{cyc}\ntilt dist. (km)',
            loc='upper left',
            bbox_to_anchor=(1.0, 1.25),
            frameon=False
        )

    cbar = fig.colorbar(bath, ax=axs, orientation='vertical',
                        fraction=.02, pad=.01)
    cbar.set_label('Depth (km)')

    for ax, label in zip(axs, ['a)', 'b)', 'c)', 'd)']):
        ax.text(-0.1, 1.01, label,
                transform=ax.transAxes,
                ha='left', va='top',
                fontsize=12, fontweight='bold')


    return fig


def _relative_mean(values):
    theta = np.deg2rad(np.asarray(values))
    vector = np.mean(np.exp(1j*theta))
    return np.rad2deg(np.angle(vector)) % 360, abs(vector)


def plot_topographic(data, grid, settings):
    """Maps and signed PV-relative tilt roses for shallow topographic eddy-days.

    Zero follows each day's PV gradient, clockwise is positive. Roses show
    percentages within each polarity, stacked by tilt distance; magenta arrows
    are circular mean relative bearings, with a fixed display length.
    """
    sample = data.copy()
    sample['relative_bearing_deg'] = (sample.TiltDir-sample.PV_grad_theta+180)%360-180
    edges = [settings['direction_min_tilt_km'],10,20,30,40,np.inf]
    width = settings['topographic_angular_bin_width_deg']
    counts = {}
    for cyc in ['AE','CE']:
        d = sample[sample.Cyc.eq(cyc)]
        counts[cyc] = tilt.windrose_counts(d.relative_bearing_deg,d.TiltDis,
            mag_bins=edges,dir_bins=np.arange(0,360+width,width),dir_shift=width/2)[0]
        if counts[cyc].sum() != len(d):
            raise ValueError('Topographic rose bins do not cover all selected eddy-days')
    step = settings['topographic_radial_tick_step_pct']
    limit = max(2*step,step*np.ceil(max(c.sum(axis=0).max()*100/c.sum()
                                     for c in counts.values())/step))
    fig = plt.figure(figsize=settings['figure12_size_inches'],constrained_layout=True)
    axes = [fig.add_subplot(1,4,k+1,projection='polar' if k%2 else None) for k in range(4)]
    for i,cyc in enumerate(['AE','CE']):
        ax,rose = axes[2*i:2*i+2]
        d = sample[sample.Cyc.eq(cyc)]
        cmap = 'Reds' if cyc=='AE' else 'Blues'
        colours = plt.get_cmap(cmap)(np.linspace(.15,1,len(edges)-1))
        bath = ax.contourf(grid.X_grid,grid.Y_grid,
            np.where(grid.mask_rho,grid.h/1000,np.nan),cmap='Greys_r')
        ax.scatter(d.xc,d.yc,s=3,c='r' if cyc=='AE' else 'dodgerblue',
                   alpha=.5,edgecolors='none',zorder=3)
        ax.contour(grid.X_grid,grid.Y_grid,grid.h,
                   levels=[settings['context_isobath_m']],colors='yellow',zorder=2)
        ax.set(xlabel='x (km)',ylabel='y (km)',
               xlim=(grid.X_grid.min(),grid.X_grid.max()),
               ylim=(grid.Y_grid.min(),grid.Y_grid.max()),aspect='equal')
        ax.set_axisbelow(True)
        ax.grid(True)
        # Exact overlay routine used in Figure 8c.
        tilt.lat_lon_contours(ax,grid)
        percent = counts[cyc]*100/len(d)
        assert np.isclose(percent.sum(),100)
        centres = np.deg2rad(np.arange(counts[cyc].shape[1])*width)
        bottom = np.zeros(len(centres))
        for j,(values,colour) in enumerate(zip(percent,colours)):
            hi = 'inf' if np.isinf(edges[j+1]) else f'{edges[j+1]:g}'
            rose.bar(centres,values,width=np.deg2rad(width),bottom=bottom,
                     color=colour,edgecolor=(0,0,0,.2),label=f'[{edges[j]:g}–{hi})')
            bottom += values
        rose.set_theta_zero_location('N');rose.set_theta_direction(-1)
        rose.set_xticks(np.deg2rad(np.arange(0,360,45)),
            ['0°\n'+r'$\nabla\mathrm{PV}$','45°','90°','135°','180°','−135°','−90°','−45°'])
        rose.set_yticks(np.arange(step,limit+.1,step));rose.set_ylim(0,limit)
        rose.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
        rose.set_rlabel_position(225 if cyc=='AE' else 315)
        mean,_ = _relative_mean(d.relative_bearing_deg)
        if np.isfinite(mean):
            rose.annotate('',xy=(np.deg2rad(mean),.9*limit),xytext=(np.deg2rad(mean),0),
                arrowprops=dict(arrowstyle='->',lw=2.5,color='magenta'))
        rose.legend(title=f'{cyc}\ntilt dist. (km)',loc='upper left',
                    bbox_to_anchor=(1.,1.25),frameon=False)
    colourbar = fig.colorbar(bath,ax=axes,orientation='vertical',fraction=.02,pad=.01)
    colourbar.set_label('Depth (km)')
    for ax,label in zip(axes,['a)','b)','c)','d)']):
        ax.text(-.1,1.01,label,transform=ax.transAxes,ha='left',va='top',
                fontsize=12,fontweight='bold')
    return fig
