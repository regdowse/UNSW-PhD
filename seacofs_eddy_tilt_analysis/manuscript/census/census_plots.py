"""Legacy census presentation, adapted locally; no legacy runtime imports.

Sources: census/tilt_census.ipynb cells 5,11,12 and seacofs_tilt_tools.py
at snapshot bb1d6c6. Population selection is supplied by census.py.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
from scipy.ndimage import gaussian_filter1d

AE_COLOUR = '#C44E52'
CE_COLOUR = '#2878A5'
TYPE_COLOURS = {'AE': AE_COLOUR, 'CE': CE_COLOUR}
SMOOTH_SIGMA_BINS = 1.15
LEVELS_LAT = [-40, -35, -30, -25]
LEVELS_LON = [150, 155, 160]
STYLE = {'figure.dpi':120, 'savefig.dpi':600, 'font.size':9.5,
         'axes.titlesize':10.5, 'axes.labelsize':9.5, 'legend.fontsize':8.5,
         'axes.spines.top':False, 'axes.spines.right':False, 'axes.linewidth':.8,
         'xtick.direction':'out', 'ytick.direction':'out', 'pdf.fonttype':42,
         'ps.fonttype':42}

def weighted_quantile(values,quantiles,weights=None):
    values=np.asarray(values,float); weights=np.ones(values.size) if weights is None else np.asarray(weights,float)
    ok=np.isfinite(values)&np.isfinite(weights)&(weights>0); values,weights=values[ok],weights[ok]
    if not values.size: return np.full(np.atleast_1d(quantiles).shape,np.nan)
    order=np.argsort(values); values,weights=values[order],weights[order]
    cdf=(np.cumsum(weights)-.5*weights)/weights.sum()
    return np.interp(np.atleast_1d(quantiles),cdf,values)

def gaussian_smooth(y,sigma=SMOOTH_SIGMA_BINS):
    radius=max(1,int(np.ceil(4*sigma))); x=np.arange(-radius,radius+1)
    kernel=np.exp(-.5*(x/sigma)**2); kernel/=kernel.sum()
    return np.convolve(y,kernel,mode='same')

def draw_vertical_weighted_box(ax,values,weights,x,colour,box_width):
    q05,q25,q50,q75,q95=weighted_quantile(values,[.05,.25,.5,.75,.95],weights)
    half_width=box_width/2
    ax.plot([x,x],[q05,q95],color=colour,lw=1.3,solid_capstyle='round',zorder=5)
    ax.plot([x-.55*half_width,x+.55*half_width],[q05,q05],color=colour,lw=1.1,zorder=5)
    ax.plot([x-.55*half_width,x+.55*half_width],[q95,q95],color=colour,lw=1.1,zorder=5)
    ax.add_patch(Rectangle((x-half_width,q25),box_width,q75-q25,facecolor='white',edgecolor=colour,alpha=.88,lw=1.3,zorder=6))
    ax.plot([x-half_width,x+half_width],[q50,q50],color=colour,lw=2.2,zorder=7)

def vertical_split_violin_box(ax,groups,column,title,ylabel,*,daily=False,bins=None,show_legend=False, smoothing_sigma=1.15):
    if bins is None:
        all_values=np.concatenate([g[column].dropna().to_numpy() for _,g in groups])
        bins=np.linspace(max(0,np.nanmin(all_values)),np.nanquantile(all_values,.995),38)
    centres=.5*(bins[:-1]+bins[1:]); curves={}; peak=0.0
    for group_index,(group_label,data) in enumerate(groups):
        for cyc in ['AE','CE']:
            g=data[data.Cyc==cyc].dropna(subset=[column]); values=g[column].to_numpy()
            weights=np.ones(len(g))
            hist,_=np.histogram(values,bins=bins,weights=weights)
            hist=hist/hist.sum()*100 if hist.sum() else hist.astype(float)
            smooth=gaussian_smooth(hist, sigma=smoothing_sigma); peak=max(peak,float(np.nanmax(smooth)))
            curves[(group_index,cyc)]=(smooth,values,weights)
    half_violin_width=.55
    for group_index,(group_label,data) in enumerate(groups):
        for cyc,sign in [('AE',-1),('CE',1)]:
            smooth,values,weights=curves[(group_index,cyc)]
            width=half_violin_width*smooth/peak if peak>0 else smooth
            edge=group_index+sign*width; colour=TYPE_COLOURS[cyc]
            ax.fill_betweenx(centres,group_index,edge,color=colour,alpha=.30,lw=0)
            ax.plot(edge,centres,color=colour,lw=1.8)
            draw_vertical_weighted_box(ax,values,weights,group_index+sign*.13,colour,.10)
    ax.set_xlim(-.55,len(groups)-.45); ax.set_xticks(np.arange(len(groups)),[label for label,_ in groups])
    ax.set(title=title,ylabel=ylabel); ax.grid(axis='y',alpha=.14,lw=.6)
    if show_legend:
        handles=[Rectangle((0,0),1,1,facecolor=TYPE_COLOURS[cyc],alpha=.45,label=cyc) for cyc in ['AE','CE']]
        ax.legend(handles=handles,frameon=False,ncol=2,loc='upper right')

def panel_label(ax,label):
    ax.text(-.10,1.02,label+')',transform=ax.transAxes,fontweight='bold',fontsize=11)

def mirrored_hist(ax, ae, ce, bins='fd', xlabel='', *, ylabel=None,
                  colors=('r', 'b'), alpha=.45, xlim=None,
                  normalize=False, smooth=True, sigma=1.2):
    """Mirrored AE/CE histogram with optional automatic bins and smoothing."""

    ae = np.asarray(ae, float)
    ce = np.asarray(ce, float)
    ae = ae[np.isfinite(ae)]
    ce = ce[np.isfinite(ce)]

    if xlim is not None:
        ae = ae[(ae >= xlim[0]) & (ae <= xlim[1])]
        ce = ce[(ce >= xlim[0]) & (ce <= xlim[1])]

    # Shared bin edges
    edges = np.histogram_bin_edges(np.r_[ae, ce], bins=bins, range=xlim)

    ae_counts, _ = np.histogram(ae, bins=edges)
    ce_counts, _ = np.histogram(ce, bins=edges)

    if normalize:
        ae_counts = ae_counts / ae_counts.sum()
        ce_counts = ce_counts / ce_counts.sum()

    centers = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    ax.bar(centers,  ae_counts, widths, color=colors[0], alpha=alpha,
           label='AE', edgecolor='none')
    ax.bar(centers, -ce_counts, widths, color=colors[1], alpha=alpha,
           label='CE', edgecolor='none')

    if smooth:
        ae_smooth = gaussian_filter1d(ae_counts.astype(float), sigma)
        ce_smooth = gaussian_filter1d(ce_counts.astype(float), sigma)

        ax.plot(centers,  ae_smooth, color=colors[0], lw=2)
        ax.plot(centers, -ce_smooth, color=colors[1], lw=2)

    ax.axhline(0, color='0.3', lw=.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel or ('Probability' if normalize else 'Frequency'))

    ax.spines[['top', 'right']].set_visible(False)

    ymax = max(ae_counts.max(), ce_counts.max()) * 1.08
    ax.set_ylim(-ymax, ymax)

    if xlim is not None:
        ax.set_xlim(xlim)

    return ax

def choose_dir_bins_cardinal(*dfs, col: str = "TiltDir", min_bins: int = 8, max_bins: int = 36, min_avg_per_sector: int = 8):
    """Choose circular direction bins while keeping cardinal directions aligned."""

    candidates = np.array([4, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 36, 45, 60, 72])
    candidates = candidates[(candidates >= min_bins) & (candidates <= max_bins)]
    candidates = candidates[candidates % 4 == 0]
    counts = []
    for df in dfs:
        x = df[col].to_numpy(float)
        counts.append(np.isfinite(x).sum())
    n = min(counts) if counts else 0
    if n == 0 or candidates.size == 0:
        k = 16
    else:
        k0 = int(np.ceil(2 * n ** (1 / 3)))
        k = candidates[np.argmin(np.abs(candidates - k0))]
        while k > candidates.min() and (n / k) < min_avg_per_sector:
            k = candidates[candidates < k].max()
    return np.linspace(0.0, 360.0, k + 1), 180.0 / k

def windrose_counts(directions_deg, magnitudes, *, mag_bins, dir_bins=None, dir_shift=None):
    """Count tilt directions and magnitudes for polar bar plots."""

    directions = np.asarray(directions_deg, float)
    magnitudes = np.asarray(magnitudes, float)
    ok = np.isfinite(directions) & np.isfinite(magnitudes)
    directions = directions[ok]
    magnitudes = magnitudes[ok]
    if directions.size == 0:
        return None
    if dir_bins is None or dir_shift is None:
        tmp = pd.DataFrame({"TiltDir": directions})
        dir_bins, dir_shift = choose_dir_bins_cardinal(tmp)
    k = len(dir_bins) - 1
    binw_deg = 360.0 / k
    angles = np.deg2rad(np.arange(k) * binw_deg)
    width = np.deg2rad(binw_deg)
    directions = np.mod(directions + dir_shift, 360.0)
    dir_idx = np.digitize(directions, dir_bins, right=False) - 1
    mag_idx = np.digitize(magnitudes, mag_bins, right=False) - 1
    counts = np.zeros((len(mag_bins) - 1, k), float)
    ok = (dir_idx >= 0) & (dir_idx < k) & (mag_idx >= 0) & (mag_idx < len(mag_bins) - 1)
    for d_i, m_i in zip(dir_idx[ok], mag_idx[ok]):
        counts[m_i, d_i] += 1
    return counts, angles, width

def plot_windrose(ax, df: pd.DataFrame, *, title: str = "", mag_bins=(0, 10, 20, 30, 40, np.inf),
                  colors=None, step=None, rlim=None, mag='TiltDis', theta='TiltDir'):
    """Draw one stacked tilt windrose."""

    if colors is None:
        cmap = "Reds" if df.Cyc.iloc[0] == "AE" else "Blues"
        colors = getattr(plt.cm, cmap)(np.linspace(0.15, 1, len(mag_bins) - 1))
    data = windrose_counts(df[theta], df[mag], mag_bins=mag_bins)
    if data is None:
        ax.set_axis_off()
        return ax
    counts, angles, width = data
    bottom = np.zeros(counts.shape[1])
    for i in range(counts.shape[0]):
        hi = "inf" if np.isinf(mag_bins[i + 1]) else f"{mag_bins[i + 1]:g}"
        ax.bar(angles, counts[i], width=width, bottom=bottom, color=colors[i], edgecolor=(0, 0, 0, 0.2), label=f"[{mag_bins[i]:g}-{hi})")
        bottom += counts[i]
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if step is not None:
        local_rmax = np.max(counts.sum(axis=0))
        local_rmax = 1 if local_rmax == 0 else local_rmax
        top = int(np.ceil(local_rmax / step) * step)
        rticks = np.arange(step, top + step, step)
        ax.set_rticks(rticks)
    if rlim is not None:
        ax.set_rlim(0, rlim)
    if np.any(df.Cyc.isin(["AE", "CE"])):
        ax.set_rlabel_position(225 if df.Cyc.iloc[0] == "AE" else 315)
    ax.set_title(title)
    return ax

def lat_lon_contours(ax, grid,
                     levels_lat=LEVELS_LAT,
                     levels_lon=LEVELS_LON,
                     land_mask=False, land_clr='k'):
    # Latitude contours
    c1 = ax.contour(
        grid.X_grid,
        grid.Y_grid,
        grid.lat_rho,
        levels=levels_lat,
        colors='k',
        linewidths=0.5
    )
    ax.clabel(
        c1,
        fmt=lambda v: f"{-v:.0f}°S",
        inline=True,
        colors='k'
    )

    # Longitude contours
    c2 = ax.contour(
        grid.X_grid,
        grid.Y_grid,
        grid.lon_rho,
        levels=levels_lon,
        colors='k',
        linewidths=0.5
    )
    ax.clabel(
        c2,
        fmt=lambda v: f"{v:.0f}°E",
        inline=True,
        colors='k'
    )

    if land_mask:
        ax.contourf(
            grid.X_grid, grid.Y_grid, np.where(grid.mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5], colors=land_clr, alpha=.25
        )
        ax.contourf(
            grid.X_grid,
            grid.Y_grid,
            np.where(grid.mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5],
            colors=land_clr,
            alpha=0.5,
        )
    
    return ax

def rose_plot(
    df_data: pd.DataFrame,
    grid: Grid,
    *,
    mag: str = "TiltDis",
    theta: str = "TiltDir",
    frac: float = 2.6,
    mag_bins=(0, 10, 20, 30, 40, np.inf),
    direction_offset: float = -20.0,
    shelf_offset: float = 80.0,
    lon_split: float = 157.0,
    lat_split: float = -33.5,
    legend_title: str = "tilt dist. (km)",
    show: bool = True,
    rtick_flag=False,
    cmaps=['Reds','Blues'],
    step=100,
    map_norm_factor=None,
    region_bins=None,
    figsize=(20,8)
):
    required={"Cyc",mag,theta}
    missing=required-set(df_data.columns)
    if missing:
        raise KeyError(f"df_data is missing required columns: {sorted(missing)}")

    region_map={1:"S1",2:"S2",3:"U1",4:"D1",5:"U2",6:"D2"}
    bin_grid = region_bins
    region_mask_grid = np.isin(bin_grid, [1, 2])
    df_plot = df_data.copy()
    df_plot['bin_id'] = df_plot.Region.map({v:k for k,v in region_map.items()})
    df_plot=df_plot.dropna(subset=["bin_id",mag,theta])
    df_plot["bin_id"]=df_plot["bin_id"].astype(int)

    colors_cmps=[
        plt.get_cmap(cmaps[0])(np.linspace(0,1,len(mag_bins)-1)),
        plt.get_cmap(cmaps[1])(np.linspace(0,1,len(mag_bins)-1)),
    ]

    cell_w=(grid.X_grid.max()-grid.X_grid.min())/3
    cell_h=(grid.Y_grid.max()-grid.Y_grid.min())/4
    cmap_bins=plt.cm.gist_ncar
    levels_bins=np.arange(.5,7.5,1)
    norm_bins=BoundaryNorm(levels_bins,cmap_bins.N)

    def get_bin_color(b,alpha=.25):
        color=list(cmap_bins(norm_bins(b))); color[-1]=alpha
        return color

    def windrose_counts_local(directions_deg,magnitudes):
        data=windrose_counts(directions_deg,magnitudes,mag_bins=mag_bins)
        if data is None:
            return None,None,None,None
        c,angles,width=data
        return c,angles,width,c.shape[1]

    def add_windrose(ax,x0,y0,data,colors,rmax):
        c,angles,width,k=data
        if c is None:
            return None
        size=frac*min(cell_w,cell_h)
        iax=ax.inset_axes([x0-size/2,y0-size/2,size,size],
                         transform=ax.transData,projection="polar")
        bottom=np.zeros(k)
        for i in range(len(mag_bins)-1):
            iax.bar(angles,c[i],width=width,bottom=bottom,color=colors[i],
                    edgecolor=(0,0,0,.1))
            bottom+=c[i]
        iax.set_rlim(0,rmax)
        iax.set_theta_zero_location("N"); iax.set_theta_direction(-1)
        iax.set_xticks([]); iax.set_yticks([]); iax.set_frame_on(False)
        return iax

    def plot_standalone_windrose(ax,data,colors,b,title="",tick_flag=False,
                                 rmax=None,rtick_flag=True):
        c,angles,width,k=data
        if c is None:
            ax.set_axis_off()
            return
        ax.set_facecolor(get_bin_color(b,alpha=.25))
        bottom=np.zeros(k)
        for i in range(len(mag_bins)-1):
            hi=r"$\infty$" if np.isinf(mag_bins[i+1]) else f"{mag_bins[i+1]:g}"
            label=f"[{mag_bins[i]:g}-{hi})"
            ax.bar(angles,c[i],width=width,bottom=bottom,color=colors[i],
                   edgecolor=(0,0,0,.1),label=label)
            bottom+=c[i]

        if rtick_flag:
            local_rmax=np.max(c.sum(axis=0))
            local_rmax=1 if local_rmax==0 else local_rmax
            ax.set_rlim(0,rmax if rmax is not None else local_rmax+5)
            top=int(np.ceil(local_rmax/step)*step)
            rticks=np.arange(step,top+step,step)
            ax.set_rticks(rticks)
            ax.set_yticklabels([f"{t:g}" for t in rticks],fontsize=8)

        ax.set_rlabel_position(135)
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        ax.set_xticks(np.deg2rad([340,70,160,250]))
        ax.set_xticklabels(["N","E","S","W"],fontsize=9)
        ax.set_title(title,fontsize=11)

        if tick_flag:
            ax.legend(title=f"{title[:2]}\n{legend_title}",loc="center left",
                      bbox_to_anchor=(.95,1.5),frameon=True,fontsize=9,title_fontsize=9)

    bin_centers={}
    for b in range(1,7):
        ii,jj=np.where(bin_grid==b)
        bin_centers[b]=(np.nanmean(grid.X_grid[ii,jj]),np.nanmean(grid.Y_grid[ii,jj]))

    # raw counts
    counts_raw={}
    for cyc in ["AE","CE"]:
        for b in range(1,7):
            sub=df_plot[(df_plot.Cyc==cyc)&(df_plot.bin_id==b)]
            directions=(sub[theta].to_numpy(float)+direction_offset)%360
            magnitudes=sub[mag].to_numpy(float)
            counts_raw[(cyc,b)]=windrose_counts_local(directions,magnitudes)

    # fixed map radial scale from raw counts
    rmax_map=0
    for cyc in ["AE","CE"]:
        for b in [3,4,5,6]:
            c=counts_raw[(cyc,b)][0]
            if c is not None:
                rmax_map=max(rmax_map,np.max(c.sum(axis=0)))
    rmax_map=1 if rmax_map==0 else rmax_map

    # plotting counts
    counts={}
    for cyc in ["AE","CE"]:
        for b in range(1,7):
            c,angles,width,k=counts_raw[(cyc,b)]

            if map_norm_factor is not None and b in [3,4,5,6] and c is not None:
                peak=c.sum(axis=0).max()
                if peak>0:
                    c=c*(map_norm_factor/peak)

            counts[(cyc,b)]=(c,angles,width,k)

    fig=plt.figure(figsize=figsize,constrained_layout=False)
    gs=fig.add_gridspec(1,4,width_ratios=[1.2,2.6,1.2,2.6],wspace=.45)
    gs_AE_small=gs[0,0].subgridspec(2,1)
    gs_CE_small=gs[0,2].subgridspec(2,1)

    small_axes={
        ("AE",1):fig.add_subplot(gs_AE_small[0,0],projection="polar"),
        ("AE",2):fig.add_subplot(gs_AE_small[1,0],projection="polar"),
        ("CE",1):fig.add_subplot(gs_CE_small[0,0],projection="polar"),
        ("CE",2):fig.add_subplot(gs_CE_small[1,0],projection="polar"),
    }
    axs={"AE":fig.add_subplot(gs[0,1]),"CE":fig.add_subplot(gs[0,3])}

    # S1/S2 unchanged
    for p,cyc in enumerate(["AE","CE"]):
        for b in [1,2]:
            plot_standalone_windrose(
                small_axes[(cyc,b)],counts[(cyc,b)],colors_cmps[p],b,
                title=f"{cyc}, S{b}",tick_flag=b==2,rtick_flag=rtick_flag
            )

    for p,cyc in enumerate(["AE","CE"]):
        ax=axs[cyc]

        ax.contourf(grid.X_grid,grid.Y_grid,
                    np.where(grid.mask_rho==0,1,np.nan),
                    levels=[.5,1.5],colors=["k"],alpha=.5)

        ax.contourf(grid.X_grid,grid.Y_grid,bin_grid,
                    levels=levels_bins,cmap=cmap_bins,norm=norm_bins,alpha=.25)

        ax=lat_lon_contours(ax,grid)

        ax.contour(grid.X_grid,grid.Y_grid,grid.h,
                   levels=[4000],colors="k",linewidths=1)

        ax.contour(grid.X_grid,grid.Y_grid,region_mask_grid.astype(float),
                   levels=[.5],colors="magenta",linewidths=2,linestyles="-")

        ax.contour(grid.X_grid,grid.Y_grid,grid.lon_rho,
                   levels=[lon_split],colors="magenta",linewidths=2,linestyles="-")

        ax.contour(grid.X_grid,grid.Y_grid,
                   np.where(grid.mask_rho==1,grid.lat_rho,np.nan),
                   levels=[lat_split],colors="magenta",linewidths=2,linestyles="-")

        for b in [3,4,5,6]:
            x0,y0=bin_centers[b]
            add_windrose(ax,x0,y0,counts[(cyc,b)],colors_cmps[p],rmax_map)

        for label,x,y in [
            ("S1",220,1300),("S2",120,50),
            ("U1",400,1450),("U2",800,1450),
            ("D1",620,750),("D2",850,750),
        ]:
            ax.text(x,y,label,ha="center",va="center",fontsize=12,fontweight="bold")

        ax.set_aspect("equal")
        ax.set_xlim(grid.X_grid.min(),grid.X_grid.max())
        ax.set_ylim(grid.Y_grid.min(),grid.Y_grid.max())
        ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
        ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")

    for key,label in [
        (("AE",1),"a)"),(("AE",2),"b)"),
        (("CE",1),"d)"),(("CE",2),"e)")
    ]:
        small_axes[key].text(-.15,1.02,label,transform=small_axes[key].transAxes,
                             fontsize=12,fontweight="bold",va="top",ha="left")

    axs["AE"].text(-.08,1.02,"c)",transform=axs["AE"].transAxes,
                   fontsize=12,fontweight="bold",va="top",ha="left")
    axs["CE"].text(-.08,1.02,"f)",transform=axs["CE"].transAxes,
                   fontsize=12,fontweight="bold",va="top",ha="left")

    if show:
        plt.show()

    return fig,{"map_axes":axs,"shelf_axes":small_axes,
                "counts":counts,"counts_raw":counts_raw,
                "bin_grid":bin_grid,"rmax_map":rmax_map}
