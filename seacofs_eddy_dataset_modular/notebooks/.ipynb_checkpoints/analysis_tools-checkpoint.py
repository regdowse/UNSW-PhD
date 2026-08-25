""" Toolkit for analysis of SEACOFS eddy dataset"""

import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
import sys
import matplotlib.cm as cm


TILT_ROOT = Path("~/UNSW-MRes/MRes/seacofs_eddy_tilt_analysis").expanduser()
if str(TILT_ROOT) not in sys.path:
    sys.path.insert(0, str(TILT_ROOT))
import seacofs_tilt_tools as tilt

def distribution_plot(df_data, grid):
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    
    # --- Panel 1: Map with eddy tracks and insets ---
    ax = axs[0]
    cs = ax.contourf(grid.lon_rho, grid.lat_rho,
                     np.where(grid.mask_rho, grid.h/1e3, np.nan), cmap='grey')
    fig.colorbar(cs, ax=ax, label='Bathymetry (km)')
    
    for eddy in df_data.Eddy.unique():
        d = df_data[df_data.Eddy == eddy]
        cyc = d.iloc[0].Cyc
        ax.plot(d.lon, d.lat, color='r' if cyc=='AE' else 'b', lw=1, alpha=1)
    
    ax.axis('equal')
    ax.set_xlim(148, 160)
    ax.set_ylim(-40, -26)
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    
    # Inset 1: Eddy counts
    df_unique = df_data[['Eddy', 'Cyc']].drop_duplicates()
    in_ax1 = inset_axes(ax, width=.4, height=1, loc='center',
                       bbox_to_anchor=(152, -27.5, 0, 0),
                       bbox_transform=ax.transData, borderpad=0)
    counts = df_unique['Cyc'].value_counts().reindex(['AE','CE']).fillna(0)
    counts.plot(kind='bar', color=['r','b'], ax=in_ax1)
    in_ax1.set_xlabel('')
    in_ax1.set_ylabel('No. eddies')
    in_ax1.tick_params(axis='x', rotation=0)
    in_ax1.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
    in_ax1.set_facecolor('none')
    
    # Inset 2: Eddy-day counts
    in_ax2 = inset_axes(ax, width=.4, height=1, loc='center',
                       bbox_to_anchor=(149.9, -32.3, 0, 0),
                       bbox_transform=ax.transData, borderpad=0)
    counts = df_data['Cyc'].value_counts().reindex(['AE','CE']).fillna(0)
    counts.plot(kind='bar', color=['r','b'], ax=in_ax2)
    in_ax2.set_xlabel('')
    in_ax2.set_ylabel('No. eddy-days')
    in_ax2.tick_params(axis='x', rotation=0)
    in_ax2.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
    in_ax2.set_facecolor('none')
    
    # --- Panel 2 & 3: Eddy-day distributions ---
    for d, subset in enumerate([df_data[df_data.Cyc == 'AE'], df_data[df_data.Cyc == 'CE']]):
        eddy_x = subset.xc.to_numpy()
        eddy_y = subset.yc.to_numpy()
    
        x_edges = tilt.bin_edges_fd(eddy_x, grid.X_grid, rule='fd')  # or rule='scott'/'fd'
        y_edges = tilt.bin_edges_fd(eddy_y, grid.Y_grid, rule='fd')
    
        H, _, _ = np.histogram2d(eddy_y, eddy_x, bins=[y_edges, x_edges])
    
        pcm = axs[d+1].pcolormesh(x_edges, y_edges, H, cmap='Reds' if d==0 else 'Blues')
        fig.colorbar(pcm, ax=axs[d+1], label='AE-day frequency' if d==0 else 'CE-day frequency')
    
        axs[d+1].contourf(grid.X_grid, grid.Y_grid, np.where(grid.mask_rho == 0, 1, np.nan),
                          levels=[0.5, 1.5], colors=['k'], alpha=.5)
    
        c1 = axs[d+1].contour(grid.X_grid, grid.Y_grid, grid.lat_rho,
                              levels=[-40, -35, -30, -25], colors='k', linewidths=.5)
        axs[d+1].clabel(c1, fmt=lambda v: f"{np.abs(v):.0f}°S", inline=True, colors='k')
        c2 = axs[d+1].contour(grid.X_grid, grid.Y_grid, grid.lon_rho,
                              levels=[150, 155, 160], colors='k', linewidths=.5)
        axs[d+1].clabel(c2, fmt=lambda v: f"{v:.0f}°E", inline=True, colors='k')
    
        axs[d+1].set_xlabel('x (km)')
        axs[d+1].set_ylabel('y (km)')
        axs[d+1].set_aspect('equal')
        axs[d+1].contour(grid.X_grid, grid.Y_grid, grid.h, levels=[4e3], colors='k')
    
    plt.show()
    return fig, axs


# def distribution_plot(
#     df_data,
#     grid,
#     *,
#     metric='Rc',
#     vmin=0,
#     vmax=120,
#     units='km',
#     rule='fd',
#     cbar_loc='right'
# ):

#     fig, axs = plt.subplots(
#         1, 3,
#         figsize=(15, 5),
#         constrained_layout=True
#     )

#     # ============================================================
#     # Panel 1: Eddy tracks + counts
#     # ============================================================

#     ax = axs[0]

#     cs = ax.contourf(
#         grid.lon_rho,
#         grid.lat_rho,
#         np.where(grid.mask_rho, grid.h / 1e3, np.nan),
#         cmap='grey'
#     )

#     fig.colorbar(
#         cs,
#         ax=ax,
#         label='Bathymetry (km)'
#     )

#     for eddy in df_data.Eddy.unique():

#         d = df_data[df_data.Eddy == eddy]
#         cyc = d.iloc[0].Cyc

#         ax.plot(
#             d.lon,
#             d.lat,
#             color='r' if cyc == 'AE' else 'b',
#             lw=1,
#             alpha=1
#         )

#     ax.axis('equal')
#     ax.set_xlim(148, 160)
#     ax.set_ylim(-40, -26)

#     ax.set_xlabel('Longitude (°E)')
#     ax.set_ylabel('Latitude (°N)')

#     # ------------------------------------------------------------
#     # Eddy counts
#     # ------------------------------------------------------------

#     df_unique = df_data[['Eddy', 'Cyc']].drop_duplicates()

#     in_ax1 = inset_axes(
#         ax,
#         width=.4,
#         height=1,
#         loc='center',
#         bbox_to_anchor=(152, -27.5, 0, 0),
#         bbox_transform=ax.transData,
#         borderpad=0
#     )

#     counts = (
#         df_unique['Cyc']
#         .value_counts()
#         .reindex(['AE', 'CE'])
#         .fillna(0)
#     )

#     counts.plot(
#         kind='bar',
#         color=['r', 'b'],
#         ax=in_ax1
#     )

#     in_ax1.set_xlabel('')
#     in_ax1.set_ylabel('No. eddies')
#     in_ax1.tick_params(axis='x', rotation=0)
#     in_ax1.ticklabel_format(
#         axis='y',
#         style='sci',
#         scilimits=(3, 3)
#     )
#     in_ax1.set_facecolor('none')

#     # ------------------------------------------------------------
#     # Eddy-day counts
#     # ------------------------------------------------------------

#     in_ax2 = inset_axes(
#         ax,
#         width=.4,
#         height=1,
#         loc='center',
#         bbox_to_anchor=(149.9, -32.3, 0, 0),
#         bbox_transform=ax.transData,
#         borderpad=0
#     )

#     counts = (
#         df_data['Cyc']
#         .value_counts()
#         .reindex(['AE', 'CE'])
#         .fillna(0)
#     )

#     counts.plot(
#         kind='bar',
#         color=['r', 'b'],
#         ax=in_ax2
#     )

#     in_ax2.set_xlabel('')
#     in_ax2.set_ylabel('No. eddy-days')
#     in_ax2.tick_params(axis='x', rotation=0)
#     in_ax2.ticklabel_format(
#         axis='y',
#         style='sci',
#         scilimits=(3, 3)
#     )
#     in_ax2.set_facecolor('none')

#     # ============================================================
#     # Panels 2–3: Binned median metric
#     # ============================================================

#     tilt.plot_binned_median_map(
#         df_data,
#         grid,
#         metric=metric,
#         vmin=vmin,
#         vmax=vmax,
#         units=units,
#         rule=rule,
#         fig=fig,
#         axs=axs[1:],
#         show=False,
#         cbar_loc=cbar_loc
#     )

#     plt.show()

#     return fig, axs



def summary_plot(df_data):
    
    fig, axs = plt.subplots(8, 1, figsize=(10, 18))
    
    def clean(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        return a[np.isfinite(a)], b[np.isfinite(b)]
    
    # --- Plot 0: Age ---
    ae = df_data.loc[df_data.Cyc == 'AE', 'Age']
    ce = df_data.loc[df_data.Cyc == 'CE', 'Age']
    bins0 = tilt.shared_bins(ae, ce, min_bins=20, max_bins=50)
    tilt.mirrored_hist(
        axs[0], ae, ce, bins0,
        xlabel='Lifespan (days)', ylabel='Number of eddies',
        # ylim=(-25000, 25000)
    )
    
    # --- Plot 1: Lat ---
    ae = df_data.loc[df_data.Cyc == 'AE', 'lat']
    ce = df_data.loc[df_data.Cyc == 'CE', 'lat']
    bins1 = tilt.shared_bins(ae, ce, min_bins=20, max_bins=50)
    tilt.mirrored_hist(
        axs[1], ae, ce, bins1,
        xlabel='Latitude (°N)', ylabel='Number of eddy-days',
        # ylim=(-4000, 4000)
    )
    
    # --- Plot 2: Propagation Distance ---
    eddy_props_AE, eddy_props_CE = [], []
    for eddy in df_data.Eddy.unique():
        df = df_data[df_data.Eddy == eddy]
        dist = np.nansum(np.hypot(df.xc.diff(), df.yc.diff()))
        if df.iloc[0].Cyc == 'AE':
            eddy_props_AE.append(dist)
        else:
            eddy_props_CE.append(dist)
    
    bins2 = tilt.shared_bins(eddy_props_AE, eddy_props_CE, min_bins=15, max_bins=40)
    tilt.mirrored_hist(
        axs[2], eddy_props_AE, eddy_props_CE, bins2,
        xlabel='Propagation distance (km)', ylabel='Number of eddies',
        # ylim=(-450, 450)
    )
    
    # --- Plot 3: Propagation Speed ---

    eddy_props = []
    for eddy in df_data.Eddy.unique():
        df = df_data[df_data.Eddy==eddy].copy()
        eddy_props.extend(np.array(np.hypot(df.xc.diff(), df.yc.diff()) * 0.011574))
    df_data['EddyProp'] = eddy_props
    
    dfAE = df_data[df_data.Cyc == 'AE'].copy().dropna()
    dfAE = dfAE[dfAE.EddyProp <= 1]
    dfCE = df_data[df_data.Cyc == 'CE'].copy().dropna()
    dfCE = dfCE[dfCE.EddyProp <= 1]
    
    bins3 = tilt.shared_bins(dfAE.EddyProp, dfCE.EddyProp, min_bins=15, max_bins=40)
    tilt.mirrored_hist(
        axs[3], dfAE.EddyProp, dfCE.EddyProp, bins3,
        xlabel=r'Propagation speed (ms$^{-1}$)', ylabel='Number of eddies',
        # ylim=(-13000, 13000)
    )
    
    # --- Plot 4: Vorticity ---
    df_v = df_data[df_data.w.abs() <= 8e-5]
    ae = df_v.loc[df_v.Cyc == 'AE', 'w']
    ce = df_v.loc[df_v.Cyc == 'CE', 'w']
    bins4 = tilt.shared_bins(ae, ce, min_bins=20, max_bins=50)
    tilt.mirrored_hist(
        axs[4], ae, ce, bins4,
        xlabel=r'Surface vorticity (s$^{-1}$)', ylabel='Number of eddy-days',
        xlim=(-8e-5, 8e-5),
        # mirror_flag=False
    )
    
    # --- Plot 5: Radius ---
    ae = df_data.loc[df_data.Cyc == 'AE', 'Rc']
    ce = df_data.loc[df_data.Cyc == 'CE', 'Rc']
    bins5 = tilt.shared_bins(ae, ce, min_bins=20, max_bins=50)
    tilt.mirrored_hist(
        axs[5], ae, ce, bins5,
        xlabel=r'Surface core radius (km)', ylabel='Number of eddy-days',
        # ylim=(-7500, 7500)
    )
    
    # --- Plot 6: Depth ---
    bin_int = 1
    max_depth_km = 5
    edges = np.arange(0, max_depth_km + bin_int, bin_int)
    
    def _fmt(x):
        return f"{x:.0f}" if float(bin_int).is_integer() else f"{x:.1f}"
    
    bin_labels = (
        ["Depth = 0", f"({_fmt(edges[0])}–{_fmt(edges[0] + bin_int)})"]
        + [f"[{_fmt(b)}–{_fmt(b+bin_int)})" for b in edges[1:-1]]
    )
    bin_labels[-1] = f"[{_fmt(edges[-2])}–{_fmt(edges[-1])}]"
    
    ae_raw = df_data.loc[df_data.Cyc == 'AE',
        'max_profile_depth'].dropna().to_numpy()
    ce_raw = df_data.loc[df_data.Cyc == 'CE',
        'max_profile_depth'].dropna().to_numpy()
    
    ae = ae_raw[(ae_raw >= edges[0]) & (ae_raw <= edges[-1])]
    ce = ce_raw[(ce_raw >= edges[0]) & (ce_raw <= edges[-1])]
    
    ae0 = int(np.isclose(ae, 0).sum())
    ce0 = int(np.isclose(ce, 0).sum())
    
    ae_nz = ae[~np.isclose(ae, 0)]
    ce_nz = ce[~np.isclose(ce, 0)]
    
    ae_counts, _ = np.histogram(ae_nz, bins=edges)
    ce_counts, _ = np.histogram(ce_nz, bins=edges)
    
    ae_counts = np.insert(ae_counts, 0, ae0)
    ce_counts = np.insert(ce_counts, 0, ce0)
    
    ax = axs[6]
    xpos = np.arange(len(bin_labels))
    ax.bar(xpos,  ae_counts, color='r')
    ax.bar(xpos, -ce_counts, color='b')
    
    ax.axhline(0, color='k', linewidth=1)
    m = max(1, int(max(ae_counts.max(), ce_counts.max()) * 1.1))
    ax.set_xticks(xpos)
    ax.set_xticklabels(bin_labels, rotation=45, ha='right')
    ax.set_xlabel('Depth bin (km)')
    ax.set_ylabel('Number of eddy-days')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels([abs(int(t)) for t in ax.get_yticks()])
    ax.set_ylim(-m, m)
    
    # --- Plot 7: Aspect ratio ---
    df_a = df_data[df_data.AR <= 5]
    ae = df_a.loc[df_a.Cyc == 'AE', 'AR']
    ce = df_a.loc[df_a.Cyc == 'CE', 'AR']
    bins7 = tilt.shared_bins(ae, ce, min_bins=15, max_bins=40)
    tilt.mirrored_hist(
        axs[7], ae, ce, bins7,
        xlabel='Surface aspect ratio', ylabel='Number of eddy-days',
        # ylim=(-11000, 11000)
    )
    
    plt.tight_layout()
    plt.show()
    return fig, axs

def eddy_trajectories(df_data):
    fig, axs = plt.subplots(1, 2, figsize=(10,5), sharey=True)
    for ax in axs:
        ax.axhline(y=0, color='k', ls='--', alpha=.5, lw=.5)
        ax.axvline(x=0, color='k', ls='--', alpha=.5, lw=.5)
    # separate groups
    ae = [e for e in df_data.Eddy.unique() if df_data[df_data.Eddy==e].iloc[0].Cyc == 'AE']
    ce = [e for e in df_data.Eddy.unique() if df_data[df_data.Eddy==e].iloc[0].Cyc != 'AE']
    
    # colour ranges
    colors_ae = cm.Reds(np.linspace(0.4, 1, len(ae)))
    colors_ce = cm.Blues(np.linspace(0.4, 1, len(ce)))
    
    for c, eddy in zip(colors_ae, ae):
        df = df_data[df_data.Eddy==eddy]
        axs[0].plot(df.lon-df.iloc[0].lon, df.lat-df.iloc[0].lat, color=c, alpha=0.8)
    
    for c, eddy in zip(colors_ce, ce):
        df = df_data[df_data.Eddy==eddy]
        axs[1].plot(df.lon-df.iloc[0].lon, df.lat-df.iloc[0].lat, color=c, alpha=0.8)
    
    axs[0].axis('equal'); axs[1].axis('equal')
    axs[0].set_title('AE trajectories')
    axs[1].set_title('CE trajectories')
    axs[0].set_xlabel('(Degrees East)'); axs[0].set_ylabel('(Degrees North)')
    axs[1].set_xlabel('(Degrees East)')
    
    plt.tight_layout()
    plt.show()
    return fig, axs