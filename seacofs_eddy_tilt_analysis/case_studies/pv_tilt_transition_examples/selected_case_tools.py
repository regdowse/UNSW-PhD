"""Plotting helpers for selected fixed-core PV-gradient/tilt case studies."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


SELECTED_EDDIES = {
    "AE": [356, 484, 1260, 1015, 1268, 2255, 604, 739],
    "CE": [1927, 1194, 1552, 1105],
}
DEFAULT_CACHE_ROOT = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/"
    "pv_gradient_selected_case_studies"
)
DEFAULT_CACHE_NAME = "selected_eddy_surface_esp_gaussian_frac1.parquet"
REGIME_COLOURS = {
    "planetary": "#2ca25f",
    "mixed": "#8c8c8c",
    "topographic": "#f28e2b",
}


def cache_path(cache_root=DEFAULT_CACHE_ROOT):
    return Path(cache_root) / DEFAULT_CACHE_NAME


def select_tracks(df, selected=SELECTED_EDDIES):
    """Select exact polarity/eddy pairs and retain the requested order."""
    pieces = []
    for cyc, eddies in selected.items():
        for order, eddy in enumerate(eddies):
            part = df[df.Cyc.eq(cyc) & df.Eddy.eq(eddy)].copy()
            if not part.empty:
                part["selection_order"] = order
                pieces.append(part)
    if not pieces:
        return df.iloc[0:0].copy()
    return pd.concat(pieces, ignore_index=True).sort_values(
        ["Cyc", "selection_order", "Day"]
    ).reset_index(drop=True)


def add_regimes(df, dominance_factor=2.0, smooth_window=7):
    """Classify sustained planetary, mixed and topographic conditions."""
    if dominance_factor <= 1:
        raise ValueError("dominance_factor must exceed one")
    if smooth_window < 1 or smooth_window % 2 == 0:
        raise ValueError("smooth_window must be a positive odd integer")
    out = df.sort_values(["Cyc", "Eddy", "Day"]).copy()
    ratio = pd.to_numeric(out.topo_plan_ratio, errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    )
    out["topo_plan_ratio"] = ratio
    out["topo_plan_ratio_smooth"] = out.groupby(
        ["Cyc", "Eddy"], sort=False
    ).topo_plan_ratio.transform(
        lambda x: x.rolling(smooth_window, center=True, min_periods=1).median()
    )
    threshold = np.log(dominance_factor)
    out["regime"] = np.select(
        [out.topo_plan_ratio_smooth <= -threshold,
         out.topo_plan_ratio_smooth >= threshold],
        ["planetary", "topographic"], default="mixed",
    )
    return out


def build_cache(eddies, grid, selected=SELECTED_EDDIES):
    """Calculate the fixed-core ESP-Gaussian terms for selected eddies."""
    import seacofs_tilt_tools as tilt

    subset = select_tracks(eddies, selected)
    if subset.empty:
        raise ValueError("None of the selected eddies was found")
    result = tilt.add_pv_gradient_terms(
        subset, grid, core_mean=True, frac=1.0,
        surface_method="esp_gaussian", averaging="nonlinear",
    )
    return add_regimes(result)


def load_cache(cache_root=DEFAULT_CACHE_ROOT):
    import seacofs_tilt_tools as tilt
    return add_regimes(tilt.read_table(cache_path(cache_root)))


def _regime_spans(ax, age, regime):
    for x, label in zip(age, regime):
        ax.axvspan(x-.5, x+.5, color=REGIME_COLOURS[label], alpha=.16, lw=0)


def _vector_endpoint(x, y, magnitude, bearing, scale=1.0, reverse=False):
    """Endpoint from a compass bearing; reverse supports surface-to-deep tilt."""
    sign = -1.0 if reverse else 1.0
    angle = np.deg2rad(float(bearing))
    return (float(x) + sign*scale*float(magnitude)*np.sin(angle),
            float(y) + sign*scale*float(magnitude)*np.cos(angle))


def plot_eddy_overview(track, grid, *, arrow_count=14):
    """Compact time series plus bathymetric track/tilt/PV-vector map."""
    import seacofs_tilt_tools as tilt

    df = add_regimes(track).sort_values("Day").copy()
    if df.empty:
        raise ValueError("The requested eddy track is empty")
    age = df.Day - df.Day.iloc[0]
    cyc, eddy = df.Cyc.iloc[0], int(df.Eddy.iloc[0])
    fig = plt.figure(figsize=(15, 10), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[1.45, 2.2])
    axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]
    axm = fig.add_subplot(gs[:, 1])
    for ax in axes:
        _regime_spans(ax, age, df.regime)

    axes[0].plot(age, df.TiltDis, color="tab:purple", lw=1.8)
    axes[0].set_ylabel("Tilt distance (km)")
    theta = np.deg2rad(df.TiltDir)
    axes[1].plot(age, df.TiltDis*np.sin(theta), color="tab:red", label="Zonal")
    axes[1].plot(age, df.TiltDis*np.cos(theta), color="tab:blue", label="Meridional")
    axes[1].axhline(0, color=".3", lw=.8)
    axes[1].set_ylabel("Tilt component (km)")
    axes[1].legend(ncol=2, frameon=False)
    axes[2].semilogy(age, df.PV_grad_mag, color="tab:green", lw=1.7, label="Net")
    axes[2].semilogy(age, df.PV_grad_mean_local_mag, color="tab:green", alpha=.35,
                    label="Mean local")
    axes[2].set_ylabel(r"$|\nabla_h q|$")
    axes[2].legend(ncol=2, frameon=False)
    axes[3].plot(age, df.h/1e3, color="saddlebrown", lw=1.6)
    axes[3].set_ylabel("Depth (km)", color="saddlebrown")
    axes[3].invert_yaxis()
    ax_pv = axes[3].twinx()
    ax_pv.plot(age, df.PV, color="tab:blue", lw=1.3)
    ax_pv.set_ylabel("PV", color="tab:blue")
    axes[3].set_xlabel("Eddy age (days)")
    for ax in axes:
        ax.grid(alpha=.2); ax.margins(x=0)

    deep_x, deep_y = _vector_endpoint(
        df.xc, df.yc, df.TiltDis, df.TiltDir, reverse=True
    ) if len(df) == 1 else (
        df.xc-df.TiltDis*np.sin(theta), df.yc-df.TiltDis*np.cos(theta)
    )
    pad = max(25, float(np.nanmedian(df.Rc)))
    xmin=min(df.xc.min(),np.nanmin(deep_x))-pad; xmax=max(df.xc.max(),np.nanmax(deep_x))+pad
    ymin=min(df.yc.min(),np.nanmin(deep_y))-pad; ymax=max(df.yc.max(),np.nanmax(deep_y))+pad
    inside=((grid.X_grid>=xmin)&(grid.X_grid<=xmax)&(grid.Y_grid>=ymin)&(grid.Y_grid<=ymax))
    bathy=np.where(grid.mask_rho.astype(bool)&inside,grid.h/1e3,np.nan)
    cf=axm.contourf(grid.X_grid,grid.Y_grid,bathy,levels=24,cmap="Greys_r")
    fig.colorbar(cf,ax=axm,orientation="horizontal",location="bottom",
                 label="Depth (km)",shrink=.72,pad=.06)

    points=np.column_stack([df.xc,df.yc])
    if len(points)>1:
        segments=np.stack([points[:-1],points[1:]],axis=1)
        colours=[REGIME_COLOURS[r] for r in df.regime.iloc[:-1]]
        axm.add_collection(LineCollection(segments,colors=colours,linewidths=3,zorder=5))
    positions=np.unique(np.linspace(0,len(df)-1,min(arrow_count,len(df)),dtype=int))
    pv_scale=.65*float(np.nanmedian(df.Rc))
    for pos in positions:
        row=df.iloc[pos]
        tx,ty=_vector_endpoint(row.xc,row.yc,row.TiltDis,row.TiltDir,reverse=True)
        axm.annotate("",xy=(tx,ty),xytext=(row.xc,row.yc),
                    arrowprops=dict(arrowstyle="-|>",color="tab:blue",lw=1.5,alpha=.75),zorder=8)
        if np.isfinite(row.PV_grad_theta):
            px,py=_vector_endpoint(row.xc,row.yc,1,row.PV_grad_theta,scale=pv_scale)
            coherence=np.clip(row.PV_grad_coherence,0,1) if np.isfinite(row.PV_grad_coherence) else 0
            axm.annotate("",xy=(px,py),xytext=(row.xc,row.yc),
                        arrowprops=dict(arrowstyle="-|>",color="magenta",lw=1.7,
                                        alpha=.25+.75*coherence),zorder=9)
        tilt.plot_ellipse(axm,row,grid,frac=1,color=REGIME_COLOURS[row.regime],
                          lw=.7,alpha=.4,zorder=6)
    axm.scatter(df.xc.iloc[0],df.yc.iloc[0],facecolor="white",edgecolor="black",s=55,zorder=10)
    axm.scatter(df.xc.iloc[-1],df.yc.iloc[-1],marker="x",color="black",s=55,zorder=10)
    legend=[Line2D([0],[0],color=REGIME_COLOURS[k],lw=3,label=k.title()) for k in REGIME_COLOURS]
    legend += [Line2D([0],[0],color="tab:blue",lw=2,label="Surface-to-deep tilt"),
               Line2D([0],[0],color="magenta",lw=2,label="Mean PV gradient")]
    axm.legend(handles=legend,frameon=False,loc="best")
    axm.set(xlim=(xmin,xmax),ylim=(ymin,ymax),aspect="equal",xlabel="x (km)",ylabel="y (km)",
            title=f"{cyc}{eddy}: surface track, tilt and PV-gradient directions")
    fig.suptitle(f"{cyc}{eddy}: fixed-core ESP-Gaussian surface PV gradient",fontsize=15)
    return fig, axes, axm


def suggested_days(track):
    """Choose strongest planetary, closest transition and strongest topo days."""
    df=add_regimes(track).dropna(subset=["topo_plan_ratio_smooth"]).sort_values("Day")
    if df.empty:
        return []
    candidates=[df.topo_plan_ratio_smooth.idxmin(),
                df.topo_plan_ratio_smooth.abs().idxmin(),
                df.topo_plan_ratio_smooth.idxmax()]
    return [int(df.loc[i,"Day"]) for i in dict.fromkeys(candidates)]


def plot_selected_day(track, vertical, grid, day):
    """Map local core gradients, weighted mean vector and vertical eddy spine."""
    import seacofs_tilt_tools as tilt
    from surface_pv_esp_gaussian import esp_pv_tools as ept

    match=track[track.Day.eq(day)]
    if match.empty:
        raise ValueError(f"Day {day} is not available for this eddy")
    row=match.iloc[0]
    local=ept.local_esp_fields(row,grid,frac=1)
    profile=vertical[(vertical.Eddy.eq(row.Eddy))&(vertical.Day.eq(row.Day))].copy()
    profile=profile.dropna(subset=["Depth","xc","yc"]).sort_values("Depth")
    pad=max(35,1.35*float(row.Rc))
    xmin=row.xc-pad; xmax=row.xc+pad; ymin=row.yc-pad; ymax=row.yc+pad
    inside=((grid.X_grid>=xmin)&(grid.X_grid<=xmax)&(grid.Y_grid>=ymin)&(grid.Y_grid<=ymax))
    bathy=np.where(grid.mask_rho.astype(bool)&inside,grid.h/1e3,np.nan)
    fig,ax=plt.subplots(figsize=(7.2,6.4),constrained_layout=True)
    cf=ax.contourf(grid.X_grid,grid.Y_grid,bathy,levels=25,cmap="terrain_r")
    tilt.plot_ellipse(ax,row,grid,frac=1,color="cyan",lw=2.5,zorder=8)
    step=max(1,len(local)//130); arrows=local.iloc[::step]
    mag=np.hypot(arrows.environment_east,arrows.environment_north); valid=mag.gt(0)
    ax.quiver(arrows.x[valid],arrows.y[valid],arrows.environment_east[valid]/mag[valid],
              arrows.environment_north[valid]/mag[valid],color="white",alpha=.7,
              scale=24,width=.003,zorder=7)
    net=np.hypot(row.PV_grad_x,row.PV_grad_y)
    if np.isfinite(net) and net>0:
        px,py=_vector_endpoint(row.xc,row.yc,.75*row.Rc,row.PV_grad_theta)
        ax.annotate("",xy=(px,py),xytext=(row.xc,row.yc),
                    arrowprops=dict(arrowstyle="-|>",color="magenta",lw=3),zorder=11)
    if not profile.empty:
        sc=ax.scatter(profile.xc,profile.yc,c=profile.Depth,cmap="viridis_r",s=30,
                      edgecolor="black",linewidth=.2,zorder=10)
        ax.plot(profile.xc,profile.yc,color="black",lw=1.2,alpha=.75,zorder=9)
        fig.colorbar(sc,ax=ax,label="Eddy-centre depth (m)",shrink=.8,pad=.02)
    else:
        ax.scatter(row.xc,row.yc,c="black",s=35,zorder=10)
    ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax),aspect="equal",xlabel="x (km)",ylabel="y (km)",
           title=(f"{row.Cyc}{int(row.Eddy)}, day {int(day)} | {row.regime}; "
                  f"coherence={row.PV_grad_coherence:.2f}"))
    fig.colorbar(cf,ax=ax,label="Water-column depth (km)",shrink=.8,pad=.08)
    return fig, ax
