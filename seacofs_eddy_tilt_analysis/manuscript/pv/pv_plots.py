"""Original Figure 2 plotting style, adapted locally from tagged legacy helpers."""
from __future__ import annotations
import numpy as np
import pandas as pd
import pv_kernels as kernels

def geographic_to_grid(east, north, angle):
    """Inverse ROMS rotation; angle is radians from geographic east to grid +x.

    Geographic = [[cos(a), -sin(a)], [sin(a), cos(a)]] @ grid.
    Apply its transpose for display on model x/y axes; magnitudes are unchanged.
    """
    if not np.isscalar(angle) or not np.isfinite(angle):
        raise ValueError('Need a finite scalar grid angle in radians')
    c, s = np.cos(angle), np.sin(angle)
    east, north = np.asarray(east), np.asarray(north)
    return c*east+s*north, -s*east+c*north

def plot_ellipse(ax, row, grid, frac=1, color='k', lw=1, zorder=None, alpha=1):
    Q = np.array([[row.q11, row.q12], [row.q12, row.q22]], dtype=float)
    dx, dy = grid.X_grid - row.xc, grid.Y_grid - row.yc
    rho2 = Q[0, 0]*dx**2 + 2*Q[0, 1]*dx*dy + Q[1, 1]*dy**2
    ax.contour(grid.X_grid, grid.Y_grid, rho2, levels=[(row.Rc**2/2)*frac**2],
               colors=[color], linewidths=lw, zorder=zorder, alpha=alpha)
    return

def environmental_grid_gradients(grid):
    """Return bathymetric slopes and beta in true east/north coordinates."""

    dhdx, dhdy = kernels.phys_grad(
        grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
    )
    dfdx, dfdy = kernels.phys_grad(
        grid.f, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
    )
    return {
        "dh_e": np.cos(grid.angle) * dhdx - np.sin(grid.angle) * dhdy,
        "dh_n": np.sin(grid.angle) * dhdx + np.cos(grid.angle) * dhdy,
        "beta": np.sin(grid.angle) * dfdx + np.cos(grid.angle) * dfdy,
    }

def local_esp_fields(row, grid, frac=1.0, grid_gradients=None):
    """Return Laplacian-derived local fields for one eddy snapshot."""

    ii, jj = kernels.core_grid_indices(row, grid, frac=frac)
    dx = grid.x_grid[ii] - float(row.xc)
    dy = grid.y_grid[jj] - float(row.yc)
    reconstruction = kernels._esp_gaussian_vorticity(
        row.w, row.Rc, row.q11, row.q12, row.q22, dx, dy
    )
    rho2 = reconstruction["rho2"]
    shape = reconstruction["gaussian"]
    zeta = reconstruction["zeta"]
    gradients = (
        environmental_grid_gradients(grid)
        if grid_gradients is None else grid_gradients
    )
    dh_e, dh_n, beta = (
        gradients["dh_e"], gradients["dh_n"], gradients["beta"]
    )
    h, f = grid.h[ii, jj], grid.f[ii, jj]
    plan_e, plan_n = np.zeros(len(ii)), beta[ii, jj]/h
    topo_e = -(f+zeta)*dh_e[ii, jj]/h**2
    topo_n = -(f+zeta)*dh_n[ii, jj]/h**2
    eddy_e = reconstruction["dzeta_dx"]/h
    eddy_n = reconstruction["dzeta_dy"]/h
    return pd.DataFrame({
        "i": ii, "j": jj, "x": grid.X_grid[ii,jj], "y": grid.Y_grid[ii,jj],
        "rho_over_Rc": np.sqrt(rho2)/float(row.Rc), "weight": shape,
        "vorticity_shape": reconstruction["laplacian_shape"],
        "zeta": zeta, "PV": (f+zeta)/h,
        "environment_east": plan_e+topo_e, "environment_north": plan_n+topo_n,
        "environment_mag": np.hypot(plan_e+topo_e, plan_n+topo_n),
        "eddy_east": eddy_e, "eddy_north": eddy_n,
        "eddy_mag": np.hypot(eddy_e, eddy_n),
        "full_east": plan_e+topo_e+eddy_e, "full_north": plan_n+topo_n+eddy_n,
        "full_mag": np.hypot(plan_e+topo_e+eddy_e, plan_n+topo_n+eddy_n),
    })

def plot_gradient_examples(
    rows,
    grid,
    *,
    vertical=None,
    vert_lim=1e3,
    min_pad_km=40.0,
    pad_rc=1.4,
    max_local_arrows=120,
    reference_percentile=90.0,
    local_arrow_length_inches=0.38,
    net_arrow_fraction=0.38,
    figsize=None,
    legend=True,
    constrained_layout=True,
    xq=.02,
    yq=1.02,
    title=True,
    quiv_clr='limegreen'
):
    """Plot geographic vectors rotated into model x/y, retaining original style."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    rows = rows.reset_index(drop=True)
    if rows.empty:
        raise ValueError("No candidate rows were supplied")
    gradients = environmental_grid_gradients(grid)
    prepared = []
    all_depths = []
    for _, row in rows.iterrows():
        local = local_esp_fields(
            row, grid, frac=1.0, grid_gradients=gradients
        )
        pad = max(float(min_pad_km), float(pad_rc) * float(row.Rc))
        inside = (
            (grid.X_grid >= row.xc - pad) & (grid.X_grid <= row.xc + pad)
            & (grid.Y_grid >= row.yc - pad) & (grid.Y_grid <= row.yc + pad)
            & grid.mask_rho.astype(bool)
        )
        bathy = np.where(inside, grid.h / 1000.0, np.nan)
        local_magnitude = local.environment_mag.to_numpy(float)
        local_magnitude = local_magnitude[
            np.isfinite(local_magnitude) & (local_magnitude > 0)
        ]
        if not len(local_magnitude):
            raise ValueError(
                f"{row.Cyc}{int(row.Eddy)} day {int(row.Day)} has no finite "
                "local PV-gradient vectors"
            )
        local_reference = float(
            np.nanpercentile(local_magnitude, reference_percentile)
        )
        prepared.append((row, local, pad, bathy, local_reference))
        all_depths.append(bathy[np.isfinite(bathy)])
    depths = np.concatenate([values for values in all_depths if len(values)])
    depth_norm = Normalize(vmin=float(np.nanmin(depths)), vmax=float(np.nanmax(depths)))

    n = len(rows)
    # figsize = figsize or (5.8 * n, 5.3)
    # fig, axes = plt.subplots(1, n, figsize=figsize, squeeze=False,
    #                          constrained_layout=True)
    ncols = min(n, 5)
    nrows = int(np.ceil(n / ncols))
    
    figsize = figsize or (5.8 * ncols, 5.3 * nrows)
    
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=figsize,
        squeeze=False,
        constrained_layout=constrained_layout
    )

    axes = axes.flatten()
    bathy_artist = None
    local_references = []
    for ax, (row, local, pad, bathy, local_reference) in zip(axes, prepared):
        bathy_artist = ax.pcolormesh(
            grid.X_grid, grid.Y_grid, bathy, cmap="Greys_r",
            norm=depth_norm, shading="auto",
        )
        ax.contour(grid.X_grid, grid.Y_grid, bathy, colors="0.35",
                   levels=8, linewidths=0.45, alpha=0.55)
        plot_ellipse(ax, row, grid, frac=1, color="red" if row.Cyc == "AE" else "dodgerblue", lw=2.0, zorder=8)
        step = max(1, int(np.ceil(len(local) / max_local_arrows)))
        arrows = local.iloc[::step]
        valid = np.isfinite(arrows.environment_east) & np.isfinite(arrows.environment_north)
        local_scale = local_reference / float(local_arrow_length_inches)
        arrow_x, arrow_y = geographic_to_grid(
            arrows.environment_east[valid], arrows.environment_north[valid], grid.angle)
        local_quiver = ax.quiver(
            arrows.x[valid], arrows.y[valid],
            arrow_x, arrow_y,
            color=quiv_clr,
            alpha=0.68, angles="xy", scale_units="inches",
            scale=local_scale, width=0.004, zorder=7,
        )
        ax.quiverkey(
            local_quiver, xq, yq, local_reference,
            rf"${local_reference:.1e}$ m$^{{-2}}$ s$^{{-1}}$",
            coordinates="axes", labelpos="E",
        )
        local_references.append(local_reference)
        net_east, net_north = float(row.PV_grad_x), float(row.PV_grad_y)
        net = float(np.hypot(net_east, net_north))
        net_x, net_y = geographic_to_grid(net_east, net_north, grid.angle)
        if np.isfinite(net) and net > 0:
            fixed_length = float(net_arrow_fraction) * pad
            ax.annotate(
                "",
                xy=(row.xc + fixed_length * net_x / net,
                    row.yc + fixed_length * net_y / net),
                xytext=(row.xc, row.yc),
                arrowprops=dict(arrowstyle="-|>", color="limegreen", lw=3),
                zorder=10,
            )
        ax.scatter(row.xc, row.yc, c="limegreen", s=24, zorder=11)
        if vertical is not None:
            spine = vertical[
                vertical.Eddy.eq(row.Eddy) & vertical.Day.eq(row.Day) & (vertical.Depth<vert_lim)
            ].sort_values("Depth")
            ax.plot(spine.xc, spine.yc, color="red" if row.Cyc == "AE" else "dodgerblue", lw=2.0, zorder=9)
        local_mean = float(row.PV_grad_mean_local_mag)
        ax.set(
            xlim=(row.xc-pad, row.xc+pad), ylim=(row.yc-pad, row.yc+pad),
            aspect="equal", xlabel="x (km)", ylabel="y (km)",
        )
        if title:
            ax.set_title(f"{row.Cyc}{int(row.Eddy)}, day {int(row.Day)}\n"
                   f"net/local = {net/local_mean:.2f}" if local_mean > 0
                   else f"{row.Cyc}{int(row.Eddy)}, day {int(row.Day)}")
        ax.plot([], [], color="limegreen", lw=3,
                label=rf"mean $|\nabla q|={net:.2e}$ m$^{{-2}}$ s$^{{-1}}$")
        ax.set_axisbelow(True)
        ax.grid(True)
        if legend:
            ax.legend(loc="lower right", frameon=True, fontsize=9)
    fig.colorbar(bathy_artist, ax=axes.tolist(), label="Bathymetry (km)", shrink=0.82)
    return fig, axes, np.asarray(local_references)
