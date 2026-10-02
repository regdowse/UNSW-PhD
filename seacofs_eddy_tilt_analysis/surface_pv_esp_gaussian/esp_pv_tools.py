"""Analysis helpers for ESP-Gaussian surface PV-gradient reconstruction."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_CACHE_ROOT = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/"
    "pv_gradient_surface_esp_gaussian"
)
DEFAULT_CACHE_NAME = "surface_pv_esp_gaussian_laplacian_v2.parquet"


def method_specs():
    """Primary fixed-core method and its equal-weight control."""
    return pd.DataFrame([
        {"method": "uniform_1", "surface_method": "uniform", "frac": 1.0},
        {"method": "esp_gaussian_1", "surface_method": "esp_gaussian", "frac": 1.0},
    ])


def cache_path(cache_root=DEFAULT_CACHE_ROOT):
    return Path(cache_root) / DEFAULT_CACHE_NAME


def build_cache(eddies, grid, specs=None):
    """Calculate the control and ESP-Gaussian surface methods."""
    import seacofs_tilt_tools as tilt

    specs = method_specs() if specs is None else specs.copy()
    tables = []
    for number, spec in enumerate(specs.itertuples(index=False), 1):
        print(f"[{number}/{len(specs)}] {spec.method}")
        table = tilt.add_pv_gradient_terms(
            eddies, grid, core_mean=True, averaging="nonlinear",
            frac=spec.frac, surface_method=spec.surface_method,
        )
        table["method"] = spec.method
        tables.append(table)
    return add_analysis_terms(pd.concat(tables, ignore_index=True))


def load_cache(cache_root=DEFAULT_CACHE_ROOT):
    import seacofs_tilt_tools as tilt
    return add_analysis_terms(tilt.read_table(cache_path(cache_root)))


def add_analysis_terms(df):
    """Add polarity-aware environmental and full-PV directional errors."""
    out = df.copy()
    for label, direction in (("environment", "PV_grad_theta"), ("full", "PV_grad_full_theta")):
        raw = np.abs((out.TiltDir - out[direction] + 180) % 360 - 180)
        out[f"{label}_preference_error_deg"] = np.where(out.Cyc.eq("CE"), np.abs(180 - raw), raw)
        out[f"{label}_preference_alignment"] = np.cos(
            np.deg2rad(out[f"{label}_preference_error_deg"])
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        out["internal_environment_ratio"] = out.PV_grad_eddy_mag / out.PV_grad_mag
    return out


def eddy_equal_scorecard(df, min_tilt_km=5.0):
    """Eddy-equal method comparison for the two candidate gradient vectors."""
    use = df[df.TiltDis.ge(min_tilt_km)].copy()
    eddy = (use.groupby(["Cyc", "method", "pv_surface_method", "ellipse_frac", "Eddy"], observed=True)
            .agg(env_error=("environment_preference_error_deg", "median"),
                 env_alignment=("environment_preference_alignment", "mean"),
                 full_error=("full_preference_error_deg", "median"),
                 full_alignment=("full_preference_alignment", "mean"),
                 env_mag=("PV_grad_mag", "median"),
                 full_mag=("PV_grad_full_mag", "median"),
                 coherence=("PV_grad_coherence", "median"),
                 effective_n=("PV_weight_effective_n", "median"))
            .reset_index())
    return (eddy.groupby(["Cyc", "method", "pv_surface_method", "ellipse_frac"], observed=True)
            .agg(eddies=("Eddy", "nunique"),
                 environmental_error=("env_error", "median"),
                 environmental_alignment=("env_alignment", "mean"),
                 environmental_within_45=("env_error", lambda x: np.mean(x <= 45)),
                 full_error=("full_error", "median"), full_alignment=("full_alignment", "mean"),
                 environmental_magnitude=("env_mag", "median"),
                 full_magnitude=("full_mag", "median"),
                 coherence=("coherence", "median"), effective_n=("effective_n", "median"))
            .reset_index())


def paired_to(df, reference="esp_gaussian_1"):
    """Attach a reference run to each matched Eddy-Day method row."""
    columns = ["PV", "PV_grad_mag", "PV_grad_theta", "PV_grad_topo_mag",
               "PV_grad_coherence", "PV_grad_full_mag", "PV_grad_full_theta"]
    ref = df[df.method.eq(reference)][["Eddy", "Day", *columns]].rename(
        columns={column: f"reference_{column}" for column in columns}
    )
    return df.merge(ref, on=["Eddy", "Day"], how="inner", validate="many_to_one")


def environmental_grid_gradients(grid):
    """Return bathymetric slopes and beta in true east/north coordinates."""
    import seacofs_tilt_tools as tilt

    dhdx, dhdy = tilt.phys_grad(
        grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
    )
    dfdx, dfdy = tilt.phys_grad(
        grid.f, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
    )
    return {
        "dh_e": np.cos(grid.angle) * dhdx - np.sin(grid.angle) * dhdy,
        "dh_n": np.sin(grid.angle) * dhdx + np.cos(grid.angle) * dhdy,
        "beta": np.sin(grid.angle) * dfdx + np.cos(grid.angle) * dfdy,
    }


def local_esp_fields(row, grid, frac=1.0, grid_gradients=None):
    """Return Laplacian-derived local fields for one eddy snapshot."""
    import seacofs_tilt_tools as tilt

    ii, jj = tilt.core_grid_indices(row, grid, frac=frac)
    dx = grid.x_grid[ii] - float(row.xc)
    dy = grid.y_grid[jj] - float(row.yc)
    reconstruction = tilt._esp_gaussian_vorticity(
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


def select_gradient_example_candidates(
    df,
    *,
    n_per_group=5,
    strong_quantile=0.90,
    seamount_max_coherence=0.35,
    seamount_min_depth_m=1500.0,
    slope_min_coherence=0.70,
    slope_depth_range_m=(200.0, 3500.0),
):
    """Select editable seamount-like and continental-slope candidate pools.

    These are screening rules, not an automatic geomorphic classification.
    Seamount-like cases combine strong local exposure with cancellation in
    deep water. Slope cases combine strong exposure with a coherent direction
    at intermediate depth. At most one day is retained from each eddy.
    """
    use = df.copy()
    if "method" in use:
        use = use[use.method.eq("esp_gaussian_1")].copy()
    required = {
        "Cyc", "Eddy", "Day", "h", "PV_grad_topo_p90_local_mag",
        "PV_grad_topo_mean_local_mag", "PV_grad_topo_coherence",
    }
    missing = required - set(use.columns)
    if missing:
        raise ValueError(f"Candidate selection requires columns: {sorted(missing)}")
    finite = np.isfinite(use[list(required - {"Cyc"})]).all(axis=1)
    use = use[finite & use.Cyc.isin(["AE", "CE"])].copy()
    use["strong_threshold"] = use.groupby("Cyc")[
        "PV_grad_topo_p90_local_mag"
    ].transform(lambda values: values.quantile(strong_quantile))
    strong = use.PV_grad_topo_p90_local_mag.ge(use.strong_threshold)

    seamount = use[
        use.PV_grad_topo_coherence.le(seamount_max_coherence)
        & use.h.ge(seamount_min_depth_m)
    ].copy()
    seamount["candidate_type"] = "seamount"
    seamount["strong_exposure"] = strong.reindex(seamount.index)
    seamount["candidate_score"] = (
        seamount.PV_grad_topo_p90_local_mag
        * (1.0 - seamount.PV_grad_topo_coherence.clip(0.0, 1.0))
    )

    low_depth, high_depth = map(float, slope_depth_range_m)
    slope = use[
        use.PV_grad_topo_coherence.ge(slope_min_coherence)
        & use.h.between(low_depth, high_depth)
    ].copy()
    slope["candidate_type"] = "slope"
    slope["strong_exposure"] = strong.reindex(slope.index)
    slope["candidate_score"] = (
        slope.PV_grad_topo_mean_local_mag
        * slope.PV_grad_topo_coherence.clip(0.0, 1.0)
    )

    pools = []
    for pool in (seamount, slope):
        pool = pool.sort_values(
            ["strong_exposure", "candidate_score"], ascending=[False, False]
        )
        pool = pool.drop_duplicates(["Cyc", "Eddy"])
        pool = pool.groupby("Cyc", group_keys=False).head(int(n_per_group))
        pool["rank"] = pool.groupby("Cyc").cumcount() + 1
        pools.append(pool)
    if not pools:
        return pd.DataFrame()
    return (pd.concat(pools, ignore_index=True)
            .sort_values(["candidate_type", "Cyc", "rank"])
            .reset_index(drop=True))


def plot_gradient_examples(
    rows,
    grid,
    *,
    vertical=None,
    min_pad_km=40.0,
    pad_rc=1.4,
    max_local_arrows=120,
    reference_percentile=90.0,
    arrow_length_inches=0.38,
    figsize=None,
):
    """Plot local and net PV-gradient vectors with one scale across panels."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    rows = rows.reset_index(drop=True)
    if rows.empty:
        raise ValueError("No candidate rows were supplied")
    gradients = environmental_grid_gradients(grid)
    prepared = []
    all_magnitudes = []
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
        prepared.append((row, local, pad, bathy))
        all_magnitudes.append(local.environment_mag.to_numpy(float))
        all_magnitudes.append(np.array([np.hypot(row.PV_grad_x, row.PV_grad_y)]))
        all_depths.append(bathy[np.isfinite(bathy)])

    magnitudes = np.concatenate(all_magnitudes)
    magnitudes = magnitudes[np.isfinite(magnitudes) & (magnitudes > 0)]
    if not len(magnitudes):
        raise ValueError("Selected examples contain no finite PV-gradient vectors")
    reference = float(np.nanpercentile(magnitudes, reference_percentile))
    quiver_scale = reference / float(arrow_length_inches)
    depths = np.concatenate([values for values in all_depths if len(values)])
    depth_norm = Normalize(vmin=float(np.nanmin(depths)), vmax=float(np.nanmax(depths)))

    n = len(rows)
    figsize = figsize or (5.8 * n, 5.3)
    fig, axes = plt.subplots(1, n, figsize=figsize, squeeze=False,
                             constrained_layout=True)
    axes = axes[0]
    bathy_artist = None
    local_quiver = None
    for ax, (row, local, pad, bathy) in zip(axes, prepared):
        bathy_artist = ax.pcolormesh(
            grid.X_grid, grid.Y_grid, bathy, cmap="Greys_r",
            norm=depth_norm, shading="auto",
        )
        ax.contour(grid.X_grid, grid.Y_grid, bathy, colors="0.35",
                   levels=8, linewidths=0.45, alpha=0.55)
        import seacofs_tilt_tools as tilt
        tilt.plot_ellipse(ax, row, grid, frac=1, color="red", lw=2.0, zorder=8)
        step = max(1, int(np.ceil(len(local) / max_local_arrows)))
        arrows = local.iloc[::step]
        valid = np.isfinite(arrows.environment_east) & np.isfinite(arrows.environment_north)
        local_quiver = ax.quiver(
            arrows.x[valid], arrows.y[valid],
            arrows.environment_east[valid], arrows.environment_north[valid],
            color="red", alpha=0.68, angles="xy", scale_units="inches",
            scale=quiver_scale, width=0.004, zorder=7,
        )
        ax.quiver(
            [row.xc], [row.yc], [row.PV_grad_x], [row.PV_grad_y],
            color="magenta", angles="xy", scale_units="inches",
            scale=quiver_scale, width=0.010, zorder=10,
        )
        ax.scatter(row.xc, row.yc, c="magenta", s=24, zorder=11)
        if vertical is not None:
            spine = vertical[
                vertical.Eddy.eq(row.Eddy) & vertical.Day.eq(row.Day)
            ].sort_values("Depth")
            ax.plot(spine.xc, spine.yc, color="limegreen", lw=2.0, zorder=9)
        local_mean = float(row.PV_grad_mean_local_mag)
        net = float(row.PV_grad_mag)
        ax.set(
            xlim=(row.xc-pad, row.xc+pad), ylim=(row.yc-pad, row.yc+pad),
            aspect="equal", xlabel="x (km)", ylabel="y (km)",
            title=(f"{row.Cyc}{int(row.Eddy)}, day {int(row.Day)}\n"
                   f"net/local = {net/local_mean:.2f}" if local_mean > 0
                   else f"{row.Cyc}{int(row.Eddy)}, day {int(row.Day)}"),
        )
    fig.colorbar(bathy_artist, ax=axes.tolist(), label="Water depth (km)", shrink=0.82)
    axes[0].quiverkey(
        local_quiver, 0.02, 1.08, reference,
        rf"${reference:.1e}\ \mathrm{{m^{{-2}}\,s^{{-1}}}}$",
        coordinates="axes", labelpos="E",
    )
    return fig, axes, reference
