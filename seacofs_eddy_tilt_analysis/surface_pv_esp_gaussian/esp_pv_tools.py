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


def local_esp_fields(row, grid, frac=1.0):
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
    dhdx, dhdy = tilt.phys_grad(grid.h, grid.X_grid*1e3, grid.Y_grid*1e3, grid.mask_rho)
    dh_n = np.sin(grid.angle)*dhdx + np.cos(grid.angle)*dhdy
    dh_e = np.cos(grid.angle)*dhdx - np.sin(grid.angle)*dhdy
    dfdx, dfdy = tilt.phys_grad(grid.f, grid.X_grid*1e3, grid.Y_grid*1e3, grid.mask_rho)
    beta = np.sin(grid.angle)*dfdx + np.cos(grid.angle)*dfdy
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
