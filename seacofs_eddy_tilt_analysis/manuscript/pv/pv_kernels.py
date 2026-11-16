"""Local scientific kernels adapted from legacy snapshot bb1d6c6.
Surface ESP geometry and amplitude are used together. No legacy runtime imports.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

def distance_km(lat1, lon1, lat2, lon2, earth_radius_km: float = 6357.0):
    """Great-circle distance in kilometres."""

    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return earth_radius_km * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

def phys_grad(F, X, Y, mask=None):
    """Gradient of ``F`` with respect to physical coordinates ``X`` and ``Y``.

    ``X`` and ``Y`` must be 2-D coordinate arrays in the same units. The
    returned gradients are in ``F`` per ``X``/``Y`` unit. A ROMS-style ocean
    mask can be supplied where 1 is ocean and 0 is land.
    """

    x_i, x_j = np.gradient(X)
    y_i, y_j = np.gradient(Y)
    F_i, F_j = np.gradient(F)
    jacobian = x_i * y_j - x_j * y_i
    dFdx = (F_i * y_j - F_j * y_i) / jacobian
    dFdy = (-F_i * x_j + F_j * x_i) / jacobian

    bad = np.isclose(jacobian, 0)
    dFdx[bad] = np.nan
    dFdy[bad] = np.nan

    if mask is not None:
        ocean = mask.astype(bool)
        dFdx = np.where(ocean, dFdx, np.nan)
        dFdy = np.where(ocean, dFdy, np.nan)
    return dFdx, dFdy

def bearing_from_xy(x, y):
    """Bearing in degrees, where 0 is north and 90 is east."""

    return (90.0 - np.degrees(np.arctan2(y, x))) % 360.0

def _validate_ellipse_fractions(frac, inner_frac=0.0):
    """Validate linear ellipse scale factors for a filled core or annulus."""
    frac = float(frac)
    inner_frac = float(inner_frac)
    if not np.isfinite(frac) or frac <= 0:
        raise ValueError("frac must be positive and finite")
    if not np.isfinite(inner_frac) or inner_frac < 0:
        raise ValueError("inner_frac must be non-negative and finite")
    if inner_frac >= frac:
        raise ValueError("inner_frac must be smaller than frac")
    return frac, inner_frac

def _weighted_percentile(values, weights, percentile):
    """Weighted percentile for finite one-dimensional samples."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return np.nan
    values, weights = values[valid], weights[valid]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) - 0.5 * weights
    cumulative /= weights.sum()
    return float(np.interp(float(percentile) / 100.0, cumulative, values))

def _esp_gaussian_vorticity(w, Rc, q11, q12, q22, dx, dy):
    """Reconstruct vorticity and its gradient from the Gaussian streamfunction.

    ``w`` is the fitted central relative vorticity, not the streamfunction
    amplitude. Coordinates and ``Rc`` are in kilometres. Returned vorticity
    gradients are per metre, matching the bathymetric and Coriolis gradients.
    The positive Gaussian envelope is returned separately for spatial
    weighting; the Laplacian-derived vorticity shape may change sign.
    """
    Rc = float(Rc)
    q11, q12, q22 = map(float, (q11, q12, q22))
    dx = np.asarray(dx, dtype=float)
    dy = np.asarray(dy, dtype=float)
    trace_q = q11 + q22
    if not np.isfinite(Rc) or Rc <= 0:
        raise ValueError("Rc must be positive and finite")
    if not np.isfinite(trace_q) or trace_q <= 0:
        raise ValueError("The ESP shape matrix must have positive finite trace")

    qr_x = q11 * dx + q12 * dy
    qr_y = q12 * dx + q22 * dy
    rho2 = q11 * dx**2 + 2.0 * q12 * dx * dy + q22 * dy**2
    gaussian = np.exp(-rho2 / Rc**2)

    # r.T @ Q**2 @ r == |Q @ r|**2 for symmetric Q.
    r_q2_r = qr_x**2 + qr_y**2
    laplacian_shape = 1.0 - 2.0 * r_q2_r / (Rc**2 * trace_q)
    zeta = float(w) * gaussian * laplacian_shape

    q2r_x = q11 * qr_x + q12 * qr_y
    q2r_y = q12 * qr_x + q22 * qr_y
    common = float(w) * gaussian
    dzeta_dx_km = common * (
        -2.0 * laplacian_shape * qr_x / Rc**2
        -4.0 * q2r_x / (Rc**2 * trace_q)
    )
    dzeta_dy_km = common * (
        -2.0 * laplacian_shape * qr_y / Rc**2
        -4.0 * q2r_y / (Rc**2 * trace_q)
    )
    return {
        "rho2": rho2,
        "gaussian": gaussian,
        "laplacian_shape": laplacian_shape,
        "zeta": zeta,
        "dzeta_dx": dzeta_dx_km / 1000.0,
        "dzeta_dy": dzeta_dy_km / 1000.0,
    }

def _surface_pv_footprint_statistics(
    df, grid, *, frac, inner_frac, surface_method="uniform",
    dh_dE=None, dh_dN=None, beta=None,
):
    """Average completed local PV-gradient fields within surface footprints.

    Existing ``*_mag`` columns are the magnitude of the spatially averaged
    vector. ``*_mean_local_mag`` and ``*_rms_local_mag`` quantify exposure
    without cancellation. ``*_coherence`` is net magnitude divided by mean
    local magnitude and lies in [0, 1]. ``surface_method='esp_gaussian'``
    reconstructs relative vorticity as the Laplacian of the fitted Gaussian
    streamfunction, with ``w`` as its central value. The positive Gaussian
    envelope remains the spatial weight.
    """
    if dh_dE is None or dh_dN is None or beta is None:
        dhdx, dhdy = phys_grad(
            grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
        )
        dh_dN = np.sin(grid.angle) * dhdx + np.cos(grid.angle) * dhdy
        dh_dE = np.cos(grid.angle) * dhdx - np.sin(grid.angle) * dhdy
        dfdx, dfdy = phys_grad(
            grid.f, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho
        )
        beta = np.sin(grid.angle) * dfdx + np.cos(grid.angle) * dfdy
    h = np.asarray(grid.h, dtype=float)
    f = np.asarray(grid.f, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        fields = {
            "h": h, "f": f, "beta": beta, "dhdx": dh_dE, "dhdy": dh_dN,
            "inv_h": 1.0 / h, "f_over_h": f / h,
            "beta_over_h": beta / h,
            "dhdx_over_h2": dh_dE / h**2,
            "dhdy_over_h2": dh_dN / h**2,
            "f_dhdx_over_h2": f * dh_dE / h**2,
            "f_dhdy_over_h2": f * dh_dN / h**2,
        }
    names = list(fields)
    stack = np.stack([
        np.where(grid.mask_rho, np.asarray(field, dtype=float), np.nan)
        for field in fields.values()
    ])
    records = []
    for row in df.itertuples(index=False):
        ii, jj = core_grid_indices(row, grid, frac=frac, inner_frac=inner_frac)
        record = {"PV_footprint_n": 0, "PV_weight_effective_n": 0.0}
        if not len(ii):
            records.append(record)
            continue
        sample = stack[:, ii, jj]
        if surface_method == "esp_gaussian":
            dx = grid.x_grid[ii] - float(row.xc)
            dy = grid.y_grid[jj] - float(row.yc)
            reconstruction = _esp_gaussian_vorticity(
                row.w, row.Rc, row.q11, row.q12, row.q22, dx, dy
            )
            zeta = reconstruction["zeta"]
            weights = reconstruction["gaussian"]
        else:
            zeta = np.full(len(ii), float(row.w))
            weights = np.ones(len(ii), dtype=float)

        def weighted_mean(values):
            values = np.asarray(values, dtype=float)
            valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
            return float(np.average(values[valid], weights=weights[valid])) if valid.any() else np.nan

        mean = {name: weighted_mean(sample[pos]) for pos, name in enumerate(names)}
        record.update({name: mean[name] for name in ("h", "f", "beta", "dhdx", "dhdy")})
        plan_x = np.zeros(sample.shape[1], dtype=float)
        plan_y = sample[names.index("beta_over_h")]
        topo_x = -(zeta * sample[names.index("dhdx_over_h2")] + sample[names.index("f_dhdx_over_h2")])
        topo_y = -(zeta * sample[names.index("dhdy_over_h2")] + sample[names.index("f_dhdy_over_h2")])
        if surface_method == "esp_gaussian":
            eddy_x = reconstruction["dzeta_dx"] * sample[names.index("inv_h")]
            eddy_y = reconstruction["dzeta_dy"] * sample[names.index("inv_h")]
        else:
            eddy_x = np.zeros(sample.shape[1], dtype=float)
            eddy_y = np.zeros(sample.shape[1], dtype=float)
        vectors = {
            "PV_grad_plan": (plan_x, plan_y),
            "PV_grad_topo": (topo_x, topo_y),
            "PV_grad": (plan_x + topo_x, plan_y + topo_y),
            "PV_grad_eddy": (eddy_x, eddy_y),
            "PV_grad_full": (plan_x + topo_x + eddy_x, plan_y + topo_y + eddy_y),
        }
        record["zeta_mean"] = weighted_mean(zeta)
        record["abs_vort"] = weighted_mean(zeta + sample[names.index("f")])
        record["PV"] = weighted_mean(
            zeta * sample[names.index("inv_h")] + sample[names.index("f_over_h")]
        )
        positive_weights = weights[np.isfinite(weights) & (weights > 0)]
        record["PV_weight_effective_n"] = float(
            positive_weights.sum() ** 2 / np.sum(positive_weights**2)
        ) if len(positive_weights) else 0.0
        for prefix, (east, north) in vectors.items():
            valid = np.isfinite(east) & np.isfinite(north) & np.isfinite(weights) & (weights > 0)
            if not valid.any():
                continue
            east_mean = float(np.average(east[valid], weights=weights[valid]))
            north_mean = float(np.average(north[valid], weights=weights[valid]))
            local_mag = np.hypot(east[valid], north[valid])
            net_mag = float(np.hypot(east_mean, north_mean))
            mean_local_mag = float(np.average(local_mag, weights=weights[valid]))
            record[f"{prefix}_x"] = east_mean
            record[f"{prefix}_y"] = north_mean
            record[f"{prefix}_mean_local_mag"] = mean_local_mag
            record[f"{prefix}_rms_local_mag"] = float(
                np.sqrt(np.average(local_mag**2, weights=weights[valid]))
            )
            record[f"{prefix}_p90_local_mag"] = _weighted_percentile(
                local_mag, weights[valid], 90
            )
            record[f"{prefix}_coherence"] = net_mag / mean_local_mag if mean_local_mag > 0 else np.nan
            if prefix == "PV_grad":
                record["PV_footprint_n"] = int(valid.sum())
        records.append(record)
    result = pd.DataFrame(records, index=df.index)
    expected = [
        "h", "f", "beta", "dhdx", "dhdy", "zeta_mean", "abs_vort", "PV",
        "PV_footprint_n", "PV_weight_effective_n",
    ]
    for prefix in ("PV_grad_plan", "PV_grad_topo", "PV_grad", "PV_grad_eddy", "PV_grad_full"):
        expected.extend([
            f"{prefix}_x", f"{prefix}_y", f"{prefix}_mean_local_mag",
            f"{prefix}_rms_local_mag", f"{prefix}_p90_local_mag",
            f"{prefix}_coherence",
        ])
    return result.reindex(columns=expected)

def _select_profile_peak_w(surface, vertical, max_depth_m):
    """Select a signed profile amplitude, preserving surface row order/index."""
    keys = ["Eddy", "Day"]
    for table, required in ((surface, keys + ["w"]),
                            (vertical, keys + ["Depth", "w"])):
        missing = set(required) - set(table.columns)
        if missing:
            raise ValueError(f"Peak-w selection requires columns: {sorted(missing)}")
    candidates = vertical[keys + ["Depth", "w"]].copy()
    candidates = candidates.loc[
        candidates[keys].notna().all(axis=1)
        & np.isfinite(candidates["Depth"]) & np.isfinite(candidates["w"])
        & candidates["Depth"].between(0, max_depth_m)
    ]
    candidates["magnitude"] = candidates["w"].abs()
    counts = candidates.groupby(keys).size()
    peaks = (candidates.sort_values(["magnitude", "Depth"], ascending=[False, True],
                                   kind="stable")
             .drop_duplicates(keys).set_index(keys))
    lookup = pd.MultiIndex.from_frame(surface[keys])
    out = surface.copy()
    if "w_surface" not in out:
        out["w_surface"] = out["w"]
    out["w"] = peaks["w"].reindex(lookup).to_numpy()
    out["w_selected_depth_m"] = peaks["Depth"].reindex(lookup).to_numpy()
    out["w_profile_n"] = counts.reindex(lookup, fill_value=0).to_numpy()
    return out

def core_grid_indices(
    row, grid: Grid, circle_region_flag: bool = False,
    frac=1, inner_frac=0.0,
):
    """Return ocean-grid indices inside a scaled core ellipse or annulus."""

    frac, inner_frac = _validate_ellipse_fractions(frac, inner_frac)

    if circle_region_flag:
        if not (hasattr(row, "rmax") and np.isfinite(row.rmax) and row.rmax > 0):
            return np.array([], dtype=int), np.array([], dtype=int)
        q = np.eye(2)
        base_threshold = float(row.rmax) ** 2
    else:
        if hasattr(row, "q11") and np.isfinite(row.q11):
            q = np.array([[row.q11, row.q12], [row.q12, row.q22]], dtype=float)
        elif hasattr(row, "Q"):
            q = np.asarray(row.Q, dtype=float)
        else:
            return np.array([], dtype=int), np.array([], dtype=int)
        if q.shape != (2, 2) or not np.isfinite(q).all() or not np.isfinite(row.Rc) or row.Rc <= 0:
            return np.array([], dtype=int), np.array([], dtype=int)
        base_threshold = float(row.Rc) ** 2 / 2.0
    threshold = base_threshold * frac**2
    inner_threshold = base_threshold * inner_frac**2
    eigenvalues = np.linalg.eigvalsh(q)
    if not np.isfinite(eigenvalues).all() or eigenvalues.min() <= 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    radius = np.sqrt(threshold / eigenvalues.min())
    i0 = max(0, int(np.searchsorted(grid.x_grid, row.xc - radius, side="left")))
    i1 = min(len(grid.x_grid), int(np.searchsorted(grid.x_grid, row.xc + radius, side="right")))
    j0 = max(0, int(np.searchsorted(grid.y_grid, row.yc - radius, side="left")))
    j1 = min(len(grid.y_grid), int(np.searchsorted(grid.y_grid, row.yc + radius, side="right")))
    if i0 >= i1 or j0 >= j1:
        return np.array([], dtype=int), np.array([], dtype=int)
    ii, jj = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1), indexing="ij")
    dx = grid.x_grid[ii] - float(row.xc)
    dy = grid.y_grid[jj] - float(row.yc)
    rho2 = q[0, 0] * dx**2 + 2.0 * q[0, 1] * dx * dy + q[1, 1] * dy**2
    use = (
        (rho2 <= threshold) & (rho2 >= inner_threshold)
        & grid.mask_rho[ii, jj].astype(bool)
    )
    return ii[use].astype(int), jj[use].astype(int)
