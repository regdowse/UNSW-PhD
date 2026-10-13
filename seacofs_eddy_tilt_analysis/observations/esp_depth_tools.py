"""Matched ADCP depth applications of the ESP SOLO and LATTE methods.

Inner estimators default to the user's ESP backend. An explicit SOLO option
uses the scaled local crossing initialisation from the ESP revision experiments.
The Gaussian outer fit uses the backend objective/sign filter, but reports failed
optimisation instead of silently returning its seed. Positions/radii are metres.
"""
from pathlib import Path
import hashlib
import importlib.util
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit, OptimizeWarning


def load_backend(path):
    path = Path(path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location("observation_esp_backend", path)
    backend = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(backend)
    return backend, {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def load_observations(path, start, end, origin, min_attendance=None):
    """Keep original depth bins; no time averaging, spatial interpolation or QC inference."""
    raw = pd.read_csv(path)
    raw["time"] = pd.to_datetime(raw["Start (UTC)"], utc=True)
    selected = raw.loc[raw.time.between(pd.Timestamp(start), pd.Timestamp(end))].copy()
    rename = {"Mean Lon": "lon", "Mean Lat": "lat", "Depth": "depth",
              "U (m/sec)": "u", "V (m/sec)": "v"}
    selected = selected.rename(columns=rename)
    for c in rename.values():
        selected[c] = pd.to_numeric(selected[c], errors="coerce")
    selected["attendance"] = pd.to_numeric(selected["Bin Attendance %"], errors="coerce")
    finite = np.isfinite(selected[["lon", "lat", "depth", "u", "v"]]).all(axis=1)
    keep = finite.copy()
    if min_attendance is not None:
        keep &= selected.attendance.ge(min_attendance)
    data = selected.loc[keep].sort_values(["depth", "time"]).copy()
    if data.empty:
        raise ValueError("No observations pass the selected time window and filters.")
    if data.duplicated(["time", "depth"]).any():
        raise ValueError("Duplicate time/depth observations: select the instrument before fitting.")
    lon0, lat0 = origin
    data["x"] = (data.lon - lon0) * 111000 * np.cos(np.deg2rad(lat0))
    data["y"] = (data.lat - lat0) * 111000
    coverage = selected.groupby("depth").size().rename("n_selected").to_frame()
    coverage["n_used"] = data.groupby("depth").size().reindex(coverage.index, fill_value=0)
    coverage["n_removed"] = coverage.n_selected - coverage.n_used
    audit = {"selected_rows": len(selected), "used_rows": len(data),
             "unique_times": data.time.nunique(), "start_utc": str(data.time.min()),
             "end_utc": str(data.time.max()), "depth_range_m": [data.depth.min(), data.depth.max()],
             "missing_avg_qc": int(selected["Avg QC"].isna().sum()),
             "zero_attendance_rows": int(selected.attendance.eq(0).sum()),
             "min_attendance": min_attendance, "origin_lon_lat": origin}
    return data, coverage, audit


def _rho(x, y, centre, q):
    r = np.column_stack([x - centre[0], y - centre[1]])
    return np.sqrt(np.maximum(np.einsum("ni,ij,nj->n", r, q, r), 0))


def _ellipse(q):
    if not np.isfinite(q).all() or not np.allclose(q, q.T):
        raise ValueError("invalid_ellipse")
    values, vectors = np.linalg.eigh(q)
    if values.min() <= 0 or not np.isclose(np.linalg.det(q), 1, rtol=1e-4):
        raise ValueError("invalid_ellipse")
    major = vectors[:, 0]
    return np.sqrt(values[-1] / values[0]), np.degrees(np.arctan2(major[1], major[0])) % 180


def _spread(x, y):
    xy = np.column_stack([x, y])
    s = np.linalg.svd(xy - xy.mean(axis=0), compute_uv=False)
    return float(s[1] / s[0]) if len(s) > 1 and s[0] > 0 else 0.


def _solo_local(l, vt, vn, centre_guess, search_half_width, core_thresh):
    """SOLO equations with a scaled, bounded crossing fit for short transects.

    Same local-initialisation variant used in UNSW-PhD/esp/local_estimators.py.
    The zero crossing must be observed and lie inside the search window.
    """
    l, vt, vn = [np.asarray(a, float) for a in (l, vt, vn)]
    order = np.argsort(l)
    l, vt, vn = [a[order] for a in (l, vt, vn)]
    mask = abs(l-centre_guess) <= search_half_width
    x, y = l[mask], vn[mask]
    if len(x) < 4:
        raise ValueError("insufficient_search")
    crossings = np.flatnonzero((y[:-1]*y[1:] <= 0) & (y[:-1] != y[1:]))
    if not len(crossings):
        raise ValueError("no_supported_solo_crossing")
    roots = x[crossings] - y[crossings]*np.diff(x)[crossings]/np.diff(y)[crossings]
    guess = roots[np.argmin(abs(roots-centre_guess))]
    scale = search_half_width

    def fit(mask, guess):
        if mask.sum() < 4 or np.unique(l[mask]).size < 4:
            raise ValueError("insufficient_core")
        xx, yy = (l[mask]-centre_guess)/scale, vn[mask]
        lower, upper = max(xx.min(), -1.), min(xx.max(), 1.)
        if upper <= lower:
            raise ValueError("no_search_overlap")
        z0 = np.clip((guess-centre_guess)/scale, lower+1e-10, upper-1e-10)
        dz = xx-z0
        c, d = np.linalg.lstsq(np.c_[dz, dz**3], yy, rcond=None)[0]
        def model(z, root, c, d):
            delta = z-root
            return c*delta+d*delta**3
        p, _ = curve_fit(model, xx, yy, p0=[z0,c,d],
                         bounds=([lower,-np.inf,-np.inf],[upper,np.inf,np.inf]), maxfev=10000)
        if min(p[0]-lower, upper-p[0]) <= 1e-7:
            raise ValueError("solo_crossing_at_bound")
        return centre_guess+scale*p[0], p[1]/scale

    initial, _ = fit(mask, guess)
    core = abs(l-initial) <= core_thresh
    l0, c = fit(core, initial)
    if not np.isfinite(c) or abs(c) < np.finfo(float).eps:
        raise ValueError("solo_zero_gradient")
    delta = (l[core]-l0)/scale
    a, _ = np.linalg.lstsq(np.c_[np.ones(core.sum()), delta**2], vt[core], rcond=None)[0]
    return l0, a/c, 2*c, np.eye(2), c


def _outer(backend, x, y, u, v, centre, q, omega, radius_limit, min_points):
    rho = _rho(x, y, centre, q)
    r = np.column_stack([x - centre[0], y - centre[1]])
    qr = np.linalg.norm(r @ q.T, axis=1)
    vt = backend.tangential_velocity(x, y, u, v, *centre, q)
    use = np.isfinite(vt) & (vt * omega >= 0) & (qr > 0)
    rho, qr, vt = rho[use], qr[use], vt[use]
    record = {"n_outer_retained": int(use.sum())}
    if len(rho) < min_points or np.unique(rho).size < 3:
        return record, "insufficient_outer"
    target = vt * rho / qr
    seed = max(float(rho[np.argmax(abs(target))] * np.sqrt(2)), 1e-6)

    def gaussian(rho2, om, rc):
        return om * np.sqrt(rho2) * np.exp(-rho2 / rc**2)

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", OptimizeWarning)
            p, _ = curve_fit(gaussian, rho**2, target, p0=[omega, seed],
                             bounds=([-np.inf, 1e-8], [np.inf, np.inf]), maxfev=10000)
        om, rc = p
        record.update(Omega=om, Rc_m=rc, psi0=-.5 * om * rc**2,
                      outer_warning="; ".join(str(w.message) for w in caught),
                      rho_min_m=rho.min(), rho_max_m=rho.max(),
                      peak_bracketed=bool(rho.min() < rc / np.sqrt(2) < rho.max()))
        if not np.isfinite(p).all() or rc <= 0:
            return record, "nonfinite_outer"
        if rc > radius_limit:
            return record, "outer_radius_limit"
        if om * omega <= 0:
            return record, "outer_rotation_reversal"
        pred = gaussian(rho**2, om, rc)
        record["outer_rmse_m_s"] = float(np.sqrt(np.mean((target - pred)**2)))
        return record, "ok"
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        record["outer_warning"] = str(exc)
        return record, "outer_optimisation_failed"


def fit_depths(data, backend, method, *, core_radius=35000., outer_radius=100000.,
               search_half_width=60000., max_centre_step=35000., radius_limit=100000.,
               min_core_points=8, min_outer_points=12, expected_sign=1,
               require_peak_bracket=True, depth_levels=None, solo_initialisation="backend"):
    """Follow accepted inner centres downwards; failed depths do not update the seed.

    Both methods start at the specified geographic origin (x=y=0). SOLO uses a
    fixed transect direction derived from the shallowest selected bin. LATTE
    repeats its core selection twice, using the last accepted ellipse as seed.
    Reconstruction acceptance is independent of inner-centre acceptance.
    """
    method = method.upper()
    if method not in ("LATTE", "SOLO"):
        raise ValueError("method must be LATTE or SOLO")
    if solo_initialisation not in ("backend", "local"):
        raise ValueError("solo_initialisation must be backend or local")
    if min(core_radius, outer_radius, search_half_width, max_centre_step, radius_limit) <= 0:
        raise ValueError("All radii and search distances must be positive.")
    if expected_sign not in (-1, 1):
        raise ValueError("expected_sign must be +1 (AE here) or -1 (CE here).")
    surface = data.loc[data.depth.eq(data.depth.min())]
    xy = surface[["x", "y"]].to_numpy()
    _, _, vh = np.linalg.svd(xy - xy.mean(axis=0), full_matrices=False)
    tangent = vh[0]
    if tangent @ (xy[-1] - xy[0]) < 0:
        tangent = -tangent
    normal = np.array([-tangent[1], tangent[0]])
    centre, q_previous = np.zeros(2), np.eye(2)
    rows, samples = [], {}
    levels = sorted(data.depth.unique() if depth_levels is None else depth_levels)
    for depth in levels:
        d = data.loc[data.depth.eq(depth)]
        x, y, u, v = [d[c].to_numpy(float) for c in ("x", "y", "u", "v")]
        rec = dict(method=method, depth_m=depth, n_obs=len(d), n_core=0,
                   inner_valid=False, valid=False, status="unattempted", warnings="",
                   xc_m=np.nan, yc_m=np.nan, q11=np.nan, q12=np.nan, q22=np.nan,
                   inner_w=np.nan, inner_Omega=np.nan, inner_R2=np.nan,
                   candidate_min_eigenvalue=np.nan, Omega=np.nan, w=np.nan,
                   Rc_m=np.nan, alpha=np.nan, angle_deg=np.nan, vector_R2=np.nan,
                   vector_rmse_m_s=np.nan, core_spread_ratio=np.nan,
                   n_outer_retained=0, peak_bracketed=False, centre_step_m=np.nan)
        sample = d.copy()
        sample["inner"] = False
        sample["outer"] = False
        samples[depth] = sample
        try:
            if len(d) < min_core_points:
                raise ValueError("insufficient_observations")
            if method == "LATTE":
                c, q = centre.copy(), q_previous.copy()
                for _ in range(2):
                    core = _rho(x, y, c, q) <= core_radius
                    rec["n_core"] = int(core.sum())
                    if core.sum() < min_core_points:
                        raise ValueError("insufficient_core")
                    rec["core_spread_ratio"] = _spread(x[core], y[core])
                    if rec["core_spread_ratio"] < 1e-6:
                        raise ValueError("collinear_latte_sampling")
                    xc, yc, w, q, om, inner_r2 = backend.latte(x[core], y[core], u[core], v[core])
                    rec["inner_R2"] = inner_r2
                    if np.isfinite(q).all():
                        rec["candidate_min_eigenvalue"] = np.linalg.eigvalsh(q).min()
                    c = np.array([xc, yc])
                    _ellipse(q)
                    if not np.isfinite(c).all() or np.linalg.norm(c - centre) > max_centre_step:
                        raise ValueError("centre_jump")
            else:
                # Use actual points (no interpolation) in a fixed local frame.
                # As in the reference, normal is a +90 degree rotation of tangent.
                xy = np.column_stack([x, y])
                along = xy @ tangent
                uv = np.column_stack([u, v])
                vt, vn = uv @ tangent, uv @ normal
                seed = centre @ tangent
                search = abs(along - seed) <= search_half_width
                if search.sum() < min_core_points:
                    raise ValueError("insufficient_search")
                order = np.argsort(along[search])
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", OptimizeWarning)
                    if solo_initialisation == "local":
                        l0, r0, w, q, om = _solo_local(
                            along[search][order], vt[search][order], vn[search][order],
                            seed, search_half_width, core_radius)
                    else:
                        l0, r0, w, q, om = backend.solo(
                            along[search][order], vt[search][order], vn[search][order],
                            core_thresh=core_radius)
                rec["warnings"] = "; ".join(str(a.message) for a in caught)
                if not np.isfinite([l0, r0, w, om]).all():
                    raise ValueError("nonfinite_inner")
                if not along[search].min() < l0 < along[search].max():
                    raise ValueError("unsupported_solo_crossing")
                # Cross-track coordinates are measured relative to the same line.
                # The line offset is fixed from the shallowest sampling track.
                offset = surface[["x", "y"]].to_numpy().mean(axis=0) @ normal
                c = l0 * tangent + (offset + r0) * normal
                core = search & (abs(along - l0) <= core_radius)
                rec["n_core"] = int(core.sum())
                rec["core_spread_ratio"] = _spread(x[core], y[core]) if core.sum() > 1 else 0.
                if core.sum() < min_core_points:
                    raise ValueError("insufficient_core")
                uh = -om * (y[core] - c[1])
                vv = om * (x[core] - c[0])
                rec["inner_R2"] = backend.vector_R2(u[core], v[core], uh, vv)
                rec["candidate_min_eigenvalue"] = 1.
            alpha, angle = _ellipse(q)
            if not np.isfinite([*c, w, om]).all() or om == 0:
                raise ValueError("nonfinite_inner")
            if w * expected_sign <= 0:
                raise ValueError("inner_rotation_reversal")
            step = np.linalg.norm(c - centre)
            if step > max_centre_step:
                raise ValueError("centre_jump")
            if rec["core_spread_ratio"] < .05:
                rec["warnings"] += "; narrow track: elliptical geometry weakly constrained" if method == "LATTE" else "; circular geometry imposed"
            rec.update(inner_valid=True, xc_m=c[0], yc_m=c[1], inner_w=w,
                       inner_Omega=om, q11=q[0, 0], q12=q[0, 1], q22=q[1, 1],
                       alpha=alpha, angle_deg=angle if alpha > 1.0001 else np.nan,
                       centre_step_m=step)
            centre, q_previous = c.copy(), q.copy()
            outer = _rho(x, y, c, q) <= outer_radius
            sample.loc[:, "inner"], sample.loc[:, "outer"] = core, outer
            out, status = _outer(backend, x[outer], y[outer], u[outer], v[outer],
                                 c, q, om, radius_limit, min_outer_points)
            rec.update(out)
            rec["status"] = status
            if status == "ok" and require_peak_bracket and not rec["peak_bracketed"]:
                rec["status"] = "peak_not_bracketed"
            if rec["status"] == "ok":
                rec["valid"] = True
                rec["w"] = rec["Omega"] * np.trace(q)
                uh, vv = backend.model_uv_at_xy(x[outer], y[outer], *c, q, rec["Omega"], rec["Rc_m"])
                rec["vector_R2"] = backend.vector_R2(u[outer], v[outer], uh, vv)
                rec["vector_rmse_m_s"] = np.sqrt(np.mean((u[outer]-uh)**2 + (v[outer]-vv)**2))
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            rec["status"] = str(exc)
        rows.append(rec)
    results = pd.DataFrame(rows)
    good = results.inner_valid
    results["offset_east_km"] = np.nan
    results["offset_north_km"] = np.nan
    if good.any():
        first = results.loc[good].iloc[0]
        results.loc[good, "offset_east_km"] = (results.loc[good, "xc_m"] - first.xc_m) / 1000
        results.loc[good, "offset_north_km"] = (results.loc[good, "yc_m"] - first.yc_m) / 1000
        results["reference_depth_m"] = first.depth_m
    results["displacement_km"] = np.hypot(results.offset_east_km, results.offset_north_km)
    return results, samples


def report(results):
    print(f"{results.method.iloc[0]}: {results.inner_valid.sum()}/{len(results)} accepted inner fits; "
          f"{results.valid.sum()}/{len(results)} accepted Gaussian reconstructions.")
    print(results.groupby("status").size().to_string())
    if results.inner_valid.any():
        inner = results.loc[results.inner_valid]
        print(f"Centre offsets reference the shallowest accepted fit: {inner.depth_m.min():g} m.")
    if results.valid.any():
        depths = results.loc[results.valid, "depth_m"]
        print(f"Reconstructed measured levels: {depths.min():g}–{depths.max():g} m (gaps retained).")
    print("The ADCP sampling limit is not an estimate of the eddy bottom. R² is in-sample.")


def plot_sampling(data, coverage, half_width_km=120):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
    d = data.loc[data.depth.eq(data.depth.min())]
    hours = (d.time - d.time.min()).dt.total_seconds() / 3600
    q = axes[0].quiver(d.x/1000, d.y/1000, d.u, d.v, hours, cmap="viridis", scale=12)
    axes[0].scatter(0, 0, marker="x", color="red", label="Initial centre")
    axes[0].set(xlabel="East (km)", ylabel="North (km)", title=f"Observed velocities at {d.depth.iloc[0]:g} m", aspect="equal")
    axes[0].legend(fontsize=8)
    fig.colorbar(q, ax=axes[0], label="Hours since transect start")
    axes[1].plot(coverage.n_selected, coverage.index, label="Selected")
    axes[1].plot(coverage.n_used, coverage.index, label="After filters")
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Observations", ylabel="Depth (m)", title="Depth coverage")
    axes[1].legend()
    axes[2].plot(hours, d.u, label="Eastward u")
    axes[2].plot(hours, d.v, label="Northward v")
    axes[2].axhline(0, color="grey", lw=.5)
    axes[2].set(xlabel="Hours since transect start", ylabel="Velocity (m/s)", title="Surface-bin time series")
    axes[2].legend()
    return fig


def plot_diagnostics(results):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True, layout="constrained")
    z = results.depth_m
    axes[0].plot(results.n_obs, z, label="All")
    axes[0].plot(results.n_core, z, label="Last core selection")
    axes[0].set(xlabel="Observation count", ylabel="Depth (m)")
    axes[0].legend(fontsize=8)
    axes[1].semilogx(results.core_spread_ratio.where(results.core_spread_ratio > 0), z, ".-")
    axes[1].axvline(.05, color="grey", ls="--", label="Narrow-track flag")
    axes[1].set_xlabel("Across/along spatial spread")
    axes[1].legend(fontsize=8)
    axes[2].plot(results.inner_R2, z, ".-", label="Attempted inner fit")
    axes[2].plot(results.vector_R2, z, ".-", label="Accepted full fit")
    axes[2].set_xlabel("Vector R²")
    axes[2].legend(fontsize=8)
    axes[0].invert_yaxis()
    fig.suptitle(f"{results.method.iloc[0]} · high R² alone does not establish valid ellipse geometry")
    return fig


def plot_profiles(results):
    if not results.inner_valid.any():
        return _no_fit(f"{results.method.iloc[0]} · depth profiles",
                       "No accepted inner-core fits.\nInspect the sampling and ellipse diagnostics above.")
    fig, axes = plt.subplots(1, 5, figsize=(15, 5), sharey=True, layout="constrained")
    z = results.depth_m
    axes[0].plot(results.offset_east_km, z, ".-", label="East")
    axes[0].plot(results.offset_north_km, z, ".-", label="North")
    axes[0].set(xlabel="Centre offset (km)", ylabel="Depth (m)")
    axes[0].legend()
    axes[1].plot(results.Rc_m.where(results.valid)/1000, z, ".-")
    axes[1].set_xlabel("Gaussian Rc (km)")
    axes[2].plot(results.alpha.where(results.inner_valid), z, ".-")
    axes[2].set_xlabel("Axis ratio")
    axes[3].plot(results.inner_w*1e5, z, ".-", label="Inner")
    axes[3].plot(results.w.where(results.valid)*1e5, z, ".-", label="Gaussian")
    axes[3].set_xlabel("Vorticity (10⁻⁵ s⁻¹)")
    axes[3].legend()
    axes[4].plot(results.vector_R2.where(results.valid), z, ".-")
    axes[4].set_xlabel("Vector R²")
    axes[0].invert_yaxis()
    for ax in axes:
        ax.grid(alpha=.2)
    ref = results.reference_depth_m.iloc[0] if "reference_depth_m" in results else np.nan
    fig.suptitle(f"{results.method.iloc[0]} · centre offsets from {ref:g} m · gaps are rejected fits")
    return fig


def _q(row):
    return np.array([[row.q11, row.q12], [row.q12, row.q22]])


def _contour(row):
    values, vectors = np.linalg.eigh(_q(row))
    angle = np.linspace(0, 2*np.pi, 150)
    ring = np.column_stack([np.cos(angle), np.sin(angle)])
    ring = (ring * (row.Rc_m / np.sqrt(2*values))) @ vectors.T
    return (ring + [row.xc_m, row.yc_m]) / 1000


def _no_fit(title, message="No accepted Gaussian reconstruction.\nInspect the depth-fit diagnostics and sampling."):
    fig, ax = plt.subplots(figsize=(8, 3), layout="constrained")
    ax.text(.5, .5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_title(title)
    ax.set_axis_off()
    return fig


def plot_maps(results, samples, backend, depths=(18, 100, 200, 300), half_width_km=120, speed_max=1.5):
    """Nearest measured levels to requests; never substitute a deeper valid fit."""
    actual = sorted(set(results.depth_m.iloc[np.argmin(abs(results.depth_m-z))] for z in depths))
    fig, axs = plt.subplots(1, len(actual), figsize=(4*len(actual), 4.5), squeeze=False, layout="constrained")
    axis = np.linspace(-half_width_km, half_width_km, 90)
    xx, yy = np.meshgrid(axis*1000, axis*1000)
    image = None
    for ax, z in zip(axs.ravel(), actual):
        row = results.loc[results.depth_m.eq(z)].iloc[0]
        d = samples[z]
        if row.valid:
            u, v = backend.model_uv_at_xy(xx, yy, row.xc_m, row.yc_m, _q(row), row.Omega, row.Rc_m)
            image = ax.pcolormesh(axis, axis, np.hypot(u, v), vmin=0, vmax=speed_max, cmap="Blues", shading="auto")
            ax.quiver(xx[::9,::9]/1000, yy[::9,::9]/1000, u[::9,::9], v[::9,::9], scale=12, color="grey", alpha=.6)
            ring = _contour(row)
            ax.plot(*ring.T, "r--", lw=1)
        else:
            ax.text(.02, .02, row.status.replace("_", " "), transform=ax.transAxes, fontsize=8, wrap=True)
        ax.quiver(d.x/1000, d.y/1000, d.u, d.v, scale=12, color="black", width=.005)
        if row.inner_valid:
            ax.scatter(row.xc_m/1000, row.yc_m/1000, marker="x", color="red")
        ax.set(xlim=(-half_width_km, half_width_km), ylim=(-half_width_km, half_width_km), aspect="equal", xlabel="East (km)", ylabel="North (km)", title=f"{z:g} m")
    if image is not None:
        fig.colorbar(image, ax=list(axs.ravel()), label="Reconstructed speed (m/s)", shrink=.8)
    fig.suptitle(f"{results.method.iloc[0]} · black: observed · grey: reconstructed · red: centre/core contour")
    return fig


def plot_3d(results, backend, half_width_km=120, speed_max=1.5):
    valid = results.loc[results.valid]
    if not results.inner_valid.any():
        return _no_fit(f"{results.method.iloc[0]} · 3D reconstruction")
    fig = plt.figure(figsize=(9, 7), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    norm = Normalize(results.depth_m.min(), results.depth_m.max())
    cmap = plt.colormaps["viridis"]
    selected = valid.iloc[np.unique(np.linspace(0, len(valid)-1, min(9, len(valid))).astype(int))] if len(valid) else valid
    axis = np.linspace(-half_width_km, half_width_km, 11)
    xx, yy = np.meshgrid(axis*1000, axis*1000)
    for _, row in selected.iterrows():
        ring = _contour(row)
        color = cmap(norm(row.depth_m))
        ax.plot(ring[:,0], ring[:,1], np.full(len(ring), row.depth_m), color=color, lw=2)
        u, v = backend.model_uv_at_xy(xx, yy, row.xc_m, row.yc_m, _q(row), row.Omega, row.Rc_m)
        # Arrow length is a common visual scale: 12 km per 1 m/s.
        inside = _rho(xx.ravel(), yy.ravel(), (row.xc_m, row.yc_m), _q(row)).reshape(xx.shape) <= min(outer_extent(row), half_width_km*1000)
        ax.quiver(xx[inside]/1000, yy[inside]/1000, np.full(inside.sum(), row.depth_m),
                  12*u[inside], 12*v[inside], np.zeros(inside.sum()), color=color, alpha=.45, arrow_length_ratio=.3)
    ax.plot(results.xc_m.where(results.inner_valid)/1000, results.yc_m.where(results.inner_valid)/1000, results.depth_m, "k.-", label="Inner-fit centreline")
    ax.set(xlim=(-half_width_km,half_width_km), ylim=(-half_width_km,half_width_km),
           zlim=(results.depth_m.max(), results.depth_m.min()), xlabel="East (km)", ylabel="North (km)", zlabel="Depth (m)")
    ax.set_box_aspect((1,1,1.1))
    ax.view_init(elev=24, azim=-55)
    detail = "No accepted Gaussian contours; inner centreline only" if valid.empty else "Contours only at accepted Gaussian fits; arrows: 12 km per 1 m/s"
    ax.set_title(f"{results.method.iloc[0]} · depth structure (vertical scale exaggerated)\n{detail}", fontsize=11)
    ax.legend()
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap), ax=ax, shrink=.55, label="Depth (m)", pad=.1)
    return fig


# def load_study(info_path, eddy):
#     """Read the single source of transect times/centres from the info.md table."""
#     text = Path(info_path).read_text()
#     marker = "Data path: `"
#     if marker not in text:
#         raise ValueError("info.md needs a Data path line")
#     data_path = Path(text.split(marker, 1)[1].split("`", 1)[0])
#     header = "| ID | Eddy | Start UTC | End UTC | Centre lon | Centre lat |"
#     if header not in text:
#         raise ValueError("Transect table header not found in info.md")
#     rows = []
#     for line in text.split(header, 1)[1].splitlines():
#         if not line.strip():
#             if rows:
#                 break
#             continue
#         if not line.startswith("|"):
#             break
#         cells = [v.strip() for v in line.strip().strip("|").split("|")]
#         if set("".join(cells)) <= {"-", ":", " "}:
#             continue
#         if len(cells) != 6:
#             raise ValueError("Each transect row must have six fields")
#         ident, name, start, end, lon, lat = cells
#         if pd.Timestamp(start).tzinfo is None or pd.Timestamp(end).tzinfo is None:
#             raise ValueError("Transect times must include UTC timezone")
#         if pd.Timestamp(end) <= pd.Timestamp(start):
#             raise ValueError("Invalid transect interval")
#         rows.append(dict(id=ident, eddy=name, start=start, end=end,
#                          lon=float(lon), lat=float(lat)))
#     all_rows = pd.DataFrame(rows)
#     if all_rows.empty or all_rows.id.duplicated().any():
#         raise ValueError("Transect IDs must be unique")
#     catalog = all_rows.loc[all_rows.eddy.str.lower().eq(eddy.lower())].set_index("id")
#     if catalog.empty:
#         raise ValueError(f"No transects for {eddy}")
#     return data_path, catalog
def load_study(info_path, eddy):
    """Load transect metadata and the full, unfiltered ADCP dataset."""
    text = Path(info_path).read_text()
    marker = "Data path: `"
    if marker not in text:
        raise ValueError("info.md needs a Data path line")

    data_path = Path(text.split(marker, 1)[1].split("`", 1)[0])

    header = "| ID | Eddy | Start UTC | End UTC | Centre lon | Centre lat |"
    if header not in text:
        raise ValueError("Transect table header not found in info.md")

    rows = []
    for line in text.split(header, 1)[1].splitlines():
        if not line.strip():
            if rows:
                break
            continue
        if not line.startswith("|"):
            break

        cells = [v.strip() for v in line.strip().strip("|").split("|")]
        if set("".join(cells)) <= {"-", ":", " "}:
            continue
        if len(cells) != 6:
            raise ValueError("Each transect row must have six fields")

        ident, name, start, end, lon, lat = cells

        if pd.Timestamp(start).tzinfo is None or pd.Timestamp(end).tzinfo is None:
            raise ValueError("Transect times must include UTC timezone")
        if pd.Timestamp(end) <= pd.Timestamp(start):
            raise ValueError("Invalid transect interval")

        rows.append(dict(
            id=ident, eddy=name, start=start, end=end,
            lon=float(lon), lat=float(lat)
        ))

    all_rows = pd.DataFrame(rows)

    if all_rows.empty or all_rows.id.duplicated().any():
        raise ValueError("Transect IDs must be unique")

    catalog = all_rows.loc[
        all_rows.eddy.str.lower().eq(eddy.lower())
    ].set_index("id")

    if catalog.empty:
        raise ValueError(f"No transects for {eddy}")

    # Load full ADCP dataset without filtering by transect times
    adcp = pd.read_csv(data_path)
    adcp = adcp.rename(columns={
        'U Ship (m/sec)': 'u',
        'V Ship (m/sec)': 'v',
        'Mean Lon': 'lon',
        'Mean Lat': 'lat'
    })
    if eddy == 'Freddy':
        adcp = adcp[(adcp.lat > -32.9)&(adcp.lat < -32.3)&(adcp.lon > 152.8)&(adcp.lon < 153.5)].copy()
    elif eddy == 'Murphy':
        adcp = adcp[(adcp.lat > -29.8)&(adcp.lat < -27.6)&(adcp.lon > 153.5)&(adcp.lon < 156.0)].copy()
    
    return data_path, catalog, adcp
    

def fit_transects(data_path, catalog, backend, fit_settings, min_attendance=None):
    """Independent depth fits per crossing: never carry a centre between times."""
    runs, summary, tables = {}, [], []
    for ident, spec in catalog.iterrows():
        data, coverage, audit = load_observations(data_path, spec.start, spec.end,
                                                 (spec.lon, spec.lat), min_attendance)
        results, samples = fit_depths(data, backend, "SOLO", depth_levels=coverage.index,
                                     **fit_settings)
        results["transect_id"] = ident
        results["xc_lon"] = spec.lon + results.xc_m/(111000*np.cos(np.deg2rad(spec.lat)))
        results["yc_lat"] = spec.lat + results.yc_m/111000
        results["start_utc"], results["end_utc"] = audit["start_utc"], audit["end_utc"]
        runs[ident] = dict(data=data, coverage=coverage, audit=audit, results=results, samples=samples)
        good = results.inner_valid
        full = results.valid
        summary.append(dict(transect=ident, times=data.time.nunique(), bins=len(results),
                            inner_fits=int(good.sum()), gaussian_fits=int(full.sum()),
                            inner_min_m=results.loc[good, "depth_m"].min(),
                            inner_max_m=results.loc[good, "depth_m"].max(),
                            gaussian_max_m=results.loc[full, "depth_m"].max(),
                            surface_status=results.iloc[0].status))
        tables.append(results)
    return runs, pd.DataFrame(summary).set_index("transect"), pd.concat(tables, ignore_index=True)


def plot_transect_comparison(runs, eddy):
    """Show depth structure per crossing relative to its own shallowest inner fit."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True, layout="constrained")
    for ident, run in runs.items():
        r = run["results"]
        ref = r.reference_depth_m.iloc[0] if "reference_depth_m" in r else np.nan
        axes[0].plot(r.offset_east_km, r.depth_m, ".-", ms=3, label=f"{ident} (ref {ref:g} m)")
        axes[1].plot(r.offset_north_km, r.depth_m, ".-", ms=3)
        axes[2].plot(r.Rc_m.where(r.valid)/1000, r.depth_m, ".-", ms=3)
    axes[0].set(ylabel="Depth (m)", xlabel="East centre offset (km)")
    axes[1].set_xlabel("North centre offset (km)")
    axes[2].set_xlabel("Accepted Gaussian Rc (km)")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle(f"{eddy} · independent SOLO profiles for each crossing\nOffsets reference each crossing's shallowest accepted centre; no averaging across times")
    return fig


def outer_extent(row):
    return 2 * row.Rc_m


def plot_sections(results, backend, half_width_km=120, speed_max=1.5):
    """Fixed vertical planes through the shallowest accepted inner centre.

    Centre trajectories are projections onto each section, not necessarily
    intersections with the fixed plane. White strips retain invalid depth bins.
    """
    if not results.valid.any():
        return _no_fit(f"{results.method.iloc[0]} · vertical sections")
    reference = results.loc[results.inner_valid].iloc[0]
    axis = np.linspace(-half_width_km, half_width_km, 180)
    zonal, meridional = [np.full((len(results),len(axis)),np.nan) for _ in range(2)]
    for j, row in results.iterrows():
        if not row.valid:
            continue
        _, zonal[j] = backend.model_uv_at_xy(axis*1000, np.full(len(axis),reference.yc_m), row.xc_m,row.yc_m,_q(row),row.Omega,row.Rc_m)
        meridional[j], _ = backend.model_uv_at_xy(np.full(len(axis),reference.xc_m), axis*1000, row.xc_m,row.yc_m,_q(row),row.Omega,row.Rc_m)
    fig, axes = plt.subplots(1,2,figsize=(13,5),sharey=True,layout="constrained")
    for ax, field, offset, label, title in zip(axes, [zonal,meridional], [results.xc_m,results.yc_m], ["East (km)","North (km)"],
        [f"Zonal slice · northward v · y={reference.yc_m/1000:.1f} km", f"Meridional slice · eastward u · x={reference.xc_m/1000:.1f} km"]):
        image=ax.pcolormesh(axis,results.depth_m,np.ma.masked_invalid(field),cmap="RdBu_r",vmin=-speed_max,vmax=speed_max,shading="nearest")
        ax.plot(offset.where(results.valid)/1000,results.depth_m,"k.--",label="Projected centre")
        ax.set(xlabel=label,title=title,ylim=(results.depth_m.max()+4,results.depth_m.min()-4))
        ax.legend(fontsize=9)
    axes[0].set_ylabel("Depth (m)")
    fig.colorbar(image,ax=axes,label="Reconstructed velocity (m/s)")
    fig.suptitle(f"{results.method.iloc[0]} · fixed planes through centre at {reference.depth_m:g} m\nCentre paths are projected; white rows are rejected fits; no depth interpolation")
    return fig
