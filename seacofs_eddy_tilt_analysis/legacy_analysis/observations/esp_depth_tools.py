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

STD_FIGSIZE = 2.5

def load_backend(path):
    path = Path(path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location("observation_esp_backend", path)
    backend = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(backend)
    return backend, {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def load_observations(path, start, end, origin, min_attendance=None, bounds=None,
                      instrument=None):
    """Keep original depth bins; optional UTC/region selection, no interpolation or QC inference.

    Pass start=end=None to use every observation in the supplied spatial region.
    bounds is (west, east, south, north) in degrees.
    instrument selects the CE23 OS150 or WH75 source; None preserves all sources.
    """
    raw = pd.read_csv(path)
    if instrument is not None:
        instrument_files = {
            "OS150": "IN2023_V06_ADCP_os150nb.nc",
            "WH75": "in2023_v06_ADCP_wh75bb.nc",
        }
        if instrument not in instrument_files:
            raise ValueError("instrument must be OS150, WH75 or None (all instruments)")
        raw = raw.loc[raw["Datafile"].eq(instrument_files[instrument])].copy()
    raw["time"] = pd.to_datetime(raw["Start (UTC)"], utc=True)
    if (start is None) != (end is None):
        raise ValueError("start and end must both be set or both be None")
    selected = (raw if start is None else
                raw.loc[raw.time.between(pd.Timestamp(start), pd.Timestamp(end))]).copy()
    rename = {"Mean Lon": "lon", "Mean Lat": "lat", "Depth": "depth",
              "U (m/sec)": "u", "V (m/sec)": "v"}
    selected = selected.rename(columns=rename)
    for c in rename.values():
        selected[c] = pd.to_numeric(selected[c], errors="coerce")
    if bounds is not None:
        west, east, south, north = bounds
        if not (west < east and south < north):
            raise ValueError("bounds must be (west, east, south, north)")
        selected = selected.loc[selected.lon.between(west, east) &
                                selected.lat.between(south, north)].copy()
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
             "min_attendance": min_attendance, "origin_lon_lat": origin,
             "instrument": instrument, "datafiles": sorted(data["Datafile"].unique().tolist()),
             "bounds_lon_lat": bounds}
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
        return record, "insufficient_outer", use
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
            return record, "nonfinite_outer", use
        if rc > radius_limit:
            return record, "outer_radius_limit", use
        if om * omega <= 0:
            return record, "outer_rotation_reversal", use
        pred = gaussian(rho**2, om, rc)
        record["outer_rmse_m_s"] = float(np.sqrt(np.mean((target - pred)**2)))
        return record, "ok", use
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        record["outer_warning"] = str(exc)
        return record, "outer_optimisation_failed", use


def fit_depths(data, backend, method, *, core_radius=35000., outer_radius=100000.,
               search_half_width=60000., max_centre_step=35000., radius_limit=100000.,
               min_core_points=8, min_outer_points=12, expected_sign=1,
               require_peak_bracket=True, depth_levels=None, solo_initialisation="backend",
               radius_tracking=False, core_radius_factor=None, outer_radius_factor=None,
               initial_Rc=None):
    """Follow accepted inner centres downwards; failed depths do not update the seed.

    Both methods start at the specified geographic origin (x=y=0). SOLO uses a
    fixed transect direction derived from the shallowest selected bin. LATTE
    repeats its core selection twice, using the last accepted ellipse as seed.
    Reconstruction acceptance is independent of inner-centre acceptance.
    With radius_tracking enabled, non-null factors multiply the last accepted
    Gaussian Rc. Until one exists, use initial_Rc or the explicit fixed radii.
    Both thresholds are frozen for each depth; rejected fits never update Rc.
    """
    method = method.upper()
    if method not in ("LATTE", "SOLO"):
        raise ValueError("method must be LATTE or SOLO")
    if solo_initialisation not in ("backend", "local"):
        raise ValueError("solo_initialisation must be backend or local")
    if min(core_radius, outer_radius, search_half_width, max_centre_step, radius_limit) <= 0:
        raise ValueError("All radii and search distances must be positive.")
    if not isinstance(radius_tracking, (bool, np.bool_)):
        raise ValueError("radius_tracking must be a boolean")
    for name, value in (("core_radius_factor", core_radius_factor),
                        ("outer_radius_factor", outer_radius_factor), ("initial_Rc", initial_Rc)):
        if value is not None and (isinstance(value, bool) or not np.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be null or a finite positive number")
    if not np.isfinite([core_radius, outer_radius]).all():
        raise ValueError("Core and outer radii must be finite")
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
    radius_reference = initial_Rc if radius_tracking else None
    radius_reference_depth = np.nan
    rows, samples = [], {}
    levels = sorted(data.depth.unique() if depth_levels is None else depth_levels)
    for depth in levels:
        core_threshold, outer_threshold = core_radius, outer_radius
        if radius_tracking and radius_reference is not None:
            if core_radius_factor is not None:
                core_threshold = core_radius_factor * radius_reference
            if outer_radius_factor is not None:
                outer_threshold = outer_radius_factor * radius_reference
        d = data.loc[data.depth.eq(depth)]
        x, y, u, v = [d[c].to_numpy(float) for c in ("x", "y", "u", "v")]
        rec = dict(method=method, depth_m=depth, n_obs=len(d), n_core=0,
                   core_radius_used_m=core_threshold, outer_radius_used_m=outer_threshold,
                   radius_reference_m=radius_reference if radius_reference is not None else np.nan,
                   radius_reference_depth_m=radius_reference_depth,
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
        sample["outer_fit"] = False
        samples[depth] = sample
        try:
            if len(d) < min_core_points:
                raise ValueError("insufficient_observations")
            if method == "LATTE":
                c, q = centre.copy(), q_previous.copy()
                for _ in range(2):
                    core = _rho(x, y, c, q) <= core_threshold
                    sample["inner"] = core
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
                            seed, search_half_width, core_threshold)
                    else:
                        l0, r0, w, q, om = backend.solo(
                            along[search][order], vt[search][order], vn[search][order],
                            core_thresh=core_threshold)
                rec["warnings"] = "; ".join(str(a.message) for a in caught)
                if not np.isfinite([l0, r0, w, om]).all():
                    raise ValueError("nonfinite_inner")
                if not along[search].min() < l0 < along[search].max():
                    raise ValueError("unsupported_solo_crossing")
                # Cross-track coordinates are measured relative to the same line.
                # The line offset is fixed from the shallowest sampling track.
                offset = surface[["x", "y"]].to_numpy().mean(axis=0) @ normal
                c = l0 * tangent + (offset + r0) * normal
                core = search & (abs(along - l0) <= core_threshold)
                sample["inner"] = core
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
            outer = _rho(x, y, c, q) <= outer_threshold
            sample["outer"] = outer
            out, status, outer_use = _outer(backend, x[outer], y[outer], u[outer], v[outer],
                                            c, q, om, radius_limit, min_outer_points)
            outer_fit = np.zeros(len(d), dtype=bool)
            outer_fit[np.flatnonzero(outer)] = outer_use
            sample["outer_fit"] = outer_fit
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
        if (radius_tracking and rec["valid"] and rec["status"] == "ok"
                and np.isfinite(rec["Rc_m"]) and rec["Rc_m"] > 0):
            radius_reference = rec["Rc_m"]
            radius_reference_depth = depth
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


def plot_buffered_track(data, crossing, seed, transect_id):
    """Show the buffered ADCP track, original SOLO window and starting centre."""
    track = data.sort_values("time").drop_duplicates("time")
    original = track.time.between(pd.Timestamp(crossing.start), pd.Timestamp(crossing.end))
    fig, ax = plt.subplots(figsize=(1.5*STD_FIGSIZE, 1.5*STD_FIGSIZE), layout="constrained")
    ax.plot(track.lon, track.lat, color="0.6", lw=1.5, label="Buffered ADCP window")
    ax.plot(track.loc[original, "lon"], track.loc[original, "lat"],
            color="magenta", lw=2.5, label="Original SOLO crossing")
    ax.scatter(seed.xc_lon, seed.yc_lat, marker="x", s=80, color="red", label="SOLO starting centre")
    ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)",
           title=f"{transect_id}: original crossing and added ADCP track")
    ax.set_aspect(1/np.cos(np.deg2rad(track.lat.median())))
    ax.legend(fontsize=8)
    return fig


def plot_regional_track(data, transects, solo, bounds, *, title,
                        shown_depth=None, transect_ids=None, colours=None,
                        figsize=(9, 6), surface_tracks=False,
                        velocity_label="All regional ADCP", track_width=1.4,
                        centre_size=70, centre_width=2, legend_kwargs=None):
    """Map regional ADCP velocities, SOLO crossing tracks and shallowest valid centres.

    ``shown_depth`` must be an observed bin; defaults to the shallowest bin.
    ``surface_tracks`` restricts crossing tracks to that bin when True.
    Returns the figure without displaying it.
    """
    shown_depth = data.depth.min() if shown_depth is None else shown_depth
    transect_ids = transects.index if transect_ids is None else transect_ids
    if colours is None:
        colours = plt.cm.tab10(np.linspace(0, 1, len(transect_ids)))
    if legend_kwargs is None:
        legend_kwargs = dict(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    surface = data.loc[data.depth.eq(shown_depth)].sort_values("time")
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    shown = surface.iloc[::max(1, len(surface)//120)]
    ax.quiver(shown.lon, shown.lat, shown.u, shown.v, color="0.65", width=.0025,
              scale=17, label=f"{velocity_label}, {shown_depth:g} m (sampled arrows)")
    track_data = surface if surface_tracks else data
    for colour, ident in zip(colours, transect_ids):
        crossing = transects.loc[ident]
        track = track_data.loc[track_data.time.between(pd.Timestamp(crossing.start), pd.Timestamp(crossing.end))]
        track = track.sort_values("time").drop_duplicates("time")
        ax.plot(track.lon, track.lat, color=colour, lw=track_width, label=f"{ident} track")
        fitted = solo.loc[solo.transect_id.eq(ident) & solo.inner_valid].sort_values("depth_m")
        if not fitted.empty:
            first = fitted.iloc[0]
            ax.scatter(first.xc_lon, first.yc_lat, marker="x", s=centre_size, lw=centre_width,
                       color=colour, zorder=5, label=f"{ident} SOLO centre ({first.depth_m:g} m)")
    ax.add_patch(plt.Rectangle((bounds[0], bounds[2]),
                           bounds[1]-bounds[0],
                           bounds[3]-bounds[2],
                           fill=False, edgecolor="firebrick", lw=1.2))
    ax.set(xlabel="Longitude (°E)", ylabel="Latitude (°N)",
           title=title)
    ax.set_aspect(1/np.cos(np.deg2rad(surface.lat.mean())))
    ax.legend(**legend_kwargs)
    return fig


def plot_sampling(data, coverage, half_width_km=120, quiv_scale=12):
    fig, axes = plt.subplots(1, 3, figsize=(3*STD_FIGSIZE, STD_FIGSIZE), layout="constrained")
    d = data.loc[data.depth.eq(data.depth.min())]
    hours = (d.time - d.time.min()).dt.total_seconds() / 3600
    q = axes[0].quiver(d.x/1000, d.y/1000, d.u, d.v, hours, cmap="viridis", scale=quiv_scale)
    axes[0].scatter(0, 0, marker="x", color="red", label="Initial centre")
    axes[0].set(xlabel="East (km)", ylabel="North (km)", title=f"Observed velocities at {d.depth.iloc[0]:g} m", aspect="equal")
    axes[0].legend()#fontsize=8)
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
    fig, axes = plt.subplots(1, 3, figsize=(3*STD_FIGSIZE, STD_FIGSIZE), sharey=True, layout="constrained")
    z = results.depth_m
    axes[0].plot(results.n_obs, z, label="All")
    axes[0].plot(results.n_core, z, label="Last core selection")
    axes[0].set(xlabel="Observation count", ylabel="Depth (m)")
    axes[0].legend()#fontsize=8)
    axes[1].semilogx(results.core_spread_ratio.where(results.core_spread_ratio > 0), z, ".-")
    axes[1].axvline(.05, color="grey", ls="--", label="Narrow-track flag")
    axes[1].set_xlabel("Across/along spatial spread")
    axes[1].legend()#fontsize=8)
    axes[2].plot(results.inner_R2, z, ".-", label="Attempted inner fit")
    axes[2].plot(results.vector_R2, z, ".-", label="Accepted full fit")
    axes[2].set_xlabel("Vector R²")
    axes[2].legend()#fontsize=8)
    axes[0].invert_yaxis()
    fig.suptitle(f"{results.method.iloc[0]} · high R² alone does not establish valid ellipse geometry")
    return fig


def plot_profiles(results):
    if not results.inner_valid.any():
        return _no_fit(f"{results.method.iloc[0]} · depth profiles",
                       "No accepted inner-core fits.\nInspect the sampling and ellipse diagnostics above.")
    fig, axes = plt.subplots(1, 5, figsize=(5*STD_FIGSIZE, STD_FIGSIZE), sharey=True, layout="constrained")
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
    fig, ax = plt.subplots(figsize=(STD_FIGSIZE, STD_FIGSIZE), layout="constrained")
    ax.text(.5, .5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_title(title)
    ax.set_axis_off()
    return fig


def plot_maps(results, samples, backend, depths=(18, 100, 200, 300), half_width_km=120,
              speed_max=1.0, quiv_width=.002, quiv_scale=None):
    """Nearest measured levels to requests; never substitute a deeper valid fit."""
    actual = sorted(set(results.depth_m.iloc[np.argmin(abs(results.depth_m-z))] for z in depths))
    ncols = 4
    nrows = int(np.ceil(len(actual) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows),
                            squeeze=False, layout="constrained")
    axis = np.linspace(-half_width_km, half_width_km, 90)
    xx, yy = np.meshgrid(axis*1000, axis*1000)
    image = None

    xy_diff = []
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
        ax.quiver(d.x/1000, d.y/1000, d.u, d.v, scale=quiv_scale, color="black", width=quiv_width, zorder=7)
        xy_diff.append(np.hypot(np.diff(d.x/1000), np.diff(d.y/1000)).mean())
        
        inner = d["inner"].to_numpy(bool)
        outer_fit = d["outer_fit"].to_numpy(bool) if "outer_fit" in d else d["outer"].to_numpy(bool)
        outer_only = outer_fit & ~inner
        ax.scatter(d.x[outer_only]/1000, d.y[outer_only]/1000, s=18, color="green", alpha=.85, zorder=5)
        ax.scatter(d.x[inner]/1000, d.y[inner]/1000, s=18, color="magenta", alpha=.85, zorder=6)
        if row.inner_valid:
            ax.scatter(row.xc_m/1000, row.yc_m/1000, marker="x", color="red")
        ax.set(xlim=(-half_width_km, half_width_km), ylim=(-half_width_km, half_width_km),
               aspect="equal", xlabel="East (km)", ylabel="North (km)", title=f"{z:g} m")
    print(f"{np.array(xy_diff).mean():.2f} km mean spacing")
    for ax in axs.ravel()[len(actual):]:
        ax.set_axis_off()

    if image is not None:
        fig.colorbar(image, ax=list(axs.ravel()), label=r"Reconstructed speed (m s$^{-1}$)", shrink=.8)
    fig.suptitle(f"{results.method.iloc[0]} · magenta: inner fit · green: additional outer-fit points · "
                 "black: observed flow · grey: reconstructed flow · red: centre/core contour")
    return fig


def plot_solo_satellite_overlay(data, results, samples, map_folder, transect_id,
                                half_width_km=120, speed_max=1.0, target_depth=None, legend=False):
    """Shallowest ESP fit samples over daily OceanCurrent velocity at the middle ADCP time.

    Missing daily DM03 files are cached beside the ADCP CSV. The colours show
    samples used by the inner fit and additional sign-retained outer-fit samples.
    """
    from datetime import datetime, timezone
    from urllib.request import urlopen
    import shutil
    import tempfile
    from netCDF4 import Dataset, num2date

    times = pd.DatetimeIndex(data.time.dropna().unique()).sort_values()
    if len(times) == 0:
        raise ValueError("Transect has no valid ADCP timestamps.")
    middle = pd.Timestamp(times[len(times)//2]).tz_convert("UTC")
    day = middle.floor("D")
    days = [day] if middle == day else [day, day + pd.Timedelta(days=1)]
    folder = Path(map_folder)
    folder.mkdir(parents=True, exist_ok=True)

    def daily_map(date):
        name = f"IMOS_OceanCurrent_HV_{date:%Y%m%d}T000000Z_GSLA_FV02_DM03.nc"
        path = folder / name
        if not path.exists():
            url = ("https://thredds.aodn.org.au/thredds/fileServer/"
                   f"IMOS/OceanCurrent/GSLA/DM/{date:%Y}/{name}")
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=folder, suffix=".nc", delete=False) as target:
                    temporary = Path(target.name)
                    with urlopen(url, timeout=60) as response:
                        shutil.copyfileobj(response, target)
                # Inspect before caching, so a failed response is never mistaken for a map.
                with Dataset(temporary) as ds:
                    if not all(key in ds.variables for key in ("TIME", "LATITUDE", "LONGITUDE", "UCUR", "VCUR")):
                        raise ValueError(f"Incomplete OceanCurrent download: {url}")
                temporary.replace(path)
                print(f"Cached OceanCurrent map: {path.name}")
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        with Dataset(path) as ds:
            time = ds.variables["TIME"]
            if len(time) != 1:
                raise ValueError(f"Expected one daily map in {path}")
            stamp = num2date(time[0], time.units, getattr(time, "calendar", "standard"))
            actual = pd.Timestamp(datetime(stamp.year, stamp.month, stamp.day,
                                           stamp.hour, stamp.minute, stamp.second,
                                           tzinfo=timezone.utc))
            if actual != date:
                raise ValueError(f"OceanCurrent time {actual} does not match {path.name}")
            lon = np.asarray(ds.variables["LONGITUDE"][:], dtype=float)
            lat = np.asarray(ds.variables["LATITUDE"][:], dtype=float)
            def field(name):
                values = ds.variables[name][0, :, :]
                return np.asarray(np.ma.filled(values, np.nan), dtype=float)
            return lon, lat, field("UCUR"), field("VCUR")

    maps = [daily_map(date) for date in days]
    # The user's 3 June subset and later full-domain daily files can have
    # different grids. Bracketing files must share a grid for linear interpolation.
    lon, lat, u, v = maps[0]
    if len(maps) == 2:
        lon_next, lat_next, u_next, v_next = maps[1]
        if not (np.array_equal(lon, lon_next) and np.array_equal(lat, lat_next)):
            raise ValueError("Bracketing OceanCurrent maps use different grids.")
        weight = (middle - day) / pd.Timedelta(days=1)
        u = (1-weight)*u + weight*u_next
        v = (1-weight)*v + weight*v_next

    # surface_depth = results.depth_m.min()
    surface_depth = (
        results.depth_m.min()
        if target_depth is None
        else results.depth_m.iloc[(results.depth_m - target_depth).abs().argmin()]
    )
    surface = samples[surface_depth].sort_values("time")
    mid_lon, mid_lat = float(surface.lon.median()), float(surface.lat.median())
    lat_half = half_width_km / 111
    lon_half = lat_half / np.cos(np.deg2rad(mid_lat))
    ix = (lon >= mid_lon-lon_half) & (lon <= mid_lon+lon_half)
    iy = (lat >= mid_lat-lat_half) & (lat <= mid_lat+lat_half)
    if ix.sum() < 2 or iy.sum() < 2:
        raise ValueError("Satellite grid does not cover the selected transect map.")
    lon, lat = lon[ix], lat[iy]
    u, v = u[np.ix_(iy, ix)], v[np.ix_(iy, ix)]
    xx, yy = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(2*STD_FIGSIZE, 2*STD_FIGSIZE), layout="constrained")
    speed = ax.pcolormesh(xx, yy, np.hypot(u, v), cmap="Blues", vmin=0,
                          vmax=speed_max, shading="auto")
    fig.colorbar(speed, ax=ax, label="OceanCurrent surface geostrophic speed (m/s)")
    ax.quiver(xx, yy, u/np.cos(np.deg2rad(yy)), v, angles="xy",
              scale_units="xy", scale=1/0.25, color="0.35", alpha=.7)
    ax.plot(surface.lon, surface.lat, color="black", lw=1, label="ADCP transect", zorder=4)
    inner = surface["inner"].to_numpy(bool)
    outer_fit = surface["outer_fit"].to_numpy(bool) if "outer_fit" in surface else surface["outer"].to_numpy(bool)
    outer_only = outer_fit & ~inner
    ax.scatter(surface.lon[outer_only], surface.lat[outer_only], s=32, color="green",
               label="Additional outer-fit samples", zorder=5)
    ax.scatter(surface.lon[inner], surface.lat[inner], s=32, color="magenta",
               label="Inner-fit samples", zorder=6)
    row = results.loc[results.depth_m.eq(surface_depth)].iloc[0]
    method = str(row.method)
    if row.inner_valid:
        ax.scatter(row.xc_lon, row.yc_lat, marker="x", s=80, color="red",
                   label=f"{method} centre", zorder=7)
    if row.valid:  # Draw only an accepted Gaussian contour
        ring = _contour(row)  # East/north coordinates in km
        origin_lat = row.yc_lat - row.yc_m / 111_000
        ring_lon = row.xc_lon + (ring[:, 0] * 1000 - row.xc_m) / (
            111_000 * np.cos(np.deg2rad(origin_lat))
        )
        ring_lat = row.yc_lat + (ring[:, 1] * 1000 - row.yc_m) / 111_000
        ax.plot(ring_lon, ring_lat, "r--", lw=1, zorder=7,
                label=f"{method} maximum-speed contour")
    ax.set(xlim=(mid_lon-lon_half, mid_lon+lon_half),
           ylim=(mid_lat-lat_half, mid_lat+lat_half),
           xlabel="Longitude (°E)", ylabel="Latitude (°N)",
           title=(f"{transect_id} · {surface_depth:g} m ADCP over OceanCurrent\n"
                  f"Satellite interpolated to {middle:%Y-%m-%d %H:%M} UTC; "
                  f"{method}: {row.status}"))
    ax.set_aspect(1/np.cos(np.deg2rad(mid_lat)))
    if legend:
        ax.legend(loc="best")#, fontsize=8)
    ax.grid(alpha=.2)
    return fig, ax

def plot_3d(results, backend, half_width_km=120, speed_max=1.5):
    valid = results.loc[results.valid]
    if not results.inner_valid.any():
        return _no_fit(f"{results.method.iloc[0]} · 3D reconstruction")
    fig = plt.figure(figsize=(2*STD_FIGSIZE, 2*STD_FIGSIZE), layout="constrained")
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



def plot_spaced_depth_structures(runs, catalog, transect_ids, *, gap_km=20.,
                                 elev=22, azim=-80, title=False):
    """Chronological centre-relative structures; layout spacing is not geography.

    Draw every accepted Gaussian rho=Rc/sqrt(2) contour at its measured depth.
    Keep rejected inner-centre rows as gaps. Equal east/north display scales and
    extent-based spacing preserve relative sizes without overlapping structures.
    Dates label the original crossing start in UTC, not the buffered window.
    """
    if not np.isfinite(gap_km) or gap_km <= 0:
        raise ValueError("gap_km must be finite and positive")
    ids = list(transect_ids)
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Provide unique transect IDs")
    missing = set(ids) - (set(runs) & set(catalog.index))
    if missing:
        raise ValueError(f"Missing crossings: {sorted(missing)}; fit the corresponding visit first")
    ids.sort(key=lambda ident: pd.to_datetime(catalog.loc[ident, "start"], utc=True))
    structures = []
    for ident in ids:
        result = runs[ident]["results"].sort_values("depth_m")
        accepted = result.loc[result.inner_valid]
        if accepted.empty:
            structures.append((ident, result, None, [], (-10., 10.)))
            continue
        reference = accepted.iloc[0]
        origin = np.array([reference.xc_m, reference.yc_m]) / 1000
        centres = result[["xc_m", "yc_m"]].where(result.inner_valid, axis=0).to_numpy() / 1000 - origin
        rings = [(row.depth_m, _contour(row) - origin)
                 for _, row in result.loc[result.valid & result.inner_valid].iterrows()]
        points = np.concatenate([centres[result.inner_valid.to_numpy()], *[ring for _, ring in rings]])
        structures.append((ident, result, centres, rings,
                           (float(points[:, 0].min()), float(points[:, 0].max()))))
    fig = plt.figure(figsize=(19, 7))
    ax = fig.add_subplot(111, projection="3d")
    ticks, labels, north, depths = [], [], [], []
    right = None
    for i, (ident, result, centres, rings, (left_extent, right_extent)) in enumerate(structures):
        shift = -left_extent if right is None else right + gap_km - left_extent
        right = shift + right_extent
        color = plt.cm.tab10(i % 10)
        ticks.append(shift)
        stamp = pd.to_datetime(catalog.loc[ident, "start"], utc=True)
        labels.append(f"{ident}\n{stamp:%d %b %H:%M}")
        depths.extend(result.depth_m.to_numpy())
        if centres is None:
            ax.text(shift, 0, result.depth_m.min(), f"{ident}: no accepted centre", color=color)
            continue
        ax.plot(centres[:, 0] + shift, centres[:, 1], result.depth_m,
                color=color, marker=".", markersize=3, lw=2, clip_on=False)
        north.extend(centres[np.isfinite(centres[:, 1]), 1])
        for depth, ring in rings:
            ax.plot(ring[:, 0] + shift, ring[:, 1], np.full(len(ring), depth),
                    color=color, lw=.65, alpha=.45, clip_on=False)
            north.extend(ring[:, 1])
    yhalf = max(10., max(abs(np.asarray(north)), default=0.) * 1.08)
    xspan = max(right, 1.)
    zmin, zmax = min(depths), max(depths)
    ax.set(xlim=(-gap_km/2, right + gap_km/2), ylim=(-yhalf, yhalf),
           zlim=(zmax + 1, zmin - 1), ylabel="North offset (km)", zlabel="Depth (m)",
           xlabel="")
    ax.set_xticks(ticks, labels, fontsize=8)
    ax.set_box_aspect((xspan + gap_km, 2*yhalf, max(2*yhalf, .22*xspan)), zoom=1.65)
    ax.set_position([0, .06, 1, .80])
    ax.set_yticks([-round(yhalf*.8), 0, round(yhalf*.8)])
    ax.view_init(elev=elev, azim=azim)
    if title:
        fig.suptitle("Freddy · buffered LATTE depth structures in chronological order\n"
                    "Each centred on its shallowest accepted inner fit; vertical scale exaggerated\n"
                    "Ellipses: accepted Gaussian ρ=Rc/√2 contours at every measured depth; centreline gaps: rejected fits", fontsize=12)
        ig.text(.5, .025, "East offsets with artificial spacing · crossing start (UTC)", ha="center")
    return fig, ax


def load_study(config_path, eddy, overview_bounds=None, load_adcp=True):
    """Read active ADCP crossing definitions from JSON, never from study notes."""
    import json

    config = json.loads(Path(config_path).read_text())
    data_path = Path(config["data_path"])
    rows = config["transects"]
    if not rows:
        raise ValueError("Transect config is empty")
    catalog = pd.DataFrame(rows)
    required = {"id", "eddy", "start", "end", "lon", "lat"}
    if not required.issubset(catalog.columns) or catalog.id.duplicated().any():
        raise ValueError("Transects need unique IDs and eddy/start/end/lon/lat fields")
    for row in catalog.itertuples():
        start, end = pd.Timestamp(row.start), pd.Timestamp(row.end)
        if start.tzinfo is None or end.tzinfo is None or start.utcoffset().total_seconds() != 0 or end.utcoffset().total_seconds() != 0 or end <= start:
            raise ValueError(f"{row.id}: crossing must have increasing UTC start/end")
        if not np.isfinite([row.lon, row.lat]).all() or not (-180 <= row.lon <= 180 and -90 <= row.lat <= 90):
            raise ValueError(f"{row.id}: invalid coordinate origin")
    catalog = catalog.loc[catalog.eddy.str.lower().eq(eddy.lower())].set_index("id")
    if catalog.empty:
        raise ValueError(f"No transects for {eddy}")
    if not load_adcp:
        return data_path, catalog, None

    adcp = pd.read_csv(data_path).rename(columns={
        "U (m/sec)": "u", "V (m/sec)": "v", "Mean Lon": "lon",
        "Mean Lat": "lat", "Start (UTC)": "time"})
    if overview_bounds is not None:
        west, east, south, north = overview_bounds
        if not (west < east and south < north):
            raise ValueError("overview_bounds must be west, east, south, north")
        adcp = adcp.loc[adcp.lon.between(west, east) & adcp.lat.between(south, north)].copy()
    return data_path, catalog, adcp


def _validate_fit_by_transect(catalog, settings):
    """Require exactly one settings block for each selected crossing."""
    if settings is not None and set(settings) != set(catalog.index):
        missing = set(catalog.index) - set(settings)
        unknown = set(settings) - set(catalog.index)
        raise ValueError(f"Invalid fit_by_transect: missing={sorted(missing)}, unknown={sorted(unknown)}")


def fit_transects(data_path, catalog, backend, fit_settings=None, min_attendance=None,
                  *, fit_by_transect=None, instrument=None):
    """Independent depth fits per crossing: never carry a centre between times."""
    _validate_fit_by_transect(catalog, fit_by_transect)
    if fit_by_transect is not None and fit_settings is not None:
        raise ValueError("Use fit_by_transect without shared settings")
    runs, summary, tables = {}, [], []
    for ident, spec in catalog.iterrows():
        data, coverage, audit = load_observations(data_path, spec.start, spec.end,
                                                 (spec.lon, spec.lat), min_attendance,
                                                 instrument=instrument)
        fit = dict(fit_settings if fit_by_transect is None else fit_by_transect[ident])
        audit["fit_settings"] = fit
        results, samples = fit_depths(data, backend, "SOLO", depth_levels=coverage.index,
                                     **fit)
        results["transect_id"] = ident
        results["xc_lon"] = spec.lon + results.xc_m/(111000*np.cos(np.deg2rad(spec.lat)))
        results["yc_lat"] = spec.lat + results.yc_m/111000
        results["start_utc"], results["end_utc"] = audit["start_utc"], audit["end_utc"]
        runs[ident] = dict(data=data, coverage=coverage, audit=audit, results=results, samples=samples, fit_settings=fit)
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


def bin_solo_along_track(data, width_km):
    """Vector-average one SOLO crossing in equal along-track cells at each depth.

    Use the shallowest-bin track direction for every depth. Each occupied cell
    contributes one observed-position/velocity vector; no depths are interpolated.
    Apply to one crossing at a time so separate visits are never averaged.
    """
    width_km = float(width_km)
    if not np.isfinite(width_km) or width_km <= 0:
        raise ValueError("Along-track bin width must be a positive number of km")
    if data.empty:
        raise ValueError("Cannot bin an empty crossing")
    surface = data.loc[data.depth.eq(data.depth.min())]
    xy = surface[["x", "y"]].to_numpy(float)
    if len(xy) < 2 or not np.any(np.ptp(xy, axis=0)):
        raise ValueError("Cannot infer the along-track direction")
    _, _, vh = np.linalg.svd(xy - xy.mean(axis=0), full_matrices=False)
    along = data[["x", "y"]].to_numpy(float) @ vh[0]
    cells = np.floor((along - along.min()) / (width_km * 1000)).astype(int)
    binned = data.assign(along_track_cell=cells).groupby(
        ["depth", "along_track_cell"], as_index=False, sort=True
    ).agg(x=("x", "mean"), y=("y", "mean"),
          u=("u", "mean"), v=("v", "mean"), n_raw=("u", "size"))
    return binned


def plot_solo_binning_comparison(fits, transect_id):
    """Compare SOLO centres, radius and inner vorticity against one raw fit."""
    raw = fits["Original"]
    accepted = raw.loc[raw.inner_valid]
    if accepted.empty:
        return _no_fit(f"{transect_id} · along-track binning", "No accepted original centre")
    reference = accepted.iloc[0]
    fig, axes = plt.subplots(1, 4, figsize=(15, 5), sharey=True, layout="constrained")
    for label, result in fits.items():
        z = result.depth_m
        good = result.inner_valid
        style = dict(marker=".", markersize=3, linewidth=1.3 if label == "Original" else 1)
        axes[0].plot((result.xc_m.where(good) - reference.xc_m) / 1000, z,
                     label=label, **style)
        axes[1].plot((result.yc_m.where(good) - reference.yc_m) / 1000, z, **style)
        axes[2].plot(result.Rc_m.where(result.valid) / 1000, z, **style)
        axes[3].plot(result.inner_w.where(good) * 1e5, z, **style)
    for ax, label in zip(axes, ["East centre offset (km)", "North centre offset (km)",
                                "Gaussian Rc (km)", "Inner vorticity (10⁻⁵ s⁻¹)"]):
        ax.set_xlabel(label)
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Depth (m)")
    axes[0].invert_yaxis()
    axes[0].legend(title="Along-track bin", fontsize=8)
    fig.suptitle(f"{transect_id} · offsets from original centre at {reference.depth_m:g} m; gaps are rejected fits")
    return fig


def fit_buffered_latte_transects(data_path, catalog, solo_results, backend, fit_settings=None,
                                 before_hours=1.0, after_hours=1.0,
                                 buffer_overrides=None, min_attendance=None,
                                 fit_overrides=None, fit_by_transect=None):
    """Fit LATTE independently around each SOLO crossing, including time buffers.

    Each crossing starts at its own shallowest accepted SOLO inner centre. All
    finite ADCP observations in the expanded time window are used, without a
    spatial or straight-line restriction. Buffers are (hours before, hours after).
    """
    _validate_fit_by_transect(catalog, fit_by_transect)
    if fit_by_transect is not None and (fit_settings is not None or fit_overrides):
        raise ValueError("Use fit_by_transect without shared settings or fit_overrides")
    overrides = {} if buffer_overrides is None else dict(buffer_overrides)
    fit_overrides = {} if fit_overrides is None else dict(fit_overrides)
    unknown_fit = set(fit_overrides) - set(catalog.index)
    if unknown_fit:
        raise ValueError(f"Unknown fit override transects: {sorted(unknown_fit)}")
    unknown = set(overrides) - set(catalog.index)
    if unknown:
        raise ValueError(f"Unknown buffer override transects: {sorted(unknown)}")
    required = {"transect_id", "depth_m", "inner_valid", "xc_lon", "yc_lat"}
    if not required.issubset(solo_results.columns):
        raise ValueError(f"Missing SOLO result columns: {required - set(solo_results.columns)}")

    runs, summary, tables = {}, [], []
    for ident, crossing in catalog.iterrows():
        before, after = overrides.get(ident, (before_hours, after_hours))
        if not np.isfinite([before, after]).all() or min(before, after) < 0:
            raise ValueError(f"{ident}: buffers must be finite, nonnegative hours")
        accepted = solo_results.loc[solo_results.transect_id.eq(ident) &
                                    solo_results.inner_valid.astype(str).str.lower().eq("true")]
        accepted = accepted.sort_values("depth_m")
        if accepted.empty:
            raise ValueError(f"{ident}: no accepted SOLO inner centre; rerun SOLO first")
        seed = accepted.iloc[0]
        origin = (float(seed.xc_lon), float(seed.yc_lat))
        start = pd.Timestamp(crossing.start) - pd.Timedelta(hours=float(before))
        end = pd.Timestamp(crossing.end) + pd.Timedelta(hours=float(after))
        data, coverage, audit = load_observations(data_path, start, end, origin,
                                                  min_attendance=min_attendance)
        fit = dict(fit_settings if fit_by_transect is None else fit_by_transect[ident])
        extra = fit_overrides.get(ident, {})
        unknown_keys = set(extra) - set(fit)
        if unknown_keys:
            raise ValueError(f"{ident}: unknown fit override settings: {sorted(unknown_keys)}")
        fit.update(extra)
        results, samples = fit_depths(data, backend, "LATTE", depth_levels=coverage.index,
                                      **fit)
        results["transect_id"] = ident
        results["xc_lon"] = origin[0] + results.xc_m/(111000*np.cos(np.deg2rad(origin[1])))
        results["yc_lat"] = origin[1] + results.yc_m/111000
        results["window_start_utc"] = str(start)
        results["window_end_utc"] = str(end)
        results["buffer_before_hours"] = float(before)
        results["buffer_after_hours"] = float(after)
        audit.update(window_start_utc=str(start), window_end_utc=str(end),
                     buffer_before_hours=float(before), buffer_after_hours=float(after),
                     solo_seed_depth_m=float(seed.depth_m), solo_seed_lon_lat=origin,
                     fit_settings=fit)
        runs[ident] = dict(data=data, coverage=coverage, audit=audit,
                           results=results, samples=samples, solo_seed=seed, fit_settings=fit)
        good, full = results.inner_valid, results.valid
        summary.append(dict(transect=ident, before_h=before, after_h=after,
                            start_utc=start, end_utc=end, times=data.time.nunique(),
                            bins=len(results), inner_fits=int(good.sum()),
                            gaussian_fits=int(full.sum()),
                            surface_status=results.iloc[0].status,
                            solo_seed_depth_m=seed.depth_m))
        tables.append(results)
    return runs, pd.DataFrame(summary).set_index("transect"), pd.concat(tables, ignore_index=True)


def plot_transect_comparison(runs, eddy):
    """Show depth structure per crossing relative to its own shallowest inner fit."""
    method = next(iter(runs.values()))["results"].method.iloc[0]
    fig, axes = plt.subplots(1, 3, figsize=(3*STD_FIGSIZE, STD_FIGSIZE), sharey=True, layout="constrained")
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
    axes[0].legend()#fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle(f"{eddy} · independent {method} profiles for each crossing\nOffsets reference each crossing's shallowest accepted centre; no averaging across times")
    return fig


def plot_buffered_latte_comparison(runs, eddy):
    """Compare independent LATTE depth profiles; never average across visits."""
    if not runs:
        raise ValueError("No buffered LATTE runs")
    fig, axes = plt.subplots(2, 4, figsize=(4*STD_FIGSIZE, 2*STD_FIGSIZE), sharey=True, layout="constrained")
    specs = [
        ("offset_east_km", "East offset (km)", "inner_valid", 1),
        ("offset_north_km", "North offset (km)", "inner_valid", 1),
        ("displacement_km", "Centre displacement (km)", "inner_valid", 1),
        ("Rc_m", "Gaussian Rc (km)", "valid", 1/1000),
        ("alpha", "Axis ratio", "inner_valid", 1),
        ("angle_deg", "Major-axis angle (° from east)", "inner_valid", 1),
        ("inner_w", "Inner vorticity (10⁻⁵ s⁻¹)", "inner_valid", 1e5),
        ("vector_R2", "Gaussian vector R²", "valid", 1),
    ]
    colours = plt.cm.tab10(np.linspace(0, 1, max(len(runs), 2)))
    for colour, (ident, run) in zip(colours, runs.items()):
        r = run["results"]
        style = "--" if ident.startswith("F") and int(ident[1:]) >= 8 else "-"
        for ax, (column, label, mask, factor) in zip(axes.flat, specs):
            value = r[column].where(r[mask]) * factor
            ax.plot(value, r.depth_m, marker=".", markersize=3, lw=1,
                    linestyle=style, color=colour, label=ident)
            ax.set_xlabel(label)
            ax.grid(alpha=.2)
    axes[0, 0].set_ylabel("Depth (m)")
    axes[1, 0].set_ylabel("Depth (m)")
    axes[0, 0].invert_yaxis()
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside right upper", title="Crossing")
    fig.suptitle(f"{eddy} buffered LATTE · independent crossings · gaps are rejected fits\n"
                 "Offsets reference each crossing's shallowest accepted centre; dashed lines mark Freddy visit 2")
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
    fig, axes = plt.subplots(1,2,figsize=(2.5*STD_FIGSIZE, 1.5*STD_FIGSIZE),sharey=True,layout="constrained")
    for ax, field, offset, label, title in zip(axes, [zonal,meridional], [results.xc_m,results.yc_m], ["East (km)","North (km)"],
        [f"Zonal slice\nnorthward v\ny={reference.yc_m/1000:.1f} km", f"Meridional slice\neastward u\nx={reference.xc_m/1000:.1f} km"]):
        image=ax.pcolormesh(axis,results.depth_m,np.ma.masked_invalid(field),cmap="RdBu_r",vmin=-speed_max,vmax=speed_max,shading="nearest")
        ax.plot(offset.where(results.valid)/1000,results.depth_m,"k.--",label="Projected centre")
        ax.set(xlabel=label,title=title,ylim=(results.depth_m.max()+4,results.depth_m.min()-4))
        # ax.legend()#(fontsize=9)
    axes[0].set_ylabel("Depth (m)")
    fig.colorbar(image,ax=axes,label="Reconstructed velocity (m/s)")
    fig.suptitle(f"{results.method.iloc[0]} · fixed planes through centre at {reference.depth_m:g} m")
    return fig
