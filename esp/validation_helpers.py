"""Small synthetic experiments for notebook 01; never edits the ESP backend.

Coordinates/radii are metres, velocity m/s, rotation/vorticity s^-1.
Angles describe the major axis anticlockwise from east, modulo 180 degrees.
"""
from __future__ import annotations

import hashlib
import importlib.util
import platform
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def load_backend(root):
    """Import the exact requested file under a unique name, without sys.path edits."""
    path = (Path(root).expanduser() / "functions.py").resolve()
    if not path.is_file():
        raise FileNotFoundError(f"ESP backend not found: {path}. Set ESP_ROOT on Katana.")
    spec = importlib.util.spec_from_file_location("esp_validation_backend", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in ("solo", "doppio", "latte", "out_core_param_fit"):
        if not callable(getattr(module, name, None)):
            raise AttributeError(f"{path} must expose {name}().")
    import scipy
    import matplotlib
    metadata = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    python=platform.python_version(), numpy=np.__version__,
                    pandas=pd.__version__, scipy=scipy.__version__,
                    matplotlib=matplotlib.__version__)
    try:
        metadata["git_commit"] = subprocess.check_output(
            ["git", "-C", str(path.parent), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
        metadata["git_status"] = subprocess.check_output(
            ["git", "-C", str(path.parent), "status", "--porcelain"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        metadata["git_commit"] = None
    return module, metadata


def shape_matrix(alpha, angle_deg):
    a = np.deg2rad(angle_deg)
    rotation = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return rotation @ np.diag([1 / alpha, alpha]) @ rotation.T


def velocity(xy, centre, omega, rc, q):
    r = np.asarray(xy) - centre
    qr = r @ q
    factor = omega * np.exp(-np.sum(r * qr, axis=1) / rc**2)
    return np.column_stack((-qr[:, 1], qr[:, 0])) * factor[:, None]


def geometry(q):
    values, vectors = np.linalg.eigh(q)
    major = vectors[:, 0]  # smallest Q eigenvalue -> longest ellipse axis
    return np.sqrt(values[1] / values[0]), np.rad2deg(np.arctan2(major[1], major[0])) % 180


def axial_error(estimate, truth):
    return (estimate - truth + 90) % 180 - 90


def make_cases(rc=85000., omega=-6.9e-5, alpha=1.527525,
               angle_deg=35., spacing=2500., outer_extent=2., seed=1729):
    """Fixed representative geometries, not a reproduction of manuscript tables.

    SOLO: eastward transect at y=20 km. DOPPIO: orthogonal transects
    intersecting at (10,10) km. LATTE: 240 points in a disk of radius .65 Rc
    and 160 in the surrounding annulus to outer_extent Rc, centred at (10,10) km.
    All observed points are available to the outer fit; its domain is fixed.
    """
    if rc <= 0 or spacing <= 0 or outer_extent <= .65 or alpha < 1 or omega == 0:
        raise ValueError("Invalid synthetic controls.")
    centre = np.zeros(2)
    n = int(np.ceil(outer_extent * rc / spacing))
    axis = np.arange(-n, n + 1) * spacing
    cases = {}
    solo_xy = np.column_stack((axis, np.full(axis.size, 20000.)))
    cases["SOLO"] = dict(xy=solo_xy, q=np.eye(2), selection_centre=np.array([0., 20000.]))
    cross = np.array([10000., 10000.])
    horizontal = np.column_stack((axis + cross[0], np.full(axis.size, cross[1])))
    vertical = np.column_stack((np.full(axis.size, cross[0]), axis + cross[1]))
    # One observation at the shared intersection, reused by both transects.
    xy = np.unique(np.vstack((horizontal, vertical)), axis=0)
    cases["DOPPIO"] = dict(xy=xy, q=shape_matrix(alpha, angle_deg), selection_centre=cross)
    rng = np.random.default_rng(seed)
    radii = rc * np.r_[.65 * np.sqrt(rng.random(240)),
                      np.sqrt(.65**2 + (outer_extent**2 - .65**2) * rng.random(160))]
    angles = rng.uniform(0, 2 * np.pi, len(radii))
    xy = cross + np.column_stack((radii * np.cos(angles), radii * np.sin(angles)))
    cases["LATTE"] = dict(xy=xy, q=shape_matrix(alpha, angle_deg), selection_centre=cross)
    for name, case in cases.items():
        case.update(name=name, centre=centre.copy(), omega=omega, rc=rc)
        case["uv"] = velocity(case["xy"], centre, omega, rc, case["q"])
        case["alpha"], case["angle"] = geometry(case["q"])
        case["vorticity"] = omega * np.trace(case["q"])
        # Peak absolute NORMALISED tangential speed, common to both eddy shapes.
        case["speed_ref"] = abs(omega) * rc / np.sqrt(2 * np.e)
    return cases


def noise_factor(xy, length):
    """Gaussian spatial covariance; independent u/v, unit marginal variance.

    C_ij = exp(-distance_ij^2 / (2 length^2)); length is in metres.
    A tiny diagonal jitter only stabilises the Cholesky decomposition.
    """
    if length <= 0:
        raise ValueError("Correlation length must be positive.")
    delta = xy[:, None, :] - xy[None, :, :]
    cov = np.exp(-np.sum(delta**2, axis=2) / (2 * length**2))
    return np.linalg.cholesky((cov + 1e-10 * np.eye(len(xy))) / (1 + 1e-10))


def inner_indices(case, radius):
    xy, c = case["xy"], case["selection_centre"]
    if case["name"] == "DOPPIO":
        i = np.flatnonzero(np.isclose(xy[:, 1], c[1]) & (abs(xy[:, 0] - c[0]) <= radius))
        j = np.flatnonzero(np.isclose(xy[:, 0], c[0]) & (abs(xy[:, 1] - c[1]) <= radius))
        return i[np.argsort(xy[i, 0])], j[np.argsort(xy[j, 1])]
    if case["name"] == "SOLO":
        return np.flatnonzero(abs(xy[:, 0] - c[0]) <= radius)
    return np.flatnonzero(np.linalg.norm(xy - c, axis=1) <= radius)


def fit_case(backend, case, uv, radius, rc_limit=100000., *,
             diagnostic_outer=False, solo_local=False, centre_guess=10000., search_half_width=60000.):
    """Original estimators by default; explicit opt-in local variants for notebook 01.

    The original outer fit does not return optimiser status and may return its
    initial values on failure. Pass an explicit data-derived seed and flag an
    unchanged pair conservatively as outer_seed_return (not confirmed success).
    This is a diagnostic heuristic, not a replacement optimiser status API.
    """
    xy = case["xy"]
    result = dict(status="unattempted", valid=False, warnings="", n_outer=len(xy),
                  n_core=0, n_outer_retained=0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            index = inner_indices(case, radius)
            if case["name"] == "SOLO":
                # Pass full transect: SOLO itself finds closest approach then masks.
                if solo_local:
                    from local_estimators import solo_local_initialisation
                    x0, y0, _, q, om = solo_local_initialisation(
                        xy[:, 0], uv[:, 0], uv[:, 1], core_thresh=radius,
                        centre_guess=centre_guess, search_half_width=search_half_width)
                else:
                    x0, y0, _, q, om = backend.solo(xy[:, 0], uv[:, 0], uv[:, 1], core_thresh=radius)
                result["n_core"] = int(np.sum(abs(xy[:, 0] - x0) <= radius))
                yc = y0 + xy[0, 1]
                xc = x0
            elif case["name"] == "DOPPIO":
                i, j = index
                result["n_core"] = len(np.union1d(i, j))
                if min(len(i), len(j)) < 4:
                    result["status"] = "insufficient_core"
                    return result
                xc, yc, _, q, om = backend.doppio(
                    xy[i, 0], xy[i, 1], uv[i, 0], uv[i, 1],
                    xy[j, 0], xy[j, 1], uv[j, 0], uv[j, 1])
            else:
                result["n_core"] = len(index)
                if len(index) < 3 or np.linalg.matrix_rank(np.c_[np.ones(len(index)), xy[index]]) < 3:
                    result["status"] = "insufficient_core"
                    return result
                xc, yc, _, q, om, _ = backend.latte(xy[index, 0], xy[index, 1], uv[index, 0], uv[index, 1])
            q = np.asarray(q)
            if not np.all(np.isfinite([xc, yc, om])) or not np.all(np.isfinite(q)):
                result["status"] = "nonfinite_inner"
                return result
            if (q.shape != (2, 2) or not np.allclose(q, q.T) or
                    np.min(np.linalg.eigvalsh(q)) <= 0 or not np.isclose(np.linalg.det(q), 1., rtol=1e-5)):
                result["status"] = "invalid_ellipse"
                return result
            if om == 0:
                result["status"] = "zero_rotation"
                return result
            # Match the backend's tangential projection, geometry normalisation,
            # rotation-sign filter and peak-based Rc initialisation.
            r = xy - [xc, yc]
            qr = r @ q
            rho = np.sqrt(np.maximum(np.sum(r * qr, axis=1), 0))
            norm2 = np.sum(qr**2, axis=1)
            safe = norm2 > 0
            vstar = np.full(len(xy), np.nan)
            vstar[safe] = ((-qr[safe, 1] * uv[safe, 0] + qr[safe, 0] * uv[safe, 1])
                          * rho[safe] / norm2[safe])
            keep = safe & np.isfinite(vstar) & ((vstar <= 0) if om < 0 else (vstar >= 0))
            result["n_outer_retained"] = int(keep.sum())
            if keep.sum() < 3 or np.unique(rho[keep]).size < 3:
                result["status"] = "insufficient_outer"
                return result
            rc0 = max(float(rho[keep][np.argmax(abs(vstar[keep]))] * np.sqrt(2)), 1e-6)
            result.update(xc_m=xc, yc_m=yc, inner_Omega=om, Rc_seed_m=rc0)
            if diagnostic_outer:
                from local_estimators import outer_fit_diagnostic
                outer = outer_fit_diagnostic(xy[:,0],xy[:,1],uv[:,0],uv[:,1],xc,yc,q,
                    Omega0=om,Rc0=rc0,Rc_max=rc_limit)
                result.update(outer)
                if outer['status'] != 'ok':
                    return result
                rc, psi, om_final = outer['Rc_m'], outer['psi0'], outer['Omega']
            else:
                rc, psi, om_final = backend.out_core_param_fit(
                    xy[:, 0], xy[:, 1], uv[:, 0], uv[:, 1], xc, yc, q,
                    Omega0=om, Rc0=rc0, Rc_max=rc_limit, plot=False)
                result.update(Rc_m=rc, Omega=om_final, psi0=psi)
                if not np.all(np.isfinite([rc, psi, om_final])) or rc <= 0:
                    result['status'] = 'nonfinite_outer'
                    return result
                if np.isclose(rc, rc0, rtol=1e-12, atol=0) and np.isclose(om_final, om, rtol=1e-12, atol=0):
                    result['status'] = 'outer_seed_return'
                    return result
                if rc > rc_limit:
                    result['status'] = 'outer_radius_limit'
                    return result
            if np.sign(om_final) != np.sign(case["omega"]):
                result["status"] = "wrong_rotation_sign"
                return result
            alpha, angle = geometry(q)
            w = om_final * np.trace(q)
            distance = np.linalg.norm(np.array([xc, yc]) - case["centre"])
            result.update(status="ok", valid=True, centre_error_km=distance / 1000,
                          centre_error_pct=100 * distance / case["rc"],
                          Rc_error_pct=100 * (rc / case["rc"] - 1),
                          vorticity_error_pct=100 * (w / case["vorticity"] - 1),
                          Omega_error_pct=100 * (om_final / case["omega"] - 1),
                          alpha=alpha, vorticity=w, angle_deg=angle,
                          alpha_error_pct=100 * (alpha / case["alpha"] - 1) if case["name"] != "SOLO" else np.nan,
                          angle_error_deg=axial_error(angle, case["angle"]) if case["name"] != "SOLO" else np.nan)
        except Exception as exc:
            from local_estimators import FitFailure
            result.update(status=str(exc) if isinstance(exc, FitFailure) else "exception",
                          exception=f"{type(exc).__name__}: {exc}")
        finally:
            result["warnings"] = " | ".join(sorted({str(w.message) for w in caught}))
    return result


def run_experiment(backend, cases, *, experiment, levels, repeats, seed,
                   baseline_radius=30000., fixed_noise=.05,
                   correlation_fraction=.15, rc_limit=100000.,
                   diagnostic_outer=False, solo_local=False, centre_guess=10000., search_half_width=60000.):
    """Paired realisations across levels; separate seeds for noise/window studies.

    levels: noise/speed_ref for 'noise', radius/Rc for 'window'. Zero-noise
    control is evaluated once per method. Correlated fields span ALL observations.
    """
    if experiment not in ("noise", "window") or repeats < 1 or min(levels) < 0:
        raise ValueError("Invalid experiment settings.")
    rows = []
    exp_id = {"noise": 1, "window": 2}[experiment]
    for method_id, (name, case) in enumerate(cases.items()):
        factor = noise_factor(case["xy"], correlation_fraction * case["rc"])
        for noise_id, kind in enumerate(("independent", "correlated")):
            for repeat in range(repeats):
                rng = np.random.default_rng(np.random.SeedSequence([seed, exp_id, method_id, noise_id, repeat]))
                error = rng.standard_normal(case["uv"].shape)
                if kind == "correlated":
                    error = factor @ error
                for level in levels:
                    fraction = float(level) if experiment == "noise" else fixed_noise
                    if fraction == 0 and (repeat > 0 or kind != "independent"):
                        continue
                    radius = baseline_radius if experiment == "noise" else float(level) * case["rc"]
                    sigma = fraction * case["speed_ref"]
                    row = fit_case(backend, case, case["uv"] + sigma * error, radius, rc_limit,
                        diagnostic_outer=diagnostic_outer, solo_local=solo_local,
                        centre_guess=centre_guess, search_half_width=search_half_width)
                    row.update(method="SOLO local" if solo_local and name == "SOLO" else name, experiment=experiment, noise_kind=kind if fraction else "none",
                               repeat=repeat, level=float(level), noise_fraction=fraction,
                               noise_sigma_ms=sigma, core_radius_m=radius,
                               core_radius_ratio=radius / case["rc"],
                               correlation_length_m=correlation_fraction * case["rc"] if kind == "correlated" else 0.)
                    rows.append(row)
        print(f"{experiment}: finished {name}", flush=True)
    return pd.DataFrame(rows)


METRICS = dict(centre_error_pct="Centre error (% Rc)", Rc_error_pct="Rc error (%)",
               vorticity_error_pct="Core vorticity error (%)", alpha_error_pct="Aspect-ratio error (%)",
               angle_error_deg="Major-axis error (degrees)")


def summarise(results):
    rows = []
    for keys, group in results.groupby(["experiment", "method", "noise_kind", "level"], sort=False):
        row = dict(zip(["experiment", "method", "noise_kind", "level"], keys))
        row.update(n_attempted=len(group), n_valid=int(group.valid.sum()), valid_fraction=group.valid.mean())
        for metric in (*METRICS, "centre_error_km", "Omega_error_pct"):
            values = group.loc[group.valid, metric].dropna() if metric in group else pd.Series(dtype=float)
            for suffix, quantile in (("p05", .05), ("median", .5), ("p95", .95)):
                row[f"{metric}_{suffix}"] = values.quantile(quantile) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot_results(results):
    import matplotlib.pyplot as plt
    summary = summarise(results)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    colours = {"SOLO":"tab:blue", "SOLO local":"tab:purple", "DOPPIO":"tab:orange", "LATTE":"tab:green"}
    noise = results.experiment.iloc[0] == "noise"
    for (method, kind), group in summary.groupby(["method", "noise_kind"]):
        if kind == "none":
            continue
        if noise:
            baseline = summary[(summary.method == method) & (summary.noise_kind == "none")]
            group = pd.concat([baseline, group])
        group = group.sort_values("level")
        x = group.level.to_numpy() * (100 if noise else 1)
        style = "-" if kind == "independent" else "--"
        label = f"{method}, {kind}"
        for ax, metric in zip(axes.flat, METRICS):
            if method.startswith("SOLO") and metric in ("alpha_error_pct", "angle_error_deg"):
                continue
            ax.plot(x, group[f"{metric}_median"], style, color=colours[method], label=label)
            ax.fill_between(x, group[f"{metric}_p05"], group[f"{metric}_p95"], color=colours[method], alpha=.10)
        axes.flat[5].plot(x, group.valid_fraction * 100, style, color=colours[method], label=label)
    for ax, title in zip(axes.flat, [*METRICS.values(), "Valid fits (%)"]):
        ax.set_title(title)
        ax.set_xlabel("Noise SD / reference speed (%)" if noise else "Inner window / Rc")
        ax.grid(alpha=.2)
        if title != "Valid fits (%)":
            ax.axhline(0, color="grey", lw=.6)
    axes.flat[5].set_ylim(-3, 103)
    axes.flat[5].legend(fontsize=7)
    fig.suptitle("Median and central 90% simulation spread among valid fits; failures retained in valid-fit fraction", fontsize=11)
    return fig


def plot_sampling(cases, radius):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, (name, case) in zip(axes, cases.items()):
        xy = case["xy"]
        index = inner_indices(case, radius)
        index = np.union1d(*index) if isinstance(index, tuple) else index
        ax.scatter(*(xy / 1000).T, s=5, color=".7", label="Outer input (all points)")
        ax.scatter(*(xy[index] / 1000).T, s=8, color="purple", label="Nominal inner input")
        values, vectors = np.linalg.eigh(case["q"])
        t = np.linspace(0, 2 * np.pi, 200)
        ellipse = (np.column_stack((np.cos(t), np.sin(t))) / np.sqrt(values)) @ vectors.T
        ellipse *= case["rc"] / np.sqrt(2)
        ax.plot(*(ellipse / 1000).T, "k--", lw=1, label="Peak normalised speed")
        ax.scatter(0, 0, marker="x", color="black")
        ax.set(title=name, xlabel="East (km)", ylabel="North (km)", aspect="equal")
    axes[0].legend(fontsize=7)
    return fig
