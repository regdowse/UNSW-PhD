"""Controlled RC3 experiments using the manuscript's unchanged ESP functions.

Public coordinates, Rc and streamfunction use SI. Only the inversion adapter
converts to the original notebooks' km and m/s, preserving numerical scaling,
seeds, bounds and stopping rules. Functions are selected by AST from the source
files; notebook imports, data loading and plotting cells are never executed.
"""
from __future__ import annotations

import ast
from contextlib import contextmanager, redirect_stdout
from dataclasses import asdict, dataclass
import hashlib
import io
import json
from pathlib import Path
import platform
from types import SimpleNamespace
from unittest.mock import patch
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy import optimize


@dataclass(frozen=True)
class Config:
    source_root: str = "/srv/scratch/z5297792/UNSW-MRes/MRes"
    seed: int = 20261027
    sample_counts: tuple = (11, 21, 51, 101)
    noise_fractions: tuple = (0., .025, .05, .10)
    replicates: int = 100
    background_fractions: tuple = (0., .05, .10, .20, .30, .50)
    background_n: int = 101
    methods: tuple = ("SOLO", "DOPPIO", "LATTE")
    all_circular: bool = False  # Additional comparison; original setup is default.


NOTEBOOKS = {
    "SOLO": "solo_idealised_tests_updated.ipynb",
    "DOPPIO": "doppio_idealised_tests_updated.ipynb",
    "LATTE": "latte_idealised_tests.ipynb",
}
METRICS = ("centre_error_km", "rc_error_pct", "alpha_error_pct")


def _definitions(text, names, filename, namespace):
    """Compile the exact original top-level definitions, without rewriting ASTs."""
    tree = ast.parse(text, filename=filename)
    tree.body = [node for node in tree.body
                 if isinstance(node, ast.FunctionDef) and node.name in names]
    exec(compile(tree, filename, "exec"), namespace)
    return {node.name for node in tree.body}


def load_backend(config):
    root = Path(config.source_root).expanduser().resolve()
    utils_path = root / "modules/utils.py"
    utils = dict(np=np, pd=pd)
    required = {"unit_det", "axis_ratio_from_Q", "tangential_velocity",
                "psi_params", "gaussian_vel_reconstruction"}
    found = _definitions(utils_path.read_text(), required, str(utils_path), utils)
    if found != required:
        raise RuntimeError(f"Missing original utilities: {required - found}")
    sources = {str(utils_path): hashlib.sha256(utils_path.read_bytes()).hexdigest()}
    methods = {}
    for method in config.methods:
        path = root / "MOCA_PAPER/final_ESP_paper_notebooks" / NOTEBOOKS[method]
        namespace = dict(np=np, pd=pd)
        required = {method.lower(), "out_core_param_fit"}
        found = set()
        for cell in json.loads(path.read_text())["cells"]:
            if cell["cell_type"] == "code":
                found |= _definitions("".join(cell["source"]), required,
                                      str(path), namespace)
        if found != required:
            raise RuntimeError(f"Missing original functions in {path}: {required - found}")
        methods[method] = SimpleNamespace(**{k: namespace[k] for k in required})
        sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return SimpleNamespace(**{k: utils[k] for k in utils if k in {
        "unit_det", "axis_ratio_from_Q", "tangential_velocity",
        "psi_params", "gaussian_vel_reconstruction"}}, methods=methods, sources=sources)


def truth_parameters(backend, method="DOPPIO", all_circular=False):
    """Manuscript circular control for SOLO; elliptical control otherwise.

    Original psi0=250 (km m/s), Rc=85 km; SI psi0=250,000 m²/s.
    u=+dpsi/dy, v=-dpsi/dx, Omega=-2 psi0/Rc²; x east, y north.
    Alpha is major/minor semiaxis length. Rc is the Gaussian scale, with
    normalised tangential-speed peak at elliptical radius Rc/sqrt(2).
    """
    off_diagonal = 0. if method == "SOLO" or all_circular else .4
    q, _ = backend.unit_det(np.array([[1., off_diagonal], [off_diagonal, 1.]]))
    rc, psi = 85_000., 250_000.
    eigenvalues, eigenvectors = np.linalg.eigh(q)
    truth = dict(xc_m=0., yc_m=0., rc_m=rc, psi0_m2_s=psi,
                 omega_s=-2 * psi / rc**2, q=q,
                 alpha=float(backend.axis_ratio_from_Q(q)),
                 major_axis_deg=(None if off_diagonal == 0 else float(np.degrees(np.arctan2(
                     eigenvectors[1, 0], eigenvectors[0, 0])) % 180)))
    truth["vorticity_s"] = truth["omega_s"] * np.trace(q)
    # The full-field maximum is on the eigenvector of Q's largest eigenvalue.
    # Evaluate the EXISTING model there, independently of any degraded samples.
    peak_xy = eigenvectors[:, -1] * rc / np.sqrt(2 * eigenvalues[-1])
    truth["u_ref_m_s"] = float(np.linalg.norm(velocity(backend, peak_xy[None], truth)))
    return truth


def velocity(backend, xy_m, truth):
    q = truth["q"]
    u, v, _, _ = backend.gaussian_vel_reconstruction(
        truth["xc_m"], truth["yc_m"], q[0, 0], q[0, 1], q[1, 1],
        truth["rc_m"], truth["psi0_m2_s"], X=xy_m[:, 0], Y=xy_m[:, 1])
    return np.column_stack((u, v))


def make_geometry(method, n, seed):
    """Fixed coverage; no geometry is regenerated between noise replicates.

    SOLO: manuscript P1/L2 eastward transect, x=-70..70 km, y=20 km.
    N is total count; original default core_thresh=30000 in km is retained.
    DOPPIO: manuscript P1/L2 cross at (10,10) km, half-length 30 km.
    N counts points PER transect; one intersection is shared (2N-1 unique).
    LATTE: scattered extension of the idealised [-50,50] km footprint.
    Four corner anchors make the convex hull identical for every N; the
    remaining nested points are uniform iid in the square. N is total count.
    """
    if int(n) != n or n < 5:
        raise ValueError("N must be an integer >= 5.")
    if method == "SOLO":
        # Manuscript P1/L2: x=-70..70 km, y=20 km; N is total.
        return dict(xy=np.column_stack((np.linspace(-70_000., 70_000., n),
                                        np.full(n, 20_000.))))
    if method == "DOPPIO":
        if n % 2 != 1:
            raise ValueError("DOPPIO requires odd N to retain the exact intersection.")
        axis = np.linspace(-30_000., 30_000., n)
        horizontal = np.column_stack((axis + 10_000., np.full(n, 10_000.)))
        vertical = np.column_stack((np.full(n, 10_000.), axis + 10_000.))
        xy, inverse = np.unique(np.vstack((horizontal, vertical)), axis=0,
                                return_inverse=True)
        return dict(xy=xy, horizontal=inverse[:n], vertical=inverse[n:])
    if method != "LATTE":
        raise ValueError(f"Unsupported method: {method}")
    rng = np.random.default_rng(np.random.SeedSequence([seed, 17]))
    corners = np.array([[-1., -1.], [-1., 1.], [1., -1.], [1., 1.]]) * 50_000.
    xy = np.vstack((corners, rng.uniform(-50_000., 50_000., size=(n - 4, 2))))
    return dict(xy=xy)


@contextmanager
def _observe_optimizer(name):
    """Observe actual SciPy returns/exceptions; forward all arguments unchanged.

    The original LATTE discards convergence status and the outer fit silently
    returns seeds after errors/cap rejection. A scoped observer exposes these
    events without replacing an objective or solver. Runs are sequential:
    patching SciPy's dispatcher is process-global and is NOT thread safe.
    """
    original = getattr(optimize, name)
    calls = []

    def observed(*args, **kwargs):
        entry = {}
        calls.append(entry)
        try:
            result = original(*args, **kwargs)
        except Exception as exc:
            entry.update(success=False, exception=f"{type(exc).__name__}: {exc}")
            raise
        if name == "least_squares":
            entry.update(success=bool(result.success), status=int(result.status),
                         message=str(result.message), nfev=int(result.nfev))
        else:
            entry.update(success=True, parameters=np.asarray(result[0]).copy())
        return result

    with patch.object(optimize, name, observed):
        yield calls


def fit_observations(backend, method, geometry, uv, truth):
    """Original km-scale inversion, with explicit convergence/physical statuses."""
    xy = geometry["xy"] / 1000.
    result = dict(success=False, inner_success=False, status="unattempted",
                  inner_optimizer_success=False, outer_optimizer_success=False,
                  inner_nfev=np.nan, n_unique=len(xy), warning_messages="",
                  shape_estimated=method != "SOLO",
                  optimizer_message="", n_outer=0,
                  **{key: np.nan for key in (*METRICS, "xc_m", "yc_m", "rc_m",
                     "alpha", "q11", "q12", "q22", "q_det", "q_min_eigenvalue",
                     "omega_inner_s", "omega_s", "psi0_m2_s", "vorticity_s",
                     "dx_km", "dy_km", "attempted_rc_m")})
    functions = backend.methods[method]
    with warnings.catch_warnings(record=True) as caught, redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        try:
            optimizer = "curve_fit" if method == "SOLO" else "least_squares"
            with _observe_optimizer(optimizer) as calls:
                if method == "SOLO":
                    fitted = functions.solo(xy[:, 0] - xy[0, 0], uv[:, 0], uv[:, 1], plot=False)
                elif method == "DOPPIO":
                    i, j = geometry["horizontal"], geometry["vertical"]
                    if min(len(i), len(j)) < 4:
                        result["status"] = "insufficient_transect_samples"
                        return result
                    fitted = functions.doppio(
                        xy[i, 0], xy[i, 1], uv[i, 0], uv[i, 1],
                        xy[j, 0], xy[j, 1], uv[j, 0], uv[j, 1], plot=False)
                else:
                    if np.linalg.matrix_rank(np.c_[np.ones(len(xy)), xy]) < 3:
                        result["status"] = "unidentifiable_geometry"
                        return result
                    fitted = functions.latte(xy[:, 0], xy[:, 1], uv[:, 0], uv[:, 1])
            if calls:
                result.update(inner_optimizer_success=all(call["success"] for call in calls),
                              inner_nfev=calls[-1].get("nfev", np.nan),
                              optimizer_message=calls[-1].get("message", calls[-1].get("exception", "")))
            xc, yc, _, q, omega = fitted[:5]
            if method == "SOLO":
                xc, yc = xc + xy[0, 0], yc + xy[0, 1]
            result.update(xc_m=1000 * xc, yc_m=1000 * yc, omega_inner_s=omega / 1000.)
            if not calls or not all(call["success"] for call in calls):
                result["status"] = "inner_optimizer_failed"
                return result
            q = np.asarray(q)
            if q.shape != (2, 2) or not np.all(np.isfinite(q)) or not np.all(
                    np.isfinite([xc, yc, omega])):
                result["status"] = "nonfinite_inner"
                return result
            result.update(q11=q[0, 0], q12=q[0, 1], q22=q[1, 1],
                          q_det=np.linalg.det(q), q_min_eigenvalue=np.linalg.eigvalsh(q)[0])
            if (not np.allclose(q, q.T, atol=1e-10) or result["q_min_eigenvalue"] <= 0
                    or not np.isclose(result["q_det"], 1., rtol=1e-8)):
                result["status"] = "invalid_shape_matrix"
                return result
            result["inner_success"] = True
            result["alpha"] = backend.axis_ratio_from_Q(q)
            result["dx_km"] = xc - truth["xc_m"] / 1000.
            result["dy_km"] = yc - truth["yc_m"] / 1000.
            result["centre_error_km"] = np.hypot(result["dx_km"], result["dy_km"])
            if method != "SOLO":
                result["alpha_error_pct"] = 100 * abs(result["alpha"] / truth["alpha"] - 1)
            frame = backend.psi_params(xc, yc, q, xy[:, 0], xy[:, 1], uv[:, 0], uv[:, 1])
            good = (np.isfinite(frame).all(axis=1) & (frame.rho2 >= 0) & (frame.Qr != 0))
            result["n_outer"] = int(good.sum())
            if good.sum() < 3 or frame.loc[good, "rho2"].nunique() < 3:
                result["status"] = "insufficient_outer_samples"
                return result
            # Preserve the exact original call choices (limits are in km).
            kwargs = {"Rc_max": 1e6} if method in ("SOLO", "DOPPIO") else {}
            with _observe_optimizer("curve_fit") as calls:
                rc, psi, omega_outer = functions.out_core_param_fit(
                    frame.rho2, frame.Qr, frame.vt, Omega0=omega, plot=False, **kwargs)
            result.update(rc_m=1000 * rc, psi0_m2_s=1000 * psi,
                          omega_s=omega_outer / 1000.,
                          vorticity_s=omega_outer / 1000. * np.trace(q))
            if not calls or not calls[-1]["success"]:
                result["status"] = "outer_optimizer_failed"
                result["optimizer_message"] = calls[-1].get("exception", "") if calls else "not called"
                return result
            result["outer_optimizer_success"] = True
            attempted = calls[-1]["parameters"]
            result["attempted_rc_m"] = 1000 * attempted[1]
            if not np.isfinite(attempted[1]) or attempted[1] > kwargs.get("Rc_max", 1e5):
                result["status"] = "outer_radius_limit_fallback"
                return result
            if not np.all(np.isfinite([rc, psi, omega_outer])) or rc <= 0:
                result["status"] = "invalid_outer_parameters"
                return result
            result["rc_error_pct"] = 100 * abs(rc * 1000 / truth["rc_m"] - 1)
            result.update(success=True, status="ok")
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            result.update(status="numerical_exception", optimizer_message=f"{type(exc).__name__}: {exc}")
        finally:
            result["warning_messages"] = " | ".join(sorted({str(w.message) for w in caught}))
    return result


def _validate(config):
    if config.replicates < 1 or not config.methods or len(set(config.methods)) != len(config.methods):
        raise ValueError("Require positive replicates and distinct methods.")
    if not set(config.methods) <= set(NOTEBOOKS):
        raise ValueError("Only SOLO, DOPPIO and LATTE are supported.")
    if not config.sample_counts or len(set(config.sample_counts)) != len(config.sample_counts):
        raise ValueError("Require distinct sample counts.")
    for fractions in (config.noise_fractions, config.background_fractions):
        if 0 not in fractions or any(not np.isfinite(f) or f < 0 for f in fractions):
            raise ValueError("Fractions must be finite/nonnegative and include zero.")


def _attach_metadata(frame, config, backend, truths, experiment):
    serial_truth = {method: {k: v.tolist() if isinstance(v, np.ndarray) else
                           (None if v is None else float(v)) for k, v in truth.items()}
                    for method, truth in truths.items()}
    frame.attrs.update(experiment=experiment, config=asdict(config), truth_by_method=serial_truth,
                       source_sha256=backend.sources,
                       versions=dict(python=platform.python_version(), numpy=np.__version__,
                                     pandas=pd.__version__, scipy=scipy.__version__))
    return frame


def run_noise_sampling(config=Config()):
    _validate(config)
    backend = load_backend(config)
    truths = {method: truth_parameters(backend, method, config.all_circular) for method in config.methods}
    rows = []
    for method in config.methods:
        truth = truths[method]
        for n in config.sample_counts:
            geometry = make_geometry(method, n, config.seed)
            clean = velocity(backend, geometry["xy"], truth)
            for replicate in range(config.replicates):
                # Method/count/replicate keys make streams invariant to loop order.
                rng = np.random.default_rng(np.random.SeedSequence(
                    [config.seed, 31, ("DOPPIO", "LATTE", "SOLO").index(method), n, replicate]))
                noise = rng.standard_normal(clean.shape)
                for fraction in config.noise_fractions:
                    if fraction == 0 and replicate != 0:
                        continue
                    uv = clean + fraction * truth["u_ref_m_s"] * noise
                    fitted = fit_observations(backend, method, geometry, uv, truth)
                    rows.append(dict(method=method, n_samples=n, noise_fraction=fraction,
                                     truth_shape="circular" if method == "SOLO" or config.all_circular else "elliptical",
                                     u_ref_m_s=truth["u_ref_m_s"],
                                     background_fraction=0., background_direction="none",
                                     replicate=replicate, **fitted))
    return _attach_metadata(pd.DataFrame(rows), config, backend, truths, "noise_sampling")


def run_background_flow(config=Config()):
    _validate(config)
    backend = load_backend(config)
    truths = {method: truth_parameters(backend, method, config.all_circular) for method in config.methods}
    rows = []
    for method in config.methods:
        truth = truths[method]
        geometry = make_geometry(method, config.background_n, config.seed)
        clean = velocity(backend, geometry["xy"], truth)
        for direction, unit in (("eastward", np.array([1., 0.])),
                                ("northward", np.array([0., 1.]))):
            for fraction in config.background_fractions:
                background = fraction * truth["u_ref_m_s"] * unit
                fitted = fit_observations(backend, method, geometry, clean + background, truth)
                rows.append(dict(method=method, n_samples=config.background_n,
                                 truth_shape="circular" if method == "SOLO" or config.all_circular else "elliptical",
                                 u_ref_m_s=truth["u_ref_m_s"],
                                 noise_fraction=0., background_fraction=fraction,
                                 background_direction=direction, background_u_m_s=background[0],
                                 background_v_m_s=background[1], replicate=0, **fitted))
    frame = pd.DataFrame(rows)
    for (_, _), group in frame.groupby(["method", "background_direction"]):
        baseline = group.loc[group.background_fraction == 0].iloc[0]
        frame.loc[group.index, "baseline_success"] = bool(baseline.success)
        for index, row in group.iterrows():
            good = row.success and baseline.success
            frame.loc[index, "incremental_centre_km"] = (
                np.hypot(row.xc_m - baseline.xc_m, row.yc_m - baseline.yc_m) / 1000 if good else np.nan)
            for component in ("xc_m", "yc_m"):
                frame.loc[index, f"incremental_{component}"] = row[component] - baseline[component] if good else np.nan
            for parameter in ("rc", "alpha"):
                column = "rc_m" if parameter == "rc" else "alpha"
                frame.loc[index, f"incremental_{parameter}_pct"] = (
                    100 * (row[column] / baseline[column] - 1)
                    if good and not (row.method == "SOLO" and parameter == "alpha") else np.nan)
    return _attach_metadata(frame, config, backend, truths, "background_flow")


def summarize(results):
    keys = ["method", "n_samples", "noise_fraction"] if results.attrs["experiment"] == "noise_sampling" else [
        "method", "background_direction", "background_fraction"]
    counts = results.groupby(keys, sort=False).agg(
        attempted=("success", "size"), successful=("success", "sum"),
        inner_converged=("inner_optimizer_success", "sum"),
        outer_converged=("outer_optimizer_success", "sum"))
    valid = results.loc[results.success].groupby(keys)[list(METRICS)].median()
    return counts.join(valid).reset_index()


def interpretation(results):
    """Generate factual Markdown from this execution, including baseline bias."""
    lines = []
    for method, group in results.groupby("method", sort=False):
        condition = "noise_fraction" if results.attrs["experiment"] == "noise_sampling" else "background_fraction"
        baseline = group.loc[(group[condition] == 0) & (group.n_samples == group.n_samples.max())]
        if len(baseline) and baseline.success.all():
            b = baseline.iloc[0]
            lines.append(f"- {method}, N={int(b.n_samples)} clean baseline: centre {b.centre_error_km:.3f} km; "
                         f"radius {b.rc_error_pct:.2f}% absolute error" + (
                             "; shape fixed circular (not estimated)." if method == "SOLO" else
                             f"; aspect ratio {b.alpha_error_pct:.2f}% absolute error."))
        else:
            lines.append(f"- {method}: the adequately sampled clean baseline failed; incremental biases require caution.")
    failures = results.loc[~results.success, "status"].value_counts().to_dict()
    lines.append(f"- Accepted {int(results.success.sum())}/{len(results)} fits; failure statuses: {failures or 'none'}. "
                 "Acceptance requires solver convergence and physical Q, with no accuracy threshold.")
    if results.attrs["experiment"] == "noise_sampling":
        high = results.loc[results.noise_fraction == results.noise_fraction.max()]
        details = []
        for method, group in high.groupby("method", sort=False):
            group = group.loc[(group.n_samples == group.n_samples.max()) & group.success]
            if len(group):
                details.append(f"{method} median centre {group.centre_error_km.median():.2f} km, "
                               f"radius {group.rc_error_pct.median():.2f}%")
        lines.append(f"- At the largest tested noise and N, among accepted fits: {'; '.join(details) or 'no accepted fits'}. "
                     "Bands show 10th–90th simulation percentiles, not confidence intervals.")
        limitation = ("These are fixed-geometry, independent component-noise perturbations; sampling-location variability, "
                      "correlated errors and observational error statistics were not tested.")
    else:
        details = []
        high = results.loc[results.background_fraction == results.background_fraction.max()]
        for _, row in high.iterrows():
            details.append(f"{row.method}/{row.background_direction}: " + (
                f"{row.centre_error_km:.2f} km total, {row.incremental_centre_km:.2f} km incremental centre displacement"
                if row.success else row.status))
        lines.append("- At the largest tested background: " + "; ".join(details) + ".")
        limitation = "These deterministic uniform-current tests do not test shear, noise, time dependence or real-data robustness."
    extreme_note = ""
    if results.attrs["config"].get("all_circular"):
        accepted = results.loc[results.success]
        if len(accepted):
            largest = accepted.loc[accepted.rc_m.idxmax()]
            true_rc = results.attrs["truth_by_method"][largest.method]["rc_m"]
            if largest.rc_m > 10 * true_rc:
                extreme_note = (f"\n\nThe largest circular-case fitted radius is {largest.rc_m/1000:,.0f} km "
                                f"({largest.method}; background fraction {largest.background_fraction:g}, "
                                f"{largest.background_direction}; noise fraction {largest.noise_fraction:g}). "
                                "This extreme estimate is retained: solver convergence does not establish "
                                "useful radius recovery. Radius plots use a logarithmic scale beyond 1% error.")
    return "\n".join(lines) + "\n\n" + limitation + (
        " The original inner-core approximations cause nonzero clean bias, especially LATTE shape bias; "
        "the simulations are sensitivity tests, not exact Gaussian-parameter estimators.") + extreme_note


def save_results(results, directory, stem):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    results.to_csv(directory / f"{stem}_trials.csv", index=False)
    summarize(results).to_csv(directory / f"{stem}_summary.csv", index=False)
    metadata = dict(results.attrs)
    metadata["helper_sha256"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__), Path(__file__).with_name("rc3_plotting.py"))}
    (directory / f"{stem}_provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
