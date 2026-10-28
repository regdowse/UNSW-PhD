"""Circular RC3 experiments; source inversions are loaded unchanged.

Public quantities use SI. The adapter alone converts to the manuscript notebooks'
km-coordinate / m/s numerical conventions. Configuration lives in rc3_config.json.
"""
from __future__ import annotations
import ast
from contextlib import contextmanager, redirect_stdout
from dataclasses import dataclass
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

METRICS = ("centre_error_km", "omega_error_pct", "alpha_error_pct", "rc_error_pct")


@dataclass(frozen=True)
class Config:
    values: dict
    path: Path

    @property
    def source_root(self):
        return self.values["source_root"]

    @property
    def methods(self):
        return tuple(self.values["methods"])

    @property
    def regimes(self):
        return tuple(self.values["sampling"]["regimes"])

    @property
    def seed(self):
        return self.values["seed"]

    @property
    def output(self):
        return self.path.parent / self.values["output_directory"]

    @property
    def supplementary_output(self):
        return self.path.parent / self.values["supplementary_output_directory"]


def load_config(path=None):
    path = Path(path or Path(__file__).with_name("rc3_config.json")).resolve()
    config = Config(json.loads(path.read_text()), path)
    v = config.values
    if config.methods != ("SOLO", "DOPPIO", "LATTE"):
        raise ValueError("This design requires SOLO, DOPPIO and LATTE in that order.")
    if config.regimes != ("LOW", "MEDIUM", "HIGH"):
        raise ValueError("Expected LOW, MEDIUM, HIGH sampling regimes.")
    if v["noise"]["replicates"] < 1:
        raise ValueError("Require a positive Monte Carlo replicate count.")
    if not np.array_equal(v["truth"]["q"], np.eye(2)) or not np.array_equal(v["truth"]["centre_m"], [0, 0]):
        raise ValueError("The revised design requires circular truth centred at zero.")
    for experiment in ("noise", "background"):
        fractions = v[experiment]["fractions"]
        if 0 not in fractions or len(set(fractions)) != len(fractions) or any(f < 0 or not np.isfinite(f) for f in fractions):
            raise ValueError("Perturbation fractions must be distinct, finite, nonnegative and include zero.")
    return config


def _definitions(text, names, filename, namespace):
    """Compile the exact original top-level definitions, without rewriting ASTs."""
    tree = ast.parse(text, filename=filename)
    tree.body = [node for node in tree.body
                 if isinstance(node, ast.FunctionDef) and node.name in names]
    exec(compile(tree, filename, "exec"), namespace)
    return {node.name for node in tree.body}


def load_backend(config):
    root = Path(config.source_root).expanduser().resolve()
    utils_path = root / config.values["source_utilities"]
    utils = dict(np=np, pd=pd)
    required = {"unit_det", "axis_ratio_from_Q", "tangential_velocity",
                "psi_params", "gaussian_vel_reconstruction"}
    found = _definitions(utils_path.read_text(), required, str(utils_path), utils)
    if found != required:
        raise RuntimeError(f"Missing original utilities: {required - found}")
    sources = {str(utils_path): hashlib.sha256(utils_path.read_bytes()).hexdigest()}
    methods = {}
    for method in config.methods:
        path = root / "MOCA_PAPER/final_ESP_paper_notebooks" / config.values["source_notebooks"][method]
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


def truth_parameters(backend, config):
    spec = config.values["truth"]
    q = np.asarray(spec["q"], float)
    rc, psi = float(spec["rc_m"]), float(spec["psi0_m2_s"])
    if rc <= 0 or not np.isfinite(psi):
        raise ValueError("Invalid Gaussian truth parameters.")
    truth = dict(xc_m=float(spec["centre_m"][0]), yc_m=float(spec["centre_m"][1]),
                 rc_m=rc, psi0_m2_s=psi, q=q, alpha=float(backend.axis_ratio_from_Q(q)),
                 omega_s=-2 * psi / rc**2, major_axis_deg=None)
    truth["vorticity_s"] = truth["omega_s"] * np.trace(q)
    eigenvalues, eigenvectors = np.linalg.eigh(q)
    peak = eigenvectors[:, -1] * rc / np.sqrt(2 * eigenvalues[-1])
    peak += [truth["xc_m"], truth["yc_m"]]
    truth["u_ref_m_s"] = float(np.linalg.norm(velocity(backend, peak[None], truth)))
    return truth


def velocity(backend, xy_m, truth):
    q = truth["q"]
    u, v, _, _ = backend.gaussian_vel_reconstruction(
        truth["xc_m"], truth["yc_m"], q[0, 0], q[0, 1], q[1, 1],
        truth["rc_m"], truth["psi0_m2_s"], X=xy_m[:, 0], Y=xy_m[:, 1])
    return np.column_stack((u, v))


def _sampling_signature(config):
    v = config.values
    return hashlib.sha256(json.dumps({k: v[k] for k in ("seed", "truth", "sampling")}, sort_keys=True).encode()).hexdigest()


def _bank_to_json(bank):
    return {method: {regime: {key: value.tolist() for key, value in geometry.items()}
                     for regime, geometry in regimes.items()} for method, regimes in bank.items()}


def _bank_from_json(data):
    return {method: {regime: {key: np.asarray(value, dtype=float if key == "xy" else int)
                             for key, value in geometry.items()}
                     for regime, geometry in regimes.items()} for method, regimes in data.items()}


def outer_count(geometry, truth):
    distance = np.linalg.norm(geometry["xy"] - [truth["xc_m"], truth["yc_m"]], axis=1)
    return int(np.count_nonzero(distance > truth["rc_m"] / np.sqrt(2)))


def _point_set(xy):
    return set(map(tuple, np.asarray(xy)))


def _assert_d4_symmetry(xy):
    """Sign/swap transforms preserve coordinate values exactly, even after JSON."""
    points = _point_set(xy)
    assert _point_set(np.column_stack((-xy[:, 1], xy[:, 0]))) == points, "Missing 90-degree counterpart"
    assert _point_set(xy * [-1, 1]) == points, "Missing reflection across y axis"
    assert _point_set(xy * [1, -1]) == points, "Missing reflection across x axis"


def _draw_open_radii(rng, n, lower, upper):
    """Uniform by stratum area, excluding floating-point boundary values."""
    radius = np.empty(n)
    pending = np.arange(n)
    while len(pending):
        draws = np.sqrt(lower**2 + rng.random(len(pending)) * (upper**2-lower**2))
        good = (draws > lower) & (draws < upper)
        radius[pending[good]] = draws[good]
        pending = pending[~good]
    return radius


def _symmetric_latte_clouds(config, truth):
    """Nested D4 orbits with exact stratum counts; no boundary anchors."""
    spec = config.values["sampling"]["latte"]
    if spec["distribution"] != "d4_stratified_uniform_area":
        raise ValueError("Expected the specified D4 stratified LATTE construction.")
    peak = truth["rc_m"] / np.sqrt(2)
    extent = config.values["sampling"]["extent_rc"] * truth["rc_m"]
    pools = {}
    for stratum, lower, upper, stream in (("inner", 0., peak, 0), ("outer", peak, extent, 1)):
        rng = np.random.default_rng(np.random.SeedSequence([config.seed, spec["stream_key"], stream]))
        n_eight = max(spec["groups"][r][stratum+"_eight"] for r in config.regimes)
        n_axis = max(spec["groups"][r][stratum+"_axis"] for r in config.regimes)
        radii = _draw_open_radii(rng, n_eight+n_axis, lower, upper)
        theta = rng.uniform(0, np.pi/4, n_eight)
        # Reject only degenerate axes/diagonal seeds, never fitted outcomes.
        degenerate = (theta <= 0) | (theta >= np.pi/4)
        while degenerate.any():
            theta[degenerate] = rng.uniform(0, np.pi/4, degenerate.sum())
            degenerate = (theta <= 0) | (theta >= np.pi/4)
        groups = []
        for radius, angle in zip(radii[:n_eight], theta):
            x, y = radius*np.cos(angle), radius*np.sin(angle)
            assert x > y > 0, "Generic eight-point seed must avoid axes and diagonals"
            groups.append(np.array([[x,y],[-x,-y],[y,-x],[-y,x],[-x,y],[x,-y],[-y,-x],[y,x]]))
        axis_groups = [np.array([[r,0],[0,r],[-r,0],[0,-r]]) for r in radii[n_eight:]]
        pools[stratum] = (groups, axis_groups)
    clouds = {}
    for regime in config.regimes:
        counts = spec["groups"][regime]
        pieces = [np.zeros((1, 2))]
        for stratum in ("inner", "outer"):
            eight, axis = pools[stratum]
            pieces.extend(eight[:counts[stratum+"_eight"]])
            pieces.extend(axis[:counts[stratum+"_axis"]])
        clouds[regime] = np.vstack(pieces)
    return clouds


def validate_geometry(bank, config, truth):
    sampling = config.values["sampling"]
    extent = sampling["extent_rc"] * truth["rc_m"]
    for method in config.methods:
        for regime in config.regimes:
            g = bank[method][regime]
            xy = g["xy"]
            n = sampling["regimes"][regime][method]
            expected = 2 * n - 1 if method == "DOPPIO" else n
            assert len(xy) == expected == len(np.unique(xy, axis=0)), (method, regime, "unique count")
            assert np.isfinite(xy).all()
            assert outer_count(g, truth) > 0
            if method == "LATTE":
                radius = np.linalg.norm(xy, axis=1)
                peak = truth["rc_m"] / np.sqrt(2)
                assert np.count_nonzero(np.all(xy == 0, axis=1)) == 1
                assert np.all(radius < extent), "LATTE points must be strictly inside Rmax"
                assert not np.any(radius == peak), "LATTE points cannot be on the stratum boundary"
                _assert_d4_symmetry(xy)
                counts = sampling["latte"]["groups"][regime]
                outer = outer_count(g, truth)
                assert (len(xy)-outer, outer) == (counts["inner_count"], counts["outer_count"])
                assert outer == outer_count(bank["DOPPIO"][regime], truth)
            else:
                paths = (g["horizontal"], g["vertical"]) if method == "DOPPIO" else (np.arange(n),)
                for axis, path in enumerate(paths):
                    np.testing.assert_allclose(xy[path, axis], np.linspace(-extent, extent, n), atol=1e-9)
                    np.testing.assert_equal(xy[path, 1-axis], 0)
                if method == "DOPPIO":
                    common = np.intersect1d(*paths)
                    assert len(common) == 1
                    np.testing.assert_equal(xy[common[0]], [0, 0])
                expected_hash = sampling["reference_transect_sha256"][method][regime]
                assert hashlib.sha256(xy.tobytes()).hexdigest() == expected_hash, "Transect changed from preceding revision"
    for small, large in zip(config.regimes[:-1], config.regimes[1:]):
        a, b = bank["LATTE"][small]["xy"], bank["LATTE"][large]["xy"]
        assert _point_set(a).issubset(_point_set(b)), "LATTE density sets must be nested"


def load_geometry_bank(config, truth):
    """One immutable coordinate manifest shared by both experiments and plots.

    D4 seed-pool prefixes give nested point sets with exact inner/outer counts.
    Reflections and quarter turns are applied only to positions, never noise.
    """
    config.output.mkdir(parents=True, exist_ok=True)
    signature = _sampling_signature(config)
    path = config.output / f"rc3_sampling_coordinates_{signature[:12]}.json"
    if path.exists():
        saved = json.loads(path.read_text())
        assert saved["sampling_signature"] == signature
        bank = _bank_from_json(saved["coordinates_m"])
    else:
        sampling = config.values["sampling"]
        extent = sampling["extent_rc"] * truth["rc_m"]
        clouds = _symmetric_latte_clouds(config, truth)
        bank = {method: {} for method in config.methods}
        for regime in config.regimes:
            for method in config.methods:
                n = sampling["regimes"][regime][method]
                if method == "LATTE":
                    geometry = dict(xy=clouds[regime].copy())
                else:
                    if n < 5 or n % 2 != 1:
                        raise ValueError("Transect counts must be odd integers >= 5.")
                    axis = np.linspace(-extent, extent, n)
                    horizontal = np.column_stack((axis, np.zeros(n)))
                    if method == "SOLO":
                        geometry = dict(xy=horizontal)
                    else:
                        vertical = np.column_stack((np.zeros(n), axis))
                        xy, inverse = np.unique(np.vstack((horizontal, vertical)), axis=0, return_inverse=True)
                        geometry = dict(xy=xy, horizontal=inverse[:n], vertical=inverse[n:])
                bank[method][regime] = geometry
        validate_geometry(bank, config, truth)
        radii = dict(Rc=truth["rc_m"], Rpeak=truth["rc_m"]/np.sqrt(2), Rmax=extent)
        path.write_text(json.dumps(dict(sampling_signature=signature, sampling_seed=config.seed,
            reference_radii_m=radii, coordinates_m=_bank_to_json(bank)), indent=2) + "\n")
    validate_geometry(bank, config, truth)
    return bank, path


def prepare_experiment(config):
    backend = load_backend(config)
    truth = truth_parameters(backend, config)
    bank, manifest = load_geometry_bank(config, truth)
    return SimpleNamespace(config=config, backend=backend, truth=truth, geometries=bank,
                           manifest=manifest, clean_fits={})


def sampling_summary(context):
    c, t = context.config, context.truth
    return pd.DataFrame([dict(method=m, density=r,
        n_per_transect=(c.values["sampling"]["regimes"][r][m] if m == "DOPPIO" else np.nan),
        n_unique=len(context.geometries[m][r]["xy"]),
        n_inner_true=len(context.geometries[m][r]["xy"])-outer_count(context.geometries[m][r], t),
        n_outer_true=outer_count(context.geometries[m][r], t))
        for m in c.methods for r in c.regimes])


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


def fit_observations(backend, method, geometry, uv, truth, config):
    """Original km-scale inversion, with explicit convergence/physical statuses."""
    settings = config.values["inversion"]
    unit = settings["metres_per_backend_length_unit"]
    xy = geometry["xy"] / unit
    result = dict(success=False, inner_success=False, status="unattempted",
                  inner_optimizer_success=False, outer_optimizer_success=False,
                  inner_nfev=np.nan, n_unique=len(xy), warning_messages="",
                  shape_estimated=method != "SOLO",
                  optimizer_message="", n_outer_fit=0,
                  **{key: np.nan for key in (*METRICS, "xc_m", "yc_m", "rc_m",
                     "alpha", "q11", "q12", "q22", "q_det", "q_min_eigenvalue",
                     "omega_inner_s", "omega_s", "psi0_m2_s", "vorticity_s",
                     "dx_km", "dy_km", "attempted_rc_m", "attempted_omega_s", "outer_omega_backend_units")})
    functions = backend.methods[method]
    with warnings.catch_warnings(record=True) as caught, redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        try:
            optimizer = "curve_fit" if method == "SOLO" else "least_squares"
            with _observe_optimizer(optimizer) as calls:
                if method == "SOLO":
                    fitted = functions.solo(xy[:, 0] - xy[0, 0], uv[:, 0], uv[:, 1], core_thresh=settings["solo_core_thresh_km"], plot=False)
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
            result.update(xc_m=unit * xc, yc_m=unit * yc, omega_inner_s=omega / unit)
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
            if (not np.allclose(q, q.T, atol=settings["q_symmetry_atol"]) or result["q_min_eigenvalue"] <= 0
                    or not np.isclose(result["q_det"], 1., rtol=settings["q_det_rtol"], atol=settings["q_det_atol"])):
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
            result["n_outer_fit"] = int(good.sum())
            if good.sum() < 3 or frame.loc[good, "rho2"].nunique() < 3:
                result["status"] = "insufficient_outer_samples"
                return result
            # Preserve the exact original call choices (limits are in km).
            kwargs = dict(Rc_max=settings["outer_rc_max_km"][method], maxfev=settings["outer_maxfev"])
            with _observe_optimizer("curve_fit") as calls:
                rc, psi, omega_outer = functions.out_core_param_fit(
                    frame.rho2, frame.Qr, frame.vt, Omega0=omega, plot=False, **kwargs)
            result.update(rc_m=unit * rc, psi0_m2_s=unit * psi,
                          omega_s=omega_outer / unit, outer_omega_backend_units=omega_outer,
                          vorticity_s=omega_outer / unit * np.trace(q))
            if not calls or not calls[-1]["success"]:
                result["status"] = "outer_optimizer_failed"
                result["optimizer_message"] = calls[-1].get("exception", "") if calls else "not called"
                return result
            result["outer_optimizer_success"] = True
            attempted = calls[-1]["parameters"]
            result["attempted_rc_m"] = unit * attempted[1]
            result["attempted_omega_s"] = attempted[0] / unit
            if not np.isfinite(attempted[1]) or attempted[1] > kwargs.get("Rc_max", 1e5):
                result["status"] = "outer_radius_limit_fallback"
                return result
            if not np.all(np.isfinite([rc, psi, omega_outer])) or rc <= 0:
                result["status"] = "invalid_outer_parameters"
                return result
            result["rc_error_pct"] = 100 * abs(rc * unit / truth["rc_m"] - 1)
            result.update(recovery_metrics(result, truth, method))
            result.update(success=True, status="ok")
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            result.update(status="numerical_exception", optimizer_message=f"{type(exc).__name__}: {exc}")
        finally:
            result["warning_messages"] = " | ".join(sorted({str(w.message) for w in caught}))
    return result


def recovery_metrics(fitted, truth, method):
    """Truth errors use FINAL outer-fit omega_s, never omega_inner_s."""
    centre = np.array([fitted["xc_m"], fitted["yc_m"]])
    return dict(centre_error_km=float(np.linalg.norm(centre - [truth["xc_m"], truth["yc_m"]]) / 1000),
                omega_error_pct=100 * abs((fitted["omega_s"] - truth["omega_s"]) / truth["omega_s"]),
                alpha_error_pct=(np.nan if method == "SOLO" else 100 * abs(fitted["alpha"] / truth["alpha"] - 1)),
                rc_error_pct=100 * abs(fitted["rc_m"] / truth["rc_m"] - 1))


def _clean_fit(context, method, regime):
    """Self-tests and experiment records share one deterministic fit per case."""
    key = (method, regime)
    if key not in context.clean_fits:
        geometry = context.geometries[method][regime]
        context.clean_fits[key] = fit_observations(context.backend, method, geometry,
            velocity(context.backend, geometry["xy"], context.truth), context.truth, context.config)
    return context.clean_fits[key].copy()


def verify_high_geometry(context):
    """Compare the actual background coordinates with saved noise coordinates.

    Standalone background execution can use the common immutable manifest;
    after the noise run, its trial provenance supplies an additional comparison.
    """
    c = context.config
    manifest_bank = _bank_from_json(json.loads(context.manifest.read_text())["coordinates_m"])
    sources = [(str(context.manifest), manifest_bank)]
    noise_path = c.output / "rc3_noise_sampling_provenance.json"
    if noise_path.exists():
        noise = json.loads(noise_path.read_text())
        if noise["sampling_signature"] != _sampling_signature(c):
            raise ValueError("Saved revised noise geometry differs from this configuration; run the noise notebook first or select a new output directory.")
        sources.append((str(noise_path), _bank_from_json(noise["sampling_coordinates_m"])))
    high = c.values["sampling"]["background_regime"]
    for _, bank in sources:
        for method in c.methods:
            for key, coordinates in context.geometries[method][high].items():
                np.testing.assert_array_equal(coordinates, bank[method][high][key])
    return [source for source, _ in sources]


def run_self_tests(context):
    c, t = context.config, context.truth
    validate_geometry(context.geometries, c, t)
    verified = verify_high_geometry(context)
    np.testing.assert_allclose(t["u_ref_m_s"], abs(t["omega_s"]) * t["rc_m"] / np.sqrt(2*np.e), rtol=1e-14)
    # Independent finite differences verify signs against the positive Gaussian psi.
    point = np.array([[.3*t["rc_m"], .2*t["rc_m"]]])
    step = 1.
    def psi(xy):
        delta = xy - [t["xc_m"], t["yc_m"]]
        return t["psi0_m2_s"] * np.exp(-np.einsum("ni,ij,nj->n", delta, t["q"], delta)/t["rc_m"]**2)
    expected_u = -(psi(point+[0, step])-psi(point-[0, step]))/(2*step)
    expected_v = (psi(point+[step, 0])-psi(point-[step, 0]))/(2*step)
    np.testing.assert_allclose(velocity(context.backend, point, t), np.column_stack((expected_u, expected_v)), rtol=1e-8)
    # Deliberately different preliminary/final Omegas expose scoring the wrong field.
    probe = dict(xc_m=t["xc_m"]+3000, yc_m=t["yc_m"]+4000,
                 omega_inner_s=9*t["omega_s"], omega_s=1.1*t["omega_s"],
                 alpha=1.2, rc_m=.9*t["rc_m"])
    errors = recovery_metrics(probe, t, "LATTE")
    np.testing.assert_allclose(list(errors.values()), [5., 10., 20., 10.], atol=1e-12)
    assert np.isnan(recovery_metrics(probe, t, "SOLO")["alpha_error_pct"])
    baselines = []
    high = c.values["sampling"]["background_regime"]
    for method in c.methods:
        fit = _clean_fit(context, method, high)
        assert fit["success"], (method, fit["status"], fit["optimizer_message"])
        assert np.isfinite([fit["centre_error_km"], fit["omega_error_pct"], fit["rc_error_pct"]]).all()
        np.testing.assert_equal(fit["omega_s"], fit["outer_omega_backend_units"] / c.values["inversion"]["metres_per_backend_length_unit"])
        baselines.append(dict(method=method, **{key: fit[key] for key in METRICS}))
    return dict(passed=True, high_geometry_checked_against=verified, clean_high_baselines=baselines,
                checks=["centred extents", "unique/outer counts", "one DOPPIO intersection",
                        "nested D4 LATTE strata", "exact LATTE inner/outer counts", "90-degree rotations and axis reflections",
                        "unchanged preceding-revision transects", "identical HIGH coordinates", "error definitions",
                        "final outer Omega and units", "streamfunction derivative signs", "common true U_ref", "finite accepted clean fits"])


def _record(context, method, regime, uv, **condition):
    c, t = context.config, context.truth
    g = context.geometries[method][regime]
    clean = condition["noise_fraction"] == 0 and condition["background_fraction"] == 0
    fit = (_clean_fit(context, method, regime) if clean else
           fit_observations(context.backend, method, g, uv, t, c))
    return dict(method=method, density=regime,
                n_samples=c.values["sampling"]["regimes"][regime][method],
                n_per_transect=(len(g["horizontal"]) if method == "DOPPIO" else np.nan),
                n_outer_true=outer_count(g, t), truth_shape="circular", u_ref_m_s=t["u_ref_m_s"],
                geometry_sha256=hashlib.sha256(g["xy"].tobytes()).hexdigest(), **condition, **fit)


def validate_trial_records(frame, context, experiment):
    """Check complete attempts, not just the subset accepted for quantiles."""
    c = context.config
    assert frame.status.notna().all()
    assert frame.success.eq(frame.status.eq("ok")).all()
    if experiment == "noise_sampling":
        keys = ["method", "density", "noise_fraction"]
        assert not frame.duplicated(keys+["replicate"]).any()
        expected = {(m, r, f): (1 if f == 0 else c.values["noise"]["replicates"])
                    for m in c.methods for r in c.regimes for f in c.values["noise"]["fractions"]}
    else:
        keys = ["method", "background_direction", "background_fraction"]
        expected = {(m, "baseline", 0.): 1 for m in c.methods}
        expected.update({(m, d, f): 1 for m in c.methods
                         for d in c.values["background"]["directions"][m]
                         for f in c.values["background"]["fractions"] if f > 0})
    assert frame.groupby(keys).size().to_dict() == expected, "Missing or duplicated trial conditions"
    assert frame.groupby(["method", "density"]).geometry_sha256.nunique().eq(1).all()
    accepted = frame.loc[frame.success]
    assert np.isfinite(accepted[["centre_error_km", "omega_error_pct", "rc_error_pct"]]).all().all()
    assert np.isfinite(accepted.loc[accepted.method != "SOLO", "alpha_error_pct"]).all()
    assert frame.loc[frame.method == "SOLO", "alpha_error_pct"].isna().all()


def _attach_metadata(frame, context, experiment, checks):
    validate_trial_records(frame, context, experiment)
    c = context.config
    truth = {key: value.tolist() if isinstance(value, np.ndarray) else value
             for key, value in context.truth.items()}
    import matplotlib
    frame.attrs.update(experiment=experiment, config=c.values, config_path=str(c.path),
                       config_sha256=hashlib.sha256(c.path.read_bytes()).hexdigest(), truth=truth,
                       sampling_signature=_sampling_signature(c),
                       sampling_manifest=str(context.manifest),
                       sampling_coordinates_m=_bank_to_json(context.geometries),
                       sampling_reference_radii_m=dict(Rc=context.truth["rc_m"],
                           Rpeak=context.truth["rc_m"]/np.sqrt(2),
                           Rmax=c.values["sampling"]["extent_rc"]*context.truth["rc_m"]),
                       sampling_seed=c.seed,
                       trial_validation="complete condition/replicate counts; explicit statuses; accepted full-fit metric validity; fixed coordinates",
                       sampling_counts=sampling_summary(context).to_dict(orient="records"),
                       source_sha256=context.backend.sources, self_tests=checks,
                       versions=dict(python=platform.python_version(), numpy=np.__version__,
                                     pandas=pd.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
                       definitions=dict(centre_error_km="norm(fitted centre - truth centre), km",
                           omega_error_pct="100*abs((FINAL outer Omega - true Omega)/true Omega)",
                           alpha_error_pct="100*abs(fitted alpha/true alpha-1); SOLO not estimated",
                           rc_error_pct="100*abs(fitted Rc/true Rc-1)",
                           n_samples="DOPPIO per transect; SOLO/LATTE unique total",
                           n_outer_true="unique observed locations with true rho > true Rc/sqrt(2)",
                           n_outer_fit="observations eligible for original outer fit; not N_outer_true",
                           units="locations/radius m; velocities m/s; Omega s^-1; psi m^2/s",
                           velocity_convention="u=-dpsi/dy; v=dpsi/dx; Omega=-2*psi0/Rc^2",
                           acceptance="source optimiser convergence, finite parameters, positive Rc, physical unit-det Q; no accuracy rejection"))
    return frame


def run_noise_sampling(context, checks=None):
    c, t = context.config, context.truth
    noise_config = c.values["noise"]
    rows = []
    for method in c.methods:
        for regime in c.regimes:
            g = context.geometries[method][regime]
            clean = velocity(context.backend, g["xy"], t)
            n = c.values["sampling"]["regimes"][regime][method]
            for replicate in range(noise_config["replicates"]):
                rng = np.random.default_rng(np.random.SeedSequence([
                    c.seed, noise_config["stream_key"], noise_config["method_stream_keys"][method], n, replicate]))
                normal = rng.standard_normal(clean.shape)
                for fraction in noise_config["fractions"]:
                    if fraction == 0 and replicate != 0:
                        continue
                    uv = clean + fraction * t["u_ref_m_s"] * normal
                    rows.append(_record(context, method, regime, uv,
                        noise_fraction=fraction, background_fraction=0., background_direction="none",
                        replicate=replicate))
    return _attach_metadata(pd.DataFrame(rows), context, "noise_sampling", checks)


def run_background_flow(context, checks=None):
    c, t = context.config, context.truth
    verified = verify_high_geometry(context)
    high = c.values["sampling"]["background_regime"]
    rows = []
    for method in c.methods:
        clean = velocity(context.backend, context.geometries[method][high]["xy"], t)
        # One baseline is shared by both SOLO directional series in the plot.
        rows.append(_record(context, method, high, clean, noise_fraction=0., background_fraction=0.,
                            background_direction="baseline", background_u_m_s=0., background_v_m_s=0., replicate=0))
        for direction in c.values["background"]["directions"][method]:
            unit = np.asarray(c.values["background"]["direction_vectors"][direction])
            np.testing.assert_allclose(np.linalg.norm(unit), 1)
            for fraction in c.values["background"]["fractions"]:
                if fraction == 0:
                    continue
                background = fraction * t["u_ref_m_s"] * unit
                rows.append(_record(context, method, high, clean + background,
                    noise_fraction=0., background_fraction=fraction, background_direction=direction,
                    background_u_m_s=background[0], background_v_m_s=background[1], replicate=0))
    frame = _attach_metadata(pd.DataFrame(rows), context, "background_flow", checks)
    frame.attrs["high_geometry_verified_against"] = verified
    return frame


def summarize(results):
    keys = (["method", "density", "noise_fraction"] if results.attrs["experiment"] == "noise_sampling"
            else ["method", "background_direction", "background_fraction"])
    counts = results.groupby(keys, sort=False).agg(attempted=("success", "size"), accepted=("success", "sum"))
    counts["acceptance_pct"] = 100 * counts.accepted / counts.attempted
    values = results.loc[results.success].groupby(keys)[list(METRICS)]
    summary = counts.join(values.median().add_suffix("_median"))
    if results.attrs["experiment"] == "noise_sampling":
        for percentile in results.attrs["config"]["plotting"]["percentiles"]:
            summary = summary.join(values.quantile(percentile/100).add_suffix(f"_p{percentile:g}"))
    return summary.reset_index()


def acceptance_table(results):
    if results.attrs["experiment"] != "noise_sampling":
        return results.groupby("method", sort=False).success.agg(accepted="sum", attempted="size")
    summary = summarize(results)
    summary["accepted/attempted"] = summary.accepted.astype(str) + "/" + summary.attempted.astype(str)
    table = summary.pivot(index=["method", "density"], columns="noise_fraction", values="accepted/attempted")
    order = pd.MultiIndex.from_product([results.attrs["config"]["methods"],
                                       results.attrs["config"]["sampling"]["regimes"]], names=["method", "density"])
    return table.reindex(index=order, columns=results.attrs["config"]["noise"]["fractions"])


def interpretation(results):
    clean_column = "noise_fraction" if results.attrs["experiment"] == "noise_sampling" else "background_fraction"
    clean = results.loc[results[clean_column] == 0]
    if "density" in clean:
        clean = clean.loc[clean.density == results.attrs["config"]["sampling"]["background_regime"]]
    lines = [f"- {row.method} clean HIGH baseline: centre {row.centre_error_km:.3f} km; "
             f"final Ω {row.omega_error_pct:.3g}%; Rc {row.rc_error_pct:.3g}% absolute error."
             if row.success else f"- {row.method} clean HIGH baseline failed: {row.status}."
             for _, row in clean.iterrows()]
    failures = results.loc[~results.success, "status"].value_counts().to_dict()
    lines.append(f"- Accepted {int(results.success.sum())}/{len(results)} full fits. Failure statuses: {failures or 'none'}. "
                 "Acceptance is numerical validity, not an accuracy or identifiability guarantee.")
    maximum = results.loc[results.success, "rc_error_pct"].max()
    if np.isfinite(maximum):
        lines.append(f"- Largest retained accepted-fit Rc error: {maximum:.4g}%; no error-based trimming is applied.")
    if results.attrs["experiment"] == "background_flow":
        alpha = results.loc[(results.method == "LATTE") & results.success, "alpha"]
        if len(alpha):
            lines.append(f"- LATTE fitted aspect ratio across the tested northward backgrounds: "
                         f"{alpha.min():.12g}–{alpha.max():.12g} (range {alpha.max()-alpha.min():.3g}). "
                         "This is an observation for these conditions, not a claim of general background-flow immunity.")
    return "\n".join(lines) + "\n\nThese fixed-geometry circular-Gaussian tests do not establish an unconditional method ranking or robustness to correlated errors, shear, temporal evolution, non-Gaussian eddies or platform-specific observation errors."


def save_results(results, context, figure_paths):
    c = context.config
    stem = "rc3_" + results.attrs["experiment"]
    results.to_csv(c.output / f"{stem}_trials.csv", index=False)
    summarize(results).to_csv(c.output / f"{stem}_summary.csv", index=False)
    sampling_summary(context).to_csv(c.output / f"{stem}_sampling.csv", index=False)
    metadata = dict(results.attrs)
    metadata["figure_paths"] = {key: str(value) for key, value in figure_paths.items()}
    metadata["helper_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__), Path(__file__).with_name("rc3_plotting.py"))}
    metadata["geometry_manifest_sha256"] = hashlib.sha256(context.manifest.read_bytes()).hexdigest()
    def json_safe(value):
        if isinstance(value, dict):
            return {k: json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(v) for v in value]
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return None
        return value
    (c.output / f"{stem}_provenance.json").write_text(json.dumps(json_safe(metadata), indent=2, allow_nan=False) + "\n")
