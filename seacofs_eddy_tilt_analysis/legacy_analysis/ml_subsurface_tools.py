"""Leakage-aware tools for predicting eddy tilt from surface information."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


STRUCTURE_FEATURES = ["Rc", "Omega"]
ENVIRONMENT_FEATURES = ["beta", "h"]
PROPAGATION_FEATURES = ["prop_east_km_day", "prop_north_km_day"]
ELLIPSE_FEATURES = ["ellipse_major_cos2", "ellipse_major_sin2"]
PV_MAG_FEATURES = ["PV_grad_mag"]
PV_DIRECTION_FEATURES = ["PV_grad_unit_east", "PV_grad_unit_north"]
FEATURES = (STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + PROPAGATION_FEATURES
            + ELLIPSE_FEATURES + PV_MAG_FEATURES + PV_DIRECTION_FEATURES)

# Each alternative removes or isolates a scientifically meaningful group.
FEATURE_SETS = {
    "full": FEATURES,
    "without_PV": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + PROPAGATION_FEATURES + ELLIPSE_FEATURES,
    "without_propagation": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + ELLIPSE_FEATURES + PV_MAG_FEATURES + PV_DIRECTION_FEATURES,
    "without_ellipse": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + PROPAGATION_FEATURES + PV_MAG_FEATURES + PV_DIRECTION_FEATURES,
    "without_beta": STRUCTURE_FEATURES + ["h"] + PROPAGATION_FEATURES + ELLIPSE_FEATURES + PV_MAG_FEATURES + PV_DIRECTION_FEATURES,
    "without_PV_magnitude": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + PROPAGATION_FEATURES + ELLIPSE_FEATURES + PV_DIRECTION_FEATURES,
    "without_PV_direction": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES + PROPAGATION_FEATURES + ELLIPSE_FEATURES + PV_MAG_FEATURES,
    "structure_environment": STRUCTURE_FEATURES + ENVIRONMENT_FEATURES,
    "beta_only": ["beta"],
}

MAGNITUDE_TARGET = ["LogTiltDis"]
DIRECTION_TARGETS = ["TiltUnitEast", "TiltUnitNorth"]


@dataclass(frozen=True)
class Configuration:
    family: str
    params_key: str
    params: Mapping
    feature_set: str


def angular_error_deg(observed, predicted):
    return np.abs((np.asarray(predicted) - np.asarray(observed) + 180.0) % 360.0 - 180.0)


def _major_axis_encoding(q11, q12, q22):
    q11, q12, q22 = map(lambda x: np.asarray(x, dtype=float), (q11, q12, q22))
    scale = np.hypot(q11 - q22, 2.0 * q12)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos2 = -(q11 - q22) / scale
        sin2 = -(2.0 * q12) / scale
    cos2[scale == 0] = np.nan
    sin2[scale == 0] = np.nan
    return sin2, cos2


def prepare_modelling_table(df, *, grid_rotation_deg=20.0):
    """Engineer predictors and separate magnitude/circular-direction targets."""

    required = {"Eddy", "Day", "Cyc", "TiltDis", "TiltDir", "beta", "h", "Omega",
                "Rc", "xc", "yc", "q11", "q12", "q22", "PV_grad_mag", "PV_grad_theta"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required modelling columns: {missing}")
    out = df.sort_values(["Eddy", "Day"]).copy()
    out["Cyc"] = out["Cyc"].astype("string")

    # Propagation requires the current and preceding surface positions. It is
    # therefore a trajectory feature rather than a single-snapshot feature.
    dt = out.groupby("Eddy")["Day"].diff().astype(float).where(lambda x: x > 0)
    out["prop_east_km_day"] = out.groupby("Eddy")["xc"].diff() / dt
    out["prop_north_km_day"] = out.groupby("Eddy")["yc"].diff() / dt

    sin2_grid, cos2_grid = _major_axis_encoding(out["q11"], out["q12"], out["q22"])
    rotation = np.deg2rad(2.0 * grid_rotation_deg)
    out["ellipse_major_sin2"] = sin2_grid * np.cos(rotation) + cos2_grid * np.sin(rotation)
    out["ellipse_major_cos2"] = cos2_grid * np.cos(rotation) - sin2_grid * np.sin(rotation)

    pv_theta = np.deg2rad(out["PV_grad_theta"].astype(float))
    out["PV_grad_unit_east"] = np.sin(pv_theta)
    out["PV_grad_unit_north"] = np.cos(pv_theta)
    out["PV_grad_east"] = out["PV_grad_mag"] * out["PV_grad_unit_east"]
    out["PV_grad_north"] = out["PV_grad_mag"] * out["PV_grad_unit_north"]
    out["PV_reference_theta"] = np.where(
        out["Cyc"].eq("AE"), (out["PV_grad_theta"] + 180.0) % 360.0, out["PV_grad_theta"]
    )

    theta = np.deg2rad(out["TiltDir"].astype(float))
    out["TiltUnitEast"] = np.sin(theta)
    out["TiltUnitNorth"] = np.cos(theta)
    out["LogTiltDis"] = np.log1p(out["TiltDis"].clip(lower=0))
    targets = out[["TiltDis", "TiltDir", "LogTiltDis", *DIRECTION_TARGETS]].astype(float)
    keep = np.isfinite(targets).all(axis=1) & out[["Eddy", "Cyc"]].notna().all(axis=1)
    out = out.loc[keep].copy()
    out["Eddy"] = out["Eddy"].astype(int)
    return out


def target_availability_summary(df):
    valid = df[["TiltDis", "TiltDir"]].notna().all(axis=1)
    rows = []
    for label, part in [("All", df), *list(df.groupby("Cyc", dropna=False))]:
        v = valid.loc[part.index]
        rows.append({"group": str(label), "rows": len(part), "valid_tilt_rows": int(v.sum()),
                     "valid_fraction": float(v.mean()), "eddies": int(part["Eddy"].nunique()),
                     "eddies_with_tilt": int(part.loc[v, "Eddy"].nunique())})
    return pd.DataFrame(rows).set_index("group")


def eddy_equal_weights(groups: Sequence):
    groups = pd.Series(groups)
    weights = 1.0 / groups.map(groups.value_counts()).to_numpy(dtype=float)
    return weights / weights.mean()


def _preprocessor(features, scale):
    steps = [("imputer", SimpleImputer(strategy="median", add_indicator=True))]
    if scale:
        steps.append(("scale", StandardScaler()))
    return ColumnTransformer([("numeric", Pipeline(steps), list(features))], verbose_feature_names_out=False)


def candidate_models():
    """Restrained candidates to keep nested CV computationally practical."""
    return [
        ("Ridge", "ridge_1", {"alpha": 1.0}),
        ("Ridge", "ridge_100", {"alpha": 100.0}),
        ("Gradient boosting", "boost_15", {"learning_rate": 0.05, "max_iter": 300,
         "max_leaf_nodes": 15, "min_samples_leaf": 75, "l2_regularization": 1.0}),
        ("Gradient boosting", "boost_regularised", {"learning_rate": 0.05, "max_iter": 300,
         "max_leaf_nodes": 15, "min_samples_leaf": 125, "l2_regularization": 10.0}),
    ]


def build_model(family, features, params, *, random_state=42):
    if family == "Ridge":
        return Pipeline([("preprocess", _preprocessor(features, True)), ("model", Ridge(**dict(params)))])
    if family == "Gradient boosting":
        estimator = MultiOutputRegressor(HistGradientBoostingRegressor(
            early_stopping=True, random_state=random_state, **dict(params)))
        return Pipeline([("preprocess", _preprocessor(features, False)), ("model", estimator)])
    raise ValueError(f"Unknown family: {family}")


def _targets(task):
    if task == "magnitude":
        return MAGNITUDE_TARGET
    if task == "direction":
        return DIRECTION_TARGETS
    raise ValueError(f"Unknown task: {task}")


def _fit(model, data, features, task):
    return model.fit(data[list(features)], data[_targets(task)],
                     model__sample_weight=eddy_equal_weights(data["Eddy"]))


def _predicted_direction(raw):
    raw = np.asarray(raw, dtype=float)
    return (np.degrees(np.arctan2(raw[:, 0], raw[:, 1])) + 360.0) % 360.0


def score_task(task, data, raw, *, minimum_tilt_km=5.0):
    if task == "magnitude":
        # Ridge returns shape (n_samples,) for a single target, whereas the
        # MultiOutputRegressor used by gradient boosting returns
        # (n_samples, 1). Flattening safely supports both representations.
        predicted = np.expm1(np.clip(np.asarray(raw, dtype=float).reshape(-1), 0.0, None))
        observed = data["TiltDis"].to_numpy(dtype=float)
        weights = eddy_equal_weights(data["Eddy"])
        scores = {
            "magnitude_MAE_km": float(np.mean(np.abs(predicted - observed))),
            "eddy_weighted_magnitude_MAE_km": float(np.average(np.abs(predicted - observed), weights=weights)),
            "magnitude_RMSE_km": float(np.sqrt(mean_squared_error(observed, predicted))),
            "magnitude_R2": float(r2_score(observed, predicted)),
            "magnitude_bias_km": float(np.average(predicted - observed, weights=weights)),
        }
        return scores, predicted
    predicted = _predicted_direction(raw)
    use = data["TiltDis"].to_numpy(dtype=float) >= minimum_tilt_km
    errors = angular_error_deg(data.loc[use, "TiltDir"], predicted[use])
    weights = eddy_equal_weights(data.loc[use, "Eddy"])
    scores = {"direction_rows": int(use.sum()), "mean_angular_error_deg": float(errors.mean()),
              "eddy_weighted_mean_angular_error_deg": float(np.average(errors, weights=weights)),
              "median_angular_error_deg": float(np.median(errors)),
              "within_30deg_fraction": float(np.mean(errors <= 30.0))}
    return scores, predicted


def _metric(task):
    return "eddy_weighted_magnitude_MAE_km" if task == "magnitude" else "eddy_weighted_mean_angular_error_deg"


def _configurations(feature_sets):
    for (family, key, params), feature_set in product(candidate_models(), feature_sets):
        yield Configuration(family, key, params, feature_set)


def select_configuration_grouped(data, task, *, feature_sets=FEATURE_SETS, n_splits=3, random_state=42):
    splits = list(GroupKFold(n_splits=n_splits).split(data, groups=data["Eddy"]))
    configs = list(_configurations(feature_sets))
    rows = []
    for config_id, config in enumerate(configs):
        features = feature_sets[config.feature_set]
        for fold, (fit_pos, val_pos) in enumerate(splits, 1):
            fit_data, val_data = data.iloc[fit_pos], data.iloc[val_pos]
            model = build_model(config.family, features, config.params, random_state=random_state + fold)
            fitted = _fit(model, fit_data, features, task)
            scores, _ = score_task(task, val_data, fitted.predict(val_data[features]))
            rows.append({"config_id": config_id, "fold": fold, "family": config.family,
                         "params_key": config.params_key, "feature_set": config.feature_set, **scores})
    results = pd.DataFrame(rows)
    winner = int(results.groupby("config_id")[_metric(task)].mean().idxmin())
    return configs[winner], results


def _baseline_predictions(task, train, validation):
    if task == "magnitude":
        value = float(train.groupby("Eddy")["TiltDis"].median().median())
        return {"Eddy-median baseline": np.full(len(validation), value)}
    per_eddy = train.groupby("Eddy")[DIRECTION_TARGETS].mean().mean()
    mean_direction = (np.degrees(np.arctan2(per_eddy.iloc[0], per_eddy.iloc[1])) + 360.0) % 360.0
    propagation = (np.degrees(np.arctan2(validation["prop_east_km_day"],
                                          validation["prop_north_km_day"])) + 360.0) % 360.0
    return {"Mean-direction baseline": np.full(len(validation), mean_direction),
            "Polarity-aligned PV baseline": validation["PV_reference_theta"].to_numpy(float),
            "Propagation-direction baseline": propagation.to_numpy(float)}


def _score_baseline(task, data, predicted, minimum_tilt_km):
    if task == "magnitude":
        observed = data["TiltDis"].to_numpy(float)
        weights = eddy_equal_weights(data["Eddy"])
        return {"magnitude_MAE_km": float(np.mean(np.abs(predicted - observed))),
                "eddy_weighted_magnitude_MAE_km": float(np.average(np.abs(predicted - observed), weights=weights)),
                "magnitude_RMSE_km": float(np.sqrt(np.mean((predicted - observed) ** 2))),
                "magnitude_R2": float(r2_score(observed, predicted)),
                "magnitude_bias_km": float(np.average(predicted - observed, weights=weights))}
    use = data["TiltDis"].to_numpy(float) >= minimum_tilt_km
    errors = angular_error_deg(data.loc[use, "TiltDir"], np.asarray(predicted)[use])
    weights = eddy_equal_weights(data.loc[use, "Eddy"])
    return {"direction_rows": int(use.sum()), "mean_angular_error_deg": float(errors.mean()),
            "eddy_weighted_mean_angular_error_deg": float(np.average(errors, weights=weights)),
            "median_angular_error_deg": float(np.median(errors)),
            "within_30deg_fraction": float(np.mean(errors <= 30.0))}


def nested_grouped_evaluation(data, task, *, feature_sets=FEATURE_SETS, outer_splits=5,
                              inner_splits=3, minimum_tilt_km=5.0, random_state=42):
    """Repeat model and feature selection inside every held-out eddy fold."""
    outer = GroupKFold(n_splits=outer_splits)
    scores, selections, predictions = [], [], []
    for outer_fold, (train_pos, test_pos) in enumerate(outer.split(data, groups=data["Eddy"]), 1):
        train, test = data.iloc[train_pos], data.iloc[test_pos]
        config, inner = select_configuration_grouped(
            train, task, feature_sets=feature_sets, n_splits=inner_splits,
            random_state=random_state + outer_fold * 100)
        features = feature_sets[config.feature_set]
        fitted = _fit(build_model(config.family, features, config.params,
                                  random_state=random_state + outer_fold), train, features, task)
        task_scores, predicted = score_task(
            task, test, fitted.predict(test[features]), minimum_tilt_km=minimum_tilt_km)
        scores.append({"outer_fold": outer_fold, "model": "Selected model", **task_scores})
        selections.append({"outer_fold": outer_fold, "family": config.family,
                           "params_key": config.params_key, "feature_set": config.feature_set,
                           "inner_score": inner.groupby("config_id")[_metric(task)].mean().min()})
        frame = test[["Eddy", "TiltDis", "TiltDir", "beta", "Rc", "PV_reference_theta"]].copy()
        frame["outer_fold"] = outer_fold
        if task == "magnitude":
            frame["predicted_magnitude"] = predicted
        else:
            frame["predicted_direction"] = predicted
            frame["angular_error"] = angular_error_deg(frame["TiltDir"], predicted)
        predictions.append(frame)
        for name, baseline in _baseline_predictions(task, train, test).items():
            valid_baseline = np.isfinite(np.asarray(baseline, dtype=float))
            if valid_baseline.any():
                baseline_data = test.iloc[np.flatnonzero(valid_baseline)]
                scores.append({"outer_fold": outer_fold, "model": name,
                               **_score_baseline(task, baseline_data,
                                                 np.asarray(baseline)[valid_baseline],
                                                 minimum_tilt_km)})
    return pd.DataFrame(scores), pd.DataFrame(selections), pd.concat(predictions)


def summarise_outer_scores(scores, task):
    secondary = "magnitude_R2" if task == "magnitude" else "within_30deg_fraction"
    return scores.groupby("model").agg(mean_score=(_metric(task), "mean"),
                                        fold_SD=(_metric(task), "std"),
                                        secondary_mean=(secondary, "mean")).sort_values("mean_score")


def consensus_configuration(selections):
    family, key, feature_set = Counter(zip(selections["family"], selections["params_key"],
                                            selections["feature_set"])).most_common(1)[0][0]
    params = next(p for fam, name, p in candidate_models() if fam == family and name == key)
    return Configuration(family, key, params, feature_set)


def feature_set_comparison(inner_results, task):
    metric = _metric(task)
    by_config = inner_results.groupby(["feature_set", "family", "params_key"])[metric].mean().reset_index()
    return by_config.loc[by_config.groupby("feature_set")[metric].idxmin()].sort_values(metric)


def assign_eddy_spatial_blocks(data, *, bins=4):
    centres = data.groupby("Eddy")[["xc", "yc"]].median()
    xbin = pd.qcut(centres["xc"], bins, labels=False, duplicates="drop")
    ybin = pd.qcut(centres["yc"], bins, labels=False, duplicates="drop")
    return data.join((xbin.astype(str) + "_" + ybin.astype(str)).rename("spatial_block"), on="Eddy")


def spatial_block_evaluation(data, task, config, *, n_splits=4, minimum_tilt_km=5.0):
    """Stress-test a fixed nested-CV consensus configuration across regions."""
    blocked = assign_eddy_spatial_blocks(data)
    features = FEATURE_SETS[config.feature_set]
    rows = []
    splitter = GroupKFold(n_splits=n_splits)
    for fold, (train_pos, test_pos) in enumerate(splitter.split(blocked, groups=blocked["spatial_block"]), 1):
        train, test = blocked.iloc[train_pos], blocked.iloc[test_pos]
        fitted = _fit(build_model(config.family, features, config.params, random_state=500 + fold),
                      train, features, task)
        task_scores, _ = score_task(task, test, fitted.predict(test[features]),
                                    minimum_tilt_km=minimum_tilt_km)
        rows.append({"spatial_fold": fold, **task_scores})
    return pd.DataFrame(rows)


def beta_magnitude_association(data, *, n_boot=1000, random_state=42):
    """Descriptive eddy-level beta/magnitude association with eddy bootstrap CI."""
    eddies = data.groupby("Eddy", as_index=False).agg(
        beta=("beta", "median"), tilt_magnitude=("TiltDis", "median"), Rc=("Rc", "median"),
        Omega=("Omega", "median"), h=("h", "median")).dropna()
    rho = float(spearmanr(eddies["beta"], eddies["tilt_magnitude"]).statistic)
    rng = np.random.default_rng(random_state)
    boot = []
    for _ in range(n_boot):
        sample = eddies.iloc[rng.integers(0, len(eddies), len(eddies))]
        boot.append(spearmanr(sample["beta"], sample["tilt_magnitude"]).statistic)
    return {"eddies": len(eddies), "spearman_rho": rho,
            "CI_low": float(np.nanpercentile(boot, 2.5)),
            "CI_high": float(np.nanpercentile(boot, 97.5))}, eddies


def propagation_confounding_summary(data):
    """Between-track versus within-track propagation diagnostic (not predictors)."""
    out = data.copy()
    for axis in ("east", "north"):
        source = f"prop_{axis}_km_day"
        out[f"track_mean_prop_{axis}"] = out.groupby("Eddy")[source].transform("mean")
        out[f"prop_{axis}_anomaly"] = out[source] - out[f"track_mean_prop_{axis}"]
    rows = []
    for name in ["track_mean_prop_east", "track_mean_prop_north",
                 "prop_east_anomaly", "prop_north_anomaly"]:
        valid = out[[name, "TiltDis"]].dropna()
        rows.append({"quantity": name,
                     "spearman_with_tilt_magnitude": spearmanr(valid[name], valid["TiltDis"]).statistic})
    return pd.DataFrame(rows).set_index("quantity"), out


def direction_performance_by_tilt(predictions, thresholds=(0, 5, 10, 20)):
    rows = []
    for threshold in thresholds:
        part = predictions[predictions["TiltDis"] >= threshold]
        errors = angular_error_deg(part["TiltDir"], part["predicted_direction"])
        rows.append({"minimum_tilt_km": threshold, "rows": len(part),
                     "median_angular_error_deg": np.median(errors),
                     "within_30deg_fraction": np.mean(errors <= 30)})
    return pd.DataFrame(rows).set_index("minimum_tilt_km")


def plot_magnitude_oof(predictions, *, title):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hexbin(predictions["TiltDis"], predictions["predicted_magnitude"],
                   gridsize=45, mincnt=1, cmap="viridis")
    limit = np.nanpercentile(np.r_[predictions["TiltDis"], predictions["predicted_magnitude"]], 99)
    axes[0].plot([0, limit], [0, limit], "k--", lw=1)
    axes[0].set(xlim=(0, limit), ylim=(0, limit), xlabel="Observed magnitude (km)",
                ylabel="Out-of-fold predicted magnitude (km)")
    residual = predictions["predicted_magnitude"] - predictions["TiltDis"]
    axes[1].scatter(predictions["TiltDis"], residual, s=4, alpha=0.15)
    axes[1].axhline(0, color="k", ls="--", lw=1)
    axes[1].set(xlabel="Observed magnitude (km)", ylabel="Residual (km)")
    fig.suptitle(title); fig.tight_layout(); return fig, axes


def plot_direction_oof(predictions, *, title):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(predictions["angular_error"], bins=np.arange(0, 185, 5), color="slateblue")
    axes[0].axvline(90, color="k", ls="--", lw=1)
    axes[0].set(xlim=(0, 180), xlabel="Angular error (degrees)", ylabel="Count")
    axes[1].scatter(predictions["TiltDis"], predictions["angular_error"], s=4, alpha=0.15)
    axes[1].set(xlabel="Observed magnitude (km)", ylabel="Angular error (degrees)", ylim=(0, 180))
    fig.suptitle(title); fig.tight_layout(); return fig, axes


def plot_beta_relationship(eddy_table, *, title, ylim=None):
    ordered = eddy_table.sort_values("beta").copy()
    ordered["beta_bin"] = pd.qcut(ordered["beta"], 10, duplicates="drop")
    binned = ordered.groupby("beta_bin", observed=True).agg(beta=("beta", "median"),
                                                             tilt=("tilt_magnitude", "median"))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(ordered["beta"], ordered["tilt_magnitude"], s=8, alpha=0.15, label="Eddy median")
    ax.plot(binned["beta"], binned["tilt"], "o-", color="crimson", label="Decile median")
    ax.set(xlabel="Beta", ylabel="Median tilt magnitude (km)", title=title)
    ax.legend(); fig.tight_layout(); 
    if ylim is not None:
        ax.set_ylim(ylim)
    return fig, ax


# ---------------------------------------------------------------------------
# Statistical association workflow
# ---------------------------------------------------------------------------

ASSOCIATION_DECOMPOSE_FEATURES = [
    "abs_beta", "Rc", "Omega", "h", "PV_grad_mag",
    "prop_east_km_day", "prop_north_km_day",
]

ASSOCIATION_GROUP_PATTERNS = {
    "structure": ["Rc_", "Omega_", "ellipse_major_"],
    "environment": ["abs_beta_", "h_"],
    "beta": ["abs_beta_"],
    "PV": ["PV_grad_mag_", "PV_grad_unit_"],
    "propagation": ["prop_east_", "prop_north_"],
    "ellipse": ["ellipse_major_"],
}


def prepare_association_table(data, *, features=ASSOCIATION_DECOMPOSE_FEATURES):
    """Create within/between covariates and standardise them within polarity.

    The returned scaling table converts standardised effect-curve axes back to
    physical units.  Between terms describe differences among eddies; within
    terms describe day-to-day departures within the same eddy.
    """

    out = data.copy()
    out["abs_beta"] = out["beta"].abs()
    out["log_tilt"] = np.log1p(out["TiltDis"].clip(lower=0))
    theta = np.deg2rad(out["TiltDir"].astype(float))
    out["tilt_unit_east"] = np.sin(theta)
    out["tilt_unit_north"] = np.cos(theta)

    for feature in features:
        mean = out.groupby("Eddy")[feature].transform("mean")
        out[f"{feature}_between"] = mean
        out[f"{feature}_within"] = out[feature] - mean

    scale_rows = []
    columns = [f"{feature}_{level}" for feature in features for level in ("between", "within")]
    for column in columns:
        mean = float(out[column].mean(skipna=True))
        sd = float(out[column].std(skipna=True))
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        out[f"{column}_z"] = (out[column] - mean) / sd
        scale_rows.append({"column": column, "mean": mean, "sd": sd})
    return out, pd.DataFrame(scale_rows).set_index("column")


def eddy_level_association_table(data):
    """One robust descriptive row per eddy."""

    work = data.copy()
    work["abs_beta"] = work["beta"].abs()
    theta = np.deg2rad(work["TiltDir"].astype(float))
    work["tilt_sin"] = np.sin(theta)
    work["tilt_cos"] = np.cos(theta)
    summary = work.groupby(["Cyc", "Eddy"], as_index=False).agg(
        tilt_magnitude=("TiltDis", "median"),
        tilt_sin=("tilt_sin", "mean"), tilt_cos=("tilt_cos", "mean"),
        abs_beta=("abs_beta", "median"), Rc=("Rc", "median"),
        Omega=("Omega", "median"), h=("h", "median"),
        PV_grad_mag=("PV_grad_mag", "median"),
        prop_east=("prop_east_km_day", "mean"),
        prop_north=("prop_north_km_day", "mean"),
        ellipse_cos2=("ellipse_major_cos2", "mean"),
        ellipse_sin2=("ellipse_major_sin2", "mean"),
    )
    summary["tilt_direction"] = (
        np.degrees(np.arctan2(summary["tilt_sin"], summary["tilt_cos"])) + 360.0
    ) % 360.0
    summary["direction_concentration"] = np.hypot(summary["tilt_sin"], summary["tilt_cos"])
    return summary


def cluster_bootstrap_spearman(data, x, y, *, cluster="Eddy", n_boot=1000, random_state=42):
    """Spearman association with complete-cluster bootstrap confidence interval."""

    clean = data[[cluster, x, y]].dropna()
    point = float(spearmanr(clean[x], clean[y]).statistic)
    clusters = clean[cluster].drop_duplicates().to_numpy()
    by_cluster = {key: part for key, part in clean.groupby(cluster, sort=False)}
    rng = np.random.default_rng(random_state)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        sampled = rng.choice(clusters, size=len(clusters), replace=True)
        sample = pd.concat([by_cluster[key] for key in sampled], ignore_index=True)
        boot[i] = spearmanr(sample[x], sample[y]).statistic
    return {
        "x": x, "y": y, "clusters": len(clusters), "spearman_rho": point,
        "CI_low": float(np.nanpercentile(boot, 2.5)),
        "CI_high": float(np.nanpercentile(boot, 97.5)),
    }


def eddy_level_spearman_table(eddy_table, predictors, *, response="tilt_magnitude",
                              n_boot=1000, random_state=42):
    rows = []
    for offset, predictor in enumerate(predictors):
        rows.append(cluster_bootstrap_spearman(
            eddy_table, predictor, response, cluster="Eddy", n_boot=n_boot,
            random_state=random_state + offset,
        ))
    return pd.DataFrame(rows).set_index("x").sort_values("spearman_rho")


def circular_mean_deg(angles):
    theta = np.deg2rad(np.asarray(angles, dtype=float))
    return float((np.degrees(np.arctan2(np.nanmean(np.sin(theta)), np.nanmean(np.cos(theta)))) + 360.0) % 360.0)


def circular_resultant_length(angles):
    theta = np.deg2rad(np.asarray(angles, dtype=float))
    return float(np.hypot(np.nanmean(np.sin(theta)), np.nanmean(np.cos(theta))))


def pv_relative_direction_summary(data):
    """Describe tilt offset from the polarity-specific PV reference by eddy."""

    out = data.copy()
    out["delta_PV_deg"] = (out["TiltDir"] - out["PV_reference_theta"] + 180.0) % 360.0 - 180.0
    rows = []
    for eddy, part in out.groupby("Eddy"):
        mean_delta = circular_mean_deg(part["delta_PV_deg"])
        rows.append({
            "Eddy": eddy, "mean_delta_PV_deg": (mean_delta + 180.0) % 360.0 - 180.0,
            "PV_alignment_concentration": circular_resultant_length(part["delta_PV_deg"]),
            "median_tilt_km": part["TiltDis"].median(),
        })
    return pd.DataFrame(rows)


def _association_group_terms():
    """Formula terms shared by magnitude and direction GEE models."""

    return {
        "structure": [
            "bs(Rc_between_z, df=5)", "bs(Rc_within_z, df=5)",
            "bs(Omega_between_z, df=5)", "bs(Omega_within_z, df=5)",
            "ellipse_major_cos2", "ellipse_major_sin2",
        ],
        "environment": [
            "bs(abs_beta_between_z, df=5)", "bs(abs_beta_within_z, df=5)",
            "bs(h_between_z, df=5)", "bs(h_within_z, df=5)",
        ],
        "PV": [
            "bs(PV_grad_mag_between_z, df=5)", "bs(PV_grad_mag_within_z, df=5)",
            "PV_grad_unit_east", "PV_grad_unit_north",
        ],
        "propagation": [
            "prop_east_km_day_between_z", "prop_east_km_day_within_z",
            "prop_north_km_day_between_z", "prop_north_km_day_within_z",
        ],
    }


def association_gee_formula(response, *, exclude_groups=()):
    terms = _association_group_terms()
    rhs = [term for group, group_terms in terms.items()
           if group not in set(exclude_groups) for term in group_terms]
    return f"{response} ~ " + " + ".join(rhs)


def _formula_columns(formula):
    """Return raw column names needed after removing formula functions."""

    import re
    response, rhs = formula.split("~", 1)
    columns = [response.strip()]
    for term in rhs.split("+"):
        term = term.strip()
        match = re.match(r"bs\(([^,]+)", term)
        columns.append(match.group(1).strip() if match else term)
    return list(dict.fromkeys(columns))


def fit_clustered_gee(data, response, *, exclude_groups=(), covariance="exchangeable",
                      maxiter=100):
    """Fit a spline GEE with eddy-clustered robust uncertainty.

    ``statsmodels`` is imported lazily so the predictive workflow does not
    require it.  Exchangeable covariance is the default; independence is a
    useful sensitivity analysis.
    """

    import statsmodels.formula.api as smf
    from statsmodels.genmod.cov_struct import Exchangeable, Independence
    from statsmodels.genmod.families import Gaussian

    formula = association_gee_formula(response, exclude_groups=exclude_groups)
    required = ["Eddy", *_formula_columns(formula)]
    fit_data = data[required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    cov = Exchangeable() if covariance == "exchangeable" else Independence()
    model = smf.gee(formula, groups="Eddy", data=fit_data,
                    cov_struct=cov, family=Gaussian())
    return model.fit(maxiter=maxiter), fit_data


def fit_magnitude_gee(data, **kwargs):
    return fit_clustered_gee(data, "log_tilt", **kwargs)


def fit_direction_gee(data, **kwargs):
    east, east_data = fit_clustered_gee(data, "tilt_unit_east", **kwargs)
    north, north_data = fit_clustered_gee(data, "tilt_unit_north", **kwargs)
    return {"east": east, "north": north}, east_data.index.intersection(north_data.index)


def gee_group_wald_table(result, *, group_patterns=ASSOCIATION_GROUP_PATTERNS):
    """Joint robust Wald tests for predeclared feature groups."""

    names = list(result.model.exog_names)
    rows = []
    for group, patterns in group_patterns.items():
        indices = [i for i, name in enumerate(names)
                   if any(pattern in name for pattern in patterns)]
        if not indices:
            continue
        restriction = np.zeros((len(indices), len(names)))
        restriction[np.arange(len(indices)), indices] = 1.0
        test = result.wald_test(restriction, scalar=True)
        rows.append({
            "group": group, "parameters_tested": len(indices),
            "wald_statistic": float(np.asarray(test.statistic).squeeze()),
            "p_value": float(np.asarray(test.pvalue).squeeze()),
        })
    return pd.DataFrame(rows).set_index("group").sort_values("p_value")


def direction_group_wald_table(results):
    east = gee_group_wald_table(results["east"])[["parameters_tested", "wald_statistic", "p_value"]]
    north = gee_group_wald_table(results["north"])[["parameters_tested", "wald_statistic", "p_value"]]
    return east.add_prefix("east_").join(north.add_prefix("north_"), how="outer")


def gee_effect_curve(result, data, feature_z, *, scaling=None, points=100,
                     response_transform=None):
    """Marginal association curve and robust pointwise 95% confidence band."""

    from patsy import build_design_matrices

    needed = [name for name in _formula_columns(result.model.formula) if name != result.model.endog_names]
    baseline = {column: float(data[column].median()) for column in needed}
    low, high = data[feature_z].quantile([0.02, 0.98])
    curve = pd.DataFrame([baseline] * points)
    curve[feature_z] = np.linspace(low, high, points)
    design = np.asarray(build_design_matrices(
        [result.model.data.design_info], curve, return_type="dataframe"
    )[0], dtype=float)
    predicted = design @ np.asarray(result.params, dtype=float)
    covariance = np.asarray(result.cov_params(), dtype=float)
    standard_error = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", design, covariance, design), 0.0))
    lower, upper = predicted - 1.96 * standard_error, predicted + 1.96 * standard_error
    if response_transform:
        predicted, lower, upper = map(response_transform, (predicted, lower, upper))
    curve["prediction"] = predicted
    curve["CI_low"] = lower
    curve["CI_high"] = upper
    curve["feature_z"] = curve[feature_z]
    original = feature_z.removesuffix("_z")
    if scaling is not None and original in scaling.index:
        curve["feature_value"] = curve["feature_z"] * scaling.loc[original, "sd"] + scaling.loc[original, "mean"]
    else:
        curve["feature_value"] = curve["feature_z"]
    return curve[["feature_value", "feature_z", "prediction", "CI_low", "CI_high"]]


def direction_effect_curve(results, data, feature_z, *, scaling=None, points=100):
    east = gee_effect_curve(results["east"], data, feature_z, scaling=scaling, points=points)
    north = gee_effect_curve(results["north"], data, feature_z, scaling=scaling, points=points)
    out = east[["feature_value", "feature_z"]].copy()
    out["predicted_east"] = east["prediction"]
    out["predicted_north"] = north["prediction"]
    out["predicted_direction"] = (
        np.degrees(np.arctan2(out["predicted_east"], out["predicted_north"])) + 360.0
    ) % 360.0
    out["predicted_concentration"] = np.hypot(out["predicted_east"], out["predicted_north"])
    return out


def plot_association_curve(curve, *, xlabel, ylabel, title):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(curve["feature_value"], curve["prediction"], color="crimson", lw=2)
    ax.fill_between(curve["feature_value"], curve["CI_low"], curve["CI_high"],
                    color="crimson", alpha=0.2, linewidth=0)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    fig.tight_layout(); return fig, ax


def plot_direction_effect_curve(curve, *, xlabel, title):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(curve["feature_value"], curve["predicted_direction"], color="darkorange", lw=2)
    axes[0].set(xlabel=xlabel, ylabel="Predicted direction (degrees)", ylim=(0, 360))
    axes[1].plot(curve["feature_value"], curve["predicted_concentration"], color="teal", lw=2)
    axes[1].set(xlabel=xlabel, ylabel="Predicted concentration")
    fig.suptitle(title); fig.tight_layout(); return fig, axes
