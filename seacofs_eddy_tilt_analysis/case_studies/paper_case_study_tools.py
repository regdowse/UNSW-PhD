"""Selection and plotting tools for paper-level eddy tilt case studies.

The signed-PV convention used here matches the current population notebooks:
AE tilt is expected along signed grad(PV), while CE tilt is expected opposite it.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PaperCaseConfig:
    smooth_window: int = 7
    min_periods: int = 5
    dominance_factor: float = 2.0
    min_planetary_depth_m: float = 3000.0
    max_topographic_depth_m: float = 2000.0
    min_tilt_distance_km: float = 0.0
    min_lifetime_days: float = 40.0
    min_directional_observations: int = 20
    min_regime_observations: int = 10
    min_sustained_run: int = 5
    min_reference_regime_fraction: float = 0.5
    min_topographic_regime_fraction: float = 0.5
    angle_tolerance_deg: float = 45.0
    min_transition_rotation_deg: float = 30.0


def signed_angle_difference(angle, reference):
    """Signed shortest difference ``angle - reference`` in [-180, 180)."""
    return (
        np.asarray(angle, dtype=float)
        - np.asarray(reference, dtype=float)
        + 180.0
    ) % 360.0 - 180.0


def circular_mean_deg(values) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    radians = np.deg2rad(values)
    return float(np.rad2deg(np.arctan2(
        np.mean(np.sin(radians)), np.mean(np.cos(radians))
    )) % 360.0)


def circular_resultant(values) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    radians = np.deg2rad(values)
    return float(np.hypot(
        np.mean(np.cos(radians)), np.mean(np.sin(radians))
    ))


def longest_true_run(values) -> int:
    longest = current = 0
    for value in np.asarray(values, dtype=bool):
        current = current + 1 if value else 0
        longest = max(longest, current)
    return int(longest)


def add_paper_case_diagnostics(
    df: pd.DataFrame,
    config: PaperCaseConfig = PaperCaseConfig(),
) -> pd.DataFrame:
    """Add regime and corrected polarity-aware directional diagnostics."""
    required = {
        "Eddy", "Day", "Cyc", "h", "TiltDis", "TiltDir",
        "PV_grad_theta", "topo_plan_ratio",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")
    if config.smooth_window < 1 or config.smooth_window % 2 == 0:
        raise ValueError("smooth_window must be a positive odd integer")
    if config.dominance_factor <= 1:
        raise ValueError("dominance_factor must be greater than one")

    out = df.sort_values(["Eddy", "Day"]).copy()
    out["topo_plan_ratio"] = pd.to_numeric(
        out["topo_plan_ratio"], errors="coerce"
    ).replace([np.inf, -np.inf], np.nan)
    out["topo_plan_ratio_smooth"] = out.groupby(
        "Eddy", sort=False
    )["topo_plan_ratio"].transform(
        lambda x: x.rolling(
            config.smooth_window,
            center=True,
            min_periods=min(config.min_periods, config.smooth_window),
        ).median()
    )
    threshold = float(np.log(config.dominance_factor))
    depth = pd.to_numeric(out["h"], errors="coerce")
    out["planetary_regime"] = (
        (out["topo_plan_ratio_smooth"] <= -threshold)
        & (depth >= config.min_planetary_depth_m)
    )
    out["topographic_regime"] = (
        (out["topo_plan_ratio_smooth"] >= threshold)
        & (depth <= config.max_topographic_depth_m)
    )
    out["mixed_regime"] = ~(
        out["planetary_regime"] | out["topographic_regime"]
    )

    raw_difference = signed_angle_difference(
        out["TiltDir"], out["PV_grad_theta"]
    )
    out["signed_dtheta_PV_grad"] = raw_difference
    out["dtheta_PV_grad"] = np.abs(raw_difference)
    # Correct signed-PV convention: AE along grad(PV), CE opposite grad(PV).
    out["expected_tilt_theta"] = np.where(
        out["Cyc"].eq("AE"),
        out["PV_grad_theta"] % 360.0,
        (out["PV_grad_theta"] + 180.0) % 360.0,
    )
    out["preference_error_deg"] = np.abs(signed_angle_difference(
        out["TiltDir"], out["expected_tilt_theta"]
    ))
    out["direction_valid"] = (
        out[["TiltDir", "PV_grad_theta"]].notna().all(axis=1)
        & (pd.to_numeric(out["TiltDis"], errors="coerce")
           >= config.min_tilt_distance_km)
    )
    out["expected_match"] = (
        out["direction_valid"]
        & (out["preference_error_deg"] <= config.angle_tolerance_deg)
    )
    return out


def _base_summary(part: pd.DataFrame, regime: pd.Series) -> dict:
    valid = regime & part["direction_valid"]
    match = valid & part["expected_match"]
    days = pd.to_numeric(part["Day"], errors="coerce")
    lifetime = float(days.max() - days.min()) if days.notna().any() else np.nan
    error = part.loc[valid, "preference_error_deg"]
    signed_target_error = signed_angle_difference(
        part.loc[valid, "TiltDir"], part.loc[valid, "expected_tilt_theta"]
    )
    return {
        "lifetime_days": lifetime,
        "n_observations": len(part),
        "directional_observations": int(part["direction_valid"].sum()),
        "directional_coverage": float(part["direction_valid"].mean()),
        "regime_observations": int(valid.sum()),
        "regime_fraction": float(valid.mean()),
        "longest_regime_run": longest_true_run(valid),
        "matching_fraction": float(match.sum() / valid.sum()) if valid.sum() else 0.0,
        "longest_matching_run": longest_true_run(match),
        "median_preference_error_deg": float(error.median()) if len(error) else np.nan,
        "preference_resultant": circular_resultant(signed_target_error),
        "median_tilt_km": float(part.loc[valid, "TiltDis"].median()) if valid.any() else np.nan,
        "median_depth_m": float(part.loc[valid, "h"].median()) if valid.any() else np.nan,
    }


def _identity(part: pd.DataFrame) -> dict:
    return {
        "Eddy": part["Eddy"].iloc[0],
        "Cyc": part["Cyc"].iloc[0],
        "Region": part["Region"].mode().iloc[0]
        if "Region" in part and not part["Region"].mode().empty else np.nan,
    }


def _weighted_percentile_score(
    ranking: pd.DataFrame,
    group_col: str,
    terms: dict[str, tuple[float, bool]],
) -> pd.Series:
    score = pd.Series(0.0, index=ranking.index)
    for column, (weight, higher_is_better) in terms.items():
        values = ranking[column] if higher_is_better else -ranking[column]
        score += weight * values.groupby(ranking[group_col]).rank(
            pct=True, method="average", na_option="bottom"
        )
    return score


def rank_planetary_reference_cases(
    df: pd.DataFrame,
    config: PaperCaseConfig = PaperCaseConfig(),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank stable deep-water planetary reference cases by polarity."""
    diagnosed = add_paper_case_diagnostics(df, config)
    rows = []
    for _, part in diagnosed.groupby("Eddy", sort=False):
        summary = _base_summary(part, part["planetary_regime"])
        summary.update(_identity(part))
        summary["eligible"] = (
            summary["lifetime_days"] >= config.min_lifetime_days
            and summary["directional_observations"] >= config.min_directional_observations
            and summary["regime_observations"] >= config.min_regime_observations
            and summary["longest_regime_run"] >= config.min_sustained_run
            and summary["longest_matching_run"] >= config.min_sustained_run
            and summary["regime_fraction"] >= config.min_reference_regime_fraction
        )
        rows.append(summary)
    ranking = pd.DataFrame(rows)
    terms = {
        "matching_fraction": (0.25, True),
        "median_preference_error_deg": (0.20, False),
        "longest_matching_run": (0.15, True),
        "regime_fraction": (0.10, True),
        "longest_regime_run": (0.10, True),
        "directional_coverage": (0.10, True),
        "lifetime_days": (0.05, True),
        "median_tilt_km": (0.05, True),
    }
    ranking["case_score"] = _weighted_percentile_score(ranking, "Cyc", terms)
    ranking.loc[~ranking["eligible"], "case_score"] = np.nan
    ranking = ranking.sort_values(
        ["Cyc", "eligible", "case_score"], ascending=[True, False, False]
    ).reset_index(drop=True)
    return diagnosed, ranking


def rank_topographic_response_cases(
    df: pd.DataFrame,
    config: PaperCaseConfig = PaperCaseConfig(),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank coherent CE opposition and directionally disrupted topographic AEs."""
    diagnosed = add_paper_case_diagnostics(df, config)
    rows = []
    for _, part in diagnosed.groupby("Eddy", sort=False):
        summary = _base_summary(part, part["topographic_regime"])
        summary.update(_identity(part))
        valid = part["topographic_regime"] & part["direction_valid"]
        relative_resultant = circular_resultant(
            part.loc[valid, "signed_dtheta_PV_grad"]
        )
        summary["relative_direction_resultant"] = relative_resultant
        summary["response_type"] = (
            "CE coherent opposition" if summary["Cyc"] == "CE"
            else "AE disrupted/variable"
        )
        if summary["Cyc"] == "CE":
            response_ok = summary["longest_matching_run"] >= config.min_sustained_run
            summary["response_quality"] = summary["matching_fraction"]
        else:
            response_ok = np.isfinite(relative_resultant)
            summary["response_quality"] = 1.0 - relative_resultant \
                if np.isfinite(relative_resultant) else 0.0
        summary["eligible"] = (
            summary["lifetime_days"] >= config.min_lifetime_days
            and summary["directional_observations"] >= config.min_directional_observations
            and summary["regime_observations"] >= config.min_regime_observations
            and summary["longest_regime_run"] >= config.min_sustained_run
            and summary["regime_fraction"] >= config.min_topographic_regime_fraction
            and response_ok
        )
        rows.append(summary)
    ranking = pd.DataFrame(rows)
    terms = {
        "response_quality": (0.30, True),
        "longest_regime_run": (0.20, True),
        "regime_fraction": (0.15, True),
        "directional_coverage": (0.10, True),
        "lifetime_days": (0.10, True),
        "median_tilt_km": (0.10, True),
        "regime_observations": (0.05, True),
    }
    ranking["case_score"] = _weighted_percentile_score(
        ranking, "response_type", terms
    )
    ranking.loc[~ranking["eligible"], "case_score"] = np.nan
    ranking = ranking.sort_values(
        ["response_type", "eligible", "case_score"],
        ascending=[True, False, False],
    ).reset_index(drop=True)
    return diagnosed, ranking


def _regime_direction_summary(part: pd.DataFrame, regime_col: str) -> dict:
    use = part[part[regime_col] & part["direction_valid"]]
    return {
        "n": len(use),
        "tilt_theta": circular_mean_deg(use["TiltDir"]),
        "pv_theta": circular_mean_deg(use["PV_grad_theta"]),
        "expected_theta": circular_mean_deg(use["expected_tilt_theta"]),
        "preference_error": float(use["preference_error_deg"].median())
        if len(use) else np.nan,
        "relative_resultant": circular_resultant(use["signed_dtheta_PV_grad"]),
        "tilt_km": float(use["TiltDis"].median()) if len(use) else np.nan,
    }


def rank_regime_transition_cases(
    df: pd.DataFrame,
    config: PaperCaseConfig = PaperCaseConfig(),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank within-eddy planetary-to-topographic transitions by polarity."""
    diagnosed = add_paper_case_diagnostics(df, config)
    rows = []
    for _, part in diagnosed.groupby("Eddy", sort=False):
        part = part.sort_values("Day")
        plan = _regime_direction_summary(part, "planetary_regime")
        topo = _regime_direction_summary(part, "topographic_regime")
        plan_run = longest_true_run(part["planetary_regime"] & part["direction_valid"])
        topo_run = longest_true_run(part["topographic_regime"] & part["direction_valid"])
        identity = _identity(part)
        days = pd.to_numeric(part["Day"], errors="coerce")
        lifetime = float(days.max() - days.min()) if days.notna().any() else np.nan
        pv_rotation = abs(signed_angle_difference(topo["pv_theta"], plan["pv_theta"]))
        tilt_rotation = abs(signed_angle_difference(topo["tilt_theta"], plan["tilt_theta"]))
        target_rotation = signed_angle_difference(
            topo["expected_theta"], plan["expected_theta"]
        )
        observed_rotation = signed_angle_difference(
            topo["tilt_theta"], plan["tilt_theta"]
        )
        rotation_error = abs(signed_angle_difference(
            observed_rotation, target_rotation
        ))
        if identity["Cyc"] == "CE":
            response_quality = 1.0 - np.nanmean([
                plan["preference_error"], topo["preference_error"], rotation_error
            ]) / 180.0
            response_label = "CE follows rotating opposite-PV target"
        else:
            topo_disruption = 1.0 - topo["relative_resultant"] \
                if np.isfinite(topo["relative_resultant"]) else np.nan
            response_quality = np.nanmean([
                1.0 - plan["preference_error"] / 180.0,
                topo_disruption,
            ])
            response_label = "AE loses planetary directional coherence"
        eligible = (
            lifetime >= config.min_lifetime_days
            and plan["n"] >= config.min_regime_observations
            and topo["n"] >= config.min_regime_observations
            and plan_run >= config.min_sustained_run
            and topo_run >= config.min_sustained_run
            and pv_rotation >= config.min_transition_rotation_deg
        )
        rows.append({
            **identity,
            "response_type": response_label,
            "lifetime_days": lifetime,
            "planetary_observations": plan["n"],
            "topographic_observations": topo["n"],
            "longest_planetary_run": plan_run,
            "longest_topographic_run": topo_run,
            "regime_balance": 2 * min(plan["n"], topo["n"]) / len(part),
            "planetary_preference_error_deg": plan["preference_error"],
            "topographic_preference_error_deg": topo["preference_error"],
            "topographic_relative_resultant": topo["relative_resultant"],
            "pv_rotation_deg": float(pv_rotation),
            "tilt_rotation_deg": float(tilt_rotation),
            "rotation_coupling_error_deg": float(rotation_error),
            "response_quality": float(response_quality),
            "planetary_tilt_km": plan["tilt_km"],
            "topographic_tilt_km": topo["tilt_km"],
            "eligible": eligible,
        })
    ranking = pd.DataFrame(rows)
    terms = {
        "response_quality": (0.25, True),
        "pv_rotation_deg": (0.15, True),
        "regime_balance": (0.15, True),
        "longest_planetary_run": (0.10, True),
        "longest_topographic_run": (0.10, True),
        "planetary_observations": (0.05, True),
        "topographic_observations": (0.05, True),
        "lifetime_days": (0.10, True),
        "rotation_coupling_error_deg": (0.05, False),
    }
    ranking["case_score"] = _weighted_percentile_score(
        ranking, "response_type", terms
    )
    ranking.loc[~ranking["eligible"], "case_score"] = np.nan
    ranking = ranking.sort_values(
        ["response_type", "eligible", "case_score"],
        ascending=[True, False, False],
    ).reset_index(drop=True)
    return diagnosed, ranking


def select_ranked_cases(
    ranking: pd.DataFrame,
    group_col: str,
    n_per_group: int = 3,
) -> dict[str, list]:
    """Select the top eligible eddies from every displayed ranking group."""
    return {
        group: part.loc[part["eligible"], "Eddy"].head(n_per_group).tolist()
        for group, part in ranking.groupby(group_col, sort=True)
    }


def plot_paper_case(
    track: pd.DataFrame,
    grid,
    *,
    config: PaperCaseConfig = PaperCaseConfig(),
    title: str | None = None,
):
    """Consistent six-panel case figure with local bathymetric context."""
    df = track.sort_values("Day").copy()
    required = {
        "Eddy", "Day", "Cyc", "lon", "lat", "h", "TiltDis", "TiltDir",
        "PV_grad_theta", "expected_tilt_theta", "preference_error_deg",
        "topo_plan_ratio", "topo_plan_ratio_smooth", "planetary_regime",
        "topographic_regime", "PV_grad_plan_mag", "PV_grad_topo_mag",
        "PV_grad_mag",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing plotting columns: {sorted(missing)}")
    cyc = df["Cyc"].iloc[0]
    day = df["Day"]
    fig = plt.figure(figsize=(14, 12), constrained_layout=True)
    gs = fig.add_gridspec(6, 2, width_ratios=[2.25, 1.35])
    axes = [fig.add_subplot(gs[i, 0]) for i in range(6)]
    ax_map = fig.add_subplot(gs[:, 1])

    axes[0].scatter(day, df["TiltDir"] % 360, s=22, color="tab:purple", label="Tilt")
    axes[0].scatter(day, df["PV_grad_theta"] % 360, s=18, marker="x",
                    color="black", label="Signed PV gradient")
    axes[0].scatter(day, df["expected_tilt_theta"] % 360, s=14, marker="|",
                    color="tab:green", label=f"{cyc} expected target")
    axes[0].set(ylim=(0, 360), yticks=[0, 90, 180, 270, 360],
                ylabel="Compass bearing (deg)")
    axes[0].legend(ncol=3, fontsize=8, frameon=False)

    tol = config.angle_tolerance_deg
    axes[1].plot(day, df["preference_error_deg"], "o-", ms=3, color="black")
    axes[1].axhspan(0, tol, color="tab:green", alpha=0.12,
                    label=f"Expected +/- {tol:g} deg")
    axes[1].set(ylim=(0, 180), yticks=[0, 45, 90, 135, 180],
                ylabel="Polarity-aware error (deg)")
    axes[1].legend(frameon=False, fontsize=8)

    ratio_limit = np.log(config.dominance_factor)
    axes[2].plot(day, df["topo_plan_ratio"], color="0.72", lw=1, label="Raw")
    axes[2].plot(day, df["topo_plan_ratio_smooth"], color="black", lw=2,
                 label="Smoothed")
    axes[2].axhline(-ratio_limit, color="tab:blue", ls="--", lw=1,
                    label=f"Planetary >= {config.dominance_factor:g}x")
    axes[2].axhline(ratio_limit, color="tab:orange", ls="--", lw=1,
                    label=f"Topographic >= {config.dominance_factor:g}x")
    axes[2].axhline(0, color="0.35", ls=":", lw=0.8)
    axes[2].set_ylabel("log topo/planetary")
    axes[2].legend(ncol=4, fontsize=8, frameon=False)

    axes[3].semilogy(day, df["PV_grad_plan_mag"], color="tab:blue", label="Planetary")
    axes[3].semilogy(day, df["PV_grad_topo_mag"], color="tab:orange", label="Topographic")
    axes[3].semilogy(day, df["PV_grad_mag"], color="black", lw=1.5, label="Total")
    axes[3].set_ylabel("PV-gradient magnitude")
    axes[3].legend(ncol=3, fontsize=8, frameon=False)

    axes[4].plot(day, df["TiltDis"], color="tab:purple", lw=1.5,
                 label="Tilt distance")
    axes[4].set_ylabel("Tilt distance (km)")
    ax_depth = axes[4].twinx()
    ax_depth.plot(day, df["h"], color="saddlebrown", alpha=0.65)
    ax_depth.set_ylabel("Depth (m)", color="saddlebrown")
    ax_depth.invert_yaxis()

    if "Ro" in df:
        axes[5].plot(day, df["Ro"], color="tab:red", label="|Ro|")
        axes[5].axhline(1, color="0.3", ls="--", lw=1)
        axes[5].set_ylabel("Rossby number")
    elif "w" in df:
        axes[5].plot(day, df["w"], color="tab:red", label="Relative vorticity")
        axes[5].set_ylabel("Relative vorticity (s$^{-1}$)")
    axes[5].set_xlabel("Day")
    axes[5].legend(frameon=False, fontsize=8)

    for ax in axes:
        ax.grid(alpha=0.2)
        ymin, ymax = ax.get_ylim()
        ax.fill_between(day, ymin, ymax, where=df["planetary_regime"],
                        color="tab:blue", alpha=0.035, step="mid")
        ax.fill_between(day, ymin, ymax, where=df["topographic_regime"],
                        color="tab:orange", alpha=0.035, step="mid")
        ax.set_ylim(ymin, ymax)

    ratio = df["topo_plan_ratio_smooth"]
    sc = ax_map.scatter(df["lon"], df["lat"], c=ratio, cmap="coolwarm",
                        vmin=-2, vmax=2, s=38, zorder=3)
    ax_map.plot(df["lon"], df["lat"], color="0.35", lw=1, zorder=2)
    ax_map.scatter(df["lon"].iloc[0], df["lat"].iloc[0], marker="o", s=90,
                   facecolor="none", edgecolor="black", label="Start", zorder=4)
    ax_map.scatter(df["lon"].iloc[-1], df["lat"].iloc[-1], marker="x", s=90,
                   color="black", label="End", zorder=4)
    pad_lon = max(0.25, 0.15 * np.ptp(df["lon"]))
    pad_lat = max(0.25, 0.15 * np.ptp(df["lat"]))
    ax_map.set_xlim(df["lon"].min() - pad_lon, df["lon"].max() + pad_lon)
    ax_map.set_ylim(df["lat"].min() - pad_lat, df["lat"].max() + pad_lat)
    bathy = np.where(grid.mask_rho, grid.h, np.nan)
    levels = np.unique(np.r_[0, 200, 500, 1000, 2000, 3000, 4000,
                             np.nanmax(grid.h)])
    ax_map.contourf(grid.lon_rho, grid.lat_rho, bathy, levels=levels,
                    cmap="Blues", alpha=0.72, extend="max")
    ax_map.contour(grid.lon_rho, grid.lat_rho, bathy,
                   levels=[200, 1000, 2000, 3000], colors="0.35", linewidths=0.7)
    ax_map.set(xlabel="Longitude", ylabel="Latitude", title="Track over bathymetry")
    ax_map.legend(frameon=False)
    cbar = fig.colorbar(sc, ax=ax_map, pad=0.02, shrink=0.75)
    cbar.set_label("Smoothed log topo/planetary")

    if title is None:
        title = f"{cyc} eddy {df['Eddy'].iloc[0]}"
    fig.suptitle(title, fontsize=15)
    return fig, axes, ax_map
