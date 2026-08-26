"""Selection and plotting tools for SEACOFS eddy case studies."""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class TransitionConfig:
    """Criteria used to rank planetary/topographic transition cases."""

    smooth_window: int = 7
    min_periods: int = 5
    min_lifetime_days: float = 100.0
    min_regime_observations: int = 10
    min_sustained_run: int = 5
    min_tilt_observations: int = 20
    regime_margin: float = 0.0


def _longest_true_run(values) -> int:
    longest = current = 0
    for value in np.asarray(values, dtype=bool):
        current = current + 1 if value else 0
        longest = max(longest, current)
    return int(longest)


def _sign_changes(values) -> int:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return 0
    signs = np.sign(values)
    signs = signs[signs != 0]
    return int(np.count_nonzero(signs[1:] != signs[:-1]))


def add_transition_diagnostics(
    df: pd.DataFrame,
    config: TransitionConfig = TransitionConfig(),
) -> pd.DataFrame:
    """Add smoothed PV dominance and sustained-regime diagnostics by eddy.

    ``topo_plan_ratio`` is log(|topographic PV gradient| / |planetary PV
    gradient|), so positive values are topographic-dominant and negative
    values are planetary-dominant.
    """
    required = {"Eddy", "Day", "topo_plan_ratio"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")
    if config.smooth_window < 1 or config.smooth_window % 2 == 0:
        raise ValueError("smooth_window must be a positive odd integer")

    out = df.sort_values(["Eddy", "Day"]).copy()
    out["topo_plan_ratio"] = pd.to_numeric(
        out["topo_plan_ratio"], errors="coerce"
    ).replace([np.inf, -np.inf], np.nan)
    grouped = out.groupby("Eddy", sort=False)["topo_plan_ratio"]
    out["topo_plan_ratio_smooth"] = grouped.transform(
        lambda x: x.rolling(
            config.smooth_window,
            center=True,
            min_periods=min(config.min_periods, config.smooth_window),
        ).median()
    )
    margin = float(config.regime_margin)
    out["topographic_dominant"] = out["topo_plan_ratio_smooth"] > margin
    out["planetary_dominant"] = out["topo_plan_ratio_smooth"] < -margin
    out["transition_zone"] = ~(
        out["topographic_dominant"] | out["planetary_dominant"]
    )
    return out


def _summarise_eddy(part: pd.DataFrame, config: TransitionConfig) -> dict:
    part = part.sort_values("Day")
    smooth = part["topo_plan_ratio_smooth"]
    valid = smooth.notna()
    topo = part["topographic_dominant"] & valid
    plan = part["planetary_dominant"] & valid
    n_valid = int(valid.sum())
    n_topo = int(topo.sum())
    n_plan = int(plan.sum())
    n_tilt = int(part.get("TiltDis", pd.Series(index=part.index, dtype=float)).notna().sum())
    days = pd.to_numeric(part["Day"], errors="coerce")
    lifetime = float(days.max() - days.min()) if days.notna().any() else np.nan
    signs = np.where(topo, 1.0, np.where(plan, -1.0, np.nan))
    changes = _sign_changes(signs)

    topo_values = smooth[topo]
    plan_values = smooth[plan]
    separation = (
        float(topo_values.median() - plan_values.median())
        if len(topo_values) and len(plan_values)
        else np.nan
    )
    balance = 2 * min(n_topo, n_plan) / n_valid if n_valid else 0.0
    coverage = n_tilt / len(part) if len(part) else 0.0
    longest_topo = _longest_true_run(topo)
    longest_plan = _longest_true_run(plan)
    sustained = min(longest_topo, longest_plan)

    eligible = (
        lifetime >= config.min_lifetime_days
        and n_topo >= config.min_regime_observations
        and n_plan >= config.min_regime_observations
        and longest_topo >= config.min_sustained_run
        and longest_plan >= config.min_sustained_run
        and n_tilt >= config.min_tilt_observations
        and changes >= 1
    )

    row = {
        "Eddy": part["Eddy"].iloc[0],
        "Cyc": part["Cyc"].iloc[0] if "Cyc" in part else np.nan,
        "Region": part["Region"].mode().iloc[0]
        if "Region" in part and not part["Region"].mode().empty
        else np.nan,
        "n_observations": len(part),
        "lifetime_days": lifetime,
        "n_valid_ratio": n_valid,
        "n_planetary": n_plan,
        "n_topographic": n_topo,
        "planetary_fraction": n_plan / n_valid if n_valid else np.nan,
        "topographic_fraction": n_topo / n_valid if n_valid else np.nan,
        "regime_balance": balance,
        "longest_planetary_run": longest_plan,
        "longest_topographic_run": longest_topo,
        "sustained_run_min": sustained,
        "regime_separation": separation,
        "sign_changes": changes,
        "smoothness": 1.0 / changes if changes else 0.0,
        "tilt_observations": n_tilt,
        "tilt_coverage": coverage,
        "median_tilt_km": part.get(
            "TiltDis", pd.Series(index=part.index, dtype=float)
        ).median(),
        "eligible": eligible,
    }
    return row


def rank_topographic_transitions(
    df: pd.DataFrame,
    config: TransitionConfig = TransitionConfig(),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank sustained planetary/topographic transitions within each polarity.

    Returns the row-level diagnostics and an eddy-level ranking. The score is
    a weighted combination of within-polarity percentile ranks, preventing the
    more numerous polarity from controlling the candidate list.
    """
    diagnosed = add_transition_diagnostics(df, config)
    rows = [
        _summarise_eddy(part, config)
        for _, part in diagnosed.groupby("Eddy", sort=False)
    ]
    ranking = pd.DataFrame(rows)
    if ranking.empty:
        ranking["transition_score"] = pd.Series(dtype=float)
        return diagnosed, ranking

    score_terms = {
        "lifetime_days": 0.20,
        "regime_balance": 0.25,
        "sustained_run_min": 0.20,
        "regime_separation": 0.15,
        "smoothness": 0.10,
        "tilt_coverage": 0.10,
    }
    ranking["transition_score"] = 0.0
    for column, weight in score_terms.items():
        percentiles = ranking.groupby("Cyc")[column].rank(
            pct=True, method="average", na_option="bottom"
        )
        ranking["transition_score"] += weight * percentiles

    ranking.loc[~ranking["eligible"], "transition_score"] = np.nan
    ranking = ranking.sort_values(
        ["Cyc", "eligible", "transition_score", "lifetime_days"],
        ascending=[True, False, False, False],
    ).reset_index(drop=True)
    return diagnosed, ranking


def select_top_cases(ranking: pd.DataFrame, n_per_polarity: int = 3) -> dict:
    """Return equally sized ranked candidate lists for AE and CE."""
    selected = {}
    for cyc in ("AE", "CE"):
        selected[cyc] = (
            ranking.loc[(ranking["Cyc"] == cyc) & ranking["eligible"], "Eddy"]
            .head(n_per_polarity)
            .tolist()
        )
    return selected


def transition_days(track: pd.DataFrame) -> np.ndarray:
    """Return midpoint days where the smoothed dominance sign changes."""
    part = track.sort_values("Day")
    ratio = part["topo_plan_ratio_smooth"].to_numpy(dtype=float)
    day = part["Day"].to_numpy(dtype=float)
    valid = np.isfinite(ratio) & np.isfinite(day) & (ratio != 0)
    ratio, day = ratio[valid], day[valid]
    if len(ratio) < 2:
        return np.array([], dtype=float)
    crossings = np.flatnonzero(np.sign(ratio[1:]) != np.sign(ratio[:-1]))
    return (day[crossings] + day[crossings + 1]) / 2


def plot_topographic_transition(track: pd.DataFrame, grid, *, title=None):
    """Plot the track and tilt/PV diagnostics for one transition candidate."""
    df = track.sort_values("Day").copy()
    required = {
        "Day", "lon", "lat", "h", "dhdx", "dhdy", "TiltDis", "TiltDir",
        "topo_plan_ratio", "topo_plan_ratio_smooth", "PV_grad_plan_mag",
        "PV_grad_topo_mag", "PV_grad_mag", "PV_grad_plan_theta",
        "PV_grad_topo_theta", "PV_grad_theta", "w", "PV",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing plotting columns: {sorted(missing)}")

    df["slope"] = np.hypot(df["dhdx"], df["dhdy"])
    fig = plt.figure(figsize=(15, 13), constrained_layout=True)
    gs = fig.add_gridspec(7, 2, width_ratios=[2.3, 1.35])
    axes = [fig.add_subplot(gs[i, 0]) for i in range(7)]
    ax_map = fig.add_subplot(gs[:, 1])
    day = df["Day"]

    axes[0].plot(day, df["topo_plan_ratio"], color="0.72", lw=1, label="Raw")
    axes[0].plot(day, df["topo_plan_ratio_smooth"], color="black", lw=2, label="Smoothed")
    axes[0].axhline(0, color="black", ls="--", lw=0.9)
    axes[0].fill_between(day, 0, df["topo_plan_ratio_smooth"],
                         where=df["topo_plan_ratio_smooth"] >= 0,
                         color="tab:orange", alpha=0.25, label="Topographic")
    axes[0].fill_between(day, 0, df["topo_plan_ratio_smooth"],
                         where=df["topo_plan_ratio_smooth"] < 0,
                         color="tab:blue", alpha=0.18, label="Planetary")
    axes[0].set_ylabel("log topo/planetary")
    axes[0].legend(ncol=4, fontsize=8, frameon=False)

    axes[1].plot(day, df["TiltDis"], "o-", ms=3, color="tab:purple")
    axes[1].set_ylabel("Tilt distance (km)")

    axes[2].scatter(day, df["TiltDir"] % 360, c=df["topo_plan_ratio_smooth"],
                    cmap="coolwarm", vmin=-2, vmax=2, s=22, label='tilt theta')
    axes[2].plot(day, df["PV_grad_theta"], color="black", label='PV_grad_theta')
    axes[2].set(ylim=(0, 360), yticks=[0, 90, 180, 270, 360],
                ylabel="Bearing (°)")
    axes[2].legend(ncol=2, fontsize=8, frameon=False, loc="upper left")

    axes[3].semilogy(day, df["PV_grad_plan_mag"], label="Planetary", color="tab:blue")
    axes[3].semilogy(day, df["PV_grad_topo_mag"], label="Topographic", color="tab:orange")
    axes[3].semilogy(day, df["PV_grad_mag"], label="Total", color="black", lw=1.5)
    axes[3].set_ylabel("PV-gradient magnitude")
    axes[3].legend(ncol=3, fontsize=8, frameon=False)

    # axes[4].scatter(day, df["PV_grad_plan_theta"] % 360, s=14, label="Planetary")
    # axes[4].scatter(day, df["PV_grad_topo_theta"] % 360, s=14, label="Topographic")
    # axes[4].scatter(day, df["PV_grad_theta"] % 360, s=14, label="Total", color="black")
    axes[4].scatter(day, df["dtheta_PV_grad"] % 360, s=14, label="Total", color="black")
    axes[4].set(ylim=(0, 360), yticks=[0, 90, 180, 270, 360],
                ylabel="|tilt - PV-gradient| bearing (°)")
    axes[4].legend(ncol=3, fontsize=8, frameon=False)

    axes[5].plot(day, df["h"], color="saddlebrown", label="Depth")
    axes[5].set_ylabel("Depth (m)", color="saddlebrown")
    ax_slope = axes[5].twinx()
    ax_slope.plot(day, df["slope"], color="tab:green", label="Slope")
    ax_slope.set_ylabel("|∇h|", color="tab:green")
    axes[5].invert_yaxis()

    axes[6].plot(day, df["w"], label="Relative vorticity", color="tab:red")
    # axes[6].plot(day, df["PV"], label="Potential vorticity", color="black")
    axes[6].set_ylabel("Vorticity (s$^{-1}$)")
    axes[6].set_xlabel("Day")
    ax_ro = axes[6].twinx()
    ax_ro.plot(day, df["PV"], color="tab:green", alpha=0.65, label="PV")
    ax_ro.set_ylabel("Potential vorticity", color="tab:green")
    axes[6].legend(ncol=2, fontsize=8, frameon=False, loc="upper left")

    for ax in axes:
        ax.grid(alpha=0.2)
        for crossing in transition_days(df):
            ax.axvline(crossing, color="0.35", ls=":", lw=1)

    sc = ax_map.scatter(df["lon"], df["lat"], c=df["topo_plan_ratio_smooth"],
                        cmap="coolwarm", vmin=-2, vmax=2, s=35, zorder=3)
    ax_map.plot(df["lon"], df["lat"], color="0.35", lw=1, zorder=2)
    ax_map.scatter(df["lon"].iloc[0], df["lat"].iloc[0], marker="o", s=90,
                   facecolor="none", edgecolor="black", label="Start", zorder=4)
    ax_map.scatter(df["lon"].iloc[-1], df["lat"].iloc[-1], marker="x", s=90,
                   color="black", label="End", zorder=4)
    pad_lon = max(0.25, 0.15 * np.ptp(df["lon"]))
    pad_lat = max(0.25, 0.15 * np.ptp(df["lat"]))
    xlim = (df["lon"].min() - pad_lon, df["lon"].max() + pad_lon)
    ylim = (df["lat"].min() - pad_lat, df["lat"].max() + pad_lat)
    levels = np.unique(np.r_[0, 200, 500, 1000, 2000, 3000, 4000,
                             np.nanmax(grid.h)])
    ax_map.contourf(grid.lon_rho, grid.lat_rho,
                    np.where(grid.mask_rho, grid.h, np.nan),
                    levels=levels, cmap="Blues", alpha=0.7, extend="max")
    ax_map.contour(grid.lon_rho, grid.lat_rho,
                   np.where(grid.mask_rho, grid.h, np.nan),
                   levels=[200, 1000, 2000], colors="0.35", linewidths=0.7)
    ax_map.set(xlim=xlim, ylim=ylim, xlabel="Longitude", ylabel="Latitude",
               title="Track over bathymetry")
    ax_map.legend(frameon=False)
    cbar = fig.colorbar(sc, ax=ax_map, pad=0.02, shrink=0.75)
    cbar.set_label("Smoothed log topo/planetary")

    if title is None:
        cyc = df["Cyc"].iloc[0] if "Cyc" in df else ""
        title = f"{cyc} eddy {df['Eddy'].iloc[0]}"
    fig.suptitle(title, fontsize=15)
    return fig, axes, ax_map
