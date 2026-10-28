"""Two circular-truth results figures and one shared geometry supplement."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.ticker import MaxNLocator

from rc3_validation import METRICS, outer_count

LABELS = ("Centre error (km)", "Final Ω error (%)", "Aspect-ratio error (%)", "Gaussian Rc error (%)")


def _panel(ax, letter, label=""):
    ax.text(.01, 1.02, f"({letter}) {label}", transform=ax.transAxes, va="bottom")
    ax.grid(alpha=.18)


def _scale(ax, config, metric, values):
    values = np.asarray(values, float)
    finite = np.abs(values[np.isfinite(values)])
    settings = config.values["plotting"]
    threshold = settings["centre_linthresh_km"] if metric == "centre_error_km" else settings["relative_linthresh_pct"]
    use_log = len(finite) and finite.max() >= threshold * 10**settings["symlog_min_decades"]
    if use_log:
        ax.set_yscale("symlog", linthresh=threshold)
        unit = "km" if metric == "centre_error_km" else "%"
        return f" (symlog: linear ≤{threshold:g} {unit})"
    return ""


def _acceptance_caption(results):
    text = []
    for method, group in results.groupby("method", sort=False):
        text.append(f"{method} {int(group.success.sum())}/{len(group)}")
    caption = "Accepted/attempted: " + "; ".join(text) + "."
    if not results.success.all():
        keys = ["method", "density", "noise_fraction"] if results.attrs["experiment"] == "noise_sampling" else ["method", "background_direction", "background_fraction"]
        rates = results.groupby(keys).success.mean()
        caption += f" Lowest condition acceptance: {100*rates.min():g}% ({rates.idxmin()}); full counts in notebook/CSV."
    return caption


def plot_noise_sampling(results, config):
    settings = config.values["plotting"]
    fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey="col")
    plotted = {key: [] for key in METRICS}
    for row, method in enumerate(config.methods):
        for regime in config.regimes:
            group = results.loc[(results.method == method) & (results.density == regime)]
            color = settings["regime_colors"][regime]
            fractions = sorted(group.noise_fraction.unique())
            x = 100*np.asarray(fractions)
            for col, metric in enumerate(METRICS):
                if method == "SOLO" and metric == "alpha_error_pct":
                    continue
                quantiles = []
                for fraction in fractions:
                    values = group.loc[(group.noise_fraction == fraction) & group.success, metric].dropna()
                    quantiles.append(values.quantile([settings["percentiles"][0]/100, .5,
                                                      settings["percentiles"][1]/100]).to_numpy()
                                     if len(values) else np.full(3, np.nan))
                low, median, high = np.asarray(quantiles).T
                ax = axes[row, col]
                ax.plot(x, median, color=color, marker="o", markersize=3)
                ax.fill_between(x, low, high, color=color, alpha=.15, linewidth=0)
                plotted[metric].extend(np.r_[low, median, high])
                missing = ~np.isfinite(median)
                if missing.any():
                    ax.plot(x[missing], np.full(missing.sum(), .95), "x", color=color,
                            transform=ax.get_xaxis_transform())
        for col, ax in enumerate(axes[row]):
            _panel(ax, chr(97+row*4+col), method)
            ax.set_ylabel(LABELS[col])
            ax.set_xticks(100*np.asarray(config.values["noise"]["fractions"]))
            if row == 2:
                ax.set_xlabel(r"Component noise SD / $U_{\rm ref}$ (%)")
        if method == "SOLO":
            axes[row, 2].text(.5, .5, "Not estimated", ha="center", va="center", transform=axes[row, 2].transAxes)
            axes[row, 2].set_axis_off()
    scale_notes = []
    for col, metric in enumerate(METRICS):
        note = _scale(axes[0, col], config, metric, plotted[metric])
        if note:
            scale_notes.append(LABELS[col] + note)
        axes[0, col].autoscale(enable=True, axis="y")
        axes[0, col].set_ylim(bottom=0)
        for row in range(3):
            axes[row, col].set_ylabel(LABELS[col] + ("\n(symlog)" if note else ""))
    handles = [Line2D([], [], color=settings["regime_colors"][r], marker="o", label=r) for r in config.regimes]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)
    caption = ("Median and 10th–90th percentiles among accepted full fits; bands are Monte Carlo spread, not confidence intervals.\n"
               f"{config.values['noise']['replicates']} trials per nonzero condition; one clean trial. No error-based trimming.\n" +
               _acceptance_caption(results))
    if scale_notes:
        caption += "\n" + "; ".join(scale_notes)
    fig.text(.5, .01, caption, ha="center", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, .14 if scale_notes else .12, 1, .94))
    return fig


def plot_background_flow(results, config):
    settings = config.values["plotting"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2), sharex=True)
    for method in config.methods:
        group = results.loc[results.method == method]
        baseline = group.loc[group.background_fraction == 0]
        assert len(baseline) == 1, "There must be one shared zero-current baseline per method."
        for direction in config.values["background"]["directions"][method]:
            series = pd.concat((baseline, group.loc[group.background_direction == direction])).sort_values("background_fraction")
            x = 100*series.background_fraction.to_numpy()
            for col, metric in enumerate(METRICS):
                if method == "SOLO" and metric == "alpha_error_pct":
                    continue
                axes[col].plot(x, series[metric].where(series.success),
                    color=settings["method_colors"][method], linestyle=settings["direction_styles"][direction],
                    marker="o", markersize=3)
                failed = ~series.success.to_numpy()
                if failed.any():
                    axes[col].plot(x[failed], np.full(failed.sum(), .94), "x",
                        color=settings["method_colors"][method], transform=axes[col].get_xaxis_transform())
    notes = []
    for col, (ax, metric) in enumerate(zip(axes, METRICS)):
        _panel(ax, chr(97+col))
        values = results.loc[results.success, metric].dropna()
        note = _scale(ax, config, metric, values)
        if note:
            notes.append(LABELS[col] + note)
            # largest = results.loc[results.success & results[metric].notna()].sort_values(metric).iloc[-1]
            # ax.annotate(f"max {largest[metric]:.3g}",
            #             xy=(100*largest.background_fraction, largest[metric]), xytext=(-5, -15),
            #             textcoords="offset points", ha="right", va="top",
            #             color=settings["method_colors"][largest.method])
        ax.set_ylim(bottom=0)
        ax.set_ylabel(LABELS[col] + ("\n(symlog)" if note else ""))
        ax.set_xlabel(r"$|U_{\rm bg}|/U_{\rm ref}$ (%)")
        ax.set_xticks(100*np.asarray(config.values["background"]["fractions"]))
    handles = [Line2D([], [], color=settings["method_colors"][m], label=m) for m in config.methods]
    handles += [Line2D([], [], color=".3", linestyle=settings["direction_styles"][d],
                        label="Eastward (SOLO)" if d == "eastward" else "Northward")
                for d in ("eastward", "northward")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False)
    caption = ("HIGH density; total absolute errors against the common circular truth. Deterministic fits, no noise or uncertainty bands.\n"
               "SOLO: both directions; DOPPIO/LATTE: northward only. SOLO α not estimated. One baseline per method.\n" +
               _acceptance_caption(results))
    if notes:
        caption += "\n" + "; ".join(notes)
    fig.text(.5, .01, caption, ha="center", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, .23 if notes else .19, 1, .9))
    return fig


def _geometry_panel(ax, geometry, context, letter, label, color):
    truth, config = context.truth, context.config
    rc = truth["rc_m"]/1000
    xy = geometry["xy"]/1000
    ax.scatter(xy[:, 0], xy[:, 1], s=9, color=color, zorder=3)
    ax.plot(0, 0, "kx", markersize=5, zorder=4)
    for radius, style in ((rc/np.sqrt(2), "--"), (rc, "-")):
        ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor=".45", linestyle=style, linewidth=.8))
    _panel(ax, letter, label)
    ax.grid(False)
    ax.text(.03, .97, f"N = {len(xy)}; N_outer = {outer_count(geometry, truth)}",
            ha="left", va="top", transform=ax.transAxes)
    limit = config.values["plotting"]["geometry_limit_rc"] * rc
    ax.set_xlim(-limit, limit); ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))


def _geometry_legend(fig):
    handles = [Line2D([], [], color="k", marker="x", linestyle="", label="True centre"),
               Line2D([], [], color=".45", linestyle="--", label=r"$r=R_c/\sqrt{2}$"),
               Line2D([], [], color=".45", linestyle="-", label=r"$r=R_c$")]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False)


def plot_noise_geometry(context):
    c = context.config
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), sharex=True, sharey=True)
    for row, method in enumerate(c.methods):
        for col, regime in enumerate(c.regimes):
            ax = axes[row, col]
            _geometry_panel(ax, context.geometries[method][regime], context,
                            chr(97+row*3+col), f"{method} — {regime}", c.values["plotting"]["regime_colors"][regime])
            if col == 0: ax.set_ylabel("y (km)")
            if row == 2: ax.set_xlabel("x (km)")
    _geometry_legend(fig)
    fig.text(.5, .01, "S1. Actual fixed observation locations used at each density; N counts unique velocity vectors.\n"
             "N_outer counts unique locations with true r > Rc/√2. DOPPIO's intersection is counted once.\n"
             "HIGH panels (c), (f), (i) are reused for background flow: SOLO eastward/northward; DOPPIO/LATTE northward only.", ha="center", va="bottom")
    fig.tight_layout(rect=(0, .09, 1, .95))
    return fig



def save_figure(figure, stem, config):
    stem.parent.mkdir(parents=True, exist_ok=True)
    path = stem.with_suffix(".png")
    figure.savefig(path, dpi=config.values["plotting"]["dpi"], bbox_inches="tight")
    return {"png": path}
