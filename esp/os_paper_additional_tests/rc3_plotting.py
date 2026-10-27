"""One Figure per experiment; all summaries retain failed-fit denominators."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from rc3_validation import METRICS


LABELS = ("Centre error (km)", "Radius error (%)", "Aspect-ratio error (%)")


def _panel(ax, letter, label):
    ax.text(.02, 1.02, f"({letter}) {label}", transform=ax.transAxes, va="bottom")
    ax.grid(alpha=.2)
    ax.xaxis.set_major_locator(MaxNLocator(5))


def plot_noise_sampling(results, config):
    """Medians/central 80% spread of accepted fits and acceptance of ALL fits."""
    fig, axes = plt.subplots(len(config.methods), 4, figsize=(15, 3 * len(config.methods) + .8),
                             sharex=True, sharey="col", squeeze=False)
    colors = plt.get_cmap("viridis")(np.linspace(.12, .85, len(config.sample_counts)))
    minimum_acceptance = 100 * results.groupby(
        ["method", "n_samples", "noise_fraction"]).success.mean().min()
    for row, method in enumerate(config.methods):
        for color, n in zip(colors, config.sample_counts):
            group = results.loc[(results.method == method) & (results.n_samples == n)]
            fractions = sorted(group.noise_fraction.unique())
            x = 100 * np.asarray(fractions)
            for col, metric in enumerate(METRICS):
                if method == "SOLO" and metric == "alpha_error_pct":
                    continue
                medians, low, high = [], [], []
                for fraction in fractions:
                    values = group.loc[(group.noise_fraction == fraction) & group.success, metric].dropna()
                    if len(values):
                        q = values.quantile([.1, .5, .9]).values
                    else:
                        q = [np.nan] * 3
                    low.append(q[0]); medians.append(q[1]); high.append(q[2])
                ax = axes[row, col]
                ax.plot(x, medians, color=color, marker="o", markersize=3)
                ax.fill_between(x, low, high, color=color, alpha=.13, linewidth=0)
            acceptance = group.groupby("noise_fraction").success.mean().reindex(fractions)
            axes[row, 3].plot(x, 100 * acceptance.values, color=color, marker="o", markersize=3)
        for col, ax in enumerate(axes[row]):
            shape = results.loc[results.method == method, "truth_shape"].iloc[0]
            _panel(ax, chr(97 + 4 * row + col), f"{method} ({shape})")
            ax.set_ylabel((*LABELS, "Accepted fits (%)")[col])
            ax.set_ylim(bottom=0)
            if row == len(config.methods) - 1:
                ax.set_xlabel(r"Component noise SD / $U_{\rm ref}$ (%)")
            ax.set_xticks(100 * np.asarray(config.noise_fractions))
        # Shared metric scales retain the large DOPPIO tail without flattening
        # LATTE's smaller radius errors. Include the exact zero-noise value.
        axes[row, 1].set_yscale("symlog", linthresh=1.)
        axes[row, 3].set_ylim(max(0, minimum_acceptance - 1), 100.5)
        axes[row, 3].yaxis.set_major_locator(MaxNLocator(4))
        if method == "SOLO":
            axes[row, 2].text(.5, .5, "Shape fixed circular\nAspect ratio not estimated",
                              ha="center", va="center", transform=axes[row, 2].transAxes)
            axes[row, 2].set_axis_off()
    # Set shared limits only after every method's lines AND percentile bands
    # exist; SOLO's empty shape panel must not fix the shared shape range to 0..1.
    for col in range(3):
        axes[0, col].autoscale(enable=True, axis="y")
        axes[0, col].set_ylim(bottom=0)
    handles = [Line2D([], [], color=c, marker="o", label=f"N = {n}")
               for c, n in zip(colors, config.sample_counts)]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False)
    accepted = int(results.success.sum())
    fig.text(.5, .012,
             f"N: per transect for DOPPIO (2N−1 unique); total for SOLO/LATTE. "
             f"Accepted {accepted}/{len(results)} fits.\n"
             f"Median and 10th–90th percentiles among accepted fits; {config.replicates} trials per noisy condition, "
             "one clean trial.\n"
             "Fixed locations; independent u/v noise; U_ref from each truth. Radius axes: linear through 1%, logarithmic above.",
             ha="center", va="bottom")
    fig.tight_layout(rect=(0, .12, 1, .94))
    return fig


def plot_background_flow(results, config):
    """Absolute truth errors above; changes from zero background below."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    colors = {method: plt.get_cmap("tab10")(i) for i, method in enumerate(config.methods)}
    styles = {"eastward": "-", "northward": "--"}
    metrics = (METRICS, ("incremental_centre_km", "incremental_rc_pct", "incremental_alpha_pct"))
    labels = (LABELS, ("Centre change from baseline (km)", "Radius change from baseline (%)",
                       "Aspect-ratio change from baseline (%)"))
    for (method, direction), group in results.groupby(["method", "background_direction"], sort=False):
        group = group.sort_values("background_fraction")
        x = group.background_fraction.to_numpy() * 100
        for row in range(2):
            for col, metric in enumerate(metrics[row]):
                if method == "SOLO" and col == 2:
                    continue
                ax = axes[row, col]
                values = group[metric].where(group.success).to_numpy()
                ax.plot(x, values, color=colors[method], linestyle=styles[direction],
                        marker="o", markersize=3)
                failed_x = x[~group.success.to_numpy()]
                if len(failed_x):
                    ax.plot(failed_x, np.full(len(failed_x), .91), "x", color=colors[method],
                            transform=ax.get_xaxis_transform(), clip_on=False)
    for row in range(2):
        for col, ax in enumerate(axes[row]):
            _panel(ax, chr(97 + 3 * row + col), "Total error" if row == 0 else "Incremental change")
            ax.set_ylabel(labels[row][col])
            ax.axhline(0, color=".5", linewidth=.6)
            if row == 0 or col == 0:
                ax.set_ylim(bottom=0)
            if row == 1:
                ax.set_xlabel(r"Uniform background speed / $U_{\rm ref}$ (%)")
            ax.set_xticks(100 * np.asarray(config.background_fractions))
    # Use matching scales for corresponding method/direction comparisons.
    handles = [Line2D([], [], color=c, label=f"{m} ({results.loc[results.method == m, 'truth_shape'].iloc[0]})")
               for m, c in colors.items()]
    handles += [Line2D([], [], color=".3", linestyle=s, label=d.capitalize())
                for d, s in styles.items()]
    if not results.success.all():
        handles.append(Line2D([], [], color=".3", marker="x", linestyle="", label="Failed fit"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False)
    if config.all_circular:
        # Preserve very large converged radii rather than clip or discard them.
        # This applies only to the added comparison, preserving the original plot.
        for ax in axes[:, 1]:
            ax.set_yscale("symlog", linthresh=1.)
            ax.autoscale(enable=True, axis="y")
        axes[0, 1].set_ylim(bottom=0)
    scale_note = ("\nRadius axes: linear through ±1%, logarithmic beyond; no extreme estimates removed."
                  if config.all_circular else "")
    fig.text(.5, .012,
             f"N = {config.background_n}: per transect for DOPPIO; total for SOLO/LATTE. "
             f"Accepted {int(results.success.sum())}/{len(results)} fits.\n"
             "U_ref from each truth; SOLO shape not estimated. "
             "Incremental radius/shape changes are signed percentages of their zero-background estimates." + scale_note,
             ha="center", va="bottom")
    fig.tight_layout(rect=(0, .12 if config.all_circular else .09, 1, .94))
    return fig
