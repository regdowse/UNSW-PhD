"""Original tagged Figure 3 three-panel layout and rendering, adapted locally."""
import numpy as np
import matplotlib.pyplot as plt


def plot(EDDY, DAY, profiles, fit, settings):
    NUM=settings["window_days"]
    MAX_DEPTH=settings["max_depth_m"]
    DEPTH_INT=settings["depth_interval_m"]
    window_days = np.arange(DAY - NUM // 2, DAY + NUM // 2 + 1)
    window = profiles.assign(Depth=profiles.Depth.abs())
    window = window.loc[window.Depth.le(MAX_DEPTH)].sort_values(['Day', 'Depth'])
    trace = fit['trace']
    points = trace[['x', 'y', 'Depth']].to_numpy()
    weights = trace.weight.to_numpy() / trace.weight.max()
    top, bottom = fit['ends']
    tilt_distance, tilt_direction = fit['TiltDis'], fit['TiltDir']
    axis = (bottom[:2] - top[:2]) / tilt_distance
    reference = window.loc[window.Day.eq(DAY)]
    reference_top, reference_bottom = reference.Depth.min(), reference.Depth.max()
    assert trace.Depth.min() >= reference_top
    assert trace.Depth.max() + DEPTH_INT <= reference_bottom + 1e-8


    plt.rcParams.update({'font.size': 11, 'axes.spines.top': True, 'axes.spines.right': True})
    fig, axs = plt.subplots(1, 3, figsize=settings['size_inches'], sharey=True)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.23, top=0.82, wspace=0.13)
    colors = ['#5399cf', '#ff993e', '#54ad55', '#df2429', '#a989d4', '#a77b6b', '#ef7bc9']
    ref = window.loc[window.Day.eq(DAY)].sort_values('Depth')
    ref_xy = ref[['xc', 'yc']].to_numpy()
    for j, day in enumerate(window_days):
        color = colors[(3 + int(day - DAY)) % len(colors)]
        g = window.loc[window.Day.eq(day)]
        if g.empty:
            continue
        offset = int(day - DAY)
        label = 'Day k' if offset == 0 else f'Day k{offset:+d}'
        xy0 = g[['xc', 'yc']].iloc[0].to_numpy()
        axs[0].plot(
            (g[['xc', 'yc']].to_numpy() - xy0) @ axis,
            g.Depth,
            color=color,
            lw=2.8 if offset == 0 else 1.1,
            label=label
        )
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=5, fontsize=11, frameon=False)
    axs[0].set_xlabel('Projected centre position (km)')
    axs[0].set_ylabel('Depth (m)')

    origin = points[0, :2] @ axis
    s = points[:, :2] @ axis - origin
    fit_s = np.array([top[:2] @ axis, bottom[:2] @ axis]) - origin
    fit_z = np.array([top[2], bottom[2]])
    raw_s = (ref_xy - ref_xy[0]) @ axis
    for ax in axs[1:]:
        ax.plot(s, points[:, 2], color='blue', lw=2, label='Mean-increment profile')
        ax.plot(fit_s, fit_z, 'k--', lw=2, label='Weighted line fit')
        ax.set_xlabel('Relative displacement (km)')
    width = weights * max(tilt_distance * 0.30, 1)
    axs[1].fill_betweenx(points[:, 2], s - width, s + width, color='#6baed6', alpha=0.5,
                         linewidth=0, label='Width ∝ fit weight')
    axs[1].lines[0].set_zorder(4)
    span = max(tilt_distance, np.ptp(raw_s), 1)
    arrow_x = fit_s - 0.10 * span
    axs[1].annotate('', xy=(arrow_x[0], fit_z[0] + 110),
                    xytext=(arrow_x[1], fit_z[1] - 60),
                    arrowprops={'arrowstyle': '-|>', 'color': 'black', 'lw': 2, 'mutation_scale': 19})
    axs[1].text(arrow_x.mean() - 0.09 * span, fit_z.mean(), 'Tilt direction', ha='right')
    axs[2].plot(raw_s, ref.Depth, color='#df2429', lw=2.5, label='Reference day k')

    def bracket(ax, left, right, depth, label, color, ends):
        ax.plot([left, right], [depth, depth], color=color, lw=1.4, clip_on=False)
        for x, end in zip([left, right], ends):
            ax.plot([x, x], [depth - 15, end], color=color, lw=1.3, clip_on=False)
        ax.text((left + right) / 2, depth - 18, label, ha='center', va='bottom', color=color)

    bracket(axs[1], *fit_s, -70, f'Tilt distance: {tilt_distance:.1f} km', 'black', fit_z)
    bracket(axs[2], *fit_s, -70, f'', 'black', fit_z)
        
    imin, imax = np.argmin(raw_s), np.argmax(raw_s)
    bracket(axs[2], raw_s[imin], raw_s[imax], -155,
            f'Max projected extent: {np.ptp(raw_s):.1f} km', '#df2429',
            [ref.Depth.iloc[imin], ref.Depth.iloc[imax]])
    lo = min(s.min() - width.max(), raw_s.min(), arrow_x.min()) - 0.20 * span
    hi = max(s.max() + width.max(), raw_s.max()) + 0.13 * span
    for ax in axs[1:]:
        ax.set_xlim(lo, hi)
    for label, ax in zip(['a)', 'b)', 'c)'], axs):
        ax.set_ylim(MAX_DEPTH, -35)
        ax.text(-0.07, 1.02, label, transform=ax.transAxes, fontweight='bold', fontsize=13)
    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), fontsize=11, frameon=False)
    return fig
