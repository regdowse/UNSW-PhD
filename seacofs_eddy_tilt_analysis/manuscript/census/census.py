"""Independent eddy-day census. Run: python census.py --config census.json."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import census_plots as plots

KEYS = ['Eddy', 'Day']
COLOURS = {'AE': '#C44E52', 'CE': '#2878A5'}
REGIONS = ['S1', 'S2', 'U1', 'D1', 'U2', 'D2']


def check_keys(frame, name):
    if frame[KEYS].isna().any().any() or frame.duplicated(KEYS).any():
        raise ValueError(f'{name}: missing or duplicate Eddy-Day keys')


def load_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    for name, value in config['inputs'].items():
        p = Path(value).expanduser()
        config['inputs'][name] = str(p if p.is_absolute() else path.parent / p)
    out = Path(config['output_dir']).expanduser()
    config['output_dir'] = str(out if out.is_absolute() else path.parent / out)
    figures = Path(config['figure_dir']).expanduser()
    config['figure_dir'] = str(figures if figures.is_absolute() else path.parent / figures)
    if config['direction_min_tilt_km'] < 0:
        raise ValueError('Direction threshold must be nonnegative')
    width = config['plot']['direction_bin_width_deg']
    if width <= 0 or not np.isclose(360 / width, round(360 / width)):
        raise ValueError('Direction bin width must divide 360')
    return config


def load_grid(path):
    """Retain legacy kilometre-grid construction; no vertical-depth file needed."""
    with netCDF4.Dataset(path) as ds:
        g = SimpleNamespace(**{k: np.asarray(ds[k][:]).T for k in
                               ['lon_rho', 'lat_rho', 'mask_rho', 'h']})
    def distance(lat1, lon1, lat2, lon2):
        a, b, c, d = map(np.radians, [lat1, lon1, lat2, lon2])
        term = np.sin((c-a)/2)**2 + np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2
        return 2*6357*np.arctan2(np.sqrt(term), np.sqrt(1-term))
    i, j = np.array(g.h.shape)//2
    x = np.r_[0, np.cumsum(distance(g.lat_rho[:-1, j], g.lon_rho[:-1, j],
                                    g.lat_rho[1:, j], g.lon_rho[1:, j]))]
    y = np.r_[0, np.cumsum(distance(g.lat_rho[i, :-1], g.lon_rho[i, :-1],
                                    g.lat_rho[i, 1:], g.lon_rho[i, 1:]))]
    g.X_grid, g.Y_grid = np.meshgrid(x, y, indexing='ij')
    return g


def region_grid(g, settings):
    """Legacy shelf construction, with one shared configurable latitude split."""
    shelf = ((g.h < settings['shelf_max_depth_m']) &
             (g.X_grid < settings['shelf_max_x_km']) &
             (g.lon_rho < settings['shelf_max_longitude_deg']) & (g.mask_rho == 1))
    dx = abs(np.nanmedian(np.diff(g.X_grid, axis=0)))
    expanded = shelf.copy()
    for shift in range(1, int(round(settings['shelf_extension_km']/dx))+1):
        shifted = np.roll(shelf, shift, axis=0)
        shifted[:shift] = False
        expanded |= shifted
    north = g.lat_rho >= settings['latitude_split_deg']
    west = g.lon_rho < settings['longitude_split_deg']
    ocean = g.mask_rho == 1
    bins = np.zeros(g.h.shape, dtype=int)
    for number, mask in enumerate([expanded & north, expanded & ~north,
                                   ~expanded & west & north, ~expanded & west & ~north,
                                   ~west & north, ~west & ~north], 1):
        bins[mask & ocean] = number
    return bins


def prepare(surface, tilt, grid, config):
    check_keys(surface, 'surface')
    check_keys(tilt, 'tilt')
    if not tilt.set_index(KEYS).index.isin(surface.set_index(KEYS).index).all():
        raise ValueError('Tilt contains Eddy-Day keys absent from surface')
    if surface.groupby('Eddy').Cyc.nunique().gt(1).any():
        raise ValueError('Polarity changes within an Eddy ID')
    if not surface.Cyc.isin(COLOURS).all():
        raise ValueError('Unrecognised polarity')
    d = surface.drop(columns=['TiltDis', 'TiltDir'], errors='ignore').merge(
        tilt[KEYS+['TiltDis', 'TiltDir']], on=KEYS, how='left', validate='one_to_one',
        indicator=True).sort_values(KEYS).reset_index(drop=True)
    finite_xy = np.isfinite(d[['xc', 'yc']]).all(axis=1)
    tree = cKDTree(np.c_[grid.X_grid.ravel(), grid.Y_grid.ravel()])
    _, index = tree.query(d.loc[finite_xy, ['xc', 'yc']].to_numpy())
    # Preserve wet-cell matches. For land matches, search from the matched grid
    # point (not the fitted eddy centre) to the closest wet grid point.
    land = grid.mask_rho.ravel()[index] == 0
    if land.any():
        wet = np.flatnonzero(grid.mask_rho.ravel() == 1)
        if not len(wet):
            raise ValueError('Cannot assign regions: grid has no wet cells')
        xy = np.c_[grid.X_grid.ravel(), grid.Y_grid.ravel()]
        _, nearest_wet = cKDTree(xy[wet]).query(xy[index[land]])
        index[land] = wet[nearest_wet]
    codes = np.zeros(len(d), dtype=int)
    codes[finite_xy] = region_grid(grid, config['regions']).ravel()[index]
    d['Region'] = pd.Series(codes).map(dict(enumerate(REGIONS, 1)))
    d['Sector'] = d.Region.map({'S1': 'Upstream', 'U1': 'Upstream', 'U2': 'Upstream',
                                'S2': 'Downstream', 'D1': 'Downstream', 'D2': 'Downstream'})
    d['valid_tilt'] = np.isfinite(d.TiltDis) & d.TiltDis.ge(0)
    d['valid_radius'] = np.isfinite(d.Rc) & d.Rc.gt(0)
    d['tilt_over_Rc'] = (d.TiltDis / d.Rc).where(d.valid_tilt & d.valid_radius)
    d['valid_direction'] = (d.valid_tilt & np.isfinite(d.TiltDir) &
                             d.TiltDis.ge(config['direction_min_tilt_km']))
    d['TiltDir'] = d.TiltDir % 360
    return d


def circular_mean(values):
    theta = np.deg2rad(np.asarray(values, float))
    theta = theta[np.isfinite(theta)]
    if not len(theta):
        return np.nan, np.nan
    z = np.exp(1j*theta).mean()
    return (float(np.angle(z, deg=True) % 360) if abs(z) > 1e-12 else np.nan,
            float(abs(z)))


def per_eddy(data):
    return data.groupby(['Cyc', 'Eddy'], as_index=False).agg(
        TiltDis=('TiltDis', 'mean'), tilt_over_Rc=('tilt_over_Rc', 'mean'),
        n_days=('Day', 'size'))


def summaries(data):
    rows = []
    for scope, column in [('All', None), ('sector', 'Sector'), ('region', 'Region')]:
        for cyc in COLOURS:
            d = data[data.Cyc.eq(cyc)]
            groups = [(scope, d)] if column is None else d.groupby(column)
            for name, group in groups:
                for unit, sample in [('eddy_day', group), ('per_eddy_mean', per_eddy(group))]:
                    for metric in ['TiltDis', 'tilt_over_Rc']:
                        v = sample[metric].dropna()
                        rows.append(dict(scope=scope, group=name, Cyc=cyc, unit=unit,
                                         metric=metric, n=len(v), eddies=sample.loc[
                                             sample[metric].notna(), 'Eddy'].nunique(),
                                         mean=v.mean(), median=v.median(), q25=v.quantile(.25),
                                         q75=v.quantile(.75), p95=v.quantile(.95), maximum=v.max()))
    return pd.DataFrame(rows)


def direction_summaries(data):
    rows = []
    eddies = []
    for (cyc, eddy), g in data.groupby(['Cyc', 'Eddy']):
        bearing, r = circular_mean(g.TiltDir)
        eddies.append(dict(Cyc=cyc, Eddy=eddy, TiltDir=bearing, R=r, n_days=len(g)))
    per = pd.DataFrame(eddies, columns=['Cyc', 'Eddy', 'TiltDir', 'R', 'n_days'])
    for scope, column in [('All', None), ('sector', 'Sector'), ('region', 'Region')]:
        for cyc in COLOURS:
            d = data[data.Cyc.eq(cyc)]
            groups = [(scope, d)] if column is None else d.groupby(column)
            for name, g in groups:
                bearing, r = circular_mean(g.TiltDir)
                rows.append(dict(scope=scope, group=name, Cyc=cyc, eddy_days=len(g),
                                 eddies=g.Eddy.nunique(), mean_bearing_deg=bearing, R=r,
                                 northward_pct=100*((g.TiltDir < 90) | (g.TiltDir >= 270)).mean(),
                                 eastward_pct=100*g.TiltDir.lt(180).mean()))
    return pd.DataFrame(rows), per


def rose_counts(data, settings):
    bins = [np.inf if v is None else v for v in settings['magnitude_bins_km']]
    result = plots.windrose_counts(data.TiltDir, data.TiltDis, mag_bins=bins)
    if result is None:
        return np.zeros((len(bins)-1, 16)), np.deg2rad(np.arange(16)*22.5), np.deg2rad(22.5)
    return result


def plot_magnitude(data, settings):
    groups = [(s, data[data.Sector.eq(s)]) for s in ['Upstream', 'Downstream']]
    eddy_groups = [(s, per_eddy(g)) for s, g in groups]
    q = settings['magnitude_display_quantile']
    fig, axes = plt.subplots(1, 4, figsize=settings['figure4_size_inches'], constrained_layout=True)
    for ax, metric, per, label in zip(axes,
            ['TiltDis', 'TiltDis', 'tilt_over_Rc', 'tilt_over_Rc'],
            [False, True, False, True],
            ['Tilt distance (km)', 'Mean tilt per eddy (km)',
             r'Tilt distance / $R_c$', r'Mean tilt / $R_c$ per eddy']):
        selected = eddy_groups if per else groups
        values = pd.concat([g[metric] for _, g in selected]).dropna()
        if per:
            bins = np.linspace(max(0, values.min()), values.quantile(q),
                               settings['per_eddy_density_edges'])
        else:
            bins = np.linspace(0, data[metric].quantile(q), settings['density_bins']+1)
        plots.vertical_split_violin_box(ax, selected, metric, '', label, bins=bins,
                    smoothing_sigma=settings['density_smoothing_sigma_bins'])
    for i, ax in enumerate(axes):
        plots.panel_label(ax, chr(97+i))
        ax.set_ylim(0, None)
    return fig


def plot_direction(data, per, settings):
    mag_bins = [np.inf if v is None else v for v in settings['magnitude_bins_km']]
    fig = plt.figure(figsize=settings['figure5_size_inches'], constrained_layout=True)
    gs = fig.add_gridspec(1, 4)
    axes = [fig.add_subplot(gs[0, i], projection='polar' if i < 2 else None) for i in range(4)]
    for ax, cyc in zip(axes[:2], COLOURS):
        sample = data[data.Cyc.eq(cyc)]
        colours = plt.get_cmap('Reds' if cyc == 'AE' else 'Blues')(
            np.linspace(.15, 1, len(mag_bins)-1))
        plots.plot_windrose(ax, sample, title='', mag_bins=mag_bins, colors=colours,
                           step=settings['figure5_radial_tick_step'],
                           rlim=settings['figure5_radial_limit'])
        bearing, _ = circular_mean(sample.TiltDir)
        if np.isfinite(bearing):
            th = np.deg2rad(bearing)
            ax.annotate('', xy=(th, .9*ax.get_ylim()[1]), xytext=(th, 0),
                        arrowprops=dict(arrowstyle='->', lw=2.5, color='magenta'))
        ax.legend(title=f'{cyc}\ntilt dist. (km)', loc='upper left',
                  bbox_to_anchor=(1.05, 1.05), frameon=False)
    edges = np.arange(0, 361, settings['direction_bin_width_deg'])
    for ax, sample, label in [(axes[2], data, 'Tilt direction (°)'),
                               (axes[3], per, 'Mean tilt dir. per eddy (°)')]:
        plots.mirrored_hist(ax, sample.loc[sample.Cyc.eq('AE'), 'TiltDir'],
                           sample.loc[sample.Cyc.eq('CE'), 'TiltDir'], edges, label,
                           normalize=True, colors=tuple(COLOURS.values()))
        ax.yaxis.set_major_formatter(plots.FuncFormatter(lambda y, _: f'{abs(y):g}'))
    for i, ax in enumerate(axes):
        plots.panel_label(ax, chr(97+i))
    return fig


def plot_regional(data, grid, config):
    settings = config['plot']
    regions = config['regions']
    mag_bins = [np.inf if v is None else v for v in settings['magnitude_bins_km']]
    fig, _ = plots.rose_plot(data, grid, mag_bins=mag_bins,
            frac=settings['figure6_rose_size_factor'],
            direction_offset=settings['figure6_display_rotation_deg'],
            lat_split=regions['latitude_split_deg'], lon_split=regions['longitude_split_deg'],
            shelf_offset=regions['shelf_extension_km'], show=False,
            legend_title='tilt dist. (km)', region_bins=region_grid(grid, regions),
            figsize=settings['figure6_size_inches'])
    return fig


def file_identity(path):
    path = Path(path)
    with path.open('rb') as handle:
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=digest.hexdigest())


def run(config_path=None):
    """Regenerate this topic's figures, CSV tables and concise writing handoff."""
    path = Path(config_path) if config_path else Path(__file__).with_suffix('.json')
    config = load_config(path)
    out = Path(config['output_dir'])
    out.mkdir(parents=True, exist_ok=True)
    figure_dir = Path(config['figure_dir'])
    figure_dir.mkdir(parents=True, exist_ok=True)
    grid = load_grid(config['inputs']['grid'])
    data = prepare(pd.read_parquet(config['inputs']['surface']),
                   pd.read_parquet(config['inputs']['tilt']), grid, config)
    magnitude = data[data.valid_tilt]
    direction = data[data.valid_direction]
    stats = summaries(magnitude)
    directions, eddies = direction_summaries(direction)
    counts = data.groupby('Cyc', as_index=False).agg(
        eddies=('Eddy', 'nunique'), eddy_days=('Day', 'size'), valid_tilt_days=('valid_tilt', 'sum'),
        directional_days=('valid_direction', 'sum'), valid_radius_days=('valid_radius', 'sum'))
    counts['tilt_coverage_pct'] = 100*counts.valid_tilt_days/counts.eddy_days
    audit = data[KEYS+['Cyc', '_merge', 'Region', 'Sector', 'valid_tilt', 'valid_radius', 'valid_direction']]
    tables = dict(counts=counts, magnitude=stats, direction=directions, per_eddy_direction=eddies,
                  sample_audit=audit)
    regional_counts = data.groupby(['Cyc', 'Region'], dropna=False).agg(
        eddy_days=('Day', 'size'), eddies=('Eddy', 'nunique'),
        magnitude_days=('valid_tilt', 'sum'), directional_days=('valid_direction', 'sum'))
    tables['regional_counts'] = regional_counts.reset_index()
    tails = []
    for cyc, sample in magnitude.groupby('Cyc'):
        ratios = sample.tilt_over_Rc.dropna()
        row = dict(Cyc=cyc, eddy_days=len(sample),
                   over_50_km_pct=100*sample.TiltDis.gt(50).mean(),
                   over_100_km_pct=100*sample.TiltDis.gt(100).mean(),
                   radius_days=len(ratios), over_one_Rc_pct=100*ratios.gt(1).mean())
        tails.append(row)
    tables['tails'] = pd.DataFrame(tails)
    for name, table in tables.items():
        table.to_csv(out/(name+'.csv'), index=False)
    rose_rows = []
    for scope, region in [('All', None)]+[(r, r) for r in REGIONS]:
        for cyc in COLOURS:
            sample = direction[direction.Cyc.eq(cyc)]
            if region:
                sample = sample[sample.Region.eq(region)]
            offset = config['plot']['figure6_display_rotation_deg'] if region else 0
            display_sample = sample.assign(TiltDir=(sample.TiltDir+offset) % 360)
            c, angles, width = rose_counts(display_sample, config['plot'])
            for i, row in enumerate(c):
                for j, value in enumerate(row):
                    rose_rows.append(dict(group=scope, Cyc=cyc, magnitude_bin=i,
                                          bearing_centre_deg=(np.rad2deg(angles[j])-offset) % 360,
                                          display_centre_deg=np.rad2deg(angles[j]),
                                          bin_width_deg=np.rad2deg(width), count=int(value)))
    pd.DataFrame(rose_rows).to_csv(out/'rose_counts.csv', index=False)
    with plt.rc_context(plots.STYLE):
        figures = {'tilt_census': plot_magnitude(magnitude, config['plot']),
                   'tilt_dir_census': plot_direction(direction, eddies, config['plot'])}
    # The regional rose notebook used the default Matplotlib style in its own kernel.
    with plt.rc_context(matplotlib.rcParamsDefault):
        figures['sectioned_rose'] = plot_regional(direction, grid, config)
    for name, fig in figures.items():
        fig.savefig(figure_dir/f'{name}.png', dpi=config['plot']['dpi'], bbox_inches='tight')
    lines = ['# Census results', '',
             f'Dataset: {len(data):,} eddy-days from {data.Eddy.nunique():,} eddies; '
             f'{len(magnitude):,} valid tilt estimates.', '',
             f'Directions require finite bearings and tilt ≥{config["direction_min_tilt_km"]:g} km '
             f'({len(direction):,} days). Magnitudes retain all finite nonnegative tilts.', '',
             'Primary summaries give each eddy-day equal weight. Per-eddy mean panels are '
             'comparisons; regional per-eddy means use only days within that sector.', '']
    for cyc in COLOURS:
        s = stats[(stats.scope == 'All') & (stats.Cyc == cyc) & (stats.unit == 'eddy_day')]
        a = s[s.metric == 'TiltDis'].iloc[0]
        b = s[s.metric == 'tilt_over_Rc'].iloc[0]
        d = directions[(directions.scope == 'All') & (directions.Cyc == cyc)].iloc[0]
        lines.append(f'- {cyc}: median tilt {a["median"]:.2f} km '
                     f'(IQR {a.q25:.2f}–{a.q75:.2f}; n={a.n:,}); median tilt/Rc '
                     f'{b["median"]:.3f} (IQR {b.q25:.3f}–{b.q75:.3f}; n={b.n:,}). '
                     f'Directional sample: {d.eddy_days:,} days / {d.eddies:,} eddies; '
                     f'circular mean {d.mean_bearing_deg:.1f}°, R={d.R:.3f}.')
        p = stats[(stats.scope == 'All') & (stats.Cyc == cyc) &
                  (stats.unit == 'per_eddy_mean') & (stats.metric == 'TiltDis')].iloc[0]
        lines.append(f'  Per-eddy comparison: median of mean tilts {p["median"]:.2f} km '
                     f'(IQR {p.q25:.2f}–{p.q75:.2f}; {p.eddies:,} eddies).')
    lines += ['', 'Figures: 4 `tilt_census`; 5 `tilt_dir_census`; 6 `sectioned_rose`.', '',
              'Figure 4: original split-violin styling and outlined boxes (median, IQR, 5–95%); '
              'density display bins truncate tails only; statistics include the full valid sample. ',
              'Figure 5: per-eddy directions average unit daily directions, without magnitude '
              'weighting; cancelling directions have no defined circular mean.',
              'Figure 6: original model-grid maps and inset roses; compass bearings are '
              'rotated only for display in grid axes. Offshore roses share a radial count '
              'scale; shelf roses have individual scales. Counts reflect sampling and preference.', '',
              'Rose sectors retain the legacy adaptive binning; the configured 15-degree '
              'bins apply to the mirrored direction histograms. rose_counts.csv records '
              'both geographic and displayed sector centres and their widths.', '',
              'Changes from legacy: primary magnitude distributions now weight eddy-days equally; '
              'direction figures exclude tilt below 5 km; all region assignments use the same '
              f'{config["regions"]["latitude_split_deg"]:g}° latitude boundary '
              '(legacy census used −33°, regional roses −33.5°).', '',
              'These are descriptive statistics, not independent-observation significance tests. '
              'IQR describes spread, not uncertainty. No bootstrap intervals are needed for '
              'these displayed census summaries. Coverage uses all source eddy-days; no '
              'assumed first/last-day eligibility correction is imposed.', '',
              f'Unassigned regions: {data.Region.isna().sum():,} days. '
              'Nearest land-cell matches use the wet cell closest to that land cell. '
              f'Valid tilt days without positive finite Rc: '
              f'{(data.valid_tilt & ~data.valid_radius).sum():,}. '
              'See sample_audit.csv for row-level exclusions; counts.csv for polarity counts.', '',
              'Tables: regional_counts.csv (including unassigned days), tails.csv '
              '(50/100 km and one-radius exceedances), magnitude.csv '
              '(all/sector/region, daily/per-eddy summaries), direction.csv '
              '(bearings and hemisphere fractions), per_eddy_direction.csv, rose_counts.csv. '
              'Distances are km; tilt/Rc and R are dimensionless; bearings are clockwise from '
              'north, deep-to-shallow. See manifest.json for exact settings and source identity.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    try:
        revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).parent,
                                             text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), config=config,
                    config_file=file_identity(path), code=file_identity(__file__), plotting_code=file_identity(plots.__file__), git_revision=revision,
                    inputs={k: file_identity(v) for k, v in config['inputs'].items()},
                    environment=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__,
                                     matplotlib=matplotlib.__version__, platform=platform.platform()))
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return dict(output_dir=out, figures=figures, tables=tables)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    result = run(args.config)
    print(result['output_dir'])
