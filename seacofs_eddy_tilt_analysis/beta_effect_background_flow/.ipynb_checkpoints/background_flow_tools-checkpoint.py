"""Efficient background-flow diagnostics for the EAC beta-effect analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
DATE_ORIGIN = pd.Timestamp("1990-01-01")
FIELD_NAMES = ("u_surface", "v_surface", "u_500", "v_500")
ANNULUS_COLUMNS = (
    "Eddy", "Day", "month", "ic", "jc",
    "ann_surface_east_ms", "ann_surface_north_ms",
    "ann_500_east_ms", "ann_500_north_ms",
)


@dataclass(frozen=True)
class BackgroundConfig:
    model_root: Path = Path("/srv/scratch/z3533156/26year_BRAN2020")
    cache_root: Path = Path(
        "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/background_flow_cache_v2"
    )
    depth_max_m: float = 500.0
    annulus_inner_rc: float = 1.5
    annulus_outer_rc: float = 3.0
    min_annulus_points: int = 20
    trimmed_percentiles: tuple[float, float] = (10.0, 90.0)

    @property
    def partition_root(self) -> Path:
        return self.cache_root / "file_partitions"

    @property
    def monthly_climatology_path(self) -> Path:
        return self.cache_root / "monthly_climatology.npz"

    @property
    def background_table_path(self) -> Path:
        return self.cache_root / "eddy_day_background.parquet"


def clean(values):
    values = np.asarray(values, dtype=float)
    return np.where(np.abs(values) > 1e30, np.nan, values)


def model_file_for_day(day: int, model_root: Path) -> Path:
    fnumber = 1461 + ((int(day) - 1462) // 30) * 30
    return model_root / f"outer_avg_{fnumber:05d}.nc"


def physical_sigma_depths_surface_first(z_r):
    """Return positive sigma-centre depths ordered from surface to bottom."""

    depth = np.abs(np.asarray(z_r, dtype=float))
    first = float(np.nanmedian(depth[..., 0]))
    last = float(np.nanmedian(depth[..., -1]))
    if first > last:
        depth = depth[..., ::-1]
    if float(np.nanmedian(depth[..., 0])) >= float(np.nanmedian(depth[..., -1])):
        raise ValueError("Could not establish surface-to-bottom z_r ordering.")
    return depth


def sigma_layer_weights(z_r, h, depth_max_m=500.0):
    """Approximate layer thickness inside 0-depth_max from sigma centres.

    ``z_r`` supplies physical centre depths at every x/y location. Interfaces
    are reconstructed midway between adjacent centres, with zero at the
    surface and bathymetry at the bottom. The interfaces are clipped at the
    requested physical depth before thicknesses are calculated.
    """

    depth = physical_sigma_depths_surface_first(z_r)
    h = np.asarray(h, dtype=float)
    edges = np.empty(depth.shape[:-1] + (depth.shape[-1] + 1,), dtype=float)
    edges[..., 0] = 0.0
    edges[..., 1:-1] = 0.5 * (depth[..., :-1] + depth[..., 1:])
    edges[..., -1] = h
    edges = np.maximum.accumulate(edges, axis=2)
    clipped = np.minimum(edges, float(depth_max_m))
    weights = np.maximum(np.diff(clipped, axis=2), 0.0)
    weights[~np.isfinite(depth)] = 0.0
    weights[~np.isfinite(h)] = 0.0
    return weights


def analysis_domain_mask(eddy_table, x_grid, y_grid, outer_factor=3.0):
    """Union of rectangular outer-annulus windows used by selected eddies."""

    mask = np.zeros((len(x_grid), len(y_grid)), dtype=bool)
    for row in eddy_table[["xc", "yc", "Rc"]].dropna().itertuples(index=False):
        radius = float(outer_factor * row.Rc)
        if not np.isfinite(radius) or radius <= 0:
            continue
        i0 = max(0, int(np.searchsorted(x_grid, row.xc - radius, side="left")))
        i1 = min(len(x_grid), int(np.searchsorted(x_grid, row.xc + radius, side="right")))
        j0 = max(0, int(np.searchsorted(y_grid, row.yc - radius, side="left")))
        j1 = min(len(y_grid), int(np.searchsorted(y_grid, row.yc + radius, side="right")))
        mask[i0:i1, j0:j1] = True
    return mask


def required_surface_sigma_levels(weights, domain_mask):
    """Deepest contiguous surface sigma index needed over sampled cells.

    Shallow columns may require every sigma level to represent their full
    water column above 500 m. In that case this correctly returns all levels.
    """

    active = (weights > 0) & np.asarray(domain_mask, bool)[..., None]
    needed = np.flatnonzero(active.any(axis=(0, 1)))
    if needed.size == 0:
        raise ValueError("No sigma levels intersect the analysis domain.")
    if needed[0] != 0:
        raise ValueError("Expected required sigma levels to begin at the surface.")
    return int(needed.max() + 1)


def depth_weighted_mean(values_surface_first, weights):
    values = clean(values_surface_first)
    valid = np.isfinite(values) & (weights > 0)
    numerator = np.sum(np.where(valid, values * weights, 0.0), axis=2)
    denominator = np.sum(np.where(valid, weights, 0.0), axis=2)
    return np.divide(
        numerator,
        denominator,
        out=np.full(numerator.shape, np.nan),
        where=denominator > 0,
    )


def read_velocity_fields(dataset, t, weights, k_read):
    """Read geographic surface and thickness-weighted upper-ocean velocity."""

    u_surface = clean(dataset["u_eastward"][t, -1, :, :].T)
    v_surface = clean(dataset["v_northward"][t, -1, :, :].T)

    # ROMS vertical storage is bottom-to-surface. Read only the required
    # upper levels, transpose to x/y/s, then reverse to surface-to-bottom.
    u_sigma = np.flip(clean(dataset["u_eastward"][t, -k_read:, :, :].T), axis=2)
    v_sigma = np.flip(clean(dataset["v_northward"][t, -k_read:, :, :].T), axis=2)
    local_weights = weights[..., :k_read]
    u_500 = depth_weighted_mean(u_sigma, local_weights)
    v_500 = depth_weighted_mean(v_sigma, local_weights)
    return u_surface, v_surface, u_500, v_500


def local_annulus_slices(row, x_grid, y_grid, config):
    outer = float(config.annulus_outer_rc * row.Rc)
    inner = float(config.annulus_inner_rc * row.Rc)
    if not np.isfinite(outer) or outer <= 0:
        return None
    i0 = max(0, int(np.searchsorted(x_grid, row.xc - outer, side="left")))
    i1 = min(len(x_grid), int(np.searchsorted(x_grid, row.xc + outer, side="right")))
    j0 = max(0, int(np.searchsorted(y_grid, row.yc - outer, side="left")))
    j1 = min(len(y_grid), int(np.searchsorted(y_grid, row.yc + outer, side="right")))
    if i1 <= i0 or j1 <= j0:
        return None
    return slice(i0, i1), slice(j0, j1), inner, outer


def annulus_backgrounds(row, fields, x_grid, y_grid, config):
    """Calculate one cropped annulus mask and apply it to all fields."""

    bounds = local_annulus_slices(row, x_grid, y_grid, config)
    if bounds is None:
        return [np.nan] * len(fields)
    islice, jslice, inner, outer = bounds
    x_local = x_grid[islice][:, None]
    y_local = y_grid[jslice][None, :]
    distance = np.hypot(x_local - row.xc, y_local - row.yc)
    annulus = (distance >= inner) & (distance <= outer)
    results = []
    for field in fields:
        values = field[islice, jslice]
        use = annulus & np.isfinite(values)
        if int(use.sum()) < config.min_annulus_points:
            results.append(np.nan)
            continue
        selected = values[use]
        low, high = np.nanpercentile(selected, config.trimmed_percentiles)
        selected = selected[(selected >= low) & (selected <= high)]
        results.append(float(np.nanmedian(selected)))
    return results


def _partition_paths(model_path, config):
    stem = model_path.stem
    return (
        config.partition_root / "annulus" / f"{stem}.parquet",
        config.partition_root / "climatology" / f"{stem}.npz",
    )


def _save_climatology_partial(path, sums, counts):
    payload = {}
    for (month, name), array in sums.items():
        payload[f"sum_m{month:02d}_{name}"] = array
        payload[f"count_m{month:02d}_{name}"] = counts[(month, name)]
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)


def process_model_file(model_path, eddy_rows, x_grid, y_grid, weights, k_read, config):
    """Process one model file once, producing restartable partial outputs."""

    annulus_path, climatology_path = _partition_paths(model_path, config)
    if annulus_path.exists() and climatology_path.exists():
        return str(annulus_path), str(climatology_path), "cached"

    day_groups = {int(day): part for day, part in eddy_rows.groupby("Day", sort=False)}
    sums, counts, output_rows = {}, {}, []
    with nc.Dataset(model_path) as dataset:
        days = np.rint(clean(dataset["ocean_time"][:]) / SECONDS_PER_DAY).astype(int)
        for t, day in enumerate(days):
            fields = read_velocity_fields(dataset, t, weights, k_read)
            month = (DATE_ORIGIN + pd.Timedelta(days=int(day))).month
            for name, field in zip(FIELD_NAMES, fields):
                key = (month, name)
                if key not in sums:
                    sums[key] = np.zeros(field.shape, dtype=np.float64)
                    counts[key] = np.zeros(field.shape, dtype=np.int16)
                valid = np.isfinite(field)
                sums[key][valid] += field[valid]
                counts[key][valid] += 1

            for row in day_groups.get(int(day), pd.DataFrame()).itertuples(index=False):
                values = annulus_backgrounds(row, fields, x_grid, y_grid, config)
                output_rows.append({
                    "Eddy": row.Eddy,
                    "Day": int(day),
                    "month": month,
                    "ic": int(row.ic),
                    "jc": int(row.jc),
                    "ann_surface_east_ms": values[0],
                    "ann_surface_north_ms": values[1],
                    "ann_500_east_ms": values[2],
                    "ann_500_north_ms": values[3],
                })

    annulus_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output_rows, columns=ANNULUS_COLUMNS).to_parquet(annulus_path, index=False)
    _save_climatology_partial(climatology_path, sums, counts)
    return str(annulus_path), str(climatology_path), "computed"


def reduce_monthly_climatology(partial_paths, output_path):
    totals, counts = {}, {}
    for path in partial_paths:
        with np.load(path) as partial:
            for key in partial.files:
                if not key.startswith("sum_"):
                    continue
                suffix = key.removeprefix("sum_")
                count_key = f"count_{suffix}"
                if suffix not in totals:
                    totals[suffix] = np.zeros_like(partial[key], dtype=np.float64)
                    counts[suffix] = np.zeros_like(partial[count_key], dtype=np.int64)
                totals[suffix] += partial[key]
                counts[suffix] += partial[count_key]
    climatology = {
        suffix: np.divide(total, counts[suffix], out=np.full_like(total, np.nan), where=counts[suffix] > 0)
        for suffix, total in totals.items()
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **climatology)
    return climatology


def attach_climatology(annulus_table, climatology):
    out = annulus_table.copy()
    for name in FIELD_NAMES:
        values = np.full(len(out), np.nan)
        for month in range(1, 13):
            use = out["month"].eq(month).to_numpy()
            if not use.any():
                continue
            key = f"m{month:02d}_{name}"
            field = climatology[key]
            values[use] = field[out.loc[use, "ic"].astype(int), out.loc[use, "jc"].astype(int)]
        component = "east" if name.startswith("u_") else "north"
        level = "surface" if name.endswith("surface") else "500"
        out[f"clim_{level}_{component}_ms"] = values
    return out


def build_background_cache(eddy_table, grid, config=BackgroundConfig(), workers=4):
    """Build all backgrounds in one restartable, file-parallel pass."""

    from joblib import Parallel, delayed

    selected = eddy_table.copy()
    domain = analysis_domain_mask(selected, grid.x_grid, grid.y_grid, config.annulus_outer_rc)
    weights = sigma_layer_weights(grid.z_r, grid.h, config.depth_max_m)
    k_read = required_surface_sigma_levels(weights, domain)
    print(f"Reading {k_read}/{weights.shape[2]} surface sigma levels.")
    if k_read == weights.shape[2]:
        print("All sigma levels are required because sampled shallow columns lie entirely above the depth limit.")

    grouped = {}
    for day, part in selected.groupby("Day", sort=False):
        grouped.setdefault(model_file_for_day(int(day), config.model_root), []).append(part)
    file_rows = {path: pd.concat(parts, ignore_index=True) for path, parts in grouped.items()}

    # The monthly climatology must include the complete archive, including
    # months/days on which no selected eddy is present. Empty eddy tables still
    # produce climatology partials but no annulus rows.
    model_paths = sorted(config.model_root.glob("outer_avg_*.nc"))
    if not model_paths:
        raise FileNotFoundError(f"No outer_avg_*.nc files found in {config.model_root}")
    empty_rows = selected.iloc[0:0].copy()

    results = Parallel(n_jobs=workers, prefer="processes", verbose=10)(
        delayed(process_model_file)(
            path, rows, grid.x_grid, grid.y_grid, weights, k_read, config
        )
        for path in model_paths
        for rows in [file_rows.get(path, empty_rows)]
    )
    annulus_paths = [result[0] for result in results]
    climatology_paths = [result[1] for result in results]
    climatology = reduce_monthly_climatology(climatology_paths, config.monthly_climatology_path)
    annulus = pd.concat([pd.read_parquet(path) for path in annulus_paths], ignore_index=True)
    background = attach_climatology(annulus, climatology)
    if background.duplicated(["Eddy", "Day"]).any():
        raise ValueError("Background output contains duplicate Eddy-Day rows.")
    config.cache_root.mkdir(parents=True, exist_ok=True)
    background.to_parquet(config.background_table_path, index=False)
    return background


def load_background_cache(config=BackgroundConfig()):
    return pd.read_parquet(config.background_table_path)
