"""Efficient background-flow diagnostics for the EAC beta-effect analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
DATE_ORIGIN = pd.Timestamp("1990-01-01")
FIELD_NAMES = (
    "u_surface", "v_surface",
    "u_200", "v_200",
    "u_500", "v_500",
)
@dataclass(frozen=True)
class BackgroundConfig:
    model_root: Path = Path("/srv/scratch/z3533156/26year_BRAN2020")
    cache_root: Path = Path(
        "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/background_flow_cache_all_eddies_v4"
    )
    shallow_depth_m: float = 200.0
    depth_max_m: float = 500.0
    climatology_window_days: int = 91

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


def analysis_domain_mask(eddy_table, x_grid, y_grid):
    """Grid cells required to sample the climatology at eddy centres."""

    mask = np.zeros((len(x_grid), len(y_grid)), dtype=bool)
    rows = eddy_table[["ic", "jc"]].dropna().astype(int)
    valid = rows.ic.between(0, len(x_grid) - 1) & rows.jc.between(0, len(y_grid) - 1)
    mask[rows.loc[valid, "ic"], rows.loc[valid, "jc"]] = True
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


def read_velocity_fields(dataset, t, weights_by_depth, k_read):
    """Read geographic surface and thickness-weighted upper-ocean velocity."""

    u_surface = clean(dataset["u_eastward"][t, -1, :, :].T)
    v_surface = clean(dataset["v_northward"][t, -1, :, :].T)

    # ROMS vertical storage is bottom-to-surface. Read only the required
    # upper levels, transpose to x/y/s, then reverse to surface-to-bottom.
    u_sigma = np.flip(clean(dataset["u_eastward"][t, -k_read:, :, :].T), axis=2)
    v_sigma = np.flip(clean(dataset["v_northward"][t, -k_read:, :, :].T), axis=2)
    output = [u_surface, v_surface]
    for depth in (200, 500):
        local_weights = weights_by_depth[depth][..., :k_read]
        output.extend([
            depth_weighted_mean(u_sigma, local_weights),
            depth_weighted_mean(v_sigma, local_weights),
        ])
    return tuple(output)


def _partition_paths(model_path, config):
    return config.partition_root / "climatology" / f"{model_path.stem}.npz"


def _save_climatology_partial(path, sums, counts):
    payload = {}
    for (month, name), array in sums.items():
        payload[f"sum_m{month:02d}_{name}"] = array
        payload[f"count_m{month:02d}_{name}"] = counts[(month, name)]
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)


def process_model_file(model_path, weights_by_depth, k_read, config):
    """Process one model file once, producing restartable partial outputs."""

    climatology_path = _partition_paths(model_path, config)
    if climatology_path.exists():
        return str(climatology_path), "cached"

    sums, counts = {}, {}
    with nc.Dataset(model_path) as dataset:
        days = np.rint(clean(dataset["ocean_time"][:]) / SECONDS_PER_DAY).astype(int)
        for t, day in enumerate(days):
            fields = read_velocity_fields(dataset, t, weights_by_depth, k_read)
            month = (DATE_ORIGIN + pd.Timedelta(days=int(day))).month
            for name, field in zip(FIELD_NAMES, fields):
                key = (month, name)
                if key not in sums:
                    sums[key] = np.zeros(field.shape, dtype=np.float64)
                    counts[key] = np.zeros(field.shape, dtype=np.int16)
                valid = np.isfinite(field)
                sums[key][valid] += field[valid]
                counts[key][valid] += 1

    _save_climatology_partial(climatology_path, sums, counts)
    return str(climatology_path), "computed"


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
    # Also form an exact full-archive mean from the original sums and counts.
    # This is count-weighted, rather than an unweighted mean of twelve monthly
    # means, so unequal month lengths and missing values are handled correctly.
    for name in FIELD_NAMES:
        month_keys = [f"m{month:02d}_{name}" for month in range(1, 13)]
        available = [key for key in month_keys if key in totals]
        if not available:
            continue
        annual_total = sum((totals[key] for key in available), start=np.zeros_like(totals[available[0]]))
        annual_count = sum((counts[key] for key in available), start=np.zeros_like(counts[available[0]]))
        climatology[f"full_{name}"] = np.divide(
            annual_total,
            annual_count,
            out=np.full_like(annual_total, np.nan),
            where=annual_count > 0,
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **climatology)
    return climatology


def _window_month_weights(day, window_days):
    """Number of days from each month in a centred climatological window."""

    if window_days < 1 or window_days % 2 != 1:
        raise ValueError("climatology_window_days must be a positive odd integer")
    date = DATE_ORIGIN + pd.Timedelta(days=int(day))
    # A non-leap reference year makes the climatological calendar circular.
    ref_day = min(date.day, 28) if date.month == 2 else date.day
    centre = pd.Timestamp(2001, date.month, ref_day)
    offsets = np.arange(-(window_days // 2), window_days // 2 + 1)
    months = pd.DatetimeIndex(centre + pd.to_timedelta(offsets, unit="D")).month
    return np.bincount(months, minlength=13)[1:].astype(float)


def attach_climatology(eddy_table, climatology, window_days=91):
    """Attach a moving seasonal climatology and full-archive mean.

    Monthly climatological fields are weighted by the number of days from
    each month falling inside a centred window about the eddy day. This is a
    smoothly moving monthly-resolution climatology, not a daily climatology.
    """

    out = eddy_table[["Eddy", "Day", "ic", "jc"]].copy()
    weights = np.vstack([_window_month_weights(day, window_days) for day in out.Day])
    for name in FIELD_NAMES:
        monthly = np.column_stack([
            climatology[f"m{month:02d}_{name}"][out.ic.to_numpy(int), out.jc.to_numpy(int)]
            for month in range(1, 13)
        ])
        valid = np.isfinite(monthly)
        numerator = np.nansum(monthly * weights, axis=1)
        denominator = np.sum(weights * valid, axis=1)
        values = np.divide(numerator, denominator, out=np.full(len(out), np.nan), where=denominator > 0)
        component = "east" if name.startswith("u_") else "north"
        level = name.split("_", 1)[1]
        out[f"clim_{level}_{component}_ms"] = values
        full_field = climatology[f"full_{name}"]
        out[f"full_{level}_{component}_ms"] = full_field[
            out["ic"].astype(int), out["jc"].astype(int)
        ]
    return out


def build_background_cache(eddy_table, grid, config=BackgroundConfig(), workers=4):
    """Build all backgrounds in one restartable, file-parallel pass."""

    from joblib import Parallel, delayed

    if not np.isclose(config.shallow_depth_m, 200.0) or not np.isclose(config.depth_max_m, 500.0):
        raise ValueError("This cache schema labels its fixed depth ranges as 0-200 m and 0-500 m.")

    selected = eddy_table.copy()
    domain = analysis_domain_mask(selected, grid.x_grid, grid.y_grid)
    weights_by_depth = {
        200: sigma_layer_weights(grid.z_r, grid.h, config.shallow_depth_m),
        500: sigma_layer_weights(grid.z_r, grid.h, config.depth_max_m),
    }
    levels_by_depth = {
        depth: required_surface_sigma_levels(weights, domain)
        for depth, weights in weights_by_depth.items()
    }
    k_read = max(levels_by_depth.values())
    nz = grid.z_r.shape[2]
    print(f"Required levels by depth: {levels_by_depth}; reading {k_read}/{nz} surface sigma levels.")
    if k_read == nz:
        print("All sigma levels are required because sampled shallow columns lie entirely above the depth limit.")

    # The climatology includes the complete archive, including days without eddies.
    model_paths = sorted(config.model_root.glob("outer_avg_*.nc"))
    if not model_paths:
        raise FileNotFoundError(f"No outer_avg_*.nc files found in {config.model_root}")
    results = Parallel(n_jobs=workers, prefer="processes", verbose=10)(
        delayed(process_model_file)(path, weights_by_depth, k_read, config)
        for path in model_paths
    )
    climatology_paths = [result[0] for result in results]
    climatology = reduce_monthly_climatology(climatology_paths, config.monthly_climatology_path)
    background = attach_climatology(selected, climatology, config.climatology_window_days)
    if background.duplicated(["Eddy", "Day"]).any():
        raise ValueError("Background output contains duplicate Eddy-Day rows.")
    config.cache_root.mkdir(parents=True, exist_ok=True)
    background.to_parquet(config.background_table_path, index=False)
    return background


def load_background_cache(config=BackgroundConfig()):
    return pd.read_parquet(config.background_table_path)
