"""Helpers unique to the EAC eddy-tilt mechanism notebooks.

Tilt geometry is never re-estimated here. ``TiltDis`` and ``TiltDir`` are
treated as authoritative measurements loaded by ``seacofs_tilt_tools``.
Existing project helpers remain the source of grid, PV, tracking, bootstrap,
and background-flow calculations.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
METRES_PER_KM = 1000.0


def require_tilt_measurements(df: pd.DataFrame) -> pd.DataFrame:
    """Validate, but never calculate, the existing tilt measurements."""

    missing = {"Eddy", "Day", "Cyc", "TiltDis", "TiltDir"} - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")
    out = df.copy()
    out.loc[~out["TiltDir"].between(0, 360, inclusive="left"), "TiltDir"] = np.nan
    out.loc[out["TiltDis"] < 0, "TiltDis"] = np.nan
    return out


def add_tilt_components(df: pd.DataFrame) -> pd.DataFrame:
    """Encode measured tilt as east/north components (bearing is from north)."""

    out = require_tilt_measurements(df)
    theta = np.deg2rad(out["TiltDir"])
    out["tilt_east_km"] = out["TiltDis"] * np.sin(theta)
    out["tilt_north_km"] = out["TiltDis"] * np.cos(theta)
    return out


def bearing_from_east_north(east, north):
    return np.degrees(np.arctan2(east, north)) % 360.0


def signed_angle_difference(observed, reference):
    """Return observed minus reference in [-180, 180) degrees."""

    return (np.asarray(observed) - np.asarray(reference) + 180.0) % 360.0 - 180.0


def merge_one_to_one_or_many_to_one(left, right, keys=("Eddy", "Day")):
    """Merge diagnostic caches while guarding against duplicated eddy-days."""

    keys = list(keys)
    if right.duplicated(keys).any():
        duplicated = right.loc[right.duplicated(keys, keep=False), keys].head()
        raise ValueError(f"Right table contains duplicate keys:\n{duplicated}")
    return left.merge(right, on=keys, how="left", validate="many_to_one")


BACKGROUND_PAIRS = {
    "ann_200": ("ann_surface", "ann_200"),
    "ann_500": ("ann_surface", "ann_500"),
    "clim_200": ("clim_surface", "clim_200"),
    "clim_500": ("clim_surface", "clim_500"),
    "full_200": ("full_surface", "full_200"),
    "full_500": ("full_surface", "full_500"),
}


def add_background_shear(df: pd.DataFrame, pairs=BACKGROUND_PAIRS) -> pd.DataFrame:
    """Add surface-minus-depth shear and its alignment with measured tilt."""

    out = add_tilt_components(df)
    for label, (surface, deep) in pairs.items():
        east = f"{label}_shear_east_ms"
        north = f"{label}_shear_north_ms"
        required = [
            f"{surface}_east_ms", f"{surface}_north_ms",
            f"{deep}_east_ms", f"{deep}_north_ms",
        ]
        missing = set(required) - set(out.columns)
        if missing:
            raise KeyError(f"Cannot calculate {label}; missing {sorted(missing)}")
        out[east] = out[f"{surface}_east_ms"] - out[f"{deep}_east_ms"]
        out[north] = out[f"{surface}_north_ms"] - out[f"{deep}_north_ms"]
        out[f"{label}_shear_mag_ms"] = np.hypot(out[east], out[north])
        out[f"{label}_shear_dir"] = bearing_from_east_north(out[east], out[north])
        out[f"{label}_tilt_shear_offset"] = signed_angle_difference(
            out["TiltDir"], out[f"{label}_shear_dir"]
        )
    return out


def _rolling_vector_integral(part, east, north, window_days):
    """Integrate a vector velocity over trailing windows using trapezoids."""

    part = part.sort_values("Day")
    days = part["Day"].to_numpy(float)
    ue = part[east].to_numpy(float)
    vn = part[north].to_numpy(float)
    ie = np.full(len(part), np.nan)
    inn = np.full(len(part), np.nan)
    for stop in range(len(part)):
        use = (days >= days[stop] - window_days) & (days <= days[stop])
        use &= np.isfinite(days) & np.isfinite(ue) & np.isfinite(vn)
        idx = np.flatnonzero(use)
        if idx.size < 2:
            continue
        t = days[idx] * SECONDS_PER_DAY
        ie[stop] = np.trapz(ue[idx], t) / METRES_PER_KM
        inn[stop] = np.trapz(vn[idx], t) / METRES_PER_KM
    return pd.DataFrame({"east": ie, "north": inn}, index=part.index)


def add_accumulated_shear(df, method="ann_500", windows=(5, 10, 20, 30)):
    """Add trailing integrals of surface-minus-depth background velocity."""

    out = df.sort_values(["Eddy", "Day"]).copy()
    east = f"{method}_shear_east_ms"
    north = f"{method}_shear_north_ms"
    for window in windows:
        accum = pd.concat(
            [_rolling_vector_integral(part, east, north, window)
             for _, part in out.groupby("Eddy", sort=False)]
        ).sort_index()
        prefix = f"{method}_accum_{int(window)}d"
        out[f"{prefix}_east_km"] = accum["east"]
        out[f"{prefix}_north_km"] = accum["north"]
        out[f"{prefix}_mag_km"] = np.hypot(accum["east"], accum["north"])
        out[f"{prefix}_dir"] = bearing_from_east_north(accum["east"], accum["north"])
        out[f"{prefix}_offset"] = signed_angle_difference(out["TiltDir"], out[f"{prefix}_dir"])
    return out


def circular_offset_summary(df, columns, group=("Cyc",), eddy_equal=True):
    """Summarise signed angular offsets with optional equal eddy weighting."""

    records = []
    group = list(group)
    for keys, part in df.groupby(group, dropna=False, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        for column in columns:
            values = part[["Eddy", column]].dropna()
            if eddy_equal:
                from seacofs_tilt_tools import circular_mean_deg_true_north

                values = values.groupby("Eddy")[column].apply(circular_mean_deg_true_north).reset_index()
                angles = values[column].to_numpy(float)
            else:
                angles = values[column].to_numpy(float)
            if not len(angles):
                continue
            radians = np.deg2rad(angles)
            vector = np.nanmean(np.exp(1j * radians))
            record = dict(zip(group, keys))
            record.update({
                "metric": column,
                "eddies": int(values["Eddy"].nunique()),
                "mean_offset_deg": float(np.degrees(np.angle(vector))),
                "resultant_length": float(np.abs(vector)),
                "median_abs_offset_deg": float(np.nanmedian(np.abs(angles))),
            })
            records.append(record)
    return pd.DataFrame(records)


def add_stratification_proxies(df: pd.DataFrame) -> pd.DataFrame:
    """Add N2/f2 and transparent constant-N deformation/Burger proxies."""

    out = df.copy()
    if "f" not in out:
        raise KeyError("The table needs local Coriolis parameter `f`.")
    for depth in (200, 500):
        n2 = f"N2_{depth}m_s2"
        if n2 not in out:
            continue
        out[f"N2_over_f2_{depth}m"] = out[n2] / out["f"].pow(2)
        # Constant-N, flat-depth first-mode proxy: Rd ~ N H / (pi |f|).
        out[f"Rd_proxy_{depth}m_km"] = (
            np.sqrt(out[n2].clip(lower=0)) * depth
            / (np.pi * out["f"].abs()) / METRES_PER_KM
        )
        out[f"Bu_proxy_{depth}m"] = (
            out[f"Rd_proxy_{depth}m_km"] / out["Rc"]
        ).pow(2)
    return out


def add_eddy_stratification_categories(
    df: pd.DataFrame,
    n2_col="N2_500m_s2",
    *,
    class_col="N2Class",
    latitude_col=None,
    labels=("Low", "Medium", "High"),
):
    """Assign polarity-specific N2 categories using one median value per eddy.

    When ``latitude_col`` is supplied, categories are based on residuals from
    a separate quadratic N2-latitude fit for each polarity. This distinguishes
    unusually weak/strong stratification from the basin-scale latitude trend.
    """

    required = {"Cyc", "Eddy", n2_col}
    if latitude_col is not None:
        required.add(latitude_col)
    if missing := required - set(df.columns):
        raise KeyError(f"Missing stratification-category columns: {sorted(missing)}")
    columns = ["Cyc", "Eddy", n2_col]
    if latitude_col is not None:
        columns.append(latitude_col)
    aggregations = {n2_col: "median"}
    if latitude_col is not None:
        aggregations[latitude_col] = "median"
    eddy = df[columns].groupby(["Cyc", "Eddy"], as_index=False).agg(aggregations)
    value_col = f"{class_col}_value"
    eddy[value_col] = eddy[n2_col]
    if latitude_col is not None:
        for _, index in eddy.groupby("Cyc").groups.items():
            part = eddy.loc[index, [latitude_col, n2_col]].dropna()
            if len(part) < 6 or part[latitude_col].nunique() < 3:
                eddy.loc[index, value_col] = np.nan
                continue
            coefficients = np.polyfit(part[latitude_col], part[n2_col], deg=2)
            expected = np.polyval(coefficients, eddy.loc[index, latitude_col])
            eddy.loc[index, value_col] = eddy.loc[index, n2_col] - expected
    eddy[class_col] = pd.NA
    for _, index in eddy.groupby("Cyc").groups.items():
        valid = eddy.loc[index, value_col].dropna()
        if len(valid) < len(labels):
            continue
        # Ranking makes equal-sized, deterministic groups even when N2 ties.
        ranked = valid.rank(method="first")
        eddy.loc[valid.index, class_col] = pd.qcut(ranked, len(labels), labels=labels)
    eddy[class_col] = pd.Categorical(eddy[class_col], categories=labels, ordered=True)
    rename = {n2_col: f"{class_col}_eddy_median_n2"}
    if latitude_col is not None:
        rename[latitude_col] = f"{class_col}_eddy_median_latitude"
    eddy = eddy.rename(columns=rename)
    return df.merge(eddy, on=["Cyc", "Eddy"], how="left", validate="many_to_one")


def residualise_against_controls(df, columns, controls):
    """Return complete cases with each target residualised on common controls.

    Controls are fitted with an intercept using ordinary least squares. This is
    intended for transparent partial-relationship plots; clustered GEE models
    remain the primary repeated-observation inference.
    """

    columns = list(columns)
    controls = list(controls)
    required = columns + controls
    if missing := set(required) - set(df.columns):
        raise KeyError(f"Missing residualisation columns: {sorted(missing)}")
    out = df[required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if not len(out):
        raise ValueError("No complete rows remain for residualisation.")
    design = np.column_stack([
        np.ones(len(out)),
        *[out[column].to_numpy(float) for column in controls],
    ])
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Residualisation controls are rank deficient.")
    for column in columns:
        values = out[column].to_numpy(float)
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        out[f"{column}_resid"] = values - design @ coefficients
    return out


def _day_values(time):
    values = np.asarray(time.values)
    if np.issubdtype(values.dtype, np.datetime64):
        return ((values - np.datetime64("1990-01-01")) / np.timedelta64(1, "D")).astype(float)
    units = str(time.attrs.get("units", ""))
    scale = SECONDS_PER_DAY if "second" in units.lower() else 1.0
    return values.astype(float) / scale


def _depth_means(n2_columns, z_columns, depths):
    """Depth-average point columns; input arrays have shape point x vertical."""

    n2_columns = np.asarray(n2_columns, float)
    z_columns = np.asarray(z_columns, float)
    if n2_columns.shape != z_columns.shape or n2_columns.ndim != 2:
        raise ValueError("N2 and z must be matching point-by-vertical arrays.")
    output = {int(depth): np.full(n2_columns.shape[0], np.nan) for depth in depths}
    for point, (n2, z) in enumerate(zip(n2_columns, z_columns)):
        for depth in depths:
            use = np.isfinite(n2) & np.isfinite(z) & (z <= 0) & (z >= -float(depth))
            if use.sum() < 2:
                continue
            order = np.argsort(z[use])
            zz = z[use][order]
            nn = n2[use][order]
            span = zz[-1] - zz[0]
            if span > 0:
                output[int(depth)][point] = np.trapz(nn, zz) / span
    return output


@dataclass(frozen=True)
class N2CacheConfig:
    """Runtime and restart settings for the file-parallel N2 cache."""

    model_root: Path = Path("/srv/scratch/z3533156/26year_BRAN2020")
    output_path: Path = Path(
        "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
        "tilt_mechanisms/n2_eddy_day_v3_core.parquet"
    )
    grid_path: Path = Path("/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc")
    z_r_path: Path = Path(
        "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/z_r.npy"
    )
    depths: tuple[int, ...] = (200, 500)
    rho0: float = 1025.0
    point_batch_size: int = 128
    skip_existing: bool = True

    @property
    def partition_root(self):
        return self.output_path.parent / f"{self.output_path.stem}_file_partitions"


def _n2_partition_path(model_path, config):
    return config.partition_root / f"{Path(model_path).stem}.parquet"


def _atomic_parquet(table, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".tmp-{os.getpid()}.parquet")
    table.to_parquet(temporary, index=False)
    os.replace(temporary, path)


def _n2_from_density(rho_columns, z_columns, rho0=1025.0):
    """Calculate N2 and midpoint depths for point-by-level density columns."""

    rho = np.asarray(rho_columns, float)
    z = np.asarray(z_columns, float)
    if rho.shape != z.shape or rho.ndim != 2:
        raise ValueError("Density and z must be matching point-by-level arrays.")
    dz = np.diff(z, axis=1)
    drho = np.diff(rho, axis=1)
    n2 = np.divide(
        -9.81 * drho,
        float(rho0) * dz,
        out=np.full_like(drho, np.nan),
        where=np.isfinite(dz) & (np.abs(dz) > 0),
    )
    z_mid = 0.5 * (z[:, 1:] + z[:, :-1])
    return n2, z_mid


def _align_z_to_roms_levels(z_columns, sigma_coordinate=None):
    """Match z-column order to raw ROMS tracer order and return negative z."""

    z = np.asarray(z_columns, float)
    if z.ndim != 2:
        raise ValueError("z must be a point-by-level array.")
    if np.nanmedian(z) > 0:
        z = -np.abs(z)
    z_first_is_deep = np.nanmedian(np.abs(z[:, 0])) > np.nanmedian(np.abs(z[:, -1]))
    if sigma_coordinate is None:
        tracer_first_is_deep = True  # Standard ROMS s_rho storage order.
    else:
        sigma = np.asarray(sigma_coordinate, float)
        if sigma.ndim != 1 or len(sigma) != z.shape[1]:
            raise ValueError("The sigma coordinate does not match the z levels.")
        tracer_first_is_deep = sigma[0] < sigma[-1]
    if z_first_is_deep != tracer_first_is_deep:
        z = z[:, ::-1]
    aligned_first_is_deep = np.nanmedian(np.abs(z[:, 0])) > np.nanmedian(np.abs(z[:, -1]))
    if aligned_first_is_deep != tracer_first_is_deep:
        raise ValueError("Could not align z_r with the ROMS tracer vertical order.")
    return z


def _load_core_grid(grid_path):
    """Load the separable kilometre grid and ocean mask used by core means."""

    import netCDF4 as nc
    from seacofs_tilt_tools import distance_km

    with nc.Dataset(grid_path) as dataset:
        lon = np.transpose(dataset.variables["lon_rho"][:], axes=(1, 0))
        lat = np.transpose(dataset.variables["lat_rho"][:], axes=(1, 0))
        mask = np.transpose(dataset.variables["mask_rho"][:], axes=(1, 0)).astype(bool)
    j_mid = lon.shape[1] // 2
    i_mid = lon.shape[0] // 2
    x = np.insert(np.cumsum(distance_km(
        lat[:-1, j_mid], lon[:-1, j_mid], lat[1:, j_mid], lon[1:, j_mid]
    )), 0, 0.0)
    y = np.insert(np.cumsum(distance_km(
        lat[i_mid, :-1], lon[i_mid, :-1], lat[i_mid, 1:], lon[i_mid, 1:]
    )), 0, 0.0)
    return x, y, mask


def _eddy_core_indices(row, x_grid, y_grid, ocean_mask):
    """Return grid indices inside the compute_core_mean elliptical contour."""

    from types import SimpleNamespace
    from seacofs_tilt_tools import core_grid_indices

    grid = SimpleNamespace(x_grid=x_grid, y_grid=y_grid, mask_rho=ocean_mask)
    return core_grid_indices(row, grid, circle_region_flag=False)


def process_n2_model_file(model_path, file_rows, config=N2CacheConfig()):
    """Calculate all requested eddy columns in one model file and partition."""

    import xarray as xr
    import xroms

    model_path = Path(model_path)
    partition = _n2_partition_path(model_path, config)
    if config.skip_existing and partition.exists():
        return str(partition), "cached"

    geometry = ["Eddy", "Day", "ic", "jc", "xc", "yc", "Rc", "q11", "q12", "q22"]
    rows = file_rows[geometry].copy()
    rows["Day"] = rows["Day"].round().astype(int)
    rows["ic"] = rows["ic"].astype(int)
    rows["jc"] = rows["jc"].astype(int)
    rows = rows.drop_duplicates(["Eddy", "Day"])
    output = []
    # z_r.npy uses the existing project convention x/y/sigma after transpose
    # and is memory-mapped separately by each worker rather than copied.
    # z_r = np.transpose(np.load(config.z_r_path, mmap_mode="r"), (1, 2, 0))
    z_r = np.load(config.z_r_path, mmap_mode="r")
    x_grid, y_grid, ocean_mask = _load_core_grid(config.grid_path)
    with xr.open_dataset(model_path, chunks=None, cache=False) as raw:
        for variable in ("temp", "salt"):
            if variable not in raw:
                raise KeyError(f"{model_path.name} does not contain {variable!r}.")
        vertical_dim = next((dim for dim in raw["temp"].dims if dim.startswith("s_")), None)
        if vertical_dim is None:
            raise ValueError(f"Could not identify the temperature vertical dimension: {raw['temp'].dims}")
        model_days = np.rint(_day_values(raw["ocean_time"])).astype(int)
        time_for_day = {int(day): t for t, day in enumerate(model_days)}
        sigma_coordinate = raw[vertical_dim].values if vertical_dim in raw.coords else None

        def calculate_points(day, ic, jc):
            """Return depth means for arbitrary columns on one model day."""

            ic = np.asarray(ic, dtype=int)
            jc = np.asarray(jc, dtype=int)
            values = {int(depth): np.full(len(ic), np.nan) for depth in config.depths}
            for start in range(0, len(ic), config.point_batch_size):
                stop = min(start + config.point_batch_size, len(ic))
                xi = xr.DataArray(ic[start:stop], dims="point")
                eta = xr.DataArray(jc[start:stop], dims="point")
                selector = {"ocean_time": time_for_day[int(day)], "xi_rho": xi, "eta_rho": eta}
                selected = raw[["temp", "salt"]].isel(selector).load()
                temp = selected["temp"].transpose("point", vertical_dim).values.astype(float)
                salt = selected["salt"].transpose("point", vertical_dim).values.astype(float)
                temp[np.abs(temp) > 1e30] = np.nan
                salt[np.abs(salt) > 1e30] = np.nan
                z_values = np.asarray(z_r[ic[start:stop], jc[start:stop], :], float)
                z_values = _align_z_to_roms_levels(z_values, sigma_coordinate)
                if z_values.shape != temp.shape:
                    raise ValueError(
                        f"z_r columns {z_values.shape} do not match temp columns {temp.shape}."
                    )
                rho = np.asarray(xroms.density(temp, salt, z=z_values), float)
                n2_values, z_mid = _n2_from_density(rho, z_values, config.rho0)
                batch_means = _depth_means(n2_values, z_mid, config.depths)
                for depth in config.depths:
                    values[int(depth)][start:stop] = batch_means[int(depth)]
            return values

        for day, day_rows in rows.groupby("Day", sort=False):
            if int(day) not in time_for_day:
                raise KeyError(f"Day {day} was assigned to {model_path.name} but is absent from ocean_time.")
            day_rows = day_rows.reset_index(drop=True)
            core_ic, core_jc, owners = [], [], []
            core_sizes = np.zeros(len(day_rows), dtype=int)
            for position, row in enumerate(day_rows.itertuples(index=False)):
                ii, jj = _eddy_core_indices(row, x_grid, y_grid, ocean_mask)
                core_sizes[position] = len(ii)
                core_ic.extend(ii)
                core_jc.extend(jj)
                owners.extend([position] * len(ii))

            core_values = calculate_points(day, core_ic, core_jc)
            centre_values = calculate_points(
                day, day_rows["ic"].to_numpy(int), day_rows["jc"].to_numpy(int)
            )
            owners = np.asarray(owners, dtype=int)
            for position, row in enumerate(day_rows.itertuples(index=False)):
                record = {
                    "Eddy": row.Eddy,
                    "Day": int(day),
                    "N2_core_cells": int(core_sizes[position]),
                }
                for depth in config.depths:
                    depth = int(depth)
                    core = core_values[depth][owners == position]
                    valid = np.isfinite(core)
                    record[f"N2_{depth}m_core_s2"] = float(np.nanmean(core)) if valid.any() else np.nan
                    record[f"N2_{depth}m_core_std_s2"] = float(np.nanstd(core)) if valid.any() else np.nan
                    record[f"N2_{depth}m_core_valid_cells"] = int(valid.sum())
                    record[f"N2_{depth}m_core_valid_fraction"] = (
                        float(valid.mean()) if len(valid) else np.nan
                    )
                    record[f"N2_{depth}m_centre_s2"] = centre_values[depth][position]
                output.append(record)

    columns = ["Eddy", "Day", "N2_core_cells"]
    for depth in config.depths:
        depth = int(depth)
        columns.extend([
            f"N2_{depth}m_core_s2", f"N2_{depth}m_core_std_s2",
            f"N2_{depth}m_core_valid_cells", f"N2_{depth}m_core_valid_fraction",
            f"N2_{depth}m_centre_s2",
        ])
    table = pd.DataFrame(output, columns=columns).sort_values(["Eddy", "Day"])
    _atomic_parquet(table, partition)
    return str(partition), "computed"


def build_n2_cache_xroms(
    eddies: pd.DataFrame,
    model_root: Path | str | None = None,
    output_path: Path | str | None = None,
    *,
    depths=(200, 500),
    rho0=1025.0,
    workers=4,
    point_batch_size=128,
    grid_path: Path | str | None = None,
    z_r_path: Path | str | None = None,
    skip_existing=True,
):
    """Build a restartable, model-file-parallel eddy-centre N2 cache.

    Only requested temperature/salinity columns are read. ``xroms.density``
    supplies the ROMS equation of state, then N2 is evaluated from its
    documented vertical-gradient definition. Results are written once as a
    parquet cache. ``ic`` indexes xi and ``jc`` indexes eta, matching the
    transposed arrays used by ``seacofs_tilt_tools``.
    """

    from joblib import Parallel, delayed
    from beta_effect_background_flow.background_flow_tools import model_file_for_day

    required = {"Eddy", "Day", "ic", "jc", "xc", "yc", "Rc", "q11", "q12", "q22"}
    missing = required - set(eddies.columns)
    if missing:
        raise KeyError(f"Missing N2 extraction columns: {sorted(missing)}")
    defaults = N2CacheConfig()
    if int(workers) < 1:
        raise ValueError("workers must be at least 1.")
    if int(point_batch_size) < 1:
        raise ValueError("point_batch_size must be at least 1.")
    config = N2CacheConfig(
        model_root=Path(model_root) if model_root is not None else defaults.model_root,
        output_path=Path(output_path) if output_path is not None else defaults.output_path,
        grid_path=Path(grid_path) if grid_path is not None else defaults.grid_path,
        z_r_path=Path(z_r_path) if z_r_path is not None else defaults.z_r_path,
        depths=tuple(int(depth) for depth in depths),
        rho0=float(rho0),
        point_batch_size=int(point_batch_size),
        skip_existing=bool(skip_existing),
    )
    work = eddies[sorted(required)].dropna().drop_duplicates(["Eddy", "Day"]).copy()
    work["model_file"] = work["Day"].map(
        lambda day: model_file_for_day(int(round(day)), config.model_root)
    )
    file_rows = {Path(path): part.drop(columns="model_file") for path, part in work.groupby("model_file")}
    missing_files = [path for path in file_rows if not path.exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing {len(missing_files)} model files; first is {missing_files[0]}")

    results = Parallel(n_jobs=int(workers), prefer="processes", verbose=10)(
        delayed(process_n2_model_file)(path, rows, config)
        for path, rows in sorted(file_rows.items())
    )
    partition_paths = [path for path, _ in results]
    result = pd.concat([pd.read_parquet(path) for path in partition_paths], ignore_index=True)
    result = result.sort_values(["Eddy", "Day"]).reset_index(drop=True)
    if result.duplicated(["Eddy", "Day"]).any():
        raise ValueError("N2 output contains duplicate Eddy-Day rows.")
    expected = work[["Eddy", "Day"]].copy()
    expected["Day"] = expected["Day"].round().astype(int)
    key_check = expected.merge(
        result[["Eddy", "Day"]], on=["Eddy", "Day"], how="outer", indicator=True
    )
    if not key_check["_merge"].eq("both").all():
        counts = key_check["_merge"].value_counts().to_dict()
        raise ValueError(f"Reduced cache does not match requested Eddy-Day keys: {counts}")
    config.output_path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_parquet(result, config.output_path)
    return result


def add_topographic_regimes(df: pd.DataFrame, shelf_depth=2000.0, dominance_ratio=1.0):
    """Attach readable shelf and planetary/topographic PV regimes."""

    out = df.copy()
    if "Region" in out:
        out["ShelfRegime"] = np.where(out["Region"].isin(["S1", "S2"]), "on_shelf", "off_shelf")
    else:
        out["ShelfRegime"] = np.where(out["h"] <= shelf_depth, "on_shelf", "off_shelf")
    ratio = out["PV_grad_topo_mag"] / out["PV_grad_plan_mag"]
    out["PVRegime"] = np.select(
        [ratio < 1 / dominance_ratio, ratio > dominance_ratio],
        ["planetary", "topographic"],
        default="mixed",
    )
    out["topo_plan_ratio_raw"] = ratio
    return out


def add_temporal_tilt_changes(df: pd.DataFrame, max_gap_days=2.0) -> pd.DataFrame:
    """Add age, lifetime, and measured tilt-change diagnostics.

    ``TiltDis`` and ``TiltDir`` remain authoritative. Changes are calculated
    only between consecutive observations of the same eddy; rates spanning a
    gap longer than ``max_gap_days`` are set to missing.
    """

    out = add_tilt_components(df).sort_values(["Eddy", "Day"]).copy()
    grouped = out.groupby("Eddy", sort=False)
    out["age_days"] = out["Day"] - grouped["Day"].transform("min")
    out["lifetime_days"] = grouped["Day"].transform("max") - grouped["Day"].transform("min")
    out["norm_age"] = np.divide(
        out["age_days"], out["lifetime_days"],
        out=np.full(len(out), np.nan), where=out["lifetime_days"].to_numpy(float) > 0,
    )
    out["dt_days"] = grouped["Day"].diff()
    valid_step = out["dt_days"].gt(0) & out["dt_days"].le(float(max_gap_days))
    out["dtilt_km_day"] = grouped["TiltDis"].diff() / out["dt_days"]
    for component in ("east", "north"):
        out[f"dtilt_{component}_km_day"] = (
            grouped[f"tilt_{component}_km"].diff() / out["dt_days"]
        )
    previous_direction = grouped["TiltDir"].shift()
    out["tilt_turn_deg_day"] = (
        signed_angle_difference(out["TiltDir"], previous_direction) / out["dt_days"]
    )
    change_columns = [
        "dtilt_km_day", "dtilt_east_km_day", "dtilt_north_km_day", "tilt_turn_deg_day",
    ]
    out.loc[~valid_step, change_columns] = np.nan
    return out


def add_eddy_pv_exposure_class(df: pd.DataFrame, stable_fraction=2 / 3) -> pd.DataFrame:
    """Classify each eddy by the fraction of days in each daily PV regime."""

    required = {"Eddy", "PVRegime"}
    if missing := required - set(df.columns):
        raise KeyError(f"Missing PV-exposure columns: {sorted(missing)}")
    fractions = (
        pd.crosstab(df["Eddy"], df["PVRegime"], normalize="index")
        .reindex(columns=["planetary", "mixed", "topographic"], fill_value=0.0)
    )
    exposure = np.select(
        [
            fractions["planetary"] >= float(stable_fraction),
            fractions["topographic"] >= float(stable_fraction),
        ],
        ["planetary-dominated", "topographic-dominated"],
        default="transitional/mixed",
    )
    summary = fractions.add_prefix("fraction_").reset_index()
    summary["PVExposure"] = exposure
    return df.merge(summary, on="Eddy", how="left", validate="many_to_one")


def add_ekman_transport(df: pd.DataFrame, tau_east="tau_east_pa", tau_north="tau_north_pa", rho0=1025.0):
    """Calculate depth-integrated Ekman transport from geographic wind stress."""

    out = df.copy()
    required = {tau_east, tau_north, "f"}
    if missing := required - set(out.columns):
        raise KeyError(f"Missing wind columns: {sorted(missing)}")
    out["ekman_east_m2s"] = out[tau_north] / (rho0 * out["f"])
    out["ekman_north_m2s"] = -out[tau_east] / (rho0 * out["f"])
    out["ekman_dir"] = bearing_from_east_north(out["ekman_east_m2s"], out["ekman_north_m2s"])
    out["tilt_ekman_offset"] = signed_angle_difference(out["TiltDir"], out["ekman_dir"])
    return out


def standardise_within_between(df, columns, group="Eddy"):
    """Add eddy-mean and daily-anomaly versions of time-varying predictors."""

    out = df.copy()
    for column in columns:
        mean = out.groupby(group)[column].transform("mean")
        out[f"{column}_between"] = mean
        out[f"{column}_within"] = out[column] - mean
    return out
