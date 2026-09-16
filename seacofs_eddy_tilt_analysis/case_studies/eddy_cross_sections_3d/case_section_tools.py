"""Reusable loaders, diagnostics, and plots for one-eddy 3-D case studies."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import netCDF4 as nc
import numpy as np
import pandas as pd


DEFAULT_MODEL_ROOT = Path("/srv/scratch/z3533156/26year_BRAN2020")
DEFAULT_ESP_ROOT = Path("/home/z5297792/ESP_zonodo")
DEFAULT_N2_CACHE_PATH = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset/"
    "tilt_mechanisms/n2_eddy_day_v4_potential_density_core.parquet"
)
FILL_THRESHOLD = 1e30
RHO0 = 1025.0
GRAVITY = 9.81


def choose_case(surface, vertical, eddy_id, day=None):
    """Select a requested eddy-day and return its surface row and fitted profile."""
    eddy_id = int(eddy_id)
    candidates = vertical.loc[vertical.Eddy.eq(eddy_id)].copy()
    required = ["Depth", "xc", "yc", "Rc", "Omega", "q11", "q12", "q22"]
    candidates = candidates.dropna(subset=required)
    if day is not None:
        candidates = candidates.loc[candidates.Day.eq(int(day))]
    if candidates.empty:
        suffix = f", day {day}" if day is not None else ""
        raise ValueError(f"No valid vertical ESP profile for eddy {eddy_id}{suffix}.")

    quality = (
        candidates.groupby(["Eddy", "Day"]).Depth
        .agg(levels="size", maximum_depth_m="max")
        .reset_index()
        .sort_values(["maximum_depth_m", "levels"], ascending=False)
    )
    choice = quality.iloc[0]
    profile = (
        candidates.loc[candidates.Day.eq(choice.Day)]
        .sort_values("Depth")
        .drop_duplicates("Depth")
        .reset_index(drop=True)
    )
    rows = surface.loc[
        surface.Eddy.eq(eddy_id) & surface.Day.eq(int(choice.Day))
    ]
    if rows.empty:
        raise ValueError(f"Surface table has no eddy {eddy_id}, day {int(choice.Day)}.")
    return rows.iloc[0], profile, quality


def resolve_model_file(row, model_root=DEFAULT_MODEL_ROOT):
    """Resolve the source outer_avg NetCDF named by a surface row."""
    supplied = Path(str(row.fname)).expanduser()
    candidates = [supplied, Path(model_root) / supplied.name]
    for path in candidates:
        if path.exists():
            return path
    checked = "\n".join(f"  - {path}" for path in candidates)
    raise FileNotFoundError(f"Could not locate {supplied.name}. Checked:\n{checked}")


def model_time_index(dataset, day):
    """Return the unique NetCDF time index matching an integer model day."""
    days = np.rint(
        np.asarray(dataset.variables["ocean_time"][:].data, float) / 86400.0
    ).astype(int)
    matches = np.flatnonzero(days == int(day))
    if matches.size != 1:
        raise KeyError(
            f"Expected one occurrence of day {int(day)} in "
            f"{Path(dataset.filepath()).name}; found {matches.size}."
        )
    return int(matches[0])


def _clean(values):
    values = np.asarray(values, float)
    return np.where(np.abs(values) > FILL_THRESHOLD, np.nan, values)


def _axis_for(dimensions, prefix):
    matches = [i for i, dim in enumerate(dimensions) if dim.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {prefix!r} dimension in {dimensions}; found {matches}.")
    return matches[0]


def _read_xyz(variable, time_index):
    """Read a ROMS variable and return it in x, y, vertical order."""
    dimensions = variable.dimensions
    time_axis = dimensions.index("ocean_time")
    indexer = [slice(None)] * variable.ndim
    indexer[time_axis] = int(time_index)
    values = _clean(variable[tuple(indexer)])
    remaining = tuple(dim for dim in dimensions if dim != "ocean_time")
    axes = (
        _axis_for(remaining, "xi_"),
        _axis_for(remaining, "eta_"),
        _axis_for(remaining, "s_"),
    )
    return np.transpose(values, axes), remaining[axes[2]]


def _align_z_to_sigma(z_r, sigma):
    """Align the stored z_r order with a raw ROMS sigma coordinate."""
    z = np.asarray(z_r, float)
    sigma = np.asarray(sigma, float)
    if z.shape[-1] != sigma.size:
        raise ValueError(f"z_r has {z.shape[-1]} levels but ROMS has {sigma.size}.")
    if np.nanmedian(z) > 0:
        z = -np.abs(z)
    z_first_deep = np.nanmedian(np.abs(z[..., 0])) > np.nanmedian(np.abs(z[..., -1]))
    sigma_first_deep = sigma[0] < sigma[-1]
    if z_first_deep != sigma_first_deep:
        z = z[..., ::-1]
    return z


def potential_density_and_n2(temp, salt, z, rho0=RHO0):
    """Calculate sigma0 and N2 exactly as the v4 stratification cache does."""
    try:
        import xroms
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "xroms is required for potential-density N2. Activate the same "
            "environment used to build the stratification cache."
        ) from exc

    sigma0 = np.asarray(xroms.potential_density(temp, salt, z=0.0), float)
    dz = np.diff(z, axis=-1)
    drho = np.diff(sigma0, axis=-1)
    n2 = np.divide(
        -GRAVITY * drho,
        float(rho0) * dz,
        out=np.full_like(drho, np.nan),
        where=np.isfinite(dz) & (np.abs(dz) > 0),
    )
    z_n2 = 0.5 * (z[..., 1:] + z[..., :-1])
    return sigma0, n2, z_n2


def load_native_volume(model_file, day, grid):
    """Load native temp/salt/east/north velocity and derive sigma0 and N2."""
    with nc.Dataset(model_file) as dataset:
        t = model_time_index(dataset, day)
        required = ("temp", "salt", "u_eastward", "v_northward")
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            raise KeyError(f"{Path(model_file).name} is missing {missing}.")

        temp, temp_zdim = _read_xyz(dataset["temp"], t)
        salt, salt_zdim = _read_xyz(dataset["salt"], t)
        u_native, u_zdim = _read_xyz(dataset["u_eastward"], t)
        v_native, v_zdim = _read_xyz(dataset["v_northward"], t)
        if len({temp_zdim, salt_zdim, u_zdim, v_zdim}) != 1:
            raise ValueError("Temperature, salinity and velocity use different vertical grids.")
        sigma = np.asarray(dataset.variables[temp_zdim][:], float)

    if not (temp.shape == salt.shape == u_native.shape == v_native.shape):
        raise ValueError(
            "Expected collocated rho-point temp, salt, u_eastward and v_northward; "
            f"got {temp.shape}, {salt.shape}, {u_native.shape}, {v_native.shape}."
        )
    z = _align_z_to_sigma(grid.z_r, sigma)
    if z.shape != temp.shape:
        raise ValueError(f"z_r shape {z.shape} does not match model fields {temp.shape}.")

    # Match the rotation used by the SEACOFS processing pipeline at every level.
    angle = float(grid.angle)
    u_east = v_native * np.sin(angle) + u_native * np.cos(angle)
    v_north = v_native * np.cos(angle) - u_native * np.sin(angle)
    sigma0, n2, z_n2 = potential_density_and_n2(temp, salt, z)
    return {
        "temp": temp,
        "salt": salt,
        "sigma0": sigma0,
        "N2": n2,
        "u": u_east,
        "v": v_north,
        "speed": np.hypot(u_east, v_north),
        "z": z,
        "z_N2": z_n2,
    }


def local_indices(grid, xc, yc, half_width_km):
    ii = np.flatnonzero(
        (grid.x_grid >= xc - half_width_km)
        & (grid.x_grid <= xc + half_width_km)
    )
    jj = np.flatnonzero(
        (grid.y_grid >= yc - half_width_km)
        & (grid.y_grid <= yc + half_width_km)
    )
    if ii.size < 5 or jj.size < 5:
        raise ValueError("Local case-study box contains fewer than five cells per axis.")
    return ii, jj


def crop_volume(volume, ii, jj):
    """Crop all x-y-z volume fields to the given index vectors."""
    index = np.ix_(ii, jj)
    cropped = {}
    for name, values in volume.items():
        cropped[name] = values[index[0], index[1], :]
    return cropped


def load_case(
    surface,
    profile,
    grid,
    *,
    model_root=DEFAULT_MODEL_ROOT,
    half_width_rc=2.25,
    minimum_half_width_km=80.0,
):
    """Load and crop one model eddy-day, preserving native sigma levels."""
    model_file = resolve_model_file(surface, model_root)
    half_width = max(
        float(minimum_half_width_km),
        float(half_width_rc) * float(profile.iloc[0].Rc),
    )
    ii, jj = local_indices(grid, profile.iloc[0].xc, profile.iloc[0].yc, half_width)
    volume = crop_volume(load_native_volume(model_file, int(surface.Day), grid), ii, jj)
    return SimpleNamespace(
        surface=surface,
        profile=profile,
        model_file=model_file,
        x=np.asarray(grid.x_grid[ii], float),
        y=np.asarray(grid.y_grid[jj], float),
        X=np.asarray(grid.X_grid[np.ix_(ii, jj)], float),
        Y=np.asarray(grid.Y_grid[np.ix_(ii, jj)], float),
        grid_angle=float(grid.angle),
        half_width_km=half_width,
        ii=ii,
        jj=jj,
        **volume,
    )


def load_n2_summary(eddy, day, path=DEFAULT_N2_CACHE_PATH):
    """Load the existing v4 eddy-day N2 summary (not a spatial section)."""
    from mechanism_tools import load_stratification_cache

    table = load_stratification_cache(path)
    match = table.loc[table.Eddy.eq(int(eddy)) & table.Day.eq(int(day))].copy()
    if match.empty:
        raise KeyError(f"N2 cache contains no eddy {int(eddy)}, day {int(day)}.")
    return match.iloc[0]


def _line_interval(x0, y0, dx, dy, x_bounds, y_bounds):
    """Return the parameter interval for a line contained in a rectangle."""
    lower, upper = -np.inf, np.inf
    for origin, direction, bounds in (
        (x0, dx, x_bounds), (y0, dy, y_bounds)
    ):
        if np.isclose(direction, 0):
            if not bounds[0] <= origin <= bounds[1]:
                raise ValueError("Section origin lies outside the local box.")
            continue
        crossings = sorted(((bounds[0] - origin) / direction,
                            (bounds[1] - origin) / direction))
        lower, upper = max(lower, crossings[0]), min(upper, crossings[1])
    if not np.isfinite(lower + upper) or lower >= upper:
        raise ValueError("Could not form a section through the local box.")
    return lower, upper


def transect_xy(case, orientation, points=None):
    """Return a true east-west or north-south line in model x/y coordinates."""
    orientation = orientation.lower()
    angle = float(case.grid_angle)
    if orientation in {"zonal", "east", "east-west"}:
        direction = (np.cos(angle), np.sin(angle))
    elif orientation in {"meridional", "north", "north-south"}:
        direction = (-np.sin(angle), np.cos(angle))
    else:
        raise ValueError("orientation must be 'zonal' or 'meridional'.")
    x0, y0 = float(case.profile.iloc[0].xc), float(case.profile.iloc[0].yc)
    lower, upper = _line_interval(
        x0, y0, *direction,
        (float(case.x.min()), float(case.x.max())),
        (float(case.y.min()), float(case.y.max())),
    )
    count = max(len(case.x), len(case.y)) if points is None else int(points)
    distance = np.linspace(lower, upper, max(5, count))
    return distance, x0 + distance * direction[0], y0 + distance * direction[1]


def _sample_horizontal(case, values, x_line, y_line):
    """Interpolate each native vertical level horizontally onto a transect."""
    from scipy.interpolate import RegularGridInterpolator

    sample_points = np.column_stack((x_line, y_line))
    output = np.full((values.shape[-1], len(x_line)), np.nan)
    for level in range(values.shape[-1]):
        interpolator = RegularGridInterpolator(
            (case.x, case.y), values[:, :, level], bounds_error=False,
            fill_value=np.nan,
        )
        output[level] = interpolator(sample_points)
    return output


def section_data(case, field, orientation):
    """Return true zonal/meridional section data and fitted-centre offsets."""
    if field not in {"temp", "salt", "sigma0", "N2", "u", "v", "speed"}:
        raise KeyError(f"Unsupported section field {field!r}.")
    values = getattr(case, field)
    z = case.z_N2 if field == "N2" else case.z
    orientation = orientation.lower()
    coordinate, x_line, y_line = transect_xy(case, orientation)
    section = _sample_horizontal(case, values, x_line, y_line)
    depth = -_sample_horizontal(case, z, x_line, y_line)
    dx = case.profile.xc.to_numpy(float) - float(case.profile.iloc[0].xc)
    dy = case.profile.yc.to_numpy(float) - float(case.profile.iloc[0].yc)
    angle = float(case.grid_angle)
    if orientation == "zonal":
        centre_coordinate = dx * np.cos(angle) + dy * np.sin(angle)
        label = "Eastward distance from surface centre (km)"
    elif orientation == "meridional":
        centre_coordinate = -dx * np.sin(angle) + dy * np.cos(angle)
        label = "Northward distance from surface centre (km)"
    else:
        raise ValueError("orientation must be 'zonal' or 'meridional'.")
    return coordinate, depth, section, centre_coordinate, label


def plot_sections(
    case,
    field,
    *,
    max_depth_m=1000.0,
    cmap=None,
    robust_percentiles=(2, 98),
):
    """Plot native-level zonal and meridional sections through the surface centre."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    signed = field in {"u", "v", "N2"}
    default_cmaps = {
        "N2": "RdBu_r", "u": "RdBu_r", "v": "RdBu_r",
        "temp": "turbo", "salt": "viridis", "sigma0": "cividis", "speed": "magma",
    }
    labels = {
        "N2": r"$N^2$ (s$^{-2}$)", "u": r"$u_E$ (m s$^{-1}$)",
        "v": r"$v_N$ (m s$^{-1}$)", "temp": r"Temperature ($^\circ$C)",
        "salt": "Salinity", "sigma0": r"$\sigma_0$ (kg m$^{-3}$)",
        "speed": r"Speed (m s$^{-1}$)",
    }
    fields = []
    for orientation in ("zonal", "meridional"):
        fields.append(section_data(case, field, orientation))
    finite = np.concatenate([item[2][np.isfinite(item[2])] for item in fields])
    if not finite.size:
        raise ValueError(f"No finite {field} values in the selected section.")
    lo, hi = np.nanpercentile(finite, robust_percentiles)
    if signed and lo < 0 < hi:
        limit = max(abs(lo), abs(hi))
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    else:
        norm = None

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    for ax, orientation, (coordinate, depth, values, centre, xlabel) in zip(
        axes, ("Zonal", "Meridional"), fields
    ):
        x2d = np.broadcast_to(coordinate, values.shape)
        mesh = ax.contourf(
            x2d, depth, values, levels=30, cmap=cmap or default_cmaps[field],
            norm=norm, vmin=None if norm else lo, vmax=None if norm else hi,
        )
        ax.plot(
            centre, case.profile.Depth, "k.-", lw=1.4, ms=4,
            label="Fitted ESP centre",
        )
        ax.set(
            xlabel=xlabel, ylabel="Depth (m)", ylim=(float(max_depth_m), 0),
            title=f"{orientation} section",
        )
        ax.legend(frameon=False, fontsize=8)
        fig.colorbar(mesh, ax=ax, label=labels[field], shrink=0.85)
    fig.suptitle(
        f"{case.surface.Cyc}{int(case.surface.Eddy)}, day {int(case.surface.Day)} — {field}",
        fontweight="bold",
    )
    return fig, axes


def _interpolate_columns(values, z, target_depths):
    """Interpolate a local native-z volume to fitted positive-down depths."""
    target_depths = np.asarray(target_depths, float)
    output = np.full(values.shape[:2] + (len(target_depths),), np.nan)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            depth = -np.asarray(z[i, j], float)
            column = np.asarray(values[i, j], float)
            valid = np.isfinite(depth) & np.isfinite(column)
            if valid.sum() < 2:
                continue
            order = np.argsort(depth[valid])
            output[i, j] = np.interp(
                target_depths, depth[valid][order], column[valid][order],
                left=np.nan, right=np.nan,
            )
    return np.moveaxis(output, -1, 0)


def reconstruct_esp(case, esp_root=DEFAULT_ESP_ROOT):
    """Return native velocity at fitted depths and the ESP 3-D reconstruction."""
    esp_root = Path(esp_root).expanduser()
    if str(esp_root) not in sys.path:
        sys.path.insert(0, str(esp_root))
    try:
        esp = importlib.import_module("functions")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            f"Could not import ESP functions.py from {esp_root}."
        ) from exc

    depths = case.profile.Depth.to_numpy(float)
    u_original = _interpolate_columns(case.u, case.z, depths)
    v_original = _interpolate_columns(case.v, case.z, depths)
    u_esp, v_esp = [], []
    for row in case.profile.itertuples(index=False):
        q = np.array([[row.q11, row.q12], [row.q12, row.q22]], float)
        u, v = esp.model_uv_at_xy(
            case.X * 1e3, case.Y * 1e3,
            row.xc * 1e3, row.yc * 1e3,
            q, row.Omega, row.Rc * 1e3,
        )
        u_esp.append(u)
        v_esp.append(v)
    return SimpleNamespace(
        depths=depths,
        u_original=u_original,
        v_original=v_original,
        u_esp=np.stack(u_esp),
        v_esp=np.stack(v_esp),
    )


def remove_outer_background(u, v, fraction=0.72):
    """Remove one robust outer-box velocity vector independently at each depth."""
    nx, ny = u.shape[-2:]
    xi, yi = np.ogrid[-1:1:complex(nx), -1:1:complex(ny)]
    outer = np.maximum(np.abs(xi), np.abs(yi)) >= float(fraction)
    ua, va = u.copy(), v.copy()
    background = np.full((len(u), 2), np.nan)
    for k in range(len(u)):
        background[k] = np.nanmedian(u[k][outer]), np.nanmedian(v[k][outer])
        ua[k] -= background[k, 0]
        va[k] -= background[k, 1]
    return ua, va, background


def esp_velocity_skill(comparison):
    """Return depth-resolved vector error and component correlations."""
    records = []
    for k, depth in enumerate(comparison.depths):
        uo, vo = comparison.u_original[k], comparison.v_original[k]
        ue, ve = comparison.u_esp[k], comparison.v_esp[k]
        valid = np.isfinite(uo) & np.isfinite(vo) & np.isfinite(ue) & np.isfinite(ve)
        if valid.sum() < 3:
            records.append({"Depth_m": depth, "N": int(valid.sum())})
            continue
        du, dv = ue[valid] - uo[valid], ve[valid] - vo[valid]
        rmse = np.sqrt(np.mean(du**2 + dv**2))
        rms_original = np.sqrt(np.mean(uo[valid]**2 + vo[valid]**2))

        def correlation(a, b):
            return (np.corrcoef(a, b)[0, 1]
                    if np.std(a) > 0 and np.std(b) > 0 else np.nan)

        records.append({
            "Depth_m": depth,
            "N": int(valid.sum()),
            "Vector_RMSE_m_s": rmse,
            "Vector_NRMSE": rmse / rms_original if rms_original > 0 else np.nan,
            "u_correlation": correlation(uo[valid], ue[valid]),
            "v_correlation": correlation(vo[valid], ve[valid]),
            "speed_correlation": correlation(
                np.hypot(uo[valid], vo[valid]), np.hypot(ue[valid], ve[valid])
            ),
        })
    return pd.DataFrame(records)


def plot_esp_sections(case, comparison):
    """Compare original velocity anomaly and ESP speed in both vertical planes."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    original_speed = np.hypot(comparison.u_original, comparison.v_original)
    esp_speed = np.hypot(comparison.u_esp, comparison.v_esp)
    difference = esp_speed - original_speed
    vmax = np.nanpercentile(np.r_[original_speed.ravel(), esp_speed.ravel()], 99)
    dlim = np.nanpercentile(np.abs(difference), 99)

    section_specs = [
        ("Zonal", int(np.argmin(np.abs(case.y - case.profile.iloc[0].yc))), case.x),
        ("Meridional", int(np.argmin(np.abs(case.x - case.profile.iloc[0].xc))), case.y),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for row_index, (name, index, coordinate) in enumerate(section_specs):
        if name == "Zonal":
            fields = [
                comparison.u_original[:, :, index],
                comparison.u_esp[:, :, index],
                difference[:, :, index],
            ]
            centre = case.profile.xc
            xlabel = "x (km)"
        else:
            fields = [
                comparison.v_original[:, index, :],
                comparison.v_esp[:, index, :],
                difference[:, index, :],
            ]
            centre = case.profile.yc
            xlabel = "y (km)"
        titles = (
            f"{name}: original component anomaly",
            f"{name}: ESP component",
            f"{name}: speed error",
        )
        for column, (ax, values, title) in enumerate(zip(axes[row_index], fields, titles)):
            if column < 2:
                limit = np.nanpercentile(np.abs(values), 99)
                mesh = ax.pcolormesh(
                    coordinate, comparison.depths, values,
                    cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                    shading="auto",
                )
                label = "m s$^{-1}$"
            else:
                mesh = ax.pcolormesh(
                    coordinate, comparison.depths, values,
                    cmap="RdBu_r", vmin=-dlim, vmax=dlim, shading="auto",
                )
                label = "speed error (m s$^{-1}$)"
            ax.plot(centre, case.profile.Depth, "k.-", lw=1.3)
            ax.set(title=title, xlabel=xlabel, ylabel="Depth (m)")
            ax.invert_yaxis()
            fig.colorbar(mesh, ax=ax, label=label)
    return fig, axes


def plot_velocity_3d(case, comparison, *, xy_step=7, z_step=2):
    """Show original horizontal anomaly and ESP velocity vectors in 3-D."""
    import matplotlib.pyplot as plt

    vmax = np.nanpercentile(
        np.r_[
            np.hypot(comparison.u_original, comparison.v_original).ravel(),
            np.hypot(comparison.u_esp, comparison.v_esp).ravel(),
        ],
        99,
    )

    def panel(ax, u, v, title):
        xx, yy = case.X[::xy_step, ::xy_step], case.Y[::xy_step, ::xy_step]
        for k in range(0, len(comparison.depths), int(z_step)):
            uu = u[k, ::xy_step, ::xy_step]
            vv = v[k, ::xy_step, ::xy_step]
            valid = np.isfinite(uu) & np.isfinite(vv)
            zz = np.full_like(xx, comparison.depths[k], dtype=float)
            colours = plt.cm.viridis(
                np.clip(np.hypot(uu[valid], vv[valid]) / max(vmax, 1e-12), 0, 1)
            )
            ax.quiver(
                xx[valid], yy[valid], zz[valid],
                uu[valid], vv[valid], np.zeros(valid.sum()),
                length=55, colors=colours, arrow_length_ratio=0.25, linewidth=0.6,
            )
        ax.plot(
            case.profile.xc, case.profile.yc, case.profile.Depth,
            "r.-", lw=2, label="Fitted ESP centre",
        )
        ax.set(xlabel="x (km)", ylabel="y (km)", zlabel="Depth (m)", title=title)
        ax.invert_zaxis()
        ax.view_init(elev=24, azim=-58)
        ax.legend(frameon=False)

    fig = plt.figure(figsize=(13, 6), constrained_layout=True)
    panel(
        fig.add_subplot(121, projection="3d"),
        comparison.u_original, comparison.v_original, "SEACOFS velocity anomaly",
    )
    panel(
        fig.add_subplot(122, projection="3d"),
        comparison.u_esp, comparison.v_esp, "ESP reconstruction",
    )
    return fig
