"""ESP composites in eddy-centred or locally onshore-aligned coordinates."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd


DEFAULT_ESP_ROOT = Path("/home/z5297792/ESP_zonodo")


def _load_esp(esp_root=DEFAULT_ESP_ROOT):
    esp_root = Path(esp_root).expanduser()
    if str(esp_root) not in sys.path:
        sys.path.insert(0, str(esp_root))
    try:
        return importlib.import_module("functions")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            f"Could not import ESP functions.py from {esp_root}."
        ) from exc


def signed_angle_difference(angle, reference):
    """Return angle minus reference in [-180, 180) degrees."""
    return (np.asarray(angle) - np.asarray(reference) + 180.0) % 360.0 - 180.0


def _bearing_from_grid_vector(x, y, grid_angle):
    """Convert a model-grid vector to a true-north bearing."""
    east = x * np.cos(grid_angle) + y * np.sin(grid_angle)
    north = -x * np.sin(grid_angle) + y * np.cos(grid_angle)
    return np.degrees(np.arctan2(east, north)) % 360.0


def _geographic_basis(grid_x, grid_y, grid_angle):
    """Convert one unit model-grid basis vector to east/north components."""
    east = grid_x * np.cos(grid_angle) + grid_y * np.sin(grid_angle)
    north = -grid_x * np.sin(grid_angle) + grid_y * np.cos(grid_angle)
    norm = np.hypot(east, north)
    return east / norm, north / norm


def bathymetry_gradients(grid):
    """Return physical dh/dx and dh/dy on the model kilometre grid."""
    import seacofs_tilt_tools as tilt

    return tilt.phys_grad(
        np.asarray(grid.h, float),
        np.asarray(grid.X_grid, float) * 1e3,
        np.asarray(grid.Y_grid, float) * 1e3,
        np.asarray(grid.mask_rho),
    )


def local_onshore_basis(row, grid, dhdx, dhdy):
    """Robustly estimate the onshore direction from the core-mean depth gradient."""
    import seacofs_tilt_tools as tilt

    ii, jj = tilt.core_grid_indices(row, grid)
    if len(ii):
        gx = np.nanmean(dhdx[ii, jj])
        gy = np.nanmean(dhdy[ii, jj])
    else:
        gx, gy = np.nan, np.nan
    if not np.all(np.isfinite([gx, gy])) or np.isclose(np.hypot(gx, gy), 0):
        gx = float(dhdx[int(row.ic), int(row.jc)])
        gy = float(dhdy[int(row.ic), int(row.jc)])
    magnitude = np.hypot(gx, gy)
    if not np.isfinite(magnitude) or magnitude == 0:
        raise ValueError(
            f"Cannot determine an onshore direction for eddy {row.Eddy}, day {row.Day}."
        )
    # h increases offshore, so minus grad(h) points onshore.
    onshore = np.array([-gx, -gy], float) / magnitude
    alongshore = np.array([-onshore[1], onshore[0]], float)
    return onshore, alongshore, magnitude


def select_eddy_days(
    surface,
    vertical,
    eddy_id,
    *,
    day_range=(None, None),
    minimum_profile_depth_m=0.0,
    maximum_bathymetry_m=None,
    grid=None,
):
    """Select valid days and profiles for one eddy with transparent filters."""
    eddy_id = int(eddy_id)
    rows = surface.loc[surface.Eddy.eq(eddy_id)].copy()
    profiles = vertical.loc[vertical.Eddy.eq(eddy_id)].copy()
    required = ["Depth", "xc", "yc", "Rc", "Omega", "q11", "q12", "q22"]
    profiles = profiles.dropna(subset=required)
    if rows.empty or profiles.empty:
        raise ValueError(f"Eddy {eddy_id} is absent from the surface or vertical table.")

    start, stop = day_range
    if start is not None:
        rows = rows.loc[rows.Day.ge(start)]
        profiles = profiles.loc[profiles.Day.ge(start)]
    if stop is not None:
        rows = rows.loc[rows.Day.le(stop)]
        profiles = profiles.loc[profiles.Day.le(stop)]

    coverage = profiles.groupby("Day").Depth.max()
    valid_days = coverage.loc[coverage.ge(float(minimum_profile_depth_m))].index
    rows = rows.loc[rows.Day.isin(valid_days)].copy()
    profiles = profiles.loc[profiles.Day.isin(valid_days)].copy()

    if maximum_bathymetry_m is not None:
        if grid is None:
            raise ValueError("Pass grid when maximum_bathymetry_m is used.")
        rows["centre_bathymetry_m"] = [
            float(grid.h[int(row.ic), int(row.jc)])
            for row in rows.itertuples(index=False)
        ]
        shelf_days = rows.loc[
            rows.centre_bathymetry_m.le(float(maximum_bathymetry_m)), "Day"
        ]
        rows = rows.loc[rows.Day.isin(shelf_days)]
        profiles = profiles.loc[profiles.Day.isin(shelf_days)]

    profiles = (
        profiles.sort_values(["Day", "Depth"])
        .drop_duplicates(["Day", "Depth"])
        .reset_index(drop=True)
    )
    rows = rows.sort_values("Day").drop_duplicates("Day").reset_index(drop=True)
    common = np.intersect1d(rows.Day.unique(), profiles.Day.unique())
    rows = rows.loc[rows.Day.isin(common)].reset_index(drop=True)
    profiles = profiles.loc[profiles.Day.isin(common)].reset_index(drop=True)
    if rows.empty:
        raise ValueError("No eddy-days remain after applying the composite filters.")
    return rows, profiles


def depth_catalogue(profiles, maximum_depth_m=1000.0):
    """Summarise exact fitted-depth availability across selected days."""
    selected = profiles.loc[profiles.Depth.between(0, float(maximum_depth_m))].copy()
    total_days = selected.Day.nunique()
    catalogue = (
        selected.groupby("Depth").Day.nunique().rename("days").reset_index()
        .sort_values("Depth").reset_index(drop=True)
    )
    catalogue["fraction_of_days"] = catalogue.days / total_days
    catalogue.insert(0, "depth_index", np.arange(len(catalogue)))
    return catalogue


def choose_depths(catalogue, *, minimum_day_fraction=0.7, depth_indices=None):
    """Choose exact cached depths, either by coverage or explicit catalogue indices."""
    if depth_indices is None:
        selected = catalogue.loc[
            catalogue.fraction_of_days.ge(float(minimum_day_fraction))
        ]
    else:
        depth_indices = tuple(int(index) for index in depth_indices)
        missing = sorted(set(depth_indices) - set(catalogue.depth_index))
        if missing:
            raise IndexError(f"Depth catalogue has no indices {missing}.")
        selected = catalogue.loc[catalogue.depth_index.isin(depth_indices)]
    depths = selected.Depth.to_numpy(float)
    if not len(depths):
        raise ValueError("No exact fitted depths satisfy the selection.")
    return depths


def _basis_for_day(row, grid, frame, dhdx, dhdy):
    frame = frame.lower()
    if frame == "onshore":
        axis1, axis2, gradient_magnitude = local_onshore_basis(
            row, grid, dhdx, dhdy
        )
    elif frame in {"grid", "recentered"}:
        axis1 = np.array([1.0, 0.0])
        axis2 = np.array([0.0, 1.0])
        gradient_magnitude = np.nan
    else:
        raise ValueError("frame must be 'onshore' or 'grid'.")
    return axis1, axis2, gradient_magnitude


def build_eddy_composite(
    surface_rows,
    profiles,
    grid,
    depths,
    *,
    frame="onshore",
    half_width=100.0,
    grid_step=5.0,
    horizontal_units="km",
    esp_root=DEFAULT_ESP_ROOT,
):
    """Average daily ESP velocity fields after recentering and optional rotation.

    Daily velocity fields are evaluated first and then averaged; ESP parameters
    themselves are never averaged. Every available eddy-day receives equal
    weight at a given depth/grid cell.
    """
    esp = _load_esp(esp_root)
    depths = np.asarray(depths, float)
    coordinate = np.arange(-float(half_width), float(half_width) + grid_step / 2, grid_step)
    axis1_grid, axis2_grid = np.meshgrid(coordinate, coordinate, indexing="ij")
    shape = (len(depths), len(coordinate), len(coordinate))
    sums = {name: np.zeros(shape, float) for name in ("u1", "u2")}
    sumsquares = {name: np.zeros(shape, float) for name in ("u1", "u2")}
    counts = np.zeros(shape, int)

    dhdx, dhdy = bathymetry_gradients(grid)
    daily_records = []
    centre_records = []
    surface_lookup = surface_rows.set_index("Day", drop=False)

    for day, profile in profiles.groupby("Day", sort=True):
        if day not in surface_lookup.index:
            continue
        surface_row = surface_lookup.loc[day]
        if isinstance(surface_row, pd.DataFrame):
            surface_row = surface_row.iloc[0]
        axis1, axis2, gradient_magnitude = _basis_for_day(
            surface_row, grid, frame, dhdx, dhdy
        )
        surface_profile = profile.loc[profile.Depth.idxmin()]
        if horizontal_units.lower() == "km":
            scale_km = 1.0
        elif horizontal_units.lower() == "rc":
            scale_km = float(surface_profile.Rc)
        else:
            raise ValueError("horizontal_units must be 'km' or 'Rc'.")

        x = (
            float(surface_profile.xc)
            + scale_km * (axis1_grid * axis1[0] + axis2_grid * axis2[0])
        )
        y = (
            float(surface_profile.yc)
            + scale_km * (axis1_grid * axis1[1] + axis2_grid * axis2[1])
        )
        velocity_axis1 = _geographic_basis(*axis1, float(grid.angle))
        velocity_axis2 = _geographic_basis(*axis2, float(grid.angle))
        axis1_bearing = _bearing_from_grid_vector(*axis1, float(grid.angle))

        available = profile.set_index("Depth", drop=False)
        used_depths = 0
        for k, depth in enumerate(depths):
            if depth not in available.index:
                continue
            level = available.loc[depth]
            if isinstance(level, pd.DataFrame):
                level = level.iloc[0]
            q = np.array(
                [[level.q11, level.q12], [level.q12, level.q22]], float
            )
            u_east, v_north = esp.model_uv_at_xy(
                x * 1e3, y * 1e3,
                level.xc * 1e3, level.yc * 1e3,
                q, level.Omega, level.Rc * 1e3,
            )
            u1 = u_east * velocity_axis1[0] + v_north * velocity_axis1[1]
            u2 = u_east * velocity_axis2[0] + v_north * velocity_axis2[1]
            valid = np.isfinite(u1) & np.isfinite(u2)
            sums["u1"][k][valid] += u1[valid]
            sums["u2"][k][valid] += u2[valid]
            sumsquares["u1"][k][valid] += u1[valid] ** 2
            sumsquares["u2"][k][valid] += u2[valid] ** 2
            counts[k][valid] += 1
            offset = np.array(
                [float(level.xc) - float(surface_profile.xc),
                 float(level.yc) - float(surface_profile.yc)]
            )
            centre_records.append({
                "Day": int(day),
                "Depth": float(depth),
                "axis1_offset": float(offset @ axis1) / scale_km,
                "axis2_offset": float(offset @ axis2) / scale_km,
            })
            used_depths += 1

        record = {
            "Eddy": int(surface_row.Eddy),
            "Day": int(day),
            "Cyc": surface_row.Cyc,
            "axis1_bearing": axis1_bearing,
            "bathymetry_gradient_m_per_m": gradient_magnitude,
            "levels_used": used_depths,
        }
        for name in ("TiltDir", "TiltDis", "Ro", "h", "lat", "lon"):
            if name in surface_row.index:
                record[name] = surface_row[name]
        if "TiltDir" in record:
            record["tilt_minus_axis1_deg"] = float(
                signed_angle_difference(record["TiltDir"], axis1_bearing)
            )
        daily_records.append(record)

    mean = {}
    std = {}
    for name in ("u1", "u2"):
        mean[name] = np.divide(
            sums[name], counts, out=np.full(shape, np.nan), where=counts > 0
        )
        second = np.divide(
            sumsquares[name], counts,
            out=np.full(shape, np.nan), where=counts > 0,
        )
        std[name] = np.sqrt(np.maximum(second - mean[name] ** 2, 0))
    mean["speed"] = np.hypot(mean["u1"], mean["u2"])

    centres = pd.DataFrame(centre_records)
    centre_mean = (
        centres.groupby("Depth")[["axis1_offset", "axis2_offset"]].mean()
        .reindex(depths)
        .reset_index()
    )
    centre_std = (
        centres.groupby("Depth")[["axis1_offset", "axis2_offset"]].std()
        .reindex(depths)
        .reset_index()
    )
    daily = pd.DataFrame(daily_records).sort_values("Day").reset_index(drop=True)
    if daily.empty:
        raise ValueError("No daily ESP fields were reconstructed.")

    return SimpleNamespace(
        eddy=int(daily.Eddy.iloc[0]),
        cyc=daily.Cyc.iloc[0],
        frame=frame,
        horizontal_units=horizontal_units,
        coordinate=coordinate,
        axis1_grid=axis1_grid,
        axis2_grid=axis2_grid,
        depths=depths,
        mean=SimpleNamespace(**mean),
        std=SimpleNamespace(**std),
        counts=counts,
        daily=daily,
        centres=centres,
        centre_mean=centre_mean,
        centre_std=centre_std,
    )


def plot_daily_alignment(composite):
    """Plot the daily tilt relative to the composite's first horizontal axis."""
    import matplotlib.pyplot as plt

    daily = composite.daily
    axis_name = "onshore" if composite.frame == "onshore" else "grid-x"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    if "tilt_minus_axis1_deg" in daily:
        axes[0].plot(daily.Day, daily.tilt_minus_axis1_deg, ".-", color="tab:purple")
        axes[0].axhline(0, color="k", lw=0.8)
        axes[0].set(
            xlabel="Model day", ylabel=f"Tilt minus {axis_name} direction (degrees)",
            ylim=(-180, 180), title="Daily directional consistency",
        )
        axes[1].hist(
            daily.tilt_minus_axis1_deg.dropna(), bins=np.arange(-180, 181, 15),
            color="tab:purple", alpha=0.8,
        )
        axes[1].axvline(0, color="k", lw=0.8)
        axes[1].set(
            xlabel=f"Tilt minus {axis_name} direction (degrees)",
            ylabel="Days", title="Alignment distribution",
        )
    return fig, axes


def plot_composite_sections(composite):
    """Plot cross-axis and along-axis sections of the composite velocity."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    coordinate = composite.coordinate
    centre_index = int(np.argmin(np.abs(coordinate)))
    fields = (
        ("u1", "Axis-1 velocity"), ("u2", "Axis-2 velocity"), ("speed", "Speed")
    )
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for row, section_name in enumerate(("axis-1 section", "axis-2 section")):
        for column, (name, title) in enumerate(fields):
            values = getattr(composite.mean, name)
            section = (
                values[:, :, centre_index] if row == 0
                else values[:, centre_index, :]
            )
            ax = axes[row, column]
            if name == "speed":
                vmax = np.nanpercentile(section, 99)
                mesh = ax.pcolormesh(
                    coordinate, composite.depths, section,
                    cmap="magma", vmin=0, vmax=vmax, shading="auto",
                )
            else:
                limit = np.nanpercentile(np.abs(section), 99)
                mesh = ax.pcolormesh(
                    coordinate, composite.depths, section, cmap="RdBu_r",
                    norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                    shading="auto",
                )
            centre_column = "axis1_offset" if row == 0 else "axis2_offset"
            ax.plot(
                composite.centre_mean[centre_column], composite.depths,
                "k.-", lw=1.6, label="Mean fitted centre",
            )
            ax.set(
                xlabel=f"Relative distance ({composite.horizontal_units})",
                ylabel="Depth (m)", title=f"{section_name}: {title}",
            )
            ax.invert_yaxis()
            ax.legend(frameon=False, fontsize=8)
            fig.colorbar(mesh, ax=ax, label="m s$^{-1}$")
    fig.suptitle(
        f"{composite.cyc}{composite.eddy}: {composite.frame}-aligned ESP composite",
        fontweight="bold",
    )
    return fig, axes


def plot_composite_3d(composite, *, xy_step=3, z_step=2):
    """Plot the mean reconstructed velocity and mean fitted centre in 3-D."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(9, 7), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    step = max(1, int(xy_step))
    xx = composite.axis1_grid[::step, ::step]
    yy = composite.axis2_grid[::step, ::step]
    vmax = np.nanpercentile(composite.mean.speed, 99)
    for k in range(0, len(composite.depths), max(1, int(z_step))):
        u = composite.mean.u1[k, ::step, ::step]
        v = composite.mean.u2[k, ::step, ::step]
        valid = np.isfinite(u) & np.isfinite(v)
        z = np.full_like(xx, composite.depths[k], dtype=float)
        colours = plt.cm.viridis(
            np.clip(np.hypot(u[valid], v[valid]) / max(vmax, 1e-12), 0, 1)
        )
        ax.quiver(
            xx[valid], yy[valid], z[valid],
            u[valid], v[valid], np.zeros(valid.sum()),
            length=45, colors=colours, arrow_length_ratio=0.25, linewidth=0.6,
        )
    ax.plot(
        composite.centre_mean.axis1_offset,
        composite.centre_mean.axis2_offset,
        composite.depths,
        "r.-", lw=2.2, label="Mean fitted centre",
    )
    axis1_label = "Onshore" if composite.frame == "onshore" else "Grid x"
    axis2_label = "Alongshore" if composite.frame == "onshore" else "Grid y"
    ax.set(
        xlabel=f"{axis1_label} ({composite.horizontal_units})",
        ylabel=f"{axis2_label} ({composite.horizontal_units})",
        zlabel="Depth (m)",
        title=f"{composite.cyc}{composite.eddy}: mean ESP velocity",
    )
    ax.invert_zaxis()
    ax.view_init(elev=24, azim=-58)
    ax.legend(frameon=False)
    return fig, ax
