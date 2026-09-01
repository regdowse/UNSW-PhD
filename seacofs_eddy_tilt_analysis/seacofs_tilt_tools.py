"""Shared utilities for the SEACOFS eddy tilt analysis notebooks.

The notebooks in this folder are intentionally thin: they define the question,
load the common data, then call the reusable helpers collected here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import linregress
from scipy.stats import binned_statistic_2d
from matplotlib.colors import Normalize
from scipy.ndimage import gaussian_filter1d



DEFAULT_EDDY_PATH = Path("/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/processed/eddy_dataset_processed.parquet")
DEFAULT_TILT_PATH = Path("/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/tilt/tilt_dataset.parquet")
DEFAULT_VERT_PATH = Path("/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/vertical_profiles_confirmed/profiles.parquet")
DEFAULT_GRID_PATH = Path("/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc")
DEFAULT_ZR_PATH = Path("/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/z_r.npy")

KM_PER_DAY_TO_M_PER_S = 1000.0 / 86400.0
LEVELS_LAT = [-40, -35, -30, -25]
LEVELS_LON = [150, 155, 160]


@dataclass
class Paths:
    """Centralised paths used by the notebooks."""

    eddies: Path = DEFAULT_EDDY_PATH
    tilt: Path = DEFAULT_TILT_PATH
    vert: Path = DEFAULT_VERT_PATH
    grid: Path = DEFAULT_GRID_PATH
    z_r: Path = DEFAULT_ZR_PATH


@dataclass
class Grid:
    lon_rho: np.ndarray
    lat_rho: np.ndarray
    mask_rho: np.ndarray
    h: np.ndarray
    f: np.ndarray
    angle: float
    z_r: np.ndarray
    x_grid: np.ndarray
    y_grid: np.ndarray
    X_grid: np.ndarray
    Y_grid: np.ndarray


def distance_km(lat1, lon1, lat2, lon2, earth_radius_km: float = 6357.0):
    """Great-circle distance in kilometres."""

    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return earth_radius_km * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def phys_grad(F, X, Y, mask=None):
    """Gradient of ``F`` with respect to physical coordinates ``X`` and ``Y``.

    ``X`` and ``Y`` must be 2-D coordinate arrays in the same units. The
    returned gradients are in ``F`` per ``X``/``Y`` unit. A ROMS-style ocean
    mask can be supplied where 1 is ocean and 0 is land.
    """

    x_i, x_j = np.gradient(X)
    y_i, y_j = np.gradient(Y)
    F_i, F_j = np.gradient(F)
    jacobian = x_i * y_j - x_j * y_i
    dFdx = (F_i * y_j - F_j * y_i) / jacobian
    dFdy = (-F_i * x_j + F_j * x_i) / jacobian

    bad = np.isclose(jacobian, 0)
    dFdx[bad] = np.nan
    dFdy[bad] = np.nan

    if mask is not None:
        ocean = mask.astype(bool)
        dFdx = np.where(ocean, dFdx, np.nan)
        dFdy = np.where(ocean, dFdy, np.nan)
    return dFdx, dFdy


def load_grid(grid_path: Path | str = DEFAULT_GRID_PATH, z_r_path: Path | str = DEFAULT_ZR_PATH) -> Grid:
    """Load the SEACOFS grid and build kilometre x/y coordinates."""

    dataset = nc.Dataset(grid_path)
    lon_rho = np.transpose(dataset.variables["lon_rho"], axes=(1, 0))
    lat_rho = np.transpose(dataset.variables["lat_rho"], axes=(1, 0))
    mask_rho = np.transpose(dataset.variables["mask_rho"], axes=(1, 0))
    h = np.transpose(dataset.variables["h"], axes=(1, 0))
    f = np.transpose(dataset.variables["f"], axes=(1, 0))
    angle = float(dataset.variables["angle"][0, 0])
    z_r = np.load(z_r_path) #np.transpose(np.load(z_r_path), (1, 2, 0))

    j_mid = lon_rho.shape[1] // 2
    i_mid = lon_rho.shape[0] // 2
    dx = distance_km(lat_rho[:-1, j_mid], lon_rho[:-1, j_mid], lat_rho[1:, j_mid], lon_rho[1:, j_mid])
    dy = distance_km(lat_rho[i_mid, :-1], lon_rho[i_mid, :-1], lat_rho[i_mid, 1:], lon_rho[i_mid, 1:])
    x_grid = np.insert(np.cumsum(dx), 0, 0)
    y_grid = np.insert(np.cumsum(dy), 0, 0)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid, indexing="ij")

    return Grid(lon_rho, lat_rho, mask_rho, h, f, angle, z_r, x_grid, y_grid, X_grid, Y_grid)


def read_table(path: Path | str):
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    if path.suffix.lower() in {".pkl", ".pickle"}:
        return pd.read_pickle(path)
    try:
        return pd.read_parquet(path)
    except Exception:
        return pd.read_pickle(path)


def load_tilt_tables(paths: Paths = Paths(), *, add_regions: bool = False, grid: Grid | None = None):
    """Load eddy and tilt tables, then merge ``TiltDis`` and ``TiltDir``."""

    df_eddies = read_table(paths.eddies)
    df_tilt = read_table(paths.tilt)
    df_eddies = df_eddies.merge(
        df_tilt[["Eddy", "Day", "TiltDis", "TiltDir"]],
        how="left",
        on=["Eddy", "Day"],
    )

    if add_regions:
        if grid is None:
            raise ValueError("Pass grid when add_regions=True.")
        df_eddies = add_region_labels(df_eddies, grid)

    return df_eddies, df_tilt

def load_vert(paths: Paths = Paths(), dic_form: bool = False):
    vertical = read_table(paths.vert)
    if not dic_form:
        return vertical
        
    if isinstance(vertical, dict):
        return vertical
    out = {}
    if vertical.empty:
        return out
    for (eddy, day), profile in vertical.groupby(["Eddy", "Day"], sort=False):
        out.setdefault(f"Eddy{int(eddy)}", {})[f"Day{int(day)}"] = profile.copy().reset_index(drop=True)

    return out


def add_region_labels(df: pd.DataFrame, grid: Grid) -> pd.DataFrame:
    """Attach the six SEACOFS region labels used in the original notebooks."""

    _, bin_grid = make_region_grids(grid, lat_split=-33.0, shelf_offset=80.0)
    tree = cKDTree(np.column_stack([grid.X_grid.ravel(), grid.Y_grid.ravel()]))
    _, idx = tree.query(np.column_stack([df.xc, df.yc]))
    region_map = {1: "S1", 2: "S2", 3: "U1", 4: "D1", 5: "U2", 6: "D2"}
    out = df.copy()
    out["Region"] = pd.Series(bin_grid.ravel()[idx], index=out.index).map(region_map)
    return out


def add_time_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Add day index and normalised lifetime coordinates per eddy."""

    out = df.copy()
    out["Day_idx"] = out.groupby("Eddy").cumcount()
    max_day = out.groupby("Eddy")["Day_idx"].transform("max")
    out["norm_time"] = np.where(max_day > 0, out["Day_idx"] / max_day, np.nan)
    return out


def bearing_from_xy(x, y):
    """Bearing in degrees, where 0 is north and 90 is east."""

    return (90.0 - np.degrees(np.arctan2(y, x))) % 360.0


def angle_diff_180(a, b):
    """Absolute angular difference between two bearings in degrees."""

    return np.abs((a - b + 180.0) % 360.0 - 180.0)


def add_pv_gradient_terms(df: pd.DataFrame, grid: Grid, core_mean: bool = False) -> pd.DataFrame:
    """Compute planetary, topographic, and total shallow-water PV-gradient terms."""

    out = df.copy()
    out["f"] = grid.f[out.ic, out.jc]
    if core_mean:
        out = compute_core_mean(
            out, grid,
            fixed_field=grid.h,
            colname="h"
        )
    else:
        out["h"] = grid.h[out.ic, out.jc]

    dhdx, dhdy = phys_grad(grid.h, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho)
    dh_dN = np.sin(grid.angle) * dhdx + np.cos(grid.angle) * dhdy
    dh_dE = np.cos(grid.angle) * dhdx - np.sin(grid.angle) * dhdy
    if core_mean:
        out = compute_core_mean(
            out, grid,
            fixed_field=dh_dE,
            colname="dhdx"
        )
        out = compute_core_mean(
            out, grid,
            fixed_field=dh_dN,
            colname="dhdy"
        )
    else:
        out["dhdx"] = dh_dE[out.ic, out.jc]
        out["dhdy"] = dh_dN[out.ic, out.jc]

    dfdx, dfdy = phys_grad(grid.f, grid.X_grid * 1e3, grid.Y_grid * 1e3, grid.mask_rho)
    df_dN = np.sin(grid.angle) * dfdx + np.cos(grid.angle) * dfdy
    out["beta"] = df_dN[out.ic, out.jc]

    omega_f = out["w"] + out["f"]
    out["abs_vort"] = omega_f
    out["PV"] = omega_f / out["h"]
    out["PV_grad_plan_x"] = 0.0
    out["PV_grad_plan_y"] = out["beta"] / out["h"]
    out["PV_grad_topo_x"] = -omega_f * out["dhdx"] / out["h"] ** 2
    out["PV_grad_topo_y"] = -omega_f * out["dhdy"] / out["h"] ** 2
    out["PV_grad_x"] = out["PV_grad_plan_x"] + out["PV_grad_topo_x"]
    out["PV_grad_y"] = out["PV_grad_plan_y"] + out["PV_grad_topo_y"]

    for prefix in ["PV_grad_plan", "PV_grad_topo", "PV_grad"]:
        out[f"{prefix}_mag"] = np.hypot(out[f"{prefix}_x"], out[f"{prefix}_y"])
        out[f"{prefix}_theta"] = bearing_from_xy(out[f"{prefix}_x"], out[f"{prefix}_y"])

    out["dtheta_PV_grad"] = angle_diff_180(out["TiltDir"], out["PV_grad_theta"])
    out["dtheta_PV_grad_topo"] = angle_diff_180(out["TiltDir"], out["PV_grad_topo_theta"])
    out["dtheta_PV_grad_plan"] = angle_diff_180(out["TiltDir"], out["PV_grad_plan_theta"])
    out["Ro"] = np.abs(out["w"] / out["f"])
    out["topo_plan_ratio"] = np.log(out["PV_grad_topo_mag"] / out["PV_grad_plan_mag"])
    return out


def add_top_bottom_speeds(df: pd.DataFrame, dic_vert: dict, zmax: float = 1000.0) -> pd.DataFrame:
    """Add top-centre, bottom-centre, and surface-bottom propagation diagnostics."""

    out = df.copy()
    dt = out.groupby("Eddy")["Day"].diff()

    out["dx_top"] = out.groupby("Eddy")["xc"].diff()
    out["dy_top"] = out.groupby("Eddy")["yc"].diff()
    out["EddyProp"] = np.hypot(out.dx_top, out.dy_top) / dt * KM_PER_DAY_TO_M_PER_S

    x_btm, y_btm, z_btm = [], [], []
    for row in out.itertuples():
        try:
            profile = dic_vert[f"Eddy{int(row.Eddy)}"][f"Day{int(row.Day)}"]
        except KeyError:
            profile = None

        if profile is None or len(profile) == 0:
            x_btm.append(np.nan)
            y_btm.append(np.nan)
            z_btm.append(np.nan)
            continue

        deep = profile[profile.Depth.abs() <= zmax]
        bottom = deep.iloc[-1] if len(deep) else profile.iloc[-1]
        x_btm.append(bottom.xc)
        y_btm.append(bottom.yc)
        z_btm.append(bottom.Depth)

    out["x_btm"] = x_btm
    out["y_btm"] = y_btm
    out["z_btm"] = z_btm
    out["dx_btm"] = out.groupby("Eddy")["x_btm"].diff()
    out["dy_btm"] = out.groupby("Eddy")["y_btm"].diff()
    out["btm_prop"] = np.hypot(out.dx_btm, out.dy_btm) / dt * KM_PER_DAY_TO_M_PER_S
    out["sep_km"] = np.hypot(out.x_btm - out.xc, out.y_btm - out.yc)
    out["sep_rate_ms"] = out.groupby("Eddy")["sep_km"].diff() / dt * KM_PER_DAY_TO_M_PER_S
    out["top_btm_diff"] = np.hypot(out.dx_btm - out.dx_top, out.dy_btm - out.dy_top) / dt * KM_PER_DAY_TO_M_PER_S
    return out


def binned_tilt_panel(
    ax,
    df: pd.DataFrame,
    xcol: str,
    xlabel: str,
    *,
    ylabel: str | None = None,
    xlim: tuple[float, float] | None = None,
    percentile_xlim: tuple[float, float] = (10, 90),
    styles: dict | None = None,
    scatter: bool = False,
    linfit: bool = False,
    bins: int | None = None,
):
    """Median/IQR tilt-distance panel split by AE/CE."""

    if styles is None:
        styles = {
            "AE": {"line": "darkred", "fill": "red"},
            "CE": {"line": "navy", "fill": "blue"},
        }

    data = df.copy()
    data[xcol] = np.ma.filled(np.ma.asarray(data[xcol]), np.nan)
    data["TiltDis"] = np.ma.filled(np.ma.asarray(data["TiltDis"]), np.nan)
    finite_x = data[xcol].to_numpy(float)
    finite_x = finite_x[np.isfinite(finite_x)]
    if finite_x.size == 0:
        ax.set_axis_off()
        return ax

    if bins is None:
        bins = min(30, max(8, int(np.sqrt(finite_x.size))))
    edges = np.unique(np.nanquantile(finite_x, np.linspace(0, 1, bins + 1)))
    if len(edges) < 3:
        edges = np.linspace(np.nanmin(finite_x), np.nanmax(finite_x), 9)
    centers = 0.5 * (edges[:-1] + edges[1:])

    for cyc in ["AE", "CE"]:
        sub = data.loc[data.Cyc == cyc].dropna(subset=[xcol, "TiltDis"])
        x = sub[xcol].to_numpy(float)
        y = sub["TiltDis"].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        idx = np.digitize(x, edges)
        med = np.array([np.nanmedian(y[idx == i]) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
        q25 = np.array([np.nanquantile(y[idx == i], 0.25) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
        q75 = np.array([np.nanquantile(y[idx == i], 0.75) if np.any(idx == i) else np.nan for i in range(1, len(edges))])

        if scatter:
            ax.scatter(x, y, s=1, alpha=0.08, color=styles[cyc]["fill"])
        mask = np.isfinite(med)
        ax.plot(centers[mask], med[mask], lw=2.5, color=styles[cyc]["line"], label=cyc)
        ax.fill_between(centers[mask], q25[mask], q75[mask], color=styles[cyc]["fill"], alpha=0.10)

        if linfit and ok.sum() > 3:
            lo, hi = np.nanpercentile(x, percentile_xlim)
            reg = (x >= lo) & (x <= hi)
            if reg.sum() > 2:
                m, c, *_ = linregress(x[reg], y[reg])
                xf = np.linspace(lo, hi, 200)
                ax.plot(xf, m * xf + c, "--", lw=2, color=styles[cyc]["fill"])

    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlim is None:
        ax.set_xlim(*np.nanpercentile(finite_x, percentile_xlim))
    else:
        ax.set_xlim(*xlim)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.2)
    return ax


def circular_mean_deg_true_north(deg):
    """Circular mean for true-north bearings in degrees."""

    x = np.asarray(deg, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    s = np.nanmean(np.sin(np.deg2rad(x)))
    c = np.nanmean(np.cos(np.deg2rad(x)))
    return np.rad2deg(np.arctan2(s, c)) % 360


def shared_bins(*arrays, min_bins: int = 12, max_bins: int = 500):
    """Freedman-Diaconis-like bins shared across one or more arrays."""

    vals = np.concatenate([np.asarray(a, float)[np.isfinite(a)] for a in arrays if len(np.asarray(a))])
    if vals.size == 0:
        return np.linspace(0, 1, min_bins + 1)
    q25, q75 = np.nanpercentile(vals, [25, 75])
    iqr = q75 - q25
    width = 2 * iqr / np.cbrt(vals.size) if iqr > 0 else (np.nanmax(vals) - np.nanmin(vals)) / min_bins
    if width <= 0 or not np.isfinite(width):
        width = 1.0
    n = int(np.clip(np.ceil((np.nanmax(vals) - np.nanmin(vals)) / width), min_bins, max_bins))
    return np.linspace(np.nanmin(vals), np.nanmax(vals), n + 1)


# def mirrored_hist(ax, ae, ce, bins, xlabel, *, ylabel=None,
#                   colors=("r", "b"), alpha=.8, xlim=None,
#                   normalize=False):
#     """Plot AE above zero and CE below zero using shared bins.

#     If normalize=True, each histogram is normalized independently
#     so that its bin heights sum to 1.
#     """
#     ae = np.asarray(ae, float)
#     ce = np.asarray(ce, float)

#     ae = ae[np.isfinite(ae)]
#     ce = ce[np.isfinite(ce)]

#     if xlim is not None:
#         ae = ae[(ae >= xlim[0]) & (ae <= xlim[1])]
#         ce = ce[(ce >= xlim[0]) & (ce <= xlim[1])]

#     ae_counts, edges = np.histogram(ae, bins=bins)
#     ce_counts, _ = np.histogram(ce, bins=edges)

#     if normalize:
#         ae_counts = ae_counts / ae_counts.sum()
#         ce_counts = ce_counts / ce_counts.sum()

#     centers = 0.5 * (edges[:-1] + edges[1:])
#     widths = np.diff(edges)

#     ax.bar(centers,  ae_counts, width=widths, align="center",
#            color=colors[0], alpha=alpha, label="AE")
#     ax.bar(centers, -ce_counts, width=widths, align="center",
#            color=colors[1], alpha=alpha, label="CE")

#     ax.axhline(0, color="0.3", lw=0.8)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel or ("Probability" if normalize else "Frequency"))

#     ax.spines["top"].set_visible(False)
#     ax.spines["right"].set_visible(False)

#     ylim_abs_max = max(np.abs(ax.get_ylim()))
#     ax.set_ylim(-ylim_abs_max, ylim_abs_max)

#     if xlim is not None:
#         ax.set_xlim(xlim)

#     return ax
def mirrored_hist(ax, ae, ce, bins='fd', xlabel='', *, ylabel=None,
                  colors=('r', 'b'), alpha=.45, xlim=None,
                  normalize=False, smooth=True, sigma=1.2):
    """Mirrored AE/CE histogram with optional automatic bins and smoothing."""

    ae = np.asarray(ae, float)
    ce = np.asarray(ce, float)
    ae = ae[np.isfinite(ae)]
    ce = ce[np.isfinite(ce)]

    if xlim is not None:
        ae = ae[(ae >= xlim[0]) & (ae <= xlim[1])]
        ce = ce[(ce >= xlim[0]) & (ce <= xlim[1])]

    # Shared bin edges
    edges = np.histogram_bin_edges(np.r_[ae, ce], bins=bins, range=xlim)

    ae_counts, _ = np.histogram(ae, bins=edges)
    ce_counts, _ = np.histogram(ce, bins=edges)

    if normalize:
        ae_counts = ae_counts / ae_counts.sum()
        ce_counts = ce_counts / ce_counts.sum()

    centers = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    ax.bar(centers,  ae_counts, widths, color=colors[0], alpha=alpha,
           label='AE', edgecolor='none')
    ax.bar(centers, -ce_counts, widths, color=colors[1], alpha=alpha,
           label='CE', edgecolor='none')

    if smooth:
        ae_smooth = gaussian_filter1d(ae_counts.astype(float), sigma)
        ce_smooth = gaussian_filter1d(ce_counts.astype(float), sigma)

        ax.plot(centers,  ae_smooth, color=colors[0], lw=2)
        ax.plot(centers, -ce_smooth, color=colors[1], lw=2)

    ax.axhline(0, color='0.3', lw=.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel or ('Probability' if normalize else 'Frequency'))

    ax.spines[['top', 'right']].set_visible(False)

    ymax = max(ae_counts.max(), ce_counts.max()) * 1.08
    ax.set_ylim(-ymax, ymax)

    if xlim is not None:
        ax.set_xlim(xlim)

    return ax

def choose_dir_bins_cardinal(*dfs, col: str = "TiltDir", min_bins: int = 8, max_bins: int = 36, min_avg_per_sector: int = 8):
    """Choose circular direction bins while keeping cardinal directions aligned."""

    candidates = np.array([4, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 30, 36, 45, 60, 72])
    candidates = candidates[(candidates >= min_bins) & (candidates <= max_bins)]
    candidates = candidates[candidates % 4 == 0]
    counts = []
    for df in dfs:
        x = df[col].to_numpy(float)
        counts.append(np.isfinite(x).sum())
    n = min(counts) if counts else 0
    if n == 0 or candidates.size == 0:
        k = 16
    else:
        k0 = int(np.ceil(2 * n ** (1 / 3)))
        k = candidates[np.argmin(np.abs(candidates - k0))]
        while k > candidates.min() and (n / k) < min_avg_per_sector:
            k = candidates[candidates < k].max()
    return np.linspace(0.0, 360.0, k + 1), 180.0 / k


def add_season(df: pd.DataFrame, time_col: str = "Date", season_col: str = "Season") -> pd.DataFrame:
    """Add austral seasons from a datetime-like column."""

    out = df.copy()
    month = pd.to_datetime(out[time_col]).dt.month
    out[season_col] = np.select(
        [month.isin([12, 1, 2]), month.isin([3, 4, 5]), month.isin([6, 7, 8]), month.isin([9, 10, 11])],
        ["DJF", "MAM", "JJA", "SON"],
        default=np.nan,
    )
    return out


def windrose_counts(directions_deg, magnitudes, *, mag_bins, dir_bins=None, dir_shift=None):
    """Count tilt directions and magnitudes for polar bar plots."""

    directions = np.asarray(directions_deg, float)
    magnitudes = np.asarray(magnitudes, float)
    ok = np.isfinite(directions) & np.isfinite(magnitudes)
    directions = directions[ok]
    magnitudes = magnitudes[ok]
    if directions.size == 0:
        return None
    if dir_bins is None or dir_shift is None:
        tmp = pd.DataFrame({"TiltDir": directions})
        dir_bins, dir_shift = choose_dir_bins_cardinal(tmp)
    k = len(dir_bins) - 1
    binw_deg = 360.0 / k
    angles = np.deg2rad(np.arange(k) * binw_deg)
    width = np.deg2rad(binw_deg)
    directions = np.mod(directions + dir_shift, 360.0)
    dir_idx = np.digitize(directions, dir_bins, right=False) - 1
    mag_idx = np.digitize(magnitudes, mag_bins, right=False) - 1
    counts = np.zeros((len(mag_bins) - 1, k), float)
    ok = (dir_idx >= 0) & (dir_idx < k) & (mag_idx >= 0) & (mag_idx < len(mag_bins) - 1)
    for d_i, m_i in zip(dir_idx[ok], mag_idx[ok]):
        counts[m_i, d_i] += 1
    return counts, angles, width


def plot_windrose(ax, df: pd.DataFrame, *, title: str = "", mag_bins=(0, 10, 20, 30, 40, np.inf),
                  colors=None, step=None, rlim=None, mag='TiltDis', theta='TiltDir'):
    """Draw one stacked tilt windrose."""

    if colors is None:
        cmap = "Reds" if df.Cyc.iloc[0] == "AE" else "Blues"
        colors = getattr(plt.cm, cmap)(np.linspace(0.15, 1, len(mag_bins) - 1))
    data = windrose_counts(df[theta], df[mag], mag_bins=mag_bins)
    if data is None:
        ax.set_axis_off()
        return ax
    counts, angles, width = data
    bottom = np.zeros(counts.shape[1])
    for i in range(counts.shape[0]):
        hi = "inf" if np.isinf(mag_bins[i + 1]) else f"{mag_bins[i + 1]:g}"
        ax.bar(angles, counts[i], width=width, bottom=bottom, color=colors[i], edgecolor=(0, 0, 0, 0.2), label=f"{mag_bins[i]:g}-{hi}")
        bottom += counts[i]
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    if step is not None:
        local_rmax = np.max(counts.sum(axis=0))
        local_rmax = 1 if local_rmax == 0 else local_rmax
        print(local_rmax)
        top = int(np.ceil(local_rmax / step) * step)
        rticks = np.arange(step, top + step, step)
        ax.set_rticks(rticks)
    if rlim is not None:
        ax.set_rlim(0, rlim)
    if np.any(df.Cyc.isin(["AE", "CE"])):
        ax.set_rlabel_position(225 if df.Cyc.iloc[0] == "AE" else 315)
    ax.set_title(title)
    return ax


def make_region_grids(
    grid: Grid,
    *,
    lon_split: float = 157.0,
    lat_split: float = -33.5,
    shelf_hmax: float = 4000.0,
    shelf_xmax: float = 400.0,
    shelf_lonmax: float = 154.85,
    shelf_offset: float = 80.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Create the shelf/upstream/downstream six-bin region grids."""

    region_mask = (
        (grid.h < shelf_hmax)
        & (grid.X_grid < shelf_xmax)
        & (grid.lon_rho < shelf_lonmax)
        & (grid.mask_rho == 1)
    )

    if shelf_offset:
        dx = abs(np.nanmedian(np.diff(grid.X_grid, axis=0)))
        max_shift = int(round(shelf_offset / dx)) if dx > 0 else 0
        expanded = region_mask.copy()
        for shift in range(1, max_shift + 1):
            shifted = np.roll(region_mask, shift, axis=0)
            shifted[:shift, :] = False
            expanded |= shifted
        region_mask = expanded

    bin_grid = np.full(grid.X_grid.shape, np.nan)
    bin_grid[region_mask & (grid.lat_rho >= lat_split)] = 1
    bin_grid[region_mask & (grid.lat_rho < lat_split)] = 2
    bin_grid[(~region_mask) & (grid.lon_rho < lon_split) & (grid.mask_rho == 1) & (grid.lat_rho >= lat_split)] = 3
    bin_grid[(~region_mask) & (grid.lon_rho < lon_split) & (grid.mask_rho == 1) & (grid.lat_rho < lat_split)] = 4
    bin_grid[(grid.lon_rho >= lon_split) & (grid.mask_rho == 1) & (grid.lat_rho >= lat_split)] = 5
    bin_grid[(grid.lon_rho >= lon_split) & (grid.mask_rho == 1) & (grid.lat_rho < lat_split)] = 6
    return region_mask, bin_grid


def assign_six_regions(
    df: pd.DataFrame,
    grid: Grid,
    *,
    lon_split: float = 157.0,
    lat_split: float = -33.5,
    shelf_offset: float = 80.0,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Assign eddy days to the six map bins used by the windrose map."""

    _, bin_grid = make_region_grids(grid, lon_split=lon_split, lat_split=lat_split, shelf_offset=shelf_offset)

    out = df.copy()
    tree = cKDTree(np.column_stack([grid.X_grid.ravel(), grid.Y_grid.ravel()]))
    _, idx = tree.query(np.column_stack([out.xc, out.yc]))
    out["bin_id"] = bin_grid.ravel()[idx]
    out = out.dropna(subset=["bin_id"])
    out["bin_id"] = out["bin_id"].astype(int)
    return out, bin_grid

def lat_lon_contours(ax, grid,
                     levels_lat=LEVELS_LAT,
                     levels_lon=LEVELS_LON):
    # Latitude contours
    c1 = ax.contour(
        grid.X_grid,
        grid.Y_grid,
        grid.lat_rho,
        levels=levels_lat,
        colors='k',
        linewidths=0.5
    )
    ax.clabel(
        c1,
        fmt=lambda v: f"{-v:.0f}°S",
        inline=True,
        colors='k'
    )

    # Longitude contours
    c2 = ax.contour(
        grid.X_grid,
        grid.Y_grid,
        grid.lon_rho,
        levels=levels_lon,
        colors='k',
        linewidths=0.5
    )
    ax.clabel(
        c2,
        fmt=lambda v: f"{v:.0f}°E",
        inline=True,
        colors='k'
    )
    return ax


def rose_plot(
    df_data: pd.DataFrame,
    grid: Grid,
    *,
    mag: str = "TiltDis",
    theta: str = "TiltDir",
    frac: float = 2.6,
    mag_bins=(0, 10, 20, 30, 40, np.inf),
    direction_offset: float = -20.0,
    shelf_offset: float = 80.0,
    lon_split: float = 157.0,
    lat_split: float = -33.5,
    legend_title: str = "tilt dist. (km)",
    show: bool = True,
    rtick_flag=False,
    cmaps = ['Reds', 'Blues'],
    step=100
):
    """Plot shelf windroses and regional map windrose insets for AE and CE.

    Parameters
    ----------
    df_data:
        Eddy-day table containing ``Cyc`` plus the magnitude and direction
        columns. If a ``Region`` column exists it is used for bin assignment;
        otherwise eddy centres are mapped to the grid by ``xc``/``yc``.
    grid:
        SEACOFS grid returned by :func:`load_grid`.
    mag, theta:
        Magnitude and true-north bearing columns to count in each windrose.
    direction_offset:
        Degrees added to ``theta`` before binning. The original plot used
        ``-20`` to align tilt bearings with the plotted grid orientation.
    """

    required = {"Cyc", mag, theta}
    missing = required - set(df_data.columns)
    if missing:
        raise KeyError(f"df_data is missing required columns: {sorted(missing)}")

    region_map = {1: "S1", 2: "S2", 3: "U1", 4: "D1", 5: "U2", 6: "D2"}
    bin_map = {v: k for k, v in region_map.items()}
    region_mask_grid, bin_grid = make_region_grids(
        grid,
        lon_split=lon_split,
        lat_split=lat_split,
        shelf_offset=shelf_offset,
    )

    df_plot = df_data.copy()
    # if "Region" in df_plot.columns:
    #     df_plot["bin_id"] = df_plot["Region"].map(bin_map)
    # else:
    df_plot, _ = assign_six_regions( # always rerun for regions
        df_plot,
        grid,
        lon_split=lon_split,
        lat_split=lat_split,
        shelf_offset=shelf_offset,
    )
    df_plot = df_plot.dropna(subset=["bin_id", mag, theta])
    df_plot["bin_id"] = df_plot["bin_id"].astype(int)

    colors_cmps = [
        plt.get_cmap(cmaps[0])(np.linspace(0, 1, len(mag_bins) - 1)),
        plt.get_cmap(cmaps[1])(np.linspace(0, 1, len(mag_bins) - 1)),
    ]
    
    cell_w = (grid.X_grid.max() - grid.X_grid.min()) / 3
    cell_h = (grid.Y_grid.max() - grid.Y_grid.min()) / 4
    cmap_bins = plt.cm.gist_ncar
    levels_bins = np.arange(0.5, 7.5, 1)
    norm_bins = BoundaryNorm(levels_bins, cmap_bins.N)

    def get_bin_color(b, alpha=0.25):
        color = list(cmap_bins(norm_bins(b)))
        color[-1] = alpha
        return color

    def windrose_counts_local(directions_deg, magnitudes):
        data = windrose_counts(directions_deg, magnitudes, mag_bins=mag_bins)
        if data is None:
            return None, None, None, None
        counts, angles, width = data
        return counts, angles, width, counts.shape[1]

    def add_windrose(ax, x0, y0, data, colors, rmax):
        counts, angles, width, k = data
        if counts is None:
            return None

        size = frac * min(cell_w, cell_h)
        iax = ax.inset_axes(
            [x0 - size / 2, y0 - size / 2, size, size],
            transform=ax.transData,
            projection="polar",
        )
        bottom = np.zeros(k)
        for i in range(len(mag_bins) - 1):
            iax.bar(
                angles,
                counts[i],
                width=width,
                bottom=bottom,
                color=colors[i],
                edgecolor=(0, 0, 0, 0.1),
            )
            bottom += counts[i]
        iax.set_rlim(0, rmax)
        iax.set_theta_zero_location("N")
        iax.set_theta_direction(-1)
        iax.set_xticks([])
        iax.set_yticks([])
        iax.set_frame_on(False)
        return iax

    def plot_standalone_windrose(ax, data, colors, b, title="",
                                 tick_flag=False, rmax=None, rtick_flag=True):
        counts, angles, width, k = data
        if counts is None:
            ax.set_axis_off()
            return

        ax.set_facecolor(get_bin_color(b, alpha=0.25))
        bottom = np.zeros(k)
        for i in range(len(mag_bins) - 1):
            hi = r"$\infty$" if np.isinf(mag_bins[i + 1]) else f"{mag_bins[i + 1]:g}"
            label = f"[{mag_bins[i]:g}-{hi})"
            ax.bar(
                angles,
                counts[i],
                width=width,
                bottom=bottom,
                color=colors[i],
                edgecolor=(0, 0, 0, 0.1),
                label=label,
            )
            bottom += counts[i]
        if rtick_flag:
            local_rmax = np.max(counts.sum(axis=0))
            local_rmax = 1 if local_rmax == 0 else local_rmax
            ax.set_rlim(0, rmax if rmax is not None else local_rmax + 5)
            top = int(np.ceil(local_rmax / step) * step)
            rticks = np.arange(step, top + step, step)
            ax.set_rticks(rticks)
            ax.set_yticklabels([f"{t:g}" for t in rticks], fontsize=8)
        ax.set_rlabel_position(135)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(np.deg2rad([340, 70, 160, 250]))
        ax.set_xticklabels(["N", "E", "S", "W"], fontsize=9)
        ax.set_title(title, fontsize=11)

        if tick_flag:
            ax.legend(
                title=f"{title[:2]}\n{legend_title}",
                loc="center left",
                bbox_to_anchor=(0.95, 1.5),
                frameon=True,
                fontsize=9,
                title_fontsize=9,
            )

    bin_centers = {}
    for b in range(1, 7):
        ii, jj = np.where(bin_grid == b)
        bin_centers[b] = (np.nanmean(grid.X_grid[ii, jj]), np.nanmean(grid.Y_grid[ii, jj]))

    counts = {}
    for cyc in ["AE", "CE"]:
        for b in range(1, 7):
            sub = df_plot[(df_plot.Cyc == cyc) & (df_plot.bin_id == b)]
            directions = (sub[theta].to_numpy(float) + direction_offset) % 360.0
            magnitudes = sub[mag].to_numpy(float)
            counts[(cyc, b)] = windrose_counts_local(directions, magnitudes)

    rmax_map = 0
    for cyc in ["AE", "CE"]:
        for b in [3, 4, 5, 6]:
            c = counts[(cyc, b)][0]
            if c is not None:
                rmax_map = max(rmax_map, np.max(c.sum(axis=0)))
    rmax_map = 1 if rmax_map == 0 else rmax_map

    fig = plt.figure(figsize=(20, 8), constrained_layout=False)
    gs = fig.add_gridspec(1, 4, width_ratios=[1.2, 2.6, 1.2, 2.6], wspace=0.45)
    gs_AE_small = gs[0, 0].subgridspec(2, 1)
    gs_CE_small = gs[0, 2].subgridspec(2, 1)
    small_axes = {
        ("AE", 1): fig.add_subplot(gs_AE_small[0, 0], projection="polar"),
        ("AE", 2): fig.add_subplot(gs_AE_small[1, 0], projection="polar"),
        ("CE", 1): fig.add_subplot(gs_CE_small[0, 0], projection="polar"),
        ("CE", 2): fig.add_subplot(gs_CE_small[1, 0], projection="polar"),
    }
    axs = {"AE": fig.add_subplot(gs[0, 1]), "CE": fig.add_subplot(gs[0, 3])}

    for p, cyc in enumerate(["AE", "CE"]):
        for b in [1, 2]:
            plot_standalone_windrose(
                small_axes[(cyc, b)],
                counts[(cyc, b)],
                colors_cmps[p],
                b,
                title=f"{cyc}, S{b}",
                tick_flag=b == 2,
                rtick_flag=rtick_flag
            )

    for p, cyc in enumerate(["AE", "CE"]):
        ax = axs[cyc]
        ax.contourf(
            grid.X_grid,
            grid.Y_grid,
            np.where(grid.mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5],
            colors=["k"],
            alpha=0.5,
        )
        ax.contourf(
            grid.X_grid,
            grid.Y_grid,
            bin_grid,
            levels=levels_bins,
            cmap=cmap_bins,
            norm=norm_bins,
            alpha=0.25,
        )
        # c1 = ax.contour(grid.X_grid, grid.Y_grid, grid.lat_rho, levels=LEVELS_LAT,
        #                 colors="k", linewidths=0.5, linestyles='-')
        # ax.clabel(c1, fmt=lambda v: f"{np.abs(v):.0f} °S", inline=True, colors="k")
        # c2 = ax.contour(grid.X_grid, grid.Y_grid, grid.lon_rho, levels=LEVELS_LON,
        #                 colors="k", linewidths=0.5, linestyles='-')
        # ax.clabel(c2, fmt=lambda v: f"{v:.0f} °E", inline=True, colors="k")
        ax = lat_lon_contours(ax, grid)
        ax.contour(grid.X_grid, grid.Y_grid, grid.h, levels=[4000], colors="k", linewidths=1)
        ax.contour(grid.X_grid, grid.Y_grid, region_mask_grid.astype(float), levels=[0.5],
                   colors="magenta", linewidths=2, linestyles='-')
        ax.contour(grid.X_grid, grid.Y_grid, grid.lon_rho, levels=[lon_split],
                   colors="magenta", linewidths=2, linestyles='-')
        ax.contour(
            grid.X_grid,
            grid.Y_grid,
            np.where(grid.mask_rho == 1, grid.lat_rho, np.nan),
            levels=[lat_split],
            colors="magenta",
            linewidths=2,
            linestyles='-'
        )

        for b in [3, 4, 5, 6]:
            x0, y0 = bin_centers[b]
            add_windrose(ax, x0, y0, counts[(cyc, b)], colors_cmps[p], rmax_map)

        for label, x, y in [
            ("S1", 220, 1300),
            ("S2", 120, 50),
            ("U1", 400, 1450),
            ("U2", 800, 1450),
            ("D1", 620, 750),
            ("D2", 850, 750),
        ]:
            ax.text(x, y, label, ha="center", va="center", fontsize=12, fontweight="bold")

        ax.set_aspect("equal")
        ax.set_xlim(grid.X_grid.min(), grid.X_grid.max())
        ax.set_ylim(grid.Y_grid.min(), grid.Y_grid.max())
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    for key, label in [(("AE", 1), "a)"), (("AE", 2), "b)"), (("CE", 1), "d)"), (("CE", 2), "e)")]:
        small_axes[key].text(-0.15, 1.02, label, transform=small_axes[key].transAxes,
                             fontsize=12, fontweight="bold", va="top", ha="left")
    axs["AE"].text(-0.08, 1.02, "c)", transform=axs["AE"].transAxes, fontsize=12, fontweight="bold", va="top", ha="left")
    axs["CE"].text(-0.08, 1.02, "f)", transform=axs["CE"].transAxes, fontsize=12, fontweight="bold", va="top", ha="left")
        
    if show:
        plt.show()
    return fig, {"map_axes": axs, "shelf_axes": small_axes, "counts": counts, "bin_grid": bin_grid}


def point_from_bearing(origin, distance, bearing_deg):
    """Endpoint from an origin, distance, and true-north bearing."""

    theta_rad = np.radians(bearing_deg)
    dx = distance * np.sin(theta_rad)
    dy = distance * np.cos(theta_rad)
    return origin[0] - dx, origin[1] - dy


def match_old_eddies(sample_eddies_old, df_eddies_old, df_eddies, min_overlap_frac: float = 0.5, max_mean_dist: float = np.inf):
    """Match old sample eddy IDs to the vertically checked eddy IDs."""

    matches = []
    for eddy_old in sample_eddies_old:
        old = (
            df_eddies_old.loc[df_eddies_old.Eddy == eddy_old, ["Day", "xc", "yc"]]
            .drop_duplicates("Day")
            .sort_values("Day")
        )
        old_days = set(old.Day)
        candidates = df_eddies.loc[df_eddies.Day.isin(old_days), "Eddy"].unique()
        scores = []
        for eddy_new in candidates:
            new = (
                df_eddies.loc[df_eddies.Eddy == eddy_new, ["Day", "xc", "yc"]]
                .drop_duplicates("Day")
                .sort_values("Day")
            )
            merged = old.merge(new, on="Day", suffixes=("_old", "_new"))
            if merged.empty:
                continue
            overlap_frac = len(merged) / len(old)
            mean_dist = np.hypot(merged.xc_old - merged.xc_new, merged.yc_old - merged.yc_new).mean()
            if overlap_frac >= min_overlap_frac and mean_dist <= max_mean_dist:
                scores.append((eddy_new, overlap_frac, mean_dist, len(merged)))
        if scores:
            best = sorted(scores, key=lambda x: (-x[1], x[2]))[0]
            matches.append({"old_eddy": eddy_old, "new_eddy": best[0], "overlap_frac": best[1], "mean_dist_km": best[2], "n_overlap": best[3]})
        else:
            matches.append({"old_eddy": eddy_old, "new_eddy": np.nan, "overlap_frac": 0.0, "mean_dist_km": np.nan, "n_overlap": 0})
    return pd.DataFrame(matches)

def core_grid_indices(row, grid: Grid, circle_region_flag: bool = False):
    """Return ocean-grid indices inside an eddy's core contour."""

    if circle_region_flag:
        if not (hasattr(row, "rmax") and np.isfinite(row.rmax) and row.rmax > 0):
            return np.array([], dtype=int), np.array([], dtype=int)
        q = np.eye(2)
        threshold = float(row.rmax) ** 2
    else:
        if hasattr(row, "q11") and np.isfinite(row.q11):
            q = np.array([[row.q11, row.q12], [row.q12, row.q22]], dtype=float)
        elif hasattr(row, "Q"):
            q = np.asarray(row.Q, dtype=float)
        else:
            return np.array([], dtype=int), np.array([], dtype=int)
        if q.shape != (2, 2) or not np.isfinite(q).all() or not np.isfinite(row.Rc) or row.Rc <= 0:
            return np.array([], dtype=int), np.array([], dtype=int)
        threshold = float(row.Rc) ** 2 / 2.0
    eigenvalues = np.linalg.eigvalsh(q)
    if not np.isfinite(eigenvalues).all() or eigenvalues.min() <= 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    radius = np.sqrt(threshold / eigenvalues.min())
    i0 = max(0, int(np.searchsorted(grid.x_grid, row.xc - radius, side="left")))
    i1 = min(len(grid.x_grid), int(np.searchsorted(grid.x_grid, row.xc + radius, side="right")))
    j0 = max(0, int(np.searchsorted(grid.y_grid, row.yc - radius, side="left")))
    j1 = min(len(grid.y_grid), int(np.searchsorted(grid.y_grid, row.yc + radius, side="right")))
    if i0 >= i1 or j0 >= j1:
        return np.array([], dtype=int), np.array([], dtype=int)
    ii, jj = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1), indexing="ij")
    dx = grid.x_grid[ii] - float(row.xc)
    dy = grid.y_grid[jj] - float(row.yc)
    rho2 = q[0, 0] * dx**2 + 2.0 * q[0, 1] * dx * dy + q[1, 1] * dy**2
    use = (rho2 <= threshold) & grid.mask_rho[ii, jj].astype(bool)
    return ii[use].astype(int), jj[use].astype(int)


def compute_core_mean(
    df_data: pd.DataFrame,
    grid: Grid,
    *,
    base_path=None, #"/srv/scratch/z3533156/26year_BRAN2020/
    varname=None,
    fixed_field=None,
    colname=None,
    circle_region_flag=False
):
    """
    Core-mean of either
      - a 3D field (x,y,t) loaded as <varname>_<fnumber>.npy, or
      - a fixed 2D field (x,y) passed in as fixed_field.
    """
    if fixed_field is None and (base_path is None or varname is None):
        raise ValueError("Either fixed_field OR (base_path and varname) must be provided.")
    mode_2d = fixed_field is not None
    if colname is None:
        if mode_2d:
            colname = "field_core"
        else:
            colname = f"{varname}"
    df = df_data.copy()
    chunks = []
    if mode_2d:
        field2d = np.where(grid.mask_rho, fixed_field, np.nan)
    for fname, df_loc in df.groupby("fname"):
        if not mode_2d:
            fnumber  = int(fname[-8:-3])
            base_day = fnumber + 1
            data3d = np.load(f"{base_path}/{varname}_{fnumber:05}.npy")
            data3d = np.where(grid.mask_rho[:, :, None], data3d, np.nan)
        df_loc = df_loc.copy().reset_index(drop=False)
        core_vals = np.full(len(df_loc), np.nan)
        for idx, row in enumerate(df_loc.itertuples(index=False)):
            ii, jj = core_grid_indices(row, grid, circle_region_flag=circle_region_flag)
            if not len(ii):
                continue
            if mode_2d:
                vals = field2d[ii, jj]
            else:
                t_idx = int(row.Day - base_day)
                vals = data3d[ii, jj, t_idx]
            core_vals[idx] = np.nanmean(vals)
        chunks.append(pd.DataFrame({
            "Eddy": df_loc["Eddy"].to_numpy(),
            "Day":  df_loc["Day"].to_numpy(),
            colname: core_vals
        }))
    df_core = pd.concat(chunks, ignore_index=True)
    df_out = df_data.merge(
        df_core[["Eddy", "Day", colname]],
        how="left",
        on=["Eddy", "Day"]
    )
    return df_out

def bin_edges_fd(x, xgrid, rule='fd'): # Freedman-Diaconis (fg) rationale
    n = len(x)
    if n < 2: return np.array([np.min(x), np.max(x)])
    rng = np.ptp(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    std = np.std(x, ddof=1)

    # raw width (km)
    if rule.lower() == 'fd':
        h = 2 * (iqr if iqr > 0 else 1.349*std) / (n ** (1/3))
    else:  # 'scott'
        h = 3.5 * std / (n ** (1/3))

    # fallback if degenerate
    if not np.isfinite(h) or h <= 0:
        h = rng / max(10, np.sqrt(n))

    # snap to grid spacing
    base = _grid_step(xgrid)
    h = _nice_step(h, base)

    lo = np.floor(np.min(x) / h) * h
    hi = np.ceil(np.max(x) / h) * h
    return np.arange(lo, hi + h, h)

def _grid_step(G):
    gx = np.diff(np.sort(np.unique(G.ravel())))
    return np.nanmedian(gx[gx > 0])

def _nice_step(h, base):
    s = h / base
    for k in [1, 2, 2.5, 5, 10]:
        if s <= k: return k * base
    return np.ceil(s) * base

def binned_median(x, y, v, xbins, ybins):
    ix = np.digitize(x, xbins) - 1
    iy = np.digitize(y, ybins) - 1

    nx, ny = len(xbins) - 1, len(ybins) - 1
    ok = (
        (ix >= 0) & (ix < nx) &
        (iy >= 0) & (iy < ny) &
        np.isfinite(v)
    )

    bins = {}
    for i, j, val in zip(ix[ok], iy[ok], v[ok]):
        bins.setdefault((j, i), []).append(val)

    out = np.full((ny, nx), np.nan)
    for (j, i), vals in bins.items():
        out[j, i] = np.nanmedian(vals)

    return out
   
def plot_binned_median_map(
    df_data: pd.DataFrame,
    grid: Grid,
    *,
    metric='Rc',
    vmin=0,
    vmax=120,
    rule='fd',
    levels_lat=[-40, -35, -30, -25],
    levels_lon=[150, 155, 160],
    cmaps={'AE': 'Reds', 'CE': 'Blues'},
    units='km',
    figsize=(9, 8),
    fig=None,
    axs=None,
    show=True,
    cbar_loc='top'
):
    df_data = df_data.copy()

    xbins = bin_edges_fd(
        pd.to_numeric(df_data.xc, errors='coerce').to_numpy(dtype=float),
        grid.X_grid,
        rule=rule
    )
    ybins = bin_edges_fd(
        pd.to_numeric(df_data.yc, errors='coerce').to_numpy(dtype=float),
        grid.Y_grid,
        rule=rule
    )

    norm = Normalize(vmin=vmin, vmax=vmax)

    # Create figure only if axes weren't supplied
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for ax, cyc in zip(axs, ['AE', 'CE']):

        df = (
            df_data[df_data.Cyc.eq(cyc)]
            .dropna(subset=['xc', 'yc', metric])
            .sort_values(metric, kind='mergesort', ignore_index=True)
        )

        H = binned_median(
            df.xc.to_numpy(dtype=float),
            df.yc.to_numpy(dtype=float),
            df[metric].to_numpy(dtype=float),
            xbins,
            ybins
        )

        m = ax.pcolormesh(
            xbins, ybins, H,
            cmap=cmaps[cyc],
            norm=norm,
            shading='auto',
            rasterized=True
        )

        cb = fig.colorbar(
            m,
            ax=ax,
            location=cbar_loc,
            shrink=0.9,
            pad=0.02
        )
        cb.set_label(
            fr'{cyc} median ${metric}$ ({units})',
            fontsize=12
        )
        cb.set_ticks(np.linspace(vmin, vmax, 5))

        # 4000 m contour
        ax.contour(
            grid.X_grid,
            grid.Y_grid,
            grid.h,
            levels=[4000],
            colors='k'
        )

        # Land
        ax.contourf(
            grid.X_grid,
            grid.Y_grid,
            np.where(grid.mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5],
            colors=['k'],
            alpha=0.5
        )

        ax = lat_lon_contours(ax, grid)
        # # Latitude contours
        # c1 = ax.contour(
        #     grid.X_grid,
        #     grid.Y_grid,
        #     grid.lat_rho,
        #     levels=levels_lat,
        #     colors='k',
        #     linewidths=0.5
        # )
        # ax.clabel(
        #     c1,
        #     fmt=lambda v: f"{-v:.0f}°S",
        #     inline=True,
        #     colors='k'
        # )

        # # Longitude contours
        # c2 = ax.contour(
        #     grid.X_grid,
        #     grid.Y_grid,
        #     grid.lon_rho,
        #     levels=levels_lon,
        #     colors='k',
        #     linewidths=0.5
        # )
        # ax.clabel(
        #     c2,
        #     fmt=lambda v: f"{v:.0f}°E",
        #     inline=True,
        #     colors='k'
        # )

        ax.set_xlim(15, grid.X_grid.max())
        ax.set_ylim(grid.Y_grid.min(), grid.Y_grid.max())
        ax.set_xlabel('x (km)', fontsize=11)
        ax.set_aspect('equal')

    axs[0].set_ylabel('y (km)', fontsize=11)

    if show:
        plt.tight_layout()
        plt.show()

    return fig, axs

def plot_pv_dominance(
    df: pd.DataFrame,
    grid: Grid,
    *,
    rule='fd',
    figsize=(6, 6),
    clabel='Fraction of eddy-days with $|\\nabla PV_{plan}| < |\\nabla PV_{topo}|$'
):

    xbins = bin_edges_fd(df.xc.values, grid.X_grid, rule=rule)
    ybins = bin_edges_fd(df.yc.values, grid.Y_grid, rule=rule)

    fig, axs = plt.subplots(1, 2, figsize=figsize, sharey=True, constrained_layout=True)

    for ax, cyc in zip(axs, ['AE', 'CE']):

        d = df[df.Cyc == cyc].copy()
        d['topo_dom'] = d.PV_grad_topo_mag > d.PV_grad_plan_mag

        ax.contourf(
            grid.X_grid, grid.Y_grid, np.where(grid.mask_rho == 0, 1, np.nan),
            levels=[0.5, 1.5], colors=['r' if cyc=='AE' else 'b'], alpha=.25
        )
        ax.contour(
            grid.X_grid,
            grid.Y_grid,
            grid.mask_rho,
            levels=[0.5],
            colors='r' if cyc == 'AE' else 'b',
            linewidths=2,
            zorder=15
        )

        H = binned_statistic_2d(
            d.xc, d.yc, d.topo_dom.astype(float),
            statistic='mean',
            bins=[xbins, ybins]
        ).statistic.T

        ax.contour(grid.X_grid, grid.Y_grid, grid.h, levels=[4000], colors='k')

        m = ax.pcolormesh(
            xbins, ybins, H,
            cmap='PiYG_r', #'RdBu_r',
            vmin=0, vmax=1,
            shading='auto'
        )

        if cyc == 'AE':
            # ax.set_title(cyc, fontweight='bold', color='r')
            ax.text(100, 900, cyc, fontweight='bold', color='r', fontsize=14)
        else:
            # ax.set_title(cyc, fontweight='bold', color='b')
            ax.text(100, 900, cyc, fontweight='bold', color='b', fontsize=14)
        ax.axis('equal')
        ax.set_xlabel('x (km)')
        ax.set_facecolor('lightgrey')

    axs[0].set_ylabel('y (km)')

    cb = fig.colorbar(m, ax=axs, location='top', shrink=0.85)
    cb.set_label(clabel, fontsize=13)

    return fig, axs

def tilt_t(df_data, grid, add_field='PV_grad_mag', field_label='PV grad.',
           figsize=(10,20), nlargest=10, width_ratios=[3, 1]):

    eddy_ids = df_data.groupby('Eddy').Age.max().nlargest(nlargest).index
    
    fig, axs = plt.subplots(
        nlargest, 2,
        figsize=figsize,
        gridspec_kw={'width_ratios': width_ratios}
    )
    
    for e, eddy in enumerate(eddy_ids):
        df = df_data[df_data.Eddy == eddy].copy()
        t = df.Day - df.Day.iloc[0]
    
        # Tilt + PV gradient
        ax = axs[e, 0]
        ax.plot(t, df.TiltDis)
        ax.set_ylabel('Tilt (km)')
    
        ax2 = ax.twinx()
        ax2.plot(t, df[add_field], 'r')
        ax2.set_ylabel(field_label, color='r')
        ax2.tick_params(axis='y', labelcolor='r')
    
        # Eddy trajectory
        ax = axs[e, 1]
        ax.contourf(grid.X_grid, grid.Y_grid, np.where(grid.mask_rho, grid.h, np.nan), cmap='Grays_r')
        ax.plot(df.xc, df.yc, color='magenta')
        ax.axis('equal')
        ax.set_title(f'{df.iloc[0].Cyc}{eddy} — {df.Age.max():.0f} days')
    
    axs[-1, 0].set_xlabel('Eddy age (days)')
    axs[-1, 1].set_xlabel('x')
    
    plt.tight_layout()
    return fig, axs
