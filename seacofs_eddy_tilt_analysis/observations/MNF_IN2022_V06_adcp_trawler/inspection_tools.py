"""Full-record ADCP inspection and OceanCurrent overlays for IN2022_V06.

esp_depth_tools.load_observations prepares finite, unique samples for ESP fitting
and adds local coordinates; it cannot preserve this full multi-instrument export.
Its fit/diagnostic plot functions do not implement geographic GSLA overlays.
Task-specific loading and plotting therefore live here; no ESP fitting or local
coordinate calculations are duplicated.
"""
from pathlib import Path
from io import StringIO
import shutil
import tempfile
from urllib.request import urlopen
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr


def load_full_record(path):
    """Read all exported rows, including literal escaped CRLF line breaks."""
    text = Path(path).read_text()
    text = text.replace(r'\r\n', '\n')
    data = pd.read_csv(StringIO(text)).rename(columns={
        'Start (UTC)': 'time', 'Mean Lat': 'lat', 'Mean Lon': 'lon',
        'Depth': 'depth', 'U (m/sec)': 'u', 'V (m/sec)': 'v',
        'Bin Attendance %': 'attendance',
    })
    data['time'] = pd.to_datetime(data['time'], utc=True, errors='coerce')
    for column in ['lat', 'lon', 'depth', 'u', 'v', 'attendance']:
        data[column] = pd.to_numeric(data[column], errors='coerce')
    data['instrument'] = data['Datafile'].str.extract(r'(os\d+\w*)\.nc$', expand=False).fillna(data['Datafile'])
    data = data.sort_values(['time', 'instrument', 'depth']).reset_index(drop=True)
    return data


def plot_record(frame, title, min_attendance=None, velocity_limit=1.5):
    track = (frame.dropna(subset=['time', 'lon', 'lat'])
             .drop_duplicates(['instrument', 'time']))
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for instrument, group in track.groupby('instrument'):
        ax.scatter(group.lon, group.lat, s=4, label=instrument, alpha=0.6)
    ax.set(xlabel='Longitude (°E)', ylabel='Latitude (°N)', title=title + ' — ship positions')
    if len(track):
        ax.set_aspect(1 / np.cos(np.deg2rad(track.lat.mean())))
        ax.legend()
    ax.grid(alpha=0.25)
    plt.show()

    instruments = sorted(frame.instrument.dropna().unique())
    if not instruments:
        print('No data in this window; adjust WINDOWS.')
        return
    fig, axes = plt.subplots(len(instruments), 2, figsize=(14, 3 * len(instruments)),
                             squeeze=False, constrained_layout=True)
    for row, instrument in enumerate(instruments):
        part = frame.loc[frame.instrument == instrument]
        if min_attendance is not None:
            part = part.loc[part.attendance >= min_attendance]
        for col, component in enumerate(['u', 'v']):
            ax = axes[row, col]
            valid = part.dropna(subset=['time', 'depth', component])
            points = ax.scatter(valid.time, valid.depth, c=valid[component], s=3,
                                cmap='RdBu_r', vmin=-velocity_limit, vmax=velocity_limit,
                                rasterized=True)
            ax.invert_yaxis()
            ax.set(title=f'{instrument}: {component.upper()}', xlabel='Time (UTC)', ylabel='Depth (m)')
            locator = mdates.AutoDateLocator()
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
            fig.colorbar(points, ax=ax, label=f'{component.upper()} (m/s)', extend='both')
    fig.suptitle(title)
    plt.show()


def load_oceancurrent(paths):
    """Load small downloaded daily subsets without requiring dask or open handles."""
    paths = list(paths)
    if not paths:
        raise FileNotFoundError('No OceanCurrent NetCDF files match GSLA_FILES.')
    parts = []
    for path in paths:
        with xr.open_dataset(path) as ds:
            parts.append(ds[['GSLA', 'UCUR', 'VCUR']].load())
    result = xr.concat(parts, dim='TIME').sortby('TIME')
    if pd.Index(result.TIME.values).has_duplicates:
        raise ValueError('Duplicate OceanCurrent map times; select one file per time.')
    return result


def select_oceancurrent(ocean, map_time, *, method='linear', tolerance_hours=12,
                        max_gap_hours=48):
    """Return a map and provenance label at a UTC time, without extrapolation.

    Linear interpolation requires bracketing maps no more than max_gap_hours
    apart. Exact timestamps use the original map. Nearest selection is bounded
    by tolerance_hours. All three fields use the same temporal weights.
    """
    if method not in ('linear', 'nearest'):
        raise ValueError("method must be 'linear' or 'nearest'.")
    if tolerance_hours < 0 or max_gap_hours <= 0:
        raise ValueError('Tolerance must be nonnegative and maximum gap positive.')
    requested = pd.Timestamp(map_time)
    requested = requested.tz_localize('UTC') if requested.tzinfo is None else requested.tz_convert('UTC')
    if pd.isna(requested):
        raise ValueError('Provide a valid map time.')
    ocean = ocean.sortby('TIME')
    times = pd.DatetimeIndex(ocean.TIME.values).tz_localize('UTC')
    if not len(times) or times.hasnans or times.has_duplicates:
        raise ValueError('OceanCurrent times must be nonempty, valid and unique.')
    index = abs(times - requested).argmin()
    actual = times[index]
    if actual == requested or method == 'nearest':
        if abs(actual - requested) > pd.Timedelta(hours=tolerance_hours):
            raise ValueError(f'No GSLA map within {tolerance_hours} h of {requested}. '
                             f'Nearest downloaded map: {actual}. Download the required day.')
        return ocean.isel(TIME=index), f'map {actual:%Y-%m-%d %H:%M} UTC'
    right = times.searchsorted(requested)
    if right == 0 or right == len(times):
        raise ValueError(f'Cannot interpolate {requested}: available maps span '
                         f'{times[0]} to {times[-1]}. Download maps bracketing the requested time.')
    before, after = times[right - 1], times[right]
    if after - before > pd.Timedelta(hours=max_gap_hours):
        raise ValueError(f'Maps bracketing {requested} are {after - before} apart; '
                         f'max_gap_hours={max_gap_hours}. Download the missing daily maps.')
    weight = float((requested - before) / (after - before))
    field = ((1 - weight) * ocean.isel(TIME=right - 1, drop=True)
             + weight * ocean.isel(TIME=right, drop=True))
    field = field.assign_coords(TIME=requested.tz_localize(None).to_datetime64())
    return field, (f'linear interpolation {requested:%Y-%m-%d %H:%M} UTC '
                   f'from {before:%Y-%m-%d %H:%M} / {after:%Y-%m-%d %H:%M} UTC')


def download_oceancurrent(map_folder, time, *, method='linear', timeout=60):
    """Cache the exact daily map or the two days needed for interpolation.

    Source: AODN THREDDS IMOS/OceanCurrent/GSLA/DM. Existing files are reused.
    Downloads are validated before being atomically moved into the cache.
    """
    requested = pd.Timestamp(time)
    requested = requested.tz_localize('UTC') if requested.tzinfo is None else requested.tz_convert('UTC')
    if pd.isna(requested) or method not in ('linear', 'nearest'):
        raise ValueError("Provide a valid UTC time and method='linear' or 'nearest'.")
    day = requested.floor('D')
    days = [day]
    if method == 'linear' and requested != day:
        days.append(day + pd.Timedelta(days=1))
    elif method == 'nearest' and requested - day > pd.Timedelta(hours=12):
        days = [day + pd.Timedelta(days=1)]
    folder = Path(map_folder)
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for date in days:
        name = f'IMOS_OceanCurrent_HV_{date:%Y%m%d}T000000Z_GSLA_FV02_DM03.nc'
        path = folder / name
        if not path.is_file():
            url = ('https://thredds.aodn.org.au/thredds/fileServer/'
                   f'IMOS/OceanCurrent/GSLA/DM/{date:%Y}/{name}')
            print(f'Downloading OceanCurrent {date:%Y-%m-%d} ...')
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=folder, suffix='.nc', delete=False) as target:
                    temporary = Path(target.name)
                    with urlopen(url, timeout=timeout) as response:
                        shutil.copyfileobj(response, target)
                check = load_oceancurrent([temporary])
                if check.sizes['TIME'] != 1 or pd.Timestamp(check.TIME.values[0]) != date.tz_localize(None):
                    raise ValueError('Downloaded file does not contain the expected daily map.')
                temporary.replace(path)
            except Exception as exc:
                raise RuntimeError(f'Could not download OceanCurrent {date:%Y-%m-%d} from {url}: {exc}') from exc
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        paths.append(path)
    return paths


def inspect_time(data, map_folder, time, *, window_hours=24, method='linear',
                 max_gap_hours=48, tolerance_hours=12, quiv_width=None, quiv_scale=None, **plot_options):
    """Automatically download/cache daily maps and plot a centred ADCP window.

    Naive times are interpreted as UTC. Returns every ADCP row in the window.
    Only the one or two maps needed at the requested time are loaded.
    """
    if not np.isfinite(window_hours) or window_hours <= 0:
        raise ValueError('window_hours must be finite and positive.')
    paths = download_oceancurrent(map_folder, time, method=method)
    ocean = load_oceancurrent(paths)
    return plot_oceancurrent_overlay(
        ocean, data, time, window_hours=window_hours, method=method,
        max_gap_hours=max_gap_hours, tolerance_hours=tolerance_hours, 
        quiv_width=quiv_width, quiv_scale=quiv_scale, **plot_options)


def plot_oceancurrent_overlay(ocean, data, map_time, *, window_hours=24,
                              tolerance_hours=12, bounds=(152, 157, -30, -25),
                              vector_stride=2, sla_limit=0.5, method='nearest',
                              max_gap_hours=48, quiv_width=None, quiv_scale=None):
    """Plot GSLA/total geostrophic velocity and return all ADCP rows in the window.

    Selection is [requested time - half width, requested time + half width).
    Bounds crop the map only; the returned observations retain the full time window.
    """
    if window_hours <= 0 or tolerance_hours < 0 or sla_limit <= 0:
        raise ValueError('Window and SLA limit must be positive; tolerance must be nonnegative.')
    if not isinstance(vector_stride, int) or vector_stride < 1:
        raise ValueError('vector_stride must be a positive integer.')
    west, east, south, north = bounds
    if not (west < east and -90 < south < north < 90):
        raise ValueError('Invalid (west, east, south, north) bounds.')
    requested = pd.Timestamp(map_time)
    requested = requested.tz_localize('UTC') if requested.tzinfo is None else requested.tz_convert('UTC')
    field, map_label = select_oceancurrent(
        ocean, requested, method=method, tolerance_hours=tolerance_hours,
        max_gap_hours=max_gap_hours)
    field = field.sortby('LONGITUDE').sortby('LATITUDE')
    field = field.sel(LONGITUDE=slice(west, east), LATITUDE=slice(south, north))
    if field.sizes['LONGITUDE'] < 2 or field.sizes['LATITUDE'] < 2:
        raise ValueError('Map bounds contain fewer than two grid points per axis.')
    start = requested - pd.Timedelta(hours=window_hours / 2)
    end = requested + pd.Timedelta(hours=window_hours / 2)
    selected = data.loc[(data.time >= start) & (data.time < end)].copy()
    track = selected.dropna(subset=['lon', 'lat', 'time']).drop_duplicates(['instrument', 'time'])
    full = data.dropna(subset=['lon', 'lat']).drop_duplicates(['instrument', 'time'])

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    lon, lat = np.meshgrid(field.LONGITUDE.values, field.LATITUDE.values)
    sla = field.GSLA.transpose('LATITUDE', 'LONGITUDE').values
    colour = ax.pcolormesh(lon, lat, sla, cmap='RdBu_r', vmin=-sla_limit,
                          vmax=sla_limit, shading='auto')
    fig.colorbar(colour, ax=ax, label='SLA (m)', extend='both', pad=0.02)
    if np.isfinite(sla).any() and np.nanmax(sla) > np.nanmin(sla):
        contours = ax.contour(lon, lat, sla, levels=np.linspace(-sla_limit, sla_limit, 11),
                              colors='0.3', linewidths=0.5)
        ax.clabel(contours, fontsize=7, fmt='%.1f')
    step = (slice(None, None, vector_stride), slice(None, None, vector_stride))
    u = field.UCUR.transpose('LATITUDE', 'LONGITUDE').values
    v = field.VCUR.transpose('LATITUDE', 'LONGITUDE').values
    # Convert east/north components to map-coordinate directions. Scale chosen
    # so 1 m/s spans ~0.35 latitude degrees; key uses the same conversion.
    cosine = np.cos(np.deg2rad(lat))
    arrows = ax.quiver(lon[step], lat[step], (u / cosine)[step], v[step],
                       angles='xy', scale_units='xy', scale=1 / 0.35,
                       color='black', width=0.003)
    ax.quiverkey(arrows, 0.80, 0.975, 1 / np.cos(np.deg2rad((south + north) / 2)),
                 '1 m/s eastward', labelpos='W', coordinates='axes')
    ax.scatter(full.lon, full.lat, s=3, color='0.65', alpha=0.5, label='Full ADCP track')
    if len(track):
        norm = plt.Normalize(mdates.date2num(start), mdates.date2num(end))
        points = ax.scatter(track.lon, track.lat, c=mdates.date2num(track.time),
                            cmap='plasma', norm=norm, s=14, #edgecolors='white',
                            linewidths=0.2, label='Selected ADCP positions', zorder=4)
        ax.quiver(track.lon, track.lat, track.u, track.v, zorder=4, scale=quiv_scale, width=quiv_width)
        bar = fig.colorbar(points, ax=ax, pad=0.10, label='ADCP time (UTC)')
        bar.formatter = mdates.DateFormatter('%d Jul %H:%M')
        bar.update_ticks()
    else:
        ax.text(0.02, 0.02, 'No ADCP data in the selected time window', transform=ax.transAxes)
    ax.set(xlim=(west, east), ylim=(south, north), xlabel='Longitude (°E)',
           ylabel='Latitude (°N)', title=(f'SLA + total geostrophic velocity: {map_label}\n'
           f'Requested: {requested:%Y-%m-%d %H:%M} UTC; ADCP {start:%d Jul %H:%M}–{end:%d Jul %H:%M} UTC'))
    ax.set_aspect(1 / np.cos(np.deg2rad((south + north) / 2)))
    ax.legend(loc='lower right')
    ax.grid(alpha=0.2)
    plt.show()
    print(f'{len(selected):,} ADCP rows in [{start}, {end}); {map_label}')
    return ax, selected
