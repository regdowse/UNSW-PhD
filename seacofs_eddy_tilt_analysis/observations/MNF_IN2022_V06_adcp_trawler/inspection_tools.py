"""Full-record ADCP inspection and OceanCurrent overlays for IN2022_V06.

esp_depth_tools.load_observations prepares finite, unique samples for ESP fitting
and adds local coordinates; it cannot preserve this full multi-instrument export.
Its fit/diagnostic plot functions do not implement geographic GSLA overlays.
Task-specific loading and plotting therefore live here; no ESP fitting or local
coordinate calculations are duplicated.
"""
from pathlib import Path
from io import StringIO
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


def plot_oceancurrent_overlay(ocean, data, map_time, *, window_hours=24,
                              tolerance_hours=12, bounds=(152, 157, -30, -25),
                              vector_stride=2, sla_limit=0.5):
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
    times = pd.DatetimeIndex(ocean.TIME.values).tz_localize('UTC')
    index = abs(times - requested).argmin()
    actual = times[index]
    if abs(actual - requested) > pd.Timedelta(hours=tolerance_hours):
        raise ValueError(f'No GSLA map within {tolerance_hours} h of {requested}. '
                         f'Nearest downloaded map: {actual}. Download the required day.')
    field = ocean.isel(TIME=index).sortby('LONGITUDE').sortby('LATITUDE')
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
    ax.quiverkey(arrows, 0.80, 1.025, 1 / np.cos(np.deg2rad((south + north) / 2)),
                 '1 m/s eastward', labelpos='W', coordinates='axes')
    ax.scatter(full.lon, full.lat, s=3, color='0.65', alpha=0.5, label='Full ADCP track')
    if len(track):
        norm = plt.Normalize(mdates.date2num(start), mdates.date2num(end))
        points = ax.scatter(track.lon, track.lat, c=mdates.date2num(track.time),
                            cmap='plasma', norm=norm, s=14, edgecolors='white',
                            linewidths=0.2, label='Selected ADCP positions', zorder=4)
        bar = fig.colorbar(points, ax=ax, pad=0.10, label='ADCP time (UTC)')
        bar.formatter = mdates.DateFormatter('%d Jul %H:%M')
        bar.update_ticks()
    else:
        ax.text(0.02, 0.02, 'No ADCP data in the selected time window', transform=ax.transAxes)
    ax.set(xlim=(west, east), ylim=(south, north), xlabel='Longitude (°E)',
           ylabel='Latitude (°N)', title=(f'SLA + total geostrophic velocity: {actual:%Y-%m-%d %H:%M} UTC\n'
           f'Requested: {requested:%Y-%m-%d %H:%M} UTC; ADCP {start:%d Jul %H:%M}–{end:%d Jul %H:%M} UTC'))
    ax.set_aspect(1 / np.cos(np.deg2rad((south + north) / 2)))
    ax.legend(loc='lower right')
    ax.grid(alpha=0.2)
    plt.show()
    print(f'{len(selected):,} ADCP rows in [{start}, {end}); map: {actual}')
    return selected
