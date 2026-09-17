"""Planetary composite analysis, using the breakdown reconstruction conventions."""
import numpy as np
import pandas as pd
import composite_breakdown_tools as cbt

def valid_profiles(vertical):
    needed = ['Eddy', 'Day', 'Depth', 'xc', 'yc', 'Rc', 'Omega', 'q11', 'q12', 'q22']
    valid = vertical.loc[np.isfinite(vertical[needed].to_numpy(float)).all(axis=1)].copy()
    valid = valid.loc[(valid.Depth >= 0) & (valid.Rc > 0)
                      & (valid.q11 > 0) & (valid.q22 > 0)
                      & (valid.q11 * valid.q22 - valid.q12**2 > 0)]
    if valid.duplicated(['Eddy', 'Day', 'Depth']).any():
        raise ValueError('Duplicate Eddy-Day-Depth rows; resolve these before averaging.')
    return valid


def choose_case(surface, vertical, eddy_id=None, day=None, reference_depth=0.):
    # Match BOTH keys to the chosen sample, before choosing an example.
    keys = surface[['Eddy', 'Day']]
    if keys.duplicated().any():
        raise ValueError('Duplicate Eddy-Day rows in the selected sample.')
    valid = valid_profiles(vertical).merge(keys, on=['Eddy', 'Day'], validate='many_to_one')
    if eddy_id is not None:
        valid = valid.loc[valid.Eddy.eq(eddy_id)]
    if day is not None:
        valid = valid.loc[valid.Day.eq(day)]
    references = valid.loc[valid.Depth.eq(reference_depth), ['Eddy', 'Day']]
    if references.empty:
        raise ValueError('No selected eddy-day has a valid profile at the reference depth.')
    # Deterministic first available day (not the deepest or most tilted example).
    choice = references.sort_values(['Eddy', 'Day']).iloc[0]
    profile = valid.loc[valid.Eddy.eq(choice.Eddy) & valid.Day.eq(choice.Day)].sort_values('Depth').copy()
    row = surface.loc[surface.Eddy.eq(choice.Eddy) & surface.Day.eq(choice.Day)].iloc[0]
    return row, profile.reset_index(drop=True)


def recenter(profile, reference_depth=0.):
    profile = profile.copy()
    reference = profile.loc[profile.Depth.eq(reference_depth)]
    if len(reference) != 1:
        raise ValueError('Each day needs exactly one row at the common reference depth.')
    profile['xc'] -= float(reference.xc.iloc[0])
    profile['yc'] -= float(reference.yc.iloc[0])
    return profile


def reconstruct(profile, X, Y, esp):
    # Input positions are kilometres; ESP expects metres.
    U, V = [], []
    for row in profile.itertuples(index=False):
        Q = np.array([[row.q11, row.q12], [row.q12, row.q22]])
        u, v = esp.model_uv_at_xy(X*1e3, Y*1e3, row.xc*1e3, row.yc*1e3,
                                  Q, row.Omega, row.Rc*1e3)
        U.append(u); V.append(v)
    if not U:
        raise ValueError('No fitted levels to reconstruct.')
    return np.stack(U, axis=-1), np.stack(V, axis=-1)


def composite_days(sample, vertical, X, Y, reference_depth=0., max_depth=None, return_centres=False, *, esp):
    keys = sample[['Eddy', 'Day']].copy()
    if keys.duplicated().any():
        raise ValueError('Duplicate Eddy-Day rows in sample would double-count days.')
    valid = valid_profiles(vertical).merge(keys, on=['Eddy', 'Day'], validate='many_to_one')
    if max_depth is not None:
        valid = valid.loc[valid.Depth.le(max_depth)]
    ref_keys = valid.loc[valid.Depth.eq(reference_depth), ['Eddy', 'Day']]
    usable = valid.merge(ref_keys, on=['Eddy', 'Day'], validate='many_to_one')
    if usable.empty:
        raise ValueError('No selected days have a valid profile at the reference depth.')
    depths = np.sort(usable.Depth.unique())
    shape = (*X.shape, len(depths))  # y, x, exact fitted depth
    u_tot, v_tot = np.zeros(shape), np.zeros(shape)
    counts = np.zeros(shape, dtype=np.int64)
    depth_days = np.zeros(len(depths), dtype=np.int64)
    included = []
    centre_records = []
    for (eddy, day), profile in usable.groupby(['Eddy', 'Day'], sort=True):
        profile = recenter(profile.sort_values('Depth'), reference_depth)
        u, v = reconstruct(profile, X, Y, esp)
        contributed = False
        for j, depth in enumerate(profile.Depth):
            k = np.searchsorted(depths, depth)
            finite = np.isfinite(u[..., j]) & np.isfinite(v[..., j])
            if not finite.any():
                continue
            u_tot[..., k] += np.where(finite, u[..., j], 0.)
            v_tot[..., k] += np.where(finite, v[..., j], 0.)
            counts[..., k] += finite
            depth_days[k] += 1
            # Record precisely the fitted levels that contributed finite velocity.
            centre_records.append(dict(Eddy=eddy, Day=day, Depth=float(depth),
                                       xc=float(profile.iloc[j].xc), yc=float(profile.iloc[j].yc)))
            contributed = True
        if contributed:
            included.append((eddy, day))
        if len(included) and len(included) % 100 == 0:
            print(f'Reconstructed {len(included)} eddy-days', flush=True)
    if not included:
        raise ValueError('ESP returned no finite velocity for the selected days.')
    u_ave = np.divide(u_tot, counts, out=np.full(shape, np.nan), where=counts > 0)
    v_ave = np.divide(v_tot, counts, out=np.full(shape, np.nan), where=counts > 0)
    included = pd.DataFrame(included, columns=['Eddy', 'Day'])
    audit = keys.merge(included.assign(included=True), on=['Eddy', 'Day'], how='left', validate='one_to_one')
    audit['included'] = audit.included.eq(True)
    ref_index = pd.MultiIndex.from_frame(ref_keys)
    has_reference = pd.MultiIndex.from_frame(audit[['Eddy', 'Day']]).isin(ref_index)
    audit['status'] = np.where(audit.included, 'included',
                             np.where(has_reference, 'no finite ESP velocity',
                                      'no valid profile at reference depth'))
    support = pd.DataFrame({'Depth': depths, 'eddy_days': depth_days,
                            'min_cell_days': counts.min(axis=(0, 1)),
                            'max_cell_days': counts.max(axis=(0, 1))})
    result = (u_ave, v_ave, depths, counts, support, audit)
    if return_centres:
        return (*result, pd.DataFrame(centre_records))
    return result


def prepare_profiles(sample, vertical, split_depth=1000.):
    """Classify each eddy-day by deepest valid fit before any display-depth limit."""
    keys = sample[['Eddy', 'Day']]
    if keys.duplicated().any():
        raise ValueError('Duplicate selected Eddy-Day keys')
    valid = valid_profiles(vertical).merge(keys, on=['Eddy','Day'], validate='many_to_one')
    refs = valid.loc[valid.Depth.eq(0), ['Eddy','Day','xc','yc']].rename(columns={'xc':'surface_xc','yc':'surface_yc'})
    valid = valid.merge(refs, on=['Eddy','Day'], validate='many_to_one')
    extent = valid.groupby(['Eddy','Day']).Depth.max().rename('max_fit_depth_m').reset_index()
    extent['extent_group'] = np.where(extent.max_fit_depth_m > split_depth, 'deep', 'shallow')
    audit = sample.merge(extent, on=['Eddy','Day'], how='left', validate='one_to_one')
    audit['profile_status'] = np.where(audit.max_fit_depth_m.notna(), 'usable', 'no valid surface fit')
    valid['xc'] -= valid.surface_xc
    valid['yc'] -= valid.surface_yc
    return audit, valid


def geographic_members(members, angle):
    """Surface-to-depth displacement; angle is ROMS grid rotation in radians."""
    out = members.copy()
    out['east_km'] = np.cos(angle)*out.xc - np.sin(angle)*out.yc
    out['north_km'] = np.sin(angle)*out.xc + np.cos(angle)*out.yc
    out['distance_km'] = np.hypot(out.east_km, out.north_km)
    out['bearing_deg'] = np.where(out.distance_km > 1e-10,
        np.degrees(np.arctan2(out.east_km, out.north_km)) % 360, np.nan)
    return out


def summarise(members, n_boot=500, seed=731):
    """Day-weighted mean vector and eddy-cluster bootstrap in geographic axes."""
    geo = members[['Eddy','Day','Depth','east_km','north_km']].rename(columns={'east_km':'xc','north_km':'yc'})
    result = cbt.centre_statistics(geo, np.sort(geo.Depth.unique()), n_boot, seed)
    result = result.rename(columns={c:c.replace('xc','east').replace('yc','north') for c in result if 'xc' in c or 'yc' in c})
    result = result.drop(columns=['mean_tilt_x','mean_tilt_y'])
    result['distance_km'] = np.hypot(result.mean_east, result.mean_north)
    result['bearing_deg'] = np.where(result.distance_km > 1e-10,
        np.degrees(np.arctan2(result.mean_east,result.mean_north)) % 360, np.nan)
    ids = np.sort(members.Eddy.unique())
    draws = np.random.default_rng(seed).multinomial(len(ids), np.ones(len(ids))/len(ids), size=n_boot)
    boot_rows=[]
    for i,row in result.iterrows():
        d=members.loc[members.Depth.eq(row.Depth)]
        a=d.groupby('Eddy').agg(e=('east_km','sum'), n=('north_km','sum'), count=('Day','size')).reindex(ids,fill_value=0)
        denom=draws@a['count'].to_numpy(float)
        xy=np.divide(draws@a[['e','n']].to_numpy(),denom[:,None],out=np.full((n_boot,2),np.nan),where=denom[:,None]>0)
        dist=np.hypot(xy[:,0],xy[:,1])
        result.loc[i,'mean_member_distance_km']=d.distance_km.mean()
        result.loc[i,'coherence']=row.distance_km/d.distance_km.mean() if d.distance_km.mean()>1e-10 else np.nan
        result.loc[i,['distance_ci_low','distance_ci_high']]=np.nanquantile(dist,[.025,.975]) if d.Eddy.nunique()>1 else [np.nan,np.nan]
        if d.Eddy.nunique()>1:
            boot_rows.append(pd.DataFrame({'Depth':row.Depth,'draw':np.arange(n_boot),'east_km':xy[:,0],'north_km':xy[:,1],
                'distance_km':dist,'bearing_deg':np.where(dist>1e-10,np.degrees(np.arctan2(xy[:,0],xy[:,1]))%360,np.nan)}))
    return result, pd.concat(boot_rows,ignore_index=True) if boot_rows else pd.DataFrame(columns=['Depth','draw','east_km','north_km','distance_km','bearing_deg'])


def matched_members(members, depths):
    """Keep only days with every requested EXACT depth; no interpolation."""
    required=np.unique(depths)
    hits=members.loc[members.Depth.isin(required)].groupby(['Eddy','Day']).Depth.nunique()
    keys=hits.loc[hits.eq(len(required))].reset_index()[['Eddy','Day']]
    return members.merge(keys,on=['Eddy','Day'],validate='many_to_one')


def rose(ax, members, bins=(0,10,20,30,40,np.inf), min_distance=0., cmap='Blues'):
    """Fixed 16 sectors, percentages of retained finite nonzero vectors."""
    import matplotlib.pyplot as plt
    good=members.loc[np.isfinite(members.bearing_deg)&np.isfinite(members.distance_km)&members.distance_km.ge(min_distance)]
    width=2*np.pi/16
    if len(good):
        sector=((good.bearing_deg.to_numpy()+11.25)%360//22.5).astype(int)
        counts=np.zeros((len(bins)-1,16))
        for i,(lo,hi) in enumerate(zip(bins[:-1],bins[1:])):
            counts[i]=np.bincount(sector[(good.distance_km>=lo)&(good.distance_km<hi)],minlength=16)*100/len(good)
        bottom=np.zeros(16)
        for i,color in enumerate(plt.get_cmap(cmap)(np.linspace(.25,.9,len(counts)))):
            ax.bar(np.arange(16)*width,counts[i],width=width,bottom=bottom,color=color,edgecolor='white',linewidth=.3,label=f'{bins[i]:g}–{bins[i+1]:g}')
            bottom+=counts[i]
    else:
        ax.text(.5,.5,'No nonzero supported vectors',ha='center',transform=ax.transAxes)
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
    return len(good)
