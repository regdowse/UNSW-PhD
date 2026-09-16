"""Audited, eddy-equal composites of existing depth-dependent ESP fits.

ESP's model_uv_at_xy uses the velocity convention of the existing SEACOFS
reconstruction notebooks (east/north output on model-grid x/y locations).
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

KEYS = ['Eddy', 'Day']
PARAMS = ['xc', 'yc', 'Rc', 'Omega', 'q11', 'q12', 'q22']


def unique(frame, keys):
    if frame[keys].isna().any().any() or frame.duplicated(keys).any():
        raise ValueError(f'Missing or duplicate keys: {keys}')


def select_depths(vertical, targets=(0, 200, 500), max_distance=60):
    """Resolve requested depths to actual levels; never interpolate."""
    available = np.sort(vertical.Depth.dropna().unique())
    if not len(available):
        raise ValueError('No fitted depths')
    chosen = np.array([available[np.argmin(abs(available-t))] for t in targets])
    if np.any(abs(chosen-np.asarray(targets)) > max_distance):
        raise ValueError(f'Requested {targets}, nearest cached levels {chosen}; choose explicitly')
    if len(np.unique(chosen)) != len(chosen) or np.any(np.diff(chosen) <= 0):
        raise ValueError('Choose distinct increasing depths')
    return chosen


def audit_population(surface, vertical, pv, depths, dominance=2):
    """Basic fit/coverage QC; preserves all source rows in the audit."""
    unique(surface, KEYS); unique(vertical, KEYS+['Depth']); unique(pv, KEYS)
    depths = np.asarray(depths, float)
    if len(depths) < 2 or np.any(np.diff(depths) <= 0) or depths[0] < 0:
        raise ValueError('Need at least two increasing nonnegative exact depths')
    if dominance <= 1:
        raise ValueError('dominance must exceed 1')
    p = vertical.loc[vertical.Depth.isin(depths)].copy()
    finite = np.isfinite(p[PARAMS].to_numpy(float)).all(axis=1)
    p['valid_fit'] = (finite & (p.Rc > 0) & (p.q11 > 0) & (p.q22 > 0)
                      & (p.q11*p.q22-p.q12**2 > 0) & p.Omega.ne(0))
    p['fit_sign'] = np.sign(p.Omega)
    stats = p.groupby(KEYS).agg(levels=('Depth','nunique'),
                                valid_levels=('valid_fit','sum'),
                                signs=('fit_sign','nunique'))
    cols = ['topo_plan_ratio'] + [c for c in ['Ro','h','PV_grad_coherence'] if c in pv]
    a = surface.drop(columns=cols, errors='ignore').merge(pv[KEYS+cols], on=KEYS,
                                                         how='left', validate='one_to_one')
    a = a.merge(stats, on=KEYS, how='left', validate='one_to_one')
    surface_cols = ['xc','yc','Rc','q11','q12','q22']
    a['invalid_surface'] = ~(np.isfinite(a[surface_cols].to_numpy(float)).all(axis=1)
                              & a.Rc.gt(0) & a.q11.gt(0) & a.q22.gt(0)
                              & (a.q11*a.q22-a.q12**2).gt(0))
    a['missing_depth'] = a.levels.fillna(0).ne(len(depths))
    a['invalid_fit'] = a.valid_levels.fillna(0).ne(len(depths))
    a['polarity_inconsistent'] = a.signs.fillna(0).ne(1)
    # Southern Hemisphere convention: positive fitted rotation AE, negative CE.
    sign = p.groupby(KEYS).fit_sign.first().rename('fit_sign')
    a = a.merge(sign, on=KEYS, how='left', validate='one_to_one')
    a['polarity_inconsistent'] |= ~((a.Cyc.eq('AE') & a.fit_sign.gt(0)) |
                                    (a.Cyc.eq('CE') & a.fit_sign.lt(0)))
    ratio = pd.to_numeric(a.topo_plan_ratio, errors='coerce')
    a['missing_regime'] = ~np.isfinite(ratio)
    a['regime'] = np.select([ratio <= -np.log(dominance), ratio >= np.log(dominance)],
                            ['Planetary','Topographic'], default='Mixed')
    a.loc[a.missing_regime, 'regime'] = 'Unknown'
    reasons = ['invalid_surface','missing_depth','invalid_fit','polarity_inconsistent','missing_regime']
    a['eligible'] = ~a[reasons].any(axis=1)
    a['exclusion_reason'] = a[reasons].apply(lambda r: '|'.join(r.index[r]), axis=1)
    a['surface_fit_provenance'] = 'unknown: processed table can contain interpolated parameters'
    a['reconstruction_skill_status'] = 'not measured: geometry QC only'
    return a, p


def slope_basis(row, grid, gx, gy):
    """Return model-grid onshore basis plus strength/coherence/support."""
    import seacofs_tilt_tools as tilt
    ii, jj = tilt.core_grid_indices(row, grid)
    if not len(ii):
        return np.full(2,np.nan), 0., 0., 0.
    x, y = gx[ii,jj], gy[ii,jj]
    ok = np.isfinite(x) & np.isfinite(y)
    fraction = float(ok.mean())
    if not ok.any():
        return np.full(2,np.nan), 0., 0., fraction
    mean = np.array([x[ok].mean(), y[ok].mean()])
    strength = np.linalg.norm(mean)
    local = np.hypot(x[ok],y[ok]).mean()
    coherence = strength/local if local > 0 else 0.
    return (-mean/strength if strength > 0 else np.full(2,np.nan)), strength, coherence, fraction


def add_frames(audit, grid, *, min_slope=1e-4, min_coherence=.5, min_fraction=.8):
    """One rigid slope frame per day; weak slope days retain geographic eligibility."""
    import seacofs_tilt_tools as tilt
    gx,gy = tilt.phys_grad(grid.h, grid.X_grid*1000, grid.Y_grid*1000, grid.mask_rho)
    a = audit.copy()
    for c in ['axis_x','axis_y','slope_strength','slope_coherence','slope_valid_fraction']:
        a[c] = np.nan
    for i,row in a.loc[a.eligible & a.regime.eq('Topographic')].iterrows():
        axis,s,c,f = slope_basis(row,grid,gx,gy)
        a.loc[i,['axis_x','axis_y','slope_strength','slope_coherence','slope_valid_fraction']] = [*axis,s,c,f]
    a['slope_eligible'] = ((a.slope_strength >= min_slope) & (a.slope_coherence >= min_coherence)
                           & (a.slope_valid_fraction >= min_fraction))
    a['selected'] = a.eligible & (a.regime.eq('Planetary') |
                                   (a.regime.eq('Topographic') & a.slope_eligible))
    bad = a.eligible & a.regime.eq('Topographic') & ~a.slope_eligible
    a.loc[bad,'exclusion_reason'] = 'unreliable_slope_frame'
    a.loc[a.eligible & a.regime.eq('Mixed'),'exclusion_reason'] = 'mixed_regime_not_in_initial_four'
    a['group'] = a.Cyc.astype(str)+'_'+a.regime
    return a


def rotation_basis(row, angle, frame):
    """Rows of B are composite axes in model coordinates; G axes in east/north."""
    R = np.array([[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]])
    if frame == 'geographic':
        B = R.copy()
    elif frame == 'onshore':
        a = np.array([row.axis_x,row.axis_y],float)
        if not np.isfinite(a).all() or not np.isclose(np.linalg.norm(a),1):
            raise ValueError('Invalid onshore axis')
        B = np.array([a,[-a[1],a[0]]])
    else:
        raise ValueError('Unknown frame')
    return B, B @ R.T


def reconstruct_day(row, profile, depths, coordinate, angle, model_uv, *, frame,
                    units='Rc', wet_mask=None):
    """Evaluate one entire column about one shallow reference without vertical interpolation."""
    unique(profile, ['Depth'])
    p = profile.set_index('Depth').loc[depths]
    ref = p.iloc[0]
    if units not in ('Rc','km'):
        raise ValueError('units must be Rc or km')
    scale = float(ref.Rc) if units == 'Rc' else 1.
    B,G = rotation_basis(row,angle,frame)
    xx,yy = np.meshgrid(coordinate,coordinate,indexing='ij')
    x = ref.xc+scale*(xx*B[0,0]+yy*B[1,0])
    y = ref.yc+scale*(xx*B[0,1]+yy*B[1,1])
    fields=[]; centres=[]
    for depth,level in p.iterrows():
        q=np.array([[level.q11,level.q12],[level.q12,level.q22]])
        u,v=model_uv(x*1000,y*1000,level.xc*1000,level.yc*1000,q,level.Omega,level.Rc*1000)
        field=np.stack([G[0,0]*u+G[0,1]*v,G[1,0]*u+G[1,1]*v])
        if wet_mask is not None:
            field[:,~wet_mask(x,y,depth)] = np.nan
        if not np.isfinite(field).any():
            raise ValueError(f'No valid velocity at depth {depth}')
        fields.append(field)
        centres.append(B @ np.array([level.xc-ref.xc,level.yc-ref.yc])/scale)
    return np.stack(fields),np.stack(centres)


def ocean_mask(grid):
    """Nearest native-cell support; outside grid/land/below seabed are missing."""
    from scipy.interpolate import RegularGridInterpolator
    coords=(grid.x_grid,grid.y_grid)
    h=RegularGridInterpolator(coords,grid.h,method='nearest',bounds_error=False,fill_value=np.nan)
    wet=RegularGridInterpolator(coords,grid.mask_rho,method='nearest',bounds_error=False,fill_value=0)
    def mask(x,y,z):
        points=np.stack([x,y],axis=-1)
        return (wet(points) > .5) & (h(points) > z)
    return mask


def build_members(selected, profiles, depths, coordinate, grid, model_uv, output, *,
                  units='Rc', mask_ocean=True):
    """Stream days into per-eddy/group means, with persisted counts and centrelines."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    unique(selected,KEYS)
    profile_lookup={key:p for key,p in profiles.groupby(KEYS,sort=False)}
    masker=ocean_mask(grid) if mask_ocean else None
    records=[]
    for (group,eddy),days in selected.groupby(['group','Eddy'],sort=True):
        frame='onshore' if group.endswith('Topographic') else 'geographic'
        sums=counts=centre_sum=None
        for row in days.itertuples(index=False):
            field,centre=reconstruct_day(row,profile_lookup[(row.Eddy,row.Day)],depths,coordinate,
                                        grid.angle,model_uv,frame=frame,units=units,wet_mask=masker)
            valid=np.isfinite(field)
            if sums is None:
                sums=np.zeros_like(field);counts=np.zeros_like(field,dtype=np.int32)
                centre_sum=np.zeros_like(centre)
            sums+=np.where(valid,field,0);counts+=valid;centre_sum+=centre
        mean=np.divide(sums,counts,out=np.full_like(sums,np.nan),where=counts>0)
        file=f'{group}_eddy_{int(eddy)}.npz'
        np.savez_compressed(output/file,mean=mean,day_counts=counts,
                            centres=centre_sum/len(days),days=days.Day.to_numpy(),
                            depths=depths,coordinate=coordinate)
        records.append(dict(group=group,Eddy=eddy,days=len(days),frame=frame,file=file))
        if len(records) % 25 == 0:
            print(f'Reconstructed {len(records)} eddy/group members', flush=True)
    manifest=pd.DataFrame(records,columns=['group','Eddy','days','frame','file'])
    manifest.to_csv(output/'members.csv',index=False)
    return manifest


def bootstrap_mean(values, draws, *, min_members=20, batch=256):
    """Pointwise percentile intervals over cluster draws, bounded by pixel batches."""
    values=np.asarray(values,float)
    flat=values.reshape(len(values),-1)
    shape=values.shape[1:]
    support=np.isfinite(flat).sum(axis=0)
    mean=np.full(flat.shape[1],np.nan);lo=mean.copy();hi=mean.copy()
    for start in range(0,flat.shape[1],batch):
        v=flat[:,start:start+batch];finite=np.isfinite(v)
        total=np.where(finite,v,0)
        denom=finite.sum(axis=0)
        m=np.divide(total.sum(axis=0),denom,out=np.full(v.shape[1],np.nan),where=denom>0)
        # draws contains cluster multiplicities, not day multiplicities.
        den=np.dot(draws.astype(float), finite.astype(float))
        boot=np.divide(np.dot(draws.astype(float), total),den,out=np.full_like(den,np.nan),where=den>0)
        good=support[start:start+batch]>=min_members
        mean[start:start+batch]=np.where(good,m,np.nan)
        for j in np.flatnonzero(good):
            samples=boot[:,j];samples=samples[np.isfinite(samples)]
            if len(samples):
                lo[start+j],hi[start+j]=np.quantile(samples,[.025,.975])
    return dict(mean=mean.reshape(shape),low=lo.reshape(shape),high=hi.reshape(shape),
                support=support.reshape(shape))


def summarise_members(manifest, member_dir, *, n_boot=500, seed=731, min_members=20):
    """Global cluster resampling retains covariance of shared eddies across groups."""
    if n_boot < 100 or min_members < 2:
        raise ValueError('Use >=100 bootstrap draws and >=2 members')
    ids=np.sort(manifest.Eddy.unique())
    if not len(ids):
        raise ValueError('No selected eddies')
    rng=np.random.default_rng(seed)
    draws=rng.multinomial(len(ids),np.full(len(ids),1/len(ids)),size=n_boot)
    results={}
    for group,part in manifest.groupby('group',sort=True):
        fields=[];centres=[]
        for f in part.file:
            with np.load(Path(member_dir)/f) as member:
                fields.append(member['mean']);centres.append(member['centres'])
        weights=draws[:,np.searchsorted(ids,part.Eddy.to_numpy())]
        results[group]={'field':bootstrap_mean(np.stack(fields),weights,min_members=min_members),
                        'centre':bootstrap_mean(np.stack(centres),weights,min_members=min_members),
                        'members':np.stack(centres),'eddies':len(part),'days':int(part.days.sum()),
                        'frame':part.frame.iloc[0]}
    return results


def plot_summary(results, depths, coordinate, units='Rc'):
    """Member spread, pointwise CI, mean velocity and depth-dependent support."""
    import matplotlib.pyplot as plt
    figures=[]
    for group,r in results.items():
        fig,axes=plt.subplots(1,4,figsize=(16,5),constrained_layout=True)
        names=('Onshore','Alongshore') if r['frame']=='onshore' else ('East','North')
        for j in range(2):
            for member in r['members']:
                axes[j].plot(member[:,j],depths,color='.6',alpha=.12,lw=.6)
            c=r['centre'];axes[j].plot(c['mean'][:,j],depths,color='k',lw=2)
            axes[j].fill_betweenx(depths,c['low'][:,j],c['high'][:,j],alpha=.3)
            axes[j].axvline(0,color='.4',ls=':');axes[j].invert_yaxis()
            axes[j].set(xlabel=f'{names[j]} centre displacement ({units})',ylabel='Depth (m)')
        k=len(depths)-1
        f=r['field'];u,v=f['mean'][k]
        xx,yy=np.meshgrid(coordinate,coordinate,indexing='ij')
        im=axes[2].pcolormesh(xx,yy,np.hypot(u,v),shading='auto',cmap='viridis')
        sl=(slice(None,None,3),slice(None,None,3))
        axes[2].quiver(xx[sl],yy[sl],u[sl],v[sl],color='white')
        axes[2].set(title=f'Mean velocity at {depths[k]:.1f} m',xlabel=f'{names[0]} ({units})',ylabel=f'{names[1]} ({units})',aspect='equal')
        fig.colorbar(im,ax=axes[2],label='Speed of mean velocity (m/s)')
        im=axes[3].pcolormesh(xx,yy,f['support'][k,0],shading='auto')
        axes[3].set(title='Contributing eddies',aspect='equal')
        fig.colorbar(im,ax=axes[3])
        fig.suptitle(f'{group}: {r["eddies"]} eddies, {r["days"]} days; shading = pointwise 95% eddy-bootstrap CI')
        figures.append((group,fig))
    return figures


def fingerprint(config, paths):
    """Record parameters and source metadata (not a hash of multi-GB source contents)."""
    payload={'config':config,'sources':{str(p):{'size':Path(p).stat().st_size,
             'mtime_ns':Path(p).stat().st_mtime_ns} for p in paths},
             'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text=json.dumps(payload,sort_keys=True,indent=2)
    return hashlib.sha256(text.encode()).hexdigest()[:16],payload


def plot_sections_3d(result, depths, coordinate, units='Rc'):
    """Two central vertical sections and a decimated 3-D horizontal-vector view."""
    import matplotlib.pyplot as plt
    f=result['field']; centre=result['centre']['mean']
    mid=int(np.argmin(abs(coordinate)))
    names=('Onshore','Alongshore') if result['frame']=='onshore' else ('East','North')
    fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
    sections=[f['mean'][:,1,:,mid],f['mean'][:,0,mid,:]]
    finite=np.concatenate([s[np.isfinite(s)] for s in sections])
    limit=np.max(abs(finite)) if len(finite) else 1.
    limit=max(limit,1e-12)
    for j,(ax,s) in enumerate(zip(axes,sections)):
        im=ax.pcolormesh(coordinate,depths,s,shading='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
        ax.plot(centre[:,j],depths,'k.-')
        ax.invert_yaxis()
        ax.set(xlabel=f'{names[j]} ({units})',ylabel='Depth (m)',
               title=f'{names[1-j]} velocity on central section')
        fig.colorbar(im,ax=ax,label='m/s')
    fig3=plt.figure(figsize=(9,7),constrained_layout=True)
    ax=fig3.add_subplot(111,projection='3d')
    stride=max(1,len(coordinate)//9)
    c=coordinate[::stride]
    xx,yy=np.meshgrid(c,c,indexing='ij')
    fields=f['mean'][:,:,::stride,::stride]
    maximum=np.nanmax(abs(fields)) if np.isfinite(fields).any() else 1.
    arrow_scale=(coordinate[-1]-coordinate[0])*.08/max(maximum,1e-12)
    for k,z in enumerate(depths):
        u,v=fields[k]
        ok=np.isfinite(u)&np.isfinite(v)
        ax.quiver(xx[ok],yy[ok],np.full(ok.sum(),z),u[ok]*arrow_scale,
                  v[ok]*arrow_scale,np.zeros(ok.sum()),color='tab:blue',alpha=.55)
    ax.plot(centre[:,0],centre[:,1],depths,'ko-',label='Mean fitted centreline')
    ax.set(xlabel=f'{names[0]} ({units})',ylabel=f'{names[1]} ({units})',zlabel='Depth (m)',
           xlim=(coordinate[0],coordinate[-1]),ylim=(coordinate[0],coordinate[-1]),
           zlim=(depths[-1],depths[0]))
    ax.legend()
    return [('sections',fig),('3d_horizontal_velocity',fig3)]
