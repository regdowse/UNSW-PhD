"""Surface topographic-PV-aligned composites using the breakdown day weighting.

x = positive signed surface topographic PV gradient; y = 90 degrees CCW
from that direction. This is a proper rotation, preserving polarity.
"""
import numpy as np
import pandas as pd
import planetary_composite_tools as pct


def reference_frames(sample, grid_angle):
    """Audit surface frame validity; retain signed PV direction without flipping.

    Bearings are geographic clockwise from north; grid_angle is radians.
    Magnitude validates the direction but never scales displacement or velocity.
    """
    if sample[['Eddy','Day']].duplicated().any():
        raise ValueError('Duplicate surface Eddy-Day keys')
    if not np.isfinite(grid_angle):
        raise ValueError('Need a finite grid angle in radians')
    out=sample.copy()
    mag=out.PV_grad_topo_mag.to_numpy(float)
    bearing=out.PV_grad_topo_theta.to_numpy(float)
    valid=np.isfinite(mag)&(mag>0)&np.isfinite(bearing)
    out['frame_valid']=valid
    out['frame_status']=np.where(valid,'usable','nonfinite bearing or nonpositive/nonfinite gradient magnitude')
    # Angle of the gradient measured anticlockwise from model-grid +x.
    out['frame_angle_rad']=np.where(valid,np.pi/2-np.deg2rad(bearing)-grid_angle,np.nan)
    # Nonlinear signed footprint means need not point exactly along mean grad(h).
    out['gradient_dot_downslope_cos']=np.nan
    if {'dhdx','dhdy'}.issubset(out):
        norm=np.hypot(out.dhdx,out.dhdy)
        dot=out.dhdx*np.sin(np.deg2rad(bearing))+out.dhdy*np.cos(np.deg2rad(bearing))
        out['gradient_dot_downslope_cos']=np.where(valid & np.isfinite(norm) & (norm>0),dot/norm,np.nan)
    return out


def rotate_xy(x,y,angle):
    """Model-grid components to along-gradient / CCW-perpendicular components."""
    c,s=np.cos(angle),np.sin(angle)
    return c*x+s*y,-s*x+c*y


def rotate_members(members, frames):
    """Rotate already surface-referenced centres in each day's own fixed frame."""
    keys=frames[['Eddy','Day','frame_angle_rad']]
    out=members.merge(keys,on=['Eddy','Day'],validate='many_to_one')
    if len(out)!=len(members) or not np.isfinite(out.frame_angle_rad).all():
        raise ValueError('Every member must have a valid surface reference frame')
    out['along_km'],out['perp_km']=rotate_xy(out.xc,out.yc,out.frame_angle_rad)
    return add_direction(out)


def add_direction(members):
    out=members.copy()
    out['distance_km']=np.hypot(out.along_km,out.perp_km)
    # Mathematical angle in the composite plane: 0 = +gradient, 90 = +perp.
    out['relative_angle_deg']=np.where(out.distance_km>1e-10,
        np.degrees(np.arctan2(out.perp_km,out.along_km))%360,np.nan)
    return out


def summarise(members,n_boot=500,seed=731):
    """Use the same cluster bootstrap estimator as the planetary comparison.

    Internal adapter uses the existing two-component implementation; returned
    columns are explicitly frame-relative, never geographic east/north.
    """
    adapter=members.copy()
    adapter['east_km']=adapter.along_km
    adapter['north_km']=adapter.perp_km
    stats,boot=pct.summarise(adapter,n_boot,seed)
    stats=stats.rename(columns={c:c.replace('east','along').replace('north','perp') for c in stats})
    stats=stats.drop(columns=['bearing_deg'])
    stats['relative_angle_deg']=np.where(stats.distance_km>1e-10,
        np.degrees(np.arctan2(stats.mean_perp,stats.mean_along))%360,np.nan)
    stats=stats.rename(columns={'cov_xy':'cov_along_perp'})
    boot=boot.rename(columns={'east_km':'along_km','north_km':'perp_km'}).drop(columns=['bearing_deg'])
    boot['relative_angle_deg']=np.where(boot.distance_km>1e-10,
        np.degrees(np.arctan2(boot.perp_km,boot.along_km))%360,np.nan)
    return stats,boot


def composite_days(sample,vertical,X,Y,*,esp):
    """Rotate each day's sampling grid AND velocities before averaging.

    X,Y are along-gradient / perpendicular km, with meshgrid [y,x] layout.
    Inverse-rotate grid points into the day's original model frame, evaluate
    its original ESP ellipse there, then rotate vector components back. This
    avoids interpolating velocity fields or averaging ellipse parameters.
    """
    if sample[['Eddy','Day']].duplicated().any():
        raise ValueError('Duplicate selected Eddy-Day keys')
    if 'frame_angle_rad' not in sample or not np.isfinite(sample.frame_angle_rad).all():
        raise ValueError('Need a finite surface frame for every selected day')
    if X.shape!=Y.shape or X.ndim!=2:
        raise ValueError('Need matching 2D X,Y arrays')
    keys=sample[['Eddy','Day','frame_angle_rad']]
    valid=pct.valid_profiles(vertical).merge(keys,on=['Eddy','Day'],validate='many_to_one')
    refs=valid.loc[valid.Depth.eq(0),['Eddy','Day']]
    usable=valid.merge(refs,on=['Eddy','Day'],validate='many_to_one')
    if usable.empty:
        raise ValueError('No selected days have a valid surface fit')
    depths=np.sort(usable.Depth.unique())
    shape=(*X.shape,len(depths))
    ut,vt=np.zeros(shape),np.zeros(shape)
    counts=np.zeros(shape,dtype=np.int64)
    depth_days=np.zeros(len(depths),dtype=np.int64)
    records=[];included=[]
    for n,((eddy,day),profile) in enumerate(usable.groupby(['Eddy','Day'],sort=True),1):
        profile=pct.recenter(profile.sort_values('Depth'),0.)
        angle=float(profile.frame_angle_rad.iloc[0])
        # Inverse rotation maps common composite coordinates to model coords.
        xm,ym=rotate_xy(X,Y,-angle)
        um,vm=pct.reconstruct(profile,xm,ym,esp)
        u,v=rotate_xy(um,vm,angle)
        contributed=False
        along,perp=rotate_xy(profile.xc.to_numpy(),profile.yc.to_numpy(),angle)
        for j,z in enumerate(profile.Depth):
            k=np.searchsorted(depths,z)
            finite=np.isfinite(u[...,j])&np.isfinite(v[...,j])
            if not finite.any():continue
            ut[...,k]+=np.where(finite,u[...,j],0.)
            vt[...,k]+=np.where(finite,v[...,j],0.)
            counts[...,k]+=finite;depth_days[k]+=1
            records.append(dict(Eddy=eddy,Day=day,Depth=float(z),
                xc=float(profile.xc.iloc[j]),yc=float(profile.yc.iloc[j]),
                along_km=along[j],perp_km=perp[j],frame_angle_rad=angle))
            contributed=True
        if contributed:included.append((eddy,day))
        if n%100==0:print(f'Reconstructed {n} eddy-days in their surface PV frames',flush=True)
    if not records:raise ValueError('No finite reconstructed velocity contributions')
    u=np.divide(ut,counts,out=np.full(shape,np.nan),where=counts>0)
    v=np.divide(vt,counts,out=np.full(shape,np.nan),where=counts>0)
    members=add_direction(pd.DataFrame(records))
    support=pd.DataFrame(dict(Depth=depths,eddy_days=depth_days,
        min_cell_days=counts.min(axis=(0,1)),max_cell_days=counts.max(axis=(0,1))))
    included=pd.DataFrame(included,columns=['Eddy','Day']).assign(included=True)
    audit=sample[['Eddy','Day']].merge(included,on=['Eddy','Day'],how='left',validate='one_to_one')
    audit['included']=audit.included.eq(True)
    reference_index=pd.MultiIndex.from_frame(refs)
    has_ref=pd.MultiIndex.from_frame(audit[['Eddy','Day']]).isin(reference_index)
    audit['status']=np.where(audit.included,'included',np.where(has_ref,'no finite ESP velocity','no valid surface fit'))
    return u,v,depths,counts,support,audit,members
