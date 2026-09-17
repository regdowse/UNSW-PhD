"""Centre statistics and inner/outer ESP refits for the simple breakdown notebook."""
from __future__ import annotations

import numpy as np
import pandas as pd


def centre_statistics(members, depths, n_boot=500, seed=731):
    """Equal-day centre mean and sample covariance; whole-eddy bootstrap of that mean.

    Variance describes member spread, NOT error in the mean. Resampling clusters
    preserves unequal track lengths and within-eddy dependence; it does not change
    the estimator to an equal-eddy mean. No fit-parameter uncertainty is included.
    """
    required=['Eddy','Day','Depth','xc','yc']
    if members.empty or members[required].isna().any().any():
        raise ValueError('Need finite contributing centres')
    if members.duplicated(['Eddy','Day','Depth']).any():
        raise ValueError('Duplicate contributing centre keys')
    if not np.isfinite(members[['Depth','xc','yc']].to_numpy(float)).all():
        raise ValueError('Nonfinite centre coordinates')
    if n_boot<100:
        raise ValueError('Use at least 100 bootstrap draws')
    ids=np.sort(members.Eddy.unique())
    draws=np.random.default_rng(seed).multinomial(len(ids),np.full(len(ids),1/len(ids)),size=n_boot)
    rows=[]
    for depth in depths:
        d=members.loc[members.Depth.eq(depth)]
        n=len(d);m=d[['xc','yc']].mean().to_numpy(float)
        covariance=d[['xc','yc']].cov().to_numpy() if n>1 else np.full((2,2),np.nan)
        summary=d.groupby('Eddy').agg(xsum=('xc','sum'),ysum=('yc','sum'),n=('Day','size')).reindex(ids,fill_value=0)
        weights=draws.astype(float)
        denominator=np.dot(weights,summary.n.to_numpy(float))
        means=np.divide(np.dot(weights,summary[['xsum','ysum']].to_numpy(float)),denominator[:,None],
                        out=np.full((n_boot,2),np.nan),where=denominator[:,None]>0)
        finite=np.isfinite(means).all(axis=1)
        lo,hi=(np.quantile(means[finite],[.025,.975],axis=0)
               if d.Eddy.nunique()>=2 and finite.any() else (np.full(2,np.nan),np.full(2,np.nan)))
        rows.append(dict(Depth=float(depth),n_eddy_days=n,n_eddies=d.Eddy.nunique(),
                         mean_xc=m[0],mean_yc=m[1],var_xc=covariance[0,0],var_yc=covariance[1,1],
                         cov_xy=covariance[0,1],std_xc=np.sqrt(covariance[0,0]),std_yc=np.sqrt(covariance[1,1]),
                         mean_xc_ci_low=lo[0],mean_xc_ci_high=hi[0],mean_yc_ci_low=lo[1],mean_yc_ci_high=hi[1],
                         mean_tilt_x=-m[0],mean_tilt_y=-m[1],bootstrap_valid_draws=int(finite.sum())))
    return pd.DataFrame(rows)


def _skill(u,v,uh,vh,mask):
    valid=mask & np.isfinite(u) & np.isfinite(v) & np.isfinite(uh) & np.isfinite(vh)
    if not valid.any():return np.nan,np.nan,np.nan
    u,v,uh,vh=(a[valid] for a in (u,v,uh,vh))
    mse=np.mean((uh-u)**2+(vh-v)**2)
    rms=np.sqrt(np.mean(u**2+v**2))
    variance=np.mean((u-u.mean())**2+(v-v.mean())**2)
    return np.sqrt(mse),np.sqrt(mse)/rms if rms>0 else np.nan,1-mse/variance if variance>0 else np.nan


def fit_composite(X,Y,u_ave,v_ave,depths,counts,esp,*,transect_radius_km=30.,
                  rho_min_km=30.,rho_max_km=200.,out_core_fac=1.75,
                  max_jump_km=100.,min_cell_days=1,expected_sign=None,initial_center=(0.,0.)):
    """Fit composite fields using the external ESP routines and pipeline radii helper.

    Array layout is [y, x, depth] (np.meshgrid default). As in the pipeline,
    fits take km and m/s, yielding raw rotation in m/s/km; divide by 1000 for
    SI reconstruction. Output Rc,xc,yc,R are km, w/Omega0/Omega are s^-1 and
    psi0 is m^2/s. Return every requested depth, including failed fits.
    """
    from seacofs_eddy_dataset.core.doppio import find_directional_radii
    X,Y=np.asarray(X,float),np.asarray(Y,float)
    depths=np.asarray(depths,float)
    if X.shape!=Y.shape or X.ndim!=2 or u_ave.shape!=(*X.shape,len(depths)) or v_ave.shape!=u_ave.shape or counts.shape!=u_ave.shape:
        raise ValueError('Expected X,Y [y,x] and velocity/counts [y,x,depth]')
    x,y=X[0,:],Y[:,0]
    if (len(x)<4 or len(y)<4 or not np.allclose(X,x[None,:]) or not np.allclose(Y,y[:,None])
        or np.any(np.diff(x)<=0) or np.any(np.diff(y)<=0) or np.any(np.diff(depths)<=0)):
        raise ValueError('Use increasing meshgrid(x,y) coordinates and depths')
    if min_cell_days<1 or not 0<rho_min_km<=rho_max_km or transect_radius_km<=0 or out_core_fac<=0 or max_jump_km<=0:
        raise ValueError('Invalid fit controls')
    if expected_sign not in (None,-1,1):raise ValueError('expected_sign must be +1, -1 or None')
    for name in ('doppio','out_core_param_fit','model_uv_at_xy'):
        if not callable(getattr(esp,name,None)):raise AttributeError(f'ESP backend missing {name}')
    previous=np.asarray(initial_center,float)
    if previous.shape!=(2,) or not np.isfinite(previous).all():raise ValueError('Invalid initial centre')
    u_fit=np.full_like(u_ave,np.nan,dtype=float);v_fit=np.full_like(v_ave,np.nan,dtype=float)
    rows=[]
    for k,depth in enumerate(depths):
        names=['xc','yc','w','q11','q12','q22','Omega0','Omega','Rc','psi0','R','rho_limit_km',
               'vector_rmse_ms','vector_nrmse','vector_r2','domain_rmse_ms','domain_nrmse','domain_r2',
               'w_full_fit','psi0_raw_km_ms','n_outer_points','outer_min_days','n_boundary_radii']
        record=dict.fromkeys(names,np.nan)
        record.update(Depth=float(depth),status='failed',inner_ok=False,fit_ok=False,
                      seed_xc=float(previous[0]),seed_yc=float(previous[1]),reason='')
        try:
            valid=np.isfinite(u_ave[...,k]) & np.isfinite(v_ave[...,k]) & (counts[...,k]>=min_cell_days)
            u=np.where(valid,u_ave[...,k],np.nan);v=np.where(valid,v_ave[...,k],np.nan)
            iy=int(np.argmin(abs(y-previous[1])));ix=int(np.argmin(abs(x-previous[0])))
            x0,y0=x[ix],y[iy]
            # Require a complete in-domain cross, as in the dataset's vertical fit.
            if min(x0-x[0],x[-1]-x0,y0-y[0],y[-1]-y0)<transect_radius_km:
                raise ValueError('Seed too close to grid boundary for requested transects')
            xs=np.flatnonzero(abs(x-x0)<=transect_radius_km)
            ys=np.flatnonzero(abs(y-y0)<=transect_radius_km)
            if (len(xs)<4 or len(ys)<4 or np.isfinite(u[iy,xs]+v[iy,xs]).sum()<4
                or np.isfinite(u[ys,ix]+v[ys,ix]).sum()<4 or not valid[iy,ix]):
                raise ValueError('Fewer than four supported points on a transect or unsupported intersection')
            xc,yc,w,Q,omega0=esp.doppio(x[xs],np.full(len(xs),y0),u[iy,xs],v[iy,xs],
                                      np.full(len(ys),x0),y[ys],u[ys,ix],v[ys,ix])
            Q=np.asarray(Q,float)
            if not np.isfinite([xc,yc,w,omega0]).all() or Q.shape!=(2,2) or not np.isfinite(Q).all():
                raise ValueError('Nonfinite inner fit')
            if not np.allclose(Q,Q.T) or np.linalg.eigvalsh(Q).min()<=0:
                raise ValueError('Inner ellipse is not positive definite')
            if w==0 or omega0==0 or np.sign(w)!=np.sign(omega0):raise ValueError('Invalid inner rotation')
            if expected_sign is not None and np.sign(w)!=expected_sign:raise ValueError('Inner fit has wrong polarity')
            if not (x[0]<=xc<=x[-1] and y[0]<=yc<=y[-1]):raise ValueError('Fitted centre outside grid')
            if np.hypot(xc-previous[0],yc-previous[1])>max_jump_km:raise ValueError('Centre jump exceeds control')
            record.update(xc=xc,yc=yc,w=w/1000.,Omega0=omega0/1000.,q11=Q[0,0],q12=Q[0,1],q22=Q[1,1],inner_ok=True)
            radii=find_directional_radii(u,v,X,Y,xc,yc,Q)
            finite_radii=np.array(list(radii.values()),float)
            finite_radii=finite_radii[np.isfinite(finite_radii)&(finite_radii>0)]
            if not len(finite_radii):raise ValueError('No supported tangential-speed maximum found')
            radius=float(finite_radii.mean())
            rho_limit=float(np.clip(radius*out_core_fac,rho_min_km,rho_max_km))
            record.update(R=radius,rho_limit_km=rho_limit,n_boundary_radii=len(finite_radii))
            dx,dy=X-xc,Y-yc
            rho2=Q[0,0]*dx**2+2*Q[0,1]*dx*dy+Q[1,1]*dy**2
            outer=valid & (rho2<=rho_limit**2)
            record['n_outer_points']=int(outer.sum())
            if outer.sum()<10:raise ValueError('Fewer than 10 supported outer-fit points')
            record['outer_min_days']=int(counts[...,k][outer].min())
            rc,psi0,omega=esp.out_core_param_fit(X[outer],Y[outer],u[outer],v[outer],xc,yc,Q,w)
            if not np.isfinite([rc,psi0,omega]).all() or rc<=0 or omega==0:
                raise ValueError('Invalid outer ESP fit')
            if np.sign(omega)!=np.sign(w):raise ValueError('Outer fit reverses inner polarity')
            uh,vh=esp.model_uv_at_xy(X*1000.,Y*1000.,xc*1000.,yc*1000.,Q,omega/1000.,rc*1000.)
            if not np.isfinite(uh[outer]).all() or not np.isfinite(vh[outer]).all():
                raise ValueError('Refit reconstruction nonfinite inside fitting region')
            rmse,nrmse,r2=_skill(u,v,uh,vh,outer)
            drmse,dnrmse,dr2=_skill(u,v,uh,vh,valid)
            record.update(Rc=float(rc),psi0=float(psi0)*1000.,psi0_raw_km_ms=float(psi0),Omega=float(omega)/1000.,
                          w_full_fit=float(omega)*np.trace(Q)/1000.,vector_rmse_ms=rmse,vector_nrmse=nrmse,vector_r2=r2,
                          domain_rmse_ms=drmse,domain_nrmse=dnrmse,domain_r2=dr2,status='ok',fit_ok=True)
            u_fit[...,k]=np.where(valid,uh,np.nan);v_fit[...,k]=np.where(valid,vh,np.nan)
            previous=np.array([xc,yc])
        except Exception as exc:
            # Report numerical/backend failures per depth rather than silently
            # retaining a previous fit or interpolating through missing fits.
            record['reason']=f'{type(exc).__name__}: {exc}'
        rows.append(record)
    fits=pd.DataFrame(rows)
    # Keep absolute composite centres. Provide an additional surface-referenced
    # displacement for comparing tilt if the reference-depth fit succeeded.
    fits['xc_from_fit_reference']=np.nan;fits['yc_from_fit_reference']=np.nan
    if fits.fit_ok.iloc[0]:
        ok=fits.fit_ok
        fits.loc[ok,'xc_from_fit_reference']=fits.loc[ok,'xc']-fits.xc.iloc[0]
        fits.loc[ok,'yc_from_fit_reference']=fits.loc[ok,'yc']-fits.yc.iloc[0]
    return fits,u_fit,v_fit


def plot_centrelines(stats,fits=None):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,9),constrained_layout=True)
    for j,c in enumerate(('xc','yc')):
        mean=stats[f'mean_{c}'].to_numpy();sd=stats[f'std_{c}'].to_numpy();z=stats.Depth.to_numpy()
        for i in range(2):
            ax=axes[i,j];ax.plot(mean,z,'k.-',label='Mean constituent centre')
            if i==0:
                ax.fill_betweenx(z,mean-sd,mean+sd,alpha=.2,label='±1 SD across eddy-days')
            else:
                ax.fill_betweenx(z,stats[f'mean_{c}_ci_low'],stats[f'mean_{c}_ci_high'],alpha=.3,label='95% eddy-bootstrap CI of mean')
            if fits is not None:
                ax.plot(fits[c].where(fits.fit_ok),fits.Depth,'r.--',label='Centre of full composite fit')
            ax.axvline(0,color='.5',ls=':');ax.set_ylim(z[-1],z[0])
            ax.set(xlabel=f'Relative model-grid {c[0]} (km)',ylabel='Depth (m)')
            ax.legend(fontsize=8)
    fig.suptitle('Constituent mean and fitted composite centrelines\nSpread of eddy-days (top) and uncertainty of the mean (bottom)')
    return fig,axes
