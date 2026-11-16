"""Local fitted-profile QC and equal-day, whole-eddy bootstrap statistics."""
import numpy as np
import pandas as pd

def valid_profiles(vertical):
    needed = ['Eddy', 'Day', 'Depth', 'xc', 'yc', 'Rc', 'Omega', 'q11', 'q12', 'q22']
    valid = vertical.loc[np.isfinite(vertical[needed].to_numpy(float)).all(axis=1)].copy()
    valid = valid.loc[(valid.Depth >= 0) & (valid.Rc > 0)
                      & (valid.q11 > 0) & (valid.q22 > 0)
                      & (valid.q11 * valid.q22 - valid.q12**2 > 0)]
    if valid.duplicated(['Eddy', 'Day', 'Depth']).any():
        raise ValueError('Duplicate Eddy-Day-Depth rows; resolve these before averaging.')
    return valid

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
