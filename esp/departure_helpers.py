"""Noise-free departures for notebook 02; SI units, explicit synthetic truth."""
import numpy as np
import matplotlib.pyplot as plt
import validation_helpers as vh


def vortex_velocity(xy, case, tail_power=np.inf):
    """Algebraic envelope with the SAME centre gradient and peak radius as Gaussian.

    s=rho²/Rc²; g_p(s)=(1+2s/(2p-1))^-p, p>1.
    g_inf=exp(-s). u=Omega J Q r g(s) is nondivergent. For every p,
    g(0)=1 and normalised speed peaks at rho=Rc/sqrt(2).
    The peak SPEED is allowed to differ. Rc is the generating length scale;
    only the Gaussian limit has a Gaussian characteristic radius.
    """
    if not (tail_power > 1):
        raise ValueError('tail_power must exceed 1 or be infinity')
    r=np.asarray(xy)-case['centre']; qr=r@case['q']
    s=np.sum(r*qr,axis=1)/case['rc']**2
    g=np.exp(-s) if np.isinf(tail_power) else (1+2*s/(2*tail_power-1))**(-tail_power)
    return case['omega']*np.c_[-qr[:,1],qr[:,0]]*g[:,None]


def background_velocity(xy, case, kind, strength, angle_deg):
    """Signed strengths: current/Vref or shear*Rc/Vref. Angle CCW from east.

    Shear in axes rotated by angle: u'=S y', v'=0. Trace=0, curl=-S.
    Its origin is the fixed synthetic centre. No background is removed/fitted.
    """
    a=np.deg2rad(angle_deg)
    rotation=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
    if kind=='current':
        return np.tile(strength*case['speed_ref']*rotation[:,0],(len(xy),1))
    if kind=='shear':
        r=(np.asarray(xy)-case['centre'])@rotation
        shear=strength*case['speed_ref']/case['rc']
        return np.c_[shear*r[:,1],np.zeros(len(r))]@rotation.T
    if kind=='profile':
        return np.zeros((len(xy),2))
    raise ValueError(kind)


def evaluation_points(case, extent=2., n=45):
    """Regular cell-centre lattice within rho<=extent Rc, excluding fit locations."""
    width=extent*case['rc']*np.sqrt(case['alpha'])
    axis=(np.arange(n)+.5)*2*width/n-width
    x,y=np.meshgrid(axis,axis)
    xy=np.c_[x.ravel(),y.ravel()]+case['centre']
    r=xy-case['centre']
    xy=xy[np.sum((r@case['q'])*r,axis=1)<=(extent*case['rc'])**2]
    # Do not score training points, even when two regular lattices happen to align.
    distance2=np.sum((xy[:,None,:]-case['xy'][None,:,:])**2,axis=2)
    return xy[np.min(distance2,axis=1)>1e-6]


def run_case(backend, case, *, kind, strength=0., angle_deg=0., tail_power=np.inf,
             radius=30000., rc_limit=100000., solo_local=True,
             centre_guess=10000., search_half_width=60000.):
    eddy=vortex_velocity(case['xy'],case,tail_power)
    background=background_velocity(case['xy'],case,kind,strength,angle_deg)
    fit=vh.fit_case(backend,case,eddy+background,radius,rc_limit,
        diagnostic_outer=True,solo_local=solo_local,
        centre_guess=centre_guess,search_half_width=search_half_width)
    # Do not report an error against a nonexistent 'true Gaussian Rc'.
    fit.pop('Rc_error_pct',None)
    fit.update(method='SOLO local' if case['name']=='SOLO' and solo_local else case['name'],
        test=kind,strength=strength,angle_deg_background=angle_deg,
        tail_power=tail_power,departure=0. if np.isinf(tail_power) else 1/tail_power,
        background_speed_ms=strength*case['speed_ref'] if kind=='current' else 0.,
        shear_s_inv=strength*case['speed_ref']/case['rc'] if kind=='shear' else 0.,
        peak_radius_truth_km=case['rc']/np.sqrt(2)/1000)
    if fit['valid']:
        # det(Q)=1: rho_peak is the equal-area radius of the peak-speed ellipse.
        fit['peak_radius_error_pct']=100*(fit['Rc_m']/case['rc']-1)
        q=vh.shape_matrix(fit['alpha'],fit['angle_deg'])
        xy=evaluation_points(case)
        prediction=vh.velocity(xy,[fit['xc_m'],fit['yc_m']],fit['Omega'],fit['Rc_m'],q)
        truth_eddy=vortex_velocity(xy,case,tail_power)
        truth_total=truth_eddy+background_velocity(xy,case,kind,strength,angle_deg)
        for label,truth in [('eddy',truth_eddy),('total',truth_total)]:
            rmse=np.sqrt(np.mean(np.sum((prediction-truth)**2,axis=1)))
            fit[f'{label}_velocity_rmse_ms']=rmse
            fit[f'{label}_velocity_error_pct']=100*rmse/case['speed_ref']
        fit['n_evaluation']=len(xy)
    return fit


def plot_test(results):
    """One small figure per test; no ensemble bands for deterministic experiments."""
    profile=results.test.iloc[0]=='profile'
    orientations=[None] if profile else sorted(results.angle_deg_background.unique())
    metrics=([('centre_error_km','Centre error (km)'),
              ('peak_radius_error_pct','Peak-radius error (%)'),
              ('eddy_velocity_error_pct','Eddy velocity RMSE (% Vref)')] if profile else
             [('centre_error_km','Centre error (km)'),
              ('vorticity_error_pct','Eddy core-vorticity error (%)'),
              ('eddy_velocity_error_pct','Eddy velocity RMSE (% Vref)')])
    fig,axes=plt.subplots(len(orientations),3,figsize=(12,3.3*len(orientations)),
                           squeeze=False,constrained_layout=True)
    colours={'SOLO local':'tab:purple','SOLO':'tab:blue','DOPPIO':'tab:orange','LATTE':'tab:green'}
    for row,orientation in enumerate(orientations):
        subset=results if profile else results[results.angle_deg_background==orientation]
        for method,group in subset.groupby('method'):
            group=group.sort_values('departure' if profile else 'strength')
            x=group['departure' if profile else 'strength']
            for ax,(metric,title) in zip(axes[row],metrics):
                # Invalid attempts remain at their x values as gaps, not connected across.
                y=group[metric].where(group.valid) if metric in group else np.full(len(group),np.nan)
                ax.plot(x,y,'o-',ms=4,color=colours[method],label=method)
                ax.set_title(title); ax.grid(alpha=.2)
                ax.set_xlabel('Profile departure 1/p (0 = Gaussian)' if profile else
                    ('Signed current / Vref' if results.test.iloc[0]=='current' else 'Signed shear Rc / Vref'))
        if not profile:
            axes[row,0].set_ylabel(f'Orientation {orientation:g}°')
    axes[0,-1].legend(fontsize=8)
    failures=int((~results.valid).sum())
    fig.suptitle(f"{results.test.iloc[0].capitalize()}: {failures}/{len(results)} fits failed diagnostics (see table); noise-free",fontsize=11)
    return fig
