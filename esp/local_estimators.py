"""Explicit local variants. Never modifies ESP_zonodo or its imported functions."""
import numpy as np
from scipy.optimize import curve_fit


class FitFailure(ValueError):
    """An expected, explicitly classified fitting failure."""


def solo_local_initialisation(l, vt, vn, *, centre_guess, search_half_width,
                              core_thresh=30000.):
    """Bounded local cubic crossing; same SOLO parameter equations.

    Initial fit uses only the fixed search interval. Refinement uses the chosen
    inner window. Both fit an odd cubic with crossing constrained to observed
    support inside the search interval. Coordinates are scaled for conditioning.
    No true eddy parameters are read. Returns the original SOLO tuple.
    """
    l, vt, vn = [np.asarray(a, float) for a in (l, vt, vn)]
    good = np.isfinite(l) & np.isfinite(vt) & np.isfinite(vn)
    order = np.argsort(l[good])
    l, vt, vn = [a[good][order] for a in (l, vt, vn)]
    if search_half_width <= 0 or core_thresh <= 0:
        raise ValueError('Search/window widths must be positive.')
    mask = abs(l - centre_guess) <= search_half_width
    if mask.sum() < 4:
        raise FitFailure('solo_insufficient_search')
    x, y = l[mask], vn[mask]
    crossings = np.flatnonzero((y[:-1] * y[1:] <= 0) & (y[:-1] != y[1:]))
    if not len(crossings):
        raise FitFailure('solo_no_supported_crossing')
    roots = x[crossings] - y[crossings] * np.diff(x)[crossings] / np.diff(y)[crossings]
    guess = roots[np.argmin(abs(roots - centre_guess))]
    scale = search_half_width

    def fit(mask, guess):
        if mask.sum() < 4 or np.unique(l[mask]).size < 4:
            raise FitFailure('solo_insufficient_core')
        xx = (l[mask] - centre_guess) / scale
        yy = vn[mask]
        lower = max(xx.min(), -1.)
        upper = min(xx.max(), 1.)
        if lower >= upper:
            raise FitFailure('solo_no_search_overlap')
        z0 = np.clip((guess - centre_guess) / scale, lower + 1e-10, upper - 1e-10)
        dz = xx-z0
        c, d = np.linalg.lstsq(np.c_[dz, dz**3], yy, rcond=None)[0]
        def model(z, root, c, d):
            delta = z-root
            return c*delta+d*delta**3
        try:
            p, _ = curve_fit(model, xx, yy, p0=[z0,c,d],
                bounds=([lower,-np.inf,-np.inf],[upper,np.inf,np.inf]), maxfev=10000)
        except (RuntimeError, ValueError, FloatingPointError) as exc:
            raise FitFailure('solo_optimisation_failed') from exc
        if min(p[0]-lower, upper-p[0]) <= 1e-7:
            raise FitFailure('solo_crossing_at_bound')
        return centre_guess + scale*p[0], p[1]/scale

    initial, _ = fit(mask, guess)
    core = abs(l-initial) <= core_thresh
    x0, c = fit(core, initial)
    if not np.isfinite(c) or abs(c) < np.finfo(float).eps:
        raise FitFailure('solo_zero_gradient')
    dl = (l[core]-x0)/scale
    a, _ = np.linalg.lstsq(np.c_[np.ones(core.sum()), dl**2], vt[core], rcond=None)[0]
    return x0, a/c, 2*c, np.eye(2), c


def outer_fit_diagnostic(xi, yi, ui, vi, xc, yc, q, *, Omega0, Rc0=None, Rc_max=1e5):
    """Original Gaussian objective/filter/bounds, explicit status instead of fallback.

    Matches the public ESP out_core_param_fit numerical formulation. No radius
    cap is imposed during optimisation; the original post-fit cap is classified.
    Returns attempted estimates even when that cap is exceeded.
    """
    xi, yi, ui, vi = [np.asarray(a, float) for a in (xi,yi,ui,vi)]
    dx, dy = xi-xc, yi-yc
    q11,q12,q22 = q[0,0],q[1,0],q[1,1]
    rho2 = q11*dx**2 + 2*q12*dx*dy + q22*dy**2
    a,b = q11*dx+q12*dy, q12*dx+q22*dy
    norm = np.sqrt(a*a+b*b)
    with np.errstate(divide='ignore', invalid='ignore'):
        vt = (-b*ui+a*vi)/norm
    keep = ((vt<=0) if Omega0<0 else (vt>=0)) & np.isfinite(rho2) & np.isfinite(norm) & np.isfinite(vt) & (rho2>=0) & (norm!=0)
    r2 = rho2[keep]
    rho = np.sqrt(r2)
    target = vt[keep]*rho/norm[keep]
    out = dict(status='insufficient_outer', n_outer_retained=int(keep.sum()),
               Rc_m=np.nan, Omega=np.nan, psi0=np.nan)
    if len(target)<3 or np.unique(rho).size<3:
        return out
    if Rc0 is None:
        Rc0=max(float(rho[np.argmax(abs(target))]*np.sqrt(2)),1e-6)
    out['Rc_seed_m']=Rc0
    def model(r2, omega, rc):
        return omega*np.sqrt(r2)*np.exp(-r2/rc**2)
    try:
        p,_=curve_fit(model,r2,target,p0=[Omega0,Rc0],
            bounds=([-np.inf,1e-8],[np.inf,np.inf]),maxfev=10000)
    except (RuntimeError, ValueError, FloatingPointError) as exc:
        out.update(status='outer_optimisation_failed', outer_message=str(exc))
        return out
    omega,rc=p
    out.update(Rc_m=rc,Omega=omega,psi0=-.5*omega*rc**2)
    if not np.all(np.isfinite(p)) or rc<=0:
        out['status']='nonfinite_outer'
    elif rc>Rc_max:
        out['status']='outer_radius_limit'
    else:
        out['status']='ok'
    return out
