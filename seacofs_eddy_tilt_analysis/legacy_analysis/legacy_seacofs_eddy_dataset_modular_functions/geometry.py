from __future__ import annotations

import numpy as np


def tangential_velocity(xp, yp, up, vp, xc: float, yc: float, Q, det1: bool = False):
    Q = np.asarray(Q, float)
    if Q.shape == (3,):
        q11, q12, q22 = Q
        Q = np.array([[q11, q12], [q12, q22]], float)
    if det1:
        det = np.linalg.det(Q)
        if det != 0:
            Q /= np.sqrt(det)

    xp, yp, up, vp = (np.asarray(a, float) for a in (xp, yp, up, vp))
    r = np.stack((xp - xc, yp - yc), axis=-1)
    gradient = 2.0 * (r @ Q.T)
    rotation = np.array([[0.0, -1.0], [1.0, 0.0]])
    tangent = gradient @ rotation.T
    norm = np.linalg.norm(tangent, axis=-1, keepdims=True)
    t_hat = np.divide(tangent, norm, out=np.zeros_like(tangent), where=norm > 0)
    velocity = np.stack((up, vp), axis=-1)
    vt = np.sum(velocity * t_hat, axis=-1)
    return np.where(norm.squeeze() > 0, vt, np.nan)


def find_directional_radii(u, v, X, Y, xc: float, yc: float, Q, return_index: bool = False):
    distance = np.hypot(X - xc, Y - yc)
    distance_f = np.where(np.isfinite(distance), distance, np.inf)
    nic, njc = np.unravel_index(np.argmin(distance_f), distance.shape)

    def walk(di: int, dj: int):
        prev = np.nan
        for step in range(1, max(X.shape) + 1):
            i = nic + di * step
            j = njc + dj * step
            if i < 0 or j < 0 or i >= X.shape[0] or j >= X.shape[1]:
                return np.nan
            vt = tangential_velocity(X[i, j], Y[i, j], u[i, j], v[i, j], xc, yc, Q)
            if not np.isfinite(vt):
                return np.nan
            if step > 1 and abs(vt) < abs(prev):
                if return_index:
                    return step - 1
                ip = nic + di * (step - 1)
                jp = njc + dj * (step - 1)
                return float(np.hypot(X[ip, jp] - xc, Y[ip, jp] - yc))
            prev = vt
        return np.nan

    return {
        "up": walk(0, 1),
        "right": walk(1, 0),
        "down": walk(0, -1),
        "left": walk(-1, 0),
    }
