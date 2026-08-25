from __future__ import annotations

import numpy as np
import pandas as pd


def transect_indexer(ic: int, jc: int, X, Y, r: float = 30.0):
    x = np.asarray(X[:, 0])
    y = np.asarray(Y[0, :])
    x0 = x[ic]
    y0 = y[jc]
    i0 = np.searchsorted(x, x0 - r, side="left")
    i1 = np.searchsorted(x, x0 + r, side="right")
    j0 = np.searchsorted(y, y0 - r, side="left")
    j1 = np.searchsorted(y, y0 + r, side="right")
    ii = np.arange(i0, i1)
    jj = np.arange(j0, j1)
    return x[ii], np.full(ii.size, y0), np.full(jj.size, x0), y[jj], ii, jj


def nearest_ij(xc: float, yc: float, X, Y) -> tuple[int, int]:
    x = np.asarray(X[:, 0])
    y = np.asarray(Y[0, :])
    ic = np.abs(x - xc).argmin()
    jc = np.abs(y - yc).argmin()
    return int(ic), int(jc)


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


def local_eddy_arrays(X, Y, u2d, v2d, xc: float, yc: float, Q, rho_max: float):
    local = (X >= xc - rho_max) & (X <= xc + rho_max) & (Y >= yc - rho_max) & (Y <= yc + rho_max)
    Xloc = X[local]
    Yloc = Y[local]
    uloc = u2d[local]
    vloc = v2d[local]
    dx = Xloc - xc
    dy = Yloc - yc
    rho2 = Q[0, 0] * dx**2 + 2 * Q[0, 1] * dx * dy + Q[1, 1] * dy**2
    rho = np.sqrt(np.where(rho2 >= 0, rho2, np.nan))
    return Xloc, Yloc, uloc, vloc, rho


def bad_doppio_row(row, fnumber: int) -> dict:
    return {
        "Day": row.Day,
        "fnumber": fnumber,
        "nxc": getattr(row, "nxc", np.nan),
        "nyc": getattr(row, "nyc", np.nan),
        "nCyc": getattr(row, "Cyc", np.nan),
        "xc": np.nan,
        "yc": np.nan,
        "w": np.nan,
        "Q": np.nan,
        "Omega0": np.nan,
        "Omega": np.nan,
        "Rc": np.nan,
        "psi0": np.nan,
        "R": np.nan,
    }


def q_columns_to_matrix(row) -> np.ndarray:
    return np.array([[row.q11, row.q12], [row.q12, row.q22]], dtype=float)


def matrix_to_q_columns(df: pd.DataFrame, column: str = "Q") -> pd.DataFrame:
    out = df.copy()
    q = np.stack(out[column].to_numpy())
    out["q11"] = q[:, 0, 0]
    out["q12"] = q[:, 0, 1]
    out["q22"] = q[:, 1, 1]
    return out.drop(columns=[column])

