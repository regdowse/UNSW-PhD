from __future__ import annotations

import numpy as np


def nencioli(u, v, X, Y, a: int, b: int, flip_sign: bool = False):
    """Nencioli eddy detection on a Cartesian grid with X/Y in km."""
    borders = max(a, b) + 1
    vel = np.sqrt(u**2 + v**2)
    bound = vel.shape

    eddy_uv = np.empty((0, 2))
    eddy_c = np.empty((0, 2))
    eddy = np.empty((0, 3))

    for i in range(borders, bound[0] - borders + 1):
        wrk = v[i, :]
        sign_change = np.diff(np.sign(wrk))
        indices = np.where((sign_change != 0) & ~np.isnan(sign_change))[0]
        indices = indices[(indices >= borders) & (indices <= bound[1] - borders - 1)]

        for ii in indices:
            var = 0
            if wrk[ii] >= 0:
                if wrk[ii - a] > wrk[ii] and wrk[ii + 1 + a] < wrk[ii + 1]:
                    var = -1
            else:
                if wrk[ii - a] < wrk[ii] and wrk[ii + 1 + a] > wrk[ii + 1]:
                    var = 1

            if var == -1:
                ok1 = (
                    u[i - a, ii] <= 0
                    and u[i - a, ii] <= u[i - 1, ii]
                    and u[i + a, ii] >= 0
                    and u[i + a, ii] >= u[i + 1, ii]
                )
                ok2 = (
                    u[i - a, ii + 1] <= 0
                    and u[i - a, ii + 1] <= u[i - 1, ii + 1]
                    and u[i + a, ii + 1] >= 0
                    and u[i + a, ii + 1] >= u[i + 1, ii + 1]
                )
                if ok1 or ok2:
                    eddy_uv = np.vstack((eddy_uv, [X[i, ii], Y[i, ii]], [X[i, ii + 1], Y[i, ii + 1]]))
                else:
                    var = 0
            elif var == 1:
                ok1 = (
                    u[i - a, ii] >= 0
                    and u[i - a, ii] >= u[i - 1, ii]
                    and u[i + a, ii] <= 0
                    and u[i + a, ii] <= u[i + 1, ii]
                )
                ok2 = (
                    u[i - a, ii + 1] >= 0
                    and u[i - a, ii + 1] >= u[i - 1, ii + 1]
                    and u[i + a, ii + 1] <= 0
                    and u[i + a, ii + 1] <= u[i + 1, ii + 1]
                )
                if ok1 or ok2:
                    eddy_uv = np.vstack((eddy_uv, [X[i, ii], Y[i, ii]], [X[i, ii + 1], Y[i, ii + 1]]))
                else:
                    var = 0

            if var != 0:
                srch = vel[i - b : i + b + 1, ii - b : ii + b + 2]
                sx = X[i - b : i + b + 1, ii - b : ii + b + 2]
                sy = Y[i - b : i + b + 1, ii - b : ii + b + 2]
                xind, yind = np.unravel_index(np.nanargmin(srch), srch.shape)
                imin = i - b + xind
                jmin = ii - b + yind
                srch2 = vel[
                    max(imin - b, 0) : min(imin + b + 1, bound[0]),
                    max(jmin - b, 0) : min(jmin + b + 1, bound[1]),
                ]
                if np.nanmin(srch2) == np.nanmin(srch):
                    eddy_c = np.vstack((eddy_c, [sx[xind, yind], sy[xind, yind]]))
                else:
                    var = 0

            if var != 0:
                d = a - 1
                u_small = u[max(imin - d, 0) : min(imin + d + 1, bound[0]), max(jmin - d, 0) : min(jmin + d + 1, bound[1])]
                v_small = v[max(imin - d, 0) : min(imin + d + 1, bound[0]), max(jmin - d, 0) : min(jmin + d + 1, bound[1])]
                if not np.isnan(u_small).any() and not np.isnan(v_small).any():
                    u_bound = np.concatenate((u_small[0, :], u_small[1:, -1], u_small[-1, -2::-1], u_small[-2:0:-1, 0]))
                    v_bound = np.concatenate((v_small[0, :], v_small[1:, -1], v_small[-1, -2::-1], v_small[-2:0:-1, 0]))
                    quadrants = np.zeros_like(u_bound, dtype=int)
                    quadrants[(u_bound >= 0) & (v_bound >= 0)] = 1
                    quadrants[(u_bound < 0) & (v_bound >= 0)] = 2
                    quadrants[(u_bound < 0) & (v_bound < 0)] = 3
                    quadrants[(u_bound >= 0) & (v_bound < 0)] = 4
                    spin = np.where(quadrants == 4)[0]
                    if spin.size > 0 and spin.size != quadrants.size:
                        if spin[0] == 0:
                            spin = np.where(quadrants != 4)[0]
                            spin = np.array([spin[0] - 1])
                        quadrants[spin[-1] + 1 :] += 4
                        if not np.any(np.diff(quadrants) > 1) and not np.any(np.diff(quadrants) < 0):
                            eddy = np.vstack((eddy, [sx[xind, yind], sy[xind, yind], var]))

    eddy_uv = np.unique(eddy_uv, axis=0)
    eddy_c = np.unique(eddy_c, axis=0)
    eddy = np.unique(eddy, axis=0)
    if flip_sign and eddy.size:
        eddy[:, 2] *= -1
    return eddy_uv, eddy_c, eddy

