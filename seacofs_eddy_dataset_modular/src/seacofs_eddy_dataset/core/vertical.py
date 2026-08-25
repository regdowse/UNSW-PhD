from __future__ import annotations

import numpy as np
import pandas as pd


def interp_3d_to_reference_depths(var3d, z3d, target_depths) -> np.ndarray:
    target_depths = np.asarray(target_depths)
    out = np.full(var3d.shape[:2] + (target_depths.size,), np.nan)
    for i in range(var3d.shape[0]):
        for j in range(var3d.shape[1]):
            z = np.abs(np.asarray(z3d[i, j, :], dtype=float))
            values = np.asarray(var3d[i, j, :], dtype=float)
            order = np.argsort(z)
            out[i, j, :] = np.interp(
                target_depths,
                z[order],
                values[order],
                left=np.nan,
                right=np.nan,
            )
    return out


def surface_profile_row(row) -> dict:
    return {
        "z": 0,
        "Depth": 0,
        "xc": row.xc,
        "yc": row.yc,
        "ic": row.ic,
        "jc": row.jc,
        "w": row.w,
        "Q": np.array([[row.q11, row.q12], [row.q12, row.q22]]),
        "Omega0": np.nan,
        "Omega": row.Omega,
        "Rc": row.Rc,
        "psi0": row.psi0,
        "R": row.R,
    }


def profile_dict_to_frame(profile: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(profile)


def completed_eddy_days(profile_dict: dict) -> pd.DataFrame:
    rows = []
    for eddy_key, day_dict in profile_dict.items():
        eddy = int(str(eddy_key).replace("Eddy", ""))
        for day_key in day_dict:
            day = int(str(day_key).replace("Day", ""))
            rows.append({"Eddy": eddy, "Day": day})
    return pd.DataFrame(rows, columns=["Eddy", "Day"])
