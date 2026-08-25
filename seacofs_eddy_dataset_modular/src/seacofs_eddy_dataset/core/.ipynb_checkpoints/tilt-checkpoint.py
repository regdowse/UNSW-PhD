from __future__ import annotations

import numpy as np
import pandas as pd


def add_tilt_displacement(profile: pd.DataFrame) -> pd.DataFrame:
    out = profile.sort_values("Depth").copy()
    if out.empty:
        return out
    x0 = out.iloc[0]["xc"]
    y0 = out.iloc[0]["yc"]
    out["tilt_dx"] = out["xc"] - x0
    out["tilt_dy"] = out["yc"] - y0
    out["tilt_distance"] = np.hypot(out["tilt_dx"], out["tilt_dy"])
    out["tilt_angle_deg"] = np.degrees(np.arctan2(out["tilt_dy"], out["tilt_dx"]))
    return out


def summarise_tilt(profile: pd.DataFrame) -> dict:
    profile = add_tilt_displacement(profile)
    valid = profile.dropna(subset=["Depth", "tilt_distance"])
    if len(valid) < 2:
        return {"tilt_slope": np.nan, "max_tilt_distance": np.nan, "max_depth": np.nan}
    depth = valid["Depth"].to_numpy(dtype=float)
    distance = valid["tilt_distance"].to_numpy(dtype=float)
    slope = np.polyfit(depth, distance, deg=1)[0]
    return {
        "tilt_slope": slope,
        "max_tilt_distance": np.nanmax(distance),
        "max_depth": np.nanmax(depth),
    }


def bearing(a, b) -> float:
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    return float((np.degrees(np.arctan2(dx, dy)) + 360) % 360)


def compute_weighted_tilt(
    profiles: pd.DataFrame,
    eddy: int,
    num: int = 6,
    depth_int: int = 10,
    max_depth: int = 1000,
    min_depth_range: int = 200,
    var: str = "Depth",
    eps: float = 1e-10,
    min_points: int = 5,
    bearing_offset: float = 20.0,
) -> pd.DataFrame:
    rows = []
    df_eddy = profiles.loc[profiles["Eddy"].eq(eddy)].copy()
    if df_eddy.empty:
        return pd.DataFrame(columns=["Eddy", "Day", "TiltDis", "TiltDir"])

    day_nums = np.arange(int(df_eddy["Day"].min()), int(df_eddy["Day"].max()) + 1)
    target_depths = np.arange(0, max_depth + depth_int, depth_int)
    full_idx = target_depths[:-1]
    diffs_xc: dict[str, pd.Series] = {}
    diffs_yc: dict[str, pd.Series] = {}

    for d, day in enumerate(day_nums):
        key = f"t_{d}"
        df_day = df_eddy.loc[df_eddy["Day"].eq(day)].copy()
        if len(df_day) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            continue

        if var == "Depth":
            df_day[var] = np.abs(df_day[var])
        df_day = df_day.loc[df_day[var] <= max_depth].set_index(var).sort_index()
        depths = df_day.index.values
        valid_depths = target_depths[(target_depths >= depths.min()) & (target_depths <= depths.max())]
        if len(valid_depths) < 2:
            diffs_xc[key] = pd.Series(np.nan, index=full_idx)
            diffs_yc[key] = pd.Series(np.nan, index=full_idx)
            continue

        xc_i = np.interp(valid_depths, depths, df_day["xc"].values, left=np.nan, right=np.nan)
        yc_i = np.interp(valid_depths, depths, df_day["yc"].values, left=np.nan, right=np.nan)
        diffs_xc[key] = pd.Series(np.diff(xc_i), index=valid_depths[:-1])
        diffs_yc[key] = pd.Series(np.diff(yc_i), index=valid_depths[:-1])

    df_x_all = pd.DataFrame(diffs_xc)
    df_y_all = pd.DataFrame(diffs_yc)

    for ref_idx in range(num // 2, len(day_nums) - num // 2):
        df_x = df_x_all.iloc[:, ref_idx - num // 2 : ref_idx + num // 2 + 1]
        df_y = df_y_all.iloc[:, ref_idx - num // 2 : ref_idx + num // 2 + 1]
        df_data = pd.DataFrame(index=df_x.index)
        df_data["dxc"] = df_x.mean(axis=1)
        df_data["dyc"] = df_y.mean(axis=1)
        df_data["sum_dxc"] = df_data["dxc"].cumsum()
        df_data["sum_dyc"] = df_data["dyc"].cumsum()
        df_data["total_var"] = df_x.var(axis=1) + df_y.var(axis=1)
        df_data["weight"] = 1 / (df_data["total_var"] + eps)
        df_data[var] = df_data.index
        df_data = df_data.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["sum_dxc", "sum_dyc", var, "weight"]
        )

        day = int(day_nums[ref_idx])
        if df_data.empty or df_data[var].max() - df_data[var].min() < min_depth_range or len(df_data) < min_points:
            rows.append({"Eddy": eddy, "Day": day, "TiltDis": np.nan, "TiltDir": np.nan})
            continue

        x = df_data["sum_dxc"].values
        y = df_data["sum_dyc"].values
        z = df_data[var].values
        weights = df_data["weight"].values
        weights = weights / np.nanmax(weights)
        total_weight = np.sum(weights)
        if not np.isfinite(total_weight) or total_weight <= 0:
            rows.append({"Eddy": eddy, "Day": day, "TiltDis": np.nan, "TiltDir": np.nan})
            continue

        mean = np.array([np.dot(weights, x), np.dot(weights, y), np.dot(weights, z)]) / total_weight
        points = np.vstack((x, y, z)).T
        try:
            _, _, vt = np.linalg.svd((points - mean) * np.sqrt(weights)[:, None], full_matrices=False)
        except np.linalg.LinAlgError:
            rows.append({"Eddy": eddy, "Day": day, "TiltDis": np.nan, "TiltDir": np.nan})
            continue

        direction = vt[0]
        if not np.all(np.isfinite(direction)) or abs(direction[2]) < eps:
            rows.append({"Eddy": eddy, "Day": day, "TiltDis": np.nan, "TiltDir": np.nan})
            continue

        z_top = np.nanmin(z)
        z_bottom = np.nanmax(z)
        p_top = mean + ((z_top - mean[2]) / direction[2]) * direction
        p_bottom = mean + ((z_bottom - mean[2]) / direction[2]) * direction
        rows.append(
            {
                "Eddy": eddy,
                "Day": day,
                "TiltDis": float(np.hypot(p_top[0] - p_bottom[0], p_top[1] - p_bottom[1])),
                "TiltDir": float((bearing(p_bottom, p_top) + bearing_offset) % 360),
            }
        )

    return pd.DataFrame(rows, columns=["Eddy", "Day", "TiltDis", "TiltDir"])
