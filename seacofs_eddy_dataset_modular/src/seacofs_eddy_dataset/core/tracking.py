from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def clean_surface_eddies(df: pd.DataFrame, X_grid, Y_grid, threshold: float = 4.0, omega0_thresh: float = 5e-5):
    required = ["nCyc", "xc", "yc", "nxc", "nyc", "Omega0"]
    clean = df.dropna(subset=required).copy()

    nenc_cyc = np.where(clean["nCyc"].eq("AE"), 1, -1)
    doppio_cyc = np.sign(clean["Omega0"].values)
    clean = clean.loc[nenc_cyc == doppio_cyc].copy()

    err = np.hypot(clean["xc"] - clean["nxc"], clean["yc"] - clean["nyc"])
    clean["err"] = err
    boundary_mask = clean["xc"].between(0, X_grid.max(), inclusive="neither") & clean["yc"].between(
        0, Y_grid.max(), inclusive="neither"
    )

    log_err = np.log1p(err)
    err_med = np.median(log_err)
    err_mad = np.median(np.abs(log_err - err_med))
    upper_err = np.inf if err_mad == 0 else np.expm1(err_med + threshold * err_mad / 0.6745)

    clean["Omega0_abs"] = clean["Omega0"].abs()
    final_mask = (
        boundary_mask
        & (err <= upper_err)
        & (clean["Omega0_abs"] <= omega0_thresh)
        & np.isfinite(err)
        & np.isfinite(clean["Omega0_abs"])
        & (clean["Omega0_abs"] > 0)
    )
    return clean.loc[final_mask].reset_index(drop=True)


def tracking_kdtree_with_omega0(
    df_data: pd.DataFrame,
    start_ids,
    next_num: int,
    length_scale: float = 50.0,
    omega0_scale: float = 1e-5,
    radius_threshold: float = 1.0,
    lookback: int = 4,
) -> pd.DataFrame:
    df = df_data.dropna(subset=["xc", "yc", "Omega0"]).copy()
    if df.empty:
        return df.assign(Eddy=pd.Series(dtype="Int64"), next_num=next_num)

    min_day = df["Day"].min()
    first_mask = df["Day"].eq(min_day)
    first_count = int(first_mask.sum())
    if np.isscalar(start_ids):
        start_ids = np.arange(int(start_ids), int(start_ids) + first_count)
    start_ids = np.asarray(start_ids)
    if start_ids.size != first_count:
        raise ValueError("start_ids must be scalar or match first-day eddy count")

    df["Eddy"] = -1
    df.loc[first_mask, "Eddy"] = start_ids
    df["Eddy"] = df["Eddy"].astype("Int64")

    unique_days = np.sort(df["Day"].unique())
    for day in unique_days[1:]:
        present = df.loc[df["Day"].eq(day)].groupby("eddy_idx", as_index=False, sort=False).first()
        assigned = set()

        for _, present_eddy in present.iterrows():
            best_id = None
            query = np.array(
                [
                    present_eddy["xc"] / length_scale,
                    present_eddy["yc"] / length_scale,
                    present_eddy["Omega0"] / omega0_scale,
                ]
            )

            for delta in range(1, lookback + 1):
                previous = (
                    df.loc[df["Day"].eq(day - delta) & df["Eddy"].ne(-1)]
                    .groupby("eddy_idx", as_index=False, sort=False)
                    .first()
                )
                if previous.empty:
                    continue

                coords = np.column_stack(
                    [
                        previous["xc"].values / length_scale,
                        previous["yc"].values / length_scale,
                        previous["Omega0"].values / omega0_scale,
                    ]
                )
                idxs = cKDTree(coords).query_ball_point(query, r=radius_threshold)
                if not idxs:
                    continue

                dists = np.linalg.norm(coords[idxs] - query, axis=1)
                for j in np.asarray(idxs)[np.argsort(dists)]:
                    previous_eddy = previous.iloc[j]
                    if present_eddy["Cyc"] == previous_eddy["Cyc"] and previous_eddy["Eddy"] not in assigned:
                        best_id = previous_eddy["Eddy"]
                        break
                if best_id is not None:
                    break

            mask = df["Day"].eq(day) & df["eddy_idx"].eq(present_eddy["eddy_idx"])
            if best_id is None:
                best_id = next_num
                next_num += 1
            df.loc[mask, "Eddy"] = best_id
            assigned.add(best_id)

    df = df.loc[df["Eddy"].ne(-1)].copy()
    if df.duplicated(subset=["Eddy", "Day"]).any():
        raise ValueError("Duplicate (Eddy, Day) pairs found after tracking")
    df["next_num"] = next_num
    return df.reset_index(drop=True)
