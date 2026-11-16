"""Local fitted ellipse geometry; no interpolation or legacy imports."""
from __future__ import annotations
import numpy as np
import pandas as pd
KEY = ['Eddy','Day']
QCOLS = ['q11','q12','q22']
COLOURS = {'AE':'#b74646','CE':'#326d9b'}

def axial_difference(a, b):
    """Signed axis difference in [-90, 90); invariant to 180-degree flips."""
    return (np.asarray(a) - np.asarray(b) + 90.0) % 180.0 - 90.0

def ellipse_geometry(frame, bearing_offset_deg=20.0):
    """Derive geometry only for finite, positive-definite Q matrices."""
    out = frame.copy()
    q = out[QCOLS].to_numpy(float)
    good = np.isfinite(q).all(axis=1)
    matrices = np.zeros((len(out), 2, 2))
    matrices[:, 0, 0], matrices[:, 0, 1] = q[:, 0], q[:, 1]
    matrices[:, 1, 0], matrices[:, 1, 1] = q[:, 1], q[:, 2]
    # Replace missing matrices only for eigensolver safety; keep them invalid.
    matrices[~good] = np.eye(2)
    vals, vecs = np.linalg.eigh(matrices)
    good &= (vals[:, 0] > 0) & np.isfinite(vals).all(axis=1)
    ratio = np.full(len(out), np.nan)
    ratio[good] = np.sqrt(vals[good, 1] / vals[good, 0])
    major = vecs[:, :, 0]
    bearing = (np.degrees(np.arctan2(major[:, 0], major[:, 1]))
               + bearing_offset_deg) % 180.0
    bearing[~good | np.isclose(ratio, 1.0, rtol=1e-8)] = np.nan
    out["AxisRatio"] = ratio
    out["MajorBearing"] = bearing
    out["Q_valid"] = good
    return out

def build_analysis_table(surface, sampled, bearing_offset_deg=20.0):
    if surface[KEY].isna().any().any() or surface.duplicated(KEY).any():
        raise ValueError("Surface Eddy/Day keys must be complete and unique")
    if sampled.duplicated(KEY + ["ShapeDepth"]).any():
        raise ValueError("Duplicate sampled geometry keys")
    base = surface.copy()
    if "Region" in base:
        base["Sector"] = base.Region.map({"S1": "Upstream", "U1": "Upstream",
            "U2": "Upstream", "S2": "Downstream", "D1": "Downstream", "D2": "Downstream"})
    surf = base.assign(ShapeDepth=0.0, DepthLower=0.0, DepthUpper=0.0, DepthMethod="surface")
    deep = sampled.merge(base.drop(columns=QCOLS), on=KEY, how="inner", validate="many_to_one")
    combined = pd.concat([surf, deep], ignore_index=True) if len(deep) else surf.reset_index(drop=True)
    out = ellipse_geometry(combined, bearing_offset_deg)
    out["SignedAlignment"] = axial_difference(out.TiltDir, out.MajorBearing)
    out["AlignmentDeg"] = out.SignedAlignment.abs()
    out["AlignmentCos2"] = np.cos(np.deg2rad(2 * out.SignedAlignment))
    out["AlignmentSin2"] = np.sin(np.deg2rad(2 * out.SignedAlignment))
    return out.sort_values(["ShapeDepth", "Eddy", "Day"]).reset_index(drop=True)

def select_rows(table, *, directional=False, min_ar=1.1, min_tilt=5.0, max_ar=5.0):
    keep = (table.Q_valid & np.isfinite(table.AxisRatio) & table.AxisRatio.between(1, max_ar)
            & np.isfinite(table.TiltDis) & table.TiltDis.ge(0))
    if directional:
        keep &= (table.AxisRatio.ge(min_ar) & table.TiltDis.ge(min_tilt)
                 & np.isfinite(table.AlignmentDeg))
    return table.loc[keep].copy()

def matched_depth_sample(table, depths):
    """Identical eddy-days at every requested depth, after outcome-specific QC."""
    selected = table.loc[table.ShapeDepth.isin(depths)]
    counts = selected.groupby(KEY).ShapeDepth.nunique()
    keys = counts[counts == len(depths)].reset_index()[KEY]
    return selected.merge(keys, on=KEY, how="inner", validate="many_to_one")
