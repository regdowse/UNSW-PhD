"""Helpers for surface PV-gradient footprint-scale sensitivity analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_CACHE_ROOT = Path(
    "/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/"
    "pv_gradient_surface_footprints"
)
DEFAULT_CACHE_NAME = "surface_pv_footprint_sensitivity.parquet"


def default_footprints():
    """Prespecified filled ellipses, annuli, and legacy comparisons."""
    rows = []
    for frac in (0.25, 0.5, 0.75, 1.0, 1.25):
        rows.append({"footprint": f"filled_{frac:g}", "kind": "Filled", "frac": frac,
                     "inner_frac": 0.0, "averaging": "nonlinear"})
    for inner, outer in ((0.25, 0.5), (0.5, 0.75), (0.75, 1.0), (0.5, 1.0)):
        rows.append({"footprint": f"annulus_{inner:g}_{outer:g}", "kind": "Annulus",
                     "frac": outer, "inner_frac": inner, "averaging": "nonlinear"})
    for frac in (0.5, 1.0):
        rows.append({"footprint": f"legacy_{frac:g}", "kind": "Legacy", "frac": frac,
                     "inner_frac": 0.0, "averaging": "legacy"})
    return pd.DataFrame(rows)


def build_cache(eddies, grid, specs=None, *, progress=True):
    """Calculate every prespecified surface footprint and concatenate results."""
    import seacofs_tilt_tools as tilt

    specs = default_footprints() if specs is None else specs.copy()
    tables = []
    for number, spec in enumerate(specs.itertuples(index=False), start=1):
        if progress:
            print(f"[{number}/{len(specs)}] {spec.footprint}")
        table = tilt.add_pv_gradient_terms(
            eddies, grid, core_mean=True, frac=spec.frac,
            inner_frac=spec.inner_frac, averaging=spec.averaging,
        )
        table["footprint"] = spec.footprint
        table["footprint_kind"] = spec.kind
        tables.append(table)
    return add_analysis_terms(pd.concat(tables, ignore_index=True))


def cache_path(cache_root=DEFAULT_CACHE_ROOT):
    return Path(cache_root) / DEFAULT_CACHE_NAME


def load_cache(cache_root=DEFAULT_CACHE_ROOT):
    import seacofs_tilt_tools as tilt
    return add_analysis_terms(tilt.read_table(cache_path(cache_root)))


def add_analysis_terms(df, dominance_factor=2.0):
    """Add exposure, cancellation, regime, and polarity-aware tilt diagnostics."""
    out = df.copy()
    threshold = np.log(float(dominance_factor))
    with np.errstate(divide="ignore", invalid="ignore"):
        out["topo_plan_exposure_ratio"] = np.log(
            out.PV_grad_topo_mean_local_mag / out.PV_grad_plan_mean_local_mag
        )
        out["topo_cancellation_gap"] = np.log(
            out.PV_grad_topo_mean_local_mag / out.PV_grad_topo_mag
        )
    out["net_regime"] = pd.Categorical(
        np.select(
            [out.topo_plan_ratio <= -threshold, out.topo_plan_ratio >= threshold],
            ["Planetary", "Topographic"], default="Mixed",
        ), categories=["Planetary", "Mixed", "Topographic"], ordered=True,
    )
    raw_error = np.abs((out.TiltDir - out.PV_grad_theta + 180) % 360 - 180)
    out["preference_error_deg"] = np.where(out.Cyc.eq("CE"), np.abs(180 - raw_error), raw_error)
    out["preference_alignment"] = np.cos(np.deg2rad(out.preference_error_deg))
    out["low_coherence"] = out.PV_grad_topo_coherence < 0.35
    return out


def eddy_equal_summary(df, value, groups=("Cyc", "footprint")):
    """Summarise daily values after reducing each eddy to one median."""
    group = list(groups)
    eddy = df.groupby([*group, "Eddy"], observed=True)[value].median().reset_index()
    return (eddy.groupby(group, observed=True)[value]
            .agg(median="median", q25=lambda x: x.quantile(.25),
                 q75=lambda x: x.quantile(.75), eddies="count").reset_index())


def footprint_scorecard(df, min_tilt_km=5.0):
    """Eddy-equal detection and directional metrics for every footprint."""
    use = df[df.TiltDis.ge(min_tilt_km)].copy()
    eddy = (use.groupby(["Cyc", "footprint_kind", "footprint", "ellipse_inner_frac", "ellipse_frac", "Eddy"], observed=True)
            .agg(error=("preference_error_deg", "median"),
                 alignment=("preference_alignment", "mean"),
                 exposure=("PV_grad_topo_mean_local_mag", "median"),
                 net=("PV_grad_topo_mag", "median"),
                 coherence=("PV_grad_topo_coherence", "median"),
                 p90=("PV_grad_topo_p90_local_mag", "median"))
            .reset_index())
    return (eddy.groupby(["Cyc", "footprint_kind", "footprint", "ellipse_inner_frac", "ellipse_frac"], observed=True)
            .agg(eddies=("Eddy", "nunique"), median_error_deg=("error", "median"),
                 mean_alignment=("alignment", "mean"),
                 fraction_within_45=("error", lambda x: np.mean(x <= 45)),
                 median_exposure=("exposure", "median"), median_net=("net", "median"),
                 median_coherence=("coherence", "median"), median_p90=("p90", "median"))
            .reset_index())


def paired_with_reference(df, reference="filled_1"):
    """Attach full-ellipse reference values to every matched Eddy-Day footprint."""
    columns = ["PV_grad_topo_mag", "PV_grad_topo_mean_local_mag",
               "PV_grad_topo_p90_local_mag", "PV_grad_topo_coherence",
               "PV_grad_theta", "topo_plan_ratio", "preference_error_deg"]
    ref = df[df.footprint.eq(reference)][["Eddy", "Day", *columns]].copy()
    ref = ref.rename(columns={column: f"reference_{column}" for column in columns})
    return df.merge(ref, on=["Eddy", "Day"], how="inner", validate="many_to_one")
