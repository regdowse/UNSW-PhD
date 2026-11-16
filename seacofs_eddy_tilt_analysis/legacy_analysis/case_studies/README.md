# Paper case-study series

These notebooks select representative eddies for the tilt paper. They share
`paper_case_study_tools.py` and use the corrected signed-PV convention:

- AE tilt is expected along signed `grad(PV)`;
- CE tilt is expected opposite signed `grad(PV)`.

Run the notebooks independently:

1. `01_planetary_reference_cases.ipynb` - stable deep-water background cases.
2. `02_topographic_response_cases.ipynb` - coherent CE opposition and disrupted
   AE directional response under strong topographic control.
3. `03_regime_transition_cases.ipynb` - long-lived eddies that move between the
   planetary and topographic regimes.

The defaults use a centred seven-observation median, a 2:1 dominance threshold,
at least 5 km of tilt for directional analysis, `h >= 3000 m` for planetary
reference days, and `h <= 2000 m` for topographic days. Rankings are screening
tools, not statistical evidence by themselves. Inspect every selected track for
coverage, boundary effects, and physical coherence before using it in a paper.

The older `case_study_tools.py` supports the exploratory notebooks. Its
directional diagnostics have also been updated to the corrected signed-PV
convention, although the historical output identifier
`topographic_ce_alignment` is retained for notebook compatibility and now means
coherent CE opposition to signed `grad(PV)`. New paper notebooks use the more
focused ranking workflow in `paper_case_study_tools.py`.
