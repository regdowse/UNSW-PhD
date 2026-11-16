# Planetary- and topographic-dominated eddies

This folder separates the population-level analysis that grew out of
`case_studies/pv_grad_tilt_theta.ipynb`.

- `planetary_dominated_eddies.ipynb` estimates the background or "natural"
  tilt of eddies where the planetary PV-gradient contribution dominates.
- `topographic_dominated_eddies.ipynb` tests whether an increasingly strong
  topographic PV gradient overrides those background directions and aligns
  eddy tilt with the total PV gradient.
- `rossby_number_pv_gradient.ipynb` isolates how signed Rossby number modifies
  the topographic PV-gradient term through the factor `1 + Ro`, including a
  no-relative-vorticity counterfactual and latitude sensitivity.

Both notebooks use a smoothed 2:1 dominance criterion, require at least 5 km
of tilt displacement for directional analyses, and report eddy-level as well
as observation-level summaries. Daily observations along a track are not
treated as independent replicates.
