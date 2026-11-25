# Supplement scientific review

Reviewed 2026-11-25. Baseline supplement `dc9b106e66062d90ac6ac2df931cdc85301579fc`.
User requests discussion of scientific issues before any proposed replacement
text; subsequently only a text-only latexdiff PDF is needed, not a clean PDF.
Keep proposed sources outside manuscript checkouts. No push authorised.

## Sources checked

- ESP methods preprint https://doi.org/10.5194/egusphere-2026-3599,
  especially equations 2, 6–7, 10–13 and the outer-fit discussion.
- ESP v1.1 Zenodo archive https://doi.org/10.5281/zenodo.24182671.
  Download checksum verified: MD5 94d44b1651389ba6c2950d9940b65289.
  Archived functions.py SHA256
  1c6b574f1084e9ee4eadfd01553464ff6fc5ae514e60180e7e804d09f212da9d;
  exactly matches current external implementation. Temporary archive/functions
  in /tmp/esp-v1.1-review.zip and /tmp/esp-v1.1-review-functions.py.
- Pipeline archived commit 143b9a31ec0586206cf440df3da129fb886bc4b8;
  current HEAD 6bbb222 has no src/config differences from that commit.
- Source tracked and processed Parquet SHA256 values match the release audit.
  Verification JSON in ignored outputs/supplement-scientific-audit-20261125/.

## Scientific distinctions to discuss

1. Directional radii and outer fitting use two different velocity quantities.
   core/doppio.py projects observed velocity onto J Q r / ||Q r||. Searches start
   at neighbours of the nearest original grid point in four model-grid directions.
   The first decrease in |v_t| selects the preceding sample; its Euclidean
   distance from the fitted centre is the directional radius. This is a local
   sample-based peak, not a global maximum or a fitted R_c/sqrt(2). Mean available
   directional radii gives Rbar. Outer fitting additionally multiplies v_t by
   rho/||Q r||, then fits Omega*rho*exp(-rho^2/Rc^2). Need both definitions.
   Rbar is not an elliptical rho coordinate. Consequently clipped
   rho_lim=max(30,min(1.75 Rbar,200)) km targets outer observations but does not
   guarantee coverage beyond every true/fitted peak-speed contour. Data are
   also clipped by a rectangular extraction half-width 3*rho_lim.

2. Original DOPPIO Q is SPD with determinant one, but component-wise temporal
   smoothing/interpolation does not preserve unit determinant. In 133,038 released
   processed surface records, determinant quantiles min/median/95%/max are
   1.0000059845 / 1.0223646861 / 1.2199734886 / 6.8859273578. None is close to one
   at rtol=1e-6, atol=1e-8. Describe original fits versus processed matrices;
   do not claim area-equivalent rho for every processed surface record, or
   silently renormalise Q (this also affects the Omega/Rc representation).

3. Vorticity provenance differs across stages. Surface original fit w is
   Omega0*trace(Q); processed surface w is overwritten with Omega*trace(Q).
   Subsurface recorded w remains the DOPPIO inner-core vorticity, even when the
   outer fit yields a different refined Omega. Surface profile row inherits
   processed surface values. Thus a universal released-profile statement
   w=Omega*trace(Q) would be wrong. Document this if describing vorticity.

4. Preliminary sampling tests are not effective fit eligibility. Pipeline guard
   is >=3 geometric samples per transect, but archived DOPPIO cleans jointly
   finite coordinates and both velocity components, requires >=4 finite samples
   per transect, exactly one common intersection, successful optimisation and
   valid elliptical geometry. Rectangular/elliptical >=10 point checks are
   geometric, prior to outer finite/sign filtering. Outer fit requires >=3
   retained observations, distinct effective radii and identifiable Gaussian
   curvature/covariance. No numerical universal R2 acceptance threshold.

## Text corrections that do not need new scientific choices

- Confirm inner outputs centre/Q/Omega0 and outer outputs Rc/refined Omega.
  Outer fit holds centre and Q fixed, filters the selected velocities by
  circulation sign, and refines both parameters jointly. Pipeline passes core
  vorticity as the initial angular-velocity argument; both have the same sign,
  but do not claim the initial guess is literally Omega0. No reason to change
  the published calculation during this text review.
- Keep psi0 out of explanatory text as explicitly requested. Code does have an
  amplitude threshold, |psi0|>300 stored units. Replaying the archived pre-fill
  duration/shape screens on tracked rows finds 1,445 amplitude exceedances and
  ZERO for which amplitude alone triggers outer screening. Thus visible Rc and
  Omega thresholds reproduce the same outer-screen mask for this dataset;
  the full amplitude rule stays in the dataset's technical documentation.
- Surface outer screening: Rc>300 km OR Rc>1.75 Rbar OR |Omega|>5e-5 s^-1
  (plus the amplitude rule above). These are post-fit processing limits, distinct
  from the 200 km observation-selection cap. Failed finite fits trigger fallback
  Omega=Omega0 and missing outer scale; a routine returning NaNs is distinct
  and subsequently subject to filling. Do not imply every missing fit uses
  the fallback.
- Gap filling covers centres, rotation, Q and Rbar as well as outer parameters;
  endpoint fills are used, with no maximum outer gap in this configuration.
- Shallow surface QC: fixed original VG grid indices, not grid nearest refined
  DOPPIO surface centre. Starts comparison from refined surface centre and checks
  successive levels through 60.17 m; 100 km applies per step, not net displacement
  from surface. User rationale: avoid progressive migration into another eddy.
  Full vertical reconstruction instead advances nearest-grid transect location.
- Full reconstruction stops for support/boundary/geometry/radius/continuity
  failures or raised exceptions. Outer function returning NaNs without raising
  does not by itself stop continuation or remove an accepted centre.
- Count-QC is >=3 recorded levels including surface; it does not require a
  successful outer fit or certify Delta eligibility. Keep this distinction.
- Preserve existing references, 21 elapsed-day duration and >15 km mean-radius
  selection. Avoid presenting 15 km as a universal mesoscale definition.
- Stage-numbered table near end should mirror prose definitions and thresholds.
  Keep outer fit and subsequent surface-fit screening as distinct stages.

No manuscript, pipeline, ESP, dataset or published metadata changed. No revised
manuscript source or PDF has been produced in this audit.

## Follow-up on Q and subsurface rotation

User accepts items 1/4 and explains retaining core vorticity beneath surface as
protection against increasingly unavailable outer fits. They ask whether Q changes
are negligible and whether raw matrices survive in the published dataset.

Raw Q survives in surface_eddies/ and tracked/; processed/ overwrites its own Q
columns with smoothed/interpolated values. No paired raw/processed Q fields in
the same final row. Final IDs differ from tracked IDs; surface rows in released
profiles copy processed Q, while new subsurface fits use fresh DOPPIO Q.

Diagnostic results in outputs/supplement-scientific-audit-20261125/
q_processing_diagnostics.json. Matched 117,954 processed records to original
tracked Q by exact Day/xc/yc after deduplicating identical raw coordinate/Q rows;
no ambiguous differing-Q matches remained. Median/p95 absolute AR changes are
0.06947/0.52837; relative changes 4.5686%/25.2609%. These describe full processing,
including interpolation of screened fits, not determinant normalisation alone.
Pure scalar determinant normalisation preserves AR and eigenvectors exactly.

Holding Rc fixed, compare processed Q to Qhat=Q/sqrt(detQ), with the same shape
and orientation. Actual semiaxes are detQ^(-1/4) times those for Qhat: median/p95
shrink 0.5514%/4.8492%, maximum38.2681%. Area shrink median1.0998%, p959.4633%,
maximum61.8918%. This is a normalisation diagnostic, NOT raw-to-processed geometry
change or a PV sensitivity result.

Important: manuscript pv/pv_kernels.py uses processed surface Q for the Gaussian
envelope, footprint selection and Gaussian-Laplacian vorticity reconstruction.
Thus current use extends beyond AR. Non-unit determinant alone does not make
that reconstruction mathematically invalid: formulas hold for SPD Q; its
normalisation/scale convention differs from canonical unit-det ESP. No PV
sensitivity has been run. Do not claim negligible effects on all conclusions.
Replacing processed Q alone with raw/unit-det matrices changes the footprint
and reconstructed field; do not recommend this as a universal correction.

Subsurface core rationale explains stored w=Omega0*trace(Q). It does not mean
released subsurface Omega is Omega0: Omega is the outer-refined fit, may be NaN,
and controls/controls.py + control_tools.py use that refined Omega for structural
angular-velocity comparisons with fit QC and matched days. Clarify in discussion.

## Agreed outcome and review artifact

User explicitly chose to RETAIN existing Q smoothing and interpolation after
discussion. Explain processed non-unit determinant (median 1.02), retain the
SPD ellipse interpretation, avoid claiming universally negligible effects, and
distinguish fitted versus processed parameter sets. No dataset rebuild or code
change. Other audit corrections accepted. User requested only a text-only
latexdiff PDF after the scientific discussion, not a clean final PDF.

Prepared a separate proposed.tex and baseline.tex under ignored
outputs/supplement-review-20261125/. Baseline exactly matches dc9b106 and
references.bib is unchanged. Revised definitions, geometry/normalised velocity
equations, outer observation selection, acceptance/termination enumerations,
gap-filling description, vorticity provenance and Q caveat. Numbered 11-stage
summary table moved to end and mirrors prose. No psi0 variable or figures added.

Delivered supplement-latexdiff.pdf (7 pages), blue underlined additions and red
struck deletions, against dc9b106. Compiled with latexdiff/latexmk/BibTeX using
texlive/20220321; no unresolved citations/references or overfull boxes. One minor
underfull spacing warning in table; selected equation/table pages visually checked
and all pages checked for out-of-page text. verification.json records checks.
Only diff PDF generated; no clean PDF. Manuscript checkout remains unchanged;
no commit/push, pipeline change, dataset change or regeneration.

User refinement: omit rectangular preselection and its 3*rho_lim half-width
from the proposed prose, continuation list and table as implementation detail.
Retain elliptical rho selection and its ten-geometric-point requirement.
Technical documentation/code remain unchanged. Rebuilt same latexdiff PDF against
dc9b106; compile/page-bound checks pass; manuscript checkout remains clean.
