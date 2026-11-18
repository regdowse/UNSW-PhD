# Coordinated manuscript revision

Status: coordinated revision drafts prepared on 2026-11-17, following user approval.
Both text repositories use `revision/analysis-refresh-2026-11-17`; publication and
validation status are recorded below. This is the single working checklist for the
main paper, supplement and supervisor feedback.
Repository locations and starting revisions: [resources.md](../resources.md).

## Working rules

Describe the dataset-producing code/configuration in the supplement. Check
provenance rather than assuming current public code and the local producing
version are identical. Align the main paper with agreed analyses and current
outputs. Bring substantive code/text discrepancies and proposed structural
changes back for discussion before implementation.

The user approved the proposed flow and text corrections, accepted all pending
tracked suggestions, and reserved figure presentation changes for their own work.
Preserve original exports and backups. Edit review branches only; merging and
Overleaf synchronisation remain separate actions.

## Read-only review

- [x] Confirm paper/supplement revisions and match supervisor anchors to the paper.
- [ ] Record/locate the original source ZIP and verify the recovery baseline.
- [x] Establish current pipeline/config/ESP identities and compare with public code; historical exact source revision remains qualified below.
- [x] Check supplement descriptions against detection, fitting, QC, tracking,
      processing, vertical reconstruction and tilt code/settings.
- [x] Check main-paper methods, results, captions and claims against current topic reports/tables.
- [x] Check shared terminology, units, coordinates, thresholds, depth support,
      temporal treatment, weighting and uncertainty across both manuscripts.
- [x] Propose a concise paper/supplement structure, including what to move or shorten.
- [x] Agree the discrepancy resolutions and revision scope before editing either text.

## Findings and proposed changes

Evidence baselines: paper `67ddd6d`, supplement `44623fb`, public pipeline `f95e42b`;
full hashes/locations are in [resources](../resources.md). P/S line numbers below
refer to `agujournaltemplate.tex` / supplementary `main.tex` at those revisions.
The following findings were recorded in the read-only review; their status column
now records the draft response. Analysis definitions and original feedback are unchanged.
No literature-wide novelty audit or journal-format compliance check was performed.

### Main paper

| ID | Text / section | Finding and code/result/comment evidence | Proposed action | Decision / status |
|---|---|---|---|---|
| M01 | Key points, abstract, Results; P48,75,436–508 | Old 2,982 eddies / 127,426 days / 101,656 valid tilts. [Current census](../../census/outputs/report.md): 3,097 / 133,038 / 106,865 (80.33% of all days), AE/CE median 18.62/15.48 km and median normalised tilt 0.228/0.214. Regional counts, tails, directional fractions and per-eddy comparisons also need replacement. | Update every occurrence from the current CSVs, with explicit denominators; retain per-eddy comparison panels. | Addressed in paper draft; analysis unchanged. |
| M02 | Tilt availability; P419,436; C010 | A failed/unavailable tilt estimate is not a zero-tilt eddy. The paper also quotes an edge-adjusted 'possible estimates' denominator absent from current census outputs. | Separate estimator coverage from tilt prevalence; use all-day coverage, or explicitly derive and justify a different denominator before quoting it. | Addressed in paper draft; analysis unchanged. |
| M03 | Regions; P440,475,484; C011 | Regional labels are not defined in methods. `census.region_grid` uses the shared -33°/157° boundaries, a restricted <4000 m/<400 km/<154.85°E shelf mask plus 80 km expansion. Upstream pools S1/U1/U2; downstream S2/D1/D2. | Give a short operational definition and map reference before regional Results. | Addressed in paper draft; analysis unchanged. |
| M04 | Direction methods/caption; P468 | Caption implies magnitude-weighted mean tilt vectors; `census.per_eddy` uses circular means of unit daily directions. Directional population also excludes tilt <5 km. | State equal-day unit-vector circular mean and directional QC, consistently in methods and caption. | Addressed in paper draft; analysis unchanged. |
| M05 | Delta methods; P362–419; C005–C009 | Core five-day Gaussian/SVD description agrees with code, including its horizontal-km/depth-m numerical convention. Missing: reference-day interval support restriction before averaging/variance/cumsum. Neighbours cannot extend a reference profile. | Add that restriction; explain SVD extracts a representative axis and endpoint distance; streamline notation and explicitly define sample variances. Preserve the existing algorithm and units. | Addressed in paper draft; analysis unchanged. |
| M06 | Figure 3 caption; P424 | Editorial REG NOTE and ±5 km characterisation remain. Updated pinned example 2191/7923 gives 3.593093 km tilt / 14.498783 km projected extent; [report](../../tilt_method/outputs/report.md). Ribbon is fit weight, not a CI. | Replace image and caption; retain the definition of projected extent and avoid claiming one selected noisy profile establishes robustness. | Addressed in paper draft; analysis unchanged. |
| M07 | PV methods; P145–357; C004 | Surface geometry/amplitude, Laplacian reconstruction and weighted completed vectors broadly agree with current code. Q is not sufficiently defined. Omission of internal ∇ζ/h is a diagnostic choice despite reconstructing spatially varying ζ; it should not imply that field is literally uniform. 'Local grid angle' is actually one reference angle; core-mean water depth drives target selection. | Define Q and the environmental proxy succinctly; retain surface amplitude, raw magnitude classification, local-before-mean calculation and distinction from full Ertel PV. | Addressed in paper draft; analysis unchanged. |
| M08 | Target populations; P528,533,546,751 | Text conflates magnitude category with a depth-selected population and uses strict >/< wording. Current code uses magnitude-only categories, then footprint-mean h≥3000 m / h≤3000 m and tilt≥5 km; excluded other-depth dominated days are not mixed. | State category and target filters separately; retain all-depth Figures 7–9; give sample counts from current reports. | Addressed in paper draft; analysis unchanged. |
| M09 | Figure 12 and interpretation; P549,553 | Caption defines on/off shelf by the 3 km isobath; actual Figure 12 split is longitude 154.75°E. It also writes a non-wrapped absolute angle difference and calls opposition 'alignment'. | Describe the actual longitude groups (or agree a later analysis change); use shortest circular separation 0–180°; call ≥135° opposition. Keep footprint depth selection distinct from the local plotted contour. | Addressed in paper draft; analysis unchanged. |
| M10 | Background methods/results; P536,540,700 | 'Centred three-month window' is imprecise; code uses a centred 91-day seasonal window of monthly climatological fields. The 200–500 m expression needs thickness weighting. Current daily contrast medians are AE 0.0473 / CE 0.0485 m/s, not ~0.06; [background report](../../background/outputs/report.md). | Define monthly climatology, layer coverage, (500U0–500−200U0–200)/300 and surface-minus-deep contrast in m/s. Separate median contrast from difference of medians. Keep 31/61/91 checks internal as agreed. | Addressed in paper draft; analysis unchanged. |
| M11 | Controls and ellipse methods; P565–602,588,707 | Current matched samples, exact depths, support thresholds, day weighting and track bootstrap are largely missing from methods. Ω is angular velocity (s⁻¹), not rotational speed. Five ellipse levels end at 327.44 m, not 300 m. Zero mean cos(2φ) is not an environment-matched independence test. | Add concise analysis protocol; update R²/medians from [controls tables](../../controls/outputs/), label Ω and exact analysed depth range, distinguish IQR from uncertainty. | Addressed in paper draft; analysis unchanged. |
| M12 | Latitude/β interpretation; P565,568,694 | Latitude and β are intrinsically linked; the binned-median regressions do not identify β independently of geography. Current fits are descriptive, unadjusted associations. | Retain observed trends, soften the claim that their correspondence rules out latitude being a geographic descriptor; do not turn binned R² into individual-day evidence or causal attribution. | Addressed in paper draft; analysis unchanged. |
| M13 | Composite methods, caption, discussion; P751,758 | [Composites](../../composites/outputs/report.md) are surface-to-depth offsets, opposite to saved deep-to-shallow tilt. At 858.92 m planetary AE/CE along-PV means are -6.63/+11.41 km; topographic CE +11.90 km. Current prose describes them as AE with / CE against without defining the reverse convention. Also shares the saved tilt≥5 km selection. | Explain own-surface referencing, signed PV frame, depthwise day means/support and track CIs. Translate signs before comparing to tilt; describe complementary geometry evidence, not an independent validation. | Addressed in paper draft; analysis unchanged. |
| M14 | Slope mechanism; P521,553,713–724,766 | Opposition to the total proxy gradient is not universally equivalent to motion up a continental slope. Bathymetric direction, signed f+ζ, cancellation and the chosen geographic groups matter. AE 'loss' of preference is descriptive, not a tested transition/competition mechanism. | Keep slope competition explicitly as a hypothesis; restrict slope-direction statements to supported environments and avoid implying a causal or track-transition test. | Addressed in paper draft; analysis unchanged. |
| M15 | Figure placement, refs, placeholders | Figures 15/16 appear in Discussion after Results use them. `\ref{fig:composites}` is undefined because caption uses `\label{fig}`. Plain-language summary, hypothesis caption, dataset subsection, Open Research Statement and appendix contain placeholders/internal notes. | Move result figures beside their results, repair the label, replace placeholders, add formal supplement/code/data references; retain excluded Figure 1 but finish its caption. | Addressed in paper draft; public product archive remains pending. |
| M16 | Figure 7 caption; P490 | The caption says "same as" the tilt rose figure. Figure 6 requires directional tilt QC; Figure 7 uses valid positive PV vectors across all depths without requiring a valid tilt. | Explicitly distinguish populations/scales and display-only grid rotation; do not imply identical samples. | Addressed in paper draft; analysis unchanged. |
| M17 | Time coverage; P75,108,129 | The inspected archive contains 9,189 consecutive daily records, 1994-01-02 to 2019-02-28. "26-year" denotes inclusive calendar-year labels, not 26 complete years of analysed daily output. | State the actual analysed dates and distinguish simulation coverage from available archive coverage. | Addressed in paper draft; analysis unchanged. |

### Supplement and producing stack

The sequence and most configured thresholds are already broadly faithful to the
pipeline. This is a compression/precision task, not a reason to redesign the code.

| ID | Text / section | Finding and code/result evidence | Proposed action | Decision / status |
|---|---|---|---|---|
| S01 | Source/config provenance; S22,169–171 | Public/current checkout is pipeline `f95e42b` (release 0.1.1). New checkout config uses `paths.esp_root`; older workspace copy uses `esp_zonodo`. Nineteen shared Python files are byte-identical; the ESP loader, config, package metadata and runner notebook differ. User intended the new ESP directory; PBS log `seacofs_eddy.o9801291` confirms both functions loaded there. Current functions match ESP v1.1 blob `2c75ec9` (SHA-256 `1c6b574f…212da9d`). | Cite public pipeline and version-specific ESP, pin used settings. Current configuration/hash comparison and historical path evidence are strong; logs/Parquet do not independently record the complete historical commit/config hashes. Do not manufacture them. | Addressed; historical source-hash limitation retained in resources. |
| S02 | Coordinates/depths; S31,136 | Code rotates geographic velocity into the model-aligned analysis frame using the single reference grid angle. Actual default depths are abs(z_r[150,150,1:]), shared physical targets, not each cell's native depth labels. | State these transformations and target levels precisely; distinguish the selected fitting depths from native sigma levels and the later 10 m Delta grid. | Addressed in supplementary draft; pipeline unchanged. |
| S03 | Detection/fitting/tracking; S29–113 | 1 km VG grid, a=4/b=3, ±30 km transects, bounded observation-selection radius and normalised KD-tree match largely agree with code. 200 km caps the selected fitting domain, not necessarily the optimised Rc. Look-back prioritises the nearest earlier available day before distance within that day. | Keep concise operational descriptions and one threshold table; do not call tracking a global optimum or a strict daily-continuity guarantee. | Addressed in supplementary draft; pipeline unchanged. |
| S04 | Post-processing; S118–127 | Only Q components have a centred three-observation mean (min_periods=1), not all ellipse/radius variables or necessarily consecutive calendar days. Duration filter is max(Day)−min(Day)≥21, not 21 observed samples. AR failures invalidate parameters before interpolation. | Describe the actual smoothing, elapsed-span filter, gap filling and final recalculation order. Distinguish upstream preprocessing from the unsmoothed daily PV classification. | Addressed in supplementary draft; pipeline unchanged. |
| S05 | Processed vorticity; absent after S127 | `stages/processing.py` finally recomputes w=Ω(q11+q22) after smoothing/interpolation and recomputes AR. Surface vorticity is not simply an unchanged raw DOPPIO output. The seeded surface profile uses these processed values; deeper fits are fitted separately. | Add one explicit sentence defining the final fields used by paper PV/controls; keep w, Ω0 and Ω distinct. | Addressed in supplementary draft; pipeline unchanged. |
| S06 | QC and failure semantics; S85–95,149–162 | Surface cleaning explicitly requires finite Ω0/sign/error, not finite every outer-fit parameter. ESP v1.1 can return NaN outer fits without raising; a row/centre can therefore survive. Full vertical reconstruction stops on `None` from continuity/core/fitting checks, but NaN returned outer parameters do not automatically stop it. Ten-point guards count geometric domain points before ESP finite/rotation filtering. Confirmed-profile QC is ≥3 recorded levels; tilt eligibility is stricter. | Separate centre-following acceptance, optional outer-fit availability, profile-count QC and downstream analysis QC. Do not claim every missing outer fit removes an eddy/day. | Addressed in supplementary draft; pipeline unchanged. |
| S07 | Physical parameters and units; S42,65–80,125 | R is the mean available directional peak-speed radius, not Rc. Pipeline coordinates/radii are km; angular quantities receive 10⁻³ scaling. ψ0 is retained in the implementation's mixed horizontal units; cutoff 300 must not be labelled SI without checking conversion. Poor outer fits may fall back to Ω0 and interpolate missing Rc/ψ0 later. | Define symbols/units and the fallback/interpolation sequence; distinguish observation-domain bounds from physical radius estimates. | Addressed in supplementary draft; pipeline unchanged. |
| S08 | Reproducibility/link; S22,171 | Intro links the older UNSW-PhD folder. Resume-by-existence is conditional, not input/config validation; supplied config has skip_existing=false and QC/consolidation regenerate. | Link/cite ocean-eddy-pipeline and pinned configuration/ESP; describe rerun behaviour accurately without shell/path clutter in the scientific narrative. | Addressed in supplementary draft; pipeline unchanged. |
| S09 | Structure/figures; S12–38,128–132,165–166 | Rough working title, three dot-only figure captions, no established updated supplementary-figure provenance, and 'Detailed in paper' placeholder. These figures were not included in the main Figure 2–17 refresh. | Convert into supporting information with one compact parameter table and focused stage descriptions; retain supplementary figures only if useful and verified against current dataset. Cross-reference the main Delta method. | Text rewritten; three unverified rough figures omitted, files retained. |

## Approved flow

Main paper: concise motivation/gap/objectives → model and dataset overview with
supplement reference → Delta estimator → PV/regions/populations and concise analysis
protocol → census → polarity/PV regimes and background-flow association → eddy
structure/ellipse associations → discussion, limits and conclusion. Put each result
figure with its results; consider introducing the composite beside the PV-regime
results and referring back in limitations. Keep Figure 1 explicitly conceptual.
Remove the duplicated Knowledge gap subsection rather than maintaining both versions.

Supplement: scope/data coordinates and settings → detection/surface fitting/QC →
tracking/post-processing → vertical reconstruction/profile QC → reproducibility and
versioned resources. Keep Delta derivation in the main paper and cross-reference it;
do not duplicate the long PV derivation. No new analysis is implied by this proposal.

## Agreed revision scope and remaining work

- Paper: update methods, numbers, captions, interpretation and flow from current
  outputs. Keep the existing Figure 12 longitude partition and describe it explicitly.
- Supplement: concise code-faithful Text S1 and settings table. The three rough
  supplementary images remain in the repository but are omitted from the active
  draft because their updated provenance has not been established. Revisit if needed.
- Keep all 16 refreshed main-paper assets unchanged. The user will handle C012–C015
  presentation requests. Introducing Delta before PV swaps their automatic numbers
  to Figures 2 and 3 respectively; other figure numbers and all filenames are stable.
- Pending submission work: persistent public archive/DOI for processed data and
  analysis products, final journal formatting, figure presentation, and locating
  the downloaded original ZIP. No new analyses were introduced.

## Draft validation and publication

- Recovered the completed edits from persistent session records after a compute-node
  change made the original node-local `/tmp` checkouts inaccessible. The editing
  repositories now live in shared `manuscript/text/`; see [resources](../resources.md).
- All 16 main-paper PNG files match current analysis exports byte-for-byte.
- Numerical anchors, author/affiliation wording, citation keys, labels and image
  paths checked. Both repository diffs pass whitespace checks.
- Both branches published and remote refs verified: paper `1f34525ebe316a49bd8d56bb1357a8fd4e73abcf`;
  supplement `db362ad1be810394e4b2fb6f2e33d4be64d5bee3`.
- PDFs rebuilt successfully: [paper](outputs/main-paper.pdf), 18 pages;
  [supplement](outputs/supplement.pdf), 3 pages. No unresolved citations/references
  or missing glyphs. One 2.91 pt overfull box remains in a model paragraph; final
  layout review is separate. [Validation record](outputs/validation.json) includes
  commits, artifact hashes and unchanged feedback fingerprints. Outputs are ignored.
- Original `main` baselines and Overleaf are unchanged. Original export files are
  retained verbatim; addressing text does not resolve a thread in Overleaf.

## Supervisor feedback

Export inventory: 14 open comment threads and eight tracked changes. Original
files are in [reviews/overleaf-comments](../reviews/overleaf-comments/).
Export stable IDs identify entries; do not rename the original IDs or edit exports.

| Export ID | Passage / section | Relevant finding or revision | Disposition |
|---|---|---|---|
| C001 | Introduction: systematic quantification claim | M12/M14 and proposed introduction; distinguish EAC-specific gap from global literature. Literature verification still required. | Addressed in draft; original Overleaf thread retained. |
| C002 | Duplicated Knowledge gap | Proposed flow removes repetitive subsection, retains one gap/objective passage. | Addressed in draft; original Overleaf thread retained. |
| C004 | Definition of Q | M07; define shape matrix/eigenstructure in the paper, not only via a previous publication. | Addressed in draft; original Overleaf thread retained. |
| C005 | 'Vertical eddy-center positions' | M05; use horizontal centre coordinates as a function of depth and time. | Addressed in draft; original Overleaf thread retained. |
| C006 | Delta notation | M05; propose x(d,t), y(d,t), with explicit increment intervals. | Addressed in draft; original Overleaf thread retained. |
| C007 | Translation explanation | M05; shorten to one purpose clause if retained. | Addressed in draft; original Overleaf thread retained. |
| C008 | Undefined s | M05; explicit sample-variance definition and available-day support. | Addressed in draft; original Overleaf thread retained. |
| C009 | Why SVD/PCA | M05; explain representative line, deep-to-shallow direction and endpoint distance, with current unit convention. | Addressed in draft; original Overleaf thread retained. |
| C010 | Missing estimates versus untilted days | M02; estimator coverage is not a physical zero-tilt classification. | Addressed in draft; original Overleaf thread retained. |
| C011 | Upstream/downstream definitions | M03; operational regions before Results, with map reference. | Addressed in draft; original Overleaf thread retained. |
| C012 | Enlarge planetary figure | M15; layout/readability decision when placing Figure 10. | Deferred to user figure presentation work. |
| C013 | Equivalent topographic roses | M09; absolute PV direction varies, so a relative-angle rose is possible; no extra panel produced. | Deferred to user figure presentation work. |
| C014 | Relative-angle wind rose | M09; current figure uses absolute separation. A directional relative rose would require agreeing signed versus folded angles and preserving numerical meaning. | Deferred to user figure presentation work. |
| C015 | Positive South latitude | M11/M12; display-axis relabelling can preserve calculations, but slope interpretation/sign must remain explicit. | Deferred to user figure presentation work. |
| T001 | Insert Christopher C. | Present in current P18; preserve during editing. Review state is not encoded in plain Git text. | Accepted by user; intent retained in revised wording. |
| T002 | Insert affil{4} | Present in current P18; preserve. | Accepted by user; intent retained in revised wording. |
| T003 | Insert Canberra affiliation | Present in current P31; preserve. | Accepted by user; intent retained in revised wording. |
| T004 | Insert 'well validated' | Reflected in P75; retain intent, verify supporting validation citations during final wording. | Accepted by user; intent retained in revised wording. |
| T005 | Insert 'of the water column' | Reflected in P75; preserve meaning while correcting sampled depth/support. | Accepted by user; intent retained in revised wording. |
| T006 | Delete 'strongly' | Target abstract clause currently omits it; retain moderated wording. | Accepted by user; intent retained in revised wording. |
| T007 | Insert 'the' | Reflected in P105; preserve when consolidating introduction. | Accepted by user; intent retained in revised wording. |
| T008 | Delete 'its' | Target clause in P105 reflects the deletion. | Accepted by user; intent retained in revised wording. |

Use `pending`, `discuss`, `addressed` or `retained with reason`. An addressed
checklist entry does not resolve the original Overleaf thread. All 14 anchors occur
exactly once in the paper snapshot; C003 is absent from the export and has not been
invented. All eight tracked suggestions are reflected in current source wording,
and the user has accepted all eight for this revision. Their actual Overleaf UI
state is separate. Original threads and exported files remain unchanged; no item
has been marked resolved in Overleaf.

## Text comparison (2026-11-17)

- User pushed annotated original tags and merged the paper revision into `main`.
  Fetched paper `main` is `067fbea233add3f91ece28979e2cb38f8c0b257e`; local
  paper checkout fast-forwarded to it without editing source.
- [Text comparison PDF](outputs/main-paper-latexdiff.pdf): original tagged
  `67ddd6d` versus GitHub `main` `067fbea`, 29 pages. Figure images omitted as
  requested; captions, equations, tables and key points retained. Blue underline
  denotes additions; red strikeout denotes deletions. Moved passages can appear
  deleted and added. Review layout and revised bibliography are used.
- Build script, source snapshots, logs and commit/hash metadata are retained in
  `outputs/main-paper-latexdiff/`. Compilation passed with no unresolved references
  or citations; no images embedded. Two modest overfull boxes remain within page
  margins. Original source repositories and feedback exports were not edited.
- Latest Overleaf push was not confirmed for this task: this compares the latest
  fetched GitHub version. Before future manuscript tasks, follow the push reminder
  in [context README](../README.md).
