# Resources

Recorded 2026-11-17 (Australia/Sydney). Repository snapshots below were checked
read-only through GitHub. They are review starting points, not a claim that those
revisions produced the current dataset. Update this map when a review baseline changes.

## Text and code repositories

| Resource | Repository / entry point | Checked `main` revision |
|---|---|---|
| Main paper | [Tilting-of-Mesoscale-Ocean-Eddies-in-The-Tasman-Sea](https://github.com/regdowse/Tilting-of-Mesoscale-Ocean-Eddies-in-The-Tasman-Sea), `agujournaltemplate.tex` | `67ddd6d4e31866d1afb40003cb8b673b502012bc` |
| Supplementary manuscript | [seacofs_eddy_dataset_modular](https://github.com/regdowse/seacofs_eddy_dataset_modular), `main.tex`, `references.bib`, `Figures/` | `44623fb2404e47a1b8d6f498f5c7c88f6bcb7bc7` |
| Dataset pipeline | [ocean-eddy-pipeline](https://github.com/regdowse/ocean-eddy-pipeline), `README.md`, `config/local.yaml`, `src/seacofs_eddy_dataset/` | `f95e42b33a8e30318f119a8bff74f041c1a2a780` |
| Local analysis/workspace | [UNSW-PhD](https://github.com/regdowse/UNSW-PhD), `seacofs_eddy_tilt_analysis/manuscript/` | `ded6021` at context organisation; working-tree changes are separate |

Local workspace: `/home/z5297792/UNSW-PhD`.

- [Manuscript analysis](..): `/home/z5297792/UNSW-PhD/seacofs_eddy_tilt_analysis/manuscript/`.
- Current public-pipeline checkout: `/home/z5297792/ocean-eddy-pipeline/`, clean
  `main` at `f95e42b33a8e30318f119a8bff74f041c1a2a780`. Configuration:
  `/home/z5297792/ocean-eddy-pipeline/config/local.yaml`.
- [Older local pipeline source](../../../seacofs_eddy_dataset_modular/):
  `/home/z5297792/UNSW-PhD/seacofs_eddy_dataset_modular/`.
  Its [configuration](../../../seacofs_eddy_dataset_modular/config/local.yaml) and
  [core code](../../../seacofs_eddy_dataset_modular/src/seacofs_eddy_dataset/core/)
  are available for comparison. This folder is part of the UNSW-PhD workspace,
  not a separate Git checkout of the supplementary manuscript.
- Main-paper editing checkout: `manuscript/text/main-paper/`.
- Supplementary-text editing checkout: `manuscript/text/supplement/`.
  These are independent Git repositories; the paper checkout is clean on `main`
  at `61bc963a442745dc5a91344cad62ad661679f3c2`, fetched and fast-forwarded
  on user request on 2026-11-23. It is a read-only reference under the current
  manuscript workflow; previous review baselines are historical records.
  The supplement checkout is on `main` at
  `3c211760c70a61c11b0876a3cb8c66fd963639dc` (last fetched 2026-11-18).
  The supplement was not updated in the 2026-11-23 paper-only pull.
  `manuscript/text/` is ignored by the
  enclosing analysis repository. Both live on shared storage. Never keep the
  sole editing checkout in node-local `/tmp`.

Review branches: [paper](https://github.com/regdowse/Tilting-of-Mesoscale-Ocean-Eddies-in-The-Tasman-Sea/tree/revision/analysis-refresh-2026-11-17),
[supplement](https://github.com/regdowse/seacofs_eddy_dataset_modular/tree/revision/analysis-refresh-2026-11-17).
See the [checklist](revision/checklist.md) for publication status and commit IDs.
Compiled review PDFs and validation records belong in `context/revision/outputs/`
(ignored). Paper compilation uses XeLaTeX; the supplement uses pdfLaTeX,
with Katana module `texlive/20220321`. The local font configuration adds TeX's
OpenType directory for the paper class's STIX fonts.

The supplementary **text repository** and local **pipeline-code folder** share
the name `seacofs_eddy_dataset_modular`; always identify which one is intended.
Current remote HEAD alone is insufficient to establish historical producing code.

### Producing-stack evidence checked during review

The user confirmed `/home/z5297792/elliptical_streamfunction_parametrisation` was
the intended ESP dependency. The current public-checkout configuration points there.
More directly, `/home/z5297792/ocean-eddy-pipeline/seacofs_eddy.o9801291` records both
`doppio` and `out_core_param_fit` loading from that directory's `functions.py`;
the job reports exit status 0. The PBS script prepends the checkout's `src/` and
checks both imported function paths before running the pipeline.

A separate read-only import-path check resolves both functions to the same file.
Current config SHA-256: `678bdb34df2f478503310c2fbb55698bf566617a302c5d70ee68ce8e3fd98efe`.
Job-log SHA-256: `59db69e8ed3fb228df7baa47d079e392dcb242ec5ea25d2d594e1f9afbede72b`.

The ESP file's current SHA-256 is `1c6b574f1084e9ee4eadfd01553464ff6fc5ae514e60180e7e804d09f212da9d`;
Git blob SHA `2c75ec9466c43e2c277aafedab60a9b5d14470a7` matches
[published ESP v1.1](https://github.com/regdowse/elliptical_streamfunction_parametrisation/tree/v1.1).
The public pipeline README names release 0.1.1 and ESP v1.1.

Nineteen shared Python files match the older workspace copy byte-for-byte.
Differences include the ESP loader (`esp_root` versus `esp_zonodo`), config,
package/release metadata, README and runner notebook. The environment also contains
stale editable-package metadata pointing at an older MRes path; that metadata is
not proof the job imported the older code, since the PBS script controls `PYTHONPATH`.

The processed Parquet has no embedded producing commit/config record. The log
establishes the dependency path, not a complete historical source/config hash.
Use this evidence with that limitation; do not assert undocumented historical
identities. Versioned current code/settings are sufficient to audit the described
operations, while any final provenance statement should retain the qualification.

## Analysis and writing evidence

All paths in this table resolve to local topic folders. Each has a notebook,
JSON settings, Python code and an `outputs/` folder containing numerical tables,
`report.md` and a provenance manifest. Final PNG figures are centralised in
`manuscript/figures/` (Git-ignored); topic runs no longer export figure PDFs.
Manuscript/latexdiff PDFs and the text repositories retain their separate locations.

| Figures | Topic | Writing report |
|---|---|---|
| 2 | [Surface PV](../pv/) | [report](../pv/outputs/report.md) |
| 3 | [Tilt method](../tilt_method/) | [report](../tilt_method/outputs/report.md) |
| 4–6 | [Census](../census/) | [report](../census/outputs/report.md) |
| 7–10, 12 | [PV dominance](../pv_dominance/) | [report](../pv_dominance/outputs/report.md) |
| 11 | [Background flow](../background/) | [report](../background/outputs/report.md) |
| 13–16 | [Controls](../controls/) | [report](../controls/outputs/report.md) |
| 17 | [Composites](../composites/) | [report](../composites/outputs/report.md) |

Figure numbers above identify the analysis workflows and original paper. The
revision introduces Delta before PV, swapping their automatic paper numbers
(Delta becomes Figure 2; PV becomes Figure 3). Stable figure labels and filenames
are unchanged; Figures 4–17 retain their sequence.

The seven notebooks were refreshed on 2026-11-17 after the source dataset update.
The output manifests record source/config/code hashes and run revisions; those
identities, not the current workspace HEAD alone, establish what each run used.
Outputs are Git-ignored local artifacts. Figures 1 and the extra post-composite
figure were excluded from recreation as agreed.

Source dataset: `/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular/`;
processed surface table, saved tilt and confirmed profiles are configured in the
topic JSON files. Independent derived products:

- PV: `/srv/scratch/z5297792/eddy_tilt_manuscript/surface_pv/surface_fit_v1/`.
- Background: `/srv/scratch/z5297792/eddy_tilt_manuscript/background_flow/monthly_v1/`.
- Model archive: `/srv/scratch/z3533156/26year_BRAN2020/`.

## Supervisor feedback and backups

Original export: [2026-11-17T02-58-08Z](reviews/overleaf-comments/Mesoscale%20Eddy%20Vertical%20Tilt%20in%20the%20East%20Australian%20Current-%20Climatology%20and%20Polarity-Depe/2026-11-17T02-58-08Z/).
Its Markdown and `comments.json` contain 14 open threads and eight tracked
changes; the export reports zero stale anchors. All 14 selected passages were
verified to occur exactly once in the checked paper source. The eight tracked
suggestions are reflected in current source wording and accepted by the user for
the revision; their Overleaf review state is not recoverable from plain Git text.
Retain the original export files unchanged; track
decisions separately in [revision/checklist.md](revision/checklist.md).

The user reports saving the Overleaf version to GitHub and downloading the original
source ZIP. The checked paper revision above is a recovery reference; its alignment
with the feedback export and downloaded ZIP remains to be verified. The ZIP has
not been located in `context/`; see [backups status](backups/README.md).
