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
- [Local pipeline source](../../../seacofs_eddy_dataset_modular/):
  `/home/z5297792/UNSW-PhD/seacofs_eddy_dataset_modular/`.
  Its [configuration](../../../seacofs_eddy_dataset_modular/config/local.yaml) and
  [core code](../../../seacofs_eddy_dataset_modular/src/seacofs_eddy_dataset/core/)
  are available for comparison. This folder is part of the UNSW-PhD workspace,
  not a separate Git checkout of the supplementary manuscript.
- The two text repositories have been accessed remotely. No dedicated local
  editing checkout has been prepared in this coordination step.

The supplementary **text repository** and local **pipeline-code folder** share
the name `seacofs_eddy_dataset_modular`; always identify which one is intended.
Confirm local/public pipeline differences and dataset-producing code, configuration
and external ESP version during the audit. Current remote HEAD alone is insufficient.

## Analysis and writing evidence

All paths in this table resolve to local topic folders. Each has a notebook,
JSON settings, Python code and an `outputs/` folder containing figures, numerical
tables, `report.md` and a provenance manifest.

| Figures | Topic | Writing report |
|---|---|---|
| 2 | [Surface PV](../pv/) | [report](../pv/outputs/report.md) |
| 3 | [Tilt method](../tilt_method/) | [report](../tilt_method/outputs/report.md) |
| 4–6 | [Census](../census/) | [report](../census/outputs/report.md) |
| 7–10, 12 | [PV dominance](../pv_dominance/) | [report](../pv_dominance/outputs/report.md) |
| 11 | [Background flow](../background/) | [report](../background/outputs/report.md) |
| 13–16 | [Controls](../controls/) | [report](../controls/outputs/report.md) |
| 17 | [Composites](../composites/) | [report](../composites/outputs/report.md) |

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
changes; the export reports zero stale anchors. Match these to the chosen text
baseline during review. Retain the original export files unchanged; track
decisions separately in [revision/checklist.md](revision/checklist.md).

The user reports saving the Overleaf version to GitHub and downloading the original
source ZIP. The checked paper revision above is a recovery reference; its alignment
with the feedback export and downloaded ZIP remains to be verified. The ZIP has
not been located in `context/`; see [backups status](backups/README.md).
