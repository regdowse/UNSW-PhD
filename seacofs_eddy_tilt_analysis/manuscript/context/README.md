# Manuscript coordination

Current stage: context organised; coordinated read-only review is next. Neither
the main paper nor supplementary manuscript has been edited in this stage.

- [Resources](resources.md): text/code repositories, revision snapshots and analysis outputs.
- [Revision checklist](revision/checklist.md): one checklist for both manuscripts and supervisor feedback.
- [Analysis usage](workflow/README.md): notebooks, configurations, caches and generated reports.
- [Analysis decisions/history](workflow/workflow_review.md), [implementation plan](workflow/implementation_plan.md)
  and [figure source map](workflow/info.md).
- [Supervisor exports](reviews/overleaf-comments/): original files, retained unchanged and Git-ignored.
- [Backups](backups/README.md): original source/PDF/ZIP snapshots and their status.

The supplementary methods should describe the code/configuration that produced
the dataset. Check that provenance before treating the latest repository as the
dataset-producing version. The main paper should match the agreed analysis
definitions and refreshed figures/numerical outputs. Discuss substantive
discrepancies before changing science or interpretation.

Preserve supervisor comments and tracked suggestions; do not silently accept,
reject or resolve them. Keep originals unchanged and record dispositions in the
revision checklist. Propose structural changes before rewriting the manuscripts.

Generated figures, numerical tables, reports and manifests stay with their topic
in `manuscript/<topic>/outputs/`. Use those current reports rather than copying
historical numbers from the workflow notes. No additional progress log is needed.
