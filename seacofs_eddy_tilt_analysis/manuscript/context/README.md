# Manuscript coordination

Current stage: coordinated paper and supplementary-text drafts prepared on
`revision/analysis-refresh-2026-11-17` in each text repository. Publication and
validation status are recorded in the revision checklist. The user subsequently
merged the paper revision into `main`; originals are retained under the
`pre-analysis-refresh-2026-11-17` tags and original feedback exports.

- [Resources](resources.md): text/code repositories, revision snapshots and analysis outputs.
- [Revision checklist](revision/checklist.md): one checklist for both manuscripts and supervisor feedback.
- [Analysis usage](workflow/README.md): notebooks, configurations, caches and generated reports.
- [Analysis decisions/history](workflow/workflow_review.md), [implementation plan](workflow/implementation_plan.md)
  and [figure source map](workflow/info.md).
- [Supervisor exports](reviews/overleaf-comments/): original files, retained unchanged and Git-ignored.
- [Backups](backups/README.md): original source/PDF/ZIP snapshots and their status.

## Manuscript workflow

Follow [Overleaf and local manuscript workflow](revision/overleaf_workflow.md).
Push the latest relevant Overleaf source to GitHub before local work. Review uses
cumulative Git diffs against [the pinned baseline](revision/review_baselines.json).
Generate latexdiff PDFs only on explicit request. Keep the detailed operating rules
in that workflow note to avoid duplication.

The supplementary methods should describe the code/configuration that produced
the dataset. Check that provenance before treating the latest repository as the
dataset-producing version. The main paper should match the agreed analysis
definitions and refreshed figures/numerical outputs. Discuss substantive
discrepancies before changing science or interpretation.

The user approved the proposed revision scope and accepted all pending tracked
suggestions. Preserve the original exports and record dispositions in the revision
checklist; an addressed item does not resolve its Overleaf thread. The user will
handle figure presentation changes. Discuss further substantive scientific changes.

Final figure PNGs go to the Git-ignored `manuscript/figures/` folder. Numerical
tables, reports and manifests stay in `manuscript/<topic>/outputs/`. Use those current reports rather than copying
historical numbers from the workflow notes. No additional progress log is needed.
