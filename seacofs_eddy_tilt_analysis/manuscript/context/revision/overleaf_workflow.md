# Overleaf and local manuscript workflow

Applies separately to the main paper and supplement. Repository locations are in
[resources.md](../resources.md).

## Current agreement (2026-11-20)

The user implements all changes to the main manuscript and supplementary text.
Codex provides analysis, checks, explanations and proposed wording in conversation;
Codex must not directly edit, commit, merge or push manuscript text or assets in
these text repositories. Earlier authorisations to revise/push manuscripts are
historical and do not apply to future work.

## Synchronise on request

1. When the user wants Codex to see updated text, they push Overleaf to GitHub
   and ask Codex to pull locally.
2. Codex checks for local changes, fetches GitHub and updates the relevant checkout
   in `manuscript/text/`, preserving any pending work. Record the fetched commit.
3. Treat the checkout as a read-only reference between requested updates. Do not
   assume it includes later Overleaf edits or fetch routinely without a request.

## Suggestions and comparisons

- Present discrepancies and exact suggested replacements for the user to implement.
- Analysis/code work does not authorise changing manuscript text or copying updated
  figures into the text repositories. Report outputs and affected passages instead.
- The former automatic cumulative-diff workflow is retired. Preserve existing
  baselines and review artifacts as historical records; do not regenerate routinely.
- Provide Git comparisons when requested. Generate text-only latexdiff PDFs only
  on explicit request; do not routinely compile review PDFs.
- Original pre-refresh tags remain recovery points.

## Preserve feedback

Keep exported supervisor comments and original backups unchanged. Record comment
dispositions in [checklist.md](checklist.md); an addressed item does not resolve its
Overleaf thread. Keep routine workflow notes here rather than adding more logs.
