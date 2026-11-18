# Overleaf and local manuscript workflow

Applies separately to the main paper and supplement. Repository locations are in
[resources.md](../resources.md).

## Synchronise before work

1. User pushes the latest relevant Overleaf source to GitHub and confirms it.
2. Codex checks for local edits, fetches GitHub and updates the local checkout in
   `manuscript/text/`, preserving any pending work.
3. Record the source commit. Unconfirmed Overleaf changes are not assumed present
   in GitHub or the local checkout.

## Edit and review

- Keep changes focused; discuss substantive scientific changes before applying them.
- After manuscript edits, provide a cumulative, readable Git diff against the latest
  confirmed pull, pinned in [review_baselines.json](review_baselines.json).
- Include all pending changes, whether committed or uncommitted, and new manuscript
  source files. Describe figure/binary changes separately.
- Save review diffs under `context/revision/outputs/review/`, identifying the baseline,
  current commit and working-tree status. Link the relevant diff in the handoff.
- Generate a text-only latexdiff PDF **only when explicitly requested**. Do not
  routinely compile PDFs for review; compile when requested or needed to validate
  a change that could affect compilation.
- GitHub updates do not create native Overleaf Reviewing suggestions. Use the
  cumulative diff to review changes before pulling them into Overleaf.

## Advance the baseline

Do not reset it after individual edits, commits, pushes or routine fetches. Advance
it after the user confirms the reviewed changes were pulled into Overleaf and the
next synchronised source has been pulled locally. Retain the previous review diff.
If remote changes arrive while local revisions are pending, reconcile them explicitly.
The original pre-refresh tags remain separate recovery points.

## Preserve feedback

Keep exported supervisor comments and original backups unchanged. Record comment
dispositions in [checklist.md](checklist.md); an addressed item does not resolve its
Overleaf thread. Keep routine workflow notes here rather than adding more logs.
