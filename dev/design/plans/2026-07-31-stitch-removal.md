# Plan: remove the stitch bundle, enforce min_assembly_length

Branch `stitch-removal`, cut from `main` (PR #10 merged). Answers doc 08 open
question 2 and stale-document item 3 from `HANDOFF.md`.

## Decisions taken (maintainer, 2026-07-31)

- The `stitched` table is dropped in BOTH a fresh schema and a migration, so the
  two cannot diverge. This is the second value-destroying migration in the
  project's history; the data it removes is stitch-era output that nothing reads.
- `annotate_project()` goes with it. It is exported, has zero callers, stitches
  the auto-selected path, and writes `annotations` rows without `path`/`scaffold`,
  which the per-scaffold schema (v7 onward) requires. It could not have worked.
- `min_assembly_length` is enforced as a FILTER at load time, not as a gate.
  Short junk contigs are noise; they are not a reason to fail a sample the way
  too many paths or scaffolds is.

## Tasks

1. Delete `R/stitch_scaffolds.R` (`order_scaffolds_by_content`,
   `stitch_scaffolds`), `R/annotate_project.R`, `inst/nextflow/bin/rp_stitch.R`,
   `dev/prototype_stitch_viz.R`.
2. Remove `stitch_one()` and `stitch_selected()` from `R/pipeline_steps.R`.
3. ~~Remove `compute_feature_alignment()`.~~ Already gone on `main`; doc 08
   section 5 records it. Checked rather than assumed.
4. Remove the `stitch` process from `inst/nextflow/modules/pipeline_r.nf`.
5. Schema: stop creating `stitched` in `R/init_db.R`, bump the version, and drop
   the table in `R/migrate_db.R`.
6. `NAMESPACE`: drop `stitch_scaffolds`, `stitch_one`, `annotate_project`.
   Delete the matching `man/*.Rd` by hand. Never run `devtools::document()`.
7. Tests: replace the "no workflow invokes stitch" check with one that the
   machinery is gone, and add a fresh-vs-migrated check that `stitched` does not
   exist either way.
8. `min_assembly_length`: filter scaffolds under the limit in `load_assemblies()`,
   before the gate check, with a note when anything is dropped.
9. Docs: doc 07 status section and doc 08 question 2. `HANDOFF.md` is left
   alone on purpose: it is being rewritten on branch `assembler-swap` (PR #11)
   and editing it here would only manufacture a conflict.

Mutation-test every new test: revert the change, confirm the test fails, restore.

## Outcome

Suite 1138 tests, 0 failures. `R CMD check` clean apart from the pre-existing
`crayon`/`dbplyr` unused-Imports NOTE.

Mutation testing, all verified to have actually applied:

- Parking short scaffolds removed -> 3 tests fail. (An earlier, redundant copy of
  the same assignment did NOT bite, which is how it was found and deleted; the
  surviving one runs after the ignore carry-over, where it matters.)
- Gate counting short scaffolds again -> 3 tests fail.
- `DROP TABLE stitched` removed -> the v18 test fails.

## Follow-up, blocked on a merge

Doc 08 section 1.9's `min_assembly_length` row (added on branch `assembler-swap`,
PR #11) still says "never enforced". Once PR #11 merges, update it: enforced as a
load-time filter, parked not discarded, counted out of the gates.
