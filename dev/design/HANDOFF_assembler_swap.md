# Handoff: folding the assembly fallback into WF1

Written 2026-07-30, mid-task. Repo `/home/dmacguig/Documents/GitHub/RiboPilot`,
remote `dmacguigan/RiboPilot`. This picks up where a paused session stopped.
Read `HANDOFF_phase3_complete.md` first for the wider project state; it is still
accurate except where noted at the bottom.

## Stop and read this first

**The work is functionally complete but NOT verified.** Tasks 1 through 6 of the
plan are done. Task 7 (whole-branch verification) was interrupted partway. The
full suite has not been run since task 1. `devtools::check()` has not been run at
all on this branch. Do not open a PR or claim completion until section
"What to do next" is finished.

## Git state, exactly

Two branches are in play. They are independent and can merge in either order.

**Branch `design-debt-stitching`** - pushed, **PR 9 is OPEN** at
`https://github.com/dmacguigan/RiboPilot/pull/9`. Two commits: deleted
`compute_feature_alignment()`, and gave `07_SCAFFOLD_SELECTION.md` a status
section recording that stitching was dropped by `b0f1b12`. That branch is
finished and verified (979 pass, check clean bar a transient clock NOTE). Leave
it alone; it just needs review and merge.

**Branch `assembler-swap`** - the current branch, cut from `main`, NOT pushed.

```
15a8e48 feat(db): schema v16, per-sample assembler column
f2769da docs(design): spec and plan for folding the fallback into WF1
```

Everything else is **uncommitted in the working tree**:

```
 M R/assemble_fallback.R
 M R/nf_command.R
 M R/run_fallback.R
 M dev/design/02_ASSEMBLY.md
 M dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md
 M inst/nextflow/main.nf
D  inst/nextflow/main_fallback.nf     (staged deletion)
 M inst/nextflow/modules/pipeline_r.nf
 M tests/testthat/test-nextflow.R
 M vignettes/RiboPilot.Rmd
?? tests/testthat/test-run-fallback.R  (new, untracked)
```

Commit these before doing anything else, so a failed verification is
distinguishable from a dirty tree.

## What this work is, and why

Spec: `dev/design/14_ASSEMBLER_SWAP.md`. Plan:
`dev/design/plans/2026-07-30-assembler-swap.md`. Both are committed; read the
spec, it is the authority and this handoff does not repeat it.

One paragraph of context. Doc 08 open question 1 asked whether the fallback
workflow was still supported. Investigation found `main_fallback.nf` was a copy
of the pre-pivot `main.nf` tail that had silently drifted out of step: it ran
assemble, score, stitch, BLAST and coverage, and never ran annotation,
screening or triage at all. It is documented in the vignette as Step 5, so users
following the docs were running it. It also stranded samples at
`assemble_switch = 5` with no reset on failure, invisible to every ordinary
query. The maintainer chose to fold it into `main.nf` as a per-sample assembler
swap rather than repair the copy, on the grounds that the duplication was the
defect and repairing it would preserve the condition that produced it.

## What is done

**Task 1 - schema v16 (COMMITTED, verified).** `assemble_opts.assembler` column,
default `getorganelle`. `RIBOPILOT_SCHEMA_VERSION` is now `16L` at
`R/migrate_db.R:8`. The migration also does the project's only value rewrite:
`UPDATE assemble SET assemble_switch = 1 WHERE assemble_switch = 5`, returning
stranded fallback samples to the ordinary queue. Three tests added to
`test-migrate.R`. One pre-existing test at what is now `test-migrate.R:289`
hardcoded `15L` as the post-migration version and was changed to track
`RIBOPILOT_SCHEMA_VERSION`. Both mutations bite (remove the column add: 3
failures; neuter the switch repair: 1 failure). Suite was 991 pass at this
commit.

**Task 2 - `query_samples` returns the assembler (UNCOMMITTED).**
`inst/nextflow/modules/pipeline_r.nf`: added `COALESCE(ao.assembler,
'getorganelle') AS assembler` and `assembler` to the hand-written header printf.
Both had to change together; the header is hand-written because `sqlite3
-header` emits nothing for zero rows.

**Task 3 - branch and merge in `main.nf` (UNCOMMITTED).** Includes
`fallback_assemble`, extends `asm_opts` with `r.assembler`, branches `asm_in`
into `go`/`fb`, and feeds `score()` from
`assemble.out.paths.mix(fallback_assemble.out.paths)`. The `fb` arm is `true`,
not an equality test, because Nextflow's `branch` silently DROPS an item
matching no condition. Header comment rewritten.

**Task 4 - deleted the duplicate (UNCOMMITTED).** `main_fallback.nf` deleted,
`query_fallback_samples` removed from `pipeline_r.nf`, and `R/nf_command.R` lost
its `"fallback"` arm, its `method` parameter and the `--method` line.

**Task 5 - `run_fallback()` reworked (UNCOMMITTED).** Same exported signature, no
breaking API change. It now repoints candidates at a derived `assemble_opts`
preset named `<source>_<method>` carrying the requested assembler, sets
`assemble_switch = 1`, and calls `run_assembly()`. New internal helper
`set_sample_assembler()`. Deriving rather than editing in place is load-bearing:
presets are shared, so editing would change the assembler for every other sample
using it. New file `tests/testthat/test-run-fallback.R`, 17 expectations.
Mutation (make `derived <- src`, i.e. edit in place) bites across 4 assertions.

**Task 6 - documentation (UNCOMMITTED).** Doc 08 open question 1 marked ANSWERED
with the two corrections found (the coverage arity bug was minor next to the
missing annotate/screen/triage half; and `assemble_fallback.R` does not write
`assemblies` at all, contrary to doc 08's claim - it writes FASTAs whose headers
were already the correct per-scaffold shape). Doc 08 open question 2 marked
NARROWED, not answered, with a note that `path_scores` should not be swept in
with the stitch bundle because it is still live. `02_ASSEMBLY.md` records that B2
is not automatic despite the design saying so, and that its baiting is a single
bowtie2 pass rather than the MITObim-style iterate loop. Vignette Step 5 reworded
so it no longer implies a separate workflow.

Verified along the way: `NAMESPACE` is unchanged (correct, no exports moved), and
every changed file is ASCII clean.

## What to do next, in order

1. **Commit the working tree.** Suggested split: one commit for the pipeline
   change (tasks 2-4, `pipeline_r.nf` + `main.nf` + deletion + `nf_command.R` +
   `test-nextflow.R`), one for `run_fallback()` (task 5), one for docs (task 6).
2. **Run the full suite.** `Rscript -e 'devtools::test()'`. Expect 991 plus the
   new nextflow and run-fallback tests, 0 failures, 0 warnings. This is the step
   that was interrupted, and it has NOT been run since task 1. The changes to
   `R/nf_command.R` are the likeliest source of a surprise, because the app's
   Update modal calls it and its signature changed.
3. **Run `devtools::check(document = FALSE, vignettes = FALSE)`.** Expect 0
   errors, 0 warnings, 2 NOTEs (crayon/dbplyr imports, `new_db.Rd` line width).
   A third NOTE reading "unable to verify current time" is a clock artifact and
   is not caused by this work. Watch specifically for a NOTE about
   `R/run_fallback.R` documentation, since its roxygen block changed.
4. **Check whether any Rd needs regenerating.** `run_fallback()`'s roxygen text
   changed but its signature did not. If `man/run_fallback.Rd` is stale, run
   `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` ONLY. Never
   `devtools::document()` - `NAMESPACE` is hand-maintained.
5. **Independent whole-branch review** before the PR, per this project's
   convention.
6. **Open the PR** against `main`.

## Things that will bite you

**No end-to-end run has happened.** The suite reads the workflow files as text;
it does not execute Nextflow. Nobody has ever run the fallback assembler
successfully, before or after this change. `bait_reads()` and `run_spades()`
shell out to `bowtie2` and `spades.py` and both degrade to SILENT NO-OPS when
the tools are absent (`R/assemble_fallback.R:22-24,47-49`), producing an
`_assembly_0` sentinel that `load_assemblies()` records as a failed assembly. So
a green suite does not mean the fallback assembles anything. Spec section 8 says
this plainly; do not overclaim in the PR.

**The `it[7]` index in `main.nf` is positional.** If anyone adds a column to
`asm_opts` before `assembler`, the branch silently misroutes every sample. The
test pins the literal `it[7] == 'getorganelle'`, which catches an edit to the
branch but not an edit to the tuple that shifts the index under it. Worth
strengthening if you touch that area.

**The stitch machinery is now unreferenced but still present.** `stitch_scaffolds`,
`stitch_one`, `stitch_selected`, `order_scaffolds_by_content`, `rp_stitch.R`,
the `stitch` process in `pipeline_r.nf`, and the `stitched` table all survive.
Two are exported public API. This was deliberate (spec section 7): doc 08 open
question 2 says decide the bundle together, and deleting exports is a breaking
change. Do not delete them as a drive-by.

**`path_scores` is NOT part of that bundle.** It is still written by the live
`score` process every WF1 run and read by the app. Only its stitch-side reader
died. Doc 08 is amended to say so.

**Doc 07 needs one sentence changed once PR 9 merges.** The status section on
`design-debt-stitching` says the only remaining stitch caller chain runs through
`main_fallback.nf`. After this branch merges that is no longer true. This
reconciliation cannot be done from `assembler-swap` because the text is not on
it. Spec section 10 and the plan's deferred list both record this. Do it in a
follow-up once both branches land.

**Migration ordering.** If someone merges a different schema change first, v16
collides. `RIBOPILOT_SCHEMA_VERSION` appears in exactly one place,
`R/migrate_db.R:8`. Nothing else on `main` is currently touching the schema, so
this is only a risk if another branch appears.

## Conventions this project enforces

Unchanged from `HANDOFF_phase3_complete.md`, but the ones this work touches:

- `NAMESPACE` is hand-maintained. Never `devtools::document()`. Regenerate Rd
  with `roxygen2::roxygenise(roclets = "rd")` only.
- ASCII only, everywhere. No em dashes, no smart quotes.
- Migrations are additive and idempotent. The v16 `assemble_switch` rewrite is
  the sole exception in the project's history and is commented as such at
  `R/migrate_db.R:441-450`; do not treat it as a precedent.
- Minimal comments, explaining why rather than what changed.
- Mutation testing is the practice to carry forward: revert the fix, confirm a
  test fails, restore, confirm the tree is clean. Every fix in this branch has
  been mutation-checked; the results are recorded above. **Check the mutation
  actually applied** - a `sed`/`perl` that matches nothing looks exactly like a
  test that correctly survived.
- Never push without being asked. `design-debt-stitching` was pushed on explicit
  instruction; `assembler-swap` has NOT been pushed and should not be without
  asking.

## Amendment to HANDOFF_phase3_complete.md

That document says "Nothing is unblocked in the contamination plan" and lists
design debt as available non-contamination work. Both branches here came out of
that list and it is otherwise still accurate. Two of its specific claims turned
out to be imprecise and are corrected in the docs themselves:

- It says doc 08 lists `stitch_scaffolds` as deprecable and "that is wrong".
  Doc 08 actually says "deprecate rather than delete". The real problem is that
  deprecating it was entangled with the unanswered fallback question, which
  this branch answers.
- It repeats doc 08's claim that `assemble_fallback.R` writes `assemblies`
  through an untraversed path. It writes no table at all.

The invertebrate clade data gathering remains the single blocker for the most
items (validation gate clade half, `min_isomer_anchor_pid`, the `human` screen
subclass) and is untouched by any of this.
