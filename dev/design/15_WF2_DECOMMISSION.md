# 15 - Decommission WF2: one workflow, one entry point

Written 2026-07-31. Branch `wf2-decommission`, cut from `main`.

## What this changes

WF1 already absorbed the first annotation pass: it annotates every scaffold as
it is assembled, gating on nothing, because nothing has been reviewed yet
(`query_sample_units`, `modules/pipeline_r.nf:59`). WF2 (`main_annotate.nf`)
survives only as the RE-annotation path for units queued again after review:
`annotate_switch IN (1,4) AND annotate_lock = 0 AND assemble_lock = 1`.

That leaves two entry points for one pipeline, and a user has to know which of
their edits needs which. This sub-project folds the re-annotation path into WF1
and deletes WF2.

## Decisions taken (maintainer, 2026-07-31)

1. **WF1 gains an annotate-only branch.** It selects two kinds of work: samples
   needing assembly, and units needing only annotation. Annotate-only units skip
   preprocess, assemble, score, the remote BLAST and coverage, and go straight to
   `annotate_tools -> reconcile -> blast_ref_align -> curate`.
2. **Annotation and curation options move to the Assemble surface, UI only.**
   The `annotate_opts` and `curate_opts` tables, their per-unit links and their
   modals are untouched. No schema change, no migration.
3. **`main_annotate.nf` and `run_annotate()` are both deleted.** Breaking change,
   consistent with how the stitch exports were handled.
4. **One Update button, on the Assemble tab.** It runs whatever is queued. The
   Annotate tab becomes purely a review surface: edit the feature model, trim,
   restore, mark, lock.

## Why Nextflow's `-resume` cannot do this instead

Asked and answered before choosing, because it is the obvious first idea.

`-resume` is already on: `nf_command()` passes it every run. `preprocess`,
`assemble` and `fallback_assemble` carry no `cache false`, so a second run
already skips re-assembly. That half is free and stays.

It cannot select the work, for two reasons:

- **Thirteen processes carry `cache false`**, including every one that touches
  the database. That is deliberate: the SQLite database is the source of truth
  and mutates between runs, so a cached task would replay against state that has
  moved.
- **`annotate_tools` reads its options from the database inside the script**
  (`modules/pipeline_r.nf:120-124`, an in-script `sqlite3` query for tool,
  kingdom and ITSx taxon). Those never enter the task's input tuple, so changing
  a ruleset changes no task hash. If annotation were cacheable, `-resume` would
  return the stale annotation as a cache hit.

The distinction: caching answers "has this task run with these inputs"; the
question here is "what work is outstanding", which lives in `annotate_switch`
and `assemble_switch`. Only a query answers that. Threading every option into
task inputs would make caching correct in principle and still would not notice a
restored scaffold or a hand-edited sequence.

## The work selection

Two queries, one workflow.

`query_samples` is unchanged: samples with `assemble_switch IN (1,4)` and
`assemble_lock = 0`.

`query_annotate_units` replaces `query_annotate_samples`. Units with
`annotate_switch IN (1,4)`, `annotate_lock = 0`, on a non-ignored scaffold,
**excluding any sample that this same run is going to assemble**. That exclusion
is the one piece of new logic and it is load-bearing: assembly re-seeds a
sample's annotate units itself, so without it a sample queued for both would be
annotated twice in one run, racing on the same rows.

The `assemble_lock = 1` gate that WF2 used is DROPPED. It existed so WF2 could
not run before the user had reviewed the assembly. With one workflow there is
nothing to run early: a unit is annotate-only work because a user queued it.

## Wiring

`annotate_tools` takes `(ID, path, scaffold, cpus, mem)` from either source, so
the two arms `mix()` into the existing chain and everything downstream is
unchanged. Two existing pieces of defensive wiring already cover the new arm:

- `blast_ref_align` waits on a `ref_done` sentinel that is `ifEmpty('ready')`,
  so an annotate-only run with no reference fetch does not stall.
- `triage` joins curation against the screen with `remainder: true`, so a sample
  that curated but never screened in this run still reaches triage.

Neither needed changing, but both are why the arm is cheap. Note them in tests
so a later simplification does not remove them.

## What is deleted

- `inst/nextflow/main_annotate.nf`
- the `query_annotate_samples` process
- `run_annotate()`, its Rd and its NAMESPACE line
- `nf_command()`'s `workflow` argument (one workflow left to build a command for)
- the Annotate tab's Update button and the run modal's annotate mode

`nf_pending_count()` stops being per-workflow and counts both kinds of
outstanding work, since one button runs both.

## Open assumption, easily changed

The maintainer said "move the assemble and curate opts into the Assemble modal".
Options today are set from link cells in the table's Options column group, not
from inside `assemble_modal()`, so this is implemented as: the Assemble tab's
Options group gains `annotate_opts` and `curate_opts`. If the intent was buttons
inside the scaffold-review modal itself, that is a small follow-up.

## Risks

**A sample queued for both kinds of work.** Covered by the exclusion above, and
it needs a test, because the failure is a race rather than an error.

**Users with `run_annotate()` in scripts.** It disappears. Nothing in the repo
calls it outside its own definition and the docs.

**The claim that `-resume` skips re-assembly is read off the process
definitions, not measured.** Verify it with a real run before relying on it in
the docs.

## Status 2026-07-31: BUILT and verified against a real project

Suite 1191 tests, 0 failures. `R CMD check` clean apart from the pre-existing
unused-Imports NOTE.

**Verified end to end.** On the three-sample baited test project with every
sample already assembled, one unit was queued for annotation
(`annotate_switch = 1`) and WF1 run. Seven tasks: `query_samples` (which found
nothing to assemble), `query_annotate_units`, then `annotate_tools`,
`reconcile`, `blast_ref_align`, `curate`, `triage` for that one unit.
`preprocess`, `assemble`, `fallback_assemble`, `score`, `blast_genbank`,
`blast_ref_fetch`, `coverage` and `query_sample_units` all sat at zero. The unit
returned to `annotate_switch = 2` with its five features, the sample settled at
2, and pending work went to 0. No NCBI round trip.

**A defect the suite missed and the run caught.** `run_assembly()` still passed
`"assemble"` positionally to `nf_command()`. With the workflow parameter gone it
bound to `resume`, and the driver died at launch. Nothing in the suite called
the driver with a project, so nothing noticed. There is now a test that checks
every `nf_command()` call site against the formals and allows only the project
path to be positional; it fails if the second positional argument comes back.

Mutation testing, all verified to have applied: dropping the assemble-queue
exclusion, dropping the `ignore` filter, restoring the Annotate tab's Update
button, restoring the annotation option columns to the Annotate query, and
making `nf_pending_count()` forget annotate-only work each break a test. One
earlier version of the call-site test did NOT bite (it allowed a second
positional argument, which is the whole bug) and was tightened until it did.

**Still unverified: the claim that `-resume` skips re-assembly.** It is read off
the process definitions, not measured. The run above proves nothing about it,
because assembly was skipped by the query rather than by the cache. Measuring it
means queueing a real assembly, which also re-runs the remote BLAST.
