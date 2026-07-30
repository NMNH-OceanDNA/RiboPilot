# Handoff: current design state

Written 2026-07-30. Repo `/home/dmacguig/Documents/GitHub/RiboPilot`, remote
`dmacguigan/RiboPilot`. This is the single current handoff; it supersedes and
replaces `HANDOFF_phase3.md`, `HANDOFF_phase3_complete.md` and
`HANDOFF_assembler_swap.md`, all deleted in the same commit. Everything from
those still worth keeping is carried below.

## Read this first

**PR #10 is open and unmerged**: "Fold the assembly fallback into WF1 as a
per-sample assembler", branch `assembler-swap`. It is verified (1122 tests,
`check` clean) but has NOT had an end-to-end run complete. Do not start design
work that assumes it merged, and do not branch new work off `assembler-swap`
without knowing that.

`main` is at the merge of PR #9. The branch was cut before that, so it is behind
by one merge. No conflicts.

## Where the code stands

Schema v17. Phases 0 through 4 of `05_ROADMAP.md` are built. The contamination
and multi-unit programme of doc 08 is implemented through sub-project 3c plus
the validation fixtures; 3d and 3e are blocked (below).

The fallback assembler is now a per-sample choice inside WF1
(`assemble_opts.assembler` = getorganelle | bait | denovo), not a separate
workflow. `main_fallback.nf` is deleted.

Test project: `~/Documents/DJMCompBio/RiboPilot_test/full_project_docker`. It has
been run through assembly on all three assemblers. Its WF1 tail past remote
BLAST has NOT completed for the fallback samples.

## The design documents, and what each is for

`00_OVERVIEW.md`, `01_ARCHITECTURE.md`, `02_ASSEMBLY.md`,
`03_ANNOTATION_CURATION.md`, `04_EXISTING_TOOLS.md` are the original design.
`05_ROADMAP.md` holds the phase plan and the eight project-level open questions.
`06_TEST_DATA.md` is the test corpus.

`07_SCAFFOLD_SELECTION.md` describes the stitching design that was later dropped;
it now opens with a status section saying so. **That status section exists only
on `main`**, not on `assembler-swap`.

`08_MULTIUNIT_UX_AND_CONTAMINATION.md` is the big one (82 KB) and the centre of
gravity for current design. Section 1 is bugs, 2 quick wins, 3 the contamination
scheme, 4 the multi-unit UX, 5 larger and optional, 6 rejected, 7 open questions.

`09` through `13` cover the merged sub-projects (triage rules, hit
classification, screen verdict, intragenomic variants, validation fixtures).
`14_ASSEMBLER_SWAP.md` covers PR #10.

Plans live in `dev/design/plans/YYYY-MM-DD-topic.md`.

**Read each document's amendments, not just its body.** In three of the five
merged sub-projects the worst defects were defects in the spec, and the
corrections are recorded in place.

## Design documents that now lag the code

Created by the last day of work. All verified 2026-07-30. This is the first
thing to fix, because every one of them will mislead a reader who trusts it.

1. **Doc 08 section 1.9 still lists `assemble_opts.max_scaffolds` as "stored and
   editable, never enforced"** (line 229). It is enforced now. Both it and
   `max_paths` are gates: a sample over either limit is failed with a note and
   its scaffolds parked, and the sequences are kept so raising the limit and
   re-running recovers it without re-assembling. Note the doc argued AGAINST a
   flat count ("degraded museum extracts legitimately fragment one repeat across
   many short contigs") and the maintainer chose the flat count anyway, with the
   reasoning "user can always go back and rerun those later". Record the decision
   and the objection it overrode.
2. **Doc 08 section 3 (line 312) reasons from "`max_scaffolds` unenforced"** when
   sizing the scaffold table. Same fix.
3. **`min_assembly_length` is still unenforced and is not recorded anywhere as
   debt.** It sits in the schema, in `new_db()` and in the options modal, and
   nothing reads it. It should either join the gate or be listed in doc 08
   section 1.9 next to where `max_scaffolds` used to be.
4. **Doc 01's `assemble_opts` field list is stale.** It names `k_list`,
   `min_length` and `expected_length`, which do not exist, and does not name
   `bowtie2`, `spades` or `min_assembly_length`, which do. Its `assembler` value
   list is also wrong: the real values are getorganelle, bait, denovo.
5. **The sample-state semantics changed and no document says so.** Assembly
   loading writes state 4 (in progress); triage, the last per-sample step, writes
   2 (done); a sample also settles when its remaining units are held or ignored.
   Consequence worth writing down: samples at 4 are picked up by the next run, so
   an interrupted run resumes, and `nf_pending_count()` reports unfinished
   samples as pending rather than zero.
6. **Doc 08 section 4.4 ("the Assemble modal becomes the evidence surface") is
   now behind the implementation.** The modal gained row-click-to-examine with a
   highlight, tick boxes, and a Copy FASTA button. The old per-row "hits" link is
   gone.
7. **Doc 07's sentence about the last stitch caller.** On `main` it says the only
   remaining stitch caller chain runs through `main_fallback.nf`. PR #10 deletes
   that file. Fix once #10 merges; it cannot be done from the branch because the
   text is not on it.

## Open design questions, and what each actually needs

From doc 08 section 7, with current status:

- **Q1 fallback workflow: ANSWERED** by doc 14 / PR #10.
- **Q2 does `path_scores` survive: NARROWED, still open.** After PR #10 nothing
  in any shipped workflow calls the stitch machinery, so the bundle
  (`stitch_scaffolds`, `stitch_one`, `stitch_selected`,
  `order_scaffolds_by_content`, `rp_stitch.R`, the `stitch` process, the
  `stitched` table) is unreferenced rather than half-referenced. It is still not
  deleted because two of those are exported public API. `path_scores` is NOT part
  of the bundle: it is written by the live `score` process every run and read by
  the app. Deciding it is a question about the app's Selected Path and Best Score
  columns.
- **Q3 how aggressive should auto-triage be by default:** needs a maintainer
  call, informed by doc 13's measurement (zero false positives on 19 clean host
  units, recall 4/4 on catchable spike-ins, fish only).
- **Q4 is `_p<path>_s<scaffold>` the right GenBank SeqID suffix:** unresolved, and
  it needs a migration story for projects that already exported bare IDs.
- **Q5 should the envelope numbers ship before validation:** unresolved.
- **Q6 run/lane/plate metadata:** blocks 3e, see below.
- **Q7 should `paralogous_variant` block export, warn, or inform:** submission
  policy call.

`05_ROADMAP.md` also carries eight older project-level questions. Numbers 1
(taxonomic scope for v1), 2 (circular vs linear output) and 8 (5S in or out) are
still genuinely open and are the ones that would shape Phase 5.

## What is blocked, and on precisely what

**Sub-project 3d, the positional and chimeric channel**, is blocked on the
identity envelope. This was deferred after measurement, not out of laziness: the
reference corpus holds 1172 hits at class distance, 2 at phylum and none beyond,
so there is no data in the band where an envelope would act, and genuine host 18S
reaches 99.6 percent identity at class distance. Unblocking needs identity data
at phylum distance and beyond. The local-BLAST machinery for that already exists
in `dev/make_screen_fixture.R`.

**Sub-project 3e, cross-sample screening and index hopping**, is blocked on
`samples.batch`, which does not exist (zero references in `R/init_db.R` and
`R/migrate_db.R`). Adding the column is trivial; the detection is meaningless
until users supply run, lane or plate metadata.

**The invertebrate clade data** is still the single blocker for the most items:
the clade half of doc 08 section 5's validation gate, `min_isomer_anchor_pid`,
and the `human` subclass of `screen_class_for()`. Verified 2026-07-30: the three
read sets are not downloaded, their clade reference DBs are not built
(`ref_dbs/getOrganelle/seeds/` holds only actinopterygii), and no invertebrate
project has ever been run. That is read downloads, three reference DBs, three
assemblies and three remote BLAST runs before any labelling starts. Long wall
clock, not a coding task. Discuss cost with the maintainer first.

**Auto-drop and the export gate** remain disabled by choice, not by blocker.

## Unblocked work available now

- The seven stale-document fixes above.
- **Per-position minor-allele fraction** (doc 08 section 5). A submitted rDNA
  scaffold is a consensus over many non-identical array copies, and that is the
  first question a systematist asks of a genome-skim submission.
  `coverage.depths` already stores the per-position track. Check whether the
  mapped reads are still available before promising it; the depth track alone may
  not be enough.
- **The corpus gaps in doc 13 section 5.1**: three named cases the fixture cannot
  catch. The congener one needs a host unit whose hits are dominated by same-genus
  relatives, which none of the ten existing samples provides.
- **Deciding the stitch bundle** (Q2). It is now a clean question.

## A latent defect found while testing, not yet fixed or recorded

`inst/nextflow/bin/rp_env.R` says "RIBOPILOT_SRC WINS over an installed package"
and exists to stop a stale install silently shadowing development code. It does
not currently work: `nf_command()` passes `--ribopilot_src <pkg>/inst`
(`R/nf_command.R:27`, via `app_sys()`), and `rp_env.R` only loads from source if
`<src>/DESCRIPTION` exists. In a development checkout `inst/` has no
`DESCRIPTION`, so the check fails and every pipeline R process silently uses the
container's installed RiboPilot instead.

Verified directly: `app_sys()` returns `.../RiboPilot/inst`,
`file.exists(".../inst/DESCRIPTION")` is FALSE, and its parent is TRUE.

This is exactly the trap the comment describes. It also means editing `R/` while
a pipeline runs is safe, and editing `inst/nextflow/bin/` is NOT, because
Nextflow stages `bin/` from the working tree per task. Anyone relying on
`RIBOPILOT_SRC` to test uncommitted package code against Docker is not testing
what they think they are testing.

## Conventions this codebase enforces

`NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd
with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only, and add export
lines by hand.

ASCII only, everywhere: R, tests, fixtures, Nextflow, documentation. No em
dashes, no smart quotes.

Comments are minimal and explain why, not what changed. Do not narrate a bugfix.

Migrations are additive and idempotent. The v16 `assemble_switch = 5` repair is
the sole value rewrite in the project's history and is commented as the exception
it is. Shared DDL constants live in `R/migrate_db.R` and are referenced by both
`R/init_db.R` and the migration, so a fresh database and a migrated one cannot
diverge. There is a test that enforces exactly that for `assemble_opts`.

`migrate_db()` now runs inside one write transaction. It has to: every pipeline
process migrates on open and annotation fans out per unit, so a dozen processes
attempt a pending migration at once, and the unserialised version killed a real
run with "duplicate column name". The test that reproduces it uses 12 forked
processes and is skipped under `R CMD check`, which caps a package at 2 cores.

The infix `%|NA|%` is **scalar only and errors on a vector**. It is defined TWICE
(`R/utils.R` and `R/app_opts.R`) and the `utils.R` one wins at load order; that
one has no NULL guard. `%||%` tests length, not NA. This has caused defects in
five separate tasks now.

`dev/` is in `.Rbuildignore`. Dev-only scripts belong there and never become
package code.

Never push without being asked. Branch from `main` and ask first for anything
sizeable.

## How the work is run

Spec (`dev/design/NN_TOPIC.md`), then implementation plan
(`dev/design/plans/YYYY-MM-DD-topic.md`), then task-by-task execution with an
independent review after each task and a whole-branch review before the PR.

**Mutation testing is the single practice to carry forward.** Revert the fix,
confirm the test fails, restore, confirm the tree is clean. It keeps finding
tests that pass while pinning nothing. In this session it caught a test that
claimed to check help-text placement in the assembler modal but would have passed
with the help text in the wrong place entirely.

Two lessons worth internalising:

- **A mutation that does not bite is information, not a nuisance.**
- **Check the mutation actually applied.** A `sed`/`perl` that matches nothing
  looks exactly like a test that correctly survived.

**Treat every claim as a hypothesis, including the design docs' own.** Confirmed
wrong so far: doc 11 named the wrong blind spot; doc 08 asked for a human spike-in
as though it were cross-kingdom; doc 08 called `assemble_fallback.R` a writer of
the `assemblies` table when it writes no table at all; doc 12 sited a params
sublist that had never been built.

**Verify with a real project, not only the suite.** The suite reads the Nextflow
files as text and never executes them. This session's two worst bugs, the
migration race and the fallback querying columns before migrating, were both
invisible to 1000+ passing tests and appeared in the first real run.

## Known open items, recorded rather than forgotten

- `samples.reference_support` is written on every `screen_sample()` run and has
  no reader. Doc 11 section 3 wants UI for it when it is 4 or worse.
- The screen's confidence weights are priors, not measurements. Any change must
  bump `SCREEN_VERSION`.
- `min_isomer_anchor_pid` at 99.5 has never been exercised as a discriminator.
- The intergenic spacer is never annotated, so doc 12's isomer message infers
  spacer copy number rather than measuring it.
- Doc 12 section 8's two salvage items are deferred with measurements.
- The `human` subclass is unreachable end to end from the fixture databases.
- All strands in the test data are `+`, so the all-minus and mixed-strand anchor
  paths in `R/screen_variants.R` have unit-test coverage only.
- The fallback assembler's two external tools degrade to SILENT NO-OPS when
  absent (`R/assemble_fallback.R`), producing an `_assembly_0` sentinel that
  loads as a failed assembly. A green suite does not mean it assembles anything.

## Practical notes that save time

**Remote BLAST is the bottleneck and it fails silently.** NCBI's BLAST compute is
a shared queue that stops answering under load; their file serving keeps working
throughout, which is the distinction to exploit. `blastn -remote` has no timeout,
so always wrap it in `timeout`. Zero CPU time with long elapsed time is the
signature.

**Prefer local BLAST against NCBI's prebuilt rRNA databases** where fidelity
allows; `dev/make_screen_fixture.R` shows the pattern. Do NOT use
`update_blastdb.pl` from the conda package: it cannot find curl, leaves `$curl`
uninitialised at its line 599, tries to exec `--user` as a program, and reports
success anyway.

**The image has RiboPilot installed but not the dev dependencies.** No `pkgload`,
no `withr`. Scripts meant to run inside it must fall back to `library(RiboPilot)`.

**Rebuild the image with `docker/deploy-local.sh test`** after changing package
code, or the pipeline keeps running the old copy. See the latent defect section
above for why that is more true than it looks.
