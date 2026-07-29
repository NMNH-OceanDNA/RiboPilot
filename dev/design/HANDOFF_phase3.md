# Handoff: RiboPilot contamination work, Phase 3

Written 2026-07-29. Repo is `/home/dmacguig/Documents/GitHub/RiboPilot`, GitHub
remote `dmacguigan/RiboPilot`.

## Where the code stands

The working branch is `contamination-phase3b`, cut from `main` at `234e97c`. It
is complete, reviewed, and **unmerged and unpushed**, because the maintainer is
testing it by hand first. The full test suite is at 705 passing with 0 failures,
`devtools::check(document = FALSE, vignettes = FALSE)` reports 0 errors, 0
warnings and 2 NOTEs (both pre-existing: crayon/dbplyr imports, and `new_db.Rd`
line width), and `nextflow lint inst/nextflow/main.nf` reports 0 errors and 6
warnings, unchanged from the branch point. A Docker image `ribopilot:test` has
been rebuilt from this branch and verified to carry schema version 14, screen
scoring version 2, and all seven triage rules.

Two earlier sub-projects are already merged into `main`. Pull request 4 delivered
a rewrite of the triage rules, and pull request 5 delivered a three-state
classification of GenBank deflines. Both are described by their own design docs,
listed below.

## The design documents, and what each is for

Everything lives in `dev/design/`, and implementation plans sit under
`dev/design/plans/`. Read them in numeric order if you need the whole story;
the later ones assume the earlier ones.

`08_MULTIUNIT_UX_AND_CONTAMINATION.md` is the original contamination plan and
remains the parent document for all of this work. It has been amended twice
during implementation, and both amendments matter. Section 3.2.1 records what
was actually measured about the identity envelope against real data, and section
3.5 channel (b) records why the tool-kingdom channel was retired rather than
built.

`09_TRIAGE_RULES_AND_DEPTH.md` covers the replacement of an anonymous if/else
chain with six named, individually controllable triage rules, and the removal of
sequencing depth from the keep/drop decision. Merged.

`10_HIT_CLASSIFICATION.md` covers the three-state classification of a GenBank
defline into rdna, other, or uninformative, plus a bitscore-weighted consensus
that replaced reading a verdict off the single top hit. Merged.

`11_SCREEN_VERDICT_LAYER.md` covers the screen verdict layer: a recomputable,
versioned `screen` table holding one contamination verdict per assembly unit,
scored from whether the unit's whole BLAST hit distribution agrees on a taxon and
whether that taxon sits on the sample's own lineage. This is what is on the
unmerged branch.

`12_INTRAGENOMIC_VARIANTS.md` is the next sub-project and has a spec but no
implementation plan yet. It separates two things that currently look identical to
a curator: assembly-graph isomers that differ only in how many intergenic-spacer
subrepeats were traversed, and genuine intragenomic paralogous rDNA copies. The
first is a decision that does not matter and the second determines what GenBank
receives.

## What is next, and what is blocked

The immediate next step is to write an implementation plan for
`12_INTRAGENOMIC_VARIANTS.md` and execute it. One caveat: that sub-project builds
directly on the `screen` table and on `SCREEN_VERSION` from the unmerged branch,
so if hand-testing of that branch forces changes, the plan will need reworking.
It may be worth waiting for the maintainer's test result before writing it.

After that, Phase 3 has nothing unblocked remaining, and this is worth
understanding rather than rediscovering.

The positional and chimeric channel is blocked on the identity envelope. That
envelope was deferred after measurement, not out of laziness. The reference
project's 1833 real BLAST hits contain 1172 at class distance, 2 at phylum
distance, and none beyond, so there is no data in the band where an envelope
would act. Worse, genuine host 18S reaches 99.6 percent identity at class
distance, so no ceiling can sit below that without flagging real host data.
Section 3.2.1 of doc 08 has the full measurement. Unblocking this needs identity
data at phylum distance and beyond, which could be generated locally by fetching
a handful of distant rDNA records by accession and running local blastn against
the stored scaffolds; that was proposed and not chosen.

Cross-sample screening and index hopping are blocked on `samples.batch`, which
does not exist in the schema. Without run, lane or plate metadata the comparison
is meaningless.

## Conventions this codebase enforces, which are easy to violate

`NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd
files with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only, and add any
new export line to `NAMESPACE` by hand.

All R code, tests, fixtures, Nextflow files and documentation are ASCII only. No
em dashes, no smart quotes, no non-ASCII bytes anywhere. The maintainer is
explicit about this.

Comments are minimal and explain why, not what changed. Do not narrate a bugfix
in a comment.

Database migrations are additive and idempotent. Nothing is ever dropped or
retyped. `RIBOPILOT_SCHEMA_VERSION` in `R/migrate_db.R` is currently 14 and moves
in exactly one place per schema change. Shared DDL constants live in
`R/migrate_db.R` and are referenced by both `R/init_db.R` and the migration so a
fresh database and a migrated one cannot diverge; the test fixtures use the same
constants, so adding a column reaches the tests automatically.

The infix operator `%|NA|%` in `R/utils.R` is scalar only and errors on a vector.
This has caused defects in three separate tasks. `%||%` tests length rather than
NA, so it does not coalesce an NA.

Auto-drop behaviour must not ship enabled before a labelled fixture set exists.
Doc 08 section 5 states this as a gate, and the fixture set does not exist. Every
rule and verdict added so far ships as a flag, with drop available for a curator
to arm deliberately.

## How the work has been run, and why it matters

Each sub-project went spec, then implementation plan, then task-by-task execution
with a fresh subagent per task, an independent review after each task, and a
whole-branch review at the end. That process caught defects that reading alone
did not, and several of the worst were defects in the spec rather than in the
implementations.

The single most valuable habit has been verifying that a test actually fails when
its target is broken. Mutation testing, in other words: revert the fix, confirm
the test fails, restore, confirm the tree is clean. Three tests on this branch
looked correct, passed, and pinned nothing at all. One compared two values
computed by the same function, so breaking that function left it green. Another
had symmetric fixtures, so removing a `WHERE` clause still produced correct rows.
If you take one practice forward, take this one.

The second most valuable habit has been treating a subagent's claim as a
hypothesis. Several were wrong: a "pre-existing test failure" that had in fact
been introduced by the branch, a claim that a DBI helper nests transactions when
it does not, and an assertion that a comment described what the code guaranteed
when it did not. Check before acting.

## Known open items, recorded rather than forgotten

Doc 10 section 10 and doc 11 section 10 each carry a short list of deliberate
deferrals with the reasoning. The notable ones are that the screen's confidence
weights are priors rather than measurements, that `samples.reference_support` is
computed and written but has no reader yet, and that the intergenic spacer is
never annotated, so the isomer test in doc 12 infers spacer length from the
unannotated remainder rather than measuring it.
