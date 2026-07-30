# Handoff: RiboPilot contamination work, after Phase 3 and the validation gate

Written 2026-07-30. Repo `/home/dmacguig/Documents/GitHub/RiboPilot`, GitHub
remote `dmacguigan/RiboPilot`. Supersedes `HANDOFF_phase3.md`, which is now
historical: everything it listed as "next" is merged.

## Where the code stands

`main` is at `fa6a4c7` and everything is merged. Pull requests 1 through 8 are
all MERGED; there is no open work and no unpushed branch. The suite on `main` is
**979 passing, 0 failing, 0 warnings**. `devtools::check(document = FALSE,
vignettes = FALSE)` gives 0 errors, 0 warnings and 2 NOTEs, both long-standing
(crayon/dbplyr imports, `new_db.Rd` line width).

`RIBOPILOT_SCHEMA_VERSION` is `15L`. `SCREEN_VERSION` is `4L`. Each appears in
exactly one place, `R/migrate_db.R:8` and `R/screen.R:16`.

The Docker image `ribopilot:test` was last built from the 3c branch and carries
schema 15 and screen version 4. The validation-gate work changed no package code,
so the image is still current, but rebuild with `docker/deploy-local.sh test`
before any hand test if you have touched `R/`.

What merged, in order: PR 4 named triage rules and removed depth from the
keep/drop decision. PR 5 (3a) three-state hit classification. PR 6 (3b) the
`screen` verdict table. PR 7 (3c) intragenomic isomers and paralogs. PR 8 the
validation fixture set.

## The design documents, and what each is for

All in `dev/design/`, implementation plans under `dev/design/plans/`. Read in
numeric order if you need the whole story; later ones assume earlier ones.

`08_MULTIUNIT_UX_AND_CONTAMINATION.md` is the parent document for all of this and
has been amended four times during implementation. The amendments matter more
than the original text in places. Section 3.2.1 records what was measured about
the identity envelope. Section 3.5 channel (b) records why the tool-kingdom
channel was retired. Section 5 records that the validation gate has SPLIT into a
delivered half and a data-blocked half.

`09_TRIAGE_RULES_AND_DEPTH.md`, `10_HIT_CLASSIFICATION.md`,
`11_SCREEN_VERDICT_LAYER.md`, `12_INTRAGENOMIC_VARIANTS.md` and
`13_VALIDATION_FIXTURES.md` cover the merged sub-projects. Each carries in-place
amendments recording where its own spec turned out to be wrong. **Read the
amendments, not just the body.** In three of the five sub-projects the worst
defects were defects in the spec rather than in the implementation.

## What is next

**Nothing is unblocked in the contamination plan.** That is the single most
important fact in this document, and it is worth confirming rather than
rediscovering. Every remaining item is blocked on data or on a schema field that
does not exist.

The concrete next piece of work, and the one that unblocks the most, is
**getting the three invertebrate clade projects to exist**. It is the same
blocker for three separate things:

1. The clade half of doc 08 section 5's validation gate: the invertebrate
   false-positive rate. The rate currently measured is actinopterygian-only,
   which is narrower than the gate asks for and narrower than where a false
   positive actually hurts.
2. `min_isomer_anchor_pid` (doc 12), still an unvalidated prior. It needs real
   multi-unit host samples carrying genuine paralogs. Spike-ins can never reach
   it: variant comparison runs only on `on_target` units and a spike-in is
   `off_target` by construction.
3. The `human` subclass of `screen_class_for()`, which is only reachable on a
   non-chordate host.

What that work requires, verified 2026-07-30 rather than assumed:

- The three read sets are NOT downloaded. `find ~/Documents -name "DRR251133*"`
  and the same for `ERR3597993` and `SRR25308289` return nothing.
- Their clade reference DBs are NOT built. `ref_dbs/getOrganelle/seeds/` holds
  only `actinopterygii_rrna_seeds.fasta`.
- No invertebrate project has ever been run. `mapping_asteroidea.csv`,
  `mapping_scyphozoa.csv` and `mapping_octocorallia.csv` in
  `~/Documents/DJMCompBio/RiboPilot_test/` each define exactly ONE sample, and
  every existing test project is actinopterygian.

So it is read downloads, three reference DBs, three assemblies and three remote
BLAST runs before any labelling can begin. That is a data-gathering task with a
long wall clock, not a coding task. Discuss the cost with the maintainer before
starting it.

## What is blocked, and on precisely what

**Sub-project 3d, the positional and chimeric channel**, is blocked on the
identity envelope, which was deferred after measurement rather than out of
laziness. Doc 08 section 3.2.1 has the numbers: the reference corpus holds 1172
hits at class distance, 2 at phylum and none beyond, so there is no data in the
band where an envelope would act; and genuine host 18S reaches 99.6 percent
identity at class distance, so no ceiling can sit below that without flagging
real host data. Unblocking it needs identity data at phylum distance and beyond,
which could be generated locally by fetching distant rDNA records by accession and
running local blastn against the stored scaffolds. That was proposed and not
chosen. **If you take this on, note that the local-BLAST machinery now exists** in
`dev/make_screen_fixture.R`, which is most of what such a measurement needs.

**Sub-project 3e, cross-sample screening and index hopping**, is blocked on
`samples.batch`, which does not exist. Verified: zero references in both
`R/init_db.R` and `R/migrate_db.R`. Adding the column is trivial; the detection
is meaningless until users actually supply run, lane or plate metadata, which doc
08 section 3.6 states plainly.

**Auto-drop and the export gate** remain disabled by choice, not by blocker. The
measurement they were waiting on now exists (doc 13): zero false positives on 19
real clean-host units and 1833 hit rows, recall 4 of 4 on catchable spike-ins.
Arming them is a maintainer decision, and the honest caveat is that the
false-positive rate is measured on fish only.

## Unblocked work that is NOT in the contamination plan

If the clade data is too expensive right now, these are real and available. All
were verified 2026-07-30.

**Design debt.** `07_SCAFFOLD_SELECTION.md` still describes stitching as the
design and says the Level-1 stitch needs integrating, which is the opposite of
what commit `b0f1b12` did. The decision to drop stitching is recorded in NO
document. `compute_feature_alignment()` has zero callers and is absent from
`NAMESPACE`, so it can be deleted. Careful: doc 08 also lists `stitch_scaffolds`
as deprecable and **that is wrong** - `R/pipeline_steps.R` and
`R/annotate_project.R` still call it.

**Per-position minor-allele fraction** (doc 08 section 5). A submitted rDNA
scaffold is a consensus over many non-identical array copies, and that is the
first question a systematist asks about a genome-skim submission.
`coverage.depths` already stores the per-position track. Check whether the mapped
reads are still available before promising this; the depth track alone may not be
enough.

**The corpus gaps in doc 13 section 5.1.** Three named cases the fixture cannot
catch. Closing the congener one needs a host unit whose hits are dominated by
same-genus relatives, which none of the ten existing samples provides.

## Conventions this codebase enforces, which are easy to violate

`NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd
with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only, and add any new
export line by hand.

All R code, tests, fixtures, Nextflow files and documentation are ASCII only. No
em dashes, no smart quotes. The maintainer is explicit about this.

Comments are minimal and explain why, not what changed. Do not narrate a bugfix
in a comment.

Database migrations are additive and idempotent. Nothing is dropped or retyped.
Shared DDL constants live in `R/migrate_db.R` and are referenced by both
`R/init_db.R` and the migration, so a fresh database and a migrated one cannot
diverge.

The infix `%|NA|%` in `R/utils.R` is **scalar only and errors on a vector**. This
has now caused defects in four separate tasks, including 3c. `%||%` tests length
rather than NA, so it does not coalesce an NA.

`dev/` is in `.Rbuildignore` and contains 0 entries in the built tarball. Dev-only
scripts belong there and never become package code.

Never push without being asked. Branch from `main` and ask first for anything
sizeable.

## How the work has been run, and why it matters

Each sub-project went spec, then implementation plan, then task-by-task
execution, with an independent review after each task and a whole-branch review
at the end. Design docs live at `dev/design/NN_TOPIC.md` and plans at
`dev/design/plans/YYYY-MM-DD-topic.md`. Follow that, not any tool's default path.

**Mutation testing is the single practice to carry forward.** Revert the fix,
confirm the test fails, restore, confirm the tree is clean. It has repeatedly
found tests that passed while pinning nothing, and in 3c it found a LIVE defect
that a passing test was hiding: the variant-clause separator contains `|`, which
as a regex is alternation, so `sub()` matched a bare space and cut a word out of
the verdict's own reason. The test meant to pin it used `force = TRUE`, which made
`screen_one()` write a fresh reason first, so the strip never ran.

Two lessons from that worth internalising:

- **A mutation that does not bite is information, not a nuisance.** In the
  validation work, two mutations that failed to bite revealed that no fixture
  unit exercised the 0.6 confidence gate at all, and that the phylum-or-deeper
  rule is enforced redundantly at two points. Both became recorded findings.
- **Check the mutation actually applied.** A `sed` that silently matches nothing
  looks exactly like a test that correctly survived.

**Treat every claim as a hypothesis, including the design docs' own.** Confirmed
wrong so far: doc 11 named the wrong blind spot (a copepod is an arthropod and a
fish a chordate, so that pair differs at phylum and IS caught); doc 08 asked for a
human spike-in as though it were cross-kingdom (it is not, on a chordate host);
doc 08 lists `stitch_scaffolds` as dead when it has live callers; doc 12 sited a
params sublist that had never been built.

**Verify with a real project, not only the suite.** The 3c hand test found a
five-unit false `isomer` group that 883 passing tests had not: two units scored
"100.00 percent identical" over an alignment covering 0.2 percent of one of them.
The test project is `~/Documents/DJMCompBio/RiboPilot_test/full_project_docker`.

## Known open items, recorded rather than forgotten

- `samples.reference_support` is computed and written on every `screen_sample()`
  run and has no reader. Doc 11 section 3 says the UI must say so when it is 4 or
  worse; that UI is deferred.
- The screen's confidence weights are priors, not measurements. Doc 13 measures
  the verdicts they produce but does not fit them, and refitting against six
  categories would overfit. Any change must bump `SCREEN_VERSION`.
- `min_isomer_anchor_pid` at 99.5 has never been exercised as a discriminator.
  In the 3c hand test every comparable pair scored exactly 100.00 and every
  rejection was decided by coverage or by disjoint element sets.
- The intergenic spacer is never annotated, so doc 12's isomer message infers
  spacer copy number from the unannotated remainder rather than measuring it.
- Doc 12 section 8's two salvage items are deferred with measurements: extending
  `parse_barrnap()`'s name map retains nothing without a store, and the ITSx
  taxonomic profile cannot reach R because both Nextflow processes declare only
  `*.itsx.positions.txt` as output.
- The `human` subclass is unreachable end to end from the curated fixture
  databases. Covered as a unit test only.
- All strands in the test data are `+`, so the all-minus and mixed-strand anchor
  paths in `R/screen_variants.R` have unit-test coverage only.

## Practical notes that will save you time

**Remote BLAST is the bottleneck and it fails silently.** NCBI runs BLAST compute
as a shared queue that stops answering under load: a capture once sat 85 minutes
at 00:00:00 CPU on a query that normally takes minutes. Their file serving keeps
working throughout, which is the distinction to exploit. `blastn -remote` has no
timeout, so always wrap it in `timeout`. Zero CPU time with long elapsed time is
the signature.

**Prefer local BLAST against NCBI's prebuilt rRNA databases** where fidelity
allows. `dev/make_screen_fixture.R` shows the pattern: curl the tarball, untar,
set `BLASTDB`, search locally. Each tarball bundles `taxdb.bt*`, which is what
supplies `staxids`. Do NOT use `update_blastdb.pl` from the conda package: it
cannot find curl, leaves `$curl` uninitialised at its line 599, tries to exec
`--user` as a program, and reports success anyway.

**The image has RiboPilot installed but not the dev dependencies.** No `pkgload`,
no `withr`. Scripts meant to run inside it must fall back to
`library(RiboPilot)`.

**Watch for missing ID predicates in SQL.** This codebase has now been bitten
twice, once in `compute_feature_alignment()` (doc 08 line 232) and once in the
fixture generator, where a readback lacking `ID = ...` pulled a copied project's
own hits and produced a contaminant fixture that was actually host data.
