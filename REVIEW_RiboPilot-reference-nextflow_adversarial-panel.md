# RiboPilot adversarial review, branch `reference-nextflow`

Twelve reviewers, one lens each: bioinformatician, Nextflow specialist, HPC admin,
geneticist, systematist, project manager, R specialist, database/SQL, GenBank
expert, and three end users (beginner, intermediate, advanced).

Findings marked **[V]** I re-checked against the code myself. Unmarked findings
come from a single reviewer and are worth confirming before acting.

Test suite passes clean, before and after the edits listed at the bottom.

---

## The regex sweep you asked for

**Clean.** The `[^\s:]` character-class bug exists in exactly two places in the
whole repo, `R/reference_nf.R:91-92`, and both already pass `perl = TRUE`. The
test at `tests/testthat/test-nextflow.R:305` pins it. **[V]**

I swept R, Nextflow, shell, and Python for three related shapes: any backslash
escape inside a bracket expression, any PCRE-only syntax (lookaround, lazy
quantifiers) without `perl = TRUE`, and any regex built by pasting a variable in.
Nothing else turned up. Note that `\s` *outside* a character class is fine in R's
default engine, so the ~127 other regex calls are not at risk.

One near-miss worth knowing about: `R/run_assembly.R:31,35` builds a pattern by
pasting `key` in unescaped. Every current caller passes a literal like `rawDir`,
so it is safe today, but a key containing `.` or `+` would match wrongly.

---

## Critical

### 1. The rank-misalignment bug is still live in the other lineage producer **[V]**

`R/taxonomy.R:163` calls `taxonkit reformat -f` with **no `-F`**, then line 168
deletes the empty rank slots with `gsub("; *;", ";", ...)`. The corrected path at
`R/taxonomy_ranks.R:340` uses `reformat -F -t`.

So for any taxon missing a rank, the family name ends up sitting in the order
slot. That string goes into `samples.lineage`, which `R/reference_groups.R:29-36`
parses **positionally** and zips against the correctly-aligned taxid vector. The
mislabelled rung becomes `clade_name` and is sent to GenBank as
`"<wrong-name>"[Organism]` when building the reference database.

This is the exact bug the header of `taxonomy_ranks.R` says was fixed. It was
fixed in one of the two producers. Two reviewers found this independently.

Downstream: `R/scaffold_scoring.R:88-100` (`taxonomic_match`) walks an unbroken
name prefix over these same strings and carries the largest single scoring weight
(0.30). A sample whose lineage lacks a class or order scores near zero against a
*congeneric* hit, so `auto_select` shows the curator the wrong assembly path.
Fails hardest on exactly the under-described invertebrates the tool targets.

### 2. `new_project(executor = "slurm")` writes an unparseable config **[V]**

`R/init_project.R:146-154` substitutes 7 tokens. `CONTAINER_ENGINE`, `QUEUE`, and
`CLUSTER_OPTIONS` are not among them; only `generate_config()` fills those
(`R/generate_config.R:207-213`). The resulting `.config` contains literal
`<<QUEUE>>`, which is not valid Groovy.

`man/new_project.Rd` lists `"slurm"` as a usable `executor` value, so the docs
lead you straight into it. `fill_config()` (`R/generate_config.R:52-60`) silently
ignores unfilled tokens rather than erroring. Same problem for `awsbatch`, which
has three tokens neither function fills.

### 3. Concurrent DB writers with no busy timeout **[V]**

`open_project_db()` (`R/pipeline_steps.R:154-159`) does the right thing: sets
`busy_timeout=10000` and picks a network-safe journal mode. Two whole categories
of writer bypass it.

- The Shiny app: `app_con()` (`R/app_utils.R:71-73`) is a bare `dbConnect` with
  no pragmas. Default busy timeout is **0**, so the first contended statement
  errors immediately.
- The Nextflow modules: 17 raw `sqlite3` CLI calls across `pipeline_r.nf`,
  `blast.nf`, `blast_ref.nf`, `coverage.nf`. The CLI default timeout is also 0.

On a network filesystem the journal mode is `delete` by design
(`R/pipeline_steps.R:117,141-145`), which means whole-file exclusion. At
`queueSize = 50` this is a global mutex.

Worse, `query_samples` (`pipeline_r.nf:7-32`) has **no `errorStrategy`**, so
Nextflow's default `terminate` applies: one locked read kills the entire run.
Same for `score` (`:94`), `reconcile` (`:150`), and `curate` (`:198`).

The one-line fix for the app half: make `app_con()` call `open_project_db()`.

### 4. Migration serialization is not actually serialized **[V]**

`R/migrate_db.R:179-180`:

```r
owns_txn <- tryCatch({ DBI::dbExecute(con, "BEGIN IMMEDIATE"); TRUE },
  error = function(e) FALSE)
```

`SQLITE_BUSY` and "already inside a transaction" are indistinguishable here. On a
busy timeout, `owns_txn` is `FALSE` and the function runs every `ALTER TABLE` in
autocommit, unserialized, which is precisely the race the block exists to
prevent. Two processes can both see a column missing and both try to add it; the
loser dies on `duplicate column name`.

This matters more than it looks because `migrate_db()` is called **per unit** by
`R/blast_genbank.R:47` and `R/curate_rrna.R:206`, from inside parallel Nextflow
tasks.

### 5. `PRAGMA user_version` is stamped after `COMMIT`, on a false premise **[V]**

`R/migrate_db.R:536-544`. The comment says "user_version is a header field, not
table data, so it is set after the transaction commits." That is wrong:
`PRAGMA user_version` *is* transactional in SQLite and can be set inside.

A crash between the COMMIT and the PRAGMA leaves a fully-migrated schema stamped
at the old version. The next run re-executes every `ver <`-gated step, including
the v9 `DELETE FROM blast_hits` at line 405, destroying real data on a project
that has since re-BLASTed.

### 6. GetOrganelle runs `-F anonym` with no size bounds

`R/init_db.R:37-41`, used at `assemble.nf:20-25`. The reviewer traced this into
the installed GetOrganelle source: with no `--target-genome-size`,
`--expected-max-size`, or `--max-extending-len`, all three are derived from the
label database's total length. The shipped labels DB is ~3.08 Mb, giving a
6.16 Mb "target genome size" and a 30.8 Mb ceiling for what should be a ~10 kb
repeat, with extension effectively unbounded.

Consequence: on a sample with a long or subrepeat-rich IGS, extension walks out
of the rDNA array into flanking genomic sequence, and nothing downstream catches
it, because the length score decays only linearly and annotation still finds one
18S and one 28S somewhere inside.

### 7. A two-repeat assembly is silently annotated as one cassette

`R/annotate.R:98`: `pick()` takes `[1, ]`, the first hit per element. Nothing
anywhere counts hits per element. `export.R:192-195`'s misassembly check tests
for out-of-order annotations, which by construction cannot fire when there is
only ever one row per element.

So a scaffold that traversed the array twice gets annotated as a single clean
cassette, reports no risks, and is submitted to GenBank as a complete repeat unit
that is actually two concatenated units.

The same `[1, ]` also discards the second fragment of a naturally split 28S
(the "hidden break" in many arthropods and crustaceans), which then trips the
truncation flag as if the assembly were bad.

### 8. `linkage_evidence` used an invalid controlled-vocabulary value **[V]**

`R/export.R:143` emitted `align_genus`. INSDC's vocabulary is space-separated
(`align genus`), and the adjacent `gap_type within scaffold` on line 142 already
uses the correct form. Fires on every export with a gap. **Fixed below.**

---

## High

**Silent data loss and mis-scoring**

- `R/blast_ref_fetch.R:54` strips IUPAC ambiguity codes (`gsub("[^ACGTNacgtn]", "", ...)`)
  from fetched reference sequences, but the feature-table coordinates parsed at
  `:22-47` are 1-based on the *original* record. Every ambiguity code before a
  feature shifts that feature's coordinates. rDNA reference records routinely
  carry them at paralog-variable sites.
- The same function is why an NCBI error body is accepted as sequence: an XML
  `<ERROR>` response degenerates into a short garbage "sequence" rather than an
  empty one, so the `nzchar` guard at `:69` never fires. `blast_ref.nf:20-28`
  only checks the file is non-empty.
- `R/pipeline_steps.R:464`: `dep[from:min(length(dep), to)]` builds a *descending*
  sequence when the depth track is shorter than the trim start, writing all-NA
  depths. Also `trim_assembly()`/`restore_assembly()` update mean/min/max but not
  `median_depth`, which is the column the UI actually reads.
- `R/pipeline_steps.R:234` parses assembly FASTA headers as `ID.path.scaffold`
  positionally, with no length check, and sample IDs are never validated to be
  dot-free. `SRR123.1` collapses every scaffold onto one key.
- `inst/nextflow/main.nf:50`: `file(r.seeds_db)` throws on an empty string, which
  is what `sqlite3` renders SQL NULL as. `seeds_db` is settable to NULL
  (`R/reference_build.R:138`). This runs for *every* row during channel
  construction, including fallback samples that never use seeds, so one bad row
  aborts the whole workflow before a single task submits. **[V]**
- `inst/nextflow/modules/coverage.nf:29-30`: Nextflow's default shell is
  `bash -ue`, not `-o pipefail`. `bowtie2 | samtools sort` takes samtools' exit
  status, and samtools exits 0 on an empty stream. A dead bowtie2 produces a
  valid empty BAM and zero-depth rows, which triage then reads as a real signal.

**Contamination can't be caught**

- `R/screen.R:199`: `keep <- hits$hit_class %in% "rdna"` discards every non-rDNA
  hit before the LCA. A bacterial 16S or fungal WGS contaminant therefore yields
  `no_reference`, which `triage_rules.R:34-38` treats as "the screen has not run"
  and only flags. The contaminant is reported as *unjudged* rather than
  off-target, which is backwards.
- `R/triage_rules.R:39-51`: only `no_feat_non_rdna` defaults to `drop`, and it
  requires zero features. Real rDNA contaminants have features, so no shipped
  rule can drop one. `R/export.R` reads neither `screen` nor `triage` beyond
  `assemblies.ignore`, so a flagged off-target unit exports normally.
- `R/curate_rrna.R:182-197`: `lineages_conflict()` is a bag-of-tokens
  intersection. Both a coral and a fish lineage contain `Metazoa`, so no
  metazoan contaminant can ever raise the flag. The rank-aware `rank_distance()`
  already exists and should be used instead.

**Operational**

- `R/annotate.R:231` passes `kingdom = "euk"` into `run_barrnap()`, whose own
  default is `"fun"`. Perl barrnap ships only `arc/bac/fun`, and `stderr = FALSE`
  at line 190 swallows the error. `pipeline_r.nf:131,138` has the same defect.
  `R/init_db.R:533-535` documents the incompatibility in a comment and nothing
  acts on it. Where pybarrnap is absent, every scaffold silently gets zero
  anchors, orientation becomes a no-op, and triage starts dropping real
  scaffolds. **[V]**
- `R/app_run.R:43-44` and `:58` paste the Nextflow argv into a shell script with
  `collapse = " "` and no quoting, while `R/run_assembly.R:73` correctly uses
  `shQuote`. `nf_command()` can return a single argv element containing spaces
  (the container bind options), so the cluster-submit path breaks. **[V]**
- `R/nf_command.R:22-23`: `%||%` is a length-0 coalesce, so `NA_character_` from
  `config_param()` passes through. Line 24 special-cases this for `api_key` but
  not `rawDir` or `minDepth`, so a config missing a key launches with
  `--rawDir NA`. **[V]**
- No `maxForks` anywhere in the repo. `blast_genbank` and `blast_ref_fetch` both
  hit NCBI, and `config.local`/`docker`/`conda` set no `queueSize` either. On a
  32-core workstation that is 30+ simultaneous remote BLAST submissions from one
  IP. `errorStrategy 'ignore'` then converts the resulting rate-limit block into
  a "successful" run with zero screening. **[V, no maxForks confirmed]**
- Memory escalation is dead code in every cluster config
  (`config.slurm:26` and five siblings): `memory = { 16.GB * task.attempt }` with
  no `errorStrategy 'retry'` anywhere to make `task.attempt` exceed 1. And the
  four processes that actually need memory override it with a static DB value and
  carry `errorStrategy 'ignore'`, so an OOM is silently dropped. **[V]**
- The NCBI API key is written to `.config` with no permission tightening
  (`R/init_project.R:51,153-156`), passed as a bare CLI arg
  (`R/nf_command.R:37`), and baked into every `.command.sh` in the work dir. On
  shared scratch, readable by anyone in the group; visible in `ps` regardless.
- No config template sets a walltime. Sites that require an explicit `-t`/
  `-l walltime=` will reject every submission, and because the assembly processes
  carry `errorStrategy 'ignore'`, the run reports success with everything dropped.
- Docker configs set no `fixOwnership` and no `-u $(id -u)`, so SQLite's `-wal`
  and `-shm` files get created root-owned in the user's project directory. The
  app then cannot open its own database, and the user cannot delete the files
  without sudo.

**Scale**

- `inst/nextflow/main.nf:111`: `blast_ref_fetch.out.collect()` is a whole-run
  barrier. No unit is aligned, curated, or triaged until every sample's remote
  NCBI fetch finishes, on the flakiest stage in the pipeline.
- `inst/nextflow/bin/rp_env.R:29` runs `pkgload::load_all()` per task, re-parsing
  the entire package for every one of ~145k tasks at 5000-sample scale.
- `export_bundle()` computes `export_summary_df()` twice per bundle
  (`R/export.R:331` and `:476`), and `feature_map_svg()` (`:304-328`) builds a
  full ggplot per unit for a five-box diagram, all embedded in one HTML file.

---

## Medium (selected)

- `R/init_db.R:107-111` creates the `samples` table with **no column types at
  all**, and column names interpolated as string literals. It only works because
  SQLite tolerates quoted strings in a column-definition position. `taxid` and
  the switch/lock columns are numeric everywhere in the code but stored with
  whatever type R bound.
- No secondary indexes anywhere in the project DB. The only `CREATE INDEX` in the
  package is on the taxdump cache (`R/tax_search.R:54`).
- No foreign keys, `PRAGMA foreign_keys` never enabled; referential integrity is
  maintained by hand-written orphan sweeps.
- Delete-then-insert without a transaction in at least six places
  (`R/pipeline_steps.R:295,335`, `R/screen.R:315`, `R/blast_genbank.R:47`,
  `R/blast_ref_fetch.R:76`, `R/reference_build.R:36-63`). `R/triage_opts.R:54`
  gets it right and documents why; nothing else follows.
- `R/migrate_db.R:504-506` runs an ungated full-table `UPDATE assemble` on every
  call, and `migrate_db()` is called per unit.
- The migration backup (`R/migrate_db.R:644-653`) uses `file.copy()`, whose
  `FALSE` return is never checked, names the destination by epoch *seconds* with
  `overwrite = FALSE`, and does not copy the `-wal`/`-shm` files.
- Remote BLAST runs with no `-evalue` and no identity or coverage floor
  (`blast.nf:61-65`), and `blast_genbank.R:118-135` applies no floor before a
  subject becomes rank-1 reference candidate.
- fastp defaults do no quality trimming at all and leave `--length_required` at
  15 bp (`R/init_db.R:161`). At rDNA's depth, bad 3' tails are exactly what
  manufacture spurious paralog groups in `screen_variants.R`.
- bowtie2 runs with no MAPQ filter and no `-k`/`-a` (`coverage.nf:29`). Reads from
  the conserved 18S/28S core map equally well to every scaffold and ties break
  pseudo-randomly, but the app invites cross-scaffold depth comparison.
- `expected_length = c(4000, 12000)` (`R/scaffold_scoring.R:120`,
  `R/select_paths.R:95`) is clade-blind and never overridden by either production
  call site, despite the per-element lengths right next door being clade-tuned.
- `min_isomer_anchor_pid = 99.5` (`R/params_rrna.R:75`) is one global threshold
  for isomer-vs-paralog, never clade-tuned, applied to the three clades the
  vignette explicitly says carry genuinely divergent rDNA copies.
- `R/taxonomy.R:128-144`: exact name resolution takes the first taxid and marks it
  `match_type = "exact"` with no `alternatives`. Cross-kingdom homonyms (*Morus*
  the mulberry vs *Morus* the gannet; *Aotus*) resolve silently to whichever
  arrives first, presented to the curator as maximum confidence.
- The GenBank reference query is keyed on a taxon **name**
  (`download_rrna.sh:99`, `R/reference_counts.R:10`) against live NCBI, while
  every other decision uses the pinned local taxdump. `txid<ID>[Organism:exp]` is
  supported and the taxid is already in hand at `reference_groups.R:130-133`.
- `SCAFFOLD_COL_HELP` at `R/constants.R:154` describes `screen_pident` as "how
  much of the scaffold matched", which is `screen_qcovs`. Wrong statistic in the
  tooltip.
- `DESCRIPTION`: `stats`, `grDevices`, and `grid` are used via `::` but not
  declared; `dbplyr` is declared and never used. Both are R CMD check failures
  waiting to happen. `NAMESPACE` is hand-maintained while every file carries
  roxygen `@export` tags.
- `inst/nextflow/bin/*.R` (all of them): `on.exit()` at script top level never
  fires, so the advertised `dbDisconnect` never happens.
- `environment.yml` pins only Nextflow and BLAST. GetOrganelle, SPAdes, fastp,
  bowtie2, barrnap, and ITSx are all unpinned, so two users a month apart can get
  different assemblies. The Dockerfile *does* pin, so the container path is
  reproducible and the documented conda path is not.

---

## What the end users hit

**Beginner** stops at `README.md:95-120`. Steps 2 and 3 require a terminal, which
is never mentioned, defined, or linked. Sea stars are not fish, and nothing early
says an invertebrate project needs a whole separate reference-database workflow
plus NCBI Entrez search syntax.

The mapping-file errors are terse and don't name the offending row:
"Duplicate IDs found", "IDs must be no more than 18 characters"
(`R/init_db.R:52-72`). Museum voucher IDs like `USNM_Asteroidea_00123_2024` are 26
characters, and the user is left guessing which to truncate.

The dangerous click: **Skip** in the Reference tab (`R/app_reference.R:340-348`).
One click, no confirmation, and the sample keeps whatever database it had, which
for a fresh project is the hardcoded Actinopterygii default
(`R/init_db.R:34-35`). Sea-star rDNA assembled against a fish-tuned seed database
half-works often enough to look fine, show "Done", and get exported.

**Intermediate** hits the SLURM config trap (critical #2), then the container:
`new_project()` warns rather than stops that `macguigand/ribopilot:0.0.0.9000` is
not published on Docker Hub, and SLURM defaults to Singularity, so the documented
cluster path dead-ends at image pull.

Then the quiet one: samples below `min_depth` are dropped by a channel filter
(`main.nf:52`) with no row written anywhere. They sit at the pending icon,
indistinguishable from "still queued", forever.

Also: launching from the app kills the run when the browser session ends
(`R/app_server.R:770-773`), and no doc anywhere says to use tmux or a batch
wrapper. And `-resume` is on by default but never documented, so nobody knows
recovery is as simple as re-running.

Not validated at all: whether the R1/R2 files exist, and single-end input (not
supported and not documented as unsupported).

**Advanced** would fork rather than adopt, for three reasons: the SQLite
write-bus concurrency model at 5000 samples; the Shiny app being mandatory
(locks gate every stage transition and **there is no exported function to set a
lock**); and the fact that the decision layer they'd most want to extend, path
scoring and auto-select, is both hardcoded and has zero tests.

Their note on the API is fair: the 55 exports are "whatever `inst/nextflow/bin/`
needed to reach", all taking an open DBI connection, while lock, state, trim, and
opts assignment are `@noRd`. And the shipped pipeline itself uses `:::`
(`rp_triage.R:11-12`, `rp_ref_db.R:11,20`), so the export list is not the
contract.

---

## Test coverage gap

Files with zero test coverage, ranked by consequence:

| File | What's untested |
|---|---|
| `R/scaffold_scoring.R` | `score_candidates`, `taxonomic_match` — decides which scaffold the curator sees first |
| `R/select_paths.R` | `score_assembly`, `blast_best_hits`, `assembly_anchors` |
| `R/score_project.R` | `score_project`, `make_blast_db` |
| `R/taxonomy.R` | all 18 functions, including the `-F` bug above |
| `R/tax_search.R` | all 9 |
| `R/resolve_project_taxonomy.R` | all 3 |

Also untested and destructive: `trim_assembly()`, `restore_assembly()` (both
rewrite stored sequences and shift annotation coordinates), and
`orient_rrna_sense()` (reverse-complements stored sequences).

`R/app_server.R` reads 72 `input$` ids; the suite sets 9. Never exercised:
`export_go`, `trim_apply`, `feat_save`, `restore_btn`, `nf_start`, `submit_job`,
`mark_reviewed`. Every destructive and every export path.

Four tests assert on `deparse()` of a function body rather than running it
(`test-app-opts.R:110-144`, `test-nf-pending.R:56-70`). A refactor that preserves
the substrings and breaks the behavior passes. `test-app-scaffold-uid.R` already
proves `testServer()` works on this app, so these can be rewritten properly.

CI is green but runs one job: ubuntu, R release, no matrix, no coverage, no lint,
and the Nextflow pipeline is never actually executed.

---

## Repo health

Better than expected. `assembly_out/`, `out/`, `work/`, `test_run/`,
`testdata_full/` (2.6 GB), and `docs/` are all present on disk and correctly
gitignored and untracked. Working tree clean. All 55 exports have man pages.

One real issue: `ref_dbs/getOrganelle/{seeds,labels}/*.fasta` is 13 MB committed
directly to git, roughly half the `.git` directory, with no LFS and no ignore
rule, while every other data blob is fetched rather than stored.

Eleven stale local branches are fully merged; `refdb-picker` is 18 ahead / 4
behind with no open PR.

Stale docs: `vignettes/RiboPilot.Rmd:290` and `RiboPilot_test.R:21` both call
`RiboPilot::run_annotate()`, removed when WF2 was decommissioned.

---

## Changes applied

Simplifications and one-token fixes only. Everything else above is left for you
to triage.

1. `R/export.R:143` — `linkage_evidence align_genus` to `align genus`. Invalid
   INSDC vocabulary value, fires on every export with a gap.
2. `R/app_opts.R:20` — deleted the duplicate `%|NA|%`. There is no `Collate:` in
   DESCRIPTION, so files load alphabetically and `R/utils.R:9` overwrites it
   anyway. Dead code, and the codebase elsewhere already documents the strict
   `utils.R` version as the real one. Note that the call sites which *don't*
   chain `%||%` first (`R/app_opts.R:76,78,92,94`) were never getting the
   NULL guard they appear to have.
3. `inst/nextflow/modules/annotate.nf` — deleted. Included by nothing
   (`main.nf:26-32` imports `annotate_tools` from `pipeline_r.nf`), and stale:
   per-sample rather than per-unit, barrnap-only, `--kingdom fun`.
4. `inst/nextflow/main.nf:55` — removed an identity `.map` that reassembled an
   8-tuple into the same 8-tuple.
5. `R/reference_nf.R:87-88` — trimmed the three-line comment narrating the fixed
   regex bug down to one line.

Test suite passes after all five.

## Status

Branch `review-fixes-2026-08`. Suite: 2028 passing, 0 failures.

**Done:** all six attack-order items below, plus every HIGH finding except the
three listed as still open. Read-file validation at project init was added on
top (not a review finding; requested separately).

**Still open, deliberately:**

- Export does not refuse an off-target unit. It flags through `submission_risks()`
  instead. Making export block is a policy call, not a bug fix.
- No rule's default action changed. A bacterial contaminant now reaches
  `off_target` and feeds triage rule 6, which still ships as Flag. Auto-drop
  wants a bigger labelled corpus than the one bacterial unit in the fixtures.
- The scale findings (the `collect()` barrier, `pkgload::load_all()` per task,
  the single-page HTML report) are architectural and untouched.
- `config.NMNH_Hydra` and `config.NOAA_SEDNA` still set no walltime. They bypass
  `generate_config()`, so they need a hardcoded site value someone who knows
  those queues has to choose.

**Two follow-ups this work created:**

- Projects that already exist hold lineage strings in the old collapsed format,
  and nothing re-resolves them. Detectable (fewer than seven `;` fields). The
  user has said not to fix this.
- `SCREEN_VERSION` went 4 -> 5, so stored screen rows recompute on next run.

## Suggested order of attack

1. `app_con()` calls `open_project_db()`. One line, removes the app half of
   critical #3.
2. Add `-F` at `R/taxonomy.R:163` and drop the slot-collapsing `gsub`, then point
   `taxonomic_match()` at `rank_distance()`. Critical #1, and the reason
   auto-select misranks invertebrates.
3. Make `new_project()` refuse a raw scheduler name, or fill the tokens.
   Critical #2, and the first thing every cluster user hits.
4. Move the `user_version` PRAGMA inside the transaction. One line, critical #5.
5. `maxForks` on the two NCBI processes, `errorStrategy` on the four processes
   that have none, and `set -o pipefail`.
6. Backfill tests for `scaffold_scoring.R` and `select_paths.R` before changing
   any of the scoring behavior above.
