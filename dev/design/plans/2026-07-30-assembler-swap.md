# Implementation plan: one workflow, two assemblers

Spec: `dev/design/14_ASSEMBLER_SWAP.md`. Branch `assembler-swap`, from `main`.
Schema v16.

Tasks are ordered so the tree is loadable and the suite green after each one.
Each task ends with the suite run and, where a task fixes something, a mutation
check: revert the fix, confirm a test fails, restore, confirm the tree is clean.

## Task 1 - schema v16: the `assembler` column

- `R/migrate_db.R`: bump `RIBOPILOT_SCHEMA_VERSION` to `16L`. Add the v16 block:
  guarded `ALTER TABLE assemble_opts ADD COLUMN assembler TEXT DEFAULT
  'getorganelle'`, then `UPDATE assemble SET assemble_switch = 1 WHERE
  assemble_switch = 5`.
- `R/init_db.R:211-223`: add `assembler TEXT DEFAULT 'getorganelle'` to the
  `CREATE TABLE assemble_opts` DDL, and to the seed row at `:225-239` if it
  enumerates columns.
- Tests (`tests/testthat/test-migrate.R` or its equivalent): a v15 database
  gains the column defaulting to `getorganelle`; a row at `assemble_switch = 5`
  moves to 1; running the migration twice changes nothing the second time; a
  fresh `new_db()` and a migrated database agree on the `assemble_opts` column
  set.

The `INSERT OR REPLACE` column-wipe trap applies to the seed row: list every
column or unlisted ones go NULL.

## Task 2 - `query_samples` returns the assembler

- `inst/nextflow/modules/pipeline_r.nf:7-32`: add `COALESCE(ao.assembler,
  'getorganelle') AS assembler` to the select list and `assembler` to the
  explicit header `printf`. The header is written by hand because `sqlite3
  -header` emits nothing for zero rows, so both must change together or
  `splitCsv(header: true)` misreads every column.
- Test: textual assertion that the header line and the select list both name
  `assembler`, and that they agree on column count.

## Task 3 - branch and merge in `main.nf`

- `inst/nextflow/main.nf`: include `fallback_assemble` from
  `./modules/fallback.nf`. Extend `asm_opts` (`:37`) with `r.assembler`. Branch
  `asm_in` into `go` / `fb`, map each to its process's input arity, and feed
  `score()` from `assemble.out.paths.mix(fallback_assemble.out.paths)`.
- The `fb` branch condition is `true`, not an equality test, so an unrecognised
  value routes to the fallback rather than being dropped. Nextflow's `branch`
  silently discards items matching no condition, which is the failure mode this
  work exists to stop repeating.
- Update the header comment (`main.nf:1-12`) to say the assembler is per-sample.
- Tests, in the style of `tests/testthat/test-nextflow.R:34-60`: `main.nf`
  branches on `assembler`; both branches reach `score` through `mix`; the
  GetOrganelle branch still receives seeds and labels; `main.nf` is the only
  workflow file naming `fallback_assemble`.

## Task 4 - delete `main_fallback.nf` and its query

- Delete `inst/nextflow/main_fallback.nf`.
- Delete `query_fallback_samples` from `inst/nextflow/modules/pipeline_r.nf:33-53`
  (its only caller goes with the file).
- `R/nf_command.R`: drop the `"fallback"` arm from `match.arg` (`:16`), the
  `main_fallback.nf` switch branch (`:20`), the `method` parameter (`:17`) and
  the `--method` line (`:38`).
- Test: no workflow file references `stitch` or `query_fallback_samples`; the
  three shipped workflow files are `main.nf`, `main_annotate.nf` and nothing
  else.

`modules/fallback.nf`, `bin/rp_fallback.R` and `R/assemble_fallback.R` all
survive untouched. The process is now called from `main.nf`.

## Task 5 - rework `run_fallback()`

- `R/run_fallback.R`: keep the exported signature. Replace the body's queueing
  half so it (a) resolves candidates as now, (b) derives a per-(source preset,
  method) `assemble_opts` preset carrying `assembler = method`, copying the
  source preset's other columns, (c) sets `assemble_switch = 1`, (d) runs the
  ordinary assembly workflow.
- Deriving rather than editing in place is load-bearing: editing the shared
  preset would change the assembler for every other sample using it, the same
  defect doc 08 records at `R/app_opts.R:347-351`.
- Fix the completion message (`:70`) to name the real WF1 tail.
- Fix the stale contract comment at `R/assemble_fallback.R:8-10`, which still
  says "load -> score -> stitch -> annotate".
- Tests against a temporary project: candidates are queued at 1 not 5; a derived
  preset is created with the right `assembler`; a sample sharing the source
  preset is not repointed; calling twice does not create a second duplicate
  preset; `ids` overrides candidate resolution.

## Task 6 - documentation reconciliation

- `dev/design/08_...md`: record that open question 1 is answered and question 2
  narrowed; correct the line 1304 claim about `assemble_fallback.R` writing
  `assemblies`.
- `dev/design/02_ASSEMBLY.md:155-159`: record that B2 is not automatic, it needs
  a `run_fallback()` call.
- `vignettes/RiboPilot.Rmd:206-216`: reword Step 5 so it no longer implies a
  separate workflow. Rebuild the vignette only if the repo tracks the built
  output.
- `NAMESPACE` is hand-maintained; no export changes in this work, so it should
  end unchanged. Regenerate Rd with `roxygen2::roxygenise(roclets = "rd")` only
  if a documented signature changed. Never `devtools::document()`.

## Task 7 - whole-branch review and verification

- Full suite, expecting 979 plus the new tests, 0 failures, 0 warnings.
- `devtools::check(document = FALSE, vignettes = FALSE)`, expecting 0 errors,
  0 warnings and the two long-standing NOTEs.
- ASCII-only check across every changed file.
- Confirm a fresh `new_db()` and a v15 database migrated to v16 produce the same
  `assemble_opts` column set.
- Independent review of the whole branch before the PR.

## Deferred to a follow-up, recorded so it is not lost

- Doc 07's stitch-caller sentence needs amending once PR 9 merges (spec section
  10). It is on the other branch, so it cannot be done here.
- The assemble-options modal (`R/app_opts.R:52`) gains no assembler control in
  this work (spec section 9).
- The fallback assembler is still unvalidated end to end; the suite cannot reach
  `bowtie2` or `spades.py` (spec section 8).
