# BLAST Stage 4a: local search and honest reporting

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Search a local rDNA database composed with NCBI decoy sets instead of remote
`core_nt`, and when that panel cannot answer a unit, say so to the user rather than
falling back to a two-hour remote search.

**Architecture:** Remote BLAST is 94 percent of pipeline CPU, almost all of it NCBI queue
time. A local search of the shipped Metazoa panel takes 0.80 s at 8 threads against the
same query. There is deliberately no automatic escalation: a unit the panel cannot answer
is recorded and surfaced, and the user decides whether a remote search is worth their time.
Most of the reporting path already exists. A no-hit unit already reaches the screen as
`verdict = "no_reference"`, `reason_code = "no_hits"`, which reads as unchecked rather than
clean, and `screen.screen_reason` is already written. It is simply never rendered.

**Tech Stack:** R (DBI, RSQLite), BLAST+ 2.17, Nextflow DSL2, Shiny, testthat 3e.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md` (Stage 4, local path only; distribution
is Stage 4b)

## Global Constraints

- R >= 4.4.0 (`DESCRIPTION`); `docker/Dockerfile` pins `r-base=4.4.*`. No package outside
  DESCRIPTION Imports.
- ASCII only in R and Nextflow sources. No non-ASCII characters, no em dashes.
- Minimal comments.
- This project ships no schema migrations before release.
- Never `git push`. Never `git pull`, `git fetch`, `git merge`, or `git rebase`: the remote
  branch is deliberately stale and pulling corrupts the history.
- Run every command in the FOREGROUND. Do not background anything and do not set up
  monitors; the full suite takes a few minutes, just let it block.
- Run tests with:
  `Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/<file>", reporter="summary")'`
- Full-suite baseline entering this plan: 2535 passed, 0 failed, 0 errors, 6 skipped.
  Measure with exactly:
  `Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); r <- testthat::test_dir("tests/testthat", reporter="silent"); d <- as.data.frame(r); cat("failed:", sum(d$failed), "error:", sum(d$error), "skipped:", sum(d$skipped), "passed:", sum(d$passed), "\n")'`
- The built panel is at `/home/dmacguig/ribopilot_refdb_build/blastdb/2026.08.1/ribopilot_rrna`:
  597,988 sequences, 688,449,026 bases, BLAST v5.
- BLAST binaries are at `/home/dmacguig/miniconda3/envs/ribopilot/bin`.

## Facts this plan rests on

- **Multi-database search works and every database must share `-parse_seqids`.** Verified:
  a two-database search returns hits from both with `staxids` resolving on each side, and
  mixing a database that carries accession info with one that does not fails outright with
  `BLAST Database error: DB contains no accession info`.
- **`params.ref_dir` is the working analogue for `params.blastdb`.** `nf_ref_dir()`
  (`R/nf_command.R:67-75`) reads the project config column directly and returns `""` when
  unset; `nf_command()` appends `--ref_dir` unconditionally (`R/nf_command.R:48`);
  `main.nf:41` declares the default; `blast_ref.nf:22` guards with `[ -n '${params.x}' ]`;
  `reference_mount_dirs()` (`R/reference_nf.R:96-123`) makes the directory visible in a
  container.
- **A no-hit unit already reports as unchecked.** `screen_channel_taxonomy()`
  (`R/screen.R:195-198`) short-circuits on zero hits with `reason_code = "no_hits"` and
  `verdict = "no_reference"`, and `screen_one()` always writes one `screen` row per unit.
- **`screen.screen_reason` is written and never displayed.** Populated at `R/screen.R:351`;
  grep across `R/app_server.R`, `R/app_modals.R`, `R/app_utils.R` finds no reader.
- **`assemblies.blast_source` currently takes only `'core_nt'` or `'seeds'`.** A third
  value falls through `rt_blast_source_badge()` (`R/utils_reactable.R:208-224`) unstyled,
  and worse, `seed_hit_block()` (`R/app_modals.R:36`) tests
  `!identical(blast_source, "core_nt")` so it would tell the user a local-panel hit came
  "from the local seed file". That is a mislabelling bug waiting for this change.
- `blast_searched` is `(ID, path, scaffold, query_md5, n_hits, time_stamp)`
  (`R/db_schema.R:47-55`), written only by `mark_unit_searched()`
  (`R/blast_genbank.R:337-345`).

---

### Task 1: Compose the database list

**Files:**
- Create: `R/blast_local.R`
- Test: `tests/testthat/test-blast-local.R` (create)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `RP_DECOY_DBS`, character vector `c("16S_ribosomal_RNA", "SSU_eukaryote_rRNA", "LSU_eukaryote_rRNA")`.
  - `local_blast_db_arg(panel, decoy_dir = NULL)` returning a single space-separated string
    suitable for `blastn -db`, or `NULL` when `panel` is blank or its BLAST volume files are
    absent. Decoys that are absent are skipped with a warning rather than failing.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-blast-local.R`:

```r
# Composing the -db argument. blastn takes a space-separated list, and every
# database in it must carry accession info or the search errors outright.

mk_fake_db <- function(dir, prefix) {
  for (ext in c("nhr", "nin", "nsq")) {
    file.create(file.path(dir, paste0(prefix, ".", ext)))
  }
  file.path(dir, prefix)
}

test_that("a panel alone composes to itself", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  expect_equal(RiboPilot:::local_blast_db_arg(p), p)
})

test_that("decoys present are appended in order", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  dd <- withr::local_tempdir()
  for (n in RiboPilot:::RP_DECOY_DBS) mk_fake_db(dd, n)
  got <- RiboPilot:::local_blast_db_arg(p, dd)
  parts <- strsplit(got, " ", fixed = TRUE)[[1]]
  expect_equal(parts[1], p)
  expect_equal(basename(parts[-1]), RiboPilot:::RP_DECOY_DBS)
})

test_that("a missing decoy is skipped with a warning, not an error", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  dd <- withr::local_tempdir()
  mk_fake_db(dd, "16S_ribosomal_RNA")
  expect_warning(got <- RiboPilot:::local_blast_db_arg(p, dd), "SSU_eukaryote_rRNA")
  parts <- strsplit(got, " ", fixed = TRUE)[[1]]
  expect_length(parts, 2L)
})

test_that("a blank or absent panel yields NULL", {
  expect_null(RiboPilot:::local_blast_db_arg(""))
  expect_null(RiboPilot:::local_blast_db_arg(NULL))
  expect_null(RiboPilot:::local_blast_db_arg(NA_character_))
  d <- withr::local_tempdir()
  expect_null(RiboPilot:::local_blast_db_arg(file.path(d, "absent")))
})

test_that("no decoy directory means panel only, no warning", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  expect_silent(got <- RiboPilot:::local_blast_db_arg(p, NULL))
  expect_equal(got, p)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, `local_blast_db_arg` not found.

- [ ] **Step 3: Implement**

Create `R/blast_local.R`:

```r
# Local BLAST against the shipped rDNA panel, composed with NCBI decoy sets so a
# non-metazoan contaminant is identified rather than merely unmatched.

#' NCBI prebuilt databases searched alongside the panel.
#'
#' The panel is metazoan nuclear rDNA, so anything else returns nothing and the
#' screen loses the class it would have reported. These restore it. They are
#' useless as a reference panel (the eukaryote sets are years stale) but a decoy
#' identifier only has to resolve to a kingdom.
#' @noRd
RP_DECOY_DBS <- c("16S_ribosomal_RNA", "SSU_eukaryote_rRNA", "LSU_eukaryote_rRNA")

#' TRUE when a BLAST database prefix has its volume files on disk
#' @noRd
blast_db_exists <- function(prefix) {
  length(prefix) == 1 && !is.na(prefix) && nzchar(prefix) &&
    length(Sys.glob(paste0(prefix, ".n*"))) > 0
}

#' Compose the -db argument for a local search
#'
#' Every database in a multi-database list must be built with -parse_seqids;
#' mixing one that carries accession info with one that does not fails the
#' search outright.
#'
#' @param panel Path prefix of the RiboPilot panel.
#' @param decoy_dir Directory holding the NCBI decoy databases, or NULL.
#' @return A single space-separated string, or NULL when the panel is unusable.
#' @noRd
local_blast_db_arg <- function(panel, decoy_dir = NULL) {
  if (!blast_db_exists(panel)) return(NULL)
  out <- panel
  if (!is.null(decoy_dir) && length(decoy_dir) == 1 && !is.na(decoy_dir) &&
      nzchar(decoy_dir)) {
    for (n in RP_DECOY_DBS) {
      p <- file.path(decoy_dir, n)
      if (blast_db_exists(p)) {
        out <- c(out, p)
      } else {
        warning("decoy database not found, skipping: ", n, call. = FALSE)
      }
    }
  }
  paste(out, collapse = " ")
}
```

- [ ] **Step 4: Run tests, then the full suite, then commit**

```bash
git add R/blast_local.R tests/testthat/test-blast-local.R
git commit -m "feat: compose the local BLAST database list with decoys"
```

---

### Task 2: Record which search produced a unit's hits

Before the pipeline can tell a user "the local panel could not answer this", it has to know
which search ran. `blast_searched` already records that a unit was searched; extend it to
record how.

**Files:**
- Modify: `R/db_schema.R` (`BLAST_SEARCHED_DDL`), `R/blast_genbank.R` (`mark_unit_searched`)
- Modify: `inst/nextflow/bin/rp_searched.R`
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: nothing.
- Produces: `blast_searched.source` TEXT, one of `"local"` or `"remote"`;
  `mark_unit_searched(con, id, path, scaffold, query_md5, n_hits, source = "remote")`.
  The default keeps every existing caller correct.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("a searched unit records which search produced it", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 3L, source = "local")
  r <- DBI::dbGetQuery(con, "SELECT source, n_hits FROM blast_searched")
  expect_equal(r$source[1], "local")
  expect_equal(r$n_hits[1], 3L)
})

test_that("the search source defaults to remote", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT source FROM blast_searched")$source[1], "remote")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, no such column `source`.

- [ ] **Step 3: Add the column and the argument**

In `R/db_schema.R`, add to `BLAST_SEARCHED_DDL` after `n_hits INTEGER,`:

```
  source TEXT,
```

In `R/blast_genbank.R`, change `mark_unit_searched()` to take `source = "remote"` and write
it, keeping the `INSERT OR REPLACE` shape and the existing column order plus the new one.

In `inst/nextflow/bin/rp_searched.R`, accept an optional 7th argument for the source,
defaulting to `"remote"` when absent, and pass it through. Keep `rp_args`' mandatory count
at 6 so the existing call keeps working.

- [ ] **Step 4: Run tests, then the full suite, then commit**

```bash
git add R/db_schema.R R/blast_genbank.R inst/nextflow/bin/rp_searched.R tests/testthat/test-blast-parse.R
git commit -m "feat: record whether a unit was searched locally or remotely"
```

---

### Task 3: Plumb the database path into the pipeline

**Files:**
- Modify: `R/nf_command.R`, `R/reference_nf.R`, `inst/nextflow/main.nf`
- Test: `tests/testthat/test-nextflow.R` (extend)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `params.blastdb`, defaulting to `""`.
  `nf_blastdb(path)` in `R/nf_command.R` returning the configured panel prefix or `""`.

The panel prefix is stored per project. Add a `blastdb` column to `reference_config`
alongside `ref_dir`, since that table already holds the project's reference settings and is
already read directly by `nf_ref_dir()`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-nextflow.R`, reusing `mk_ref_config_project()` already
defined in that file around line 343:

```r
test_that("main.nf declares a blastdb param", {
  p <- system.file("nextflow", "main.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/main.nf"
  expect_true(any(grepl("params.blastdb", readLines(p, warn = FALSE), fixed = TRUE)))
})

test_that("nf_blastdb is blank when the project configures no panel", {
  d <- withr::local_tempdir()
  proj <- mk_ref_config_project(d, NA_character_)
  expect_equal(RiboPilot:::nf_blastdb(proj), "")
})

test_that("nf_blastdb returns the configured panel prefix", {
  d <- withr::local_tempdir()
  proj <- mk_ref_config_project(d, NA_character_, blastdb = "/opt/rp/ribopilot_rrna")
  expect_equal(RiboPilot:::nf_blastdb(proj), "/opt/rp/ribopilot_rrna")
})
```

`mk_ref_config_project(dir, ref_dir_value)` is defined at `tests/testthat/test-nextflow.R:343`
and takes those two POSITIONAL arguments. Add a third, `blastdb = NA_character_`, and widen
its INSERT to include the column. Its existing callers pass two positional arguments and
must keep working unchanged.

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL on all three.

- [ ] **Step 3: Implement**

Add `blastdb TEXT` to `REFERENCE_CONFIG_DDL` in `R/db_schema.R`.

Add `nf_blastdb(path)` to `R/nf_command.R`, modelled exactly on `nf_ref_dir()`: read
`reference_config.blastdb` directly from the project sqlite, return `""` on any error or
blank. Append `"--blastdb", blastdb` to the argument vector next to `--ref_dir`.

Declare `params.blastdb = ''` in `inst/nextflow/main.nf`.

In `R/reference_nf.R`, add the panel's containing directory to what
`reference_mount_dirs()` considers, so the database is visible inside a container. Use
`dirname()` of the configured prefix, and only when it is explicitly set.

- [ ] **Step 4: Run tests, then the full suite, then commit**

```bash
git add R/db_schema.R R/nf_command.R R/reference_nf.R inst/nextflow/main.nf tests/testthat/test-nextflow.R
git commit -m "feat: plumb the local BLAST panel path into the pipeline"
```

---

### Task 4: Search locally when a panel is configured

**Files:**
- Modify: `inst/nextflow/modules/blast.nf`
- Create: `inst/nextflow/bin/rp_dbarg.R`
- Test: `tests/testthat/test-nextflow.R` (extend)

**Interfaces:**
- Consumes: `local_blast_db_arg()` (Task 1), `params.blastdb` (Task 3),
  `mark_unit_searched(source=)` (Task 2).
- Produces: `rp_dbarg.R <panel> [decoy_dir]` printing the composed `-db` string to stdout,
  or printing nothing and exiting non-zero when the panel is unusable.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-nextflow.R`:

```r
test_that("blast.nf searches locally when a panel is configured", {
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("params.blastdb", txt, fixed = TRUE))
  expect_true(grepl("rp_dbarg.R", txt, fixed = TRUE))
  expect_true(grepl("-remote", txt, fixed = TRUE),
    info = "the remote path must remain for projects with no panel")
})

test_that("blast.nf records the search source", {
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("rp_searched.R", txt, fixed = TRUE))
  expect_true(grepl("SRC", txt, fixed = TRUE))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: the first two assertions of test one FAIL.

- [ ] **Step 3: Add the composer script**

Create `inst/nextflow/bin/rp_dbarg.R`, following `inst/nextflow/bin/rp_rederive.R`:

```r
#!/usr/bin/env Rscript
# Print the composed -db argument for a local search, or exit 1 if unusable.
# Args: <panel_prefix> [decoy_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 1, "rp_dbarg.R <panel_prefix> [decoy_dir]")
local({
  # decoys ship beside the panel, so default to its own directory; the optional
  # second argument exists only to let a site keep them somewhere else
  dd <- if (length(a) >= 2 && nzchar(a[2])) a[2] else dirname(a[1])
  arg <- suppressWarnings(RiboPilot:::local_blast_db_arg(a[1], dd))
  if (is.null(arg)) quit(status = 1L)
  cat(arg, "\n", sep = "")
})
```

`chmod +x inst/nextflow/bin/rp_dbarg.R`. Check `tests/testthat/test-nf-bin-contract.R` and
satisfy every rule it enforces. Note this script does not open the project database, so the
contract's connection rules may not apply; confirm which rules do.

- [ ] **Step 4: Add the local branch to the module**

In `inst/nextflow/modules/blast.nf`, before the `for U in \$UNITS` loop, resolve the
database argument once per task:

```bash
    DBARG=""
    SRC="remote"
    if [ -n '${params.blastdb}' ]; then
        if DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}'); then
            SRC="local"
        else
            echo "configured local panel is unusable, falling back to remote for ${id}" >&2
            DBARG=""
        fi
    fi
```

Then branch the two `blastn` invocations. When `SRC` is `local`, drop `-remote`, use
`-db "\$DBARG"`, and add `-num_threads ${task.cpus}`; otherwise keep the existing
`-remote -db core_nt` call unchanged. `-entrez_query` is remote-only
(`inst/tool_help/blastn.txt:289`), so omit it on the local path and log once that it was
ignored. Keep the dc-megablast second pass on both paths; locally it is nearly free.

Add a `cpus` directive to the process, following `inst/nextflow/modules/blast_ref.nf`.

Pass the source through to the marker: `rp_searched.R ... "\$QMD5" "\$NHITS" "\$SRC"`.

Note the fallback in that snippet is to the REMOTE search, not to nothing: a misconfigured
panel must not silently stop searching. That is different from the shelved automatic
escalation, which was about a correctly-configured panel that simply had no answer.

- [ ] **Step 5: Verify against the real database**

Set up a scratch project pointing at the built panel and run one unit through the local
path by hand, confirming hits are stored and `blast_searched.source` reads `local`. Report
the wall time.

- [ ] **Step 6: Run tests, run nextflow lint if available, commit**

```bash
git add inst/nextflow/modules/blast.nf inst/nextflow/bin/rp_dbarg.R tests/testthat/test-nextflow.R
git commit -m "feat: search the local panel when one is configured"
```

---

### Task 5: Label the source honestly and stop mislabelling it

`assemblies.blast_source` has only ever held `'core_nt'` or `'seeds'`. Adding a third value
without touching the UI would tell the user a local-panel hit came "from the local seed
file", because `seed_hit_block()` tests `!= "core_nt"`.

**Files:**
- Modify: `R/blast_genbank.R` (`derive_blast_products`), `R/app_modals.R`,
  `R/utils_reactable.R`
- Test: `tests/testthat/test-blast-parse.R` (derivation) and
  `tests/testthat/test-blast-provenance.R` (the UI sites). `seed_hit_block()` is already
  covered there. `rt_blast_source_badge()` has NO test anywhere, so add one.

**Interfaces:**
- Consumes: `blast_searched.source` (Task 2).
- Produces: `assemblies.blast_source` gains `'local_rrna'`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("a locally searched unit is labelled local_rrna, not core_nt", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 1L, source = "local")
  RiboPilot:::derive_blast_products(con, "S1", 1, 1)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT blast_source FROM assemblies WHERE ID='S1'")$blast_source[1], "local_rrna")
})

test_that("a remotely searched unit is still labelled core_nt", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT blast_source FROM assemblies WHERE ID='S1'")$blast_source[1], "core_nt")
})
```

- [ ] **Step 2: Run tests to verify the first fails**

Expected: the first FAILS with `core_nt`; the second PASSES.

- [ ] **Step 3: Set the label from the search record**

In `derive_blast_products()`, look up `blast_searched.source` for the unit and write
`'local_rrna'` when it is `"local"`, otherwise `'core_nt'`. A missing row means remote,
which preserves current behaviour for every existing project.

- [ ] **Step 4: Fix the two UI sites**

In `R/app_modals.R`, `seed_hit_block()` currently decides "from the seed file" with
`!identical(r$blast_source %|NA|% "", "core_nt")`. Change it to test explicitly for
`"seeds"` so a local-panel hit is not mislabelled. Read the surrounding block and keep its
wording correct for all three values.

In `R/utils_reactable.R`, extend `rt_blast_source_badge()` so `'local_rrna'` renders with
its own label rather than falling through unstyled. Match the existing badge style.

- [ ] **Step 5: Run tests, then the full suite, then commit**

```bash
git add R/blast_genbank.R R/app_modals.R R/utils_reactable.R tests/testthat/
git commit -m "feat: label locally searched units and stop mislabelling them as seed hits"
```

---

### Task 6: Tell the user when the panel could not answer

The screen already reaches `verdict = "no_reference"` with `reason_code = "no_hits"` for a
unit with nothing stored, and already writes `screen.screen_reason`. Nothing renders it.

**Files:**
- Modify: `R/screen.R` (reason text only), `R/app_utils.R` (`unit_evidence`),
  `R/app_server.R` (the scaffolds reactable)
- Test: `tests/testthat/test-screen.R`, `tests/testthat/test-app-utils.R` (extend)

**Interfaces:**
- Consumes: `blast_searched.source` (Task 2).
- Produces: `screen.screen_reason` mentions the local panel and suggests a remote search
  when the unanswered unit was searched locally. `unit_evidence()` returns `screen_reason`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-screen.R`, using `mk_screen_db()` already defined there:

```r
test_that("an unanswered local search suggests a remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 0L, source = "local")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM screen")
  expect_equal(r$screen_verdict[1], "no_reference")
  expect_match(r$screen_reason[1], "local")
  expect_match(r$screen_reason[1], "remote")
})

test_that("an unanswered remote search does not suggest another remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 0L, source = "remote")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  expect_false(grepl("remote search",
    DBI::dbGetQuery(con, "SELECT screen_reason FROM screen")$screen_reason[1]))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: the first FAILS on the reason text.

- [ ] **Step 3: Extend the reason text**

In `R/screen.R`, where the zero-hit short circuit sets its reason, look up
`blast_searched.source` for the unit. When it is `"local"`, write a reason naming the local
panel and suggesting a remote search, for example: "the local reference panel returned no
usable hit for this scaffold; a remote BLAST search may identify it". Otherwise keep the
current wording. Wrap the lookup so a project without the table behaves as before.

- [ ] **Step 4: Surface it in the app**

Add `screen_reason` to the columns `unit_evidence()` selects (`R/app_utils.R:156-189`).

In the scaffolds reactable (`R/app_server.R:1043-1127`), render it next to
`screen_verdict`, using `rt_longtext()` as the `triage_reason` column already does at
`R/app_server.R:1118-1120`. Keep it narrow; it is an advisory, not a primary column.

- [ ] **Step 5: Run tests, then the full suite, then commit**

```bash
git add R/screen.R R/app_utils.R R/app_server.R tests/testthat/
git commit -m "feat: surface why a unit could not be identified"
```

---

### Task 7: Acceptance run on an invertebrate clade

The one remaining gate. Not a code change; a measurement that decides whether the local
path is useful on the clades it was built for.

**Files:** none. Report only.

- [ ] **Step 1: Pick a project and record its baseline**

Use one of `/home/dmacguig/ribopilot_validation/{asteroidea,scyphozoa,keratoisis}`. These
predate the BLAST schema and have no `blast_hits` table, so work on a COPY and note what
tables exist.

- [ ] **Step 2: Run the units against the local panel**

For each active unit, run the composed local search by hand with the built panel at
`/home/dmacguig/ribopilot_refdb_build/blastdb/2026.08.1/ribopilot_rrna`, using the same
`-outfmt`, `-task megablast`, `-max_target_seqs` and `-max_hsps` the module uses, and parse
the results with `parse_blast_candidates()`.

- [ ] **Step 3: Report the numbers**

For that clade report: how many units got any hit; how many got a reference candidate that
clears the existing floors; the distribution of top `pident`; how many units the panel could
not answer and by which trigger (zero hits, no candidate, or top identity below the clade's
`min_anchor_pid["18S"]`, which for an invertebrate ruleset is 93); and local query wall time
per unit.

Report the numbers honestly whatever they are. A high unanswerable rate is a legitimate
result and is precisely what this task exists to discover. Do NOT adjust thresholds, swap
projects, or subset units to improve the figure. If the panel cannot answer most
invertebrate units, that is the finding, and it belongs in the design doc.

- [ ] **Step 4: Record the result in the design doc**

Add the measurement to the Stage 4 section of `dev/design/17_REMOTE_BLAST_PERFORMANCE.md`,
including the date, the clade, and the unanswerable rate. Commit.

---

## After this plan

Stage 4a is complete when tasks 1 to 6 are committed with no new suite failures and task 7's
measurement is recorded.

Not in scope:

- **Stage 4b, distribution.** Fetching the decoys at build time, the Docker layer strategy,
  the GitHub release, the manifest, the download helper, and recording which database
  version produced a project's results. Until 4b lands, a user must build or place the panel
  themselves and set `reference_config.blastdb`.
- Still open from earlier stages: the API key reaches `.command.run`; `REF_ALN_MIN_COV`
  duplicates and silently caps the tunable `min_ref_coverage`; `compute_feature_msa` has no
  coverage floors while `compute_ref_alignment` does; the single-record fetch path is dead
  but retained.
