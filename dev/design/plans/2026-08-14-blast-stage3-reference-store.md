# BLAST Stage 3: the shared reference store

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Keep fetched and locally-annotated GenBank references in one store outside any
project, so a second project in the same clade neither re-fetches nor re-annotates them.

**Architecture:** Reference records are project-independent: the same accession means the
same sequence and the same element coordinates whoever asks. Today each project fetches
and annotates its own copy, paying NCBI round trips and about 200 s of ITSx per project
for records another project already has. This adds a SQLite store under the user-chosen
reference folder, read before fetching and written after. The read-through needs no change
to the fetch logic: the fetch process already selects only candidate accessions absent
from `blast_ref_sequences`, so priming the project database from the store makes those
accessions drop out of the fetch list on their own.

**Tech Stack:** R (DBI, RSQLite, Biostrings), SQLite, Nextflow DSL2, testthat 3e.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md` (Stage 3)

## Global Constraints

- R >= 4.4.0 (`DESCRIPTION`); `docker/Dockerfile` pins `r-base=4.4.*`. No package outside
  DESCRIPTION Imports.
- ASCII only in R and Nextflow sources. No non-ASCII characters, no em dashes.
- Minimal comments. Do not narrate fixes in comments.
- This project ships no schema migrations before release. The store is a NEW file, so it
  has no migration burden of its own; only additions to the project database do.
- Never `git push`. Never `git pull`, `git fetch`, `git merge`, or `git rebase`: the remote
  branch is deliberately stale relative to this local branch and pulling corrupts history.
- Run tests with:
  `Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/<file>", reporter="summary")'`
- Run every command in the FOREGROUND. Do not background anything and do not set up
  monitors; the full suite takes a few minutes, just let it block.
- Full suite baseline entering this plan: 2465 passed, 0 failed, 0 errors, 6 skipped.
- The Perl `barrnap` in the local conda env is broken (missing `Path::Tiny`). `pybarrnap`
  0.5.2 and `ITSx` 1.1.3 work at `/home/dmacguig/miniconda3/envs/ribopilot/bin`. Prefer the
  fake-tool helpers in `tests/testthat/test-rrna-kingdom.R` over real binaries.

## Facts this plan rests on

- `open_project_db()` (`R/pipeline_steps.R:170-175`) sets `busy_timeout=10000` then
  `journal_mode` from `sqlite_journal_mode()` (`R/pipeline_steps.R:151-161`), which already
  returns `delete` on a network filesystem and `wal` otherwise, honouring a
  `RIBOPILOT_SQLITE_JOURNAL` override. Reuse it; do not invent a second policy.
- **There is no read-only open path anywhere in the package.** A `journal_mode` pragma
  takes an exclusive lock, so it fails on a read-only file, and because every caller wraps
  the connection in `tryCatch`, the failure is swallowed and the store silently contributes
  nothing. A group-shared store that one administrator populates and everyone else reads is
  the recommended HPC deployment, so this must work.
- The fetch query at `inst/nextflow/modules/blast_ref.nf:22` already reads
  `... LEFT JOIN blast_ref_sequences s ON s.accession = c.accession WHERE ... s.accession IS NULL`.
  Anything present in the project's `blast_ref_sequences` is therefore not fetched.
- `reference_mount_dirs()` (`R/reference_nf.R:96-113`) plus `container_bind_opts()`
  (`R/reference_nf.R:57-66`) emit a per-run `-process.containerOptions` override appended by
  `nf_command()` (`R/nf_command.R:51-52`). That is the mechanism for making a directory
  visible in-container WITHOUT rewriting a saved project config, which the static
  `<<EXTRA_MOUNTS>>` path in `R/init_project.R:117-148` cannot do for an existing project.
- `params.taxdump_dir` (`R/nf_command.R:46`, consumed at `inst/nextflow/modules/blast.nf:94`)
  is the closest existing analogue for plumbing a filesystem path into a process.
- `write_blast_refs()` (`R/blast_ref_fetch.R:121-154`) persists only
  `element, type, start, stop, strand`; `annotate_records()` also returns `tool` and
  `model_cov`, which are dropped. `blast_ref_sequences.taxid` is always written `NA`.
- No software version string is captured anywhere in the package.

---

### Task 1: Close Stage 2's parked regression hole

Carried over. Stage 2 fixed a defect where ITSx reported sense-frame coordinates regardless
of input orientation, by passing `--complement F` from `annotate_records()`. That fix has no
test: the fake ITSx in `tests/testthat/test-annotate-multi.R` ignores every flag but `-i`
and `-o`, so reverting the fix passes the whole suite silently.

`tests/testthat/test-rrna-kingdom.R:188-195` already shows the pattern, asserting on the
`-t` argv token via `mk_fake_tool()` and `fake_args()`.

**Files:**
- Test: `tests/testthat/test-rrna-kingdom.R` (extend)

**Interfaces:**
- Consumes: `run_itsx(fasta, cpu, taxon, complement)` and `annotate_records()`.
- Produces: nothing.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-rrna-kingdom.R`:

```r
test_that("run_itsx omits --complement unless asked", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa))
  expect_false(grepl("--complement", fake_args(d, "ITSx")))
})

test_that("run_itsx passes --complement F when asked", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa, complement = FALSE))
  expect_match(fake_args(d, "ITSx"), "--complement F")
})

test_that("annotate_records asks ITSx not to search the complement", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "pybarrnap", status = 0L)
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- file.path(d, "refs.fa")
  writeLines(c(">ACC1.1 x", paste(rep("ACGT", 60), collapse = "")), fa)
  suppressWarnings(RiboPilot:::annotate_records(fa))
  expect_match(fake_args(d, "ITSx"), "--complement F")
})
```

The third test is the one that matters: it pins the caller, not just the wrapper. If
`annotate_records()` runs ITSx more than once, `fake_args()` may hold only the last
invocation; that is fine, every invocation must carry the flag.

- [ ] **Step 2: Run tests to verify they fail against a reverted fix**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-rrna-kingdom.R", reporter="summary")'`

All three should PASS immediately, because the fix is already in place. That is expected and
is NOT sufficient. Prove they bite: temporarily change `complement = FALSE` to
`complement = NULL` at the `run_on()` call inside `annotate_records()` in `R/annotate.R`,
re-run, and confirm the THIRD test fails. Then restore the file exactly and confirm all
three pass again. Report both runs.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-rrna-kingdom.R
git commit -m "test: pin the ITSx complement flag so the fix cannot silently revert"
```

---

### Task 2: Open a SQLite database read-only without taking a write lock

The store's whole point is that many projects, and on a cluster many users, read one copy.
A `journal_mode` pragma takes an exclusive lock and fails outright on a read-only file, so
today a shared read-only store would fail at connect and every read with it.

**Files:**
- Modify: `R/pipeline_steps.R` (append near `open_project_db()`)
- Test: `tests/testthat/test-ref-store.R` (create)

**Interfaces:**
- Consumes: `sqlite_journal_mode()` (`R/pipeline_steps.R:151`).
- Produces:
  `open_sqlite_ro(path)` returning a DBI connection opened read-only, with
  `busy_timeout` set and NO `journal_mode` pragma. Errors if the file does not exist.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-ref-store.R`:

```r
# The shared store is read by many projects and, on a cluster, by many users.
# Reading must not require write access.

mk_ro_db <- function(dir) {
  p <- file.path(dir, "ro.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "CREATE TABLE t (a TEXT)")
  DBI::dbExecute(con, "INSERT INTO t VALUES ('x')")
  DBI::dbDisconnect(con)
  Sys.chmod(p, "0444")
  p
}

test_that("open_sqlite_ro reads a write-protected database", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT a FROM t")$a, "x")
})

test_that("open_sqlite_ro refuses to write", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_error(DBI::dbExecute(con, "INSERT INTO t VALUES ('y')"))
})

test_that("open_sqlite_ro errors on a missing file rather than creating one", {
  d <- withr::local_tempdir()
  p <- file.path(d, "absent.sqlite")
  expect_error(RiboPilot:::open_sqlite_ro(p))
  expect_false(file.exists(p))
})
```

`Sys.chmod` is a no-op when the test runs as root. If any of these fail because the process
can write regardless of mode, report that rather than weakening the test; a root-only
failure is worth knowing about.

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL, `open_sqlite_ro` not found.

- [ ] **Step 3: Implement**

Append to `R/pipeline_steps.R`:

```r
#' Open a SQLite database read-only
#'
#' The shared reference store may live on a path the running user cannot write,
#' which is the recommended cluster deployment. `journal_mode` takes an
#' exclusive lock and fails there, so this deliberately sets only busy_timeout.
#'
#' @param path Path to an existing SQLite file.
#' @return An open read-only DBI connection.
#' @noRd
open_sqlite_ro <- function(path) {
  if (!file.exists(path)) {
    stop("no such database: ", path, call. = FALSE)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
  con
}
```

- [ ] **Step 4: Run tests, then the full suite, then commit**

```bash
git add R/pipeline_steps.R tests/testthat/test-ref-store.R
git commit -m "feat: read-only SQLite open for a shared store"
```

---

### Task 3: The store itself

**Files:**
- Create: `R/ref_store.R`
- Test: `tests/testthat/test-ref-store.R` (extend)

**Interfaces:**
- Consumes: `open_sqlite_ro()` (Task 2), `sqlite_journal_mode()`.
- Produces:
  - `ref_store_path(ref_dir)` returning `<ref_dir>/refstore/refstore.sqlite`.
  - `REF_STORE_SEQ_DDL`, `REF_STORE_ANN_DDL` string constants.
  - `open_ref_store(ref_dir, write = FALSE)` returning a connection, creating the file and
    tables when `write = TRUE`, opening read-only when `write = FALSE`, and returning NULL
    when the store does not exist and `write` is FALSE.
  - `ref_store_get(store, accessions)` returning
    `list(seqs = data.frame, ann = data.frame)`, empty frames when nothing matches.
  - `ref_store_put(store, seqs, ann, kingdom, taxon, tool)` returning, invisibly, the number
    of accessions written.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ref-store.R`:

```r
mk_store_dir <- function() withr::local_tempdir(.local_envir = parent.frame())

seq_fixture <- function() data.frame(
  accession = c("A.1", "B.1"), species = c("Genus a", "Genus b"),
  title = c("Genus a 18S", "Genus b 28S"), sequence = c("ACGT", "TTTT"),
  length = c(4L, 4L), stringsAsFactors = FALSE)

ann_fixture <- function() data.frame(
  accession = c("A.1", "A.1", "B.1"), element = c("18S", "ITS1", "28S"),
  type = c("rRNA", "spacer", "rRNA"), start = c(1L, 3L, 1L), stop = c(2L, 4L, 4L),
  strand = c("+", "+", "-"), stringsAsFactors = FALSE)

test_that("a store round-trips sequences and annotations", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  n <- RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(),
    kingdom = "euk", taxon = "all", tool = "pybarrnap")
  expect_equal(as.integer(n), 2L)
  got <- RiboPilot:::ref_store_get(st, c("A.1", "B.1"))
  expect_equal(sort(got$seqs$accession), c("A.1", "B.1"))
  expect_equal(nrow(got$ann), 3L)
  expect_equal(sort(got$ann$element[got$ann$accession == "A.1"]), c("18S", "ITS1"))
})

test_that("ref_store_get returns only what was asked for", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  got <- RiboPilot:::ref_store_get(st, "A.1")
  expect_equal(got$seqs$accession, "A.1")
  expect_true(all(got$ann$accession == "A.1"))
})

test_that("ref_store_get on an unknown accession returns empty frames, not an error", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  got <- RiboPilot:::ref_store_get(st, "ZZ.9")
  expect_equal(nrow(got$seqs), 0L)
  expect_equal(nrow(got$ann), 0L)
})

test_that("re-putting an accession replaces rather than duplicating", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_sequences")$n, 2L)
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_annotations")$n, 3L)
})

test_that("the store records how its coordinates were produced", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "M", "pybarrnap")
  r <- DBI::dbGetQuery(st, "SELECT kingdom, taxon, tool FROM ref_sequences WHERE accession='A.1'")
  expect_equal(r$kingdom, "euk")
  expect_equal(r$taxon, "M")
  expect_equal(r$tool, "pybarrnap")
})

test_that("open_ref_store returns NULL for a store that does not exist", {
  d <- mk_store_dir()
  expect_null(RiboPilot:::open_ref_store(d, write = FALSE))
})

test_that("an existing store opens read-only without a write lock", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)
  Sys.chmod(RiboPilot:::ref_store_path(d), "0444")
  ro <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(ro))
  expect_equal(nrow(RiboPilot:::ref_store_get(ro, "A.1")$seqs), 1L)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, `open_ref_store` not found.

- [ ] **Step 3: Implement**

Create `R/ref_store.R`:

```r
# Cross-project store of fetched reference sequences and their locally derived
# element coordinates. Lives under the user-chosen reference folder so a second
# project in the same clade neither re-fetches nor re-annotates.

#' Path to the store inside a reference folder
#' @noRd
ref_store_path <- function(ref_dir) {
  file.path(path.expand(ref_dir), "refstore", "refstore.sqlite")
}

#' @noRd
REF_STORE_SEQ_DDL <- "CREATE TABLE IF NOT EXISTS ref_sequences (
  accession TEXT NOT NULL,
  species TEXT,
  title TEXT,
  sequence TEXT,
  length INTEGER,
  kingdom TEXT,
  taxon TEXT,
  tool TEXT,
  stored_utc TEXT,
  PRIMARY KEY (accession)
);"

#' @noRd
REF_STORE_ANN_DDL <- "CREATE TABLE IF NOT EXISTS ref_annotations (
  accession TEXT NOT NULL,
  element TEXT NOT NULL,
  type TEXT,
  start INTEGER,
  stop INTEGER,
  strand TEXT,
  PRIMARY KEY (accession, element)
);"

#' Open the reference store under a reference folder
#'
#' @param ref_dir Reference folder, as resolved by [ref_config()].
#' @param write TRUE to create the store if absent and open it writable.
#' @return A DBI connection, or NULL when the store does not exist and `write`
#'   is FALSE.
#' @noRd
open_ref_store <- function(ref_dir, write = FALSE) {
  p <- ref_store_path(ref_dir)
  if (!write) {
    if (!file.exists(p)) return(NULL)
    return(open_sqlite_ro(p))
  }
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
  DBI::dbExecute(con, sprintf("PRAGMA journal_mode=%s;", sqlite_journal_mode(p)))
  DBI::dbExecute(con, REF_STORE_SEQ_DDL)
  DBI::dbExecute(con, REF_STORE_ANN_DDL)
  con
}

#' Fetch stored records for a set of accessions
#' @return list(seqs, ann); zero-row frames when nothing matches.
#' @noRd
ref_store_get <- function(store, accessions) {
  empty <- list(seqs = data.frame(), ann = data.frame())
  if (is.null(store) || length(accessions) == 0) return(empty)
  ph <- paste(rep("?", length(accessions)), collapse = ",")
  seqs <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_sequences WHERE accession IN (%s)", ph),
    params = as.list(accessions))
  if (nrow(seqs) == 0) return(empty)
  ann <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_annotations WHERE accession IN (%s)", ph),
    params = as.list(accessions))
  list(seqs = seqs, ann = ann)
}

#' Write records into the store
#'
#' One transaction per accession, so a failure partway leaves earlier accessions
#' stored and a concurrent writer contends for a short lock rather than a long
#' one.
#' @return (invisibly) the number of accessions written.
#' @noRd
ref_store_put <- function(store, seqs, ann, kingdom, taxon, tool) {
  if (is.null(store) || is.null(seqs) || nrow(seqs) == 0) return(invisible(0L))
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  n <- 0L
  for (i in seq_len(nrow(seqs))) {
    acc <- seqs$accession[i]
    rows <- if (!is.null(ann) && nrow(ann) > 0) ann[ann$accession == acc, , drop = FALSE] else NULL
    DBI::dbWithTransaction(store, {
      DBI::dbExecute(store, "DELETE FROM ref_sequences WHERE accession = ?", params = list(acc))
      DBI::dbAppendTable(store, "ref_sequences", data.frame(
        accession = acc, species = seqs$species[i], title = seqs$title[i],
        sequence = seqs$sequence[i], length = seqs$length[i],
        kingdom = kingdom, taxon = taxon, tool = tool, stored_utc = ts,
        stringsAsFactors = FALSE))
      DBI::dbExecute(store, "DELETE FROM ref_annotations WHERE accession = ?", params = list(acc))
      if (!is.null(rows) && nrow(rows) > 0) {
        DBI::dbAppendTable(store, "ref_annotations", data.frame(
          accession = acc, element = rows$element, type = rows$type,
          start = rows$start, stop = rows$stop, strand = rows$strand,
          stringsAsFactors = FALSE))
      }
    })
    n <- n + 1L
  }
  invisible(n)
}
```

- [ ] **Step 4: Run tests, then the full suite, then commit**

```bash
git add R/ref_store.R tests/testthat/test-ref-store.R
git commit -m "feat: cross-project reference store"
```

---

### Task 4: Prime a project from the store before fetching

The fetch process already excludes candidate accessions that are present in the project's
`blast_ref_sequences`. So copying store hits into the project database is the whole
read-through: the fetch list shrinks by itself.

**Files:**
- Modify: `R/ref_store.R` (append)
- Create: `inst/nextflow/bin/rp_ref_prime.R`
- Test: `tests/testthat/test-ref-store.R` (extend)

**Interfaces:**
- Consumes: `open_ref_store()`, `ref_store_get()`.
- Produces:
  `prime_from_ref_store(con, ref_dir, ids = NULL)` returning, invisibly, the number of
  accessions copied into the project database. A missing store is not an error.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ref-store.R`:

```r
mk_project_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT,
    rank INT, accession TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences(accession TEXT PRIMARY KEY,
    species TEXT, title TEXT, taxid TEXT, lineage TEXT, sequence TEXT, length INT,
    time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations(accession TEXT, element TEXT,
    type TEXT, start INT, stop INT, strand TEXT, time_stamp INT,
    PRIMARY KEY (accession, element))")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,1,'A.1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,2,'ZZ.9')")
  con
}

test_that("priming copies a stored accession into the project and skips unknown ones", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::prime_from_ref_store(con, d)
  expect_equal(as.integer(n), 1L)
  s <- DBI::dbGetQuery(con, "SELECT accession FROM blast_ref_sequences")
  expect_equal(s$accession, "A.1")
  a <- DBI::dbGetQuery(con, "SELECT element FROM blast_ref_annotations ORDER BY element")
  expect_equal(a$element, c("18S", "ITS1"))
})

test_that("priming is a no-op when no store exists", {
  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  d <- mk_store_dir()
  expect_equal(as.integer(RiboPilot:::prime_from_ref_store(con, d)), 0L)
})

test_that("priming does not overwrite an accession the project already has", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence) VALUES ('A.1','KEEPME')")
  RiboPilot:::prime_from_ref_store(con, d)
  expect_equal(DBI::dbGetQuery(con, "SELECT sequence FROM blast_ref_sequences WHERE accession='A.1'")$sequence,
    "KEEPME")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, `prime_from_ref_store` not found.

- [ ] **Step 3: Implement**

Append to `R/ref_store.R`:

```r
#' Copy stored references for a project's candidates into the project database
#'
#' The fetch step selects only candidate accessions absent from
#' `blast_ref_sequences`, so anything copied here is not fetched again. An
#' accession the project already holds is left alone.
#'
#' @param con Open project DB connection.
#' @param ref_dir Reference folder, as resolved by [ref_config()].
#' @param ids Optional sample IDs to limit to; default is every candidate.
#' @return (invisibly) the number of accessions copied.
#' @export
prime_from_ref_store <- function(con, ref_dir, ids = NULL) {
  store <- tryCatch(open_ref_store(ref_dir, write = FALSE), error = function(e) NULL)
  if (is.null(store)) return(invisible(0L))
  on.exit(DBI::dbDisconnect(store), add = TRUE)

  sql <- "SELECT DISTINCT c.accession FROM blast_ref_candidates c
          LEFT JOIN blast_ref_sequences s ON s.accession = c.accession
          WHERE s.accession IS NULL"
  want <- if (is.null(ids)) {
    DBI::dbGetQuery(con, sql)$accession
  } else {
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbGetQuery(con, sprintf("%s AND c.ID IN (%s)", sql, ph), params = as.list(ids))$accession
  }
  if (length(want) == 0) return(invisible(0L))

  got <- ref_store_get(store, want)
  if (nrow(got$seqs) == 0) return(invisible(0L))

  ts <- as.integer(Sys.time())
  n <- 0L
  for (i in seq_len(nrow(got$seqs))) {
    acc <- got$seqs$accession[i]
    rows <- if (nrow(got$ann) > 0) got$ann[got$ann$accession == acc, , drop = FALSE] else NULL
    DBI::dbWithTransaction(con, {
      DBI::dbAppendTable(con, "blast_ref_sequences", data.frame(
        accession = acc, species = got$seqs$species[i], title = got$seqs$title[i],
        taxid = NA_character_, lineage = NA_character_, sequence = got$seqs$sequence[i],
        length = got$seqs$length[i], time_stamp = ts, stringsAsFactors = FALSE))
      if (!is.null(rows) && nrow(rows) > 0) {
        DBI::dbAppendTable(con, "blast_ref_annotations", data.frame(
          accession = acc, element = rows$element, type = rows$type,
          start = rows$start, stop = rows$stop, strand = rows$strand,
          time_stamp = ts, stringsAsFactors = FALSE))
      }
    })
    n <- n + 1L
  }
  invisible(n)
}
```

- [ ] **Step 4: Add the pipeline entry point**

Create `inst/nextflow/bin/rp_ref_prime.R`, following `inst/nextflow/bin/rp_rederive.R`:

```r
#!/usr/bin/env Rscript
# Copy stored references into the project before fetching.
# Args: <db> <ref_dir> [id ...]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_prime.R <db> <ref_dir> [id ...]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  ids <- if (length(a) > 2) a[-(1:2)] else NULL
  n <- prime_from_ref_store(con, a[2], ids = ids)
  cat("primed from store:", n, "\n")
})
```

`chmod +x inst/nextflow/bin/rp_ref_prime.R`. Check `tests/testthat/test-nf-bin-contract.R`
and satisfy every rule it enforces.

- [ ] **Step 5: Run tests, then the full suite, regenerate docs, commit**

`prime_from_ref_store` is `@export`, so run `Rscript -e 'roxygen2::roxygenise()'` and include
the regenerated NAMESPACE and `man/`.

```bash
git add R/ref_store.R inst/nextflow/bin/rp_ref_prime.R NAMESPACE man/ tests/testthat/test-ref-store.R
git commit -m "feat: prime a project from the shared reference store"
```

---

### Task 5: Write through to the store after fetching

**Files:**
- Modify: `R/ref_store.R` (append), `inst/nextflow/bin/rp_ref_store.R`
- Test: `tests/testthat/test-ref-store.R` (extend)

**Interfaces:**
- Consumes: `open_ref_store()`, `ref_store_put()`, `annotate_records()` output.
- Produces:
  `push_to_ref_store(ref_dir, fasta_file, ann, kingdom, taxon, tool)` returning, invisibly,
  the number of accessions written. A ref_dir that cannot be written is not an error: the
  run continues without a store.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-ref-store.R`:

```r
test_that("pushing writes fetched records into the store", {
  d <- mk_store_dir()
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 Genus a 18S ribosomal RNA", "ACGT",
               ">B.1 Genus b 28S ribosomal RNA", "TTTT"), fa)
  ann <- data.frame(seqid = c("A.1", "B.1"), element = c("18S", "28S"),
    type = c("rRNA", "rRNA"), start = c(1L, 1L), stop = c(4L, 4L),
    strand = c("+", "-"), stringsAsFactors = FALSE)
  n <- RiboPilot:::push_to_ref_store(d, fa, ann, "euk", "all", "pybarrnap")
  expect_equal(as.integer(n), 2L)

  st <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(st))
  got <- RiboPilot:::ref_store_get(st, c("A.1", "B.1"))
  expect_equal(sort(got$seqs$accession), c("A.1", "B.1"))
  expect_equal(got$seqs$species[got$seqs$accession == "A.1"], "Genus a")
  expect_equal(nrow(got$ann), 2L)
})

test_that("pushing to an unwritable reference folder is not an error", {
  d <- mk_store_dir()
  Sys.chmod(d, "0555")
  on.exit(Sys.chmod(d, "0755"), add = TRUE)
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "ACGT"), fa)
  expect_equal(as.integer(RiboPilot:::push_to_ref_store(d, fa, NULL, "euk", "all", "pybarrnap")), 0L)
})
```

If the second test's `Sys.chmod` has no effect because the process is root, report that
rather than deleting the test.

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, `push_to_ref_store` not found.

- [ ] **Step 3: Implement**

Append to `R/ref_store.R`:

```r
#' Write fetched records and their annotations into the shared store
#'
#' Best effort. A reference folder that cannot be written leaves the run
#' unaffected: the project keeps its own copy either way.
#'
#' @param ref_dir Reference folder.
#' @param fasta_file The multi-record fasta that was fetched.
#' @param ann annotate_records() output, or NULL.
#' @param kingdom,taxon,tool Provenance for the coordinates.
#' @return (invisibly) the number of accessions written.
#' @export
push_to_ref_store <- function(ref_dir, fasta_file, ann, kingdom, taxon, tool) {
  store <- tryCatch(open_ref_store(ref_dir, write = TRUE), error = function(e) NULL)
  if (is.null(store)) return(invisible(0L))
  on.exit(DBI::dbDisconnect(store), add = TRUE)

  seqs <- tryCatch(Biostrings::readDNAStringSet(fasta_file), error = function(e) NULL)
  if (is.null(seqs) || length(seqs) == 0) return(invisible(0L))
  full <- names(seqs)
  acc <- sub("\\s.*$", "", full)
  title <- trimws(sub("^\\S+\\s*", "", full))
  sdf <- data.frame(accession = acc, species = extract_species(title), title = title,
    sequence = as.character(seqs), length = as.integer(Biostrings::width(seqs)),
    stringsAsFactors = FALSE)

  adf <- NULL
  if (!is.null(ann) && nrow(ann) > 0) {
    adf <- data.frame(accession = ann$seqid, element = ann$element, type = ann$type,
      start = ann$start, stop = ann$stop, strand = ann$strand, stringsAsFactors = FALSE)
  }
  tryCatch(ref_store_put(store, sdf, adf, kingdom, taxon, tool),
    error = function(e) invisible(0L))
}
```

- [ ] **Step 4: Call it from the store script**

In `inst/nextflow/bin/rp_ref_store.R`, after `write_blast_refs()` and the
`stop_if_nothing_stored()` guard, push to the store when a reference folder was supplied.
Add an optional trailing `ref_dir` argument; when it is absent or empty, skip the push
entirely. Report the argument order you chose.

- [ ] **Step 5: Run tests, then the full suite, regenerate docs, commit**

```bash
git add R/ref_store.R inst/nextflow/bin/rp_ref_store.R NAMESPACE man/ tests/testthat/test-ref-store.R
git commit -m "feat: write fetched references through to the shared store"
```

---

### Task 6: Plumb the reference folder into the pipeline

**Files:**
- Modify: `R/nf_command.R`, `inst/nextflow/main.nf`, `inst/nextflow/modules/blast_ref.nf`
- Test: `tests/testthat/test-nextflow.R`, `tests/testthat/test-nf-command-quoting.R` (extend)

**Interfaces:**
- Consumes: `ref_config()` (`R/reference_dir.R:15`), `prime_from_ref_store()`,
  `push_to_ref_store()`.
- Produces: `params.ref_dir`, defaulting to the empty string, consumed by
  `blast_ref_fetch`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-nextflow.R`:

```r
test_that("main.nf declares a ref_dir param", {
  p <- system.file("nextflow", "main.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/main.nf"
  expect_true(any(grepl("params.ref_dir", readLines(p, warn = FALSE), fixed = TRUE)))
})

test_that("the fetch module primes from the store before fetching", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- readLines(p, warn = FALSE)
  prime <- grep("rp_ref_prime.R", txt)
  accs <- grep("SELECT DISTINCT c.accession", txt)
  expect_length(prime, 1L)
  expect_true(prime < accs[1],
    info = "priming must run before the un-fetched accession list is computed")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: both FAIL.

- [ ] **Step 3: Pass the folder from R**

In `R/nf_command.R`, resolve the reference folder the same way the rest of the package does
and pass it as a param, alongside the existing `--taxdump_dir`. Read how `taxdump` is
resolved and follow that shape. The value comes from `ref_config(con)$ref_dir` on the
project database; if that cannot be read, pass an empty string rather than failing.

The folder must also be visible inside the container. `reference_mount_dirs()`
(`R/reference_nf.R:96-113`) already computes per-run bind mounts and drops directories that
are already covered; add the reference folder to what it considers. Do NOT touch the static
`<<EXTRA_MOUNTS>>` path in `R/init_project.R`, because that is baked in at project creation
and cannot help an existing project.

- [ ] **Step 4: Declare and consume the param**

In `inst/nextflow/main.nf`, declare `params.ref_dir = ''` beside the other param defaults.

In `inst/nextflow/modules/blast_ref.nf`, at the very top of the `blast_ref_fetch` script,
before the accession list is computed, prime from the store when a folder was given:

```bash
    if [ -n '${params.ref_dir}' ]; then
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_prime.R '${params.db}' '${params.ref_dir}' '${id}' || true
    fi
```

The `|| true` is deliberate: a store that cannot be read must not fail the run, because the
fetch below is a complete fallback.

Then pass the folder to `rp_ref_store.R` so the write-through happens, matching the argument
order you chose in Task 5.

- [ ] **Step 5: Verify end to end on real data**

This is the acceptance check for the whole stage. Work on COPIES.

Take the reference sequences already in
`/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/e2e_run2/.sqlite`, push them into a
store in a temporary reference folder together with annotations from `annotate_records()`,
then take a COPY of `/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/full_project_docker/.sqlite`,
empty its `blast_ref_sequences` and `blast_ref_annotations`, and run
`prime_from_ref_store()` against that store.

Report: how many of the project's candidate accessions the store could satisfy, how many it
could not, and how long priming took. The two projects are different fish datasets, so
overlap may be low or zero. A low number is a legitimate result and is exactly what the
design predicted for within-clade overlap; report it honestly rather than choosing projects
to make it look better.

- [ ] **Step 6: Run the full suite and commit**

If `nextflow lint` is available, run it on both changed `.nf` files.

```bash
git add R/nf_command.R inst/nextflow/main.nf inst/nextflow/modules/blast_ref.nf tests/testthat/test-nextflow.R
git commit -m "feat: plumb the reference folder into the fetch step"
```

---

## After this plan

Stage 3 is complete when all six tasks are committed and the full suite shows no new
failures.

Deliberately not in scope:

- **Migrating an existing project's references into a store.** A one-shot importer would be
  useful and is a small function over the same primitives, but nothing depends on it.
- **A store manifest or any human-readable index.** `append_ref_manifest()`
  (`R/reference_dir.R:96-117`) is a full read-modify-write with no locking and no atomic
  rename, so two concurrent users can silently lose each other's rows. The store
  deliberately uses SQLite instead, which handles concurrency through `busy_timeout` and the
  journal mode `sqlite_journal_mode()` already selects per filesystem. Do not add a TSV
  index alongside it.
- **Stage 4**: the shipped local rDNA BLAST database, the `blast_local` and
  `blast_remote_escalate` split, and the five escalation triggers. Stage 4 must not ship
  before its acceptance gate is measured on Asteroidea, Scyphozoa, or Octocorallia: local
  query time and peak memory, escalation rate, and local-versus-remote concordance.
- Still open from earlier stages: the API key reaches `.command.run`;
  `REF_ALN_MIN_COV` duplicates and silently caps the tunable `min_ref_coverage`;
  `compute_feature_msa` has no coverage floors while `compute_ref_alignment` does; the
  single-record fetch path is dead but retained.
