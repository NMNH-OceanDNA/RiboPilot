# BLAST Stage 1: correctness without NCBI traffic

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Separate BLAST result *interpretation* from BLAST *searching*, so candidate
selection and the contamination screen become re-derivable from stored hits, then fix
the query-coverage floor that rejects the most complete assemblies.

**Architecture:** `parse_blast_candidates()` currently parses a hits file and writes
`blast_hits`, `blast_ref_candidates`, and `blast_screen` in one pass, applying the
candidate filter at parse time. Because the filter runs during parsing, changing any
threshold costs a full remote re-search at roughly two hours per unit. `blast_hits`
already stores every column the filter needs, so this plan extracts the derivation into
a function over the table. Every later threshold change then costs seconds, and the
`qcovs` fix applies retroactively to existing projects with no network access.

**Tech Stack:** R (DBI, RSQLite), testthat 3e, Nextflow DSL2, sqlite3 CLI.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md` (Stage 1)

## Global Constraints

- R >= 4.4.0 (`DESCRIPTION`); `docker/Dockerfile` pins `r-base=4.4.*`. Do not use
  `tools::sha256sum()` (R 4.5.0) or any package not already in `Imports`.
- ASCII only in all R source. No non-ASCII characters.
- Minimal comments. Do not narrate bugfixes in comments.
- This project ships **no schema migrations** before release. Adding a column breaks
  projects created earlier, and that is accepted policy, not an oversight.
  `R/db_schema.R:147` is the only column-existence check in the codebase.
- Never `git push`.
- Run tests with:
  `Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/<file>", reporter="summary")'`
- Baseline: `tests/testthat/test-blast-parse.R` has 29 passing tests before this plan.
  It must still have at least 29 passing after every task.

---

### Task 1: Drop the write-only alignment strings

`blast_ref_alignment.sample_aln` and `.ref_aln` have no reader anywhere in `R/`,
`inst/`, or `tests/`. Verified: sole writer is `R/blast_ref_align.R:117`, sole schema is
`R/init_db.R:460-461`, and the only reader of the table at all is `R/curate_rrna.R:316`,
which selects `element, pid, aln_cov` by name. The export report builds its alignments
from `annotations` and `assemblies` with DECIPHER, not from this table.

Both columns are nullable and neither is in the primary key, so an old project's table
still accepts the narrower insert. This is the one task in the plan that does not break
existing projects.

**Files:**
- Modify: `R/blast_ref_align.R:111-117`
- Modify: `R/init_db.R:460-461`
- Test: `tests/testthat/test-ref-alignment.R` (create)

**Interfaces:**
- Consumes: nothing.
- Produces: `blast_ref_alignment` rows with columns
  `ID, path, scaffold, accession, element, pid, aln_cov, aln_len, time_stamp`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-ref-alignment.R`. This asserts behaviour rather than schema
text: the table is built without the two columns, so the old writer fails on
`dbAppendTable` with "table blast_ref_alignment has no column named sample_aln".

```r
# compute_ref_alignment() stores per-element identity and coverage, not the
# aligned strings themselves.

mk_aln_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_alignment (ID TEXT, path INT, scaffold INT,
    accession TEXT, element TEXT, pid REAL, aln_cov REAL, aln_len INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
    element TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT,
    rank INT, accession TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences(accession TEXT, sequence TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations(accession TEXT, element TEXT,
    start INT, stop INT, strand TEXT)")

  s <- paste(rep("ACGTACGTAC", 20), collapse = "")   # 200 bp
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S',1,200)")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,1,'ACC1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences VALUES('ACC1',?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations VALUES('ACC1','18S',1,200,'+')")
  con
}

test_that("compute_ref_alignment writes identity and coverage but no aligned strings", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 1L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_alignment")
  expect_equal(nrow(r), 1L)
  expect_equal(r$element[1], "18S")
  expect_equal(r$pid[1], 100)
  expect_equal(r$aln_len[1], 200L)
  expect_false("sample_aln" %in% names(r))
})
```

Separately, confirm by eye that `R/init_db.R` no longer declares the two columns; there is
no cheap programmatic check because that DDL is inline rather than a named constant.

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-ref-alignment.R", reporter="summary")'`

Expected: FAIL. The current writer builds a frame containing `sample_aln` and `ref_aln`,
so `dbAppendTable` errors with "table blast_ref_alignment has no column named sample_aln".

- [ ] **Step 3: Remove the columns from the schema**

In `R/init_db.R`, delete these two lines from the `blast_ref_alignment` DDL:

```r
      sample_aln TEXT,
      ref_aln TEXT,
```

- [ ] **Step 4: Remove the writer**

In `R/blast_ref_align.R`, delete line 112 entirely:

```r
    r_aln <- as.character(pwalign::subject(aln))
```

`r_aln` feeds storage only. Keep line 111 (`s_aln`): both surviving metrics derive from
it, `aln_cov` via `alignment_coverage(s_aln, nchar(s_sub))` and `aln_len` via
`nchar(s_aln)`.

Then remove `sample_aln` and `ref_aln` from the `dbAppendTable` frame so it reads:

```r
    DBI::dbAppendTable(con, "blast_ref_alignment", data.frame(ID = id, path = path, scaffold = scaffold,
      accession = acc, element = el, pid = round(pid, 1), aln_cov = round(aln_cov, 3), aln_len = nchar(s_aln),
      time_stamp = ts, stringsAsFactors = FALSE))
```

- [ ] **Step 5: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-ref-alignment.R", reporter="summary")'`
Expected: PASS.

Then confirm nothing else regressed:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-curate.R", reporter="summary")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add R/blast_ref_align.R R/init_db.R tests/testthat/test-ref-alignment.R
git commit -m "refactor: drop write-only alignment strings from blast_ref_alignment"
```

---

### Task 2: Extract candidate and screen derivation from the parser

Pure refactor. Behaviour must not change; the tests prove it. `parse_blast_candidates()`
keeps parsing the file and writing `blast_hits`, then delegates the two derived tables to
a new function that reads `blast_hits` instead of the in-memory frame.

Transaction boundaries change deliberately. Today one transaction covers all three
tables. After this task, `blast_hits` is written in its own transaction and the derived
products in a second. Each is atomic, and a failure in derivation now leaves `blast_hits`
intact and re-derivable, which is strictly better than losing the search.

**Files:**
- Modify: `R/blast_genbank.R:50-160`
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: `BLAST_HITS_DDL` (`R/db_schema.R:11`), `write_blast_screen()`
  (`R/blast_genbank.R:224`), `consensus_hit_class()` (`R/hit_class.R`).
- Produces:
  - `ref_candidate_ok(x)`: takes a data.frame with columns `hit_class`, `pident`,
    `qcovs`; returns a logical vector, one element per row. Task 3 modifies this.
  - `derive_blast_products(con, id, path = 1L, scaffold = 1L, ts = as.integer(Sys.time()))`
    reads `blast_hits` for one unit, replaces that unit's `blast_ref_candidates` and
    `blast_screen` rows, updates `assemblies.blast_*`. Returns the candidate data.frame
    invisibly, or an empty data.frame.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("derive_blast_products reproduces what the parser wrote", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(c(
    hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600),
    hit("AA2", "Danio rerio 28S ribosomal RNA gene", 98.0, 60, 900, 950, 1, 950)))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  before_c <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_candidates ORDER BY rank")
  before_s <- DBI::dbGetQuery(con, "SELECT * FROM blast_screen")

  RiboPilot:::derive_blast_products(con, "S1", 1, 1, ts = before_c$time_stamp[1])
  after_c <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_candidates ORDER BY rank")
  after_s <- DBI::dbGetQuery(con, "SELECT * FROM blast_screen")

  expect_equal(after_c, before_c)
  expect_equal(after_s, before_s)
})

test_that("derive_blast_products is idempotent and does not duplicate rows", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  RiboPilot:::derive_blast_products(con, "S1", 1, 1)
  RiboPilot:::derive_blast_products(con, "S1", 1, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 1L)
})

test_that("derive_blast_products on a unit with no hits clears its derived rows", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  DBI::dbExecute(con, "DELETE FROM blast_hits")
  RiboPilot:::derive_blast_products(con, "S1", 1, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 0L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: FAIL, `could not find function "derive_blast_products"`.

- [ ] **Step 3: Add the predicate and the derivation function**

In `R/blast_genbank.R`, immediately after the `REF_CANDIDATE_MIN_QCOVS` constant
(line 25), add:

```r
#' Does this subject qualify as a reference the sample is aligned to?
#'
#' @param x data.frame with hit_class, pident, qcovs.
#' @return logical vector, one element per row.
#' @noRd
ref_candidate_ok <- function(x) {
  !is.na(x$hit_class) & x$hit_class == "rdna" &
    !is.na(x$pident) & x$pident >= REF_CANDIDATE_MIN_PIDENT &
    !is.na(x$qcovs) & x$qcovs >= REF_CANDIDATE_MIN_QCOVS
}
```

Then add the derivation function after `parse_blast_candidates()`:

```r
#' Derive a unit's reference panel and contamination screen from stored hits
#'
#' Reads `blast_hits` rather than a hits file, so thresholds can change without
#' re-running a remote search. Replaces the unit's rows in
#' `blast_ref_candidates` and `blast_screen` and updates `assemblies.blast_*`.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param ts Timestamp to stamp the written rows with.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
derive_blast_products <- function(con, id, path = 1L, scaffold = 1L,
                                  ts = as.integer(Sys.time())) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  clear <- function() {
    for (tbl in c("blast_ref_candidates", "blast_screen")) {
      DBI::dbExecute(con,
        sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
        params = list(id, path, scaffold))
    }
  }
  h <- DBI::dbGetQuery(con,
    "SELECT * FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp",
    params = list(id, path, scaffold))
  if (nrow(h) == 0) {
    DBI::dbWithTransaction(con, clear())
    return(invisible(data.frame()))
  }

  subj <- h[!duplicated(h$hit_rank), , drop = FALSE]
  ref <- subj[ref_candidate_ok(subj), , drop = FALSE]
  df <- NULL
  if (nrow(ref) > 0) {
    ref$rank <- seq_len(nrow(ref))
    df <- data.frame(ID = id, path = path, scaffold = scaffold, rank = ref$rank,
      accession = ref$accession, species = ref$species, title = ref$title,
      pident = ref$pident, qcovs = ref$qcovs, evalue = ref$evalue,
      bitscore = ref$bitscore, time_stamp = ts, stringsAsFactors = FALSE)
  }

  DBI::dbWithTransaction(con, {
    clear()
    write_blast_screen(con, id, path, scaffold, subj[1, ], hits = h, ts = ts)
    if (!is.null(df)) {
      DBI::dbAppendTable(con, "blast_ref_candidates", df)
      top <- ref[1, ]
      DBI::dbExecute(con,
        "UPDATE assemblies SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?,
           blast_evalue = ?, blast_source = 'core_nt'
         WHERE ID = ? AND path = ? AND scaffold = ?",
        params = list(top$accession, top$species, top$pident, top$qcovs, top$evalue,
          id, path, scaffold))
    }
  })
  invisible(df %||% data.frame())
}
```

- [ ] **Step 4: Make the parser delegate**

In `parse_blast_candidates()`, replace everything from the comment
`# one row per subject, best HSP first, for the two derived tables` (line 116) through
the end of the function body with:

```r
  DBI::dbWithTransaction(con, {
    clear_unit()
    DBI::dbAppendTable(con, "blast_hits", hits)
  })
  derive_blast_products(con, id, path, scaffold, ts = ts)
}
```

Delete the now-unused `subj`, `ref`, and `df` construction from
`parse_blast_candidates()`. Keep `clear_unit()` as it is: it must still clear all three
tables, because a re-parse replaces the unit wholesale.

- [ ] **Step 5: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: PASS, at least 32 tests.

Then the screen and curation suites, which consume these tables:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); for (f in c("test-screen.R","test-screen-variants.R","test-curate.R","test-hit-class.R")) testthat::test_file(file.path("tests/testthat", f), reporter="summary")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add R/blast_genbank.R tests/testthat/test-blast-parse.R
git commit -m "refactor: derive candidates and screen from blast_hits, not the hits file"
```

---

### Task 3: Accept a subject on absolute aligned length as well as query coverage

`qcovs` for a complete rDNA assembly measures how much intergenic spacer the assembler
recovered, inverted, not how good the match is. GenBank cistron records stop at 28S, so a
15 kb full-repeat assembly matching a 7 kb reference end to end scores `qcovs` 46 and is
rejected, while a fragmented assembly scores 100. Measured: 10 of 26 units in
`full_project_docker` got zero candidates from this, all of them `SRR14433904` paths with
a 7208 bp alignment at 100.0 percent identity.

Fix is the bug report's Option A: accept on `qcovs` **or** on an absolute aligned length.

**Files:**
- Modify: `R/blast_genbank.R` (constants block and `ref_candidate_ok()`)
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: `ref_candidate_ok()` from Task 2.
- Produces: `REF_CANDIDATE_MIN_ALN_BP` constant; `ref_candidate_ok()` gains an
  `aln_len` dependency, so any caller must pass a frame carrying that column.
  `blast_hits` has it.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("a long assembly matching a short reference end to end is kept", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  # SRR14433904: 15120 bp unit, 6997 bp reference matched end to end, qcovs 46
  f <- write_hits(hit("MZ286333.1", "Gymnothorax fimbriatus small subunit ribosomal RNA gene",
    100.0, 46, 12922, 6997, 4998, 11994))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  cand <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_candidates")
  expect_equal(nrow(cand), 1)
  expect_equal(cand$accession[1], "MZ286333.1")
})

test_that("a short incidental alignment is still rejected", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(hit("AA9", "Danio rerio 18S ribosomal RNA gene", 99.0, 5, 300, 200, 1, 200))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 0L)
})

test_that("the identity floor still applies on the aligned-length branch", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(hit("AA8", "Danio rerio 28S ribosomal RNA gene", 70.0, 20, 2000, 3000, 1, 3000))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 0L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: the first new test FAILS with `nrow(cand)` equal to 0, the other two PASS.

- [ ] **Step 3: Add the constant and widen the predicate**

In `R/blast_genbank.R`, after `REF_CANDIDATE_MIN_QCOVS <- 50`, add:

```r
#' Absolute aligned length that admits a subject regardless of query coverage.
#'
#' qcovs against a full rDNA repeat measures recovered intergenic spacer, not
#' match quality, so it falls as an assembly improves. 1500 bp is most of an 18S
#' and well under a 28S. NOT yet calibrated outside Actinopterygii: confirm on a
#' clade with shorter elements before treating it as settled.
#' @noRd
REF_CANDIDATE_MIN_ALN_BP <- 1500
```

Change `ref_candidate_ok()` to:

```r
ref_candidate_ok <- function(x) {
  aln <- if ("aln_len" %in% names(x)) x$aln_len else rep(NA_integer_, nrow(x))
  !is.na(x$hit_class) & x$hit_class == "rdna" &
    !is.na(x$pident) & x$pident >= REF_CANDIDATE_MIN_PIDENT &
    ((!is.na(x$qcovs) & x$qcovs >= REF_CANDIDATE_MIN_QCOVS) |
      (!is.na(aln) & aln >= REF_CANDIDATE_MIN_ALN_BP))
}
```

The `aln_len` guard keeps the predicate working on a pre-v9 hits frame that lacks the
column.

- [ ] **Step 4: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: PASS, at least 35 tests.

- [ ] **Step 5: Commit**

```bash
git add R/blast_genbank.R tests/testthat/test-blast-parse.R
git commit -m "fix: accept reference candidates on aligned length as well as qcovs"
```

---

### Task 4: Re-derive an existing project without touching NCBI

Task 3 only helps projects searched after it ships unless the derivation can be replayed.
`blast_hits` already holds everything needed, so a project-wide re-derive recovers the
affected units at zero network cost. This is also the mechanism that makes every future
threshold change cheap.

**Files:**
- Modify: `R/blast_genbank.R` (append function)
- Create: `inst/nextflow/bin/rp_rederive.R`
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: `derive_blast_products()` from Task 2.
- Produces: `rederive_blast_products(con, ids = NULL)` returning an integer count of
  units re-derived.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("rederive_blast_products recovers units across a whole project", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',2,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  f1 <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  f2 <- write_hits(hit("AA2", "Danio rerio 28S ribosomal RNA gene", 98.0, 70, 1400, 1500, 1, 1500))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f1)
  RiboPilot:::parse_blast_candidates(con, "S1", 2, 1, f2)

  DBI::dbExecute(con, "DELETE FROM blast_ref_candidates")
  DBI::dbExecute(con, "DELETE FROM blast_screen")
  n <- RiboPilot:::rederive_blast_products(con)

  expect_equal(n, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 2L)
})

test_that("rederive_blast_products can be limited to given samples", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S2',1,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2',7955,NULL,NULL)")
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  RiboPilot:::parse_blast_candidates(con, "S2", 1, 1, f)
  expect_equal(RiboPilot:::rederive_blast_products(con, ids = "S2"), 1L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: FAIL, `could not find function "rederive_blast_products"`.

- [ ] **Step 3: Implement**

Append to `R/blast_genbank.R`:

```r
#' Re-derive stored BLAST products for every unit that has hits
#'
#' Replays [derive_blast_products()] over `blast_hits`, so a changed threshold or
#' a fixed filter reaches existing projects without re-running a remote search.
#'
#' @param con Open project DB connection.
#' @param ids Optional character vector of sample IDs; default is every sample.
#' @return (invisibly) the number of units re-derived.
#' @export
rederive_blast_products <- function(con, ids = NULL) {
  sql <- "SELECT DISTINCT ID, path, scaffold FROM blast_hits"
  units <- if (is.null(ids)) {
    DBI::dbGetQuery(con, paste(sql, "ORDER BY ID, path, scaffold"))
  } else {
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbGetQuery(con, sprintf("%s WHERE ID IN (%s) ORDER BY ID, path, scaffold", sql, ph),
      params = as.list(ids))
  }
  if (nrow(units) == 0) return(invisible(0L))
  for (i in seq_len(nrow(units))) {
    derive_blast_products(con, units$ID[i], units$path[i], units$scaffold[i])
  }
  invisible(nrow(units))
}
```

- [ ] **Step 4: Add the pipeline entry point**

Create `inst/nextflow/bin/rp_rederive.R`, matching the shape of
`inst/nextflow/bin/rp_blast.R`:

```r
#!/usr/bin/env Rscript
# Re-derive reference candidates + screen from stored blast_hits.
# Args: <db> [id ...]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 1, "rp_rederive.R <db> [id ...]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  ids <- if (length(a) > 1) a[-1] else NULL
  n <- rederive_blast_products(con, ids = ids)
  cat("units re-derived:", n, "\n")
})
```

Make it executable: `chmod +x inst/nextflow/bin/rp_rederive.R`

- [ ] **Step 5: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: PASS, at least 37 tests.

`tests/testthat/test-nf-bin-contract.R` asserts things about the scripts in
`inst/nextflow/bin/`. Run it and fix any contract it enforces (for example a required
shebang or argument-count check):
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-nf-bin-contract.R", reporter="summary")'`
Expected: PASS.

- [ ] **Step 6: Verify against the real affected project**

This is the acceptance check for Tasks 3 and 4 together. Work on a copy, never the
original:

```bash
cp /home/dmacguig/Documents/DJMCompBio/RiboPilot_test/full_project_docker/.sqlite /tmp/rederive_test.sqlite
Rscript -e '
  suppressMessages(pkgload::load_all(".", quiet=TRUE))
  con <- DBI::dbConnect(RSQLite::SQLite(), "/tmp/rederive_test.sqlite")
  before <- DBI::dbGetQuery(con, "SELECT COUNT(DISTINCT ID||path||scaffold) n FROM blast_ref_candidates")$n
  n <- RiboPilot:::rederive_blast_products(con)
  after <- DBI::dbGetQuery(con, "SELECT COUNT(DISTINCT ID||path||scaffold) n FROM blast_ref_candidates")$n
  cat("units re-derived:", n, "  units with candidates before:", before, " after:", after, "\n")
  cat("SRR14433904 candidates:",
      DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates WHERE ID=\"SRR14433904\"")$n, "\n")
  DBI::dbDisconnect(con)'
rm -f /tmp/rederive_test.sqlite
```

Expected: units with candidates rises from 16 to 26, and `SRR14433904` goes from 0
candidates to a non-zero count across its 10 paths.

If the count does not rise to 26, stop and report rather than adjusting
`REF_CANDIDATE_MIN_ALN_BP` to force it.

- [ ] **Step 7: Commit**

```bash
git add R/blast_genbank.R inst/nextflow/bin/rp_rederive.R tests/testthat/test-blast-parse.R
git commit -m "feat: re-derive BLAST products from stored hits without re-searching"
```

---

### Task 5: Capture slen and qseqid in the search output

Every `-outfmt` change costs one full remote re-run of every project, at roughly two
hours per unit. So add every column that might be wanted now, once. `slen` enables the
principled subject-coverage filter later with no second re-run; `qseqid` is insurance if
per-sample query batching is ever revisited.

`parse_blast_hits_file()` reads trailing columns through a guarded accessor, so older
hits files return NA rather than breaking.

**Files:**
- Modify: `inst/nextflow/modules/blast.nf:63`
- Modify: `R/blast_genbank.R` (`parse_blast_hits_file()`, and the `hits` frame)
- Modify: `R/db_schema.R` (`BLAST_HITS_DDL`)
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: `parse_blast_hits_file()` from the existing code.
- Produces: `blast_hits.slen` (INTEGER, nullable). `qseqid` is captured in the search
  output but deliberately **not** stored: nothing consumes it yet, and adding an unread
  column is what Task 1 just removed.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("slen is parsed when present and NA when absent", {
  new_fmt <- paste("BB1", "Danio rerio 18S ribosomal RNA gene", 99.0, 80, "1e-90",
    7955, 1500, 1600, 1, 1600, 1, 1600, 6997, "q1", sep = "\t")
  d <- RiboPilot:::parse_blast_hits_file(write_hits(new_fmt))
  expect_equal(d$slen[1], 6997L)

  old_fmt <- paste("BB2", "Danio rerio 18S ribosomal RNA gene", 99.0, 80, "1e-90",
    7955, 1500, 1600, 1, 1600, 1, 1600, sep = "\t")
  d2 <- RiboPilot:::parse_blast_hits_file(write_hits(old_fmt))
  expect_true(is.na(d2$slen[1]))
})

test_that("slen reaches blast_hits", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  f <- write_hits(paste("BB1", "Danio rerio 18S ribosomal RNA gene", 99.0, 80, "1e-90",
    7955, 1500, 1600, 1, 1600, 1, 1600, 6997, "q1", sep = "\t"))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, f)
  expect_equal(DBI::dbGetQuery(con, "SELECT slen FROM blast_hits")$slen[1], 6997L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: FAIL, `d$slen` is NULL.

- [ ] **Step 3: Add the column to the schema**

In `R/db_schema.R`, in `BLAST_HITS_DDL`, add after the `send INTEGER,` line:

```
  slen INTEGER,
```

- [ ] **Step 4: Parse it**

In `R/blast_genbank.R`, inside `parse_blast_hits_file()`, add `slen` to the row frame
after `send`:

```r
      sstart = int(g(11)), send = int(g(12)), slen = int(g(13)),
```

Then in `parse_blast_candidates()`, add `slen = h$slen,` to the `hits` data.frame,
immediately after `send = h$send,`.

- [ ] **Step 5: Extend the search output format**

In `inst/nextflow/modules/blast.nf`, change line 63 from:

```
        OUTFMT="6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send"
```

to:

```
        # slen enables a subject-coverage filter without a second remote re-run;
        # qseqid is captured for a future batched query. Appended last so the
        # positional parser in parse_blast_hits_file() needs no change.
        OUTFMT="6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"
```

- [ ] **Step 6: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: PASS, at least 39 tests.

Then the Nextflow contract tests, which read the module files:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); for (f in c("test-nextflow.R","test-nf-bin-contract.R")) testthat::test_file(file.path("tests/testthat", f), reporter="summary")'`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add R/blast_genbank.R R/db_schema.R inst/nextflow/modules/blast.nf tests/testthat/test-blast-parse.R
git commit -m "feat: capture slen and qseqid in the remote search output"
```

---

### Task 6: Stop swallowing efetch failures and stop leaking the API key

Two independent defects in the same file, both fixed by editing adjacent lines, so they
share a task.

`inst/nextflow/modules/blast_ref.nf:31-32` ends both efetch calls with `|| true`, and the
error-body grep misses NCBI's letter-spaced form. Measured consequence: `e2e_ponytail`
requested 187 accessions, stored 50, and reported success.

`export NCBI_API_KEY='${params.ncbi_api_key}'` is Groovy interpolation, written verbatim
into `.command.sh` and `.command.run` in every task work directory. Work directories on
shared HPC filesystems are typically group readable. The key is also dead code in
`blast.nf`: it is an E-utilities credential with no effect on the BLAST URL API.

**Files:**
- Modify: `inst/nextflow/modules/blast_ref.nf:20-45`
- Modify: `inst/nextflow/modules/blast.nf:41`
- Modify: `inst/nextflow/nextflow.config` (env scope)
- Test: `tests/testthat/test-nextflow.R` (extend)

**Interfaces:**
- Consumes: nothing.
- Produces: no R-level interface. The fetch loop writes a per-accession outcome line to
  stdout, and the process fails loudly if it stored nothing for a sample that had
  candidates.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-nextflow.R`:

```r
test_that("the API key is not interpolated into task scripts", {
  for (f in c("modules/blast.nf", "modules/blast_ref.nf")) {
    txt <- readLines(system.file("nextflow", f, package = "RiboPilot"), warn = FALSE)
    expect_false(any(grepl("export NCBI_API_KEY='\\$\\{params", txt, perl = TRUE)),
      info = paste(f, "interpolates the API key into the task script"))
  }
})

test_that("efetch failures are not swallowed", {
  txt <- readLines(system.file("nextflow", "modules/blast_ref.nf", package = "RiboPilot"),
    warn = FALSE)
  efetch <- grep("efetch -db nuccore", txt, value = TRUE)
  expect_true(length(efetch) > 0)
  expect_false(any(grepl("\\|\\| true", efetch)),
    info = "efetch still ends with || true")
})
```

`system.file(..., package = "RiboPilot")` resolves against the loaded source tree under
`pkgload::load_all()`. If it returns `""`, use `file.path("inst/nextflow", f)` instead.

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-nextflow.R", reporter="summary")'`
Expected: FAIL on both tests.

- [ ] **Step 3: Move the API key out of the task script**

In `inst/nextflow/nextflow.config`, add an `env` scope that reads from the launching
shell rather than from an interpolated parameter:

```groovy
env {
    NCBI_API_KEY = System.getenv('NCBI_API_KEY') ?: ''
}
```

Then delete this line from `inst/nextflow/modules/blast_ref.nf` (line 23):

```
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'
```

And delete the same line from `inst/nextflow/modules/blast.nf` (line 41). In `blast.nf`
it is dead code regardless: `NCBI_API_KEY` is an E-utilities credential and has no effect
on `blastn -remote`.

- [ ] **Step 4: Make fetch failures visible**

In `inst/nextflow/modules/blast_ref.nf`, replace the fetch loop body (lines 28-44) with:

```bash
    STORED=0
    WANTED=0
    for acc in \$ACCS; do
        WANTED=\$(( WANTED + 1 ))
        have=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT 1 FROM blast_ref_sequences WHERE accession = '\$acc';")
        if [ -n "\$have" ]; then STORED=\$(( STORED + 1 )); continue; fi

        if ! efetch -db nuccore -id "\$acc" -format fasta > "\$acc.fa"; then
            echo "efetch failed for \$acc (non-zero exit)"; continue
        fi
        if ! efetch -db nuccore -id "\$acc" -format ft > "\$acc.ft"; then
            echo "efetch ft failed for \$acc (non-zero exit)"; continue
        fi
        sleep \$PACE

        if [ ! -s "\$acc.fa" ]; then echo "efetch returned an empty body for \$acc"; continue; fi
        # NCBI answers HTTP 200 with an XML error body, and also with a
        # letter-spaced plain-text error that the XML patterns do not match.
        if grep -q -e '<ERROR>' -e '<eFetchResult>' -e 'F a i l e d' -e '^Error' "\$acc.fa" "\$acc.ft" 2>/dev/null; then
            echo "efetch returned an error body for \$acc"; continue
        fi
        if [ "\$(head -c 1 "\$acc.fa")" != ">" ]; then echo "efetch returned a non-fasta body for \$acc"; continue; fi
        if [ -s "\$acc.ft" ] && ! head -n 1 "\$acc.ft" | grep -q '^>Feature'; then
            echo "efetch returned a non-feature-table body for \$acc"; continue
        fi
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_fetch.R '${params.db}' "\$acc" "\$acc.fa" "\$acc.ft"
        STORED=\$(( STORED + 1 ))
    done

    echo "reference fetch for ${id}: stored \${STORED} of \${WANTED} candidate accession(s)"
    if [ "\${WANTED}" -gt 0 ] && [ "\${STORED}" -eq 0 ]; then
        echo "stored nothing for ${id} despite \${WANTED} candidates" >&2
        exit 1
    fi
```

The final guard turns the silent case into a retry. A sample that genuinely has no
fetchable references still exits 0 once `WANTED` is 0, which the existing early return
above already handles.

- [ ] **Step 5: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-nextflow.R", reporter="summary")'`
Expected: PASS.

Confirm the config still parses. If Nextflow is on PATH:
`nextflow config -profile standard inst/nextflow/nextflow.config >/dev/null && echo "config OK"`
If Nextflow is not installed, skip this and note it in the commit body.

- [ ] **Step 6: Commit**

```bash
git add inst/nextflow/modules/blast_ref.nf inst/nextflow/modules/blast.nf inst/nextflow/nextflow.config tests/testthat/test-nextflow.R
git commit -m "fix: surface efetch failures and keep the API key out of work dirs"
```

---

### Task 7: Record that a unit was searched, keyed on its query sequence

The resume guard at `inst/nextflow/modules/blast.nf:51` counts rows in `blast_hits`, so a
unit that returned nothing is re-searched on every retry and every re-run. Measured: 20 of
26 units in `full_project_docker` returned zero hits and cost two remote searches each,
40 of the 46 searches in that project.

Keying the marker on a hash of the query sequence fixes a second bug in the same change:
nothing currently clears `blast_hits` when an assembly is re-assembled, trimmed, or
edited, so `blast_screen` and `blast_ref_candidates` can keep describing a contig that no
longer exists.

**Files:**
- Modify: `R/db_schema.R` (add DDL constant)
- Modify: `R/init_db.R` (create the table)
- Modify: `inst/nextflow/modules/blast.nf:45-53` (guard) and after the search (write)
- Create: `inst/nextflow/bin/rp_searched.R`
- Test: `tests/testthat/test-blast-parse.R` (extend)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: table `blast_searched(ID, path, scaffold, query_md5, n_hits, time_stamp)`
  with `PRIMARY KEY (ID, path, scaffold)`; R function
  `mark_unit_searched(con, id, path, scaffold, query_md5, n_hits)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-blast-parse.R`:

```r
test_that("a searched unit is recorded even when it found nothing", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc123", 0L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM blast_searched")
  expect_equal(nrow(r), 1L)
  expect_equal(r$query_md5[1], "abc123")
  expect_equal(r$n_hits[1], 0L)
})

test_that("re-marking the same unit replaces rather than duplicates", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc123", 0L)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "def456", 7L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM blast_searched")
  expect_equal(nrow(r), 1L)
  expect_equal(r$query_md5[1], "def456")
  expect_equal(r$n_hits[1], 7L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: FAIL, object `BLAST_SEARCHED_DDL` not found.

- [ ] **Step 3: Add the DDL and the writer**

In `R/db_schema.R`, after `BLAST_HITS_DDL`, add:

```r
#' DDL for blast_searched: which units have been sent to the remote search.
#'
#' Keyed on a hash of the query sequence, so a re-assembled or trimmed unit is
#' searched again rather than inheriting the previous contig's hits. A unit that
#' found nothing still gets a row, which is what stops it being re-searched on
#' every retry and every re-run.
#' @noRd
BLAST_SEARCHED_DDL <- "CREATE TABLE IF NOT EXISTS blast_searched (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  query_md5 TEXT NOT NULL,
  n_hits INTEGER,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"
```

In `R/init_db.R`, next to the `BLAST_HITS_DDL` execution (line 349), add:

```r
  DBI::dbExecute(con, BLAST_SEARCHED_DDL)
```

In `R/blast_genbank.R`, append:

```r
#' Record that a unit's query was sent to the remote search
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param query_md5 md5 of the exact query sequence submitted.
#' @param n_hits Number of hits the search returned.
#' @return (invisibly) TRUE.
#' @export
mark_unit_searched <- function(con, id, path = 1L, scaffold = 1L, query_md5,
                               n_hits = 0L) {
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO blast_searched (ID, path, scaffold, query_md5, n_hits, time_stamp)
     VALUES (?, ?, ?, ?, ?, ?)",
    params = list(id, as.integer(path), as.integer(scaffold), query_md5,
      as.integer(n_hits), as.integer(Sys.time())))
  invisible(TRUE)
}
```

- [ ] **Step 4: Add the pipeline entry point**

Create `inst/nextflow/bin/rp_searched.R`:

```r
#!/usr/bin/env Rscript
# Record that a unit was searched. Args: <db> <id> <path> <scaffold> <md5> <n_hits>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_searched.R <db> <id> <path> <scaffold> <md5> <n_hits>")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  mark_unit_searched(con, a[2], as.integer(a[3]), as.integer(a[4]), a[5], as.integer(a[6]))
})
cat("searched:", a[2], a[3], a[4], "\n")
```

`chmod +x inst/nextflow/bin/rp_searched.R`

- [ ] **Step 5: Use the marker as the resume guard**

In `inst/nextflow/modules/blast.nf`, replace the `DONE` guard (lines 51-52) with a
hash comparison, and write the marker after the search. Inside the `for U in \$UNITS`
loop, after `query.fa` is written and the `[ -s query.fa ] || continue` check:

```bash
        QMD5=\$(grep -v '^>' query.fa | tr -d '[:space:]' | md5sum | cut -d' ' -f1)
        PREV=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT query_md5 FROM blast_searched WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" 2>/dev/null || echo "")
        if [ "\${PREV}" = "\${QMD5}" ]; then echo "already searched ${id} \$P.\$S"; continue; fi
```

Move the existing `query.fa` write above this block, and delete the old `DONE` lookup
entirely. After `rp_blast.R` runs at the end of the loop body, add:

```bash
        NHITS=\$(grep -c . hits.txt 2>/dev/null || echo 0)
        grep -q NO_SIGNIFICANT_HITS hits.txt && NHITS=0
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_searched.R '${params.db}' '${id}' "\$P" "\$S" "\$QMD5" "\$NHITS"
```

- [ ] **Step 6: Run tests to verify they pass**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-blast-parse.R", reporter="summary")'`
Expected: PASS, at least 41 tests.

Then the full suite, because this task touches `init_db.R`:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_dir("tests/testthat", reporter="summary")'`
Expected: no new failures against the pre-plan baseline. Record the baseline first if you
have not: `git stash && Rscript -e '...test_dir...' > /tmp/baseline.txt; git stash pop`.

- [ ] **Step 7: Commit**

```bash
git add R/db_schema.R R/init_db.R R/blast_genbank.R inst/nextflow/bin/rp_searched.R inst/nextflow/modules/blast.nf tests/testthat/test-blast-parse.R
git commit -m "fix: record searched units by query hash so empty results are not re-searched"
```

---

## After this plan

Stage 1 is complete when all seven tasks are committed and the full suite shows no new
failures. Stages 2 through 4 get their own plans:

- **Stage 2** (reference yield): sequence-only batched fetch, one R invocation per
  sample, local barrnap and ITSx annotation of references, and the three blockers the
  measurement found (multi-FASTA collapse in `parse_itsx()`, the ITSx 100 kb batch abort,
  and the barrnap/ITSx coordinate-frame mismatch on minus-strand records).
- **Stage 3** (shared reference store): the cross-project store under
  `ribopilot_ref_dir()`, including the read-only SQLite open path.
- **Stage 4** (local database): `params.blastdb` plumbing, the `blast_local` and
  `blast_remote_escalate` split, the five escalation triggers, and the NCBI decoy panel.

Stage 4 must not ship before its acceptance gate is measured on Asteroidea, Scyphozoa, or
Octocorallia: local query time and peak RSS, escalation rate, and local-versus-remote
concordance. The fish projects cannot answer any of the three.
