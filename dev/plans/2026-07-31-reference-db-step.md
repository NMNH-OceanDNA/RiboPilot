# Reference Database Step Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an in-app step between taxonomy assignment and the assembly workflow that suggests, builds, and wires per-clade GetOrganelle reference databases.

**Architecture:** Sample lineages already frozen by the taxonomy step are walked up rank by rank until a clade has enough GenBank rRNA records; samples landing on the same clade form a group. Groups live in a new `reference` project table and are built by detached `processx` jobs wrapping the existing `build_reference_db()`. A successful build creates an `assemble_opts` preset named after the clade and moves the group's samples onto it, so nothing in the Nextflow layer changes.

**Tech Stack:** R, Shiny, `reactable`, `shinyWidgets`, `RSQLite`/`DBI`, `processx`, `taxonkit` (via the project's pinned NCBI taxdump), Entrez Direct (via the shipped `download_rrna.sh`), `testthat` edition 3.

Spec: `dev/specs/2026-07-31-reference-db-step-design.md`.

## Global Constraints

- ASCII only in R code. No non-ASCII characters anywhere in R sources, including comments and roxygen.
- No em dashes anywhere. Use commas, semicolons, colons.
- Minimal comments. Comment the "why" where it is not obvious; do not narrate bugfixes.
- Narrow scope. Do not refactor code the task does not need to touch.
- Existing style: two-space indent, `%||%` for NULL-coalesce, `%|NA|%` for NA-coalesce (defined in `R/app_opts.R:20`), roxygen `@noRd` for internal helpers and full roxygen with `@export` for user-facing functions.
- Schema changes are additive and idempotent, guarded by `has_tbl` / `has_col` / `add_col` inside `migrate_db()`, exactly as every earlier migration does.
- Never contact NCBI from a test. Every network call is behind an injectable function argument.
- The taxonomy step already writes `samples.lineage` (7 rank-aligned names, `;`-joined) and `samples.lineage_taxids` (7 rank-aligned taxids, `;`-joined, empty at filled ranks). Rank order is `RRNA_LINEAGE_RANKS` = kingdom, phylum, class, order, family, genus, species (`R/taxonomy_ranks.R:20`). Parse them with the existing `parse_rank_names()` and `parse_rank_taxids()`. Do not add a new taxonkit call for this.
- Default record window: `min_records = 1000`, `max_records = 25000`.
- Reference-mode state codes reuse the app's existing state vocabulary: 0 hold, 1 ready, 2 done, 3 failed, 4 building, plus 5 skipped which is new.
- Commit after every task with a short conventional-commit subject. Never push.

---

### Task 1: Schema v19 (reference tables, `samples.ref_id`, `ref_dir` project argument)

**Files:**
- Modify: `R/migrate_db.R` (bump `RIBOPILOT_SCHEMA_VERSION`, add shared DDL constants, add the v19 block)
- Modify: `R/init_db.R` (create the same two tables in a fresh project, add `ref_id` to `samples`)
- Modify: `R/init_project.R` (new `ref_dir` argument, passed to `new_db()`)
- Test: `tests/testthat/test-reference-schema.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `REFERENCE_DDL`, `REFERENCE_CONFIG_DDL` character constants in `R/migrate_db.R`.
  - Table `reference(ref_id TEXT PK, clade_taxid INTEGER, clade_name TEXT, clade_rank TEXT, query TEXT, n_samples INTEGER, n_records INTEGER, reference_switch INTEGER, source TEXT, seeds_db TEXT, labels_db TEXT, log_file TEXT, pid INTEGER, fail_reason TEXT, time_stamp INTEGER)`.
  - Table `reference_config(rowid_key INTEGER PK, ref_dir TEXT, min_records INTEGER, max_records INTEGER)` holding exactly one row with `rowid_key = 1`.
  - Column `samples.ref_id TEXT`.
  - `new_db(..., ref_dir = NULL)` and `new_project(..., ref_dir = NULL)`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-schema.R`:

```r
# Schema v19: the in-app reference-database step.

mk_v18_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT PRIMARY KEY, cpus INT,
    memory INT, getOrganelle TEXT, seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INT,
    max_scaffolds INT, min_assembly_length INT, assembler TEXT, bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa;Chordata', '33208;7711', 1, 2)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1', 1, 0, 'default', 'default')")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default',4,16,'','s.fa','l.fa',NA,10,10,1000,'getorganelle','--local','--isolate')")
  DBI::dbExecute(con, "PRAGMA user_version = 18")
  con
}

test_that("migrating a v18 project adds the reference tables and ref_id", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("reference" %in% DBI::dbListTables(con))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("reference" %in% DBI::dbListTables(con))
  expect_true("reference_config" %in% DBI::dbListTables(con))
  expect_true("ref_id" %in% DBI::dbListFields(con, "samples"))
  for (f in c("ref_id", "clade_taxid", "clade_name", "clade_rank", "query", "n_samples",
              "n_records", "reference_switch", "source", "seeds_db", "labels_db",
              "log_file", "pid", "fail_reason", "time_stamp")) {
    expect_true(f %in% DBI::dbListFields(con, "reference"), info = f)
  }
})

test_that("the migration seeds exactly one reference_config row with the default window", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  cfg <- DBI::dbGetQuery(con, "SELECT * FROM reference_config")
  expect_equal(nrow(cfg), 1L)
  expect_true(is.na(cfg$ref_dir[1]))
  expect_equal(cfg$min_records[1], 1000)
  expect_equal(cfg$max_records[1], 25000)
})

test_that("the migration is idempotent and leaves existing rows alone", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  DBI::dbExecute(con, "UPDATE reference_config SET min_records = 250")
  migrate_db(con)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM reference_config")), 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT min_records FROM reference_config")$min_records[1], 250)
  expect_equal(DBI::dbGetQuery(con, "SELECT taxid FROM samples")$taxid[1], 7955)
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-schema.R")'`
Expected: FAIL, `reference` is not in `DBI::dbListTables(con)`.

- [ ] **Step 3: Add the DDL constants and the v19 migration block**

In `R/migrate_db.R`, bump the version:

```r
RIBOPILOT_SCHEMA_VERSION <- 19L
```

Add these constants next to the other DDL constants near the top of the file:

```r
#' DDL for the reference table (one row per clade reference-database group).
#' Shared by new_db() and the v19 migration so both agree on the schema.
#' @noRd
REFERENCE_DDL <- "CREATE TABLE IF NOT EXISTS reference (
  ref_id           TEXT NOT NULL,
  clade_taxid      INTEGER,
  clade_name       TEXT,
  clade_rank       TEXT,
  query            TEXT,
  n_samples        INTEGER,
  n_records        INTEGER,
  reference_switch INTEGER,
  source           TEXT,
  seeds_db         TEXT,
  labels_db        TEXT,
  log_file         TEXT,
  pid              INTEGER,
  fail_reason      TEXT,
  time_stamp       INTEGER,
  PRIMARY KEY (ref_id)
);"

#' DDL for reference_config: one row of project-level reference settings.
#'
#' `ref_dir` is the folder built databases are written to and reused from; NULL
#' means fall back to RIBOPILOT_REF_DIR and then ribopilot_ref_dir(). The record
#' window bounds the clade suggestion walk.
#' @noRd
REFERENCE_CONFIG_DDL <- "CREATE TABLE IF NOT EXISTS reference_config (
  rowid_key   INTEGER NOT NULL,
  ref_dir     TEXT,
  min_records INTEGER,
  max_records INTEGER,
  PRIMARY KEY (rowid_key)
);"
```

Add the migration block immediately before the `user_version` comment at the end of `migrate_db()`:

```r
  # --- schema v19: in-app reference-database step ---
  # Clade groups get their own table rather than columns on samples: a group is a
  # unit of work (built, failed, skipped) that outlives any one sample's membership.
  DBI::dbExecute(con, REFERENCE_DDL)
  DBI::dbExecute(con, REFERENCE_CONFIG_DDL)
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM reference_config")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
       VALUES (1, NULL, 1000, 25000)")
  }
  add_col("samples", "ref_id", "TEXT")
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-schema.R")'`
Expected: PASS, 3 tests.

- [ ] **Step 5: Create the same tables in a fresh project**

In `R/init_db.R`, add `ref_id = NA_character_` to the `mapping |> dplyr::mutate(...)` block that defines the `samples` columns (after `lineage_taxids = NA_character_`), and add this block after the `curate_opts` table is created, before the function returns:

```r
  # reference (clade DB groups) + reference_config ----
  DBI::dbExecute(con, REFERENCE_DDL)
  DBI::dbExecute(con, REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con,
    "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
     VALUES (1, ?, 1000, 25000)",
    params = list(ref_dir %||% NA_character_))
```

Add `ref_dir = NULL` to the `new_db()` signature and document it:

```r
#' @param ref_dir Folder for built reference databases (default: the
#'   `RIBOPILOT_REF_DIR` environment variable, else [ribopilot_ref_dir()]).
```

- [ ] **Step 6: Pass `ref_dir` through `new_project()`**

In `R/init_project.R`, add `ref_dir = NULL` to the `new_project()` signature next to `custom_labels_db`, document it with the same roxygen text, and pass `ref_dir = ref_dir` in the `new_db()` call alongside `seeds_db` / `labels_db`.

- [ ] **Step 7: Add a fresh-project test**

Append to `tests/testthat/test-reference-schema.R`:

```r
test_that("a fresh project database has the reference tables and stores ref_dir", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  mapping <- data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = file.path(dir, "S1_R1.fq.gz"), R2 = file.path(dir, "S1_R2.fq.gz"))
  file.create(mapping$R1, mapping$R2)
  refs <- withr::local_tempdir()
  new_db(db_path = file.path(dir, ".sqlite"), mapping = mapping, ref_dir = refs)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_true("reference" %in% DBI::dbListTables(con))
  expect_true("ref_id" %in% DBI::dbListFields(con, "samples"))
  expect_equal(DBI::dbGetQuery(con, "SELECT ref_dir FROM reference_config")$ref_dir[1], refs)
})
```

If `new_db()` rejects this minimal mapping (it validates read files exist and the seed/label databases are reachable), pass the shipped test seeds via `seeds_db` / `labels_db` pointing at `ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta` and `.../labels/actinopterygii_rrna_labels.fasta` resolved with `testthat::test_path("..", "..", ...)`, and add `skip_if_not(file.exists(...))`.

- [ ] **Step 8: Run the full suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: all tests pass, including the pre-existing `test-migrate.R`.

- [ ] **Step 9: Commit**

```bash
git add R/migrate_db.R R/init_db.R R/init_project.R tests/testthat/test-reference-schema.R
git commit -m "feat(reference): schema v19 for the reference database step"
```

---

### Task 2: Reference folder resolution and manifest

**Files:**
- Create: `R/reference_dir.R`
- Test: `tests/testthat/test-reference-dir.R`

**Interfaces:**
- Consumes: `reference_config` table from Task 1; `ribopilot_ref_dir()` from `R/build_reference_db.R:13`.
- Produces:
  - `ref_config(con)` returns `list(ref_dir = <chr>, min_records = <int>, max_records = <int>)`, `ref_dir` already resolved through the fallback chain and `path.expand`ed.
  - `set_ref_config(con, ref_dir = NULL, min_records = NULL, max_records = NULL)` returns the new config invisibly; only non-NULL arguments are written.
  - `ref_manifest_path(dir)` returns `file.path(dir, "manifest.tsv")`.
  - `read_ref_manifest(dir)` returns `data.frame(clade, taxid, rank, query, n_records, built, seeds, labels)`, zero rows when absent.
  - `append_ref_manifest(dir, clade, taxid, rank, query, n_records, seeds, labels)` returns the manifest path invisibly.
  - `validate_ref_dir(dir)` returns `list(ok = <lgl>, msg = <chr>)`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-dir.R`:

```r
# Reference-folder resolution order and the cross-project manifest.

mk_cfg_db <- function(ref_dir = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con,
    "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
     VALUES (1, ?, 1000, 25000)", params = list(ref_dir))
  con
}

test_that("the project setting wins over the environment variable", {
  d <- withr::local_tempdir()
  con <- mk_cfg_db(d); on.exit(DBI::dbDisconnect(con))
  withr::local_envvar(RIBOPILOT_REF_DIR = "/somewhere/else")
  expect_equal(RiboPilot:::ref_config(con)$ref_dir, path.expand(d))
})

test_that("an unset project setting falls back to the environment variable", {
  d <- withr::local_tempdir()
  con <- mk_cfg_db(NA_character_); on.exit(DBI::dbDisconnect(con))
  withr::local_envvar(RIBOPILOT_REF_DIR = d)
  expect_equal(RiboPilot:::ref_config(con)$ref_dir, path.expand(d))
})

test_that("set_ref_config writes only the arguments it is given", {
  con <- mk_cfg_db(NA_character_); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_ref_config(con, min_records = 200)
  cfg <- DBI::dbGetQuery(con, "SELECT * FROM reference_config")
  expect_equal(cfg$min_records[1], 200)
  expect_equal(cfg$max_records[1], 25000)
  expect_true(is.na(cfg$ref_dir[1]))
})

test_that("validate_ref_dir rejects a path that is not a writable directory", {
  f <- withr::local_tempfile()
  writeLines("x", f)
  expect_false(RiboPilot:::validate_ref_dir(f)$ok)
  expect_true(RiboPilot:::validate_ref_dir(withr::local_tempdir())$ok)
})

test_that("the manifest round-trips a built database", {
  d <- withr::local_tempdir()
  expect_equal(nrow(RiboPilot:::read_ref_manifest(d)), 0L)
  RiboPilot:::append_ref_manifest(d, clade = "asteroidea", taxid = 7588L, rank = "class",
    query = "\"Asteroidea\"[Organism]", n_records = 4102L,
    seeds = file.path(d, "seeds", "asteroidea_rrna_seeds.fasta"),
    labels = file.path(d, "labels", "asteroidea_rrna_labels.fasta"))
  m <- RiboPilot:::read_ref_manifest(d)
  expect_equal(nrow(m), 1L)
  expect_equal(m$clade[1], "asteroidea")
  expect_equal(m$taxid[1], 7588L)
  expect_equal(m$n_records[1], 4102L)
})

test_that("appending twice keeps one row per clade, the newest", {
  d <- withr::local_tempdir()
  RiboPilot:::append_ref_manifest(d, "asteroidea", 7588L, "class", "q", 10L, "s1", "l1")
  RiboPilot:::append_ref_manifest(d, "asteroidea", 7588L, "class", "q", 20L, "s2", "l2")
  m <- RiboPilot:::read_ref_manifest(d)
  expect_equal(nrow(m), 1L)
  expect_equal(m$n_records[1], 20L)
  expect_equal(m$seeds[1], "s2")
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-dir.R")'`
Expected: FAIL, `ref_config` not found.

- [ ] **Step 3: Write the implementation**

Create `R/reference_dir.R`:

```r
# Where built reference databases live, and what is already there.
#
# The folder is deliberately outside any project so one clade database serves
# every project (and, on a shared drive, every user) that needs that clade. The
# manifest is the only index of what a folder holds, so a fresh project can
# offer an existing database instead of downloading it again.

#' Project-level reference settings, with the folder resolved
#'
#' Resolution order: the project setting, then `RIBOPILOT_REF_DIR`, then
#' [ribopilot_ref_dir()].
#' @param con Open project DB connection.
#' @return list(ref_dir, min_records, max_records).
#' @noRd
ref_config <- function(con) {
  d <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM reference_config WHERE rowid_key = 1"),
    error = function(e) data.frame())
  stored <- if (nrow(d) > 0) d$ref_dir[1] else NA_character_
  dir <- if (!is.na(stored) && nzchar(trimws(stored))) trimws(stored) else ribopilot_ref_dir()
  list(
    ref_dir = path.expand(dir),
    min_records = if (nrow(d) > 0) as.integer(d$min_records[1] %|NA|% 1000L) else 1000L,
    max_records = if (nrow(d) > 0) as.integer(d$max_records[1] %|NA|% 25000L) else 25000L
  )
}

#' Write project-level reference settings (only the arguments supplied)
#' @param con Open project DB connection.
#' @param ref_dir,min_records,max_records New values, or NULL to leave as is.
#' @return (invisibly) the resulting [ref_config()].
#' @noRd
set_ref_config <- function(con, ref_dir = NULL, min_records = NULL, max_records = NULL) {
  if (!is.null(ref_dir)) {
    DBI::dbExecute(con, "UPDATE reference_config SET ref_dir = ? WHERE rowid_key = 1",
      params = list(if (nzchar(trimws(ref_dir))) path.expand(trimws(ref_dir)) else NA_character_))
  }
  if (!is.null(min_records)) {
    DBI::dbExecute(con, "UPDATE reference_config SET min_records = ? WHERE rowid_key = 1",
      params = list(as.integer(min_records)))
  }
  if (!is.null(max_records)) {
    DBI::dbExecute(con, "UPDATE reference_config SET max_records = ? WHERE rowid_key = 1",
      params = list(as.integer(max_records)))
  }
  invisible(ref_config(con))
}

#' Is `dir` usable as a reference folder? (created if missing)
#' @return list(ok, msg).
#' @noRd
validate_ref_dir <- function(dir) {
  dir <- path.expand(trimws(dir %||% ""))
  if (!nzchar(dir)) return(list(ok = FALSE, msg = "Give a folder path."))
  if (file.exists(dir) && !dir.exists(dir)) {
    return(list(ok = FALSE, msg = paste0(dir, " exists but is a file, not a folder.")))
  }
  if (!dir.exists(dir)) {
    ok <- dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) return(list(ok = FALSE, msg = paste0("Could not create ", dir, ".")))
  }
  if (unname(file.access(dir, mode = 2)) != 0) {
    return(list(ok = FALSE, msg = paste0(dir, " is not writable.")))
  }
  list(ok = TRUE, msg = dir)
}

#' Path to a reference folder's manifest
#' @noRd
ref_manifest_path <- function(dir) file.path(path.expand(dir), "manifest.tsv")

#' Manifest columns, in file order
#' @noRd
REF_MANIFEST_COLS <- c("clade", "taxid", "rank", "query", "n_records", "built", "seeds", "labels")

#' Read a reference folder's manifest
#' @return data.frame with [REF_MANIFEST_COLS]; zero rows when there is no manifest.
#' @noRd
read_ref_manifest <- function(dir) {
  empty <- data.frame(clade = character(0), taxid = integer(0), rank = character(0),
    query = character(0), n_records = integer(0), built = character(0),
    seeds = character(0), labels = character(0), stringsAsFactors = FALSE)
  p <- ref_manifest_path(dir)
  if (!file.exists(p)) return(empty)
  d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character"),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !all(REF_MANIFEST_COLS %in% names(d))) return(empty)
  d <- d[, REF_MANIFEST_COLS, drop = FALSE]
  d$taxid <- suppressWarnings(as.integer(d$taxid))
  d$n_records <- suppressWarnings(as.integer(d$n_records))
  d
}

#' Record a built database in the folder's manifest (one row per clade, newest wins)
#' @return (invisibly) the manifest path.
#' @noRd
append_ref_manifest <- function(dir, clade, taxid, rank, query, n_records, seeds, labels) {
  dir.create(path.expand(dir), recursive = TRUE, showWarnings = FALSE)
  m <- read_ref_manifest(dir)
  m <- m[m$clade != clade, , drop = FALSE]
  m <- rbind(m, data.frame(
    clade = clade, taxid = as.integer(taxid), rank = rank %||% NA_character_,
    query = query, n_records = as.integer(n_records),
    built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seeds = seeds, labels = labels, stringsAsFactors = FALSE))
  p <- ref_manifest_path(dir)
  utils::write.table(m, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  invisible(p)
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-dir.R")'`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
git add R/reference_dir.R tests/testthat/test-reference-dir.R
git commit -m "feat(reference): reference folder resolution and manifest"
```

---

### Task 3: GenBank record counts (network call, cached, injectable)

**Files:**
- Create: `R/reference_counts.R`
- Test: `tests/testthat/test-reference-counts.R`

**Interfaces:**
- Consumes: `read_ref_manifest()` conventions from Task 2 (same folder, sibling file); `app_sys()` and `resolve_conda_bin()` used by `R/build_reference_db.R`.
- Produces:
  - `clade_query(name)` returns `"\"<name>\"[Organism]"`.
  - `parse_record_count(lines)` returns an integer or `NA_integer_`.
  - `genbank_rrna_count(query, api_key = NULL, refseq = FALSE, conda_env = NULL)` returns an integer or `NA_integer_`.
  - `cached_count_fn(dir, count_fn = genbank_rrna_count, ...)` returns `function(taxid, name)` giving an integer count, caching to `<dir>/counts.tsv`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-counts.R`:

```r
# GenBank record counts: parsing the builder script's output, and the on-disk
# cache that keeps the suggestion walk down to a handful of network calls.

test_that("parse_record_count reads the builder script's count line", {
  out <- c("GenBank query:", "  (\"Asteroidea\"[Organism]) AND (...)",
           "", "Found 4102 matching nuclear rRNA record(s).")
  expect_equal(RiboPilot:::parse_record_count(out), 4102L)
})

test_that("parse_record_count returns NA when the count line is absent", {
  expect_true(is.na(RiboPilot:::parse_record_count(c("boom", "esearch: command not found"))))
})

test_that("clade_query quotes the organism term", {
  expect_equal(RiboPilot:::clade_query("Asteroidea"), "\"Asteroidea\"[Organism]")
})

test_that("cached_count_fn calls the network once per taxid", {
  d <- withr::local_tempdir()
  calls <- 0L
  fake <- function(query, ...) { calls <<- calls + 1L; 4102L }
  f <- RiboPilot:::cached_count_fn(d, count_fn = fake)
  expect_equal(f(7588L, "Asteroidea"), 4102L)
  expect_equal(f(7588L, "Asteroidea"), 4102L)
  expect_equal(calls, 1L)
})

test_that("the count cache survives a new session", {
  d <- withr::local_tempdir()
  f1 <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) 99L)
  f1(7588L, "Asteroidea")
  f2 <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) stop("no network"))
  expect_equal(f2(7588L, "Asteroidea"), 99L)
})

test_that("a failed lookup is not cached", {
  d <- withr::local_tempdir()
  n <- 0L
  f <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) {
    n <<- n + 1L
    if (n == 1L) NA_integer_ else 7L
  })
  expect_true(is.na(f(1L, "Nowhere")))
  expect_equal(f(1L, "Nowhere"), 7L)
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-counts.R")'`
Expected: FAIL, `parse_record_count` not found.

- [ ] **Step 3: Write the implementation**

Create `R/reference_counts.R`:

```r
# How many nuclear rRNA records GenBank holds for a clade.
#
# This is the evidence the clade suggestion walks on: a clade nobody has
# sequenced makes a useless seed database, and one covering half of Metazoa
# makes an unspecific one. The count comes from the same query the builder will
# run, so the number shown is the number that will be downloaded.

#' GenBank organism term for a clade name
#' @noRd
clade_query <- function(name) sprintf("\"%s\"[Organism]", trimws(name))

#' Pull the record count out of the builder script's output
#' @param lines Character vector of script stdout.
#' @return integer, or NA when the script did not report a count.
#' @noRd
parse_record_count <- function(lines) {
  hit <- grep("^Found [0-9]+ matching", lines, value = TRUE)
  if (length(hit) == 0) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Found ([0-9]+) .*$", "\\1", hit[1])))
}

#' Count GenBank nuclear rRNA records for a query
#'
#' Runs the shipped builder script in `--count-only` mode, so the count matches
#' exactly what a build of the same query would fetch. Requires Entrez Direct on
#' PATH (pass `conda_env` when the app is not running inside that environment).
#'
#' @param query Advanced GenBank search term, e.g. `'"Asteroidea"[Organism]'`.
#' @param api_key Optional NCBI API key.
#' @param refseq Restrict to RefSeq records.
#' @param conda_env Optional conda env (name or path) providing Entrez Direct.
#' @return integer record count, or `NA_integer_` if the count could not be read.
#' @noRd
genbank_rrna_count <- function(query, api_key = NULL, refseq = FALSE, conda_env = NULL) {
  script <- app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
  if (!nzchar(script) || !file.exists(script)) return(NA_integer_)
  if (!is.null(conda_env) && nzchar(conda_env)) {
    bin <- tryCatch(resolve_conda_bin(conda_env), error = function(e) NULL)
    if (!is.null(bin)) {
      old <- Sys.getenv("PATH")
      Sys.setenv(PATH = paste(bin, old, sep = .Platform$path.sep))
      on.exit(Sys.setenv(PATH = old), add = TRUE)
    }
  }
  args <- c(shQuote(script), "--count-only")
  if (isTRUE(refseq)) args <- c(args, "--refseq")
  if (!is.null(api_key) && nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))
  args <- c(args, shQuote(query), "count_probe")
  out <- suppressWarnings(tryCatch(system2("bash", args, stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)))
  parse_record_count(out)
}

#' Path to a reference folder's record-count cache
#' @noRd
ref_counts_path <- function(dir) file.path(path.expand(dir), "counts.tsv")

#' Record-count lookup with an on-disk cache
#'
#' Counts change slowly and the walk revisits the same high ranks for every
#' project, so they are cached per reference folder. Failures are not cached: a
#' missing Entrez Direct today should not poison the folder forever.
#'
#' @param dir Reference folder holding the cache.
#' @param count_fn Function(query, api_key, refseq, conda_env) -> integer. Injected
#'   so tests never reach the network.
#' @param ... Passed to `count_fn` (api_key, refseq, conda_env).
#' @return function(taxid, name) -> integer count (NA when unavailable).
#' @noRd
cached_count_fn <- function(dir, count_fn = genbank_rrna_count, ...) {
  extra <- list(...)
  p <- ref_counts_path(dir)
  read_cache <- function() {
    if (!file.exists(p)) {
      return(data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE))
    }
    d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !all(c("taxid", "n_records") %in% names(d))) {
      return(data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE))
    }
    d$taxid <- suppressWarnings(as.integer(d$taxid))
    d$n_records <- suppressWarnings(as.integer(d$n_records))
    d
  }
  function(taxid, name) {
    taxid <- suppressWarnings(as.integer(taxid))
    cache <- read_cache()
    i <- match(taxid, cache$taxid)
    if (!is.na(i) && !is.na(cache$n_records[i])) return(cache$n_records[i])
    n <- do.call(count_fn, c(list(clade_query(name)), extra))
    n <- suppressWarnings(as.integer(n))
    if (is.na(n)) return(NA_integer_)
    dir.create(path.expand(dir), recursive = TRUE, showWarnings = FALSE)
    cache <- cache[cache$taxid != taxid | is.na(cache$taxid), , drop = FALSE]
    cache <- rbind(cache, data.frame(taxid = taxid, n_records = n, stringsAsFactors = FALSE))
    utils::write.table(cache, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    n
  }
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-counts.R")'`
Expected: PASS, 6 tests.

- [ ] **Step 5: Commit**

```bash
git add R/reference_counts.R tests/testthat/test-reference-counts.R
git commit -m "feat(reference): cached GenBank record counts"
```

---

### Task 4: Lineage ladders and clade grouping

**Files:**
- Create: `R/reference_groups.R`
- Test: `tests/testthat/test-reference-groups.R`

**Interfaces:**
- Consumes: `RRNA_LINEAGE_RANKS`, `parse_rank_names()`, `parse_rank_taxids()` from `R/taxonomy_ranks.R`; `cached_count_fn()` shape from Task 3 (`function(taxid, name) -> integer`).
- Produces:
  - `sample_ladders(con)` returns `data.frame(ID, level, rank, taxid, name)`, one row per resolved rank per sample, `level` 1..7 matching `RRNA_LINEAGE_RANKS` order, unresolved ranks dropped.
  - `ref_id_from_name(name)` returns a lowercase, underscore-separated slug.
  - `suggest_reference_groups(ladders, count_fn, min_records = 1000, max_records = 25000, start_rank = "family")` returns `list(groups = data.frame(ref_id, clade_taxid, clade_name, clade_rank, n_records, n_samples), members = data.frame(ID, ref_id))`.
  - `group_ladder(ladders, ids, count_fn = NULL)` returns `data.frame(level, rank, taxid, name, n_samples, n_records)` for the shared lineage of `ids`, deepest last; `n_records` is `NA` when `count_fn` is `NULL`.
  - `ladder_coverage(ladders, taxid)` returns the character vector of sample IDs whose lineage contains `taxid`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-groups.R`:

```r
# Clade grouping: walk each sample's lineage up until the clade has enough
# GenBank records, then merge samples that land on the same clade.

# Rank-aligned kingdom..species lineages, as samples.lineage / lineage_taxids
# hold them (7 fields, ";"-joined, empty taxid at a filled rank).
FISH <- list(
  names = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
  taxids = "33208;7711;186623;7952;1002673;7954;7955")
FISH2 <- list(
  names = "Metazoa;Chordata;Actinopteri;Perciformes;Percidae;Perca;Perca flavescens",
  taxids = "33208;7711;186623;8111;8168;8166;8167")
STAR <- list(
  names = "Metazoa;Echinodermata;Asteroidea;Forcipulatida;Asteriidae;Asterias;Asterias rubens",
  taxids = "33208;7586;7588;41692;7607;7608;7610")
STAR2 <- list(
  names = "Metazoa;Echinodermata;Asteroidea;Valvatida;Asterinidae;Patiria;Patiria miniata",
  taxids = "33208;7586;7588;41695;41693;104852;46514")

mk_ladder_db <- function(rows) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT, ref_id TEXT)")
  for (nm in names(rows)) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,?,?,?,1,2,NULL)",
      params = list(nm, NA_integer_, rows[[nm]]$names, rows[[nm]]$taxids))
  }
  con
}

# counts by taxid; anything absent counts as 0 records
counter <- function(map) function(taxid, name) {
  v <- map[[as.character(taxid)]]
  if (is.null(v)) 0L else as.integer(v)
}

test_that("sample_ladders yields one row per resolved rank, deepest last", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_equal(nrow(l), 7L)
  expect_equal(l$rank[1], "kingdom")
  expect_equal(l$rank[7], "species")
  expect_equal(l$name[3], "Actinopteri")
  expect_equal(l$taxid[3], 186623L)
})

test_that("a rank with no taxid is dropped from the ladder", {
  gappy <- list(names = "Metazoa;Chordata;unclassified Chordata;Cypriniformes;Danionidae;Danio;Danio rerio",
    taxids = "33208;7711;;7952;1002673;7954;7955")
  con <- mk_ladder_db(list(F1 = gappy)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_false("class" %in% l$rank)
  expect_equal(nrow(l), 6L)
})

test_that("fish and sea stars land in two groups", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2, S1 = STAR, S2 = STAR2))
  on.exit(DBI::dbDisconnect(con))
  # families and orders are record-poor here; the two classes clear the floor
  cf <- counter(list("186623" = 30000, "7588" = 4102, "7711" = 900000, "7586" = 60000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 50000)
  expect_equal(nrow(res$groups), 2L)
  expect_setequal(res$groups$clade_name, c("Actinopteri", "Asteroidea"))
  expect_equal(res$members$ref_id[res$members$ID == "F1"], "actinopteri")
  expect_equal(res$members$ref_id[res$members$ID == "S2"], "asteroidea")
  expect_equal(res$groups$n_samples[res$groups$clade_name == "Actinopteri"], 2L)
})

test_that("a record-rich family stops the walk at family", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("1002673" = 5000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(res$groups$clade_name, "Danionidae")
  expect_equal(res$groups$clade_rank, "family")
})

test_that("the ceiling stops the walk before an over-broad clade", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  # nothing reaches the floor until kingdom, which is far over the ceiling
  cf <- counter(list("33208" = 900000, "7711" = 400000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  # the walk keeps the best clade it saw rather than one above the ceiling
  expect_equal(res$groups$clade_rank, "phylum")
  expect_equal(res$groups$clade_name, "Chordata")
})

test_that("a nested group merges into its parent when the parent is under the ceiling", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2)); on.exit(DBI::dbDisconnect(con))
  # Danionidae alone clears the floor, Perciformes only does at class level
  cf <- counter(list("1002673" = 2000, "186623" = 20000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(nrow(res$groups), 1L)
  expect_equal(res$groups$clade_name, "Actinopteri")
  expect_equal(res$groups$n_samples, 2L)
})

test_that("a nested group stays separate when the parent is over the ceiling", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("1002673" = 2000, "186623" = 90000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(nrow(res$groups), 2L)
})

test_that("group_ladder reports coverage for the shared lineage", {
  con <- mk_ladder_db(list(S1 = STAR, S2 = STAR2)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  g <- RiboPilot:::group_ladder(l, c("S1", "S2"))
  expect_equal(g$name[nrow(g)], "Asteroidea")
  expect_equal(g$n_samples[nrow(g)], 2L)
  expect_true(all(is.na(g$n_records)))
})

test_that("ladder_coverage names the samples under a taxid", {
  con <- mk_ladder_db(list(F1 = FISH, S1 = STAR, S2 = STAR2)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_setequal(RiboPilot:::ladder_coverage(l, 7588L), c("S1", "S2"))
  expect_setequal(RiboPilot:::ladder_coverage(l, 41692L), "S1")
})

test_that("samples with no usable lineage are reported, not grouped", {
  con <- mk_ladder_db(list(F1 = FISH))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO samples VALUES('X1',NULL,NULL,NULL,0,1,NULL)")
  l <- RiboPilot:::sample_ladders(con)
  expect_false("X1" %in% l$ID)
  res <- RiboPilot:::suggest_reference_groups(l, counter(list("1002673" = 5000)))
  expect_false("X1" %in% res$members$ID)
})

test_that("ref_id_from_name makes a filename-safe slug", {
  expect_equal(RiboPilot:::ref_id_from_name("Asteroidea"), "asteroidea")
  expect_equal(RiboPilot:::ref_id_from_name("Danio rerio"), "danio_rerio")
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-groups.R")'`
Expected: FAIL, `sample_ladders` not found.

- [ ] **Step 3: Write the implementation**

Create `R/reference_groups.R`:

```r
# Turning frozen sample taxonomy into clade reference-database groups.
#
# A seed database is only useful in a narrow band: too narrow and GenBank has
# nothing to seed from, too broad and the seeds stop being specific to the
# organism being assembled. The walk therefore climbs each sample's lineage
# until the clade clears a record floor, and refuses to climb past a ceiling.
# All of this is pure given a count function, which is what makes it testable
# without touching NCBI.

#' Rank-aligned ladders for every sample with a usable lineage
#'
#' One row per resolved rank per sample, kingdom first and species last. Ranks
#' the taxonomy could not resolve carry no taxid and are dropped, so a sample
#' with an unassigned class simply has no class rung.
#'
#' @param con Open project DB connection.
#' @return data.frame(ID, level, rank, taxid, name).
#' @noRd
sample_ladders <- function(con) {
  empty <- data.frame(ID = character(0), level = integer(0), rank = character(0),
    taxid = integer(0), name = character(0), stringsAsFactors = FALSE)
  s <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ID, lineage, lineage_taxids FROM samples
     WHERE lineage_taxids IS NOT NULL AND TRIM(lineage_taxids) != ''"),
    error = function(e) NULL)
  if (is.null(s) || nrow(s) == 0) return(empty)
  rows <- lapply(seq_len(nrow(s)), function(i) {
    tx <- parse_rank_taxids(s$lineage_taxids[i])
    nm <- parse_rank_names(s$lineage[i])
    keep <- which(!is.na(tx) & !is.na(nm))
    if (length(keep) == 0) return(NULL)
    data.frame(ID = s$ID[i], level = as.integer(keep), rank = RRNA_LINEAGE_RANKS[keep],
      taxid = unname(tx[keep]), name = unname(nm[keep]), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Filename-safe group key from a clade name
#' @noRd
ref_id_from_name <- function(name) {
  x <- tolower(trimws(name %||% ""))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

#' Sample IDs whose lineage passes through a taxid
#' @param ladders From [sample_ladders()]. @param taxid NCBI taxid.
#' @noRd
ladder_coverage <- function(ladders, taxid) {
  unique(ladders$ID[ladders$taxid == as.integer(taxid)])
}

#' The lineage shared by a set of samples, with coverage and optional counts
#'
#' @param ladders From [sample_ladders()].
#' @param ids Sample IDs in the group.
#' @param count_fn Optional function(taxid, name) -> integer record count.
#' @return data.frame(level, rank, taxid, name, n_samples, n_records), deepest last.
#' @noRd
group_ladder <- function(ladders, ids, count_fn = NULL) {
  l <- ladders[ladders$ID %in% ids, , drop = FALSE]
  empty <- data.frame(level = integer(0), rank = character(0), taxid = integer(0),
    name = character(0), n_samples = integer(0), n_records = integer(0),
    stringsAsFactors = FALSE)
  if (nrow(l) == 0) return(empty)
  # a rung belongs to the group's ladder when every member passes through it,
  # plus the deeper rungs of a single member so a one-sample group still has
  # somewhere to go down to
  agg <- stats::aggregate(list(n_samples = l$ID), by = list(taxid = l$taxid), FUN = function(x) length(unique(x)))
  key <- l[!duplicated(l$taxid), c("level", "rank", "taxid", "name")]
  out <- merge(key, agg, by = "taxid")
  out <- out[order(out$level), , drop = FALSE]
  out <- out[out$n_samples == max(out$n_samples) | out$n_samples == length(unique(l$ID)), , drop = FALSE]
  out$n_records <- NA_integer_
  if (!is.null(count_fn)) {
    out$n_records <- vapply(seq_len(nrow(out)),
      function(i) as.integer(count_fn(out$taxid[i], out$name[i]) %|NA|% NA_integer_), integer(1))
  }
  rownames(out) <- NULL
  out[, c("level", "rank", "taxid", "name", "n_samples", "n_records"), drop = FALSE]
}

#' Walk one sample's ladder up to the first clade clearing the record floor
#'
#' Returns the deepest rung at or above `start_rank` whose count reaches
#' `min_records`. If nothing does, the deepest rung whose count stays at or
#' below `max_records` is kept: a database that is too small is still usable,
#' one that is far too broad is not.
#' @return one-row data.frame(level, rank, taxid, name, n_records), or NULL.
#' @noRd
walk_to_clade <- function(rungs, count_fn, min_records, max_records, start_rank) {
  start <- match(start_rank, RRNA_LINEAGE_RANKS)
  cand <- rungs[rungs$level <= start, , drop = FALSE]
  cand <- cand[order(-cand$level), , drop = FALSE]   # deepest first, walking up
  if (nrow(cand) == 0) return(NULL)
  best <- NULL
  for (i in seq_len(nrow(cand))) {
    n <- suppressWarnings(as.integer(count_fn(cand$taxid[i], cand$name[i])))
    if (is.na(n)) next
    row <- data.frame(level = cand$level[i], rank = cand$rank[i], taxid = cand$taxid[i],
      name = cand$name[i], n_records = n, stringsAsFactors = FALSE)
    if (n <= max_records) best <- row
    if (n >= min_records) return(if (n <= max_records) row else (best %||% row))
  }
  best %||% NULL
}

#' Suggest clade groups for a project's samples
#'
#' @param ladders From [sample_ladders()].
#' @param count_fn function(taxid, name) -> integer GenBank record count.
#' @param min_records,max_records Record window bounding the walk.
#' @param start_rank Rank the walk starts from (default "family").
#' @return list(groups = data.frame(ref_id, clade_taxid, clade_name, clade_rank,
#'   n_records, n_samples), members = data.frame(ID, ref_id)).
#' @noRd
suggest_reference_groups <- function(ladders, count_fn, min_records = 1000,
                                     max_records = 25000, start_rank = "family") {
  ids <- unique(ladders$ID)
  empty <- list(
    groups = data.frame(ref_id = character(0), clade_taxid = integer(0),
      clade_name = character(0), clade_rank = character(0), n_records = integer(0),
      n_samples = integer(0), stringsAsFactors = FALSE),
    members = data.frame(ID = character(0), ref_id = character(0), stringsAsFactors = FALSE))
  if (length(ids) == 0) return(empty)

  picks <- lapply(ids, function(id) {
    p <- walk_to_clade(ladders[ladders$ID == id, , drop = FALSE], count_fn,
      min_records, max_records, start_rank)
    if (is.null(p)) return(NULL)
    p$ID <- id
    p
  })
  picks <- Filter(Negate(is.null), picks)
  if (length(picks) == 0) return(empty)
  picks <- do.call(rbind, picks)

  # Fold a group into an ancestor group when one exists and stays under the
  # ceiling: two fish families are better served by one class database than by
  # two family databases that each miss the other's reads.
  chosen <- unique(picks[, c("level", "rank", "taxid", "name", "n_records")])
  for (i in seq_len(nrow(picks))) {
    anc <- ladders$taxid[ladders$ID == picks$ID[i]]
    up <- chosen[chosen$taxid %in% anc & chosen$level < picks$level[i] &
        chosen$n_records <= max_records, , drop = FALSE]
    if (nrow(up) > 0) {
      up <- up[order(-up$level), , drop = FALSE][1, ]
      picks[i, c("level", "rank", "taxid", "name", "n_records")] <- up
    }
  }

  picks$ref_id <- ref_id_from_name(picks$name)
  members <- picks[, c("ID", "ref_id")]
  rownames(members) <- NULL
  g <- picks[!duplicated(picks$ref_id), , drop = FALSE]
  groups <- data.frame(
    ref_id = g$ref_id, clade_taxid = as.integer(g$taxid), clade_name = g$name,
    clade_rank = g$rank, n_records = as.integer(g$n_records),
    n_samples = as.integer(table(members$ref_id)[g$ref_id]),
    stringsAsFactors = FALSE)
  rownames(groups) <- NULL
  list(groups = groups, members = members)
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-groups.R")'`
Expected: PASS, 11 tests. If the merge step or `group_ladder` filtering does not satisfy a case, fix the implementation, not the test: the tests encode the spec's stated behavior.

- [ ] **Step 5: Commit**

```bash
git add R/reference_groups.R tests/testthat/test-reference-groups.R
git commit -m "feat(reference): lineage ladders and clade grouping"
```

---

### Task 5: Persisting groups, building, and wiring to assemble_opts

**Files:**
- Create: `R/reference_build.R`
- Test: `tests/testthat/test-reference-build.R`

**Interfaces:**
- Consumes: `reference` / `reference_config` / `samples.ref_id` from Task 1; `ref_config()`, `append_ref_manifest()` from Task 2; `clade_query()`, `cached_count_fn()` from Task 3; `sample_ladders()`, `suggest_reference_groups()`, `ref_id_from_name()` from Task 4; `build_reference_db()` from `R/build_reference_db.R`.
- Produces:
  - `write_reference_groups(con, groups, members)` replaces the suggestion, preserving `reference_switch`, `source`, `seeds_db`, `labels_db` for any `ref_id` whose `clade_taxid` is unchanged; returns the number of groups invisibly.
  - `set_reference_clade(con, ref_id, clade_taxid, clade_name, clade_rank, n_records, member_ids)` writes one group and its membership (used by the ladder modal and by splits); returns the new `ref_id` invisibly.
  - `skip_reference_group(con, ref_id)` sets `reference_switch = 5`, `source = 'skipped'`.
  - `wire_reference_db(con, ref_id, seeds, labels, source = "built")` creates or updates the `assemble_opts` preset named `ref_id`, moves the group's samples onto it, re-queues them (`assemble_switch = 1` where unlocked), and sets `reference_switch = 2`; returns the number of samples moved invisibly.
  - `start_reference_build(con, project, ref_id, api_key = NULL, refseq = FALSE, conda_env = NULL, overwrite = FALSE, launcher = processx_launcher)` writes `reference_switch = 4`, `pid`, `log_file`; returns the log path invisibly.
  - `poll_reference_builds(con, alive_fn = process_alive)` finishes any `reference_switch = 4` row whose process has exited, wiring successes and recording failures; returns the number of rows changed invisibly.
  - `reference_pending(con)` returns the `ref_id`s that are neither done nor skipped.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-build.R`:

```r
# Persisting clade groups, wiring a built database into assemble_opts, and
# finishing background builds.

mk_ref_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT, ref_id TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT PRIMARY KEY, cpus INT,
    memory INT, getOrganelle TEXT, seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INT,
    max_scaffolds INT, min_assembly_length INT, assembler TEXT, bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, NULL, 1000, 25000)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default',4,16,'--flag','old_s.fa','old_l.fa',NULL,10,10,1000,'getorganelle','--local','--isolate')")
  for (id in c("S1", "S2")) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,NULL,NULL,NULL,1,2,NULL)", params = list(id))
    DBI::dbExecute(con, "INSERT INTO assemble VALUES(?,2,0,'default','default')", params = list(id))
  }
  con
}

grp <- data.frame(ref_id = "asteroidea", clade_taxid = 7588L, clade_name = "Asteroidea",
  clade_rank = "class", n_records = 4102L, n_samples = 2L, stringsAsFactors = FALSE)
mem <- data.frame(ID = c("S1", "S2"), ref_id = "asteroidea", stringsAsFactors = FALSE)

test_that("write_reference_groups stores groups and stamps membership", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(nrow(r), 1L)
  expect_equal(r$clade_name[1], "Asteroidea")
  expect_equal(r$reference_switch[1], 1L)
  expect_equal(r$query[1], "\"Asteroidea\"[Organism]")
  expect_setequal(DBI::dbGetQuery(con, "SELECT ref_id FROM samples")$ref_id, "asteroidea")
})

test_that("re-suggesting keeps a built group's database and state", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  RiboPilot:::write_reference_groups(con, grp, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 2L)
  expect_equal(r$seeds_db[1], "s.fa")
})

test_that("re-suggesting resets a group whose clade changed", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  moved <- grp; moved$clade_taxid <- 7586L
  RiboPilot:::write_reference_groups(con, moved, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 1L)
  expect_true(is.na(r$seeds_db[1]))
})

test_that("wire_reference_db creates a preset and moves the group's samples", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  n <- RiboPilot:::wire_reference_db(con, "asteroidea", "/refs/a_seeds.fa", "/refs/a_labels.fa")
  expect_equal(n, 2L)
  o <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'asteroidea'")
  expect_equal(nrow(o), 1L)
  expect_equal(o$seeds_db[1], "/refs/a_seeds.fa")
  expect_equal(o$labels_db[1], "/refs/a_labels.fa")
  expect_equal(o$cpus[1], 4)              # copied from the samples' previous preset
  expect_equal(o$getOrganelle[1], "--flag")
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts,
    "asteroidea")
  expect_equal(DBI::dbGetQuery(con, "SELECT reference_switch FROM reference")$reference_switch[1], 2L)
})

test_that("wiring re-queues unlocked samples and leaves locked ones alone", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'S2'")
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  a <- DBI::dbGetQuery(con, "SELECT ID, assemble_switch, assemble_opts FROM assemble ORDER BY ID")
  expect_equal(a$assemble_switch[a$ID == "S1"], 1L)
  expect_equal(a$assemble_switch[a$ID == "S2"], 2L)
  expect_equal(a$assemble_opts[a$ID == "S2"], "default")
})

test_that("skip leaves the samples' preset untouched", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::skip_reference_group(con, "asteroidea")
  r <- DBI::dbGetQuery(con, "SELECT reference_switch, source FROM reference")
  expect_equal(r$reference_switch[1], 5L)
  expect_equal(r$source[1], "skipped")
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("set_reference_clade splits a subset into its own group", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::set_reference_clade(con, ref_id = NULL, clade_taxid = 41692L,
    clade_name = "Forcipulatida", clade_rank = "order", n_records = 900L,
    member_ids = "S1")
  r <- DBI::dbGetQuery(con, "SELECT ref_id, n_samples FROM reference ORDER BY ref_id")
  expect_setequal(r$ref_id, c("asteroidea", "forcipulatida"))
  expect_equal(r$n_samples[r$ref_id == "forcipulatida"], 1L)
  expect_equal(r$n_samples[r$ref_id == "asteroidea"], 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT ref_id FROM samples WHERE ID = 'S1'")$ref_id, "forcipulatida")
})

test_that("start_reference_build records the job and poll finishes it", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = refs)
  RiboPilot:::write_reference_groups(con, grp, mem)

  seeds <- file.path(refs, "seeds", "asteroidea_rrna_seeds.fasta")
  labels <- file.path(refs, "labels", "asteroidea_rrna_labels.fasta")
  launcher <- function(cmd, args, log) {
    dir.create(dirname(seeds), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(labels), recursive = TRUE, showWarnings = FALSE)
    writeLines(">x", seeds); writeLines(">18S", labels)
    writeLines("done", log)
    4242L
  }
  RiboPilot:::start_reference_build(con, proj, "asteroidea", launcher = launcher)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 4L)
  expect_equal(r$pid[1], 4242L)

  n <- RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) FALSE)
  expect_equal(n, 1L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 2L)
  expect_equal(r$seeds_db[1], seeds)
  expect_equal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble WHERE ID='S1'")$assemble_opts,
    "asteroidea")
  expect_equal(nrow(RiboPilot:::read_ref_manifest(refs)), 1L)
})

test_that("a build that produced no FASTAs is recorded as failed", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = withr::local_tempdir())
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::start_reference_build(con, proj, "asteroidea",
    launcher = function(cmd, args, log) { writeLines("esearch: command not found", log); 7L })
  RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) FALSE)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 3L)
  expect_true(nzchar(r$fail_reason[1]))
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("a still-running build is left alone", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_ref_config(con, ref_dir = withr::local_tempdir())
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::start_reference_build(con, withr::local_tempdir(), "asteroidea",
    launcher = function(cmd, args, log) { writeLines("working", log); 11L })
  expect_equal(RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) TRUE), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT reference_switch FROM reference")$reference_switch[1], 4L)
})

test_that("reference_pending lists everything not done or skipped", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  expect_equal(RiboPilot:::reference_pending(con), "asteroidea")
  RiboPilot:::skip_reference_group(con, "asteroidea")
  expect_equal(length(RiboPilot:::reference_pending(con)), 0L)
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-build.R")'`
Expected: FAIL, `write_reference_groups` not found.

- [ ] **Step 3: Write the implementation**

Create `R/reference_build.R`:

```r
# Clade groups in the project database: how they are stored, built, and wired
# into the assembly presets the pipeline actually reads.
#
# A build is a detached process rather than an in-app wait: downloading tens of
# thousands of GenBank records takes minutes, and the app has to stay usable.
# The database row is the only state that survives a closed app, so the row
# carries the pid and log path and the finishing logic is a poll, not a callback.

#' Expected seed/label paths for a clade in a reference folder
#' @return list(seeds, labels).
#' @noRd
reference_db_paths <- function(dir, ref_id) {
  prefix <- paste0(ref_id, "_rrna")
  list(
    seeds = file.path(path.expand(dir), "seeds", paste0(prefix, "_seeds.fasta")),
    labels = file.path(path.expand(dir), "labels", paste0(prefix, "_labels.fasta")))
}

#' Store a suggestion, keeping what an unchanged group already achieved
#'
#' Re-suggesting must not throw away a database somebody waited ten minutes for,
#' so a group whose clade taxid is unchanged keeps its state and paths; a group
#' pointing somewhere new starts over.
#'
#' @param con Open project DB connection.
#' @param groups data.frame from [suggest_reference_groups()]$groups.
#' @param members data.frame(ID, ref_id).
#' @return (invisibly) the number of groups written.
#' @noRd
write_reference_groups <- function(con, groups, members) {
  old <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM reference"), error = function(e) NULL)
  DBI::dbExecute(con, "DELETE FROM reference")
  now <- as.integer(Sys.time())
  for (i in seq_len(nrow(groups))) {
    keep <- NULL
    if (!is.null(old) && nrow(old) > 0) {
      j <- which(old$ref_id == groups$ref_id[i] &
          !is.na(old$clade_taxid) & old$clade_taxid == groups$clade_taxid[i])
      if (length(j) > 0) keep <- old[j[1], ]
    }
    DBI::dbExecute(con,
      "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank, query, n_samples,
         n_records, reference_switch, source, seeds_db, labels_db, log_file, pid, fail_reason,
         time_stamp)
       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
      params = list(groups$ref_id[i], as.integer(groups$clade_taxid[i]), groups$clade_name[i],
        groups$clade_rank[i], clade_query(groups$clade_name[i]),
        as.integer(groups$n_samples[i]), as.integer(groups$n_records[i]),
        if (is.null(keep)) 1L else as.integer(keep$reference_switch),
        if (is.null(keep)) NA_character_ else keep$source,
        if (is.null(keep)) NA_character_ else keep$seeds_db,
        if (is.null(keep)) NA_character_ else keep$labels_db,
        if (is.null(keep)) NA_character_ else keep$log_file,
        if (is.null(keep)) NA_integer_ else keep$pid,
        if (is.null(keep)) NA_character_ else keep$fail_reason,
        now))
  }
  DBI::dbExecute(con, "UPDATE samples SET ref_id = NULL")
  if (nrow(members) > 0) {
    upd <- DBI::dbSendStatement(con, "UPDATE samples SET ref_id = ? WHERE ID = ?")
    on.exit(DBI::dbClearResult(upd), add = TRUE)
    DBI::dbBind(upd, list(members$ref_id, members$ID))
  }
  invisible(nrow(groups))
}

#' Recount `reference.n_samples` from current membership, dropping empty groups
#' @noRd
refresh_reference_counts <- function(con) {
  DBI::dbExecute(con,
    "UPDATE reference SET n_samples =
       (SELECT COUNT(*) FROM samples s WHERE s.ref_id = reference.ref_id)")
  DBI::dbExecute(con, "DELETE FROM reference WHERE COALESCE(n_samples, 0) = 0")
  invisible(NULL)
}

#' Point a set of samples at a chosen clade, creating the group if needed
#'
#' @param con Open project DB connection.
#' @param ref_id Group key, or NULL to derive one from `clade_name`.
#' @param clade_taxid,clade_name,clade_rank,n_records The chosen NCBI node.
#' @param member_ids Sample IDs to move into the group.
#' @return (invisibly) the group key.
#' @noRd
set_reference_clade <- function(con, ref_id = NULL, clade_taxid, clade_name, clade_rank,
                               n_records = NA_integer_, member_ids = character(0)) {
  key <- ref_id %||% ref_id_from_name(clade_name)
  exists_row <- nrow(DBI::dbGetQuery(con, "SELECT 1 FROM reference WHERE ref_id = ? LIMIT 1",
    params = list(key))) > 0
  now <- as.integer(Sys.time())
  if (exists_row) {
    # a clade change invalidates whatever was built for the old one
    DBI::dbExecute(con,
      "UPDATE reference SET clade_taxid = ?, clade_name = ?, clade_rank = ?, query = ?,
         n_records = ?, reference_switch = 1, source = NULL, seeds_db = NULL, labels_db = NULL,
         fail_reason = NULL, time_stamp = ?
       WHERE ref_id = ? AND (clade_taxid IS NULL OR clade_taxid != ?)",
      params = list(as.integer(clade_taxid), clade_name, clade_rank, clade_query(clade_name),
        as.integer(n_records), now, key, as.integer(clade_taxid)))
  } else {
    DBI::dbExecute(con,
      "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank, query, n_samples,
         n_records, reference_switch, time_stamp)
       VALUES (?,?,?,?,?,0,?,1,?)",
      params = list(key, as.integer(clade_taxid), clade_name, clade_rank,
        clade_query(clade_name), as.integer(n_records), now))
  }
  if (length(member_ids) > 0) {
    upd <- DBI::dbSendStatement(con, "UPDATE samples SET ref_id = ? WHERE ID = ?")
    on.exit(DBI::dbClearResult(upd), add = TRUE)
    DBI::dbBind(upd, list(rep(key, length(member_ids)), member_ids))
  }
  refresh_reference_counts(con)
  invisible(key)
}

#' Mark a group skipped (its samples keep the database they already have)
#' @noRd
skip_reference_group <- function(con, ref_id) {
  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 5, source = 'skipped', fail_reason = NULL,
       time_stamp = ? WHERE ref_id = ?",
    params = list(as.integer(Sys.time()), ref_id))
  invisible(NULL)
}

#' Give a group's samples an assemble_opts preset pointing at its database
#'
#' The preset is named after the group and copies every other setting from the
#' preset the group's samples are on, so switching database changes the database
#' and nothing else. Locked samples keep their preset: their stored results were
#' produced with the old one.
#'
#' @param con Open project DB connection.
#' @param ref_id Group key.
#' @param seeds,labels Absolute paths to the seed and label FASTAs.
#' @param source "built" or "existing".
#' @return (invisibly) the number of samples moved.
#' @noRd
wire_reference_db <- function(con, ref_id, seeds, labels, source = "built") {
  ids <- DBI::dbGetQuery(con,
    "SELECT s.ID FROM samples s JOIN assemble a ON a.ID = s.ID
     WHERE s.ref_id = ? AND COALESCE(a.assemble_lock, 0) = 0",
    params = list(ref_id))$ID
  base <- DBI::dbGetQuery(con,
    "SELECT o.* FROM assemble_opts o JOIN assemble a ON a.assemble_opts = o.assemble_opts
     JOIN samples s ON s.ID = a.ID WHERE s.ref_id = ? LIMIT 1", params = list(ref_id))
  if (nrow(base) == 0) {
    base <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'default'")
  }
  if (nrow(base) > 0) {
    row <- base[1, , drop = FALSE]
    row$assemble_opts <- ref_id
    row$seeds_db <- seeds
    row$labels_db <- labels
    DBI::dbExecute(con, "DELETE FROM assemble_opts WHERE assemble_opts = ?", params = list(ref_id))
    DBI::dbAppendTable(con, "assemble_opts", row)
  }
  if (length(ids) > 0) {
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbExecute(con,
      sprintf("UPDATE assemble SET assemble_opts = ?, assemble_switch = 1, time_stamp = ?
               WHERE ID IN (%s)", ph),
      params = c(list(ref_id, as.integer(Sys.time())), as.list(ids)))
  }
  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 2, source = ?, seeds_db = ?, labels_db = ?,
       fail_reason = NULL, pid = NULL, time_stamp = ? WHERE ref_id = ?",
    params = list(source, seeds, labels, as.integer(Sys.time()), ref_id))
  invisible(length(ids))
}

#' Point a group at seed/label FASTAs that already exist
#' @return (invisibly) the number of samples moved.
#' @noRd
use_existing_reference <- function(con, ref_id, seeds, labels) {
  seeds <- path.expand(seeds); labels <- path.expand(labels)
  for (f in c(seeds, labels)) {
    if (!file.exists(f)) stop("Reference file not found: ", f, call. = FALSE)
  }
  wire_reference_db(con, ref_id, seeds, labels, source = "existing")
}

#' Launch a detached command, appending its output to `log`
#'
#' Split out so tests can substitute a launcher that never spawns anything.
#' @return the pid.
#' @noRd
processx_launcher <- function(cmd, args, log) {
  p <- processx::process$new(cmd, args, stdout = log, stderr = "2>&1", cleanup = FALSE)
  p$get_pid()
}

#' Is a pid still running?
#' @noRd
process_alive <- function(pid) {
  if (is.na(pid)) return(FALSE)
  isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(as.integer(pid))), error = function(e) FALSE))
}

#' Start a background build for one clade group
#'
#' @param con Open project DB connection.
#' @param project Project directory (the log lives in `<project>/.reference`).
#' @param ref_id Group key.
#' @param api_key,refseq,conda_env,overwrite Passed to [build_reference_db()].
#' @param launcher function(cmd, args, log) -> pid.
#' @return (invisibly) the log path.
#' @noRd
start_reference_build <- function(con, project, ref_id, api_key = NULL, refseq = FALSE,
                                  conda_env = NULL, overwrite = FALSE,
                                  launcher = processx_launcher) {
  row <- DBI::dbGetQuery(con, "SELECT * FROM reference WHERE ref_id = ?", params = list(ref_id))
  if (nrow(row) == 0) stop("No reference group '", ref_id, "'.", call. = FALSE)
  cfg <- ref_config(con)
  logdir <- file.path(project, ".reference")
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(logdir, paste0(ref_id, ".log"))
  unlink(log)

  call <- sprintf(
    "RiboPilot::build_reference_db(query = %s, clade = %s, out_dir = %s, api_key = %s, refseq = %s, conda_env = %s, overwrite = %s)",
    deparse(row$query[1]), deparse(ref_id), deparse(cfg$ref_dir),
    if (is.null(api_key) || !nzchar(api_key)) "NULL" else deparse(api_key),
    if (isTRUE(refseq)) "TRUE" else "FALSE",
    if (is.null(conda_env) || !nzchar(conda_env)) "NULL" else deparse(conda_env),
    if (isTRUE(overwrite)) "TRUE" else "FALSE")
  pid <- launcher(file.path(R.home("bin"), "Rscript"), c("-e", call), log)

  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 4, pid = ?, log_file = ?, fail_reason = NULL,
       time_stamp = ? WHERE ref_id = ?",
    params = list(as.integer(pid), log, as.integer(Sys.time()), ref_id))
  invisible(log)
}

#' Last non-empty lines of a build log, for the failure message
#' @noRd
log_tail <- function(log, n = 3) {
  if (is.null(log) || is.na(log) || !file.exists(log)) return("")
  x <- tryCatch(readLines(log, warn = FALSE), error = function(e) character(0))
  x <- trimws(x[nzchar(trimws(x))])
  if (length(x) == 0) return("")
  paste(utils::tail(x, n), collapse = " | ")
}

#' Finish any background build whose process has exited
#'
#' @param con Open project DB connection.
#' @param alive_fn function(pid) -> logical.
#' @return (invisibly) the number of rows changed.
#' @noRd
poll_reference_builds <- function(con, alive_fn = process_alive) {
  running <- tryCatch(DBI::dbGetQuery(con,
    "SELECT * FROM reference WHERE reference_switch = 4"), error = function(e) data.frame())
  if (nrow(running) == 0) return(invisible(0L))
  cfg <- ref_config(con)
  n <- 0L
  for (i in seq_len(nrow(running))) {
    if (isTRUE(alive_fn(running$pid[i]))) next
    ref_id <- running$ref_id[i]
    p <- reference_db_paths(cfg$ref_dir, ref_id)
    if (file.exists(p$seeds) && file.exists(p$labels)) {
      wire_reference_db(con, ref_id, p$seeds, p$labels, source = "built")
      append_ref_manifest(cfg$ref_dir, clade = ref_id, taxid = running$clade_taxid[i],
        rank = running$clade_rank[i], query = running$query[i],
        n_records = running$n_records[i], seeds = p$seeds, labels = p$labels)
    } else {
      DBI::dbExecute(con,
        "UPDATE reference SET reference_switch = 3, pid = NULL, fail_reason = ?, time_stamp = ?
         WHERE ref_id = ?",
        params = list(paste0("Build produced no seed/label FASTAs. ",
          log_tail(running$log_file[i])), as.integer(Sys.time()), ref_id))
    }
    n <- n + 1L
  }
  invisible(n)
}

#' Groups that are neither built nor skipped
#' @noRd
reference_pending <- function(con) {
  tryCatch(DBI::dbGetQuery(con,
    "SELECT ref_id FROM reference WHERE COALESCE(reference_switch, 0) NOT IN (2, 5)")$ref_id,
    error = function(e) character(0))
}
```

- [ ] **Step 4: Add `ps` to Imports**

`process_alive()` uses the `ps` package, which `processx` already depends on but the package does not import. In `DESCRIPTION`, add `ps,` to `Imports:` in alphabetical position (after `purrr`). Run `Rscript -e 'devtools::document()'` afterwards if NAMESPACE needs regenerating (these helpers are `@noRd`, so it should not).

- [ ] **Step 5: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-build.R")'`
Expected: PASS, 11 tests.

- [ ] **Step 6: Commit**

```bash
git add R/reference_build.R DESCRIPTION tests/testthat/test-reference-build.R
git commit -m "feat(reference): persist, build, and wire clade databases"
```

---

### Task 6: Reference tab UI (table, ladder modal, handlers)

**Files:**
- Create: `R/app_reference.R`
- Test: `tests/testthat/test-reference-modal.R`

**Interfaces:**
- Consumes: everything from Tasks 2 to 5; `opts_help()` from `R/app_opts.R`; `rt_link()`, `rt_dynamicIcon()` from `R/utils_reactable.R`; `tax_search()` from `R/tax_search.R`; `app_con()`, `app_project()`, `app_query()` from `R/app_utils.R`.
- Produces:
  - `reference_table(con)` returns `data.frame(reference_switch, ID, clade_name, clade_rank, n_samples, n_records, source, database, fail_reason, time_stamp)` where `ID` is the `ref_id` (so the app's shared selection helpers work unchanged) and `database` is the seed path or `""`.
  - `reference_ladder_modal(ref_id, clade_name, ladder, chosen_taxid, n_group)` returns a `modalDialog`.
  - `reference_build_modal(ref_ids, ref_dir)` returns a `modalDialog`.
  - `reference_existing_modal(manifest, ref_id)` returns a `modalDialog`.
  - `reference_handlers(input, output, session, deps)` wires every observer and returns `list(open_ladder = function(ref_id), poll = function())`.
    `deps` is `list(con = fn -> DBI con, project = fn -> chr, refresh = reactiveVal, selected_ids = fn -> chr, table_data = reactive)`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-modal.R`. The repo already tests modals by inspecting the returned tag tree (see `test-triage-modal.R`); follow that style.

```r
# The reference tab's table shape and modal content.

mk_tab_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, NULL, 1000, 25000)")
  DBI::dbExecute(con, "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank,
    query, n_samples, n_records, reference_switch, source, seeds_db, labels_db, time_stamp)
    VALUES ('asteroidea', 7588, 'Asteroidea', 'class', '\"Asteroidea\"[Organism]', 4, 4102,
      2, 'built', '/refs/a_seeds.fa', '/refs/a_labels.fa', 1)")
  con
}

test_that("reference_table exposes ref_id as ID and the seed path as database", {
  con <- mk_tab_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::reference_table(con)
  expect_equal(d$ID[1], "asteroidea")
  expect_equal(d$clade_name[1], "Asteroidea")
  expect_equal(d$database[1], "/refs/a_seeds.fa")
  expect_true("reference_switch" %in% names(d))
})

test_that("the ladder modal shows every rung with coverage and counts", {
  ladder <- data.frame(level = 2:5, rank = c("phylum", "class", "order", "family"),
    taxid = c(7586L, 7588L, 41692L, 7607L),
    name = c("Echinodermata", "Asteroidea", "Forcipulatida", "Asteriidae"),
    n_samples = c(4L, 4L, 3L, 3L), n_records = c(31500L, 4102L, 900L, 340L),
    stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_ladder_modal("asteroidea", "Asteroidea", ladder, 7588L, 4L)
  html <- as.character(m)
  expect_true(grepl("Asteriidae", html, fixed = TRUE))
  expect_true(grepl("Echinodermata", html, fixed = TRUE))
  expect_true(grepl("4102", html, fixed = TRUE))
  expect_true(grepl("ref_rung", html, fixed = TRUE))
  expect_true(grepl("ref_search", html, fixed = TRUE))
})

test_that("the ladder modal marks the current clade as chosen", {
  ladder <- data.frame(level = 3L, rank = "class", taxid = 7588L, name = "Asteroidea",
    n_samples = 4L, n_records = 4102L, stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_ladder_modal("asteroidea", "Asteroidea", ladder, 7588L, 4L)
  expect_true(grepl("current", as.character(m), fixed = TRUE))
})

test_that("the build modal names the folder it will write to", {
  m <- RiboPilot:::reference_build_modal(c("asteroidea", "actinopteri"), "/lab/refs")
  html <- as.character(m)
  expect_true(grepl("/lab/refs", html, fixed = TRUE))
  expect_true(grepl("2", html, fixed = TRUE))
  expect_true(grepl("ref_api_key", html, fixed = TRUE))
  expect_true(grepl("ref_conda", html, fixed = TRUE))
  expect_true(grepl("ref_refseq", html, fixed = TRUE))
})

test_that("the existing-database modal lists manifest entries", {
  man <- data.frame(clade = "gastropoda", taxid = 6448L, rank = "class",
    query = "q", n_records = 5000L, built = "2026-07-01 09:00:00",
    seeds = "/refs/seeds/gastropoda_rrna_seeds.fasta",
    labels = "/refs/labels/gastropoda_rrna_labels.fasta", stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("gastropoda", html, fixed = TRUE))
  expect_true(grepl("5000", html, fixed = TRUE))
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
})

test_that("the existing-database modal still offers typed paths with an empty manifest", {
  man <- RiboPilot:::read_ref_manifest(withr::local_tempdir())
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
  expect_true(grepl("ref_existing_labels", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-modal.R")'`
Expected: FAIL, `reference_table` not found.

- [ ] **Step 3: Write the table and the modals**

Create `R/app_reference.R` with the table builder and the three modals:

```r
# The Reference tab: one row per clade group, a ladder to change the clade, and
# the build / reuse / skip actions.
#
# The table's key column is named ID so the app's shared selection, CSV export
# and click handlers work on it unchanged; the value is the group key, not a
# sample.

#' Reference-tab table data
#' @param con Open project DB connection.
#' @return data.frame with the columns the tab renders.
#' @noRd
reference_table <- function(con) {
  d <- tryCatch(DBI::dbGetQuery(con,
    "SELECT reference_switch, ref_id AS ID, clade_name, clade_rank, n_samples, n_records,
       COALESCE(source, '') AS source, COALESCE(seeds_db, '') AS database,
       COALESCE(fail_reason, '') AS fail_reason, time_stamp
     FROM reference ORDER BY ref_id"), error = function(e) NULL)
  if (is.null(d)) {
    d <- data.frame(reference_switch = integer(0), ID = character(0), clade_name = character(0),
      clade_rank = character(0), n_samples = integer(0), n_records = integer(0),
      source = character(0), database = character(0), fail_reason = character(0),
      time_stamp = integer(0), stringsAsFactors = FALSE)
  }
  d
}

#' One rung of the clade ladder, as a radio choice label
#' @noRd
rung_label <- function(name, rank, n_samples, n_records, chosen) {
  recs <- if (is.na(n_records)) "records unknown" else paste0(format(n_records, big.mark = ","), " records")
  sprintf("%s [%s] - %d sample(s), %s%s", name, rank, n_samples, recs,
    if (isTRUE(chosen)) "  (current)" else "")
}

#' Clade ladder modal: pick a rung, or search the taxdump for anything else
#'
#' @param ref_id Group key. @param clade_name Current clade.
#' @param ladder data.frame from [group_ladder()].
#' @param chosen_taxid The group's current clade taxid.
#' @param n_group Samples in the group.
#' @noRd
reference_ladder_modal <- function(ref_id, clade_name, ladder, chosen_taxid, n_group) {
  ladder <- ladder[order(-ladder$level), , drop = FALSE]
  choices <- stats::setNames(as.character(ladder$taxid),
    vapply(seq_len(nrow(ladder)), function(i) rung_label(ladder$name[i], ladder$rank[i],
      ladder$n_samples[i], ladder$n_records[i],
      !is.na(chosen_taxid) && ladder$taxid[i] == chosen_taxid), character(1)))
  shiny::modalDialog(
    title = sprintf("Reference clade for '%s' (%d sample(s))", ref_id, n_group), size = "l",
    shinyWidgets::prettyRadioButtons("ref_rung", "Shared lineage (deepest first):",
      choices = choices,
      selected = if (!is.na(chosen_taxid)) as.character(chosen_taxid) else NULL,
      shape = "round", status = "primary"),
    opts_help("A clade with too few records makes a weak seed database; one with ",
      "too many stops being specific. Picking a rung that covers only some of the ",
      "group's samples offers to split the rest into their own group."),
    shiny::textInput("ref_search", "Or search NCBI for another clade:", width = "100%",
      placeholder = "e.g. Echinodermata"),
    shiny::uiOutput("ref_search_results"),
    footer = shiny::tagList(
      shiny::actionButton("ref_clade_apply", "Use this clade", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Build modal: the options a build run takes
#' @param ref_ids Groups that will be built. @param ref_dir Destination folder.
#' @noRd
reference_build_modal <- function(ref_ids, ref_dir) {
  shiny::modalDialog(
    title = sprintf("Build %d reference database(s)", length(ref_ids)), size = "m",
    shiny::tags$p(class = "text-muted", style = "font-size:12px;",
      "Downloads GenBank records and writes seed + label FASTAs to:"),
    shiny::tags$pre(style = "font-size:12px;", ref_dir),
    shiny::tags$pre(style = "max-height:120px; overflow:auto; font-size:12px;",
      paste(ref_ids, collapse = "\n")),
    shiny::textInput("ref_api_key", "NCBI API key (optional):", width = "100%",
      value = Sys.getenv("NCBI_API_KEY", "")),
    shiny::textInput("ref_conda", "Conda environment providing Entrez Direct (optional):",
      width = "100%", value = ""),
    shinyWidgets::prettyCheckbox("ref_refseq", "RefSeq records only", value = FALSE,
      status = "primary", shape = "curve"),
    shinyWidgets::prettyCheckbox("ref_overwrite", "Rebuild even if the database exists",
      value = FALSE, status = "primary", shape = "curve"),
    opts_help("Builds run in the background; the tab shows their state and log."),
    footer = shiny::tagList(
      shiny::actionButton("ref_build_go", "Start build", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Use-an-existing-database modal: manifest entries plus typed paths
#' @param manifest From [read_ref_manifest()]. @param ref_id Group key.
#' @noRd
reference_existing_modal <- function(manifest, ref_id) {
  listing <- if (nrow(manifest) == 0) {
    shiny::helpText("No databases recorded in this reference folder yet.")
  } else {
    shinyWidgets::prettyRadioButtons("ref_existing_pick", "Databases in the reference folder:",
      choices = stats::setNames(as.character(seq_len(nrow(manifest))),
        sprintf("%s [%s] - %s records, built %s", manifest$clade,
          manifest$rank %|NA|% "?", format(manifest$n_records, big.mark = ","),
          manifest$built)),
      shape = "round", status = "primary")
  }
  shiny::modalDialog(
    title = sprintf("Use an existing database for '%s'", ref_id), size = "l",
    listing,
    shiny::textInput("ref_existing_seeds", "Seed FASTA path:", width = "100%"),
    shiny::textInput("ref_existing_labels", "Label FASTA path:", width = "100%"),
    opts_help("Choosing a listed database fills these paths; you can also type ",
      "paths to a database built somewhere else."),
    footer = shiny::tagList(
      shiny::actionButton("ref_existing_apply", "Use this database", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-modal.R")'`
Expected: PASS, 6 tests.

- [ ] **Step 5: Add the handlers**

Append `reference_handlers()` to `R/app_reference.R`:

```r
#' Wire the Reference tab
#'
#' @param deps list(con = fn -> DBI con, project = fn -> chr, refresh = reactiveVal,
#'   selected_ids = fn -> chr, table_data = reactive -> data.frame).
#' @return list(open_ladder = function(ref_id), poll = function()).
#' @noRd
reference_handlers <- function(input, output, session, deps) {
  ladder_target <- shiny::reactiveVal(NULL)
  existing_target <- shiny::reactiveVal(NULL)
  build_targets <- shiny::reactiveVal(character(0))
  picked <- shiny::reactiveVal(NULL)   # taxid chosen from the search box

  with_con <- function(fn) {
    con <- deps$con(); on.exit(DBI::dbDisconnect(con)); fn(con)
  }
  bump <- function() deps$refresh(shiny::isolate(deps$refresh()) + 1)
  notify <- function(msg, type = "message") shiny::showNotification(msg, type = type, duration = 8)

  # ---- Suggest ----
  shiny::observeEvent(input$ref_suggest_btn, {
    shinyWidgets::sendSweetAlert(title = "Suggesting clades...",
      text = "Reading sample lineages and checking GenBank record counts. This can take a moment.",
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        res <- tryCatch(with_con(function(con) {
          cfg <- ref_config(con)
          l <- sample_ladders(con)
          if (nrow(l) == 0) stop("No samples have a resolved lineage yet. Assign taxonomy first.",
            call. = FALSE)
          cf <- cached_count_fn(cfg$ref_dir,
            api_key = Sys.getenv("NCBI_API_KEY", ""), conda_env = NULL)
          s <- suggest_reference_groups(l, cf, cfg$min_records, cfg$max_records)
          if (nrow(s$groups) == 0) stop("Could not count GenBank records for any clade (is Entrez Direct on PATH?).",
            call. = FALSE)
          write_reference_groups(con, s$groups, s$members)
          nrow(s$groups)
        }), error = function(e) e)
        if (inherits(res, "error")) {
          shinyWidgets::sendSweetAlert(title = "Suggestion failed",
            text = conditionMessage(res), type = "error")
        } else {
          bump()
          shinyWidgets::sendSweetAlert(title = "Clades suggested",
            text = sprintf("Proposed %d reference group(s). Review, then build or skip.", res),
            type = "success")
        }
      })
    }, delay = 0.1)
  })

  # ---- ladder modal ----
  open_ladder <- function(ref_id) {
    if (length(ref_id) == 0) { notify("Select a group first.", "warning"); return() }
    ref_id <- ref_id[1]
    picked(NULL)
    ladder_target(ref_id)
    with_con(function(con) {
      row <- DBI::dbGetQuery(con, "SELECT * FROM reference WHERE ref_id = ?", params = list(ref_id))
      if (nrow(row) == 0) { notify("Group not found.", "warning"); return() }
      ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE ref_id = ?", params = list(ref_id))$ID
      cfg <- ref_config(con)
      l <- group_ladder(sample_ladders(con), ids,
        count_fn = cached_count_fn(cfg$ref_dir, api_key = Sys.getenv("NCBI_API_KEY", "")))
      shiny::showModal(reference_ladder_modal(ref_id, row$clade_name[1], l,
        row$clade_taxid[1], length(ids)))
    })
  }
  shiny::observeEvent(input$ref_edit_btn, open_ladder(deps$selected_ids()))
  shiny::observeEvent(input$ref_open_clade, {
    i <- suppressWarnings(as.integer(input$ref_open_clade))
    d <- deps$table_data()
    if (!is.na(i) && i >= 1 && i <= nrow(d)) open_ladder(d$ID[i])
  })

  # taxdump search inside the ladder modal
  q <- shiny::debounce(shiny::reactive(input$ref_search), 300)
  output$ref_search_results <- shiny::renderUI({
    query <- trimws(q() %||% "")
    if (nchar(query) < 2) return(NULL)
    r <- tax_search(project_taxdump_dir(deps$project()), query, 20)
    if (nrow(r) == 0) return(shiny::helpText("No matches in the local taxdump."))
    shinyWidgets::prettyRadioButtons("ref_search_pick", NULL,
      choices = stats::setNames(as.character(r$taxid),
        sprintf("%s   [%s]   (taxid %s)", r$name, ifelse(is.na(r$rank), "?", r$rank), r$taxid)),
      shape = "round", status = "primary")
  })

  shiny::observeEvent(input$ref_clade_apply, {
    ref_id <- ladder_target()
    tid <- suppressWarnings(as.integer(input$ref_search_pick %||% input$ref_rung %||% NA))
    if (is.na(tid)) { notify("Pick a clade first.", "warning"); return() }
    res <- tryCatch(with_con(function(con) {
      ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE ref_id = ?",
        params = list(ref_id))$ID
      ladders <- sample_ladders(con)
      covered <- intersect(ladder_coverage(ladders, tid), ids)
      info <- ladders[ladders$taxid == tid, , drop = FALSE]
      if (nrow(info) == 0) {
        info <- with_taxdump_data(project_taxdump_dir(deps$project()), function(dd) {
          i <- tk_taxid_info(tid, dd)
          data.frame(name = i$sci_name[1], rank = i$rank[1], stringsAsFactors = FALSE)
        })
      }
      cfg <- ref_config(con)
      n <- cached_count_fn(cfg$ref_dir, api_key = Sys.getenv("NCBI_API_KEY", ""))(tid, info$name[1])
      set_reference_clade(con, ref_id = ref_id, clade_taxid = tid, clade_name = info$name[1],
        clade_rank = info$rank[1], n_records = n,
        member_ids = if (length(covered) > 0) covered else ids)
      list(name = info$name[1], left = setdiff(ids, covered))
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal()
    bump()
    if (length(res$left) > 0) {
      notify(sprintf("%d sample(s) are not inside %s and were left in their own group; open it to pick a clade.",
        length(res$left), res$name), "warning")
    } else {
      notify(sprintf("Group now uses %s.", res$name))
    }
  })

  # ---- skip ----
  shiny::observeEvent(input$ref_skip_btn, {
    ids <- deps$selected_ids()
    if (length(ids) == 0) { notify("Select a group first.", "warning"); return() }
    with_con(function(con) for (r in ids) skip_reference_group(con, r))
    bump()
    notify(sprintf("Skipped %d group(s); their samples keep their current database.", length(ids)))
  })

  # ---- use existing ----
  shiny::observeEvent(input$ref_existing_btn, {
    ids <- deps$selected_ids()
    if (length(ids) != 1) { notify("Select exactly one group.", "warning"); return() }
    existing_target(ids[1])
    with_con(function(con) {
      shiny::showModal(reference_existing_modal(read_ref_manifest(ref_config(con)$ref_dir), ids[1]))
    })
  })
  shiny::observeEvent(input$ref_existing_pick, {
    i <- suppressWarnings(as.integer(input$ref_existing_pick))
    with_con(function(con) {
      m <- read_ref_manifest(ref_config(con)$ref_dir)
      if (!is.na(i) && i >= 1 && i <= nrow(m)) {
        shiny::updateTextInput(session, "ref_existing_seeds", value = m$seeds[i])
        shiny::updateTextInput(session, "ref_existing_labels", value = m$labels[i])
      }
    })
  })
  shiny::observeEvent(input$ref_existing_apply, {
    res <- tryCatch(with_con(function(con) use_existing_reference(con, existing_target(),
      input$ref_existing_seeds %||% "", input$ref_existing_labels %||% "")),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    notify(sprintf("Pointed %d sample(s) at the existing database.", res))
  })

  # ---- build ----
  shiny::observeEvent(input$ref_build_btn, {
    ids <- deps$selected_ids()
    if (length(ids) == 0) { notify("Select the group(s) to build.", "warning"); return() }
    build_targets(ids)
    with_con(function(con) shiny::showModal(reference_build_modal(ids, ref_config(con)$ref_dir)))
  })
  shiny::observeEvent(input$ref_build_go, {
    ids <- build_targets()
    res <- tryCatch(with_con(function(con) {
      for (r in ids) {
        start_reference_build(con, deps$project(), r,
          api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
          conda_env = input$ref_conda %||% "", overwrite = isTRUE(input$ref_overwrite))
      }
      length(ids)
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    notify(sprintf("Started %d build(s); watch the State column.", res))
  })

  # ---- log tail ----
  shiny::observeEvent(input$ref_open_log, {
    i <- suppressWarnings(as.integer(input$ref_open_log))
    d <- deps$table_data()
    if (is.na(i) || i < 1 || i > nrow(d)) return()
    ref_id <- d$ID[i]
    row <- with_con(function(con) DBI::dbGetQuery(con,
      "SELECT log_file, fail_reason FROM reference WHERE ref_id = ?", params = list(ref_id)))
    log <- row$log_file[1]
    txt <- if (is.na(log) || !file.exists(log)) "No log for this group yet." else
      paste(utils::tail(readLines(log, warn = FALSE), 200), collapse = "\n")
    shiny::showModal(shiny::modalDialog(title = paste("Build log:", ref_id), size = "l",
      easyClose = TRUE,
      shiny::tags$pre(style = "max-height:60vh; overflow:auto; font-size:12px;", txt),
      footer = shiny::modalButton("Close")))
  })

  # ---- reference folder ----
  shiny::observeEvent(input$ref_dir_btn, {
    with_con(function(con) {
      shiny::showModal(shiny::modalDialog(title = "Reference database folder", size = "l",
        shiny::textInput("ref_dir_path", NULL, width = "100%", value = ref_config(con)$ref_dir),
        opts_help("Databases built here are reused by every project pointed at this ",
          "folder, including a shared drive. Changing it does not move databases ",
          "already built, and groups already wired keep working."),
        shiny::numericInput("ref_min", "Minimum GenBank records for a suggested clade:",
          value = ref_config(con)$min_records, min = 1, step = 100),
        shiny::numericInput("ref_max", "Maximum GenBank records for a suggested clade:",
          value = ref_config(con)$max_records, min = 1, step = 1000),
        footer = shiny::tagList(
          shiny::actionButton("ref_dir_apply", "Save", class = "btn-primary btn-sm"),
          shiny::modalButton("Cancel"))))
    })
  })
  shiny::observeEvent(input$ref_dir_apply, {
    v <- validate_ref_dir(input$ref_dir_path %||% "")
    if (!isTRUE(v$ok)) { notify(v$msg, "error"); return() }
    with_con(function(con) set_ref_config(con, ref_dir = v$msg,
      min_records = input$ref_min, max_records = input$ref_max))
    shiny::removeModal(); bump()
    notify(paste("Reference folder set to", v$msg))
  })

  list(open_ladder = open_ladder,
    poll = function() tryCatch(with_con(poll_reference_builds), error = function(e) invisible(0L)))
}
```

- [ ] **Step 6: Run the whole suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add R/app_reference.R tests/testthat/test-reference-modal.R
git commit -m "feat(reference): Reference tab table, modals, and handlers"
```

---

### Task 7: Wire the tab into the app

**Files:**
- Modify: `R/app_ui.R` (mode picker entry, control bar, folder header)
- Modify: `R/app_server.R` (mode metadata, table data, column definitions, handlers, build polling, run-modal warning)
- Test: `tests/testthat/test-reference-app.R`

**Interfaces:**
- Consumes: `reference_table()`, `reference_handlers()` from Task 6; `reference_pending()` from Task 5.
- Produces: no new exported names. The Reference mode string is `"Reference"`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-reference-app.R`:

```r
# The Reference mode is reachable from the UI and its rows carry the state icon.

test_that("the mode picker offers Reference between Taxonomy and Assemble", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  expect_true(grepl("Reference", html, fixed = TRUE))
  expect_lt(regexpr("Reference", html, fixed = TRUE),
    regexpr("Annotation Review", html, fixed = TRUE))
})

test_that("the Reference control bar exposes every action", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  for (id in c("ref_suggest_btn", "ref_edit_btn", "ref_existing_btn", "ref_skip_btn",
               "ref_build_btn", "ref_dir_btn")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
})

test_that("skipped rows get their own state icon", {
  expect_true("5" %in% names(RiboPilot:::REFERENCE_STATE_ICONS))
})
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-app.R")'`
Expected: FAIL, "Reference" is not in the UI.

- [ ] **Step 3: Add the mode and its control bar**

In `R/app_ui.R`, change the mode picker choices to:

```r
            choices = c(Taxonomy = "Taxonomy", Reference = "Reference", Assemble = "Assemble",
              "Annotation Review" = "Annotate", Export = "Export")),
```

and add a conditional panel after the Taxonomy one:

```r
          shiny::conditionalPanel("input.mode == 'Reference'",
            ctrls(bttn("ref_suggest_btn", "Suggest", "wand-magic-sparkles"),
              bttn("ref_edit_btn", "Edit Clade"), bttn("ref_existing_btn", "Use Existing"),
              bttn("ref_skip_btn", "Skip"), bttn("ref_build_btn", "Build", "circle-play"),
              bttn("ref_dir_btn", "Folder", "folder-open"))),
```

- [ ] **Step 4: Wire the server**

In `R/app_server.R`:

1. Add the skipped-state icon map near the top of the file (module level, outside `app_server`), so the test can see it:

```r
#' State icons for the Reference table (adds "skipped" to the shared vocabulary)
#' @noRd
REFERENCE_STATE_ICONS <- c(`0` = "fa fa-hourglass", `1` = "fa fa-person-running",
  `2` = "fa fa-circle-check", `3` = "fa fa-triangle-exclamation",
  `4` = "fa fa-circle-half-stroke", `5` = "fa fa-circle-minus")
```

2. In `mode_meta()`, add the Reference branch:

```r
      Reference = list(tbl = "reference", pfx = "reference"),
```

3. In `samples_data()`, add the Reference branch before the Assemble one:

```r
    } else if (mode == "Reference") {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      reference_table(con)
```

4. In the `output$tbl` column definitions, add Reference-specific columns next to the existing per-mode column definitions, and use the extended icon map when the mode is Reference:

```r
    if (mode == "Reference") {
      cols[["clade_name"]] <- reactable::colDef(name = "Clade", minWidth = 140, html = TRUE,
        cell = rt_link("ref_open_clade"))
      cols[["clade_rank"]] <- reactable::colDef(name = "Rank", width = 90)
      cols[["n_samples"]] <- reactable::colDef(name = "Samples", width = 90)
      cols[["n_records"]] <- reactable::colDef(name = "GenBank records", width = 140)
      cols[["source"]] <- reactable::colDef(name = "Source", width = 90)
      cols[["database"]] <- reactable::colDef(name = "Database", minWidth = 240, html = TRUE,
        cell = rt_longtext())
      cols[["fail_reason"]] <- reactable::colDef(name = "Note", minWidth = 200, html = TRUE,
        cell = rt_longtext())
      cols[[state_col]] <- reactable::colDef(name = "", sticky = "left", html = TRUE, width = 32,
        align = "center", filterable = FALSE, cell = rt_dynamicIcon(REFERENCE_STATE_ICONS))
    }
```

Place this after the block that sets the generic `state_col` column definition so it overrides it. Add `Reference = "No reference groups yet. Press Suggest to propose clades from the assigned taxonomy."` to the `nodata` switch.

5. Wire the handlers next to the existing `opts_handlers(...)` / `taxonomy_handlers(...)` calls:

```r
  ref_h <- reference_handlers(input, output, session, list(
    con = app_con, project = app_project, refresh = refresh,
    selected_ids = selected_ids, table_data = samples_data))
```

6. Poll running builds while the tab is open:

```r
  observe({
    if ((input$mode %||% "") != "Reference") return()
    invalidateLater(4000, session)
    if (isTRUE(ref_h$poll() > 0)) refresh(isolate(refresh()) + 1)
  })
```

7. Warn (do not block) when the workflow is launched with groups still pending. In `observeEvent(input$run_btn, ...)`, immediately after the Taxonomy early return:

```r
    pending <- tryCatch({
      c1 <- app_con(); on.exit(DBI::dbDisconnect(c1), add = TRUE); reference_pending(c1)
    }, error = function(e) character(0))
    if (length(pending) > 0) {
      showNotification(sprintf(
        "%d reference group(s) are not built or skipped (%s). Those samples will assemble with the database their options set already has.",
        length(pending), paste(pending, collapse = ", ")), type = "warning", duration = 12)
    }
```

8. In `observeEvent(input$open_detail, ...)`, add a Reference branch so clicking the key cell opens the ladder:

```r
    } else if (mode == "Reference") {
      ref_h$open_ladder(id)
```

- [ ] **Step 5: Run the test and confirm it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-app.R")'`
Expected: PASS, 3 tests.

- [ ] **Step 6: Check the app actually starts**

Run: `Rscript -e 'devtools::load_all("."); ui <- RiboPilot:::app_ui(NULL); stopifnot(inherits(ui, "shiny.tag.list")); cat("ui ok\n")'`
Expected: `ui ok`, no error.

Then open a real project and click through the tab manually:

```r
devtools::load_all(".")
migrate_project("<a test project>")
RiboPilot(<a test project>)
```

Confirm: the Reference tab lists nothing before Suggest, Suggest fills it, the clade cell opens the ladder, Skip marks a row, and Folder saves a new path.

- [ ] **Step 7: Run the whole suite**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: all tests pass.

- [ ] **Step 8: Commit**

```bash
git add R/app_ui.R R/app_server.R tests/testthat/test-reference-app.R
git commit -m "feat(app): add the Reference tab"
```

---

### Task 8: Documentation

**Files:**
- Modify: `vignettes/Reference-Databases.Rmd`
- Modify: `inst/ref_dbs/getOrganelle/GenBankDownload/README.md`
- Modify: `R/build_reference_db.R` (roxygen cross-reference only)

**Interfaces:**
- Consumes: the behavior built in Tasks 1 to 7.
- Produces: no code interfaces.

- [ ] **Step 1: Read the current vignette**

Run: `Rscript -e 'cat(readLines("vignettes/Reference-Databases.Rmd"), sep = "\n")'`
Note where it tells the reader to call `build_reference_db()` by hand.

- [ ] **Step 2: Rewrite the vignette's workflow section**

Lead with the in-app path and keep the console path as the alternative. The new section must state:

- The Reference tab sits between Taxonomy and Assemble and runs after taxonomy is assigned and locked.
- Suggest groups samples by walking each lineage up until the clade has at least 1,000 GenBank rRNA records, stopping before 25,000; a project with distant groups gets one database per group.
- The clade cell opens the shared lineage; each rung shows samples covered and record count; a rung covering part of the group offers to split the rest.
- Build runs in the background and writes to the reference folder; the folder is set with the Folder button or `new_project(ref_dir = )` and can be a shared drive.
- A finished build creates an `assemble_opts` preset named after the clade and moves that group's samples onto it; locked samples are left alone.
- Skip keeps whatever database the samples' preset already has; Use Existing points at a database built earlier or by someone else.
- The step is optional: Assemble is never blocked, it only warns.

- [ ] **Step 3: Update the builder README**

In `inst/ref_dbs/getOrganelle/GenBankDownload/README.md`, add a short paragraph above "### From R (recommended)" saying the app's Reference tab is now the usual way to build these databases and that the R and bash paths remain for scripted use.

- [ ] **Step 4: Cross-reference from the function docs**

In `R/build_reference_db.R`, extend the `build_reference_db()` description with one sentence: the RiboPilot app's Reference tab calls this in the background and wires the result into the project, so calling it directly is for scripted or headless use.

- [ ] **Step 5: Rebuild docs and check the package**

Run: `Rscript -e 'devtools::document()'`
Run: `Rscript -e 'devtools::check(document = FALSE, args = "--no-manual", error_on = "warning")'`
Expected: no ERRORs or WARNINGs. Notes about the vignette build needing external tools are acceptable if they predate this work.

- [ ] **Step 6: Commit**

```bash
git add vignettes/Reference-Databases.Rmd inst/ref_dbs/getOrganelle/GenBankDownload/README.md R/build_reference_db.R man
git commit -m "docs(reference): document the in-app reference database step"
```

---

## Self-Review Notes

Spec coverage check, section by section:

- Data model (`reference`, `reference_config`, `samples.ref_id`, manifest, `ref_dir` argument): Tasks 1 and 2. The spec calls the state column `state`; the implementation names it `reference_switch` so the app's existing per-mode `<pfx>_switch` machinery works without a special case. Same values, same meaning.
- Suggestion algorithm (ranked lineage, walk up from family, floor, ceiling, merge, cached counts, editable window): Tasks 3 and 4, window editable in Task 6's folder modal. The spec proposed a new `taxonkit reformat -t` helper; `samples.lineage_taxids` already carries exactly that, so Task 4 reads it instead of adding a call.
- Reference tab, ladder modal, partial-coverage split, off-ladder search, Use Existing, Skip, folder editing, no gating on Assemble: Tasks 6 and 7.
- Building, polling, conda / API key / RefSeq / rebuild options, named failures, wiring on success: Tasks 5 and 6.
- Testing (offline grouping tests, wiring against a temp SQLite project, manifest round trip, folder resolution order): Tasks 2 to 5.
