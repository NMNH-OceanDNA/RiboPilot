# Reference databases through Nextflow: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move GenBank record counting and clade-database building out of host
processes and into a Nextflow workflow, and give Build the same run modal the
Assemble tab uses.

**Architecture:** A new `inst/nextflow/main_reference.nf` with two entry points
(`COUNT`, `BUILD`), driven by a work-list TSV the app writes. R gains
`R/reference_nf.R` (manifest, command builder, blocking count runner, count
merging) and `R/nf_progress.R` (the ANSI board parser, extracted from
`app_server.R` so both run modals share it). The group state machine in
`reference_build.R` and the ladder walk in `reference_groups.R` are unchanged
apart from where their counts come from.

**Tech stack:** R (shiny, processx, DBI/RSQLite), testthat, Nextflow DSL2, bash.

Spec: `dev/design/16_REFERENCE_NEXTFLOW.md`.

## Global constraints

- ASCII only in R code, and in the Nextflow and bash sources.
- Minimal comments; comment the why, not the what.
- Every `download_rrna.sh` invocation runs with stdin detached (`stdin = "/dev/null"`
  from R, `< /dev/null` from a shell). Entrez Direct blocks forever on an open
  idle stdin.
- Suite green between tasks: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`.
- Commit after each task. Do not push.
- Branch: cut `reference-nextflow` from `refdb-picker` before Task 1.
- Count cache file stays `<ref_dir>/counts.tsv`, columns `taxid`, `n_records`,
  tab separated with a header. Existing caches must keep working.
- Seed/label output paths stay exactly what `reference_db_paths()` returns:
  `<ref_dir>/seeds/<ref_id>_rrna_seeds.fasta`, `<ref_dir>/labels/<ref_id>_rrna_labels.fasta`.

## File structure

| File | Responsibility |
| --- | --- |
| `R/reference_nf.R` (new) | Work-list TSV, Nextflow command builder, blocking count runner, count-file reading and merging |
| `R/nf_progress.R` (new) | ANSI progress-board parser and its render, shared by both run modals |
| `R/reference_counts.R` | Gains `read_counts_cache()` / `cache_count_fn()`; loses nothing |
| `R/reference_build.R` | `start_reference_build()` launches Nextflow instead of Rscript |
| `R/app_modals.R` | `reference_run_modal()` |
| `R/app_reference.R` | Suggest and apply-clade use the count runner; Build opens the run modal |
| `R/app_server.R` | Uses the extracted parser; mounts the reference run modal handlers |
| `inst/nextflow/main_reference.nf` (new) | COUNT and BUILD entry points |
| `inst/nextflow/bin/rp_ref_db.R` (new) | Resolves and runs `download_rrna.sh` inside a task |

---

### Task 1: Count cache reading, split out

**Files:**
- Modify: `R/reference_counts.R`
- Test: `tests/testthat/test-reference-counts.R`

**Interfaces:**
- Produces: `read_counts_cache(dir)` -> `data.frame(taxid = integer, n_records = integer)`,
  empty frame when the file is missing or malformed. `cache_count_fn(dir)` ->
  `function(taxid, name)` returning the cached count or `NA_integer_`, never
  touching the network.

- [ ] **Step 1: Write the failing tests**

```r
test_that("read_counts_cache returns an empty frame when there is no cache", {
  d <- withr::local_tempdir()
  x <- RiboPilot:::read_counts_cache(d)
  expect_equal(nrow(x), 0L)
  expect_equal(names(x), c("taxid", "n_records"))
})

test_that("cache_count_fn answers from the cache and never calls out", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = 7898L, n_records = 11750L),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  f <- RiboPilot:::cache_count_fn(d)
  expect_equal(f(7898L, "Actinopterygii"), 11750L)
  expect_true(is.na(f(6656L, "Arthropoda")))
})
```

- [ ] **Step 2: Run them and watch them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-counts.R")'`
Expected: FAIL, "could not find function read_counts_cache".

- [ ] **Step 3: Implement**

In `R/reference_counts.R`, lift the `read_cache()` closure out of
`cached_count_fn()` into a top-level function and have `cached_count_fn()` call
it, so there is one reader:

```r
#' Read a reference folder's record-count cache
#' @noRd
read_counts_cache <- function(dir) {
  empty <- data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE)
  p <- ref_counts_path(dir)
  if (!file.exists(p)) return(empty)
  d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !all(c("taxid", "n_records") %in% names(d))) return(empty)
  d$taxid <- suppressWarnings(as.integer(d$taxid))
  d$n_records <- suppressWarnings(as.integer(d$n_records))
  d[, c("taxid", "n_records"), drop = FALSE]
}

#' Cache-only count lookup: NA for anything not counted yet
#' @noRd
cache_count_fn <- function(dir) {
  function(taxid, name) {
    cache <- read_counts_cache(dir)
    i <- match(suppressWarnings(as.integer(taxid)), cache$taxid)
    if (is.na(i)) NA_integer_ else cache$n_records[i]
  }
}
```

- [ ] **Step 4: Run the file, then the suite**

Both green.

- [ ] **Step 5: Commit**

```bash
git add R/reference_counts.R tests/testthat/test-reference-counts.R
git commit -m "refactor(reference): split out the count-cache reader"
```

---

### Task 2: Work list and count merging

**Files:**
- Create: `R/reference_nf.R`
- Test: `tests/testthat/test-reference-nf.R` (new)

**Interfaces:**
- Consumes: `read_counts_cache()` from Task 1.
- Produces:
  - `write_clade_manifest(path, clades)` -> `path`. `clades` is a data.frame with
    `ref_id`, `taxid`, `clade_name`, `query`; writes those four columns as a
    tab-separated file with a header.
  - `read_count_files(dir)` -> `data.frame(taxid, n_records)` from every
    `<taxid>.count` file in `dir` (one line, `taxid<TAB>n_records`). Unreadable
    or empty files are skipped.
  - `merge_counts(ref_dir, counts)` -> the merged `data.frame`, also written to
    `<ref_dir>/counts.tsv`. New values replace existing rows for the same taxid.

- [ ] **Step 1: Write the failing tests**

```r
test_that("write_clade_manifest writes the four columns Nextflow reads", {
  d <- withr::local_tempdir()
  p <- RiboPilot:::write_clade_manifest(file.path(d, "clades.tsv"),
    data.frame(ref_id = "actinopteri", taxid = 186623L, clade_name = "Actinopteri",
      query = '"Actinopteri"[Organism]', stringsAsFactors = FALSE))
  x <- utils::read.delim(p, stringsAsFactors = FALSE)
  expect_equal(names(x), c("ref_id", "taxid", "clade_name", "query"))
  expect_equal(x$taxid, 186623L)
})

test_that("read_count_files reads every count file and skips junk", {
  d <- withr::local_tempdir()
  writeLines("7898\t11750", file.path(d, "7898.count"))
  writeLines("186623\t9000", file.path(d, "186623.count"))
  writeLines("", file.path(d, "6656.count"))
  x <- RiboPilot:::read_count_files(d)
  expect_equal(nrow(x), 2L)
  expect_equal(x$n_records[x$taxid == 7898L], 11750L)
})

test_that("merge_counts replaces stale rows and keeps the rest", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = c(7898L, 33208L), n_records = c(10L, 500L)),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out <- RiboPilot:::merge_counts(d, data.frame(taxid = 7898L, n_records = 11750L))
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 500L)
  expect_equal(nrow(RiboPilot:::read_counts_cache(d)), 2L)
})
```

- [ ] **Step 2: Run and watch them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-nf.R")'`

- [ ] **Step 3: Implement**

```r
# Launching the reference workflow, and getting its results back into R.
#
# Counts come back as one small file per clade rather than appended rows: the
# tasks run in parallel and would otherwise race on a single cache file.

#' Write the clade work list the reference workflow reads
#' @noRd
write_clade_manifest <- function(path, clades) {
  cols <- c("ref_id", "taxid", "clade_name", "query")
  d <- clades[, cols, drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(d, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  path
}

#' Read the per-clade count files a COUNT run produced
#' @noRd
read_count_files <- function(dir) {
  empty <- data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE)
  files <- list.files(dir, pattern = "\\.count$", full.names = TRUE)
  if (length(files) == 0) return(empty)
  rows <- lapply(files, function(f) {
    x <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    x <- trimws(x[nzchar(trimws(x))])
    if (length(x) == 0) return(NULL)
    p <- strsplit(x[1], "\t")[[1]]
    if (length(p) < 2) return(NULL)
    tx <- suppressWarnings(as.integer(p[1])); n <- suppressWarnings(as.integer(p[2]))
    if (is.na(tx) || is.na(n)) return(NULL)
    data.frame(taxid = tx, n_records = n, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Merge fresh counts into a reference folder's cache
#' @noRd
merge_counts <- function(ref_dir, counts) {
  cache <- read_counts_cache(ref_dir)
  if (nrow(counts) > 0) {
    cache <- cache[!(cache$taxid %in% counts$taxid), , drop = FALSE]
    cache <- rbind(cache, counts[, c("taxid", "n_records"), drop = FALSE])
  }
  dir.create(path.expand(ref_dir), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(cache, ref_counts_path(ref_dir), sep = "\t",
    row.names = FALSE, quote = FALSE, na = "")
  cache
}
```

- [ ] **Step 4: Run the file, then the suite**

- [ ] **Step 5: Commit**

```bash
git add R/reference_nf.R tests/testthat/test-reference-nf.R
git commit -m "feat(reference): clade work list and count merging"
```

---

### Task 3: The Nextflow command builder

**Files:**
- Modify: `R/reference_nf.R`
- Test: `tests/testthat/test-reference-nf.R`

**Interfaces:**
- Produces: `nf_reference_command(project, entry, manifest, ref_dir, out_dir,
  api_key = "", refseq = FALSE, resume = TRUE, work = file.path(project, "work_ref"))`
  -> character vector of arguments following the `nextflow` binary. `entry` is
  `"COUNT"` or `"BUILD"`.
  `ref_bind_option(cfg_lines, ref_dir)` -> character(0), or the two-element
  `c("-process.containerOptions", "<value>")` when `.config` enables docker
  (`-v dir:dir`) or singularity/apptainer (`--bind dir`).

- [ ] **Step 1: Write the failing tests**

```r
test_that("nf_reference_command names the entry point and every param", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  a <- RiboPilot:::nf_reference_command(d, "COUNT", file.path(d, "clades.tsv"),
    "/refs", file.path(d, "out"))
  expect_equal(a[1], "run")
  expect_true("-entry" %in% a)
  expect_equal(a[which(a == "-entry") + 1], "COUNT")
  expect_equal(a[which(a == "--refDir") + 1], "/refs")
  expect_true("--clades" %in% a && "--outDir" %in% a)
})

test_that("the reference folder is bind-mounted for docker and for singularity", {
  docker <- c("docker {", "  enabled = true", "}")
  sing <- c("singularity {", "  enabled = true", "}")
  expect_equal(RiboPilot:::ref_bind_option(docker, "/refs"),
    c("-process.containerOptions", "-v /refs:/refs"))
  expect_equal(RiboPilot:::ref_bind_option(sing, "/refs"),
    c("-process.containerOptions", "--bind /refs"))
  expect_length(RiboPilot:::ref_bind_option("process { executor = 'local' }", "/refs"), 0L)
})

test_that("refseq and the api key only appear when asked for", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  plain <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o")
  expect_false("--refseq" %in% plain)
  expect_false("--ncbi_api_key" %in% plain)
  keyed <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o",
    api_key = "abc", refseq = TRUE)
  expect_equal(keyed[which(keyed == "--ncbi_api_key") + 1], "abc")
  expect_equal(keyed[which(keyed == "--refseq") + 1], "true")
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Implement**

```r
#' Container bind option for the reference folder, if the config uses a container
#'
#' The folder lives outside the project, so a container task cannot see it. The
#' mount is added per run rather than written into the project config, so a user
#' who moves the folder needs no project edits.
#' @noRd
ref_bind_option <- function(cfg_lines, ref_dir) {
  txt <- paste(cfg_lines, collapse = "\n")
  d <- path.expand(ref_dir)
  if (grepl("docker\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions", paste0("-v ", d, ":", d)))
  }
  if (grepl("(singularity|apptainer)\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions", paste0("--bind ", d)))
  }
  character(0)
}

#' Nextflow arguments for the reference workflow
#' @noRd
nf_reference_command <- function(project, entry, manifest, ref_dir, out_dir,
                                 api_key = "", refseq = FALSE, resume = TRUE,
                                 work = file.path(project, "work_ref")) {
  entry <- match.arg(entry, c("COUNT", "BUILD"))
  main_nf <- app_sys("nextflow", "main_reference.nf")
  cfg <- file.path(project, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg) else character(0)
  c("run", main_nf, "-entry", entry,
    if (file.exists(cfg)) c("-c", cfg),
    "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    "--clades", path.expand(manifest),
    "--refDir", path.expand(ref_dir),
    "--outDir", path.expand(out_dir),
    "--ribopilot_src", ribopilot_src_default(),
    if (isTRUE(refseq)) c("--refseq", "true"),
    if (nzchar(api_key %||% "")) c("--ncbi_api_key", api_key),
    ref_bind_option(cfg_lines, ref_dir))
}
```

- [ ] **Step 4: Run the file, then the suite**

- [ ] **Step 5: Commit**

```bash
git add R/reference_nf.R tests/testthat/test-reference-nf.R
git commit -m "feat(reference): nextflow command builder for the reference workflow"
```

---

### Task 4: The workflow and its task wrapper

**Files:**
- Create: `inst/nextflow/main_reference.nf`, `inst/nextflow/bin/rp_ref_db.R`
- Test: `tests/testthat/test-reference-workflow.R` (new)

**Interfaces:**
- Produces: entry points `COUNT` and `BUILD`; `rp_ref_db.R` usable as
  `rp_ref_db.R count <taxid> <query> <outfile>` and
  `rp_ref_db.R build <prefix> <query> <outdir>`.

- [ ] **Step 1: Write the failing tests**

These pin the contract without running Nextflow (no cluster, no network in CI):

```r
test_that("the reference workflow ships with both entry points", {
  f <- system.file("nextflow", "main_reference.nf", package = "RiboPilot")
  expect_true(nzchar(f))
  txt <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("workflow COUNT", txt, fixed = TRUE))
  expect_true(grepl("workflow BUILD", txt, fixed = TRUE))
  # the work list is read as a header TSV, matching write_clade_manifest()
  expect_true(grepl("splitCsv(header: true, sep: '\\t')", txt, fixed = TRUE))
})

test_that("the task wrapper detaches stdin from Entrez Direct", {
  f <- system.file("nextflow", "bin", "rp_ref_db.R", package = "RiboPilot")
  expect_true(nzchar(f))
  txt <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("/dev/null", txt, fixed = TRUE))
})
```

- [ ] **Step 2: Run and watch them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reference-workflow.R")'`

- [ ] **Step 3: Write `inst/nextflow/bin/rp_ref_db.R`**

```r
#!/usr/bin/env Rscript
# Run the bundled GenBank rRNA script for one clade.
# Args: count <taxid> <query> <outfile> | build <prefix> <query> <outdir>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))

script <- app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
if (!nzchar(script) || !file.exists(script)) stop("download_rrna.sh not found")
api_key <- Sys.getenv("NCBI_API_KEY", "")
refseq <- identical(Sys.getenv("RP_REFSEQ", ""), "true")

args <- c(shQuote(script))
if (refseq) args <- c(args, "--refseq")
if (nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))

if (identical(a[1], "count")) {
  args <- c(args, "--count-only", shQuote(a[3]), "count_probe")
  out <- suppressWarnings(system2("bash", args, stdout = TRUE, stderr = TRUE,
    stdin = "/dev/null"))
  n <- parse_record_count(out)
  if (is.na(n)) { cat(paste(out, collapse = "\n"), "\n"); stop("count failed") }
  writeLines(paste(a[2], n, sep = "\t"), a[4])
} else {
  args <- c(args, "--outdir", shQuote(a[4]), shQuote(a[3]), shQuote(a[2]))
  status <- system2("bash", args, stdin = "/dev/null")
  if (status != 0) stop("build failed for ", a[2])
}
```

- [ ] **Step 4: Write `inst/nextflow/main_reference.nf`**

```groovy
// Reference clade databases: count GenBank records, or build GetOrganelle seed
// and label FASTAs. Driven by a work-list TSV written by the app
// (write_clade_manifest): ref_id, taxid, clade_name, query.
//
// Two entry points rather than two files: the processes share their inputs,
// their container, and their retry policy, and a second file drifts out of step
// (see dev/design/15_WF2_DECOMMISSION.md).

params.clades = ''
params.refDir = ''
params.outDir = ''
params.ribopilot_src = ''
params.ncbi_api_key = ''
params.refseq = false

process count_clade {
    cache false
    tag "${clade_name}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(taxid), val(clade_name), val(query)

    output:
    path "${taxid}.count"

    script:
    """
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 20 )); fi
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'
    export RP_REFSEQ='${params.refseq}'
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_db.R count '${taxid}' '${query}' '${taxid}.count'
    """
}

process build_ref_db {
    cache false
    tag "${ref_id}"
    errorStrategy 'ignore'
    publishDir "${params.refDir}/seeds", mode: 'copy', pattern: '*_seeds.fasta'
    publishDir "${params.refDir}/labels", mode: 'copy', pattern: '*_labels.fasta'

    input:
    tuple val(ref_id), val(query)

    output:
    path "*_seeds.fasta", optional: true
    path "*_labels.fasta", optional: true

    script:
    """
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'
    export RP_REFSEQ='${params.refseq}'
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_db.R build '${ref_id}_rrna' '${query}' .
    """
}

workflow COUNT {
    Channel.fromPath(params.clades)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.taxid, r.clade_name, r.query) }
        | count_clade
}

workflow BUILD {
    Channel.fromPath(params.clades)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.ref_id, r.query) }
        | build_ref_db
}
```

- [ ] **Step 5: Run the tests, then the suite**

Both green. `parse_record_count()` already exists in `R/reference_counts.R`;
confirm it is reachable from the wrapper (it is package-internal, and `rp_env.R`
loads the package).

- [ ] **Step 6: Commit**

```bash
git add inst/nextflow/main_reference.nf inst/nextflow/bin/rp_ref_db.R tests/testthat/test-reference-workflow.R
git commit -m "feat(reference): count and build entry points for the reference workflow"
```

---

### Task 5: The blocking count runner

**Files:**
- Modify: `R/reference_nf.R`
- Test: `tests/testthat/test-reference-nf.R`

**Interfaces:**
- Produces: `run_clade_counts(project, ref_dir, clades, api_key = "", refseq = FALSE,
  runner = nf_run_blocking)` -> the merged count `data.frame`. Writes the work
  list and the run output under `<project>/.reference/`, runs COUNT, merges the
  results into the cache, and returns the merged cache. Clades already in the
  cache are dropped before launching; with nothing left to count it returns the
  cache without running anything.
  `nf_run_blocking(project, args)` -> exit status integer.

- [ ] **Step 1: Write the failing tests**

```r
test_that("run_clade_counts skips cached clades and merges what the run produced", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  utils::write.table(data.frame(taxid = 33208L, n_records = 900000L),
    file.path(ref, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  seen <- new.env()
  runner <- function(project, args) {
    seen$args <- args
    out <- args[which(args == "--outDir") + 1]
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines("7898\t11750", file.path(out, "7898.count"))
    0L
  }
  clades <- data.frame(ref_id = c("metazoa", "actinopterygii"), taxid = c(33208L, 7898L),
    clade_name = c("Metazoa", "Actinopterygii"),
    query = c('"Metazoa"[Organism]', '"Actinopterygii"[Organism]'),
    stringsAsFactors = FALSE)
  out <- RiboPilot:::run_clade_counts(d, ref, clades, runner = runner)
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 900000L)
  # only the uncached clade was sent to the run
  man <- utils::read.delim(seen$args[which(seen$args == "--clades") + 1], stringsAsFactors = FALSE)
  expect_equal(man$taxid, 7898L)
})

test_that("run_clade_counts launches nothing when every clade is cached", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  utils::write.table(data.frame(taxid = 7898L, n_records = 11750L),
    file.path(ref, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  runner <- function(project, args) stop("must not launch")
  out <- RiboPilot:::run_clade_counts(d, ref,
    data.frame(ref_id = "a", taxid = 7898L, clade_name = "Actinopterygii",
      query = "q", stringsAsFactors = FALSE), runner = runner)
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
})

test_that("a failed run still merges whatever counts landed", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  runner <- function(project, args) {
    out <- args[which(args == "--outDir") + 1]
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines("7898\t11750", file.path(out, "7898.count"))
    1L
  }
  out <- RiboPilot:::run_clade_counts(d, ref,
    data.frame(ref_id = c("a", "b"), taxid = c(7898L, 6656L),
      clade_name = c("Actinopterygii", "Arthropoda"), query = c("q1", "q2"),
      stringsAsFactors = FALSE), runner = runner)
  expect_equal(nrow(out), 1L)
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Implement**

```r
#' Run a Nextflow command to completion, returning its exit status
#'
#' Split out so tests can substitute a runner that spawns nothing.
#' @noRd
nf_run_blocking <- function(project, args) {
  if (!nzchar(Sys.which("nextflow"))) {
    stop("Nextflow not found on PATH.", call. = FALSE)
  }
  p <- processx::run("nextflow", args, wd = project, error_on_status = FALSE,
    stdout = "|", stderr = "2>&1")
  p$status
}

#' Count GenBank records for a set of clades, waiting for the run
#'
#' Only uncached clades are sent to Nextflow; the walk reads the cache.
#' @noRd
run_clade_counts <- function(project, ref_dir, clades, api_key = "", refseq = FALSE,
                             runner = nf_run_blocking) {
  cache <- read_counts_cache(ref_dir)
  need <- clades[!(as.integer(clades$taxid) %in% cache$taxid), , drop = FALSE]
  need <- need[!duplicated(need$taxid), , drop = FALSE]
  if (nrow(need) == 0) return(cache)
  rundir <- file.path(project, ".reference")
  out <- file.path(rundir, "counts")
  unlink(out, recursive = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  man <- write_clade_manifest(file.path(rundir, "clades_count.tsv"), need)
  args <- nf_reference_command(project, "COUNT", man, ref_dir, out,
    api_key = api_key, refseq = refseq, resume = FALSE)
  # a non-zero status is not fatal: clades that did report are still worth
  # keeping, and the walk treats a missing count as unknown
  runner(project, args)
  merge_counts(ref_dir, read_count_files(out))
}
```

- [ ] **Step 4: Run the file, then the suite**

- [ ] **Step 5: Commit**

```bash
git add R/reference_nf.R tests/testthat/test-reference-nf.R
git commit -m "feat(reference): blocking clade-count runner"
```

---

### Task 6: Suggest and apply-clade count through Nextflow

**Files:**
- Modify: `R/app_reference.R`, `R/reference_groups.R`
- Test: `tests/testthat/test-reference-groups.R`, `tests/testthat/test-reference-modal.R`

**Interfaces:**
- Produces: `candidate_rungs(ladders, start_rank = "family")` ->
  `data.frame(ref_id, taxid, clade_name, query)`, deduplicated by taxid, holding
  every rung at or above `start_rank`. `ref_id` is `ref_id_from_name(clade_name)`,
  `query` is `clade_query(clade_name)`.
- Consumes: `run_clade_counts()`, `cache_count_fn()`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("candidate_rungs covers every rung at or above the start rank, once each", {
  con <- mk_ladder_db(); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  r <- RiboPilot:::candidate_rungs(l, start_rank = "family")
  expect_true(all(c("ref_id", "taxid", "clade_name", "query") %in% names(r)))
  expect_false(any(duplicated(r$taxid)))
  expect_true(all(r$taxid %in% l$taxid[l$rank %in%
    c("kingdom", "phylum", "class", "order", "family")]))
  expect_false(any(r$rank %in% "species"))
  expect_true(grepl("[Organism]", r$query[1], fixed = TRUE))
})
```

(`mk_ladder_db()` already exists at the top of `test-reference-groups.R`; reuse
it. If its samples stop at family, assert against the ranks it actually has.)

In `test-reference-modal.R`, add a Suggest test driven through
`reference_handlers()` with a stub runner, mirroring the existing apply tests:

```r
test_that("Suggest counts only uncached clades and then walks offline", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)          # its counts.tsv already covers 33208/7586/7588
  deps <- mk_split_deps(dbfile, dir)
  calls <- 0L
  withr::local_options(list(ribopilot.count_runner = function(project, args) {
    calls <<- calls + 1L; 0L
  }))
  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_suggest_btn = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  expect_gt(nrow(DBI::dbGetQuery(con, "SELECT * FROM reference")), 0L)
  expect_equal(calls, 0L)   # every rung was already cached
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Implement `candidate_rungs()` in `R/reference_groups.R`**

```r
#' Every clade the walk could ask about, once each
#'
#' The walk starts at `start_rank` and climbs, so anything deeper is never
#' counted. Counting the whole set up front is what lets the walk itself run
#' offline against the cache.
#' @noRd
candidate_rungs <- function(ladders, start_rank = "family") {
  top <- match(start_rank, RRNA_LINEAGE_RANKS)
  l <- ladders[ladders$level <= top, , drop = FALSE]
  l <- l[!duplicated(l$taxid), , drop = FALSE]
  if (nrow(l) == 0) {
    return(data.frame(ref_id = character(0), taxid = integer(0),
      clade_name = character(0), query = character(0), stringsAsFactors = FALSE))
  }
  data.frame(ref_id = vapply(l$name, ref_id_from_name, character(1), USE.NAMES = FALSE),
    taxid = as.integer(l$taxid), clade_name = l$name,
    query = vapply(l$name, clade_query, character(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE)
}
```

- [ ] **Step 4: Rewire Suggest in `R/app_reference.R`**

Inside the Suggest handler, replace the `cached_count_fn(...)` call with a
count run followed by a cache-only walk. The runner comes from an option so
tests can stub it:

```r
count_runner <- function() getOption("ribopilot.count_runner", nf_run_blocking)

# ... inside the Suggest tryCatch, after `l <- ladders_now(con)`:
run_clade_counts(deps$project(), cfg$ref_dir, candidate_rungs(l),
  api_key = Sys.getenv("NCBI_API_KEY", ""), runner = count_runner())
s <- suggest_reference_groups(l, cache_count_fn(cfg$ref_dir),
  cfg$min_records, cfg$max_records)
```

Keep the existing "could not count records for any clade" error: it now fires
when the walk produced no groups because every count was `NA`.

- [ ] **Step 5: Rewire apply-clade in `R/app_reference.R`**

In `apply_clade()`, replace the `cached_count_fn(...)(tid, info$name[1])` call:

```r
run_clade_counts(deps$project(), cfg$ref_dir,
  data.frame(ref_id = ref_id_from_name(info$name[1]), taxid = as.integer(tid),
    clade_name = info$name[1], query = clade_query(info$name[1]),
    stringsAsFactors = FALSE),
  api_key = Sys.getenv("NCBI_API_KEY", ""), runner = count_runner())
n <- cache_count_fn(cfg$ref_dir)(tid, info$name[1])
```

- [ ] **Step 6: Run the suite**

Green. The pre-seeded `counts.tsv` in `mk_split_db()` is what keeps these tests
off the network.

- [ ] **Step 7: Commit**

```bash
git add R/app_reference.R R/reference_groups.R tests/testthat/test-reference-groups.R tests/testthat/test-reference-modal.R
git commit -m "feat(reference): count clades through nextflow, walk offline"
```

---

### Task 7: Extract the Nextflow progress parser

**Files:**
- Create: `R/nf_progress.R`
- Modify: `R/app_server.R:723-795` (parser, `progress_update()`, `output$nf_prog_all`)
- Test: `tests/testthat/test-nf-progress.R` (new)

**Interfaces:**
- Produces:
  - `nf_progress_empty()` -> `list(header = character(0), executor = "", proc = list(), footer = character(0))`
  - `nf_progress_parse(state, lines)` -> updated state. Process rows are keyed by
    process name and overwritten in place; the newest executor line wins; lines
    before the first process row are header, later ones are footer.
  - `nf_progress_text(state, footer_n = 40)` -> one string, empty sections omitted.

- [ ] **Step 1: Write the failing tests**

```r
test_that("a redrawn process row overwrites rather than repeats", {
  s <- RiboPilot:::nf_progress_empty()
  s <- RiboPilot:::nf_progress_parse(s, c("N E X T F L O W  ~  version 25.10.6",
    "executor >  local (1)", "[a1/b2c3] preprocess (S1) | 0 of 1"))
  s <- RiboPilot:::nf_progress_parse(s, "[a1/b2c3] preprocess (S1) | 1 of 1 ✔")
  txt <- RiboPilot:::nf_progress_text(s)
  expect_equal(length(gregexpr("preprocess", txt)[[1]]), 1L)
  expect_true(grepl("1 of 1", txt, fixed = TRUE))
  expect_true(grepl("executor >  local (1)", txt, fixed = TRUE))
})

test_that("banner lines are header and post-process lines are footer", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("N E X T F L O W", "[a1/b2c3] count_clade (Metazoa) | 1 of 1", "Completed at: now"))
  expect_true(grepl("N E X T F L O W", RiboPilot:::nf_progress_text(s), fixed = TRUE))
  expect_true(grepl("Completed at: now", RiboPilot:::nf_progress_text(s), fixed = TRUE))
})

test_that("ANSI escapes are stripped", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    "\033[32m[a1/b2c3] build_ref_db (actinopteri) | 1 of 1\033[0m")
  expect_false(grepl("\033", RiboPilot:::nf_progress_text(s), fixed = TRUE))
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Move the parser into `R/nf_progress.R`**

Copy the body of `progress_update()` and the `nf_prog_all` render from
`app_server.R` verbatim into the three functions above, replacing the four
`reactiveVal`s with fields of a plain list. Do not change the regex or the
classification rules; this is a move.

- [ ] **Step 4: Rewire `app_server.R`**

Replace the four `reactiveVal`s with one `nf_state <- reactiveVal(nf_progress_empty())`,
`progress_update()` becomes `nf_state(nf_progress_parse(nf_state(), lines))`, and
`output$nf_prog_all <- renderText(nf_progress_text(nf_state()))`. Every reset
site (`input$run_btn`, `input$nf_start`) sets `nf_state(nf_progress_empty())`.

- [ ] **Step 5: Run the suite**

Green, including any existing Assemble-tab tests.

- [ ] **Step 6: Commit**

```bash
git add R/nf_progress.R R/app_server.R tests/testthat/test-nf-progress.R
git commit -m "refactor(app): extract the nextflow progress parser"
```

---

### Task 8: The reference run modal

**Files:**
- Modify: `R/app_modals.R`
- Test: `tests/testthat/test-reference-modal.R`

**Interfaces:**
- Produces: `reference_run_modal(n, ref_dir, headless = FALSE, sched = NULL, script = "")`
  -> a `modalDialog`. Inputs: `ref_api_key`, `ref_refseq`, `ref_nf_resume`.
  Outputs: `ref_nf_cmd_text`, `ref_nf_prog_all`. Buttons: `ref_nf_start`,
  `ref_nf_stop`, `ref_nf_close`, plus `ref_submit_job` / `ref_save_script` when
  headless. Progress div id `ref_nf_progress`, gears span id `ref_gears`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("the reference run modal carries the command, progress, and controls", {
  html <- as.character(RiboPilot:::reference_run_modal(2, "/lab/refs"))
  for (id in c("ref_api_key", "ref_refseq", "ref_nf_resume", "ref_nf_cmd_text",
               "ref_nf_prog_all", "ref_nf_start", "ref_nf_stop", "ref_nf_close",
               "ref_gears")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  expect_true(grepl("/lab/refs", html, fixed = TRUE))
  expect_true(grepl("2", html, fixed = TRUE))
})

test_that("the headless reference modal offers the job script instead of a local run", {
  html <- as.character(RiboPilot:::reference_run_modal(1, "/lab/refs", headless = TRUE,
    sched = list(name = "Slurm", cmd = "sbatch"), script = "#!/bin/bash\nnextflow run"))
  expect_true(grepl("ref_submit_job", html, fixed = TRUE))
  expect_true(grepl("ref_save_script", html, fixed = TRUE))
  expect_false(grepl("ref_nf_start", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Implement, mirroring `run_modal()`**

```r
#' Reference build modal: run the reference workflow with live progress
#'
#' The same shape as [run_modal()] so a database build looks like any other
#' RiboPilot run: command preview, progress board, Run / Stop / Close, and the
#' cluster-submit section when the app is served headless.
#' @noRd
reference_run_modal <- function(n, ref_dir, headless = FALSE, sched = NULL, script = "") {
  cluster_section <- if (headless) {
    tagList(
      hr(style = "margin:10px 0;"),
      h5("Submit to cluster", style = "margin:6px 0 4px;"),
      div(class = "text-muted", style = "font-size:0.82em; margin-bottom:4px;",
        if (is.null(sched)) {
          "No scheduler (sbatch/qsub/bsub) found on PATH. Save the script and submit it yourself; running from the login node is discouraged."
        } else {
          sprintf("%s detected. Edit the job script if needed, then Submit as job (or Save script).", sched$name)
        }),
      tags$textarea(id = "ref_submit_script_text", rows = 11,
        style = "width:100%; font-family:monospace; font-size:0.78em;", script))
  }
  run_progress <- if (!headless) {
    shinyjs::hidden(div(id = "ref_nf_progress",
      hr(style = "margin:12px 0;"), h5("Progress:", style = "margin:0 0 4px;"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit; white-space:pre-wrap;",
          textOutput("ref_nf_prog_all", container = tags$code)))))
  }
  footer <- if (headless) {
    tagList(
      if (!is.null(sched)) actionButton("ref_submit_job", "Submit as job", class = "btn-success btn-sm"),
      actionButton("ref_save_script", "Save script", class = "btn-sm"),
      actionButton("ref_nf_close", "Close", class = "btn-sm"))
  } else {
    tagList(
      actionButton("ref_nf_start", "Run from App", class = "btn-success btn-sm"),
      shinyjs::hidden(actionButton("ref_nf_stop", "Stop / Interrupt", class = "btn-sm")),
      actionButton("ref_nf_close", "Close", class = "btn-sm"))
  }
  modalDialog(
    title = div(style = "display:flex; justify-content:space-between; align-items:center; height:38px;",
      span(sprintf("Reference databases - %d clade(s) to build", n)),
      span(id = "ref_gears", class = "gears paused")),
    size = "l",
    div(class = "text-muted", style = "font-size:0.85em; margin-bottom:6px;",
      sprintf("Databases are written to %s", ref_dir)),
    textInput("ref_api_key", "NCBI API key (optional):", width = "100%"),
    shinyWidgets::prettyCheckbox("ref_refseq", "RefSeq records only", value = FALSE,
      status = "primary", inline = TRUE),
    shinyWidgets::prettyCheckbox("ref_nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("ref_nf_cmd_text")),
    run_progress,
    cluster_section,
    footer = footer
  )
}
```

- [ ] **Step 4: Run the file, then the suite**

- [ ] **Step 5: Commit**

```bash
git add R/app_modals.R tests/testthat/test-reference-modal.R
git commit -m "feat(app): reference build run modal"
```

---

### Task 9: Build launches the workflow

**Files:**
- Modify: `R/reference_build.R:286-326` (`start_reference_build`), `R/app_reference.R`
- Test: `tests/testthat/test-reference-build.R`

**Interfaces:**
- Consumes: `write_clade_manifest()`, `nf_reference_command()`, the run modal.
- Produces: `start_reference_build(con, project, ref_ids, api_key = NULL,
  refseq = FALSE, resume = TRUE, overwrite = FALSE, launcher = processx_launcher)`
  -> invisibly the log path, or `NULL` when every group was reused. Takes a
  VECTOR of `ref_ids` now, launches ONE run for all of them, and writes the same
  pid, log, and launch time onto every launched row.
  `poll_reference_builds()` is unchanged.

- [ ] **Step 1: Write the failing tests**

```r
test_that("one build run covers every group and shares its pid", {
  con <- mk_build_db(); on.exit(DBI::dbDisconnect(con))     # existing helper
  proj <- withr::local_tempdir()
  seen <- new.env()
  launcher <- function(cmd, args, log) { seen$cmd <- cmd; seen$args <- args; 4242L }
  RiboPilot:::start_reference_build(con, proj, c("asteroidea", "gastropoda"),
    launcher = launcher)
  r <- DBI::dbGetQuery(con, "SELECT ref_id, pid, reference_switch, log_file FROM reference ORDER BY ref_id")
  expect_true(all(r$pid == 4242L))
  expect_true(all(r$reference_switch == 4L))
  expect_equal(length(unique(r$log_file)), 1L)
  expect_equal(seen$cmd, "nextflow")
  expect_true("-entry" %in% seen$args)
  expect_equal(seen$args[which(seen$args == "-entry") + 1], "BUILD")
  man <- utils::read.delim(seen$args[which(seen$args == "--clades") + 1], stringsAsFactors = FALSE)
  expect_setequal(man$ref_id, c("asteroidea", "gastropoda"))
})

test_that("a group whose database already exists is wired without launching", {
  # keep the existing reuse test, passing ref_id as a length-1 vector
})
```

- [ ] **Step 2: Run and watch them fail**

- [ ] **Step 3: Implement**

Rewrite `start_reference_build()`: resolve reusable groups first (existing
logic, per group), write one manifest for the rest, build the BUILD command with
`nf_reference_command()`, launch it through `launcher("nextflow", args, log)`,
then stamp `reference_switch = 4`, the pid, the shared log path, and the launch
time onto every launched row. The log lives at
`<project>/.reference/build_<timestamp>.log`. Leave `poll_reference_builds()`
alone: it already judges each row by its own FASTAs.

- [ ] **Step 4: Wire the Build button in `R/app_reference.R`**

`input$ref_build_btn` opens `reference_run_modal()` instead of
`reference_build_modal()`; `output$ref_nf_cmd_text` renders the BUILD command
for the current inputs; `input$ref_nf_start` calls `start_reference_build()`
with an attached `processx` process, shows the progress pane, and streams output
through `nf_progress_parse()`; `ref_nf_stop` and `ref_nf_close` kill it exactly
as the Update modal does. Delete `reference_build_modal()` and its test, and
drop the conda-environment input.

- [ ] **Step 5: Run the suite**

- [ ] **Step 6: Commit**

```bash
git add R/reference_build.R R/app_reference.R R/app_modals.R tests/testthat/test-reference-build.R tests/testthat/test-reference-modal.R
git commit -m "feat(reference): build clade databases through nextflow"
```

---

### Task 10: Docs and dead host paths

**Files:**
- Modify: `dev/design/16_REFERENCE_NEXTFLOW.md` (status line), `R/build_reference_db.R`
  (roxygen), any vignette or README section describing the Reference tab
- Test: `tests/testthat/test-reference-build.R`

- [ ] **Step 1: Find what still describes the old behavior**

Run: `grep -rn "conda_env\|detached\|background build" R/*.R vignettes man README.md | grep -i refer`

- [ ] **Step 2: Update the prose**

`build_reference_db()` keeps working as the CLI entry point and keeps its host
tools; say so in its roxygen, and say that the app no longer calls it. Update
the Reference-tab section of the docs to describe the run modal. Flip the spec's
status line to "implemented".

- [ ] **Step 3: Regenerate docs and run the suite**

Run: `Rscript -e 'devtools::document(); devtools::load_all("."); testthat::test_dir("tests/testthat")'`

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "docs(reference): describe the nextflow build path"
```

---

## Self-review notes

- Spec coverage: workflow and entry points (Task 4), mounts (Task 3), work-list
  TSV (Task 2), Suggest and apply counting (Tasks 5 and 6), build launching and
  polling (Task 9), run modal including headless (Task 8), shared progress
  parser (Task 7), error table (Tasks 5, 9), migration note (Task 10).
- The stdin-detach rule appears in the wrapper (Task 4) and is already true of
  the R-side calls fixed earlier on this branch.
- Names used across tasks: `read_counts_cache`, `cache_count_fn`,
  `write_clade_manifest`, `read_count_files`, `merge_counts`,
  `nf_reference_command`, `ref_bind_option`, `nf_run_blocking`,
  `run_clade_counts`, `candidate_rungs`, `nf_progress_empty`,
  `nf_progress_parse`, `nf_progress_text`, `reference_run_modal`,
  `start_reference_build`.
