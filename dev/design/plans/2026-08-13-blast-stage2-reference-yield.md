# BLAST Stage 2: reference yield

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Derive reference element coordinates locally with barrnap and ITSx instead of
parsing NCBI feature tables, so every fetched reference is usable, and fetch only the
sequence.

**Architecture:** Reference coordinates currently come from `efetch -format ft`, which
yields a usable coordinate for only 41.6 percent of records because most rDNA submissions
carry a single `misc_RNA` with a descriptive `note` rather than per-element `product`
qualifiers. Running the pipeline's own annotation tools over the fetched sequences raises
that to 100 percent, is a strict superset of what GenBank provides, and makes the
sample-versus-reference comparison like-for-like for the first time. Because only the
sequence is then needed, the fetch drops from two requests per accession to one, and a
multi-FASTA response splits unambiguously on `>`.

**Tech Stack:** R (DBI, RSQLite, Biostrings), pybarrnap, ITSx, testthat 3e, Nextflow DSL2,
Entrez Direct.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md` (Stage 2)

## Global Constraints

- R >= 4.4.0 (`DESCRIPTION`); `docker/Dockerfile` pins `r-base=4.4.*`. No package outside
  DESCRIPTION Imports.
- ASCII only in R and Nextflow sources. No non-ASCII characters, no em dashes.
- Minimal comments. Do not narrate fixes in comments.
- This project ships no schema migrations before release.
- Never `git push`. Never `git pull`, `git fetch`, `git merge`, or `git rebase`: the
  remote branch is deliberately stale relative to this local branch and pulling corrupts
  the history.
- Run tests with:
  `Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/<file>", reporter="summary")'`
- Full suite baseline entering this plan: 2401 passed, 0 failed, 0 errors, 6 skipped.
  A rising passed count is expected as tests are added; failures and errors are not.
- **`barrnap` (the Perl one) is broken in the local conda env** (missing `Path::Tiny`).
  `pybarrnap` 0.5.2 and `ITSx` 1.1.3 both work at
  `/home/dmacguig/miniconda3/envs/ribopilot/bin`. Tests must not require the Perl barrnap;
  use the fake-tool helpers in `tests/testthat/test-rrna-kingdom.R` where a binary is
  needed.

## Measured facts this plan rests on

Running the pipeline's own `pybarrnap -k euk` and `ITSx -t all --preserve T` over the 517
reference records already fetched into the three test projects, then feeding the output
through RiboPilot's real parsers:

| | GenBank feature table | local barrnap + ITSx |
|---|---|---|
| records with >=1 usable element | 215 of 517 (41.6%) | 517 of 517 (100%) |
| records with all 5 elements | 44 (8.5%) | 330 (63.8%) |

Strict superset: all 371 accession/element pairs already stored were reproduced. Across
those 371 pairs no start or stop differs by more than 20 bp, worst case 10 bp, and
re-running the real `curate_rrna_core()` with local coordinates produced zero flag
changes. Cost: 232 records in 156 s at 8 cores, ITSx being 72 to 81 percent of it.

## The three blockers

Each fails silently and yields plausible, wrong output. Tasks 1 to 3 exist for these.

1. **Multi-FASTA collapses to one record.** `parse_itsx()` reads only `pos_lines[1]`
   (`R/annotate.R:62`). `parse_barrnap()` never extracts the GFF seqid (column 1), taking
   only columns 4, 5, 7, 9 (`R/annotate.R:45-52`), so rows from different records are
   indistinguishable. `reconcile_features()`'s `pick()` (`R/annotate.R:98`) then takes
   `df[df$element == el, ][1, ]`. Verified: real 235-record tool output through the real
   functions returns one 5-element model, no error, no warning.
2. **One record over 100 kb aborts the whole ITSx batch.** ITSx bundles HMMER's protein
   pipeline, which hard-fails past 100K. Tested with a 1.4 Mb record: fatal error, no
   positions file written. Tested with 20 normal records plus one large: same abort, all
   20 lost. `run_itsx()` (`R/annotate.R:259-273`) passes `stdout = FALSE, stderr = FALSE`
   and returns `data.frame()` when the file is absent, so this is completely silent.
   `pybarrnap` handled the same 1.4 Mb record fine, so the gate must be per tool.
3. **Minus-strand records mix coordinate frames.** barrnap reports the input frame; ITSx
   with `--preserve T` silently searches the other strand and reports flipped. Verified on
   reverse-complemented records: 5.8S came out about 1900 bp from its true position and
   ITS1/ITS2 were dropped by the spacer clipper.

---

### Task 1: Make the parsers record-aware

Additive only. Both parsers gain a `seqid` column; nothing else about their output
changes, so existing single-record callers (`reconcile_one()` at
`R/pipeline_steps.R:92-96`, `assembly_anchors()` at `R/select_paths.R:20-26`) keep working
byte-identically.

**Files:**
- Modify: `R/annotate.R` (`parse_barrnap()` at 21-53, `parse_itsx()` at 58-74)
- Test: `tests/testthat/test-reconcile.R` (extend)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `parse_barrnap(gff_lines)` returns `data.frame(seqid, element, start, stop, strand, model_cov)`.
  - `parse_itsx(pos_lines)` returns `data.frame(seqid, element, start, stop)`, one group of
    rows per input line.
  - `reconcile_features(bar, its)` is UNCHANGED. Callers that group by record must subset
    and drop `seqid` before calling it.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-reconcile.R`:

```r
test_that("parse_barrnap carries the GFF seqid", {
  gff <- c(
    "ACC1\tbarrnap\trRNA\t10\t1800\t.\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA",
    "ACC2\tbarrnap\trRNA\t20\t1900\t.\t-\t.\tName=18S_rRNA;product=18S ribosomal RNA")
  d <- RiboPilot:::parse_barrnap(gff)
  expect_equal(nrow(d), 2L)
  expect_equal(d$seqid, c("ACC1", "ACC2"))
  expect_equal(d$element, c("18S", "18S"))
  expect_equal(d$start, c(10L, 20L))
  expect_equal(d$strand, c("+", "-"))
})

test_that("parse_itsx reads every record, not just the first", {
  pos <- c(
    "ACC1\t7000 bp.\tSSU: 1-1800\tITS1: 1801-2100\t5.8S: 2101-2258\tITS2: 2259-2500\tLSU: 2501-6900",
    "ACC2\t4000 bp.\tSSU: Not found\tITS1: Not found\t5.8S: Not found\tITS2: Not found\tLSU: 1-3950")
  d <- RiboPilot:::parse_itsx(pos)
  expect_equal(sort(unique(d$seqid)), c("ACC1", "ACC2"))
  expect_equal(nrow(d[d$seqid == "ACC1", ]), 5L)
  expect_equal(nrow(d[d$seqid == "ACC2", ]), 1L)
  expect_equal(d$element[d$seqid == "ACC2"], "28S")
})

test_that("a single-record positions line still parses exactly as before", {
  pos <- "S1\t7000 bp.\tSSU: 1-1800\tITS1: 1801-2100\t5.8S: 2101-2258\tITS2: 2259-2500\tLSU: 2501-6900"
  d <- RiboPilot:::parse_itsx(pos)
  expect_equal(d$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
  expect_equal(d$start, c(1L, 1801L, 2101L, 2259L, 2501L))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-reconcile.R", reporter="summary")'`
Expected: the first two FAIL (`d$seqid` is NULL; `parse_itsx` returns only ACC1's rows).
The third PASSES already.

- [ ] **Step 3: Add seqid to parse_barrnap**

In `R/annotate.R`, change the returned frame so `seqid` is the first column:

```r
  data.frame(
    seqid = vapply(f, `[`, character(1), 1)[keep],
    element = el[keep],
    start = as.integer(vapply(f, `[`, character(1), 4))[keep],
    stop = as.integer(vapply(f, `[`, character(1), 5))[keep],
    strand = vapply(f, `[`, character(1), 7)[keep],
    model_cov = model_cov,
    stringsAsFactors = FALSE
  )
```

Update the roxygen `@return` above the function to list `seqid`.

- [ ] **Step 4: Make parse_itsx read every line**

Replace the body of `parse_itsx()` with:

```r
parse_itsx <- function(pos_lines) {
  if (length(pos_lines) == 0) {
    return(data.frame())
  }
  pos_lines <- pos_lines[nzchar(trimws(pos_lines))]
  labmap <- c("SSU" = "18S", "ITS1" = "ITS1", "5.8S" = "5.8S", "ITS2" = "ITS2", "LSU" = "28S")
  out <- list()
  for (ln in pos_lines) {
    fields <- strsplit(ln, "\t")[[1]]
    if (length(fields) == 0) next
    sid <- fields[1]
    for (fd in fields) {
      m <- regmatches(fd, regexec("^\\s*(SSU|ITS1|5\\.8S|ITS2|LSU):\\s*(\\d+)-(\\d+)", fd))[[1]]
      if (length(m) == 4) {
        out[[length(out) + 1]] <- data.frame(
          seqid = sid, element = unname(labmap[m[2]]),
          start = as.integer(m[3]), stop = as.integer(m[4]),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}
```

Update its roxygen `@return` to list `seqid`.

- [ ] **Step 5: Confirm reconcile_features is unaffected**

`reconcile_features()` filters on `element` and takes `[1, ]`, so an extra column changes
nothing for single-record input. Verify by running the reconcile and kingdom suites:

`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); for (f in c("test-reconcile.R","test-rrna-kingdom.R")) testthat::test_file(file.path("tests/testthat", f), reporter="summary")'`
Expected: PASS.

If `reconcile_features()` turns out to reference columns positionally anywhere, do NOT
change it; instead make the callers drop `seqid` before calling. Report if you hit this.

- [ ] **Step 6: Run the full suite and commit**

`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); r <- testthat::test_dir("tests/testthat", reporter="silent"); d <- as.data.frame(r); cat("failed:", sum(d$failed), "error:", sum(d$error), "skipped:", sum(d$skipped), "passed:", sum(d$passed), "\n")'`
Expected: 0 failed, 0 errors.

```bash
git add R/annotate.R tests/testthat/test-reconcile.R
git commit -m "feat: carry the sequence id through the barrnap and ITSx parsers"
```

---

### Task 2: Make ITSx failures visible and keep oversized records away from it

`run_itsx()` discards both streams and returns an empty frame when the positions file is
missing, so a fatal abort is indistinguishable from "no ITS found". `run_rrna_tool()`
already does this correctly for barrnap; copy that shape.

**Files:**
- Modify: `R/annotate.R` (`run_itsx()` at 259-273)
- Test: `tests/testthat/test-rrna-kingdom.R` (extend)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `ITSX_MAX_BP` constant, 100000.
  - `run_itsx(fasta, cpu = 2, taxon = "all")` unchanged in return type, but warns on
    non-zero exit or a missing positions file, and accepts an ITSx `-t` taxon.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-rrna-kingdom.R`:

```r
test_that("run_itsx warns when ITSx exits non-zero", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 1L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  expect_warning(RiboPilot:::run_itsx(fa), "ITSx")
})

test_that("run_itsx warns when no positions file is produced", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  expect_warning(RiboPilot:::run_itsx(fa), "positions")
})

test_that("run_itsx passes the requested taxon through to -t", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa, taxon = "M"))
  args <- fake_args(d, "ITSx")
  expect_true(any(args == "M"))
})

test_that("ITSX_MAX_BP is the documented HMMER ceiling", {
  expect_equal(RiboPilot:::ITSX_MAX_BP, 100000)
})
```

`mk_fake_tool()`, `fake_args()` and `mk_fasta()` are already defined at the top of that
file; reuse them, do not redefine.

- [ ] **Step 2: Run tests to verify they fail**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-rrna-kingdom.R", reporter="summary")'`
Expected: the four new tests FAIL (no warning emitted; `ITSX_MAX_BP` not found).

- [ ] **Step 3: Add the constant**

In `R/annotate.R`, above `run_itsx()`:

```r
#' Largest record ITSx will accept.
#'
#' ITSx runs HMMER's protein pipeline, which refuses a target over 100K and
#' aborts the whole run rather than skipping the record, writing no positions
#' file at all. pybarrnap has no such limit.
#' @noRd
ITSX_MAX_BP <- 100000
```

- [ ] **Step 4: Rewrite run_itsx**

```r
run_itsx <- function(fasta, cpu = 2, taxon = "all") {
  bin <- unname(Sys.which("ITSx"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  pre <- tempfile("rp_itsx_")
  err <- tempfile("rp_itsx_err_")
  on.exit(unlink(err), add = TRUE)
  st <- suppressWarnings(system2(bin, c("-i", shQuote(fasta), "-o", shQuote(pre),
    "--save_regions", "none", "--graphical", "F", "--cpu", cpu,
    "-t", taxon, "--preserve", "T"),
    stdout = FALSE, stderr = err))
  tail_msg <- function() {
    tryCatch(paste(utils::tail(readLines(err, warn = FALSE), 3L), collapse = " | "),
      error = function(e) "")
  }
  if (!identical(as.integer(st), 0L)) {
    warning(sprintf("ITSx failed on %s: %s", basename(fasta), tail_msg()), call. = FALSE)
    return(data.frame())
  }
  pf <- paste0(pre, ".positions.txt")
  if (!file.exists(pf)) {
    warning(sprintf("ITSx wrote no positions file for %s: %s", basename(fasta), tail_msg()),
      call. = FALSE)
    return(data.frame())
  }
  parse_itsx(readLines(pf, warn = FALSE))
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run the kingdom file, then the full suite as in Task 1 Step 6.
Expected: PASS, 0 failed, 0 errors.

- [ ] **Step 6: Commit**

```bash
git add R/annotate.R tests/testthat/test-rrna-kingdom.R
git commit -m "fix: surface ITSx failures instead of returning an empty frame"
```

---

### Task 3: Annotate a multi-record FASTA, handling strand and the size ceiling

The composing task. Runs each tool once over a whole multi-FASTA, demultiplexes by
sequence id, and reconciles per record. Handles all three blockers.

Strand handling: barrnap reports the input frame, ITSx with `--preserve T` reports the
flipped frame on a reverse-oriented record. Rather than reorienting stored sequences,
split the input by the strand barrnap already reported, run ITSx separately on a
reverse-complemented copy of the minus-strand records, and map those coordinates back with
`start' = L - stop + 1`, `stop' = L - start + 1`.

**Files:**
- Modify: `R/annotate.R` (append)
- Test: `tests/testthat/test-annotate-multi.R` (create)

**Interfaces:**
- Consumes: `parse_barrnap()`, `parse_itsx()`, `reconcile_features()`, `rrna_caller()`,
  `run_rrna_tool()`, `rrna_kingdom_for()`, `run_itsx()`, `ITSX_MAX_BP`.
- Produces:
  `annotate_records(fasta, kingdom = "euk", cpu = 2, taxon = "all")` returning
  `data.frame(seqid, element, type, start, stop, length, strand, model_cov, tool)`, zero
  rows when nothing is annotatable. One group of rows per input record that yielded any
  feature.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-annotate-multi.R`:

```r
# annotate_records() maps tool output back to the right record and the right
# coordinate frame. Tool output is injected rather than invoking real binaries.

test_that("annotate_records demultiplexes barrnap and ITSx output by record", {
  bar <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("18S", "18S"),
    start = c(10L, 20L), stop = c(1810L, 1920L),
    strand = c("+", "+"), model_cov = c(100, 100), stringsAsFactors = FALSE)
  its <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("ITS1", "ITS2"),
    start = c(1811L, 1921L), stop = c(2000L, 2100L), stringsAsFactors = FALSE)
  out <- RiboPilot:::reconcile_records(bar, its)
  expect_equal(sort(unique(out$seqid)), c("ACC1", "ACC2"))
  # the load-bearing property is that a record's coordinates come from its own
  # rows, not from whichever row sorted first
  expect_equal(out$start[out$seqid == "ACC1" & out$element == "18S"], 10L)
  expect_equal(out$start[out$seqid == "ACC2" & out$element == "18S"], 20L)
  expect_true("ITS1" %in% out$element[out$seqid == "ACC1"])
  expect_true("ITS2" %in% out$element[out$seqid == "ACC2"])
})

test_that("minus-strand ITSx coordinates are mapped back to the input frame", {
  # a 1000 bp record; ITSx searched the reverse complement and reported 100-200
  # there, which is 801-901 in the input frame
  flipped <- data.frame(seqid = "ACC1", element = "5.8S",
    start = 100L, stop = 200L, stringsAsFactors = FALSE)
  back <- RiboPilot:::flip_itsx_coords(flipped, c(ACC1 = 1000L))
  expect_equal(back$start, 801L)
  expect_equal(back$stop, 901L)
})

test_that("flip_itsx_coords leaves records it has no length for untouched", {
  d <- data.frame(seqid = "ACCX", element = "5.8S", start = 5L, stop = 10L,
    stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::flip_itsx_coords(d, c(ACC1 = 1000L)), d)
})

test_that("records over the ITSx ceiling are excluded from the ITSx input only", {
  lens <- c(SMALL = 5000L, HUGE = 200000L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "SMALL")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-annotate-multi.R", reporter="summary")'`
Expected: FAIL, `reconcile_records` / `flip_itsx_coords` / `itsx_eligible` not found.

- [ ] **Step 3: Add the three helpers**

Append to `R/annotate.R`:

```r
#' Sequence ids short enough for ITSx
#' @param lens Named integer vector of sequence lengths.
#' @noRd
itsx_eligible <- function(lens) {
  names(lens)[!is.na(lens) & lens <= ITSX_MAX_BP]
}

#' Map ITSx coordinates from the reverse-complement frame back to the input frame
#' @param its parse_itsx() output for reverse-complemented records.
#' @param lens Named integer vector of input-frame sequence lengths.
#' @noRd
flip_itsx_coords <- function(its, lens) {
  if (is.null(its) || nrow(its) == 0) return(its)
  L <- lens[its$seqid]
  ok <- !is.na(L)
  if (!any(ok)) return(its)
  s <- its$start; e <- its$stop
  its$start[ok] <- as.integer(L[ok] - e[ok] + 1L)
  its$stop[ok] <- as.integer(L[ok] - s[ok] + 1L)
  its
}

#' Reconcile per record
#' @param bar,its parse_barrnap()/parse_itsx() output carrying a seqid column.
#' @return reconcile_features() output with a leading seqid column.
#' @noRd
reconcile_records <- function(bar, its) {
  ids <- unique(c(
    if (!is.null(bar) && nrow(bar) > 0) bar$seqid else character(0),
    if (!is.null(its) && nrow(its) > 0) its$seqid else character(0)))
  if (length(ids) == 0) return(data.frame())
  out <- list()
  drop_seqid <- function(d) {
    if (is.null(d) || nrow(d) == 0) return(data.frame())
    d[, setdiff(names(d), "seqid"), drop = FALSE]
  }
  for (sid in ids) {
    b <- drop_seqid(if (!is.null(bar) && nrow(bar) > 0) bar[bar$seqid == sid, , drop = FALSE] else NULL)
    i <- drop_seqid(if (!is.null(its) && nrow(its) > 0) its[its$seqid == sid, , drop = FALSE] else NULL)
    r <- reconcile_features(b, i)
    if (is.null(r) || nrow(r) == 0) next
    out[[length(out) + 1]] <- cbind(seqid = sid, r, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}
```

- [ ] **Step 4: Add the driver**

Append to `R/annotate.R`:

```r
#' Annotate every record in a multi-record fasta
#'
#' Runs the rRNA anchor caller once over the whole file and ITSx once per strand
#' group, then reconciles per record. Records longer than [ITSX_MAX_BP] are kept
#' but get anchors only, because ITSx aborts the entire run on an oversized
#' target rather than skipping it.
#'
#' @param fasta Path to a multi-record fasta.
#' @param kingdom Anchor-caller kingdom.
#' @param cpu Threads for ITSx.
#' @param taxon ITSx `-t` profile set.
#' @return data.frame(seqid, element, type, start, stop, length, strand, model_cov, tool).
#' @noRd
annotate_records <- function(fasta, kingdom = "euk", cpu = 2, taxon = "all") {
  seqs <- Biostrings::readDNAStringSet(fasta)
  if (length(seqs) == 0) return(data.frame())
  names(seqs) <- sub("\\s.*$", "", names(seqs))
  lens <- stats::setNames(as.integer(Biostrings::width(seqs)), names(seqs))

  caller <- rrna_caller()
  bar <- data.frame()
  if (!is.null(caller)) {
    gff <- run_rrna_tool(caller$bin, caller$tool, fasta,
      rrna_kingdom_for(kingdom, caller$tool))
    bar <- parse_barrnap(gff)
  }

  minus <- if (nrow(bar) > 0) unique(bar$seqid[bar$strand == "-"]) else character(0)
  eligible <- itsx_eligible(lens)
  plus_ids <- setdiff(intersect(names(seqs), eligible), minus)
  minus_ids <- intersect(minus, eligible)

  run_on <- function(ids, revcomp) {
    if (length(ids) == 0) return(data.frame())
    sub <- seqs[ids]
    if (revcomp) sub <- Biostrings::reverseComplement(sub)
    f <- tempfile(fileext = ".fa")
    on.exit(unlink(f), add = TRUE)
    Biostrings::writeXStringSet(sub, f)
    run_itsx(f, cpu = cpu, taxon = taxon)
  }

  its_plus <- run_on(plus_ids, FALSE)
  its_minus <- flip_itsx_coords(run_on(minus_ids, TRUE), lens)
  its <- rbind(its_plus, its_minus)

  reconcile_records(bar, its)
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run the new file, then `test-reconcile.R`, `test-rrna-kingdom.R`, then the full suite.
Expected: PASS, 0 failed, 0 errors.

- [ ] **Step 6: Commit**

```bash
git add R/annotate.R tests/testthat/test-annotate-multi.R
git commit -m "feat: annotate a multi-record fasta per record and per strand"
```

---

### Task 4: Store references from sequence plus local annotation

`write_blast_ref()` takes one accession and a feature-table file. Add a batch sibling that
takes a multi-FASTA and the annotation frame from Task 3. Keep the old function: it is
exported and covered by `tests/testthat/test-reference-fetch.R` and
`tests/testthat/test-write-atomicity.R`.

**Files:**
- Modify: `R/blast_ref_fetch.R` (append)
- Test: `tests/testthat/test-reference-fetch.R` (extend)

**Interfaces:**
- Consumes: `annotate_records()` from Task 3.
- Produces:
  `write_blast_refs(con, fasta_file, ann = NULL)` returning, invisibly, the number of
  accessions stored. Reads the multi-FASTA itself, derives `species` from each defline,
  writes `blast_ref_sequences` and `blast_ref_annotations`, replacing any existing rows for
  each accession, all in one transaction per accession.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-reference-fetch.R`:

```r
test_that("write_blast_refs stores every record in a multi-fasta", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(
    ">ACC1.1 Danio rerio 18S ribosomal RNA gene",
    paste(rep("ACGT", 50), collapse = ""),
    ">ACC2.1 Danio rerio 28S ribosomal RNA gene",
    paste(rep("TGCA", 60), collapse = "")), f)
  ann <- data.frame(
    seqid = c("ACC1.1", "ACC2.1"), element = c("18S", "28S"),
    type = c("rRNA", "rRNA"), start = c(1L, 1L), stop = c(200L, 240L),
    length = c(200L, 240L), strand = c("+", "+"),
    model_cov = c(100, 100), tool = c("barrnap", "barrnap"),
    stringsAsFactors = FALSE)
  n <- RiboPilot:::write_blast_refs(con, f, ann)
  expect_equal(as.integer(n), 2L)
  s <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_sequences ORDER BY accession")
  expect_equal(s$accession, c("ACC1.1", "ACC2.1"))
  expect_equal(nchar(s$sequence[1]), 200L)
  expect_equal(s$species[1], "Danio rerio")
  a <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_annotations ORDER BY accession")
  expect_equal(a$element, c("18S", "28S"))
})

test_that("write_blast_refs does not concatenate records into one sequence", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "AAAA", ">B.1 y", "CCCC"), f)
  RiboPilot:::write_blast_refs(con, f, NULL)
  s <- DBI::dbGetQuery(con, "SELECT accession, sequence FROM blast_ref_sequences ORDER BY accession")
  expect_equal(s$sequence, c("AAAA", "CCCC"))
})

test_that("write_blast_refs replaces rather than duplicating on a re-run", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "AAAA"), f)
  RiboPilot:::write_blast_refs(con, f, NULL)
  RiboPilot:::write_blast_refs(con, f, NULL)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_sequences")$n, 1L)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-reference-fetch.R", reporter="summary")'`
Expected: FAIL, `write_blast_refs` not found.

- [ ] **Step 3: Implement**

Append to `R/blast_ref_fetch.R`:

```r
#' Store every record of a multi-record fasta as a reference
#'
#' Sequence comes from the fasta; element coordinates come from `ann` rather than
#' from an NCBI feature table. An accession with no rows in `ann` is still stored,
#' with no annotations.
#'
#' @param con Open project DB connection.
#' @param fasta_file Multi-record fasta, deflines `>ACC.V description`.
#' @param ann Annotation frame carrying seqid, element, type, start, stop, strand,
#'   or NULL.
#' @return (invisibly) the number of accessions stored.
#' @export
write_blast_refs <- function(con, fasta_file, ann = NULL) {
  if (!is_fasta_file(fasta_file)) return(invisible(0L))
  seqs <- Biostrings::readDNAStringSet(fasta_file)
  if (length(seqs) == 0) return(invisible(0L))
  ts <- as.integer(Sys.time())
  full <- names(seqs)
  acc <- sub("\\s.*$", "", full)
  title <- trimws(sub("^\\S+\\s*", "", full))
  n <- 0L
  for (i in seq_along(seqs)) {
    s <- as.character(seqs[[i]])
    if (!nzchar(s)) next
    a <- ann
    rows <- if (!is.null(a) && nrow(a) > 0) a[a$seqid == acc[i], , drop = FALSE] else NULL
    DBI::dbWithTransaction(con, {
      DBI::dbExecute(con, "DELETE FROM blast_ref_sequences WHERE accession = ?",
        params = list(acc[i]))
      DBI::dbAppendTable(con, "blast_ref_sequences", data.frame(
        accession = acc[i], species = extract_species(title[i]), title = title[i],
        taxid = NA_character_, lineage = NA_character_, sequence = s,
        length = nchar(s), time_stamp = ts, stringsAsFactors = FALSE))
      DBI::dbExecute(con, "DELETE FROM blast_ref_annotations WHERE accession = ?",
        params = list(acc[i]))
      if (!is.null(rows) && nrow(rows) > 0) {
        DBI::dbAppendTable(con, "blast_ref_annotations", data.frame(
          accession = acc[i], element = rows$element, type = rows$type,
          start = rows$start, stop = rows$stop, strand = rows$strand,
          time_stamp = ts, stringsAsFactors = FALSE))
      }
    })
    n <- n + 1L
  }
  invisible(n)
}
```

`Biostrings` is already declared in DESCRIPTION Imports (line 17), so no dependency change
is needed here or in Task 3.

- [ ] **Step 4: Run tests, regenerate docs, commit**

Run the reference-fetch file, then the full suite.
Then `Rscript -e 'roxygen2::roxygenise()'` and include the regenerated NAMESPACE and
`man/write_blast_refs.Rd`.

```bash
git add R/blast_ref_fetch.R NAMESPACE man/ tests/testthat/test-reference-fetch.R
git commit -m "feat: store references from a multi-fasta plus local annotation"
```

---

### Task 5: Fetch sequence only, in batches, with one R invocation per sample

**Files:**
- Modify: `inst/nextflow/modules/blast_ref.nf` (`blast_ref_fetch` process body)
- Create: `inst/nextflow/bin/rp_ref_store.R`
- Test: `tests/testthat/test-reference-fetch.R` (extend, static module checks)

**Interfaces:**
- Consumes: `write_blast_refs()` (Task 4), `annotate_records()` (Task 3).
- Produces: `rp_ref_store.R <db> <fasta> [kingdom] [cpu]`, which annotates the fasta and
  stores every record.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-reference-fetch.R`:

```r
test_that("the fetch module asks for fasta only, not feature tables", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- readLines(p, warn = FALSE)
  expect_false(any(grepl("-format ft", txt, fixed = TRUE)))
  expect_true(any(grepl("-format fasta", txt, fixed = TRUE)))
})

test_that("the fetch module batches accessions rather than looping one at a time", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("rp_ref_store.R", txt, fixed = TRUE))
  expect_false(grepl("rp_ref_fetch.R", txt, fixed = TRUE))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: both FAIL (the module still uses `-format ft` and `rp_ref_fetch.R`).

- [ ] **Step 3: Add the store script**

Create `inst/nextflow/bin/rp_ref_store.R`, following the shape of
`inst/nextflow/bin/rp_rederive.R`:

```r
#!/usr/bin/env Rscript
# Annotate a multi-record reference fasta and store it.
# Args: <db> <fasta> [kingdom] [cpu]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_store.R <db> <fasta> [kingdom] [cpu]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  kingdom <- if (length(a) >= 3 && nzchar(a[3])) a[3] else "euk"
  cpu <- if (length(a) >= 4 && nzchar(a[4])) as.integer(a[4]) else 2L
  ann <- annotate_records(a[2], kingdom = kingdom, cpu = cpu)
  n <- write_blast_refs(con, a[2], ann)
  cat("references stored:", n, "\n")
})
```

`chmod +x inst/nextflow/bin/rp_ref_store.R`

`annotate_records()` is `@noRd`, so `rp_env.R` must be able to see it. Check how
`rp_env.R` exposes internals: if it uses `library(RiboPilot)` as a fallback, mark
`annotate_records()` `@export` and regenerate docs, matching what Stage 1 did for
`derive_blast_products()`. Report which route you took.

- [ ] **Step 4: Rewrite the fetch process body**

Replace the `for acc in \$ACCS` loop in `inst/nextflow/modules/blast_ref.nf` with a
batched fetch. Keep the existing early return for a sample with no candidates, the
retry backoff, and the stored-versus-wanted guard.

```bash
    ACCS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT DISTINCT c.accession FROM blast_ref_candidates c LEFT JOIN blast_ref_sequences s ON s.accession = c.accession WHERE c.ID = '${id}' AND s.accession IS NULL;")
    if [ -z "\$ACCS" ]; then echo "nothing to fetch for ${id}"; exit 0; fi
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 30 )); fi

    printf '%s\\n' \$ACCS > wanted.txt
    WANTED=\$(grep -c . wanted.txt)

    # one request per chunk instead of two per accession
    : > refs.fa
    split -l 200 -d -a 4 wanted.txt chunk_
    for f in chunk_*; do
        IDS=\$(paste -sd, "\$f")
        if ! efetch -db nuccore -id "\$IDS" -format fasta >> refs.fa; then
            echo "efetch failed for chunk \$f" >&2
        fi
        sleep \$PACE
    done

    if [ ! -s refs.fa ]; then
        echo "no reference sequence returned for ${id}" >&2
        exit 1
    fi
    grep -q -e '<ERROR>' -e '<eFetchResult>' -e 'F a i l e d' refs.fa && { echo "efetch error body for ${id}" >&2; exit 1; }

    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_store.R '${params.db}' refs.fa euk ${task.cpus ?: 2}

    GOT=\$(grep -c '^>' refs.fa)
    echo "reference fetch for ${id}: returned \${GOT} of \${WANTED} requested accession(s)"
```

Keep the `PACE` selection line above this block. Delete the now-unused per-accession
`have=` lookup, both `efetch` calls, the `.ft` checks, and the `rp_ref_fetch.R` call. Do
NOT delete `inst/nextflow/bin/rp_ref_fetch.R` itself or `write_blast_ref()`; both stay for
the single-record path and its tests.

Map records back to accessions by the accession parsed from each defline, never by
position: efetch silently omits unresolvable ids and EDirect re-sorts. `write_blast_refs()`
already keys on the defline, so this holds as long as nothing reorders records.

- [ ] **Step 5: Verify and commit**

Run the reference-fetch tests and the full suite. If `nextflow lint` is available, run it
on the changed module and report the result.

```bash
git add inst/nextflow/modules/blast_ref.nf inst/nextflow/bin/rp_ref_store.R tests/testthat/test-reference-fetch.R
git commit -m "feat: fetch reference sequence only, batched, one store call per sample"
```

---

### Task 6: Refuse to store a degenerate reference alignment

With local annotation, ITS1 and ITS2 coordinates exist for far more references than
before, and `pairwise_nt()` is ends-free, so two hypervariable spacers from different
families can collapse to a near-zero-length alignment. Recomputing on
`full_project_docker` produced 24 new rows of which 8 were degenerate: `pid` NaN, which
lands in SQLite as NULL, or `aln_cov` 0.002.

The codebase already learned this: `R/params_rrna.R:61-67` documents a 1277 bp unit
scoring 100.00 percent against a 6310 bp one over 0.2 percent of it, and pairs
`min_isomer_anchor_pid` with `min_isomer_anchor_coverage` and `min_isomer_anchor_bp`.

**Files:**
- Modify: `R/blast_ref_align.R` (`compute_ref_alignment()`)
- Test: `tests/testthat/test-ref-alignment.R` (extend)

**Interfaces:**
- Consumes: nothing.
- Produces: `REF_ALN_MIN_COV` (0.5) and `REF_ALN_MIN_BP` (50) constants;
  `compute_ref_alignment()` skips a row below either floor.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-ref-alignment.R`, reusing `mk_aln_db()` already defined
there:

```r
test_that("a degenerate element alignment is not stored", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  # replace the reference with a short unrelated spacer so the ends-free
  # alignment covers almost none of the sample element
  DBI::dbExecute(con, "UPDATE blast_ref_sequences SET sequence = ? WHERE accession = 'ACC1'",
    params = list(paste(rep("TTTTTTTTTT", 3), collapse = "")))
  DBI::dbExecute(con, "UPDATE blast_ref_annotations SET stop = 30 WHERE accession = 'ACC1'")
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_alignment")$n, 0L)
})

test_that("a good alignment is still stored", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 1L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Expected: the first test FAILS (a row is written with tiny coverage); the second PASSES.

- [ ] **Step 3: Add the floors**

In `R/blast_ref_align.R`, above `compute_ref_alignment()`:

```r
#' Floors an element alignment must clear to be stored.
#'
#' pairwise_nt() is ends-free, so two unrelated hypervariable spacers can align
#' over a handful of bases at high identity. Storing that reports a similarity
#' that does not exist.
#' @noRd
REF_ALN_MIN_COV <- 0.5
REF_ALN_MIN_BP <- 50
```

Inside the element loop, after `aln_cov` is computed and before `dbAppendTable`, add:

```r
    if (is.na(pid) || is.na(aln_cov) ||
        aln_cov < REF_ALN_MIN_COV || nchar(gsub("-", "", s_aln)) < REF_ALN_MIN_BP) {
      next
    }
```

- [ ] **Step 4: Run tests and check the real effect**

Run the ref-alignment file and the full suite.

Then measure the effect on real data, on a COPY:

```bash
cp /home/dmacguig/Documents/DJMCompBio/RiboPilot_test/full_project_docker/.sqlite /tmp/aln_floor.sqlite
Rscript -e '
  suppressMessages(pkgload::load_all(".", quiet=TRUE))
  con <- DBI::dbConnect(RSQLite::SQLite(), "/tmp/aln_floor.sqlite")
  cat("rows before:", DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_alignment")$n, "\n")
  u <- DBI::dbGetQuery(con, "SELECT DISTINCT ID, path, scaffold FROM blast_ref_alignment")
  for (i in seq_len(nrow(u))) RiboPilot:::compute_ref_alignment(con, u$ID[i], u$path[i], u$scaffold[i])
  d <- DBI::dbGetQuery(con, "SELECT COUNT(*) n, SUM(pid IS NULL) nullpid, MIN(aln_cov) mincov FROM blast_ref_alignment")
  cat("rows after:", d$n, " null pid:", d$nullpid, " min aln_cov:", d$mincov, "\n")
  DBI::dbDisconnect(con)'
rm -f /tmp/aln_floor.sqlite
```

Expected: no NULL `pid` rows remain and `min aln_cov` is at or above 0.5. Report the
numbers. If rows with NULL `pid` survive, the guard is in the wrong place; report rather
than adjusting the constants.

- [ ] **Step 5: Commit**

```bash
git add R/blast_ref_align.R tests/testthat/test-ref-alignment.R
git commit -m "fix: refuse to store an element alignment below coverage and length floors"
```

---

## After this plan

Stage 2 is complete when all six tasks are committed and the full suite shows no new
failures.

Not in scope, and deliberately deferred:

- **The end-to-end yield claim is not verified by this plan.** The 41.6 to 100 percent
  figure was measured offline on already-fetched sequences. Confirming it through the real
  pipeline needs a WF1 run, which needs remote BLAST. Do that before Stage 3.
- **Stage 3** (shared cross-project reference store under `ribopilot_ref_dir()`, including
  the read-only SQLite open path) and **Stage 4** (the shipped local database, the
  `blast_local` / `blast_remote_escalate` split, and the five escalation triggers).
- **The API key still lands in `.command.run`.** Stage 1 removed it from `.command.sh`
  only. A real fix needs Nextflow secrets or a key file outside the work directory.
- `docker/Dockerfile` installs `pybarrnap`, `barrnap` and `itsx` unpinned while
  `environment.yml` pins them, despite a comment claiming the Dockerfile mirrors it. Worth
  reconciling, unrelated to this plan.
