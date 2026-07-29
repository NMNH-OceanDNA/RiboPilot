# Three-State Hit Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the single-hit binary nuclear-rDNA test with a three-state title classification and a bitscore-weighted consensus, and repair existing projects offline.

**Architecture:** A pure classifier (`R/hit_class.R`) turns a GenBank defline into `rdna` / `other` / `uninformative`, and a pure consensus turns a unit's hit frame into `1L` / `0L` / `NA`. `blast_hits` gains a `hit_class` column; because classification is a pure function of the stored title, a v13 migration backfills it and recomputes `blast_screen.is_nuclear_rrna` without re-running the remote BLAST. The parse path then writes `hit_class` going forward.

**Tech Stack:** R 4.4, DBI + RSQLite, testthat 3e, roxygen2 (Rd only).

**Spec:** `dev/design/10_HIT_CLASSIFICATION.md`

## Global Constraints

- ASCII only in all R code, tests, fixtures and docs. No em dashes, no smart quotes, no non-ASCII bytes.
- Comments minimal, explaining WHY rather than narrating a change.
- `NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only.
- `RIBOPILOT_SCHEMA_VERSION` moves from `12L` to `13L` in this plan and nowhere else.
- Migrations only ever add tables and columns. Never drop, never retype.
- The three class strings are exactly `rdna`, `other`, `uninformative`, lower case.
- Classification precedence is fixed and must not be reordered: organellar veto, then prokaryotic-without-eukaryotic-marker, then nuclear rDNA, then assembly-level, then `other` as the default.
- This layer answers "is this record nuclear rDNA", NOT "is this the host". A fungal or human rDNA record is correctly `rdna`. Say so in the file header comment.
- Run tests from the repo root: `Rscript -e 'testthat::test_local()'`. Baseline is 505 passing, 0 failing.
- Branch is `contamination-phase3`, already created from `main` at `f9c5889`.

---

## File Structure

**Create:**
- `R/hit_class.R` - the pure classifier and the pure consensus. No DBI.
- `tests/testthat/fixtures/hit_titles.csv` - the labelled corpus.
- `tests/testthat/test-hit-class.R` - classifier and consensus tests.

**Modify:**
- `R/constants.R` - broaden `NUCLEAR_RRNA_TITLE_RE`, add the prokaryotic, eukaryotic-marker and assembly patterns. `ORGANELLAR_TITLE_RE` is unchanged.
- `R/migrate_db.R` - `hit_class` in `BLAST_HITS_DDL`, the v13 block, `backfill_hit_class()`, version bump.
- `R/blast_genbank.R` - write `hit_class`, filter the reference panel on it, and take the screen flag from the consensus.
- `R/triage_rules.R` - the `no_screen` reason string.
- `tests/testthat/test-migrate.R`, `tests/testthat/test-blast-parse.R` - coverage for the above.
- `vignettes/Curation-Flags.Rmd` - one sentence.

---

### Task 1: The pure classifier and its corpus

**Files:**
- Create: `R/hit_class.R`, `tests/testthat/fixtures/hit_titles.csv`, `tests/testthat/test-hit-class.R`
- Modify: `R/constants.R` (the `NUCLEAR_RRNA_TITLE_RE` block at lines 50-61)

**Interfaces:**
- Consumes: `%||%` from `R/utils.R`; `ORGANELLAR_TITLE_RE` from `R/constants.R`, unchanged.
- Produces:
  - `HIT_CLASSES` - `c("rdna", "other", "uninformative")`
  - `classify_hit_title(title)` -> character vector, same length as `title`, each element one of `HIT_CLASSES`
  - `NUCLEAR_RRNA_TITLE_RE`, `PROKARYOTIC_RRNA_RE`, `EUKARYOTIC_RRNA_MARKER_RE`, `ASSEMBLY_TITLE_RE` in `R/constants.R`

- [ ] **Step 1: Write the fixture corpus**

Create `tests/testthat/fixtures/hit_titles.csv`. The 26 `uninformative` rows are real titles from a ten-sample fish project and are the regression this work exists to fix; do not edit their text, including the `unplaced-saffold` typo, which is present in the GenBank record itself.

```csv
title,expected,source,note
"Micropterus salmoides isolate LMB-PRFRI chromosome Unknown",uninformative,real,chromosome-level host assembly
"Scomber australasicus DNA, Alternate haplotype, chromosome 5, sequence",uninformative,real,haplotype-resolved host assembly
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_103, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_116, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_26, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_28, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_36, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_43, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_89, sequence",uninformative,real,
"Scomber australasicus DNA, Alternate haplotype, unplaced-scaffold_91, sequence",uninformative,real,
"Scomber australasicus DNA, Primary haplotype, chromosome 8, sequence",uninformative,real,
"Scomber australasicus DNA, Primary haplotype, unplaced-scaffold_27, sequence",uninformative,real,
"Scomber australasicus DNA, Primary haplotype, unplaced-scaffold_30, sequence",uninformative,real,
"Scomber australasicus DNA, Primary haplotype, unplaced-scaffold_32, sequence",uninformative,real,
"Scomber australasicus DNA, Primary haplotype, unplaced-scaffold_33, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, chromosome 1, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, chromosome 17, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, chromosome 20, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, unplaced-scaffold_25, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, unplaced-scaffold_28, sequence",uninformative,real,
"Scomber japonicus DNA, alternate haplotype, unplaced-scaffold_33, sequence",uninformative,real,
"Scomber japonicus DNA, primary haplotype, chromosome 11, sequence",uninformative,real,
"Scomber japonicus DNA, primary haplotype, chromosome 17, sequence",uninformative,real,
"Scomber japonicus DNA, primary haplotype, chromosome 4, sequence",uninformative,real,
"Scomber japonicus DNA, primary haplotype, unplaced-saffold_25, sequence",uninformative,real,typo present in the GenBank record
"Scomber japonicus DNA, primary haplotype, unplaced-saffold_26, sequence",uninformative,real,typo present in the GenBank record
"Cleisthenes herzensteini 18S ribosomal RNA gene, complete sequence",rdna,real,18S complete
"Gonorynchus greyi isolate 94609 18S ribosomal RNA gene, partial sequence",rdna,real,18S partial
"Kareius bicoloratus isolate QD1k 18S ribosomal RNA gene, complete sequence",rdna,real,
"Micromesistius poutassou 28S ribosomal RNA (AC3HAE_023566), rRNA",rdna,real,28S with no gene keyword
"Micromesistius poutassou 28S ribosomal RNA (AC3HAE_040873), rRNA",rdna,real,28S with no gene keyword
"Abantennarius coccineus voucher USNM:FISH:446381 small subunit ribosomal RNA gene, partial sequence; internal transcribed spacer 1, complete sequence",rdna,real,SSU phrasing carrying ITS
"Acanthocybium solandri voucher USNM:FISH:405884 small subunit ribosomal RNA gene, partial sequence; internal transcribed spacer 1, complete sequence",rdna,real,
"Antimora microlepis voucher SIO:Marine Vertebrates:05_103 small subunit ribosomal RNA gene, partial sequence; internal transcribed spacer 1, complete sequence",rdna,real,
"Danio rerio 5.8S ribosomal RNA gene, complete sequence",rdna,authored,5.8S
"Danio rerio 5S ribosomal RNA gene, complete sequence",rdna,authored,5S is nuclear rDNA by design decision
"Homo sapiens 18S ribosomal RNA gene, complete sequence",rdna,authored,rdna even though it is a contaminant; taxonomy decides host
"Aspergillus niger internal transcribed spacer 1, complete sequence",rdna,authored,fungal ITS is still rDNA
"Danio rerio large subunit ribosomal RNA gene, partial sequence",rdna,authored,LSU phrasing
"Danio rerio ribosomal DNA intergenic spacer region",rdna,authored,IGS carrying a ribosomal marker
"Escherichia coli 16S ribosomal RNA gene, partial sequence",other,authored,prokaryotic SSU
"Bacillus subtilis 23S ribosomal RNA gene, partial sequence",other,authored,prokaryotic LSU
"Danio rerio mitochondrion, complete genome",other,authored,organellar
"Danio rerio mitochondrial 16S ribosomal RNA gene, partial sequence",other,authored,organellar veto and prokaryotic-style name
"Zea mays chloroplast large subunit ribosomal RNA gene, complete sequence",other,authored,organellar veto must beat the rDNA phrase
"Danio rerio actin alpha 1 mRNA, complete cds",other,authored,protein coding
"Danio rerio recombination activating gene 1 (rag1), partial cds",other,authored,protein coding
"Danio rerio clone CH211-152C2, complete sequence",other,authored,unrecognised phrasing defaults to other
"",uninformative,authored,empty title says nothing
```

- [ ] **Step 2: Write the failing test**

Create `tests/testthat/test-hit-class.R`:

```r
# Three-state classification of a GenBank defline. This layer answers "is this
# record nuclear rDNA", NOT "is this the host": a fungal or human rDNA record is
# correctly rdna, and catching it as a contaminant is the taxonomy channel's job.

fixtures <- read.csv(testthat::test_path("fixtures", "hit_titles.csv"),
  stringsAsFactors = FALSE, colClasses = "character")

test_that("every fixture title classifies as labelled", {
  got <- RiboPilot:::classify_hit_title(fixtures$title)
  bad <- which(got != fixtures$expected)
  expect_equal(got, fixtures$expected,
    info = paste("mismatched:", paste(fixtures$title[bad], collapse = " | ")))
})

test_that("the corpus actually covers all three classes and both sources", {
  # a corpus that drifted to one class would make the test above vacuous
  expect_setequal(unique(fixtures$expected), RiboPilot:::HIT_CLASSES)
  expect_setequal(unique(fixtures$source), c("real", "authored"))
  expect_gte(sum(fixtures$expected == "uninformative" & fixtures$source == "real"), 26L)
})

test_that("the organellar veto beats an rDNA phrase in the same title", {
  expect_equal(RiboPilot:::classify_hit_title(
    "Zea mays chloroplast large subunit ribosomal RNA gene, complete sequence"), "other")
  expect_equal(RiboPilot:::classify_hit_title(
    "Danio rerio mitochondrial 18S ribosomal RNA gene"), "other")
})

test_that("a prokaryotic record does not match the broadened rRNA gene pattern", {
  expect_equal(RiboPilot:::classify_hit_title(
    "Escherichia coli 16S ribosomal RNA gene, partial sequence"), "other")
  # but a record naming both keeps its eukaryotic reading
  expect_equal(RiboPilot:::classify_hit_title(
    "Synthetic construct 16S and 18S ribosomal RNA genes"), "rdna")
})

test_that("an assembly title naming an rDNA locus is rdna, not uninformative", {
  # precedence: the title says what the locus is, so the container is irrelevant
  expect_equal(RiboPilot:::classify_hit_title(
    "Danio rerio chromosome 5, 18S ribosomal RNA gene, complete sequence"), "rdna")
})

test_that("NA and empty titles are uninformative rather than an error", {
  expect_equal(RiboPilot:::classify_hit_title(NA_character_), "uninformative")
  expect_equal(RiboPilot:::classify_hit_title(""), "uninformative")
  expect_equal(RiboPilot:::classify_hit_title(character(0)), character(0))
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "hit-class")'`
Expected: FAIL, `could not find function "classify_hit_title"`.

- [ ] **Step 4: Broaden the patterns**

In `R/constants.R`, replace the `NUCLEAR_RRNA_TITLE_RE` definition (lines 54-57) with the four patterns below, keeping the existing roxygen prose above it and leaving `ORGANELLAR_TITLE_RE` untouched.

```r
#' Titles naming a nuclear rDNA locus.
#'
#' Broadened from a literal-phrase list because GenBank writes the same locus
#' many ways: "28S ribosomal RNA (ACC), rRNA" carries no "gene", and a
#' chromosome-level host record may name the locus in passing. Safe to include
#' the bare "rRNA gene" construction only because PROKARYOTIC_RRNA_RE is tested
#' first.
#' @noRd
NUCLEAR_RRNA_TITLE_RE <- paste0(
  "internal transcribed spacer|external transcribed spacer|\\bITS[12]?\\b|",
  "\\b(18S|28S|5\\.8S|5S)\\b[^;]{0,40}?(ribosomal RNA|rRNA)|",
  "small subunit ribosomal RNA|large subunit ribosomal RNA|",
  "\\b(SSU|LSU)\\b[^;]{0,20}?rRNA|",
  "ribosomal RNA gene|\\brRNA gene\\b|ribosomal DNA|\\brDNA\\b")

#' Titles naming a prokaryotic ribosomal subunit.
#' @noRd
PROKARYOTIC_RRNA_RE <- "\\b(16S|23S)\\b"

#' Eukaryotic nuclear rDNA markers, used to spare a record that names both.
#' @noRd
EUKARYOTIC_RRNA_MARKER_RE <- paste0(
  "\\b(18S|28S|5\\.8S)\\b|internal transcribed spacer|\\bITS[12]?\\b")

#' Titles describing an assembly rather than a locus.
#'
#' A chromosome record is not evidence against rDNA, it is no evidence either
#' way. Several markers are matched because real deflines contain typos:
#' "unplaced-saffold_25" occurs in GenBank.
#' @noRd
ASSEMBLY_TITLE_RE <- paste0(
  "chromosome|unplaced|scaffold|contig|whole genome shotgun|\\bWGS\\b|haplotype")
```

- [ ] **Step 5: Write the classifier**

Create `R/hit_class.R`:

```r
# Three-state classification of a GenBank defline.
#
# This answers "does this record describe nuclear rDNA". It does NOT answer "is
# this the host": a fungal ITS record and a human 18S record are both correctly
# rdna. Deciding that either is a contaminant is the taxonomy channel's job,
# from rank distance and the identity envelope.
#
# The third state exists because a chromosome-level or whole-genome record is
# not evidence AGAINST rDNA. Counting it as such is what made genuine host
# scaffolds read as "not a nuclear rDNA record".

#' The three hit classes, in no particular order.
#' @noRd
HIT_CLASSES <- c("rdna", "other", "uninformative")

#' Classify GenBank defline(s) into rdna / other / uninformative
#'
#' Precedence is fixed: the organellar veto wins over any rDNA phrase, because a
#' mitochondrial LSU record contains "large subunit ribosomal RNA"; the
#' prokaryotic test comes next so that "16S rRNA gene" cannot match the bare
#' "rRNA gene" construction; then nuclear rDNA phrasing; then assembly-level
#' phrasing; and anything unrecognised stays "other", as it was before.
#'
#' @param title Character vector of defline titles.
#' @return Character vector of the same length, each one of [HIT_CLASSES].
#' @noRd
classify_hit_title <- function(title) {
  t <- as.character(title)
  if (length(t) == 0) return(character(0))
  t[is.na(t)] <- ""
  out <- rep("other", length(t))
  is_rdna <- grepl(NUCLEAR_RRNA_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_prok <- grepl(PROKARYOTIC_RRNA_RE, t, ignore.case = TRUE, perl = TRUE) &
    !grepl(EUKARYOTIC_RRNA_MARKER_RE, t, ignore.case = TRUE, perl = TRUE)
  is_asm <- grepl(ASSEMBLY_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_org <- grepl(ORGANELLAR_TITLE_RE, t, ignore.case = TRUE)
  # assigned low precedence first, so each later line overrides the earlier
  out[is_asm] <- "uninformative"
  out[is_rdna] <- "rdna"
  out[is_prok] <- "other"
  out[is_org] <- "other"
  out[!nzchar(t)] <- "uninformative"
  out
}
```

- [ ] **Step 6: Point the legacy predicate at the new classifier**

`is_nuclear_rrna_title()` (R/blast_genbank.R:13) reads `NUCLEAR_RRNA_TITLE_RE`, which Step 4 just broadened, so leaving it alone would silently change its behaviour: "Escherichia coli 16S ribosomal RNA gene" would start matching the bare "ribosomal RNA gene" construction, which `test-export.R:105` correctly asserts must be FALSE. Redefine it as a thin wrapper so the two cannot diverge:

```r
#' Does this title name a nuclear rDNA record?
#'
#' Thin wrapper over [classify_hit_title()] so the binary predicate and the
#' three-state classifier can never disagree. `uninformative` reads as FALSE
#' here, which is what every existing caller already assumed.
#' @param title Character vector of subject titles.
#' @noRd
is_nuclear_rrna_title <- function(title) {
  classify_hit_title(title) == "rdna"
}
```

Its two callers (R/blast_genbank.R:60 and :205) keep working unchanged, and `blast_hits.is_nuclear_rrna` keeps its existing 1/0 meaning. Do NOT modify those callers in this task; Task 4 rewires them.

- [ ] **Step 7: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "hit-class")'`
Expected: PASS.

If a fixture row fails, fix the PATTERN, not the fixture, unless the fixture's `expected` is genuinely wrong on the spec's own rules; in that case say so explicitly in your report.

- [ ] **Step 8: Check the whole suite and commit**

Run: `Rscript -e 'testthat::test_local()'`
Expected: 0 failures, with `test-export.R:100`'s existing assertions passing UNCHANGED, including that a bacterial 16S record is not nuclear rDNA. If any pre-existing assertion fails, the pattern or the wrapper is wrong; do not edit the assertion.

```bash
git add R/hit_class.R R/constants.R tests/testthat/fixtures/hit_titles.csv tests/testthat/test-hit-class.R
git commit -m "feat(blast): three-state classification of a GenBank defline"
```

---

### Task 2: The bitscore-weighted consensus

**Files:**
- Modify: `R/hit_class.R`
- Test: `tests/testthat/test-hit-class.R`

**Interfaces:**
- Consumes: `HIT_CLASSES` from Task 1.
- Produces: `consensus_hit_class(hit_class, bitscore)` -> `1L`, `0L` or `NA_integer_`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-hit-class.R`:

```r
test_that("consensus reads the informative hits and ignores the rest", {
  expect_equal(RiboPilot:::consensus_hit_class(
    c("rdna", "rdna", "rdna"), c(900, 890, 880)), 1L)
  expect_equal(RiboPilot:::consensus_hit_class(
    c("other", "other", "other"), c(900, 890, 880)), 0L)
  # the reference project's case: every top hit is a chromosome record
  expect_true(is.na(RiboPilot:::consensus_hit_class(
    c("uninformative", "uninformative", "uninformative"), c(900, 890, 880))))
})

test_that("one mislabelled top hit cannot flip the unit", {
  # the outlier holds the TOP bitscore, which is exactly the case hit #1 lost to
  cls <- c("other", "rdna", "rdna", "rdna")
  bs <- c(1000, 990, 985, 980)
  expect_equal(RiboPilot:::consensus_hit_class(cls, bs), 1L)
})

test_that("hits far below the top bitscore do not vote", {
  # 3 rdna near the top, a pile of weak other hits that must not outweigh them
  cls <- c(rep("rdna", 3), rep("other", 20))
  bs <- c(1000, 995, 990, rep(200, 20))
  expect_equal(RiboPilot:::consensus_hit_class(cls, bs), 1L)
})

test_that("a tie resolves to rdna", {
  # a false 0 costs a spurious triage flag; a false 1 only costs its absence
  expect_equal(RiboPilot:::consensus_hit_class(c("rdna", "other"), c(900, 900)), 1L)
})

test_that("fewer than three hits still vote, and empty input is NA", {
  expect_equal(RiboPilot:::consensus_hit_class("rdna", 900), 1L)
  expect_equal(RiboPilot:::consensus_hit_class(c("other", "rdna"), c(900, 100)), 0L)
  expect_true(is.na(RiboPilot:::consensus_hit_class(character(0), numeric(0))))
  expect_true(is.na(RiboPilot:::consensus_hit_class(NA_character_, NA_real_)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "hit-class")'`
Expected: FAIL, `could not find function "consensus_hit_class"`.

- [ ] **Step 3: Implement**

Append to `R/hit_class.R`:

```r
#' Decide a unit's nuclear-rDNA flag from its hit distribution
#'
#' Votes over hits within 10 percent of the top bitscore, taking at least three
#' where that many exist, which is the same defence against one mislabelled
#' GenBank accession that [bitscore_weighted_lca()] uses. Uninformative hits are
#' dropped from the vote rather than counted against rDNA.
#'
#' @param hit_class Character vector of [HIT_CLASSES].
#' @param bitscore Numeric vector, same length.
#' @return `1L` (nuclear rDNA), `0L` (not), or `NA_integer_` when no hit is
#'   informative.
#' @noRd
consensus_hit_class <- function(hit_class, bitscore) {
  # parallel vectors: a mismatch either crashes on a poisoned max() or silently
  # votes on misaligned rows, so fail loudly instead
  stopifnot(length(hit_class) == length(bitscore))
  hc <- as.character(hit_class)
  bs <- suppressWarnings(as.numeric(bitscore))
  ok <- !is.na(hc) & !is.na(bs)
  hc <- hc[ok]; bs <- bs[ok]
  if (length(bs) == 0) return(NA_integer_)
  sel <- bs >= 0.9 * max(bs)
  if (sum(sel) < 3L) {
    sel <- seq_along(bs) %in% order(bs, decreasing = TRUE)[seq_len(min(3L, length(bs)))]
  }
  hc <- hc[sel]; bs <- bs[sel]
  keep <- hc %in% c("rdna", "other")
  if (!any(keep)) return(NA_integer_)
  hc <- hc[keep]; bs <- bs[keep]
  # tie resolves to rdna: a false 0 raises a spurious triage flag, a false 1
  # only fails to raise one
  if (sum(bs[hc == "rdna"]) >= sum(bs[hc == "other"])) 1L else 0L
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "hit-class")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/hit_class.R tests/testthat/test-hit-class.R
git commit -m "feat(blast): bitscore-weighted consensus over a unit's hit classes"
```

---

### Task 3: Schema v13 and the offline backfill

**Files:**
- Modify: `R/migrate_db.R` (`RIBOPILOT_SCHEMA_VERSION` at line 8, `BLAST_HITS_DDL` at line 18, a new v13 block after the v12 block, and a new `backfill_hit_class()` next to `backfill_coverage_median()`)
- Test: `tests/testthat/test-migrate.R`

**Interfaces:**
- Consumes: `classify_hit_title()` and `consensus_hit_class()` from Tasks 1 and 2.
- Produces: `blast_hits.hit_class TEXT`; `backfill_hit_class(con)` -> invisibly the number of `blast_hits` rows classified.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-migrate.R`:

```r
test_that("v13 adds hit_class and backfills it from stored titles", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  # a v12-shaped blast_hits, i.e. with every column EXCEPT hit_class
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, hsp INT, accession TEXT, title TEXT, species TEXT, staxid INT,
    lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, aln_len INT, qstart INT, qend INT, sstart INT,
    send INT, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT,
    is_nuclear_rrna INT)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,2,0)")
  ins <- "INSERT INTO blast_hits (ID,path,scaffold,hit_rank,hsp,title,bitscore,is_nuclear_rrna)
          VALUES (?,?,?,?,1,?,?,?)"
  # unit 1: all chromosome records, the case that used to score 0
  for (i in 1:3) DBI::dbExecute(con, ins, params = list("S1", 1L, 1L, i,
    sprintf("Scomber japonicus DNA, primary haplotype, chromosome %d, sequence", i), 900 - i, 0L))
  # unit 2: genuine rDNA
  for (i in 1:3) DBI::dbExecute(con, ins, params = list("S1", 1L, 2L, i,
    "Danio rerio 18S ribosomal RNA gene, complete sequence", 900 - i, 1L))
  DBI::dbExecute(con, "PRAGMA user_version = 12")

  migrate_db(con)

  expect_true("hit_class" %in% DBI::dbListFields(con, "blast_hits"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 13L)
  cls <- DBI::dbGetQuery(con, "SELECT scaffold, hit_class FROM blast_hits ORDER BY scaffold, hit_rank")
  expect_equal(unique(cls$hit_class[cls$scaffold == 1]), "uninformative")
  expect_equal(unique(cls$hit_class[cls$scaffold == 2]), "rdna")

  # the repair: unit 1 goes from a false 0 to NA, unit 2 to 1
  scr <- DBI::dbGetQuery(con, "SELECT scaffold, is_nuclear_rrna FROM blast_screen ORDER BY scaffold")
  expect_true(is.na(scr$is_nuclear_rrna[scr$scaffold == 1]))
  expect_equal(scr$is_nuclear_rrna[scr$scaffold == 2], 1L)
})

test_that("the backfill leaves units with no hits alone and is idempotent", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, title TEXT, bitscore REAL, hit_class TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT,
    is_nuclear_rrna INT)")
  # a screened unit with NO hits rows must keep whatever it had
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,9,1)")
  RiboPilot:::backfill_hit_class(con)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT is_nuclear_rrna FROM blast_screen WHERE scaffold=9")$is_nuclear_rrna, 1L)

  DBI::dbExecute(con, "INSERT INTO blast_hits VALUES('S1',1,1,1,
    'Danio rerio 18S ribosomal RNA gene, complete sequence', 900, NULL)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,0)")
  n1 <- RiboPilot:::backfill_hit_class(con)
  expect_equal(n1, 1L)
  n2 <- RiboPilot:::backfill_hit_class(con)
  expect_equal(n2, 0L)   # nothing left to classify
  expect_equal(DBI::dbGetQuery(con,
    "SELECT is_nuclear_rrna FROM blast_screen WHERE scaffold=1")$is_nuclear_rrna, 1L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: FAIL, `hit_class` not in `dbListFields`.

- [ ] **Step 3: Implement**

In `R/migrate_db.R`, set line 8 to:

```r
RIBOPILOT_SCHEMA_VERSION <- 13L
```

Add `hit_class TEXT,` to `BLAST_HITS_DDL` immediately after the `is_nuclear_rrna INTEGER,` line, so a fresh database and a migrated one agree.

Add after the v12 block, before the `user_version` write:

```r
  # --- schema v13: three-state hit classification ---
  # is_nuclear_rrna was a literal-phrase regex over ONE hit's title, so a
  # chromosome-level host record read as "not nuclear rDNA". Classification is a
  # pure function of the stored title, so history is repairable offline: the
  # backfill costs no remote BLAST, which matters because a real run takes hours.
  add_col("blast_hits", "hit_class", "TEXT")
  if (has_tbl("blast_hits")) backfill_hit_class(con)
```

Add next to `backfill_coverage_median()`:

```r
#' Classify stored BLAST hit titles and recompute the per-unit screen flag
#'
#' Runs offline over titles already in the database, so an existing project is
#' repaired without repeating its remote BLAST. Units with no `blast_hits` rows
#' are left exactly as they are: absence of hits is not evidence.
#'
#' @param con Open project DB connection.
#' @return (invisibly) the number of `blast_hits` rows classified.
#' @noRd
backfill_hit_class <- function(con) {
  if (!("hit_class" %in% DBI::dbListFields(con, "blast_hits"))) {
    return(invisible(0L))
  }
  todo <- DBI::dbGetQuery(con,
    "SELECT rowid AS rid, title FROM blast_hits WHERE hit_class IS NULL")
  n <- nrow(todo)
  if (n > 0) {
    todo$hit_class <- classify_hit_title(todo$title)
    DBI::dbWithTransaction(con, {
      up <- DBI::dbSendStatement(con, "UPDATE blast_hits SET hit_class = ? WHERE rowid = ?")
      on.exit(DBI::dbClearResult(up), add = TRUE)
      DBI::dbBind(up, list(todo$hit_class, todo$rid))
    })
  }
  # recompute the screen flag for every unit that HAS hits
  units <- DBI::dbGetQuery(con,
    "SELECT DISTINCT ID, path, scaffold FROM blast_hits")
  if (nrow(units) > 0 && "blast_screen" %in% DBI::dbListTables(con)) {
    DBI::dbWithTransaction(con, {
      for (i in seq_len(nrow(units))) {
        h <- DBI::dbGetQuery(con,
          "SELECT hit_class, bitscore FROM blast_hits
           WHERE ID = ? AND path = ? AND scaffold = ?",
          params = list(units$ID[i], units$path[i], units$scaffold[i]))
        DBI::dbExecute(con,
          "UPDATE blast_screen SET is_nuclear_rrna = ?
           WHERE ID = ? AND path = ? AND scaffold = ?",
          params = list(consensus_hit_class(h$hit_class, h$bitscore),
            units$ID[i], units$path[i], units$scaffold[i]))
      }
    })
  }
  invisible(n)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: PASS, including the pre-existing "migrate_db is idempotent" test.

- [ ] **Step 5: Commit**

```bash
git add R/migrate_db.R tests/testthat/test-migrate.R
git commit -m "feat(db): hit_class column with an offline backfill (v13)"
```

---

### Task 4: Wire the parse path onto the new classification

**Files:**
- Modify: `R/blast_genbank.R` (line 60, the `hits` frame at 100-108, the reference filter at 117, and `write_blast_screen()` at 195-211)
- Test: `tests/testthat/test-blast-parse.R`

**Interfaces:**
- Consumes: `classify_hit_title()`, `consensus_hit_class()`.
- Produces: `write_blast_screen(con, id, path, scaffold, top, hits = NULL, ts = ...)` - `hits` is the unit's full hit frame, carrying at least `hit_class` and `bitscore`; when NULL the flag falls back to the top row's own class, which keeps the function usable on its own.

- [ ] **Step 1: Write the failing test**

The file already provides `mk_blast_db()` (line 11) and `write_hits()` (line 27); use both rather than building your own. `mk_blast_db()` creates `blast_hits` from `RiboPilot:::BLAST_HITS_DDL`, so it picks up `hit_class` automatically once Task 3 adds it to that constant. The tab-separated column order `parse_blast_hits_file()` expects is, verbatim: accession, title, pident, qcovs, evalue, staxid, bitscore, aln_len, qstart, qend, sstart, send.

```r
test_that("parse writes hit_class and takes the screen flag from consensus", {
  con <- mk_blast_db()
  on.exit(DBI::dbDisconnect(con))
  # three chromosome-level host hits: informative about nothing
  f <- write_hits(c(
    "CM001.1\tScomber japonicus DNA, primary haplotype, chromosome 1, sequence\t99\t90\t0.0\t9999\t900\t1500\t1\t1500\t1\t1500",
    "CM002.1\tScomber japonicus DNA, primary haplotype, chromosome 4, sequence\t98\t90\t0.0\t9999\t890\t1500\t1\t1500\t1\t1500",
    "CM003.1\tScomber japonicus DNA, alternate haplotype, chromosome 1, sequence\t97\t90\t0.0\t9999\t880\t1500\t1\t1500\t1\t1500"))
  RiboPilot:::parse_blast_candidates(con, "S1", 1, 1, hits_file = f)

  cls <- DBI::dbGetQuery(con, "SELECT DISTINCT hit_class FROM blast_hits")$hit_class
  expect_equal(cls, "uninformative")
  flag <- DBI::dbGetQuery(con, "SELECT is_nuclear_rrna FROM blast_screen")$is_nuclear_rrna
  expect_true(is.na(flag))
  # and an assembly-level record is not offered as a reference
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 0L)
})
```

If anything in that fixture does not line up with the file as you find it, adapt to the file rather than inventing a new helper, and say so in your report.

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "blast-parse")'`
Expected: FAIL, `hit_class` is not a column of `blast_hits` in the fixture, or the flag is 0 rather than NA.

- [ ] **Step 3: Implement**

In `R/blast_genbank.R`:

At line 60, keep the existing `is_nuclear_rrna` assignment for back-compatibility and add the class beside it:

```r
  h$is_nuclear_rrna <- as.integer(is_nuclear_rrna_title(h$title))
  h$hit_class <- classify_hit_title(h$title)
```

Add `hit_class = h$hit_class,` to the `hits` data.frame (line 100-107) next to `is_nuclear_rrna`.

Change the screen call at line 114 to pass the full frame:

```r
  write_blast_screen(con, id, path, scaffold, subj[1, ], hits = h, ts = ts)
```

Change the reference filter at line 117 from `subj$is_nuclear_rrna == 1L` to:

```r
  # assembly-level records are poor per-element references, so the panel stays
  # rdna-only rather than admitting uninformative hits
  ref <- subj[subj$hit_class == "rdna", , drop = FALSE]
```

In `write_blast_screen()`, add the `hits` argument and take the flag from the consensus:

```r
write_blast_screen <- function(con, id, path, scaffold, top, hits = NULL,
                               ts = as.integer(Sys.time())) {
```

and replace the `is_nuclear_rrna = ...` element of the row it writes with:

```r
    is_nuclear_rrna = if (!is.null(hits)) {
      consensus_hit_class(hits$hit_class, hits$bitscore)
    } else {
      consensus_hit_class(top$hit_class %||% classify_hit_title(top$title), top$bitscore)
    },
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "blast")'`
Expected: PASS.

Then the whole suite: `Rscript -e 'testthat::test_local()'`
Expected: 0 failures. Existing blast-parse tests assert `is_nuclear_rrna` on rDNA fixtures, which consensus preserves; if any now fails, report which and why before changing it.

- [ ] **Step 5: Commit**

```bash
git add R/blast_genbank.R tests/testthat/test-blast-parse.R
git commit -m "feat(blast): classify every hit and screen from the consensus"
```

---

### Task 5: The triage reason string, docs, and final verification

**Files:**
- Modify: `R/triage_rules.R` (the `no_screen` entry), `vignettes/Curation-Flags.Rmd`
- Test: `tests/testthat/test-triage-rules.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: no signature changes.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-triage-rules.R`:

```r
test_that("the no_screen reason covers an unclassifiable hit, not just a missing run", {
  # is_nuclear_rrna is now NA both when the screen did not run and when it ran
  # and every hit was assembly-level, so the reason must not assert the former
  d <- RiboPilot:::triage_rules(us(list(feats = 5, nuc = NA)))
  expect_equal(d$new_rule, "no_screen")
  expect_false(grepl("has not run for this scaffold", d$new_reason, fixed = TRUE))
  expect_match(d$new_reason, "nuclear rDNA")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage-rules")'`
Expected: FAIL, the reason still reads "the GenBank screen has not run for this scaffold".

- [ ] **Step 3: Implement**

In `R/triage_rules.R`, in the `no_screen` rule definition, change only the returned reason string to:

```r
      "no GenBank hit says whether this scaffold is nuclear rDNA"
```

Leave the rule's `key`, `label`, `meaning`, `actions`, `default` and its position in the registry exactly as they are. Its `meaning` ("No usable core_nt hit is recorded, so the scaffold cannot be judged against GenBank") is already correct for both cases.

- [ ] **Step 4: Update the vignette**

In `vignettes/Curation-Flags.Rmd`, in the triage section, add one sentence in the file's existing voice explaining that a scaffold whose only GenBank matches are whole-genome or chromosome-level assemblies is reported as unjudged rather than as not-rDNA, because such a record does not say what the aligned region is. Read the surrounding prose and match it; do not paste from the design doc.

- [ ] **Step 5: Run the full suite and check the package**

```bash
Rscript -e 'testthat::test_local()'
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'
```

Expected: 0 test failures; `devtools::check` at 0 errors and 0 warnings. Three NOTEs are pre-existing (crayon/dbplyr imports, `new_db.Rd` line width, hidden files); report the exact counts and flag any NOTE not in that list.

- [ ] **Step 6: Commit**

```bash
git add R/triage_rules.R vignettes/ tests/testthat/test-triage-rules.R man/
git commit -m "fix(triage): say no hit is classifiable, not that the screen never ran"
```

---

## Self-Review

**Spec coverage:**
- Spec section 1 (what this is and is not) -> Task 1 file header comment and fixture rows for fungal and human rDNA.
- Section 3 classification and precedence -> Task 1 Steps 4-5.
- Section 4 storage and migration -> Task 3.
- Section 5 consensus -> Task 2, wired in Task 4.
- Section 6 triage consequence -> Task 5.
- Section 7 fixture corpus -> Task 1 Step 1.
- Section 8 testing -> distributed across Tasks 1-4; the reference-candidate exclusion is Task 4 Step 1; migration idempotency is Task 3 Step 1.
- Section 9 out of scope -> nothing in this plan touches the envelope, the `screen` table, or any Nextflow process.

**Placeholders:** none. Every code step carries its code, and both places where a fixture may not match the existing helper (Task 4 Step 1) say explicitly to read the authoritative source and report rather than guess.

**Type consistency:** `classify_hit_title()` returns a character vector in every use. `consensus_hit_class()` returns `integer(1)` or `NA_integer_` in every use. `hit_class` is `TEXT` in the DDL and character in R throughout. `backfill_hit_class()` returns an integer count, asserted in Task 3.

**Assumptions verified while writing this plan, so no task depends on a guess:**
- `parse_blast_hits_file()` (R/blast_genbank.R:147) reads the columns in the order used by Task 4's fixture: accession, title, pident, qcovs, evalue, staxid, bitscore, aln_len, qstart, qend, sstart, send.
- `tests/testthat/test-blast-parse.R` already defines `mk_blast_db()` (line 11) and `write_hits()` (line 27), and `mk_blast_db()` builds `blast_hits` from `RiboPilot:::BLAST_HITS_DDL`, so Task 3's DDL change reaches the fixture with no test edit.
- `parse_blast_candidates()` calls `migrate_db(con)` at line 48 inside a `tryCatch`, so a project database reaching the parse path is migrated to v13 before any write.
- `ORGANELLAR_TITLE_RE` (R/constants.R:61) is unchanged by this plan and remains the hard veto.
