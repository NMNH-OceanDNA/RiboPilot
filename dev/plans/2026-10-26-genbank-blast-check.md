# GenBank BLAST Check Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store, per contig, the two BLAST flags GenBank's rRNA/ITS curation raises (top-hit query coverage under 75%, more than one HSP on the top hit) and show them at export so curators are ready for the reply.

**Architecture:** A new R/genbank_check.R runs two default megablasts of the contig (ITS set; SSU plus LSU set) and stores one text value, `genbank_check`, on `annotate` (pipeline) or `review_units` (recheck), exactly as the stock ITSx check does with `itsx_warnings`. Export, the Checks column, and the app tables read it through `review_units_v`. No new Nextflow process: the check runs inside the existing `feature_blast` process and inside `recheck_unit()`.

**Tech Stack:** R package (DBI/RSQLite, testthat, reactable), blastn 2.17 (default task megablast), Nextflow bin scripts, shell release scripts under dev/refdb.

**Spec:** dev/specs/2026-10-26-genbank-blast-check-design.md

## Global Constraints

- Branch `feat/genbank-blast-check` off devel. Never push. No Claude attribution in commits. Short commit messages.
- ASCII only in R code. Minimal comments. No em dashes anywhere.
- Schema change = edit `new_db()` (R/init_db.R, R/db_schema.R) AND append step 6 to `db_migrations()` (version 6 to 7). `tests/testthat/fixtures/schema_v1.sql` must not change. A step that changes a view drops it first (R/db_migrate.R header).
- Thresholds: `GENBANK_CHECK_MIN_QCOVS <- 75`, `GENBANK_CHECK_MAX_HSPS <- 1L`. blastn: default task, `-evalue 0.05 -max_target_seqs 100`, no `-max_hsps`, outfmt exactly `RECHECK_BLAST$outfmt` (R/recheck.R) so `parse_blast_hits_file()` reads it.
- Database names: `ITS_eukaryote_sequences`, `SSU_eukaryote_rRNA`, `LSU_eukaryote_rRNA`, all looked up in one `db_dir`. `RP_DECOY_DBS` is not changed. `RP_DB_VERSION` is not changed (no release in this work).
- Stored text grammar (one TEXT value):
  - clean: `""`
  - flags joined by `"; "`, in this order: `genbank_no_hit(ITS)`, `genbank_low_coverage(ITS <q>%)`, `genbank_multi_hsp(ITS <n> HSPs: <s>-<e>, <s>-<e>)`, `genbank_no_hit(rRNA)`, `genbank_multi_hsp(rRNA <n> HSPs: ...)`; `<q>` is `qcovs` as an integer; spans are the top hit's HSPs sorted by `qstart`.
  - when any flag is present and at least one search hit, append `" top hits: ITS <acc> <species> <q>%, rRNA <acc> <species> <q>%"` (a search with no hit contributes `ITS none` / `rRNA none`).
  - failure: `not_run(<reason>)`. Reasons used: `no local BLAST`, `no <dbname> database`, `no sequence`, `BLAST failed: <message>`.
  - rRNA coverage is never flagged.
- Amber only: nothing in this feature feeds `has_screen_risk()`, Lock, or the Export gate.
- User-facing name: "GenBank BLAST check"; column `genbank_check` everywhere (DB, CSV, app, report).

## Review Focus

1. An empty or NULL contig sequence must store `not_run(no sequence)`, not crash blastn (Task 2 test).
2. blastn exiting non-zero on one of the two searches must store `not_run(BLAST failed: ...)` and leave no partial text (Task 2 test).
3. A top hit whose HSPs come back in bitscore order must report spans in query order (Task 2 test).
4. A `not_run(...)` value must count as "not checked" and never as flagged in `unit_checks()`, the export counts, or the CSV (Task 4 test).
5. A review copy made after the pipeline ran must carry the original's `genbank_check`, and an edit must clear it so the recheck rewrites it (Task 1 test).

---

### Task 1: Schema column and migration step 6

**Files:**
- Modify: `R/init_db.R:469-477` (annotate DDL, add `genbank_check TEXT` as the LAST column, so new_db() and the ALTER TABLE in step 6 agree on order)
- Modify: `R/db_schema.R:372-395` (`REVIEW_UNITS_DDL`, add `genbank_check TEXT` as the LAST column), `R/db_schema.R:458-476` (`review_units_v`: add `x.genbank_check` after `x.itsx_scores` in the copy branch and `an.genbank_check` after `an.itsx_scores` in the original branch)
- Modify: `R/db_migrate.R:10-36` (append step 6)
- Modify: `R/review.R:110-118` (`review_invalidate`: add `genbank_check = NULL`)
- Test: `tests/testthat/test-review-schema.R` (append), `tests/testthat/test-backwards-compat.R` (existing parity test must pass unchanged)

**Interfaces:**
- Produces: column `genbank_check TEXT` on `annotate`, `review_units`, exposed by `review_units_v`; `current_schema_version()` returns 7.

- [ ] **Step 1: Write the failing test** (append to test-review-schema.R)

```r
test_that("genbank_check follows itsx_warnings through copy, view and invalidate", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, genbank_check) VALUES ('S1',1,1,'genbank_low_coverage(ITS 70%)')")
  expect_equal(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units_v")$genbank_check, "genbank_low_coverage(ITS 70%)")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units")$genbank_check, "genbank_low_coverage(ITS 70%)")
  RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units_v")$genbank_check))
  expect_equal(RiboPilot:::current_schema_version(), 7L)
})
```

- [ ] **Step 2: Run it, expect failure** `Rscript -e 'devtools::test(filter="review-schema")'` fails with "no such column: genbank_check".

- [ ] **Step 3: Add the column** to the three DDL places and `review_invalidate`. Step 6 in `db_migrations()`:

```r
list(note = "GenBank BLAST check (annotate.genbank_check, review_units.genbank_check)",
  fn = function(con) {
    DBI::dbExecute(con, "DROP VIEW IF EXISTS review_units_v")
    DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN genbank_check TEXT")
    DBI::dbExecute(con, "ALTER TABLE review_units ADD COLUMN genbank_check TEXT")
    DBI::dbExecute(con, REVIEW_VIEWS_DDL[["review_units_v"]])
  })
```

- [ ] **Step 4: Run** `Rscript -e 'devtools::test(filter="review-schema|backwards-compat|review-copy|review-state")'`. Expected: all pass, including the v1-plus-steps parity test.

- [ ] **Step 5: Commit** `git add -A R tests && git commit -m "schema: genbank_check column, migration step 6"`

---

### Task 2: The check (R/genbank_check.R) with the SUB16547976 fixture

**Files:**
- Create: `R/genbank_check.R`
- Create: `tests/testthat/fixtures/genbank_check/` by copying the four files from `/tmp/claude-1000/-home-dmacguig-Documents-GitHub-RiboPilot/7e6d4cd3-f648-4d90-ac10-f6edc95e8466/scratchpad/fixture/` (`SUB16547976.fasta`, `SUB16547976_ITS_eukaryote_sequences.tsv`, `SUB16547976_SSU_LSU_eukaryote_rRNA.tsv`, `SUB16547976_ncbi_flags.tsv`). The two .tsv files are blastn output in `RECHECK_BLAST$outfmt`, all ten queries in one file, `qseqid` in column 14. The flags file is `seqid, low_coverage, hsp_problem` (0/1) transcribed from NCBI's 2026-10-15 reply.
- Create: `tests/testthat/test-genbank-check.R`

**Interfaces:**
- Produces:
  - `GENBANK_CHECK_DBS <- c(ITS = "ITS_eukaryote_sequences", SSU = "SSU_eukaryote_rRNA", LSU = "LSU_eukaryote_rRNA")`
  - `genbank_check_calls(its_hits, rrna_hits)` -> `character(1)` stored text. Each argument is a `parse_blast_hits_file()` frame for ONE query (possibly 0 rows). Top hit = first `accession` in frame order; its HSPs = all rows with that accession. Species via `extract_species(title)`.
  - `genbank_check_display(text)` -> `character(length(text))` short form for the export window: `genbank_low_coverage(ITS 70%)` -> `low coverage ITS 70%`; `genbank_multi_hsp(ITS 2 HSPs: ...)` -> `ITS 2 HSPs`; `genbank_multi_hsp(rRNA ...)` -> `rRNA n HSPs`; `genbank_no_hit(X)` -> `no X hit`; joined with `", "`; the `top hits:` tail dropped; `""` and `not_run(...)` returned unchanged.
  - `genbank_blast_run(query, dbarg, out, threads = 1L, blastn = "blastn")` -> hits frame; stops on non-zero exit (mockable, like `feature_blast_run`).
  - `genbank_check_unit(con, id, path = 1L, scaffold = 1L, db_dir = NULL, workdir = tempfile("genbank_check_"), threads = 1L, review = FALSE, blastn = "blastn")` -> invisible stored text. Sequence from `assemblies` (review FALSE) or `review_units_v` (review TRUE). Writes `annotate.genbank_check` or `review_units.genbank_check`.
  - `genbank_store(con, id, path, scaffold, text, review = FALSE)` -> invisible text (the one UPDATE; a failed write on `annotate` is silent, on `review_units` raises, as `itsx_store_warnings`).
  - `genbank_check_db_dir(dbarg)` -> `dirname()` of the first whitespace-separated token of a `rp_dbarg.R` string, or NULL for NULL/empty.

- [ ] **Step 1: Write the failing tests**

```r
gc_fix <- function(f) test_path("fixtures", "genbank_check", f)
gc_hits <- function(f) RiboPilot:::parse_blast_hits_file(gc_fix(f))
gc_split <- function(h) split(h, h$qseqid)

test_that("the SUB16547976 fixture reproduces NCBI's two flag lists exactly", {
  its <- gc_split(gc_hits("SUB16547976_ITS_eukaryote_sequences.tsv"))
  rr <- gc_split(gc_hits("SUB16547976_SSU_LSU_eukaryote_rRNA.tsv"))
  want <- read.delim(gc_fix("SUB16547976_ncbi_flags.tsv"), stringsAsFactors = FALSE)
  got <- vapply(want$seqid, function(s) RiboPilot:::genbank_check_calls(its[[s]], rr[[s]]), character(1))
  expect_equal(unname(grepl("genbank_low_coverage(ITS", got, fixed = TRUE)), want$low_coverage == 1)
  expect_equal(unname(grepl("genbank_multi_hsp(ITS", got, fixed = TRUE)), want$hsp_problem == 1)
  expect_false(any(grepl("rRNA", got[want$low_coverage == 0 & want$hsp_problem == 0])))
  expect_equal(unname(got[c("07-01-61_S93", "07-02-111_S335")]), c("", ""))
  expect_match(got[["07-01-65_S97"]], "^genbank_low_coverage\\(ITS 70%\\); genbank_multi_hsp\\(ITS 2 HSPs: 1-1855, 2274-2457\\) top hits: ITS OP250806.1 Datnioides microlepis 70%, rRNA MG895428.1 Brachypleura novaezeelandiae 63%$")
})

test_that("genbank_check_calls covers the grammar", {
  mk <- function(acc, spans, qcovs, title = "Danio rerio SSU ITS1") data.frame(accession = acc, title = title,
    pident = 99, qcovs = qcovs, evalue = 0, staxid = 1L, bitscore = rev(seq_along(spans)) * 100,
    aln_len = 10L, qstart = spans, qend = spans + 9L, sstart = 1L, send = 10L, slen = 100L, qseqid = "q",
    stringsAsFactors = FALSE)
  none <- data.frame()
  expect_equal(RiboPilot:::genbank_check_calls(mk("A1", 1L, 100), mk("B1", 1L, 60)), "")
  expect_equal(RiboPilot:::genbank_check_calls(none, none), "genbank_no_hit(ITS); genbank_no_hit(rRNA)")
  expect_match(RiboPilot:::genbank_check_calls(none, mk("B1", 1L, 60)), "^genbank_no_hit\\(ITS\\) top hits: ITS none, rRNA B1 Danio rerio 60%$")
  # spans in query order even when HSP rows arrive in bitscore order
  expect_match(RiboPilot:::genbank_check_calls(mk("A1", c(2000L, 1L), 80), mk("B1", 1L, 60)),
    "genbank_multi_hsp\\(ITS 2 HSPs: 1-10, 2000-2009\\)")
  expect_match(RiboPilot:::genbank_check_calls(mk("A1", 1L, 100), mk("B1", c(1L, 50L), 60)),
    "^genbank_multi_hsp\\(rRNA 2 HSPs: 1-10, 50-59\\) top hits")
  expect_equal(RiboPilot:::genbank_check_display(c("", "not_run(x)",
    "genbank_low_coverage(ITS 70%); genbank_multi_hsp(ITS 2 HSPs: 1-2, 3-4); genbank_no_hit(rRNA) top hits: ITS A 1 2 70%, rRNA none")),
    c("", "not_run(x)", "low coverage ITS 70%, ITS 2 HSPs, no rRNA hit"))
  expect_equal(RiboPilot:::genbank_check_db_dir("/x/2026/ribopilot_rrna /x/2026/16S_ribosomal_RNA"), "/x/2026")
  expect_null(RiboPilot:::genbank_check_db_dir(""))
})

test_that("genbank_check_unit stores text, not_run reasons, and the review copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  fb_seed(con)  # from test-feature-blast.R helpers: S1 1 1 with a sequence and annotate row
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 1)")
  d <- withr::local_tempdir()
  for (n in RiboPilot:::GENBANK_CHECK_DBS) file.create(file.path(d, paste0(n, ".nin")))
  gc <- function() DBI::dbGetQuery(con, "SELECT genbank_check FROM annotate WHERE ID = 'S1'")$genbank_check
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = NULL)
  expect_equal(gc(), "not_run(no local BLAST)")
  unlink(file.path(d, "ITS_eukaryote_sequences.nin"))
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(no ITS_eukaryote_sequences database)")
  file.create(file.path(d, "ITS_eukaryote_sequences.nin"))
  testthat::local_mocked_bindings(genbank_blast_run = function(...) stop("boom"), .package = "RiboPilot")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(BLAST failed: boom)")
  its <- gc_hits("SUB16547976_ITS_eukaryote_sequences.tsv"); its <- its[its$qseqid == "07-01-65_S97", ]
  rr <- gc_hits("SUB16547976_SSU_LSU_eukaryote_rRNA.tsv"); rr <- rr[rr$qseqid == "07-01-65_S97", ]
  testthat::local_mocked_bindings(genbank_blast_run = function(query, dbarg, ...)
    if (grepl("ITS_eukaryote", dbarg)) its else rr, .package = "RiboPilot")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_match(gc(), "^genbank_low_coverage\\(ITS 70%\\)")
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = NULL")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(no sequence)")
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = 'ACGT'")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d, review = TRUE)
  expect_match(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units")$genbank_check, "^genbank_low_coverage")
  expect_equal(gc(), "not_run(no sequence)")
})
```

`fb_seed()` is file-local at tests/testthat/test-feature-blast.R:121; move it to `tests/testthat/helper-chimera.R` so both files share it.

- [ ] **Step 2: Run** `Rscript -e 'devtools::test(filter="genbank-check")'`. Expected: fails, functions not found.

- [ ] **Step 3: Implement R/genbank_check.R** with the signatures above. `genbank_check_unit` order: resolve sequence (`not_run(no sequence)` when NA/empty); `db_dir` NULL/empty -> `not_run(no local BLAST)`; each of the three prefixes must pass `blast_db_exists()` else `not_run(no <name> database)`; write `query.fa`; run ITS then rRNA (`dbarg` = the two prefixes joined by a space) inside one `tryCatch` -> `not_run(BLAST failed: <msg>)`; else `genbank_check_calls()`; `genbank_store()`. Add `extract_species` reuse; no new helpers beyond those listed.

- [ ] **Step 4: Run** the same filter. Expected: PASS. Also `Rscript -e 'devtools::document()'` if any roxygen `@export` is added (prefer `@noRd`; nothing here needs exporting).

- [ ] **Step 5: Commit** `git add R/genbank_check.R tests && git commit -m "feat: GenBank BLAST check with SUB16547976 fixture"`

---

### Task 3: Run it in the pipeline and the recheck

**Files:**
- Modify: `inst/nextflow/bin/rp_feature_blast.R:12-16`
- Modify: `R/recheck.R:173-175` (after `feature_blast_unit`)
- Test: `tests/testthat/test-recheck.R` (append), `tests/testthat/test-nf-bin-contract.R` (only if it lists rp_feature_blast.R's calls)

**Interfaces:**
- Consumes: `genbank_check_unit`, `genbank_check_db_dir` (Task 2).

- [ ] **Step 1: Write the failing test** (test-recheck.R; use the file's existing recheck_unit setup helper)

```r
# beside "an ITSx check that throws is stored as not_run and the contig still finishes" (:369),
# using that file's rc_project(), rc_rev(), rc_tools(), rc_copy() helpers
test_that("recheck_unit stores the GenBank check on the copy and survives its failure", {
  con <- rc_project()
  seen <- NULL
  testthat::local_mocked_bindings(genbank_check_unit = function(con, id, path, scaffold, db_dir = NULL, ..., review = FALSE) {
    seen <<- list(db_dir = db_dir, review = review)
    RiboPilot:::genbank_store(con, id, path, scaffold, "genbank_no_hit(ITS); genbank_no_hit(rRNA)", review = review)
  }, .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", rc_tools(), dbarg = "/p/v1/ribopilot_rrna /p/v1/16S_ribosomal_RNA"), "done")
  expect_equal(seen, list(db_dir = "/p/v1", review = TRUE))
  expect_equal(rc_copy(con)$genbank_check, "genbank_no_hit(ITS); genbank_no_hit(rRNA)")
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 1 WHERE", RC_UNIT))
  testthat::local_mocked_bindings(genbank_check_unit = function(...) stop("boom"), .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", rc_tools()), "done")
  expect_equal(rc_copy(con)$genbank_check, "not_run(GenBank check failed: boom)")
})
```

- [ ] **Step 2: Run** `Rscript -e 'devtools::test(filter="^recheck")'`. Expected: new test fails.

- [ ] **Step 3: Wire it.** In `recheck_unit()` after the feature BLAST call:

```r
tryCatch(genbank_check_unit(con, id, path, scaffold, db_dir = genbank_check_db_dir(dbarg),
  threads = threads, workdir = file.path(workdir, "genbank_check"), review = TRUE),
  error = function(e) genbank_store(con, id, path, scaffold,
    paste0("not_run(GenBank check failed: ", conditionMessage(e), ")"), review = TRUE))
```

In rp_feature_blast.R, after `feature_blast_unit()` and before `review_follow_original()`, the same call with `review = FALSE`, `db_dir = RiboPilot:::genbank_check_db_dir(opt(5))`, `workdir = file.path(getwd(), "genbank_check_work")`, errors caught the same way (store `not_run(GenBank check failed: ...)` on annotate). Print one line `genbank_check: <id> <path> <scaffold> - <display text or not_run>`.

- [ ] **Step 4: Run** `Rscript -e 'devtools::test(filter="recheck|nf-bin-contract|nf-main-run|nf-recheck-run")'`. Expected: PASS.

- [ ] **Step 5: Commit** `git add -A R inst/nextflow tests && git commit -m "feat: run GenBank check in feature_blast step and recheck"`

---

### Task 4: Export surface (CSV, window, README, HTML report) and the Checks column

**Files:**
- Modify: `R/export.R`: `unit_checks()` (:259-312), `submission_readme()` (after the ITSx section ending :408), `export_summary_df()` (:459-513), `export_report_html()` (:712)
- Modify: `R/app_modals.R:869-899` (`export_check_summary_ui`)
- Modify: `inst/export_report.Rmd:37-38`
- Test: `tests/testthat/test-export.R` (append), `tests/testthat/test-app-export.R` (append)

**Interfaces:**
- Consumes: `genbank_check_display` (Task 2).
- Produces: `export_summary_df()` column `genbank_check` after `itsx_check` (NULL stored -> `not_run(not checked)`); `unit_checks()` entry `"GenBank BLAST"` (warning tone, listed after `"ITSx"`); `export_check_summary_ui(smry, not_done)` adds one `tags$p` after the ITSx line.

- [ ] **Step 1: Write the failing tests**

```r
# test-export.R, inside or beside the unit_checks test: scaffold 1 flagged, scaffold 2 not_run
test_that("the GenBank BLAST check is amber in unit_checks and present in the summary", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  # seed two kept units as the unit_checks test does, screen on_target, itsx_warnings ''
  DBI::dbExecute(con, "UPDATE annotate SET genbank_check = 'genbank_low_coverage(ITS 70%) top hits: ITS A 1 2 70%, rRNA none' WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET genbank_check = 'not_run(no local BLAST)' WHERE scaffold = 2")
  u <- RiboPilot:::unit_checks(con, "unit")
  expect_equal(u$checks, c("GenBank BLAST", "ok"))
  expect_equal(u$checks_tone, c("warning", "success"))
  smry <- RiboPilot:::export_summary_df(con, cbind(u[1:3], seqid = c("a", "b")))
  expect_equal(smry$genbank_check, c("genbank_low_coverage(ITS 70%) top hits: ITS A 1 2 70%, rRNA none", "not_run(no local BLAST)"))
  expect_false(RiboPilot:::has_screen_risk(smry$submission_risk)[1])
  expect_match(RiboPilot:::submission_readme("g"), "GENBANK BLAST CHECK")
  expect_match(RiboPilot:::submission_readme("g"), "HSP_problem", fixed = TRUE)
})
# test-app-export.R
test_that("export_check_summary_ui counts and names GenBank BLAST flags", {
  smry <- data.frame(seqid = c("a", "b", "c"), submission_risk = "", itsx_check = "",
    genbank_check = c("genbank_low_coverage(ITS 70%); genbank_multi_hsp(ITS 2 HSPs: 1-2, 3-4) top hits: ITS A x y 70%, rRNA none", "", "not_run(not checked)"),
    stringsAsFactors = FALSE)
  html <- as.character(RiboPilot:::export_check_summary_ui(smry))
  expect_match(html, "GenBank BLAST check: 1 flagged, 1 clean, 1 not checked")
  expect_match(html, "a \\(low coverage ITS 70%, ITS 2 HSPs\\)")
})
```

- [ ] **Step 2: Run** `Rscript -e 'devtools::test(filter="^export$|app-export")'`. Expected: fail.

- [ ] **Step 3: Implement.** `unit_checks()`: fetch `an.genbank_check` in the query, add `"GenBank BLAST" = nzchar(gb) & !startsWith(gb, "not_run(")` to `hit` right after `"ITSx"` (so it is a warning, not in `h[1:4]`). `export_summary_df()`: fetch like `itx`, NA -> `not_run(not checked)`, add `genbank_check = genbank` after `itsx_check`. `export_report_html()` `hdr` and the Rmd `cols`: append `"genbank_check"`. `export_check_summary_ui()`: after the ITSx block, same three-way branch with the text `GenBank BLAST check: %d flagged, %d clean, %d not checked.`, flagged list `sprintf("%s (%s)", seqid, genbank_check_display(gb))`, suffix `See the genbank_check column in the summary CSV.` README section, placed after the ITSx section, exact content (plain text, 80 columns):

```
GENBANK BLAST CHECK
-------------------
The genbank_check column in <label>_summary.csv predicts the two BLAST flags
GenBank's rRNA/ITS curation raises before it assigns accessions. NCBI runs a
default megablast of each record against its ITS_eukaryote_sequences database
and its eukaryote SSU/LSU rRNA database. Its reply names a record as
  "Query coverage is less than 75%."   when the top hit covers under 75% of it
  "HSP_problem"                        when the top hit aligns in more than one
                                       piece (a "potential misassembly")
RiboPilot runs the same two searches on each contig in the annotation step, and
again when Update on Annotation Review rechecks an edited contig, and stores:
  genbank_low_coverage(ITS <n>%)       NCBI's first flag, from the ITS search
  genbank_multi_hsp(ITS <n> HSPs: ...) NCBI's second flag, with the aligned
                                       pieces as positions on the record
  genbank_multi_hsp(rRNA ...)          the same from the SSU/LSU search
  genbank_no_hit(ITS) / (rRNA)         no hit at all in that database
  top hits: ...                        the top hit and its coverage per search
Empty means clean. not_run(...) means the check could not run; the README's
ITSx section says how to run a contig's checks again.

Expect these flags on most fish records, and on other groups the ITS database
covers thinly: it holds a few hundred fish records among 77,000 plant and
fungal ones, so the top hit is a distant relative, the ITS spacers do not
align to it, and the record splits into an 18S piece and a 5.8S piece. That
is not evidence of a bad assembly. NCBI's copy of the database is newer than
RiboPilot's, so a record near 75% can be called differently.

When GenBank sends the submission back for "low similarity and/or possible
misassembly", reply to gb-admin that the records below were reviewed with
RiboPilot's screen, chimera, per-feature BLAST and ITSx checks and are believed
correct; NCBI then runs further verification and proceeds. Records flagged in
this export:
  <one line per flagged seqid: seqid  display text>
```

The final list is built from `smry` at export time (`submission_readme(label)` gains a `flagged = character(0)` argument; pass `sprintf("  %s  %s", seqid, genbank_check_display(gb))` for flagged rows, or `  none` when empty).

- [ ] **Step 4: Run** `Rscript -e 'devtools::test(filter="export")'`. Expected: PASS (all export test files).

- [ ] **Step 5: Commit** `git add -A R inst/export_report.Rmd tests && git commit -m "feat: GenBank BLAST check in export CSV, window, README, report"`

---

### Task 5: App columns and help text

**Files:**
- Modify: `R/app_server.R:346` (scaffold table: add `genbank_check = cd("genbank_check", minWidth = 130, html = TRUE, cell = rt_gated(rt_itsx_badge()))` after `itsx_warnings`), `R/app_server.R:2192-2193` (Annotation Review table: same colDef with name "GenBank BLAST" and `rt_help_header("GenBank BLAST", SCAFFOLD_COL_HELP$genbank_check)`), `R/app_server.R:469` (add `"genbank_check"` to `Check details` after `"itsx_warnings"`), plus the SELECTs feeding those two tables (`v.genbank_check` next to `v.itsx_warnings` at :228 and the review table query near :2192)
- Modify: `R/constants.R:145` (add `SCAFFOLD_COL_HELP$genbank_check`)
- Test: `tests/testthat/test-app-table.R` (append)

**Interfaces:**
- Consumes: `genbank_check` from `review_units_v` (Task 1).

- [ ] **Step 1: Write the failing test** (mirror whatever test-app-table.R asserts for `itsx_warnings` columns: that the column is present in the table data and in the `Check details` group; add `genbank_check` to the same expectation). Also `expect_match(RiboPilot:::SCAFFOLD_COL_HELP$genbank_check, "GenBank")`.

- [ ] **Step 2: Run** `Rscript -e 'devtools::test(filter="app-table")'`. Expected: fail.

- [ ] **Step 3: Implement.** Help text:

```r
genbank_check = paste("The two BLAST flags GenBank's rRNA/ITS curation raises (default megablast against NCBI's",
  "ITS_eukaryote_sequences and SSU/LSU eukaryote rRNA databases): top-hit query coverage under 75%,",
  "and a top hit aligned in more than one piece. Most fish records trip both because the ITS database",
  "holds few fishes; it predicts GenBank's reply, not a bad assembly. Filled in by the annotation step,",
  "and by Update on Annotation Review after a trim or an edit. Empty means clean; 'not run' means",
  "the check could not run (hover for why); 'not checked' means none is stored yet.")
```

`rt_itsx_badge()` is reused as is (it already renders empty, `not_run(`, and flagged text).

- [ ] **Step 4: Run** `Rscript -e 'devtools::test(filter="app-table|app-headless|ui-reactable")'`. Expected: PASS.

- [ ] **Step 5: Commit** `git add -A R tests && git commit -m "app: GenBank BLAST column on scaffold and review tables"`

---

### Task 6: Ship ITS_eukaryote_sequences with the next panel release

**Files:**
- Modify: `R/blast_local.R:11` (add `RP_CHECK_DBS <- "ITS_eukaryote_sequences"` below `RP_DECOY_DBS`, with a two-line comment: shipped beside the panel for the GenBank check, never searched by the screen)
- Modify: `dev/refdb/fetch_decoys.sh:15` (`SEARCHABLE_DBS` gains `ITS_eukaryote_sequences`), `dev/refdb/package_release.sh:62` (loop list), `dev/refdb/build_refdb.sh:121` (`DECOY_DBS`)
- Modify: `docker/Dockerfile:131-133` (one more `blastdbcmd -db .../ITS_eukaryote_sequences -info` line in the verify chain; it must stay commented out or conditional until a release ships the db, because the current `RP_DB_VERSION` tarball lacks it and the image build would fail: add the line with a `# enable with the first panel release that ships ITS_eukaryote_sequences` comment above it, commented out)
- Modify: `R/db_dist.R:238` roxygen (decoys now "roughly 150 MB")
- Test: `tests/testthat/test-db-release.R:235-245`

**Interfaces:**
- Produces: `RP_CHECK_DBS`.

- [ ] **Step 1: Change the test** so the shell lists are expected to equal `c(RP_DECOY_DBS, RP_CHECK_DBS)` plus taxdb; rename the test "the shell decoy lists are the shipped databases plus exactly taxdb".

- [ ] **Step 2: Run** `Rscript -e 'devtools::test(filter="db-release")'`. Expected: fails on the lists.

- [ ] **Step 3: Edit the three scripts, the constant, the Dockerfile comment, and the roxygen line.**

- [ ] **Step 4: Run** `Rscript -e 'devtools::test(filter="db-release|blast-local|db-dist")'`. Expected: PASS.

- [ ] **Step 5: Commit** `git add -A R dev/refdb docker tests && git commit -m "refdb: ship ITS_eukaryote_sequences with the decoys"`

---

### Task 7: Vignette and the opt-in real-run test

**Files:**
- Modify: `vignettes/GenBank-Submission.Rmd` (after the ITSx table, :176)
- Modify: `tests/testthat/test-genbank-check.R` (append)

- [ ] **Step 1: Vignette block** headed `**The GenBank BLAST check.**`, same content as the README section condensed to two paragraphs plus a four-row table of values (`genbank_low_coverage(ITS <n>%)`, `genbank_multi_hsp(<db> <n> HSPs: ...)`, `genbank_no_hit(<db>)`, `not_run(...)`). Name the columns `genbank_check` (CSV) and **GenBank BLAST** (app).

- [ ] **Step 2: Real-run test**

```r
test_that("genbank_check_unit reproduces NCBI on SUB16547976 with real blastn", {
  d <- Sys.getenv("RIBOPILOT_GENBANK_CHECK_DBDIR", "")
  skip_if(!nzchar(d), "RIBOPILOT_GENBANK_CHECK_DBDIR not set")
  skip_if(!nzchar(Sys.which("blastn")), "blastn not on PATH")
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  fa <- readLines(gc_fix("SUB16547976.fasta"))
  ids <- sub("^>", "", fa[startsWith(fa, ">")]); seqs <- fa[!startsWith(fa, ">")]
  want <- read.delim(gc_fix("SUB16547976_ncbi_flags.tsv"), stringsAsFactors = FALSE)
  for (i in seq_along(ids)) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES (?,1,1,?,?)", params = list(ids[i], seqs[i], nchar(seqs[i])))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES (?,1,1)", params = list(ids[i]))
    RiboPilot:::genbank_check_unit(con, ids[i], 1L, 1L, db_dir = d)
  }
  got <- DBI::dbGetQuery(con, "SELECT ID, genbank_check FROM annotate")
  got <- got$genbank_check[match(want$seqid, got$ID)]
  message("ITS db build: ", system2("blastdbcmd", c("-db", shQuote(file.path(d, "ITS_eukaryote_sequences")), "-info"), stdout = TRUE)[3])
  expect_equal(grepl("genbank_low_coverage(ITS", got, fixed = TRUE), want$low_coverage == 1)
  expect_equal(grepl("genbank_multi_hsp(ITS", got, fixed = TRUE), want$hsp_problem == 1)
})
```

Run it once now with `RIBOPILOT_GENBANK_CHECK_DBDIR=/tmp/claude-1000/-home-dmacguig-Documents-GitHub-RiboPilot/7e6d4cd3-f648-4d90-ac10-f6edc95e8466/scratchpad/ncbidb` and `PATH` including `~/miniconda3/envs/blast/bin`. Expected: PASS and the message line shows `Oct 12, 2026`.

- [ ] **Step 3: Full suite** `Rscript -e 'devtools::test()'`. Expected: no new failures (compare against devel's known failures if any; see tests/testthat/_problems).

- [ ] **Step 4: Commit** `git add -A vignettes tests && git commit -m "docs: GenBank BLAST check in the submission vignette; real-run test"`

---

## Not in this plan

- No panel release, no image rebuild, no `RP_DB_VERSION` bump (user decision 2026-10-26). Until a release ships the ITS db, every contig stores `not_run(no ITS_eukaryote_sequences database)`.
- No backfill of locked contigs.
- No NEWS.md entry; the release flow writes it.
- No sample-level roll-up column of `genbank_check` on the Samples table (the ITSx one at R/app_server.R:187 has no counterpart here).
