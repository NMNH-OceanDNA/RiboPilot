# Contig Recheck Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Edits on Annotation Review (trim, gene-boundary nudge, 5'/3' choice, State) go to per-contig review copies; a second Nextflow workflow, launched by an Update button on Annotation Review, reruns every offline check on those copies; the original tables keep the pipeline's results for the Assemble step; export reads only the review layer and no longer runs ITSx.

**Architecture:** Eight `review_*` tables hold a copy of a contig once it is first touched; five SQL views return the copy when one exists and the original otherwise, so Annotation Review and export read one layer. The shared check functions gain a `review = FALSE` argument that redirects their reads and writes to the review tables, so the main pipeline and the recheck score with the same code. `main_recheck.nf` runs `query_recheck` -> `recheck_tools` (the same `rp_rrna_tools.sh` the main pipeline uses) -> `recheck` (`recheck_unit()` in R) per contig, guarded by a per-copy `rev` counter so an edit during a run always wins.

**Tech Stack:** R package (DBI/RSQLite, shiny, reactable, testthat 3e), SQLite views and savepoints, Nextflow DSL2 (25.10.x), bash, barrnap/pybarrnap, ITSx, BLAST+.

**Spec:** `dev/specs/2026-10-08-contig-recheck-design.md` (approved 2026-10-08). Section references below look like S5.2.

## Global Constraints

- Work on branch `feat/contig-recheck` (off `devel`). One commit per task, short message, no attribution lines. Stage only the task's files: the user has uncommitted edits in `vignettes/Troubleshooting.Rmd` and `vignettes/Your-Own-Project.Rmd` that no task may stage, revert, or overwrite. Never push.
- Schema goes from version 4 to 5 through one appended step 4 in `db_migrations()`; `new_db()` and step 4 share the DDL; `tests/testthat/fixtures/schema_v1.sql` and released steps 1-3 never change.
- The originals (`assemblies`, `annotations`, `annotate` check columns, `coverage`, `screen`, `chimera`, `chimera_windows`, `blast_hits`, `blast_searched`, `blast_screen`, `blast_ref_candidates`, `blast_ref_alignment`, `path_scores`) are written only by the main pipeline and the existing in-app exceptions; Annotation Review edits, the recheck, and export never write them (S4.1).
- Annotation Review and the Export step read contig data only through `review_units_v`, `review_annotations_v`, `review_screen_v`, `review_chimera_v`, `review_chimera_windows_v`; the Assemble step reads originals (S8).
- The recheck runs every offline check plus a local BLAST and never contacts NCBI (D1).
- Every recheck process is `cache false`, runs per contig, and is tagged `ID.path.scaffold`; the workflow launches from `<project>/.recheck` with `-work-dir <project>/work_recheck` and no `-resume` (S7.2, S7.4).
- The main pipeline and the recheck run the same barrnap/pybarrnap and ITSx commands through `inst/nextflow/bin/rp_rrna_tools.sh` (S7.1).
- The recheck's last write sets `check_switch = 2` only `WHERE rev = <passed rev> AND check_switch IN (1, 4)` (S5.4).
- Lock on Annotation Review is refused for contigs whose checks are not done; Exclude & Lock is exempt (D2).
- Not-run strings exactly: `not_run(ITSx did not finish)`, `not_run(ITSx check failed: <message>)`, `not_run(not checked before upgrade)`, and `not_run(not checked)` for NA in the export summary (S6.3, S9, S11).
- ASCII only in code, tests, and docs. User-facing text says "contig" and "step" ("Assemble step", "Annotation Review", "Export step"), never "<Step> tab".
- NAMESPACE is hand-maintained: never run `devtools::document()`. Changed exported functions get their Rd regenerated with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`; new internals stay `@noRd` and bin scripts call them as `RiboPilot:::fn`.
- Until Task 22 lands, never run `backwards_compatibility()` or the app on a real project; use a fresh test project or a copied `.sqlite`.
- Real-Nextflow and other skip_on_cran tests need `NOT_CRAN=true` in the environment, or they skip.
- The shared names, signatures, and SQL in each task's Interfaces block are the contract between tasks. Do not rename them.

## Review Focus

1. **Upgrading real version-4 projects (Task 22).** One bad unit must not abort step 4, because a failed upgrade locks the user out of the project. Expected: the step finishes and lists degraded units. Task 22 tests a version-4 fixture holding a current sequence not found in its backup, a NULL backup sequence, a snapshot JSON that no longer fits, a backup row with no annotate row, and a locked unit with NULL `itsx_warnings`.
2. **Trims at and through gene ends on + and - strand contigs, then export (Tasks 4, 19).** Expected: cut ends get `partial_*_user` in the feature's own orientation, a trim exactly at a gene end marks nothing partial, a trimmed end is never an automatic assembly end, and the `<`/`>` marks in the `.tbl` and the defline follow. Task 4 tests both strands, an exact-end trim, and a second trim; Task 19 tests a minus-strand trimmed copy end to end.
3. **Tool failures inside `rp_rrna_tools.sh` (Tasks 9, 10, 11).** Expected: in `annotate_tools` a barrnap exit 1 fails the task as today and exit 137 still reaches the out-of-memory retry; in `recheck_tools` a failure is recorded and only out-of-memory codes fail the task; a stock ITSx crash gives `not_run(ITSx did not finish)`. Task 9 tests the script's exit status and status file with failing fake tools; Task 10 tests `recheck_unit` with failed tools; Task 11 tests the `recheck_tools` exit rule.
4. **No local BLAST panel at recheck time (Tasks 10, 11).** Expected: an empty `dbarg` (no panel, `rp_dbarg.R` failure, or `RIBOPILOT_NO_LOCAL_BLAST`) falls back to mapped hits and the contig still finishes, never calling NCBI. Task 10 tests `recheck_hits` with no `dbarg`; Task 11 runs the workflow with `RIBOPILOT_NO_LOCAL_BLAST=1`.
5. **Upgraded databases with ALTER-appended columns (Tasks 1, 2).** Every user project is version 1 plus the steps, so column order differs from `new_db()`. Expected: the views and `review_copy()` work there. Task 1 queries all five views on a version-1 database upgraded through every step and holding data; Task 2 runs `review_copy()`, `review_invalidate()`, and `review_drop()` on the same.

## Task overview

| Task | What it delivers |
|---|---|
| 1 | Review schema, views, ITSx columns, step 4, fixture helpers |
| 2 | Review copy lifecycle: savepoints, copy, invalidate, drop |
| 3 | Coverage window and mapped BLAST hits for copies |
| 4 | Review edit functions: trim, bounds, partials, reset, State |
| 5 | ITSx states: finished-run detection, not_run rules, stored origin and scores |
| 6 | Review mode for curation, reference alignment, and the live alignment view |
| 7 | Review mode for BLAST products and the screen |
| 8 | Review mode for chimera and the stock ITSx check |
| 9 | Shared rRNA tool script used by the main pipeline |
| 10 | recheck_unit: the per-contig recheck in R |
| 11 | Recheck Nextflow workflow and its command |
| 12 | Main pipeline follows copies; taxonomy refresh |
| 13 | Update button and run window on Annotation Review |
| 14 | Annotation Review table: State column, gated checks, Checks column |
| 15 | Trim and Restore on copies |
| 16 | Edit mode and State on copies |
| 17 | Details window readers and the Assemble chimera panel |
| 18 | Export without the ITSx run |
| 19 | Export reads the review layer |
| 20 | Export step table |
| 21 | Lock gate |
| 22 | Upgrade step 4 conversion; assembly_backup removed |
| 23 | Docs and stale comments |
| 24 | Full verification |

---

### Task 1: Review schema, views, ITSx columns, step 4, fixture helpers

**Files:**
- Create: `tests/testthat/test-review-schema.R`
- Modify: `R/db_schema.R` (after `CHIMERA_WINDOWS_DDL`, :143; after `meta_apply_schema()`, :262-264)
- Modify: `R/init_db.R` (`new_db()`: blast_screen :392-413, blast_ref_alignment :447-461, annotations :484-510, annotate :538, after :604)
- Modify: `R/db_migrate.R` (header comment :3-4; `db_migrations()` :16-20)
- Modify: `tests/testthat/helper-app-fixture.R` (end of `fixture_fill_schema()`, :59-61; three new helpers)
- Modify: `tests/testthat/test-backwards-compat.R` (after `bc_schema()`, :23-24; parity test :127-128)
- Modify: `tests/testthat/test-batch-trim.R` (`mk_batch_project()`, :11-16)
- Modify: `tests/testthat/test-app-export.R` (`add_trim_cols()`, :562-570; one test appended at the end)
- Test: `tests/testthat/test-review-schema.R`, `tests/testthat/test-backwards-compat.R`, `tests/testthat/test-app-export.R`

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `R/db_schema.R` constants (all `CREATE TABLE IF NOT EXISTS`): `ANNOTATIONS_DDL`, `BLAST_SCREEN_DDL`, `BLAST_REF_ALIGNMENT_DDL` (used by `new_db()`); `REVIEW_UNITS_DDL`; mirrors `REVIEW_ANNOTATIONS_DDL` (annotations columns plus `partial_start_user INTEGER, partial_stop_user INTEGER` before the PRIMARY KEY), `REVIEW_BLAST_HITS_DDL`, `REVIEW_BLAST_SCREEN_DDL`, `REVIEW_SCREEN_DDL`, `REVIEW_CHIMERA_DDL`, `REVIEW_CHIMERA_WINDOWS_DDL`, `REVIEW_REF_ALIGNMENT_DDL`.
  - `REVIEW_MIRRORS`: named chr, original -> mirror: `annotations = "review_annotations", blast_hits = "review_blast_hits", blast_screen = "review_blast_screen", screen = "review_screen", chimera = "review_chimera", chimera_windows = "review_chimera_windows", blast_ref_alignment = "review_ref_alignment"`.
  - `REVIEW_TABLES`: `c("review_units", unname(REVIEW_MIRRORS))` (8 names, review_units first).
  - `REVIEW_VIEWS_DDL`: named chr of 5 `CREATE VIEW IF NOT EXISTS` statements, names `review_units_v`, `review_annotations_v`, `review_screen_v`, `review_chimera_v`, `review_chimera_windows_v`.
  - `review_units_v` columns, in order: `ID, path, scaffold, has_copy, topology, orig_length, sequence, length, trim_start, trim_end, edited, check_switch, check_note, rev, checked_at, n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores, blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp`. Copy branch: review_units JOIN assemblies (topology, orig_length). Fallback: assemblies JOIN annotate, `has_copy 0`, `trim_start 1`, `trim_end = length = orig_length = assemblies.length`, `edited 0`, `check_switch = CASE WHEN an.annotate_switch = 2 AND an.itsx_warnings IS NULL THEN 1 ELSE an.annotate_switch END`, check_note/rev/checked_at NULL, check columns and time_stamp from annotate, blast_* from assemblies.
  - `review_annotations_v` columns: the 18 annotations columns, then `partial_start_user, partial_stop_user` (NULL on fallback rows). `review_screen_v`, `review_chimera_v`, `review_chimera_windows_v`: exactly their original's columns, in order. Mirror rows count only when the unit has a `review_units` row; a unit with a copy never shows original rows.
  - `review_apply_schema(con)`: the 8 tables, then the 5 views; idempotent; no transaction of its own.
  - Internal helpers in `R/db_schema.R`: `review_mirror_ddl(ddl, orig)`, `ddl_cols(ddl)`, `review_view_ddl(view, orig, cols, fallback = paste0("x.", cols))`, `REVIEW_COPY_SQL` (the `EXISTS` subquery on `review_units` that matches the unit of the row aliased `x`).
  - `annotate.itsx_origin TEXT`, `annotate.itsx_scores TEXT`, the last two annotate columns (after `export_time_stamp`) in new_db and after step 4.
  - `db_migrations()[[4]]`: note `"review copies for edited contigs; ITSx origin and scores"`, fn = the two ALTERs, then `review_apply_schema(con)`. `current_schema_version()` is 5. Task 22 appends its conversion to this fn.
  - Test helpers (`tests/testthat/helper-app-fixture.R`): `fixture_review_schema(con)`, `fixture_mark_done(con)`, `fixture_ref_schema()` (session cache). `fixture_fill_schema()` calls `fixture_review_schema(con)` last, so `mk_app_project()` and `mk_export_project()` carry every canonical column of assemblies, annotate, annotations, coverage, screen, chimera, chimera_windows, blast_hits, blast_screen, blast_ref_alignment, plus the review tables and views.

This task adds the data layer the rest of the plan builds on (S4.1-4.3, S11.1-2): eight `review_*` tables, five views that return a contig's review copy when it has one and its original rows otherwise, and the two `annotate` columns `itsx_origin` and `itsx_scores`. `new_db()` and the new upgrade step 4 run the same DDL constants, and the parity test now compares view SQL, so a version 1 database run through every step matches a fresh project. Every view branch names its columns because `UNION ALL` matches columns by position and column order differs between `new_db()`, upgraded databases, and hand-built fixtures. The hand-built app fixtures get a helper that widens them to the real columns and adds the review schema, so later tasks can point app code at the views without patching each fixture. The version 1 upgrade test below covers Review Focus 5 for this task. Nothing reads the views yet.

- [ ] **Step 1: Write the failing schema tests**

Create `tests/testthat/test-review-schema.R`:

```r
# Review copies and the views over them (SPEC S4.2, S4.3, S11).

RS_SEQ <- strrep("ACGTACGTAC", 10)

# Four 100 bp contigs of S1 in the original tables. Scaffold 1 is done without
# a stored ITSx value, 2 is done, 3 is running in the main pipeline with no
# stored ITSx value. Scaffold 5 is a dropped contig with no annotate row; the
# review layer does not list it.
rs_originals <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence,
        blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
      VALUES ('S1', 1, ?, 'linear', 100, ?, 'A1', 'Danio rerio', 99.5, 88, 1e-50, 'core_nt', 999)",
      params = list(s, RS_SEQ))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, time_stamp)
      VALUES ('S1', 1, ?, ?, 0, 2, 'ITS1', 0, 3, 'length', ?, 111)",
      params = list(s, c(2L, 2L, 4L, 2L)[s], c(NA, "", NA, "")[s]))
    for (f in list(c("18S", 1, 50), c("28S", 60, 100))) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
        VALUES ('S1', 1, ?, ?, ?, ?, '+')", params = list(s, f[1], as.integer(f[2]), as.integer(f[3])))
    }
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, ?, 'on_target')",
      params = list(s))
    DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1', 1, ?, 'clean')",
      params = list(s))
    for (w in c(1L, 51L)) {
      DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, ?, ?, ?)",
        params = list(s, w, w + 49L))
    }
  }
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 5, 'linear', 100, ?, 1)", params = list(RS_SEQ))
}

# The stored ITSx origin of every contig; a copy of scaffold 4 (trimmed to
# 11..90, one feature, no chimera row); and rows a discarded recheck left for
# scaffold 1, which has no copy.
rs_copy <- function(con) {
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M', itsx_scores = 'M:50'")
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, trim_start,
      trim_end, edited, check_switch, check_note, rev, checked_at, n_features, missing, chimeric,
      n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores, blast_accession, blast_species,
      blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
    VALUES ('S1', 1, 4, ?, 80, 11, 90, 1, 3, 'boom', CAST(? AS INTEGER), 333, 1, '28S', 1, 0, '',
      'itsx_x', 'F', 'F:10', 'B2', 'Copy sp', 90, 70, 1e-10, 'local', 222)",
    params = list(substr(RS_SEQ, 11, 90), "1760000000123"))
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, start, stop, strand,
      partial_start_user) VALUES ('S1', 1, 4, '18S', 1, 40, '+', 1)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 4, 'uncertain')")
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, 4, 1, 50)")
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, start, stop)
    VALUES ('S1', 1, 1, 'ITS1', 51, 59)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 1, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label) VALUES ('S1', 1, 1, 'chimeric')")
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, 1, 7, 8)")
}

rs_check_views <- function(con) {
  u <- DBI::dbGetQuery(con, "SELECT scaffold, has_copy, topology, orig_length, sequence, length,
      trim_start, trim_end, edited, check_switch, check_note, CAST(rev AS TEXT) AS rev, checked_at,
      n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores,
      blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp
    FROM review_units_v ORDER BY scaffold")
  expect_equal(u$scaffold, 1:4)
  expect_equal(u$has_copy, c(0L, 0L, 0L, 1L))
  expect_equal(u$topology, rep("linear", 4))
  expect_equal(u$orig_length, rep(100L, 4))
  expect_equal(u$sequence, c(rep(RS_SEQ, 3), substr(RS_SEQ, 11, 90)))
  expect_equal(u$length, c(100L, 100L, 100L, 80L))
  expect_equal(u$trim_start, c(1L, 1L, 1L, 11L))
  expect_equal(u$trim_end, c(100L, 100L, 100L, 90L))
  expect_equal(u$edited, c(0L, 0L, 0L, 1L))
  # done without a stored ITSx value reads as not run
  expect_equal(u$check_switch, c(1L, 2L, 4L, 3L))
  expect_equal(u$check_note, c(NA, NA, NA, "boom"))
  expect_equal(u$rev, c(NA, NA, NA, "1760000000123"))
  expect_equal(u$checked_at, c(NA, NA, NA, 333L))
  expect_equal(u$n_features, c(2L, 2L, 2L, 1L))
  expect_equal(u$missing, c("ITS1", "ITS1", "ITS1", "28S"))
  expect_equal(u$chimeric, c(0L, 0L, 0L, 1L))
  expect_equal(u$n_flags, c(3L, 3L, 3L, 0L))
  expect_equal(u$curate_flags, c("length", "length", "length", ""))
  expect_equal(u$itsx_warnings, c(NA, "", NA, "itsx_x"))
  expect_equal(u$itsx_origin, c("M", "M", "M", "F"))
  expect_equal(u$itsx_scores, c("M:50", "M:50", "M:50", "F:10"))
  expect_equal(u$blast_accession, c("A1", "A1", "A1", "B2"))
  expect_equal(u$blast_species, c("Danio rerio", "Danio rerio", "Danio rerio", "Copy sp"))
  expect_equal(u$blast_pident, c(99.5, 99.5, 99.5, 90))
  expect_equal(u$blast_qcovs, c(88, 88, 88, 70))
  expect_equal(u$blast_evalue, c(1e-50, 1e-50, 1e-50, 1e-10))
  expect_equal(u$blast_source, c("core_nt", "core_nt", "core_nt", "local"))
  expect_equal(u$time_stamp, c(111L, 111L, 111L, 222L))

  a <- DBI::dbGetQuery(con, "SELECT scaffold, element, start, stop, partial_start_user
    FROM review_annotations_v ORDER BY scaffold, start")
  expect_equal(a$scaffold, c(1L, 1L, 2L, 2L, 3L, 3L, 4L))
  expect_equal(a$element, c(rep(c("18S", "28S"), 3), "18S"))
  expect_equal(a$stop, c(rep(c(50L, 100L), 3), 40L))
  expect_equal(a$partial_start_user, c(rep(NA, 6), 1L))

  s <- DBI::dbGetQuery(con, "SELECT scaffold, screen_verdict FROM review_screen_v ORDER BY scaffold")
  expect_equal(s$scaffold, 1:4)
  expect_equal(s$screen_verdict, c("on_target", "on_target", "on_target", "uncertain"))
  ch <- DBI::dbGetQuery(con, "SELECT scaffold, label FROM review_chimera_v ORDER BY scaffold")
  expect_equal(ch$scaffold, 1:3)
  expect_equal(ch$label, rep("clean", 3))
  w <- DBI::dbGetQuery(con, "SELECT scaffold, start FROM review_chimera_windows_v ORDER BY scaffold, start")
  expect_equal(w$scaffold, c(1L, 1L, 2L, 2L, 3L, 3L, 4L))
  expect_equal(w$start, c(rep(c(1L, 51L), 3), 1L))
}

test_that("new_db() makes each review table like its original, and annotate gains the ITSx columns", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  info <- function(t) DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", t))[,
    c("name", "type", "notnull", "dflt_value", "pk")]
  m <- RiboPilot:::REVIEW_MIRRORS
  expect_equal(names(m), c("annotations", "blast_hits", "blast_screen", "screen", "chimera",
    "chimera_windows", "blast_ref_alignment"))
  expect_equal(RiboPilot:::REVIEW_TABLES, c("review_units", unname(m)))
  for (orig in names(m)[-1]) expect_equal(info(m[[orig]]), info(orig), info = orig)
  ra <- info("review_annotations")
  expect_equal(nrow(ra), 20L)
  expect_equal(ra[1:18, ], info("annotations"))
  expect_equal(ra$name[19:20], c("partial_start_user", "partial_stop_user"))
  expect_equal(ra$type[19:20], c("INTEGER", "INTEGER"))

  ru <- info("review_units")
  expect_equal(ru$name, c("ID", "path", "scaffold", "sequence", "length", "trim_start", "trim_end",
    "edited", "check_switch", "check_note", "rev", "checked_at", "n_features", "missing", "chimeric",
    "n_flags", "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession",
    "blast_species", "blast_pident", "blast_qcovs", "blast_evalue", "blast_source", "time_stamp"))
  expect_equal(ru$pk, c(1:3, rep(0L, 24)))
  expect_equal(ru$notnull[ru$name %in% c("ID", "path", "scaffold", "edited", "rev")], rep(1L, 5))
  expect_equal(ru$dflt_value[ru$name == "edited"], "0")

  expect_equal(utils::tail(DBI::dbListFields(con, "annotate"), 3),
    c("export_time_stamp", "itsx_origin", "itsx_scores"))

  views <- DBI::dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type = 'view'")$name
  expect_setequal(views, names(RiboPilot:::REVIEW_VIEWS_DDL))
  expect_equal(DBI::dbListFields(con, "review_units_v"), c("ID", "path", "scaffold", "has_copy",
    "topology", "orig_length", "sequence", "length", "trim_start", "trim_end", "edited",
    "check_switch", "check_note", "rev", "checked_at", "n_features", "missing", "chimeric", "n_flags",
    "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession", "blast_species",
    "blast_pident", "blast_qcovs", "blast_evalue", "blast_source", "time_stamp"))
  expect_equal(DBI::dbListFields(con, "review_annotations_v"), ra$name)
  for (t in c("screen", "chimera", "chimera_windows")) {
    expect_equal(DBI::dbListFields(con, paste0("review_", t, "_v")), DBI::dbListFields(con, t), info = t)
  }
})

test_that("review_apply_schema() is idempotent", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  before <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name")
  RiboPilot:::review_apply_schema(con)
  RiboPilot:::review_apply_schema(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name"), before)
})

test_that("the views read a contig's copy when it has one and its originals otherwise", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rs_originals(con)
  rs_copy(con)
  rs_check_views(con)

  # a copy whose features were all trimmed away shows none, not the original ones
  n_feat <- function() DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM review_annotations_v WHERE scaffold = 4")$n
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE scaffold = 4")
  expect_equal(n_feat(), 0)

  # without its review_units row the contig reads its originals again
  DBI::dbExecute(con, "DELETE FROM review_units WHERE scaffold = 4")
  expect_equal(n_feat(), 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_verdict FROM review_screen_v WHERE scaffold = 4")$screen_verdict,
    "on_target")
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v WHERE scaffold = 4"),
    data.frame(has_copy = 0L, length = 100L, check_switch = 2L))
})

test_that("a version 1 database holding data upgrades through step 4 and its views work", {
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  rs_originals(con)
  DBI::dbDisconnect(con)

  suppressMessages(backwards_compatibility(dir))
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  h <- DBI::dbGetQuery(con, "SELECT schema_version, note FROM schema_history ORDER BY rowid")
  expect_equal(h$schema_version[h$note %in% "review copies for edited contigs; ITSx origin and scores"], 5L)
  expect_equal(utils::tail(DBI::dbListFields(con, "annotate"), 3),
    c("export_time_stamp", "itsx_origin", "itsx_scores"))
  rs_copy(con)
  rs_check_views(con)
})
```

In `tests/testthat/test-backwards-compat.R`, the parity test now also compares the views (S11):

Find:

```r
  stats::setNames(lapply(tbls, function(t) sort(DBI::dbListFields(con, t))), tbls)
}
```

Replace with:

```r
  stats::setNames(lapply(tbls, function(t) sort(DBI::dbListFields(con, t))), tbls)
}

bc_views <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  v <- DBI::dbGetQuery(con, "SELECT name, sql FROM sqlite_master WHERE type = 'view' ORDER BY name")
  stats::setNames(gsub("\\s+", " ", v$sql), v$name)
}
```

Find:

```r
  expect_equal(bc_schema(old), bc_schema(file.path(dir, ".sqlite")))
})
```

Replace with:

```r
  expect_equal(bc_schema(old), bc_schema(file.path(dir, ".sqlite")))
  expect_equal(bc_views(old), bc_views(file.path(dir, ".sqlite")))
  expect_setequal(names(bc_views(old)), names(RiboPilot:::REVIEW_VIEWS_DDL))
})
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "backwards-compat")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-schema.R `[ FAIL 6 | WARN 0 | SKIP 0 | PASS 0 ]` with `object 'REVIEW_MIRRORS' not found`, `object 'review_apply_schema' not found`, `no such column: itsx_origin`, and the step 4 note missing from schema_history. test-backwards-compat.R `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 37 ]`: "a v1 database run through every step matches a fresh new_db()" errors with `object 'REVIEW_VIEWS_DDL' not found`.

- [ ] **Step 3: Add the DDL constants, the views, and review_apply_schema() to R/db_schema.R**

Find (end of `CHIMERA_WINDOWS_DDL`):

```r
  PRIMARY KEY (ID, path, scaffold, start)
);"
```

Replace with:

```r
  PRIMARY KEY (ID, path, scaffold, start)
);"

#' DDL for annotations: the nrDNA feature model per unit.
#'
#' status/flags are filled by curation (curate_one); tool/edited by annotation.
#' partial_start/partial_stop: 5'/3' end partial in the feature's own orientation.
#' @noRd
ANNOTATIONS_DDL <- "CREATE TABLE IF NOT EXISTS annotations (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  element TEXT NOT NULL,
  type TEXT,
  start INTEGER,
  stop INTEGER,
  length INTEGER,
  strand TEXT,
  model_coverage REAL,
  partial_start INTEGER,
  partial_stop INTEGER,
  tool TEXT,
  status TEXT,
  flags TEXT,
  edited INTEGER,
  notes TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, element)
);"

#' DDL for blast_screen: the unfiltered top remote hit the contaminant check reads.
#' @noRd
BLAST_SCREEN_DDL <- "CREATE TABLE IF NOT EXISTS blast_screen (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  accession TEXT,
  title TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  is_nuclear_rrna INTEGER,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for blast_ref_alignment: per-feature sample-vs-reference alignment.
#' @noRd
BLAST_REF_ALIGNMENT_DDL <- "CREATE TABLE IF NOT EXISTS blast_ref_alignment (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  accession TEXT NOT NULL,
  element TEXT NOT NULL,
  pid REAL,
  aln_cov REAL,
  aln_len INTEGER,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, accession, element)
);"
```

Find (end of `meta_apply_schema()`, the last lines of the file):

```r
  .nmnh_ensure_table(con)
  invisible()
}
```

Replace with:

```r
  .nmnh_ensure_table(con)
  invisible()
}

#' DDL for review_units: a contig's review copy; a row here means it has one.
#'
#' trim_start/trim_end: the kept window in original coordinates. check_switch:
#' 2 done, 1 not run, 4 running, 3 failed, 0 on hold. rev starts at the time in
#' milliseconds and goes up by one on every invalidation.
#' @noRd
REVIEW_UNITS_DDL <- "CREATE TABLE IF NOT EXISTS review_units (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL,
  scaffold INTEGER NOT NULL,
  sequence TEXT,
  length INTEGER,
  trim_start INTEGER,
  trim_end INTEGER,
  edited INTEGER NOT NULL DEFAULT 0,
  check_switch INTEGER,
  check_note TEXT,
  rev INTEGER NOT NULL,
  checked_at INTEGER,
  n_features INTEGER,
  missing TEXT,
  chimeric INTEGER,
  n_flags INTEGER,
  curate_flags TEXT,
  itsx_warnings TEXT,
  itsx_origin TEXT,
  itsx_scores TEXT,
  blast_accession TEXT,
  blast_species TEXT,
  blast_pident REAL,
  blast_qcovs REAL,
  blast_evalue REAL,
  blast_source TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' Review mirror of each original per-contig table
#' @noRd
REVIEW_MIRRORS <- c(annotations = "review_annotations", blast_hits = "review_blast_hits",
  blast_screen = "review_blast_screen", screen = "review_screen", chimera = "review_chimera",
  chimera_windows = "review_chimera_windows", blast_ref_alignment = "review_ref_alignment")

#' The eight review tables, review_units first
#' @noRd
REVIEW_TABLES <- c("review_units", unname(REVIEW_MIRRORS))

#' An original's DDL with the table renamed to its mirror
#' @noRd
review_mirror_ddl <- function(ddl, orig) {
  sub(paste0("CREATE TABLE IF NOT EXISTS ", orig, " ("),
    paste0("CREATE TABLE IF NOT EXISTS ", REVIEW_MIRRORS[[orig]], " ("), ddl, fixed = TRUE)
}

REVIEW_ANNOTATIONS_DDL <- sub("  PRIMARY KEY",
  "  partial_start_user INTEGER,\n  partial_stop_user INTEGER,\n  PRIMARY KEY",
  review_mirror_ddl(ANNOTATIONS_DDL, "annotations"), fixed = TRUE)
REVIEW_BLAST_HITS_DDL <- review_mirror_ddl(BLAST_HITS_DDL, "blast_hits")
REVIEW_BLAST_SCREEN_DDL <- review_mirror_ddl(BLAST_SCREEN_DDL, "blast_screen")
REVIEW_SCREEN_DDL <- review_mirror_ddl(SCREEN_DDL, "screen")
REVIEW_CHIMERA_DDL <- review_mirror_ddl(CHIMERA_DDL, "chimera")
REVIEW_CHIMERA_WINDOWS_DDL <- review_mirror_ddl(CHIMERA_WINDOWS_DDL, "chimera_windows")
REVIEW_REF_ALIGNMENT_DDL <- review_mirror_ddl(BLAST_REF_ALIGNMENT_DDL, "blast_ref_alignment")

#' Column names of a CREATE TABLE constant, in order
#' @noRd
ddl_cols <- function(ddl) {
  x <- sub("[ ,].*", "", trimws(strsplit(ddl, "\n", fixed = TRUE)[[1]][-1]))
  x[grepl("^[A-Za-z_]", x) & x != "PRIMARY"]
}

REVIEW_COPY_SQL <- "EXISTS (SELECT 1 FROM review_units r
  WHERE r.ID = x.ID AND r.path = x.path AND r.scaffold = x.scaffold)"

#' CREATE VIEW: mirror rows of contigs with a copy, original rows of the rest
#' @noRd
review_view_ddl <- function(view, orig, cols, fallback = paste0("x.", cols)) {
  paste0("CREATE VIEW IF NOT EXISTS ", view, " AS\nSELECT ", paste0("x.", cols, collapse = ", "),
    "\nFROM ", REVIEW_MIRRORS[[orig]], " x\nWHERE ", REVIEW_COPY_SQL,
    "\nUNION ALL\nSELECT ", paste(fallback, collapse = ", "),
    "\nFROM ", orig, " x\nWHERE NOT ", REVIEW_COPY_SQL)
}

#' The review layer Annotation Review and export read
#'
#' A contig nobody touched reads its originals. In review_units_v such a contig
#' reads as not run (1) when it is done (2) but has no stored ITSx value.
#' @noRd
REVIEW_VIEWS_DDL <- c(
  review_units_v = paste0("CREATE VIEW IF NOT EXISTS review_units_v AS
SELECT x.ID, x.path, x.scaffold, 1 AS has_copy, a.topology, a.length AS orig_length,
  x.sequence, x.length, x.trim_start, x.trim_end, x.edited, x.check_switch, x.check_note,
  x.rev, x.checked_at, x.n_features, x.missing, x.chimeric, x.n_flags, x.curate_flags,
  x.itsx_warnings, x.itsx_origin, x.itsx_scores, x.blast_accession, x.blast_species,
  x.blast_pident, x.blast_qcovs, x.blast_evalue, x.blast_source, x.time_stamp
FROM review_units x
JOIN assemblies a ON a.ID = x.ID AND a.path = x.path AND a.scaffold = x.scaffold
UNION ALL
SELECT x.ID, x.path, x.scaffold, 0, x.topology, x.length,
  x.sequence, x.length, 1, x.length, 0,
  CASE WHEN an.annotate_switch = 2 AND an.itsx_warnings IS NULL THEN 1 ELSE an.annotate_switch END,
  NULL, NULL, NULL, an.n_features, an.missing, an.chimeric, an.n_flags, an.curate_flags,
  an.itsx_warnings, an.itsx_origin, an.itsx_scores, x.blast_accession, x.blast_species,
  x.blast_pident, x.blast_qcovs, x.blast_evalue, x.blast_source, an.time_stamp
FROM assemblies x
JOIN annotate an ON an.ID = x.ID AND an.path = x.path AND an.scaffold = x.scaffold
WHERE NOT ", REVIEW_COPY_SQL),
  review_annotations_v = review_view_ddl("review_annotations_v", "annotations",
    c(ddl_cols(ANNOTATIONS_DDL), "partial_start_user", "partial_stop_user"),
    c(paste0("x.", ddl_cols(ANNOTATIONS_DDL)), "NULL AS partial_start_user", "NULL AS partial_stop_user")),
  review_screen_v = review_view_ddl("review_screen_v", "screen", ddl_cols(SCREEN_DDL)),
  review_chimera_v = review_view_ddl("review_chimera_v", "chimera", ddl_cols(CHIMERA_DDL)),
  review_chimera_windows_v = review_view_ddl("review_chimera_windows_v", "chimera_windows",
    ddl_cols(CHIMERA_WINDOWS_DDL)))

#' Create the review tables, then the views over them (idempotent)
#' @noRd
review_apply_schema <- function(con) {
  for (ddl in c(REVIEW_UNITS_DDL, REVIEW_ANNOTATIONS_DDL, REVIEW_BLAST_HITS_DDL,
                REVIEW_BLAST_SCREEN_DDL, REVIEW_SCREEN_DDL, REVIEW_CHIMERA_DDL,
                REVIEW_CHIMERA_WINDOWS_DDL, REVIEW_REF_ALIGNMENT_DDL, REVIEW_VIEWS_DDL)) {
    DBI::dbExecute(con, ddl)
  }
  invisible()
}
```

Then check the file parses: `Rscript -e 'invisible(parse("R/db_schema.R"))'`.

- [ ] **Step 4: Build them in new_db() (R/init_db.R)**

Find:

```r
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_screen (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      accession TEXT,
      title TEXT,
      species TEXT,
      staxid INTEGER,
      lineage TEXT,
      lineage_taxids TEXT,
      rank_distance INTEGER,
      pident REAL,
      qcovs REAL,
      evalue REAL,
      bitscore REAL,
      is_nuclear_rrna INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    );"
  )
```

Replace with:

```r
  DBI::dbExecute(con, BLAST_SCREEN_DDL)
```

Find:

```r
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_ref_alignment (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      accession TEXT NOT NULL,
      element TEXT NOT NULL,
      pid REAL,
      aln_cov REAL,
      aln_len INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, accession, element)
    );"
  )
```

Replace with:

```r
  DBI::dbExecute(con, BLAST_REF_ALIGNMENT_DDL)
```

Find:

```r
  # annotations (nrDNA feature model per assembly unit) ----
  # status/flags are filled by curation (curate_one); tool/edited by annotation.
  # partial_start/partial_stop: 5'/3' end partial in the feature's own orientation.
  DBI::dbExecute(
    con,
    "CREATE TABLE annotations (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      element TEXT NOT NULL,
      type TEXT,
      start INTEGER,
      stop INTEGER,
      length INTEGER,
      strand TEXT,
      model_coverage REAL,
      partial_start INTEGER,
      partial_stop INTEGER,
      tool TEXT,
      status TEXT,
      flags TEXT,
      edited INTEGER,
      notes TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, element)
    );"
  )
```

Replace with:

```r
  # annotations (nrDNA feature model per assembly unit; see ANNOTATIONS_DDL) ----
  DBI::dbExecute(con, ANNOTATIONS_DDL)
```

Find (inside `CREATE TABLE annotate`):

```r
      export_time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
```

Replace with:

```r
      export_time_stamp INTEGER,
      itsx_origin TEXT,
      itsx_scores TEXT,
      PRIMARY KEY (ID, path, scaffold)
```

Find:

```r
  DBI::dbExecute(con, CHIMERA_WINDOWS_DDL)
```

Replace with:

```r
  DBI::dbExecute(con, CHIMERA_WINDOWS_DDL)

  # review copies of edited contigs and the views over them (see REVIEW_VIEWS_DDL) ----
  review_apply_schema(con)
```

Every table a view reads (assemblies, annotate, annotations, screen, chimera, chimera_windows) exists by then. Check: `Rscript -e 'invisible(parse("R/init_db.R"))'`.

- [ ] **Step 5: Append upgrade step 4 (R/db_migrate.R)**

Find:

```r
# Any table or column change edits new_db() AND appends a step here. Step k
# takes a version k database to version k + 1. Released steps are never edited.
```

Replace with:

```r
# Any table or column change edits new_db() AND appends a step here. Step k
# takes a version k database to version k + 1. Released steps are never edited.
# Views (REVIEW_VIEWS_DDL) are CREATE VIEW IF NOT EXISTS: a step that changes one
# drops it first. While a view is broken every ALTER TABLE RENAME fails, so a
# step that rebuilds a table a view reads drops the views and reruns
# review_apply_schema() afterwards.
```

Find:

```r
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
      })
  )
```

Replace with:

```r
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
      }),
    list(note = "review copies for edited contigs; ITSx origin and scores",
      fn = function(con) {
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_origin TEXT")
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_scores TEXT")
        review_apply_schema(con)
      })
  )
```

The ALTERs run in the same order as the columns appear in new_db's annotate. `current_schema_version()` becomes 5 on its own (1 + 4 steps). Check: `Rscript -e 'invisible(parse("R/db_migrate.R"))'`.

- [ ] **Step 6: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "backwards-compat")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-schema.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 106 ]`, test-backwards-compat.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 38 ]`.

Then the files that build a real `new_db()` (they must stay green; the new columns and views change nothing they read):
Run: `Rscript -e 'devtools::load_all("."); for (f in c("chimera", "curate", "load-assemblies", "export", "utils", "init-assemblies", "reference-schema", "nextflow", "nf-pending", "select-paths", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0`.

- [ ] **Step 7: Write the failing fixture tests**

Append to the end of `tests/testthat/test-review-schema.R`:

```r

test_that("the app fixture has every column the review layer reads, and its views answer", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ref <- chimera_test_con()
  for (t in c("assemblies", "annotate", "annotations", "coverage", "screen", "chimera",
              "chimera_windows", "blast_hits", "blast_screen", "blast_ref_alignment")) {
    expect_true(all(DBI::dbListFields(ref, t) %in% DBI::dbListFields(con, t)), info = t)
  }
  DBI::dbDisconnect(ref)
  for (v in names(RiboPilot:::REVIEW_VIEWS_DDL)) {
    expect_no_error(DBI::dbGetQuery(con, sprintf("SELECT ID, path, scaffold FROM %s LIMIT 1", v)))
  }
  u <- DBI::dbGetQuery(con, "SELECT has_copy, trim_end, check_switch FROM review_units_v ORDER BY path, scaffold")
  expect_equal(u$has_copy, rep(0L, 4))
  expect_equal(u$trim_end, c(6120L, 1442L, 5980L, 890L))
  expect_equal(u$check_switch, rep(1L, 4))

  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_x' WHERE path = 2 AND scaffold = 2")
  fixture_mark_done(con)
  d <- DBI::dbGetQuery(con, "SELECT v.check_switch, an.itsx_warnings FROM review_units_v v
    JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold
    ORDER BY v.path, v.scaffold")
  expect_equal(d$check_switch, rep(2L, 4))
  expect_equal(d$itsx_warnings, c("", "", "", "itsx_x"))

  before <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name")
  fixture_review_schema(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name"), before)
})
```

Append to the end of `tests/testthat/test-app-export.R` (`mk_export_project()` lives in this file, so its check goes here):

```r

test_that("the export fixture answers every review view", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  for (v in names(RiboPilot:::REVIEW_VIEWS_DDL)) {
    expect_no_error(DBI::dbGetQuery(con, sprintf("SELECT ID, path, scaffold FROM %s LIMIT 1", v)))
  }
  # locked and annotated, but no stored ITSx value: not done until fixture_mark_done()
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v"),
    data.frame(has_copy = 0L, length = 600L, check_switch = 1L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_annotations_v")$n, 3L)
})
```

- [ ] **Step 8: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "app-export")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-schema.R `[ FAIL 6 | WARN 0 | SKIP 0 | PASS 110 ]`: five "Expected `all(DBI::dbListFields(ref, t) %in% DBI::dbListFields(con, t))` to be TRUE" failures (assemblies, annotate, annotations, coverage, blast_screen lack columns) and `Error: no such table: blast_ref_alignment`. test-app-export.R: "the export fixture answers every review view" fails with `no such table: review_units_v` (and the other four views); every other test in the file passes.

- [ ] **Step 9: Add the fixture helpers**

In `tests/testthat/helper-app-fixture.R`:

Find (end of `fixture_fill_schema()`):

```r
  for (col in RiboPilot:::META_SAMPLE_COLS) add_col("samples", col, "TEXT")
  invisible()
}
```

Replace with:

```r
  for (col in RiboPilot:::META_SAMPLE_COLS) add_col("samples", col, "TEXT")
  fixture_review_schema(con)
  invisible()
}

# The tables the review layer reads, as a fresh new_db() makes them; built once.
fixture_ref <- new.env()
fixture_ref_schema <- function() {
  if (is.null(fixture_ref$sql)) {
    tbls <- c("assemblies", "annotate", "annotations", "coverage", "screen", "chimera",
      "chimera_windows", "blast_hits", "blast_screen", "blast_ref_alignment")
    con <- chimera_test_con()
    fixture_ref$sql <- stats::setNames(vapply(tbls, function(t) DBI::dbGetQuery(con,
      "SELECT sql FROM sqlite_master WHERE type = 'table' AND name = ?", params = list(t))$sql,
      character(1)), tbls)
    fixture_ref$cols <- stats::setNames(lapply(tbls, function(t)
      DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", t))[, c("name", "type")]), tbls)
    DBI::dbDisconnect(con)
  }
  fixture_ref
}

# Widen a hand-built fixture to those tables (creating missing ones), then add
# the review tables and views. Idempotent; run it after positional INSERTs.
fixture_review_schema <- function(con) {
  ref <- fixture_ref_schema()
  for (t in names(ref$sql)) {
    if (!DBI::dbExistsTable(con, t)) {
      DBI::dbExecute(con, ref$sql[[t]])
      next
    }
    miss <- ref$cols[[t]][!ref$cols[[t]]$name %in% DBI::dbListFields(con, t), ]
    for (i in seq_len(nrow(miss))) {
      DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", t, miss$name[i], miss$type[i]))
    }
  }
  RiboPilot:::review_apply_schema(con)
  invisible()
}

# Every contig annotated with a stored ITSx result, so the review layer reads it as done.
fixture_mark_done <- function(con) {
  DBI::dbExecute(con,
    "UPDATE annotate SET annotate_switch = 2, itsx_warnings = COALESCE(itsx_warnings, '')")
  invisible()
}
```

`chimera_test_con()` (helper-chimera.R) registers its temp folder on `fixture_ref_schema()`'s frame, so the connection is closed before that function returns. `ASSEMBLY_BACKUP_DDL` stays in `fixture_fill_schema()` until Task 22.

In `tests/testthat/test-batch-trim.R`, `mk_app_project()` now has these columns, so the ALTERs would fail with "duplicate column name":

Find:

```r
  on.exit(DBI::dbDisconnect(con))
  for (sql in c("ALTER TABLE assemblies ADD COLUMN edited INT",
    "ALTER TABLE annotations ADD COLUMN length INT", "ALTER TABLE annotations ADD COLUMN edited INT",
    "ALTER TABLE annotations ADD COLUMN type TEXT", "ALTER TABLE coverage ADD COLUMN clips TEXT")) {
    DBI::dbExecute(con, sql)
  }
```

Replace with:

```r
  on.exit(DBI::dbDisconnect(con))
```

In `tests/testthat/test-app-export.R`, the same duplicate-column failure; `add_trim_cols()` keeps its name so its three callers stay as they are:

Find:

```r
# trim_assembly() writes these; the export fixture was made without them
add_trim_cols <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "ALTER TABLE assemblies ADD COLUMN edited INT")
  DBI::dbExecute(con, "ALTER TABLE annotations ADD COLUMN length INT")
  DBI::dbExecute(con, "ALTER TABLE annotations ADD COLUMN edited INT")
  DBI::dbExecute(con, "ALTER TABLE coverage ADD COLUMN clips TEXT")
}
```

Replace with:

```r
# trim_assembly() writes these; fixture_review_schema() adds any the fixture lacks
add_trim_cols <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  fixture_review_schema(con)
}
```

- [ ] **Step 10: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "backwards-compat", "app-export", "batch-trim")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-schema.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 127 ]`, test-backwards-compat.R `PASS 38`, test-app-export.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 142 ]`, test-batch-trim.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 51 ]`.

Then every other file built on the app fixtures (the widened columns and the views must not change what they read; test-app-export.R:896 renames `samples`, which works because no view is broken):
Run: `Rscript -e 'devtools::load_all("."); for (f in c("app-table", "blast-provenance", "app-opts", "notes-partial", "seqview", "app-run", "specimen-app", "app-scaffold-uid")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` (PASS 47, 23, 63, 9, 17, 31, 99, 123).

Then the whole suite once: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: `FAIL 0`. On a heavily loaded machine the two `open_path` tests in test-app-utils.R can time out waiting for their marker file; rerun that file alone, where they pass.

- [ ] **Step 11: Commit**

```bash
git add R/db_schema.R R/init_db.R R/db_migrate.R tests/testthat/helper-app-fixture.R tests/testthat/test-review-schema.R tests/testthat/test-backwards-compat.R tests/testthat/test-batch-trim.R tests/testthat/test-app-export.R
git commit -m "Review tables, views, and schema step 4"
```


### Task 2: Review copy lifecycle

**Files:**
- Create: `R/review.R`
- Create: `tests/testthat/test-review-copy.R`
- Modify: `R/pipeline_steps.R` (`UNIT_EVIDENCE_TABLES` and its roxygen lines, :237-241)
- Test: `tests/testthat/test-review-copy.R`

**Interfaces:**
- Consumes (Task 1): `REVIEW_UNITS_DDL`, `REVIEW_MIRRORS` (named chr, original -> mirror: annotations, blast_hits, blast_screen, screen, chimera, chimera_windows, blast_ref_alignment), `REVIEW_TABLES` (the 8 names, review_units first), `ddl_cols(ddl)` (column names of a CREATE TABLE constant, in order), all in `R/db_schema.R`; the 8 review tables and the 5 views that `new_db()` and schema step 4 create, in particular the `review_units_v` fallback row (`has_copy = 0`, `trim_start = 1`, `trim_end = length`, `edited = 0`, derived `check_switch`); `current_schema_version()` = 5; test helpers `chimera_test_con()` and `mk_app_project()` (which runs `fixture_review_schema()`).
- Produces (`R/review.R`, all `@noRd`; call them as `RiboPilot:::fn` from tests and bin scripts):
  - `REVIEW_EFF_PARTIALS`: the SQL fragment `"COALESCE(partial_start_user, partial_start) AS partial_start, COALESCE(partial_stop_user, partial_stop) AS partial_stop"`, for explicit column lists on `review_annotations` / `review_annotations_v`.
  - `with_savepoint(con, code, name = "rp_review")`: runs `code` atomically and returns its value. No open transaction: `BEGIN IMMEDIATE` ... `COMMIT` (the write lock is taken before the first read). Inside an open transaction (for example `DBI::dbWithTransaction()` or another `with_savepoint()`): `SAVEPOINT <name>` ... `RELEASE <name>`. An error rolls back (`ROLLBACK`, or `ROLLBACK TO` then `RELEASE`) and propagates. Leaving `code` early with `return()` also rolls back, so a block's result must be its last value, never a `return()`. Code inside must not call `DBI::dbWithTransaction()` or dplyr `rows_*(in_place = TRUE)`; `DBI::dbAppendTable()` is fine.
  - `review_new_rev()`: character, the current time in milliseconds with every digit (`format(round(as.numeric(Sys.time()) * 1000), scientific = FALSE)`).
  - `review_has_copy(con, id, path, scaffold)`: logical(1).
  - `review_copy(con, id, path, scaffold)`: one `with_savepoint()`. Returns `invisible(FALSE)` when a copy exists (nothing written). Otherwise deletes the unit's rows in the 7 mirrors, inserts the `review_units` row from the unit's `review_units_v` row (`has_copy = 0`) with `rev = CAST(review_new_rev() AS INTEGER)`, `time_stamp = as.integer(Sys.time())`, `check_note` and `checked_at` NULL, then copies each original's rows into its mirror by column name (`intersect(dbListFields(orig), dbListFields(mirror))`; user partials NULL); `invisible(TRUE)`. Stops with `"No sequence for contig <ID>.<path>.<scaffold>."` (for example `No sequence for contig S1.1.2.`) when the unit has no `assemblies` row or no `annotate` row, and writes nothing.
  - `review_counts(con, id, path, scaffold)`: `list(n_features = <integer row count of review_annotations>, missing = paste(setdiff(c(RRNA_GENES, RRNA_SPACERS), elements), collapse = ","))`; `missing` is `""` when nothing is missing.
  - `review_invalidate(con, id, path, scaffold)`: one `with_savepoint()`. `review_units`: `check_switch = 1`, `rev = rev + 1`, `check_note`, `checked_at`, `n_flags`, `curate_flags`, `itsx_warnings`, `itsx_origin`, `itsx_scores` NULL, `n_features` and `missing` from `review_counts()`, `time_stamp` now; `chimeric`, `blast_*`, sequence, window, and `edited` kept. Deletes the unit's rows in review_blast_hits, review_blast_screen, review_screen, review_chimera, review_chimera_windows, review_ref_alignment. On review_annotations sets `status`, `flags`, `partial_start`, `partial_stop` NULL (user partials, coordinates, `model_coverage` kept). Returns `invisible(<new rev as text>)`. Stops with `"No review copy for contig <ID>.<path>.<scaffold>."` when the unit has no copy, and writes nothing.
  - `review_drop(con, id, path, scaffold)`: one `with_savepoint()` deleting the unit's rows in all 8 review tables (no error when there are none); `invisible(NULL)`.
  - `UNIT_EVIDENCE_TABLES` (`R/pipeline_steps.R`) ends with `"assembly_backup", REVIEW_TABLES)`, so `reset_unit_evidence()` and the `load_assemblies()` orphan sweep delete a unit's review copy. Tables only, never a view name.

This task adds the functions that make, invalidate, and delete a contig's review copy (SPEC S5.1, S5.2), plus `with_savepoint()`, which every later `review_*` writer wraps itself in so it is atomic on its own and also safe inside an open transaction such as the upgrade step's (S5.1, S11). `review_copy()` first deletes rows a discarded recheck left in the mirrors (S5.4) and takes its `review_units` row from the view's fallback row, so the derived check state has a single definition (S4.3); it copies mirror rows by column name because user projects are version 1 plus ALTER-appended columns and the hand-built fixtures hold columns in other orders (Review Focus 5). `review_invalidate()` keeps `chimeric` and the `blast_*` summary, because the recheck reads them "as copied" (S7.3 steps 2 and 6). Adding the review tables to `UNIT_EVIDENCE_TABLES` makes a re-assembly that changes or removes a contig delete its copy (S5.3, re-assembly row). Nothing in the app calls these functions yet.

- [ ] **Step 1: Write the failing lifecycle tests**

Create `tests/testthat/test-review-copy.R`:

```r
# Review copy lifecycle: first touch, invalidate, drop (SPEC S5.1-S5.4).

RC_SEQ <- strrep("ACGTTGCA", 25)

# Contigs S1.1.1 and S1.1.2 (200 bp), each with rows in every original a copy
# reads. Only columns a version 1 database has, so the upgrade test can use it.
rc_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:2) {
    p <- list(s)
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore,
        blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
      VALUES ('S1', 1, ?, 'linear', 200, ?, 0, 'A1', 'Danio rerio', 99.5, 88, 1e-50, 'core_nt', 7)",
      params = list(s, RC_SEQ))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, time_stamp)
      VALUES ('S1', 1, ?, 2, 0, 3, '5.8S,ITS2', 1, 2, 'truncated', ?, 9)",
      params = list(s, c("itsx_x", NA)[s]))
    for (f in list(list("18S", 1L, 60L, 1L, NA), list("ITS1", 61L, 100L, NA, NA),
                   list("28S", 101L, 200L, NA, 1L))) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length,
          strand, model_coverage, partial_start, partial_stop, tool, status, flags, edited, notes, time_stamp)
        VALUES ('S1', 1, ?, ?, 'rRNA', ?, ?, ?, '+', 95.5, ?, ?, 'barrnap', 'ok', 'f', 0, 'n', 9)",
        params = list(s, f[[1]], f[[2]], f[[3]], f[[3]] - f[[2]] + 1L, f[[4]], f[[5]]))
    }
    for (r in 1:2) {
      DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession, species,
          pident, qcovs, evalue, bitscore, aln_len, qstart, qend, sstart, send, slen, time_stamp)
        VALUES ('S1', 1, ?, ?, 1, ?, 'Danio rerio', 99, 90, 1e-40, 300, 180, 11, 190, 1, 180, 1800, 9)",
        params = list(s, r, paste0("H", r)))
    }
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, accession, species, pident, qcovs,
        time_stamp) VALUES ('S1', 1, ?, 'H1', 'Danio rerio', 99, 90, 9)", params = p)
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version,
        variant_class, variant_group, variant_identity, time_stamp)
      VALUES ('S1', 1, ?, 'on_target', 3, 'major', 'g1', 99.1, 9)", params = p)
    DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, reason, time_stamp)
      VALUES ('S1', 1, ?, 0.1, 'clean', 'none', 9)", params = p)
    for (w in c(1L, 101L)) {
      DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, accession, pident)
        VALUES ('S1', 1, ?, ?, ?, 'H1', 98)", params = list(s, w, w + 99L))
    }
    DBI::dbExecute(con, "INSERT INTO blast_ref_alignment (ID, path, scaffold, accession, element, pid,
        aln_cov, aln_len, time_stamp) VALUES ('S1', 1, ?, 'R1', '18S', 98, 100, 60, 9)", params = p)
  }
}

# One row in each mirror, as a discarded recheck or an edit would leave them.
rc_mirror_rows <- function(con, s) {
  p <- list(s)
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, ?, '18S')",
    params = p)
  DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank) VALUES ('S1', 1, ?, 9)",
    params = p)
  for (t in c("review_blast_screen", "review_screen", "review_chimera")) {
    DBI::dbExecute(con, sprintf("INSERT INTO %s (ID, path, scaffold) VALUES ('S1', 1, ?)", t), params = p)
  }
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start) VALUES ('S1', 1, ?, 7)",
    params = p)
  DBI::dbExecute(con, "INSERT INTO review_ref_alignment (ID, path, scaffold, accession, element)
    VALUES ('S1', 1, ?, 'R9', '28S')", params = p)
}

# Rows of a table, by column name, in a stable order.
rc_rows <- function(con, tbl, cols = DBI::dbListFields(con, tbl), path = 1L, scaffold = NULL) {
  ord <- paste(cols, collapse = ", ")
  if (is.null(scaffold)) {
    return(DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", ord, tbl, ord)))
  }
  DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s WHERE ID = 'S1' AND path = ? AND scaffold = ? ORDER BY %s",
    ord, tbl, ord), params = list(as.integer(path), as.integer(scaffold)))
}

rc_counts <- function(con, scaffold, path = 1L) {
  vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = ? AND scaffold = ?", t),
    params = list(as.integer(path), as.integer(scaffold)))$n, integer(1))
}

rc_originals <- function(con) {
  lapply(stats::setNames(nm = c("assemblies", "annotate", names(RiboPilot:::REVIEW_MIRRORS))),
    function(t) rc_rows(con, t))
}

rc_unit <- function(con, scaffold) {
  DBI::dbGetQuery(con, "SELECT sequence, length, trim_start, trim_end, edited, check_switch, check_note,
      CAST(rev AS TEXT) AS rev, checked_at, n_features, missing, chimeric, n_flags, curate_flags,
      itsx_warnings, itsx_origin, itsx_scores, blast_accession, blast_species, blast_pident, blast_qcovs,
      blast_evalue, blast_source, time_stamp
    FROM review_units WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(as.integer(scaffold)))
}

rc_now_ms <- function() as.numeric(Sys.time()) * 1000

# Copy, invalidate, and drop on a seeded database (new_db or an upgraded one).
rc_lifecycle <- function(con) {
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M', itsx_scores = 'M:50'")
  orig <- rc_originals(con)
  mirrors <- RiboPilot:::REVIEW_MIRRORS

  # first touch: the unit's originals, the view's check state, a fresh rev
  expect_false(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  t0 <- rc_now_ms()
  expect_true(RiboPilot:::review_copy(con, "S1", 1L, 1L))
  t1 <- rc_now_ms()
  expect_true(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  u <- rc_unit(con, 1)
  expect_equal(u[, c("sequence", "length", "trim_start", "trim_end", "edited", "check_switch",
      "check_note", "checked_at", "n_features", "missing", "chimeric", "n_flags", "curate_flags",
      "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession", "blast_species", "blast_pident",
      "blast_qcovs", "blast_evalue", "blast_source")],
    data.frame(sequence = RC_SEQ, length = 200L, trim_start = 1L, trim_end = 200L, edited = 0L,
      check_switch = 2L, check_note = NA_character_, checked_at = NA_integer_, n_features = 3L,
      missing = "5.8S,ITS2", chimeric = 1L, n_flags = 2L, curate_flags = "truncated",
      itsx_warnings = "itsx_x", itsx_origin = "M", itsx_scores = "M:50", blast_accession = "A1",
      blast_species = "Danio rerio", blast_pident = 99.5, blast_qcovs = 88, blast_evalue = 1e-50,
      blast_source = "core_nt"))
  expect_match(u$rev, "^[0-9]{13}$")
  expect_gte(as.numeric(u$rev), floor(t0))
  expect_lte(as.numeric(u$rev), ceiling(t1))
  expect_gte(u$time_stamp, floor(t0 / 1000))
  expect_lte(u$time_stamp, ceiling(t1 / 1000))
  for (o in names(mirrors)) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, mirrors[[o]], cols, scaffold = 1), rc_rows(con, o, cols, scaffold = 1), info = o)
  }
  expect_equal(nrow(rc_rows(con, "review_annotations", scaffold = 1)), 3L)
  expect_true(all(is.na(rc_rows(con, "review_annotations", c("partial_start_user", "partial_stop_user")))))
  expect_equal(unname(rc_counts(con, 2)), rep(0L, 8))
  expect_identical(rc_originals(con), orig)

  # a second touch changes nothing
  expect_false(RiboPilot:::review_copy(con, "S1", 1L, 1L))
  expect_equal(rc_unit(con, 1), u)

  # the user's partial wins over curation's on the review layer
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 0 WHERE scaffold = 1 AND element = '18S'")
  eff <- DBI::dbGetQuery(con, paste("SELECT element,", RiboPilot:::REVIEW_EFF_PARTIALS,
    "FROM review_annotations_v WHERE scaffold = 1 ORDER BY start"))
  expect_equal(eff$partial_start, c(0L, NA, NA))
  expect_equal(eff$partial_stop, c(NA, NA, 1L))

  # rows a discarded recheck left for a contig without a copy are replaced;
  # done without a stored ITSx value starts as not run
  rc_mirror_rows(con, 2)
  expect_true(RiboPilot:::review_copy(con, "S1", 1L, 2L))
  for (o in names(mirrors)) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, mirrors[[o]], cols, scaffold = 2), rc_rows(con, o, cols, scaffold = 2), info = o)
  }
  expect_equal(rc_unit(con, 2)$check_switch, 1L)
  expect_true(is.na(rc_unit(con, 2)$itsx_warnings))

  # invalidate after a finished recheck and an edit that removed ITS1
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, check_note = 'old', checked_at = 5,
    time_stamp = 1 WHERE scaffold = 1")
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE scaffold = 1 AND element = 'ITS1'")
  other <- lapply(stats::setNames(nm = RiboPilot:::REVIEW_TABLES[-1]), function(t) rc_rows(con, t, scaffold = 2))
  other_unit <- rc_unit(con, 2)
  rev0 <- u$rev
  t0 <- rc_now_ms()
  rev1 <- RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  t1 <- rc_now_ms()
  v <- rc_unit(con, 1)
  expect_equal(rev1, v$rev)
  expect_equal(DBI::dbGetQuery(con, "SELECT CAST(rev - 1 AS TEXT) AS r FROM review_units WHERE scaffold = 1")$r,
    rev0)
  expect_equal(v[, c("check_switch", "check_note", "checked_at", "n_features", "missing", "n_flags",
      "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores")],
    data.frame(check_switch = 1L, check_note = NA_character_, checked_at = NA_integer_, n_features = 2L,
      missing = "5.8S,ITS1,ITS2", n_flags = NA_integer_, curate_flags = NA_character_,
      itsx_warnings = NA_character_, itsx_origin = NA_character_, itsx_scores = NA_character_))
  # kept as copied: sequence, window, the ITSx chimera flag, and the BLAST summary
  keep <- c("sequence", "length", "trim_start", "trim_end", "edited", "chimeric", "blast_accession",
    "blast_species", "blast_pident", "blast_qcovs", "blast_evalue", "blast_source")
  expect_equal(v[, keep], u[, keep])
  expect_gte(v$time_stamp, floor(t0 / 1000))
  expect_lte(v$time_stamp, ceiling(t1 / 1000))
  ra <- rc_rows(con, "review_annotations", c("element", "start", "stop", "model_coverage", "tool", "status",
    "flags", "partial_start", "partial_stop", "partial_start_user"), scaffold = 1)
  expect_equal(ra$element, c("18S", "28S"))
  expect_equal(ra$start, c(1L, 101L))
  expect_equal(ra$model_coverage, c(95.5, 95.5))
  expect_equal(ra$tool, c("barrnap", "barrnap"))
  expect_true(all(is.na(ra[, c("status", "flags", "partial_start", "partial_stop")])))
  expect_equal(ra$partial_start_user, c(0L, NA))
  expect_equal(unname(rc_counts(con, 1)), c(1L, 2L, rep(0L, 6)))
  # the other contig's copy and every original are untouched
  expect_identical(lapply(stats::setNames(nm = RiboPilot:::REVIEW_TABLES[-1]), function(t)
    rc_rows(con, t, scaffold = 2)), other)
  expect_equal(rc_unit(con, 2), other_unit)
  expect_identical(rc_originals(con), orig)
  rev2 <- RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT CAST(rev - 2 AS TEXT) AS r FROM review_units
    WHERE scaffold = 1 AND rev = CAST(? AS INTEGER)", params = list(rev2))$r, rev0)

  # drop: the contig reads its originals again; the other copy stays
  RiboPilot:::review_drop(con, "S1", 1L, 1L)
  expect_false(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  expect_equal(unname(rc_counts(con, 1)), rep(0L, 8))
  expect_equal(unname(rc_counts(con, 2)), c(1L, 3L, 2L, 1L, 1L, 1L, 2L, 1L))
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1"), data.frame(has_copy = 0L, length = 200L, check_switch = 2L))
  expect_identical(rc_originals(con), orig)
  expect_no_error(RiboPilot:::review_drop(con, "S1", 1L, 1L))
}

test_that("with_savepoint() commits, rolls back, and nests inside an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  ws <- RiboPilot:::with_savepoint
  DBI::dbExecute(con, "CREATE TABLE t (x INTEGER)")
  xs <- function() DBI::dbGetQuery(con, "SELECT x FROM t ORDER BY x")$x

  expect_equal(ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (1)")
    "value"
  }), "value")
  expect_error(ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (2)")
    stop("boom")
  }), "boom")
  expect_equal(xs(), 1L)
  expect_false(RSQLite::sqliteIsTransacting(con))

  # inside DBI::dbWithTransaction: a failed block undoes only its own writes
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (3)")
    expect_error(ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (4)")
      stop("inner")
    }), "inner")
    ws(con, ws(con, DBI::dbExecute(con, "INSERT INTO t VALUES (5)")))
    expect_true(RSQLite::sqliteIsTransacting(con))
  })
  expect_equal(xs(), c(1L, 3L, 5L))

  # a nested failure caught inside an outer block
  ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (6)")
    try(ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (7)")
      stop("nested")
    }), silent = TRUE)
  })
  expect_equal(xs(), c(1L, 3L, 5L, 6L))

  # leaving the block early rolls back like an error
  f <- function() {
    ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (8)")
      return("early")
    })
    "late"
  }
  expect_equal(f(), "early")
  expect_equal(xs(), c(1L, 3L, 5L, 6L))
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("with_savepoint() takes the write lock before its first statement", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con), add = TRUE)
  other <- DBI::dbConnect(RSQLite::SQLite(), con@dbname)
  on.exit(DBI::dbDisconnect(other), add = TRUE)
  DBI::dbExecute(other, "PRAGMA busy_timeout = 0")
  DBI::dbExecute(con, "CREATE TABLE t (x INTEGER)")
  got <- RiboPilot:::with_savepoint(con, {
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM t")
    tryCatch({
      DBI::dbExecute(other, "INSERT INTO t VALUES (1)")
      "written"
    }, error = conditionMessage)
  })
  expect_match(got, "locked")
  expect_equal(DBI::dbExecute(other, "INSERT INTO t VALUES (1)"), 1)
})

test_that("review_new_rev() is the time in milliseconds with every digit", {
  withr::local_options(digits = 3)
  t0 <- rc_now_ms()
  r <- RiboPilot:::review_new_rev()
  t1 <- rc_now_ms()
  expect_type(r, "character")
  expect_match(r, "^[0-9]{13}$")
  expect_gte(as.numeric(r), floor(t0))
  expect_lte(as.numeric(r), ceiling(t1))
})

test_that("review_counts() counts the copy's features as reconcile_one does", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 1L),
    list(n_features = 0L, missing = "18S,5.8S,28S,ITS1,ITS2"))
  for (e in c("28S", "ITS1", "18S", "5.8S", "ITS2")) {
    DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, 1, ?)",
      params = list(e))
  }
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, 2, '18S')")
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 1L), list(n_features = 5L, missing = ""))
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 2L),
    list(n_features = 1L, missing = "5.8S,28S,ITS1,ITS2"))
})

test_that("copy, invalidate, and drop on a new project", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  rc_lifecycle(con)
})

test_that("copy, invalidate, and drop on a version 1 database holding data, upgraded through every step", {
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  rc_seed(con)
  DBI::dbDisconnect(con)

  suppressMessages(RiboPilot::backwards_compatibility(dir))
  con <- RiboPilot:::open_project_db(db)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(max(DBI::dbGetQuery(con, "SELECT schema_version FROM schema_history")$schema_version),
    RiboPilot:::current_schema_version())
  rc_lifecycle(con)
})

test_that("review_copy() copies by column name on a hand-built database", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length,
      strand, status, flags, partial_start, edited, time_stamp)
    VALUES ('S1', 2, 1, '18S', 'rRNA', 5, 40, 36, '-', 'ok', 'f', 1, 0, 9)")
  DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, accession, species, pident, time_stamp)
    VALUES ('S1', 2, 1, 'H1', 'Danio rerio', 97.5, 9)")
  RiboPilot:::review_copy(con, "S1", 2L, 1L)
  for (o in c("annotations", "blast_screen")) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, RiboPilot:::REVIEW_MIRRORS[[o]], cols, path = 2, scaffold = 1),
      rc_rows(con, o, cols, path = 2, scaffold = 1), info = o)
  }
  expect_equal(DBI::dbGetQuery(con, "SELECT path, scaffold, sequence, length, trim_end, check_switch,
      n_features, missing, blast_accession, blast_qcovs FROM review_units"),
    data.frame(path = 2L, scaffold = 1L, sequence = "ACGT", length = 5980L, trim_end = 5980L,
      check_switch = 1L, n_features = 3L, missing = "", blast_accession = "A1", blast_qcovs = 88))

  RiboPilot:::review_invalidate(con, "S1", 2L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT n_features, missing FROM review_units"),
    data.frame(n_features = 1L, missing = "5.8S,28S,ITS1,ITS2"))
  expect_equal(unname(rc_counts(con, 1, path = 2)), c(1L, 1L, rep(0L, 6)))
  RiboPilot:::review_drop(con, "S1", 2L, 1L)
  expect_equal(unname(rc_counts(con, 1, path = 2)), rep(0L, 8))
})

test_that("review_copy() refuses a contig with no annotate or assemblies row and writes nothing", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  DBI::dbExecute(con, "DELETE FROM annotate WHERE scaffold = 2")
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 2L), "No sequence for contig S1.1.2.", fixed = TRUE)
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 9L), "No sequence for contig S1.1.9.", fixed = TRUE)
  expect_equal(unname(rc_counts(con, 2)), rep(0L, 8))
  expect_error(RiboPilot:::review_invalidate(con, "S1", 1L, 1L), "No review copy for contig S1.1.1.",
    fixed = TRUE)
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("a failure part way through review_copy() leaves no copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  DBI::dbExecute(con, "CREATE TRIGGER rc_boom BEFORE INSERT ON review_ref_alignment
    BEGIN SELECT RAISE(ABORT, 'rc_boom'); END")
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 1L), "rc_boom")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, 8))
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("review_copy, review_invalidate, and review_drop join an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  expect_error(DBI::dbWithTransaction(con, {
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    stop("outer failed")
  }), "outer failed")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, 8))

  DBI::dbWithTransaction(con, {
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    RiboPilot:::review_copy(con, "S1", 1L, 2L)
    RiboPilot:::review_drop(con, "S1", 1L, 2L)
  })
  expect_equal(rc_unit(con, 1)$check_switch, 1L)
  expect_equal(unname(rc_counts(con, 1)), c(1L, 3L, rep(0L, 6)))
  expect_equal(unname(rc_counts(con, 2)), rep(0L, 8))
})
```

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-copy.R")'`
Expected: FAIL. Every test errors with `object 'with_savepoint' not found`, `object 'review_new_rev' not found`, `object 'review_counts' not found`, `object 'review_has_copy' not found`, or `object 'review_copy' not found`, and testthat stops with "Maximum number of failures exceeded; quitting." (with `testthat::set_max_fails(Inf)` the count is `[ FAIL 10 | WARN 0 | SKIP 0 | PASS 1 ]`).

- [ ] **Step 3: Create R/review.R**

Create `R/review.R`:

```r
# Review copies: the per-contig working copy that Annotation Review edits and
# the recheck scores, so the originals keep the main pipeline's results.
# Design: dev/specs/2026-10-08-contig-recheck-design.md (S4, S5).

#' Effective partial ends on the review layer: the user's choice, else curation's
#' @noRd
REVIEW_EFF_PARTIALS <- paste("COALESCE(partial_start_user, partial_start) AS partial_start,",
  "COALESCE(partial_stop_user, partial_stop) AS partial_stop")

#' Run code atomically, inside or outside an open transaction
#'
#' With no transaction open it runs BEGIN IMMEDIATE ... COMMIT, so the write
#' lock is taken before the first read and busy_timeout applies. Inside an open
#' transaction it runs SAVEPOINT ... RELEASE. An error, or leaving code early
#' (return() inside code, an interrupt), rolls back and leaves no transaction or
#' savepoint open; the error propagates. Code inside must not call
#' DBI::dbWithTransaction(). Returns the value of code.
#' @noRd
with_savepoint <- function(con, code, name = "rp_review") {
  outer <- !RSQLite::sqliteIsTransacting(con)
  DBI::dbExecute(con, if (outer) "BEGIN IMMEDIATE" else paste0("SAVEPOINT ", name))
  ok <- FALSE
  on.exit(if (!ok) {
    if (outer) {
      try(DBI::dbExecute(con, "ROLLBACK"), silent = TRUE)
    } else {
      try(DBI::dbExecute(con, paste0("ROLLBACK TO ", name)), silent = TRUE)
      try(DBI::dbExecute(con, paste0("RELEASE ", name)), silent = TRUE)
    }
  })
  res <- force(code)
  DBI::dbExecute(con, if (outer) "COMMIT" else paste0("RELEASE ", name))
  ok <- TRUE
  res
}

#' A new review_units.rev: the current time in milliseconds, as text
#'
#' rev does not fit an R integer; carry it as text and compare it in SQL
#' (rev = CAST(? AS INTEGER)).
#' @noRd
review_new_rev <- function() format(round(as.numeric(Sys.time()) * 1000), scientific = FALSE)

#' @noRd
review_has_copy <- function(con, id, path, scaffold) {
  nrow(DBI::dbGetQuery(con,
    "SELECT 1 AS x FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))) > 0
}

#' Make a contig's review copy from its originals, unless it has one (S5.1)
#'
#' Rows a discarded recheck left in the mirrors are deleted first. The
#' review_units row comes from the review_units_v fallback row, so the copy
#' starts with the check state the view showed. Mirror rows are copied by column
#' name, since column order differs between new_db(), upgraded databases, and
#' hand-built ones.
#' @return invisible TRUE when a copy was made, FALSE when one existed.
#' @noRd
review_copy <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  cols <- paste(setdiff(ddl_cols(REVIEW_UNITS_DDL), c("check_note", "rev", "checked_at", "time_stamp")),
    collapse = ", ")
  invisible(with_savepoint(con, {
    made <- !review_has_copy(con, id, path, scaffold)
    if (made) {
      for (m in REVIEW_MIRRORS) {
        DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", m),
          params = p)
      }
      n <- DBI::dbExecute(con, sprintf(
        "INSERT INTO review_units (%s, rev, time_stamp)
         SELECT %s, CAST(? AS INTEGER), ? FROM review_units_v
         WHERE ID = ? AND path = ? AND scaffold = ? AND has_copy = 0", cols, cols),
        params = c(list(review_new_rev(), as.integer(Sys.time())), p))
      if (n == 0) stop(sprintf("No sequence for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
      for (orig in names(REVIEW_MIRRORS)) {
        mc <- paste(intersect(DBI::dbListFields(con, orig), DBI::dbListFields(con, REVIEW_MIRRORS[[orig]])),
          collapse = ", ")
        DBI::dbExecute(con, sprintf(
          "INSERT INTO %s (%s) SELECT %s FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
          REVIEW_MIRRORS[[orig]], mc, mc, orig), params = p)
      }
    }
    made
  }))
}

#' Feature count and missing elements of a copy, as reconcile_one counts them
#' @noRd
review_counts <- function(con, id, path, scaffold) {
  el <- DBI::dbGetQuery(con,
    "SELECT element FROM review_annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))$element
  list(n_features = length(el),
    missing = paste(setdiff(c(RRNA_GENES, RRNA_SPACERS), el), collapse = ","))
}

#' Mark a copy's checks not run after an edit (S5.2)
#'
#' Raises rev, clears the check results and the curation-owned feature columns,
#' and deletes the copy's check rows. Counts are recomputed, not cleared;
#' chimeric and blast_* stay as copied until the recheck rewrites them.
#' @return invisible new rev, as text.
#' @noRd
review_invalidate <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  invisible(with_savepoint(con, {
    cnt <- review_counts(con, id, path, scaffold)
    n <- DBI::dbExecute(con,
      "UPDATE review_units SET check_switch = 1, rev = rev + 1, check_note = NULL, checked_at = NULL,
         n_features = ?, missing = ?, n_flags = NULL, curate_flags = NULL, itsx_warnings = NULL,
         itsx_origin = NULL, itsx_scores = NULL, time_stamp = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(cnt$n_features, cnt$missing, as.integer(Sys.time())), p))
    if (n == 0) stop(sprintf("No review copy for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
    for (m in REVIEW_MIRRORS[c("blast_hits", "blast_screen", "screen", "chimera", "chimera_windows",
                               "blast_ref_alignment")]) {
      DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", m),
        params = p)
    }
    DBI::dbExecute(con,
      "UPDATE review_annotations SET status = NULL, flags = NULL, partial_start = NULL, partial_stop = NULL
       WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    DBI::dbGetQuery(con,
      "SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?",
      params = p)$rev
  }))
}

#' Delete a contig's review copy, so it reads its originals again (S5.2 Restore)
#' @noRd
review_drop <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  with_savepoint(con, for (t in REVIEW_TABLES) {
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", t), params = p)
  })
  invisible()
}
```

`REVIEW_UNITS_DDL`, `REVIEW_MIRRORS`, `REVIEW_TABLES`, and `ddl_cols()` come from `R/db_schema.R` (Task 1), which R sources before `R/review.R` (DESCRIPTION has no Collate field, so files load in name order). The `review_units` INSERT names every column except `check_note`, `rev`, `checked_at`, and `time_stamp`, and selects the same names from `review_units_v`, so the derived `check_switch` and `trim_end` come from the view. Check the file parses: `Rscript -e 'invisible(parse("R/review.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-copy.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 154 ]`.

- [ ] **Step 5: Write the failing re-assembly test**

Append to the end of `tests/testthat/test-review-copy.R`:

```r

test_that("a re-assembly deletes the copies of changed and vanished contigs", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  out <- withr::local_tempdir()
  sq <- lapply(1:4, function(k) withr::with_seed(k, paste(sample(c("A", "C", "G", "T"), 1200, TRUE),
    collapse = "")))
  fa <- function(x) writeLines(unlist(lapply(names(x), function(k)
    c(sprintf(">S1.1.%s linear c%s", k, k), x[[k]]))), file.path(out, "S1_assembly_1.fasta"))
  fa(list("1" = sq[[1]], "2" = sq[[2]], "3" = sq[[3]]))
  RiboPilot:::load_assemblies(con, out, "S1")
  for (s in 1:3) {
    RiboPilot:::review_copy(con, "S1", 1L, s)
    rc_mirror_rows(con, s)
  }
  expect_equal(unname(rc_counts(con, 1)), rep(1L, 8))

  # scaffold 1 comes back with another sequence, 2 unchanged, 3 gone
  fa(list("1" = sq[[4]], "2" = sq[[2]]))
  RiboPilot:::load_assemblies(con, out, "S1")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, 8))
  expect_equal(unname(rc_counts(con, 2)), rep(1L, 8))
  expect_equal(unname(rc_counts(con, 3)), rep(0L, 8))

  RiboPilot:::reset_unit_evidence(con, "S1", 1L, 2L)
  expect_equal(unname(rc_counts(con, 2)), rep(0L, 8))
})
```

- [ ] **Step 6: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-copy.R")'`
Expected: FAIL, `[ FAIL 3 | WARN 0 | SKIP 0 | PASS 156 ]`. The three failures are in "a re-assembly deletes the copies of changed and vanished contigs": `unname(rc_counts(con, 1))`, `unname(rc_counts(con, 3))`, and the last `unname(rc_counts(con, 2))` are not equal to `rep(0L, 8)` (the copies survive, because `reset_unit_evidence()` and the orphan sweep do not know the review tables yet).

- [ ] **Step 7: Add the review tables to UNIT_EVIDENCE_TABLES (R/pipeline_steps.R)**

Find:

```r
#' Per-unit tables whose rows describe one assembled sequence
#' @noRd
UNIT_EVIDENCE_TABLES <- c("annotations", "coverage", "screen", "chimera", "chimera_windows",
  "blast_hits", "blast_screen", "blast_searched", "blast_ref_candidates",
  "blast_ref_alignment", "assembly_backup")
```

Replace with:

```r
#' Per-unit tables whose rows describe one assembled sequence, review copies included
#'
#' Tables only, never a view: callers DELETE from each name.
#' @noRd
UNIT_EVIDENCE_TABLES <- c("annotations", "coverage", "screen", "chimera", "chimera_windows",
  "blast_hits", "blast_screen", "blast_searched", "blast_ref_candidates",
  "blast_ref_alignment", "assembly_backup", REVIEW_TABLES)
```

Both users, `reset_unit_evidence()` and the orphan sweep in `load_assemblies()`, intersect this list with `DBI::dbListTables(con)`, so hand-built test databases without the review tables keep working. `dbListTables()` also returns views; the intersect never picks one because no view name is in the list. Check: `Rscript -e 'invisible(parse("R/pipeline_steps.R"))'`.

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-copy.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 159 ]`.

Then the existing files that run `load_assemblies()` or `reset_unit_evidence()`, or check the schema and the sources:
Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "backwards-compat", "load-assemblies", "sample-state", "write-atomicity", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` (PASS 127, 38, 119, 24, 22, 1, 237).

Then the bin-script and Nextflow files that call `load_assemblies()` through `rp_score.R` / `rp_find_rrna.R` (they skip without `NOT_CRAN=true`):
Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("pipeline-drivers", "nextflow", "nf-main-run")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` and `SKIP 0` (PASS 17, 299, 33; test-nf-main-run.R runs the real main.nf with fake tools and takes several minutes).

- [ ] **Step 9: Commit**

```bash
git add R/review.R R/pipeline_steps.R tests/testthat/test-review-copy.R
git commit -m "Review copy lifecycle: copy, invalidate, drop"
```


### Task 3: Coverage window and mapped BLAST hits

**Files:**
- Create: `tests/testthat/test-review-window.R`
- Modify: `R/review.R` (append after `review_drop()`, the end of the file Task 2 creates)
- Test: `tests/testthat/test-review-window.R`

**Interfaces:**
- Consumes: Task 1: `REVIEW_UNITS_DDL`, `ddl_cols(ddl)` (`R/db_schema.R`), the `review_units` and `review_blast_hits` tables that `new_db()` creates. Task 2: `review_copy(con, id, path, scaffold)` (tests only). Existing: `BLAST_HITS_DDL` (`R/db_schema.R`), `parse_depths(depths_blob)` (`R/coverage.R`), `write_coverage_summary(con, id, path, scaffold, dep)` (`R/pipeline_steps.R`, tests only), `chimera_test_con()` (`tests/testthat/helper-chimera.R`).
- Produces (`R/review.R`, all `@noRd`; call them as `RiboPilot:::fn` from tests and bin scripts):
  - `review_coverage(con, id, path, scaffold)`: always `list(depths = , clips = )`, numeric vectors in copy coordinates. Reads the unit's `coverage` row with `SELECT *` (fixtures lack `clips`) and parses each track with `parse_depths()`. Without a copy: the whole original tracks. With a copy: each track cut to `[trim_start, trim_end]` of its `review_units` row; a track shorter than the window is cut, never padded. An element is `NULL` when the coverage row or the track (column or value) is missing, or nothing of the track lies in the window. Never writes.
  - `review_coverage_summary(dep)`: `list(mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2), min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)))`, the values `write_coverage_summary()` stores; `NULL` for an empty or `NULL` track.
  - `review_map_hits(con, id, path, scaffold)`: stops with `"No review copy for contig <ID>.<path>.<scaffold>."` when the unit has no copy. Otherwise a data.frame with the columns `ddl_cols(BLAST_HITS_DDL)` in that order (the `review_blast_hits` columns), ordered by `hit_rank, hsp`: the original `blast_hits` HSPs with `qend >= trim_start` and `qstart <= trim_end`, clipped to the window and shifted by `trim_start - 1`. A clipped HSP gets `aln_len` = its clipped span and `sstart = send = NA`; unclipped HSPs keep `aln_len`, `sstart`, `send`. `hit_rank`, `hsp`, and every other column keep their stored values, so a subject or HSP clipped away leaves a gap. `qcovs` per subject (`hit_rank`) = `floor(100 * <length of the union of its clipped spans> / (trim_end - trim_start + 1) + 0.5)`. HSPs without query coordinates are left out. Zero rows (same columns) when nothing overlaps or the original has no hits. Never writes; the caller (Task 10, `recheck_hits()` map mode) deletes the unit's `review_blast_hits` rows and appends this frame.

A trim only removes sequence from the contig's ends, so the original's per-position data maps onto a review copy with an offset of `trim_start - 1` (SPEC S4.5). This task adds the two readers that do that mapping: `review_coverage()` gives the copy its read-depth and soft-clip tracks (the recheck's chimera check, S7.3 step 10, and the Annotation Review viewers, S8, read it instead of `coverage`), and `review_map_hits()` gives the copy BLAST hits when the recheck does no fresh search (S7.3 step 2, `map` mode). Coverage is never copied, so a trim no longer loses the clip track. On an untrimmed copy the mapped hits equal the original rows, and the recomputed `qcovs` reproduces every stored value (checked on all 1,333 HSPs of a real test project). Nothing calls these functions yet, and they never write the originals (S4.1).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-review-window.R`:

```r
# Coverage window and mapped BLAST hits of a review copy (SPEC S4.5).

# Kept contigs S1.1.<scaffold> of `len` bp, done in the main pipeline.
rw_units <- function(con, scaffolds, len) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in scaffolds) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', ?, ?, 0)", params = list(s, len, strrep("ACGTA", len / 5)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, '')", params = list(s))
  }
}

# Make the contig's copy and set its kept window, in original coordinates.
rw_window <- function(con, scaffold, from, to) {
  RiboPilot:::review_copy(con, "S1", 1L, scaffold)
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = ?, trim_end = ?, length = ?
    WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(from, to, to - from + 1L, scaffold))
}

rw_track <- function(x) if (is.null(x)) NA_character_ else paste(x, collapse = ",")

test_that("review_coverage() cuts the original tracks to the copy's window", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:4, 10L)
  cov <- function(s, depths, clips) DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips)
    VALUES ('S1', 1, ?, ?, ?)", params = list(s, rw_track(depths), rw_track(clips)))
  cov(1, 11:20, 0:9)
  cov(2, 31:35, 40:49)
  cov(3, 51:60, NULL)
  rc <- function(s) RiboPilot:::review_coverage(con, "S1", 1L, s)

  rw_window(con, 2, 3L, 7L)
  rw_window(con, 3, 2L, 9L)
  # no copy: the whole original tracks
  expect_equal(rc(1), list(depths = as.numeric(11:20), clips = as.numeric(0:9)))
  # no coverage row
  expect_equal(rc(4), list(depths = NULL, clips = NULL))

  rw_window(con, 1, 3L, 7L)
  expect_equal(rc(1), list(depths = as.numeric(13:17), clips = as.numeric(2:6)))
  # a depth track shorter than the contig is cut, never padded
  expect_equal(rc(2), list(depths = as.numeric(33:35), clips = as.numeric(42:46)))
  # clips lost by an old trim
  expect_equal(rc(3), list(depths = as.numeric(52:59), clips = NULL))
  # a window that starts past the end of a track
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = 6, trim_end = 9 WHERE scaffold = 2")
  expect_equal(rc(2), list(depths = NULL, clips = as.numeric(45:48)))
  rw_window(con, 4, 1L, 5L)
  expect_equal(rc(4), list(depths = NULL, clips = NULL))

  expect_equal(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE scaffold = 1"),
    data.frame(depths = rw_track(11:20), clips = rw_track(0:9)))
})

test_that("review_coverage() reads a coverage table without a clips column", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage (ID TEXT, path INT, scaffold INT, mean_depth REAL, depths TEXT)")
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbExecute(con, "INSERT INTO coverage VALUES ('S1', 1, 1, 3, '1,2,3,4,5')")
  expect_equal(RiboPilot:::review_coverage(con, "S1", 1L, 1L), list(depths = as.numeric(1:5), clips = NULL))
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, trim_start, trim_end, rev)
    VALUES ('S1', 1, 1, 2, 4, 1)")
  expect_equal(RiboPilot:::review_coverage(con, "S1", 1L, 1L), list(depths = as.numeric(2:4), clips = NULL))
})

test_that("review_coverage_summary() rounds as write_coverage_summary() stores", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  dep <- c(2, 9, 4, 4, 1, 3)
  expect_identical(RiboPilot:::review_coverage_summary(dep),
    list(mean_depth = 3.83, median_depth = 3.5, min_depth = 1L, max_depth = 9L))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold) VALUES ('S1', 1, 1)")
  RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep)
  expect_equal(as.list(DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, min_depth, max_depth FROM coverage")),
    RiboPilot:::review_coverage_summary(dep))
  expect_null(RiboPilot:::review_coverage_summary(NULL))
  expect_null(RiboPilot:::review_coverage_summary(numeric(0)))
})

# HSPs of contig S1.1.1 (300 bp). qcovs is stored as BLAST computes it on the
# whole contig: the union of a subject's HSPs over 300, rounded. Subject 7 has
# no query coordinates (an old hits file).
RW_HITS <- data.frame(
  hit_rank = c(1L, 1L, 2L, 3L, 3L, 3L, 4L, 4L, 5L, 6L, 6L, 7L),
  hsp = c(1L, 2L, 1L, 1L, 2L, 3L, 1L, 2L, 1L, 1L, 2L, 1L),
  qstart = c(1L, 61L, 1L, 241L, 251L, 51L, 230L, 200L, 250L, 40L, 30L, NA),
  qend = c(300L, 90L, 50L, 300L, 280L, 65L, 260L, 250L, 300L, 60L, 51L, NA),
  sstart = c(1L, 500L, 1L, 1L, 30L, 1000L, 9L, 77L, 1L, 5L, 2L, NA),
  send = c(300L, 470L, 50L, 60L, 59L, 1014L, 39L, 127L, 51L, 25L, 23L, NA),
  aln_len = c(300L, 31L, 50L, 60L, 30L, 15L, 31L, 51L, 51L, 21L, 22L, 99L),
  qcovs = c(100, 100, 17, 25, 25, 25, 20, 20, 17, 10, 10, 0))

rw_hits <- function(con) {
  h <- cbind(data.frame(ID = "S1", path = 1L, scaffold = 1L), RW_HITS,
    accession = paste0("H", RW_HITS$hit_rank), species = "Danio rerio", hit_class = "rdna",
    pident = 90 + seq_len(nrow(RW_HITS)) / 10, bitscore = 1000 - seq_len(nrow(RW_HITS)),
    evalue = 1e-30, slen = 2000L, time_stamp = 9L)
  DBI::dbAppendTable(con, "blast_hits", h)
  # another contig's hit, outside the window of its own copy
  DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession, pident, qcovs,
      aln_len, qstart, qend, sstart, send) VALUES ('S1', 1, 2, 1, 1, 'X1', 99, 33, 100, 1, 100, 1, 100)")
}

rw_orig_hits <- function(con, scaffold) {
  DBI::dbGetQuery(con, sprintf("SELECT %s FROM blast_hits WHERE ID = 'S1' AND path = 1 AND scaffold = ?
    ORDER BY hit_rank, hsp", paste(RiboPilot:::ddl_cols(RiboPilot:::BLAST_HITS_DDL), collapse = ", ")),
    params = list(scaffold))
}

test_that("review_map_hits() refuses a contig without a copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1, 300L)
  rw_hits(con)
  expect_error(RiboPilot:::review_map_hits(con, "S1", 1L, 1L), "No review copy for contig S1.1.1.", fixed = TRUE)
})

test_that("review_map_hits() on an untrimmed copy returns the original hits", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:2, 300L)
  rw_hits(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  rw_window(con, 2, 101L, 300L)
  orig <- rw_orig_hits(con, 1)
  expect_equal(RiboPilot:::review_map_hits(con, "S1", 1L, 1L), orig[!is.na(orig$qstart), ], ignore_attr = TRUE)
})

test_that("review_map_hits() clips and shifts hits to the kept window", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:3, 300L)
  rw_hits(con)
  rw_window(con, 2, 101L, 300L)
  rw_window(con, 1, 51L, 250L)
  rw_window(con, 3, 1L, 300L)
  orig <- rw_orig_hits(con, 1)
  m <- RiboPilot:::review_map_hits(con, "S1", 1L, 1L)

  expect_named(m, DBI::dbListFields(con, "review_blast_hits"))
  # subjects 2 and 7 fall outside; the others keep their stored rank and hsp
  expect_equal(m[, c("hit_rank", "hsp", "qstart", "qend", "sstart", "send", "aln_len", "qcovs")],
    data.frame(hit_rank = c(1L, 1L, 3L, 3L, 4L, 4L, 5L, 6L, 6L), hsp = c(1L, 2L, 1L, 3L, 1L, 2L, 1L, 1L, 2L),
      qstart = c(1L, 11L, 191L, 1L, 180L, 150L, 200L, 1L, 1L), qend = c(200L, 40L, 200L, 15L, 200L, 200L, 200L, 10L, 1L),
      sstart = c(NA, 500L, NA, 1000L, NA, 77L, NA, NA, NA), send = c(NA, 470L, NA, 1014L, NA, 127L, NA, NA, NA),
      aln_len = c(200L, 31L, 10L, 15L, 21L, 51L, 1L, 10L, 1L), qcovs = c(100, 100, 13, 13, 26, 26, 1, 5, 5)))
  # everything else as stored
  kept <- match(paste(m$hit_rank, m$hsp), paste(orig$hit_rank, orig$hsp))
  rest <- setdiff(names(m), c("qstart", "qend", "sstart", "send", "aln_len", "qcovs"))
  expect_equal(m[, rest], orig[kept, rest], ignore_attr = TRUE)

  # appendable to the mirror as the recheck writes it
  DBI::dbExecute(con, "DELETE FROM review_blast_hits WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  DBI::dbAppendTable(con, "review_blast_hits", m)
  expect_equal(DBI::dbGetQuery(con, sprintf("SELECT %s FROM review_blast_hits WHERE scaffold = 1 ORDER BY hit_rank, hsp",
    paste(names(m), collapse = ", "))), m)

  # nothing overlaps the window, or the original has no hits
  for (s in 2:3) {
    z <- RiboPilot:::review_map_hits(con, "S1", 1L, s)
    expect_equal(nrow(z), 0L)
    expect_named(z, names(m))
  }
  expect_equal(rw_orig_hits(con, 1), orig)
})
```

What the window test pins (window 51..250 of a 300 bp contig, copy length 200, shift 50): HSPs ending at 50 (subject 2) or starting at 251 (subject 3, hsp 2) fall outside, while HSPs ending exactly at 51 or starting exactly at 250 overlap by 1 bp; subject 3's union of 25 bp gives 12.5, which `floor(x + 0.5)` rounds to 13 (R's `round()` gives 12); subject 4's two HSPs overlap, so their union is 51 bp (26%), not 72 bp; subjects 2 and 7 leave gaps in `hit_rank`. The other contigs' copies have other windows and S1.1.2 has its own hit, so a lookup that drops the unit key fails.

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-window.R")'`
Expected: FAIL, `[ FAIL 6 | WARN 0 | SKIP 0 | PASS 0 ]`. Every test errors with `object 'review_coverage' not found`, `object 'review_coverage_summary' not found`, or `object 'review_map_hits' not found` (the "refuses a contig without a copy" test reports that its error message does not match).

- [ ] **Step 3: Append the window readers to R/review.R**

Append to the end of `R/review.R` (after `review_drop()`):

```r

#' A contig's read coverage on the review layer (S4.5)
#'
#' The original coverage tracks, cut to the copy's kept window when the contig
#' has a copy. A short track is cut, never padded. A track is NULL when the
#' coverage row or the track is missing, or nothing of it lies in the window.
#' Coverage is read with SELECT * because older projects and test fixtures lack
#' the clips column.
#' @return list(depths, clips): numeric vectors in copy coordinates, or NULL.
#' @noRd
review_coverage <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  cv <- DBI::dbGetQuery(con, "SELECT * FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  w <- DBI::dbGetQuery(con,
    "SELECT trim_start, trim_end FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  cut <- function(col) {
    if (!nrow(cv) || !col %in% names(cv)) return(NULL)
    v <- parse_depths(cv[[col]][1])
    if (nrow(w)) v <- v[seq_along(v) >= w$trim_start & seq_along(v) <= w$trim_end]
    if (length(v)) v
  }
  list(depths = cut("depths"), clips = cut("clips"))
}

#' Coverage summary of a depth track, rounded as write_coverage_summary() rounds it
#' @return list(mean_depth, median_depth, min_depth, max_depth), or NULL for an
#'   empty track.
#' @noRd
review_coverage_summary <- function(dep) {
  if (!length(dep)) return(NULL)
  list(mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2),
    min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)))
}

#' A copy's BLAST hits mapped from the original's (S4.5)
#'
#' The original HSPs that overlap the copy's kept window, clipped to it and
#' shifted to copy coordinates. A clipped HSP gets aln_len = its clipped span
#' and NA sstart and send; pident and bitscore stay as stored. qcovs is
#' recomputed per subject as the union of its kept spans over the copy length,
#' rounded as BLAST rounds it. hit_rank and hsp keep their stored values, so a
#' subject or HSP clipped away leaves a gap.
#' @return data.frame with the blast_hits columns; zero rows when nothing overlaps.
#' @noRd
review_map_hits <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  w <- DBI::dbGetQuery(con,
    "SELECT trim_start, trim_end FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  if (!nrow(w)) stop(sprintf("No review copy for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
  h <- DBI::dbGetQuery(con, sprintf(
    "SELECT %s FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY hit_rank, hsp",
    paste(ddl_cols(BLAST_HITS_DDL), collapse = ", ")), params = p)
  h <- h[!is.na(h$qstart) & !is.na(h$qend) & h$qend >= w$trim_start & h$qstart <= w$trim_end, , drop = FALSE]
  qs <- pmax(h$qstart, w$trim_start)
  qe <- pmin(h$qend, w$trim_end)
  cut <- qs != h$qstart | qe != h$qend
  h$aln_len[cut] <- qe[cut] - qs[cut] + 1L
  h$sstart[cut] <- NA
  h$send[cut] <- NA
  h$qstart <- qs - w$trim_start + 1L
  h$qend <- qe - w$trim_start + 1L
  span <- stats::ave(seq_len(nrow(h)), h$hit_rank,
    FUN = function(i) length(unique(unlist(Map(seq, h$qstart[i], h$qend[i])))))
  h$qcovs <- floor(100 * span / (w$trim_end - w$trim_start + 1) + 0.5)
  rownames(h) <- NULL
  h
}
```

Notes for the implementer:
- `BLAST_HITS_DDL` (existing) and `ddl_cols()` (Task 1) live in `R/db_schema.R`, which loads before `R/review.R`. Reading `blast_hits` by explicit column list keeps the frame in `review_blast_hits` column order on every database (new, upgraded, or a fixture widened by `fixture_review_schema()`); `coverage` alone is read with `SELECT *`, because it is not one of the review tables and hand-built fixtures lack `clips`.
- BLAST reports query coordinates ascending (`qstart <= qend`; a minus-strand hit has `sstart > send`), so no strand handling is needed.
- Do not renumber `hit_rank` or `hsp`: `derive_blast_products()` takes the first row per `hit_rank` after `ORDER BY hit_rank, hsp` (R/blast_genbank.R:166-176), and `screen_one()` orders the same way (R/screen.R:335-337), so gaps are harmless and keep subjects traceable to the original rows.

Check the file parses: `Rscript -e 'invisible(parse("R/review.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-window.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 25 ]`.

Then the other tests of `R/review.R` and the source scans:
Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-copy", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` (PASS 159, 1, 237).

Mutation check (optional, then revert): changing `floor(100 * span / (w$trim_end - w$trim_start + 1) + 0.5)` to `round(100 * span / (w$trim_end - w$trim_start + 1))`, or `h$qend >= w$trim_start` to `h$qend > w$trim_start`, or `seq_along(v) >= w$trim_start` to `seq_along(v) > w$trim_start`, each makes test-review-window.R fail.

- [ ] **Step 5: Commit**

```bash
git add R/review.R tests/testthat/test-review-window.R
git commit -m "Review copies: coverage window and mapped BLAST hits"
```


### Task 4: Review edit functions

**Files:**
- Modify: `R/review.R` (append after the last function in the file; Tasks 2 and 3 wrote everything above it)
- Create: `tests/testthat/test-review-edits.R`
- Test: `tests/testthat/test-review-edits.R`

**Interfaces:**
- Consumes (Task 2, `R/review.R`): `with_savepoint(con, code, name = "rp_review")` (BEGIN IMMEDIATE outside a transaction, SAVEPOINT inside one; a block's result is its last value, never a `return()`), `review_has_copy(con, id, path, scaffold)`, `review_copy(con, id, path, scaffold)` (no-op when a copy exists; stops with `"No sequence for contig <ID>.<path>.<scaffold>."` when the unit has no `review_units_v` row), `review_invalidate(con, id, path, scaffold)` (check_switch 1, rev + 1, check results and curation partials cleared, counts recomputed; returns invisible new rev text), `REVIEW_EFF_PARTIALS`. (Task 1): the review tables and views, `REVIEW_TABLES`, test helpers `chimera_test_con()` and `mk_app_project()` (widened by `fixture_review_schema()`).
- Produces (`R/review.R`, all `@noRd`; tests and later tasks call them as `RiboPilot:::fn` or from package code):
  - `review_trim(con, id, path, scaffold, from, to)`: one `with_savepoint()`. `review_copy()`; `from`/`to` are copy coordinates, clamped to `1..nchar(copy sequence)`; stops with `"No sequence for contig <ID>.<path>.<scaffold>."` when the copy (or the unit) has no sequence and with `"Invalid trim range."` when `from > to` (or either is NA) after clamping; new window in original coordinates `trim_start + from - 1 .. trim_start + to - 1`; `sequence = substr(assemblies.sequence, window)`, `length = nchar`, `edited = 1`; `review_annotations` shifted by `from - 1`, features wholly outside deleted, cut features clamped with `partial_<end>_user = 1` in the feature's own orientation (left cut: start on `+`, stop on `-`; right cut the reverse), `edited = 1` on every kept feature; then `review_invalidate()`. Returns `invisible(<new length>)`. A failure anywhere leaves no copy and no change.
  - `review_set_bounds(con, id, path, scaffold, element, start, stop)`: one `with_savepoint()`: `review_copy()`; the feature's `start`, `stop`, `length = stop - start + 1`, `edited = 1`; `review_units.edited = 1`; `review_invalidate()`. Returns `invisible(<new rev text>)`. Stops with `"No feature <element> on contig <ID>.<path>.<scaffold>."` when the copy has no such feature (nothing written).
  - `review_set_partial(con, id, path, scaffold, element, end = c("start", "stop"), value)`: one `with_savepoint()`: `review_copy()`; `partial_<end>_user = value`, feature `edited = 1`; `review_units.edited = 1`; check state, rev, and results unchanged. Same "No feature" error. Returns `invisible(NULL)`.
  - `review_reset_session(con, id, path, scaffold, backup)`: `backup = data.frame(element, start, stop, partial_start_user, partial_stop_user)` as read from `review_annotations_v` when Edit was ticked. One `with_savepoint()`: no copy -> nothing written, `invisible(FALSE)`; else writes start, stop, length, and both user partials back for every backup feature the copy still has, then `review_invalidate()` only when a start or stop differed from the copy's current value. Returns `invisible(TRUE)` when it invalidated, else `invisible(FALSE)`.
  - `review_set_state(con, units, state)`: `units = data.frame(ID, path, scaffold)`, `state` 1 (Ready to run) or 0 (On hold), else stops with `"state must be 1 (Ready to run) or 0 (On hold)."`. Character vector, one per unit: `NA` where applied, `"unlock first"` when `annotate.annotate_lock = 1`, `"press Update on the Assemble step first"` for state 1 on a unit without a copy whose `annotate_switch` is 1 or 4 or whose `n_features` IS NULL. One `with_savepoint()` per unit (reads and write together): state 1 = `review_copy()` + `review_invalidate()`; state 0 = `review_copy()` + `check_switch = 0`, `time_stamp` now, `rev` unchanged. Errors (for example `review_copy()`'s) propagate; units before the failing one stay written.
  - Internal helper `review_edit_feature(con, id, path, scaffold, element, set, vals)` (used by the two functions above; no other caller).

This task adds the functions behind every Annotation Review edit that writes a contig's review copy (SPEC S5.1 first touch, S5.2 table): a trim cuts the copy's sequence from the original by its kept window and gives cut gene ends a user partial mark in the feature's own orientation (S4.4), a boundary nudge moves one feature, a 5'/3' choice sets a user partial without touching the check state, Reset and Discard put an Edit session's starting values back, and State queues (Ready to run) or holds (On hold) the recheck while refusing locked contigs and contigs the main pipeline has not finished (S5.2, D7). Each function is one `with_savepoint()` block, so the first-touch copy, the edit, and the invalidation commit together or not at all, and they also work inside an open transaction. The originals are never written (S4.1). Nothing in the app calls these yet: Task 15 wires Trim, Task 16 wires Edit mode and State. This task covers its part of Review Focus 2: both strands, a trim exactly at a gene end, and a second trim of an already trimmed copy.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-review-edits.R`:

```r
# Annotation Review edits on a contig's review copy: trim, boundary nudge,
# 5'/3' choice, Reset, and State (SPEC S5.1, S5.2).

RE_SEQ <- withr::with_seed(4, replicate(2, paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")))

# S1.1.1 has its features on the plus strand, S1.1.2 on the minus strand.
RE_FEATS <- list(
  list(element = c("18S", "ITS1", "5.8S", "ITS2", "28S"), start = c(11L, 41L, 51L, 61L, 71L),
    stop = c(40L, 50L, 60L, 70L, 95L), strand = "+"),
  list(element = c("28S", "ITS2", "5.8S", "ITS1", "18S"), start = c(6L, 31L, 41L, 51L, 61L),
    stop = c(30L, 40L, 50L, 60L, 90L), strand = "-"))

# Two done 100 bp contigs; curation marked every feature's 3' end partial.
re_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:2) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', 100, ?, 0)", params = list(s, RE_SEQ[s]))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features,
        missing, n_flags, curate_flags, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, 5, '', 1, 'truncated', '')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths) VALUES ('S1', 1, ?, ?)",
      params = list(s, paste(1:100, collapse = ",")))
    f <- RE_FEATS[[s]]
    for (i in 1:5) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length, strand,
          partial_stop, status, edited)
        VALUES ('S1', 1, ?, ?, ?, ?, ?, ?, 1, 'ok', 0)",
        params = list(s, f$element[i], f$start[i], f$stop[i], f$stop[i] - f$start[i] + 1L, f$strand))
    }
  }
}

re_ints <- function(d) {
  k <- setdiff(names(d), c("element", "sequence", "missing"))
  d[k] <- lapply(d[k], as.integer)
  d
}

re_feats <- function(con, s) {
  re_ints(DBI::dbGetQuery(con, "SELECT element, start, stop, length, edited, partial_start_user,
      partial_stop_user FROM review_annotations WHERE ID = 'S1' AND path = 1 AND scaffold = ? ORDER BY start",
    params = list(as.integer(s))))
}

re_unit <- function(con, s) {
  re_ints(DBI::dbGetQuery(con, "SELECT sequence, length, trim_start, trim_end, edited, check_switch,
      n_features, missing, n_flags FROM review_units WHERE ID = 'S1' AND path = 1 AND scaffold = ?",
    params = list(as.integer(s))))
}

re_rev <- function(con, s, minus = 0L) {
  DBI::dbGetQuery(con, "SELECT CAST(rev - ? AS TEXT) AS rev FROM review_units
    WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(as.integer(minus), as.integer(s)))$rev
}

re_rows <- function(con, s) {
  sum(vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = ?", t),
    params = list(as.integer(s)))$n, integer(1)))
}

re_originals <- function(con) {
  lapply(stats::setNames(nm = c("assemblies", "annotate", "annotations", "coverage")), function(t) {
    cols <- paste(DBI::dbListFields(con, t), collapse = ", ")
    DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", cols, t, cols))
  })
}

re_feat_df <- function(element, start, stop, edited, ps, pe) {
  data.frame(element = element, start = start, stop = stop, length = stop - start + 1L,
    edited = edited, partial_start_user = ps, partial_stop_user = pe)
}

test_that("a trim cuts the copy from the original window and marks cut ends on the plus strand", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 80), 60L)
  expect_equal(re_unit(con, 1), data.frame(sequence = substr(RE_SEQ[1], 21, 80), length = 60L,
    trim_start = 21L, trim_end = 80L, edited = 1L, check_switch = 1L, n_features = 5L, missing = "",
    n_flags = NA_))
  # 18S is cut on the left (its 5' end); 28S on the right (its 3' end)
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 60L), rep(1L, 5), c(1L, NA_, NA_, NA_, NA_), c(NA_, NA_, NA_, NA_, 1L)))
  # curation's partial ends are cleared until the recheck; the cut ends stay
  eff <- DBI::dbGetQuery(con, paste("SELECT", RiboPilot:::REVIEW_EFF_PARTIALS, "FROM review_annotations_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1 ORDER BY start"))
  expect_equal(as.integer(eff$partial_start), c(1L, NA_, NA_, NA_, NA_))
  expect_equal(as.integer(eff$partial_stop), c(NA_, NA_, NA_, NA_, 1L))
  rev1 <- re_rev(con, 1)

  # a second trim works in the copy's coordinates: 5..50 of the copy is 25..70
  # of the original; ITS2 now ends exactly at the right end, and 28S is gone
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 1L, 5, 50), 46L)
  expect_equal(re_unit(con, 1), data.frame(sequence = substr(RE_SEQ[1], 25, 70), length = 46L,
    trim_start = 25L, trim_end = 70L, edited = 1L, check_switch = 1L, n_features = 4L, missing = "28S",
    n_flags = NA_))
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element[1:4], c(1L, 17L, 27L, 37L),
    c(16L, 26L, 36L, 46L), rep(1L, 4), c(1L, NA_, NA_, NA_), rep(NA_, 4)))
  expect_equal(re_rev(con, 1, 1L), rev1)
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, orig_length, length, trim_start, trim_end
      FROM review_units_v WHERE ID = 'S1' AND path = 1 AND scaffold = 1"),
    data.frame(has_copy = 1L, orig_length = 100L, length = 46L, trim_start = 25L, trim_end = 70L))

  # the other contig and every original are untouched
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("on the minus strand a left cut marks the 3' end and a right cut the 5' end", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  # 28S is cut on the left; 18S ends exactly at the new right end, so it is not partial
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 2L, 11, 90), 80L)
  expect_equal(re_feats(con, 2), re_feat_df(RE_FEATS[[2]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 80L), rep(1L, 5), rep(NA_, 5), c(1L, NA_, NA_, NA_, NA_)))

  # a 5'/3' choice the user made stays on an uncut feature; a cut end overrides it
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 0
    WHERE scaffold = 2 AND element IN ('ITS1', '18S')")
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 2L, 1, 70), 70L)
  expect_equal(re_feats(con, 2), re_feat_df(RE_FEATS[[2]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 70L), rep(1L, 5), c(NA_, NA_, NA_, 0L, 1L), c(1L, NA_, NA_, NA_, NA_)))
  expect_equal(re_unit(con, 2)[, c("sequence", "length", "trim_start", "trim_end")],
    data.frame(sequence = substr(RE_SEQ[2], 11, 80), length = 70L, trim_start = 11L, trim_end = 80L))
  expect_equal(re_rows(con, 1), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("a trim clamps to the copy's sequence and refuses bad ranges and missing sequences", {
  dir <- mk_app_project()
  app <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(app))
  # stored length 6120, sequence 'ACGT': the sequence decides
  expect_equal(RiboPilot:::review_trim(app, "S1", 1L, 1L, 2, 6000), 3L)
  expect_equal(RiboPilot:::review_trim(app, "S1", 1L, 1L, -5, 2), 2L)
  expect_equal(DBI::dbGetQuery(app, "SELECT sequence, length, trim_start, trim_end FROM review_units"),
    data.frame(sequence = "CG", length = 2L, trim_start = 2L, trim_end = 3L))

  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con), add = TRUE)
  re_seed(con)
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, ignore) VALUES ('S1', 1, 3, 100, 0)")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, n_features)
    VALUES ('S1', 1, 3, 2, 0)")
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 9L, 1, 10), "No sequence for contig S1.1.9.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 3L, 1, 10), "No sequence for contig S1.1.3.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 1L, 50, 40), "Invalid trim range.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 1L, 150, 200), "Invalid trim range.", fixed = TRUE)
  expect_equal(vapply(c(1, 3, 9), function(s) re_rows(con, s), integer(1)), c(0L, 0L, 0L))
})

test_that("a boundary nudge writes the copy and marks its checks not run", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_
  f0 <- re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop, rep(0L, 5), rep(NA_, 5),
    rep(NA_, 5))

  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40)
  f <- f0; f[1, c("start", "length", "edited")] <- c(5L, 36L, 1L)
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_unit(con, 1), data.frame(sequence = RE_SEQ[1], length = 100L, trim_start = 1L,
    trim_end = 100L, edited = 1L, check_switch = 1L, n_features = 5L, missing = "", n_flags = NA_))
  rev1 <- re_rev(con, 1)

  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "28S", 71, 99)
  f[5, c("stop", "length", "edited")] <- c(99L, 29L, 1L)
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1, 1L), rev1)

  expect_error(RiboPilot:::review_set_bounds(con, "S1", 1L, 2L, "XX", 1, 5), "No feature XX on contig S1.1.2.",
    fixed = TRUE)
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("a 5'/3' choice writes the user column and keeps the checks", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "stop", 0)
  rev1 <- re_rev(con, 1)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "start", 1)
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop,
    c(1L, 0L, 0L, 0L, 1L), c(NA_, NA_, NA_, NA_, 1L), c(0L, NA_, NA_, NA_, NA_)))
  expect_equal(re_unit(con, 1)[, c("edited", "check_switch", "n_flags")],
    data.frame(edited = 1L, check_switch = 2L, n_flags = 1L))
  expect_equal(re_rev(con, 1), rev1)
  # the user's choice wins over curation's, which stays until the next recheck
  eff <- DBI::dbGetQuery(con, paste("SELECT", RiboPilot:::REVIEW_EFF_PARTIALS, "FROM review_annotations_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1 ORDER BY start"))
  expect_equal(as.integer(eff$partial_start), c(NA_, NA_, NA_, NA_, 1L))
  expect_equal(as.integer(eff$partial_stop), c(0L, 1L, 1L, 1L, 1L))

  expect_error(RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "middle", 1), "should be one of")
  expect_error(RiboPilot:::review_set_partial(con, "S1", 1L, 2L, "XX", "start", 1),
    "No feature XX on contig S1.1.2.", fixed = TRUE)
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("Reset writes the session's start back and invalidates only when a boundary moved", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  NA_ <- NA_integer_
  backup <- function() DBI::dbGetQuery(con, "SELECT element, start, stop, partial_start_user, partial_stop_user
    FROM review_annotations_v WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  f0 <- re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop, rep(0L, 5), rep(NA_, 5),
    rep(NA_, 5))

  # Edit ticked on a contig without a copy, nothing changed: nothing to write
  b0 <- backup()
  expect_false(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b0))
  expect_equal(re_rows(con, 1), 0L)

  # a nudge and a toggle, then Reset
  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 1)
  rev1 <- re_rev(con, 1)
  expect_true(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b0))
  f <- f0; f$edited[c(1, 5)] <- 1L
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1, 1L), rev1)
  expect_equal(re_unit(con, 1)$check_switch, 1L)

  # a session with only toggles; the backup may name a feature the copy no longer has
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, n_flags = 4 WHERE scaffold = 1")
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 1)
  b1 <- rbind(backup(), data.frame(element = "XX", start = 1L, stop = 2L, partial_start_user = NA,
    partial_stop_user = NA))
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "start", 0)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 0)
  rev2 <- re_rev(con, 1)
  expect_false(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b1))
  f$partial_stop_user[5] <- 1L
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1), rev2)
  expect_equal(re_unit(con, 1)[, c("check_switch", "n_flags")], data.frame(check_switch = 2L, n_flags = 4L))
})

test_that("State queues or holds copies and refuses locked and waiting contigs", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  # 1 done; 2 locked; 3 running and 6 queued in the main pipeline; 4 never
  # annotated; 5 queued in the main pipeline but already has a copy
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold = 2")
  for (s in 3:6) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES ('S1', 1, ?, 4, 'ACGT', 0)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features)
      VALUES ('S1', 1, ?, ?, 0, ?)", params = list(s, c(4L, 2L, 1L, 1L)[s - 2], c(0L, NA, 0L, 0L)[s - 2]))
  }
  RiboPilot:::review_copy(con, "S1", 1L, 5L)
  rev5 <- re_rev(con, 5)
  u <- data.frame(ID = "S1", path = 1L, scaffold = 1:6)
  wait <- "press Update on the Assemble step first"
  sw <- function() DBI::dbGetQuery(con, "SELECT scaffold, check_switch FROM review_units ORDER BY scaffold")

  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::review_set_state(con, u, 1), c(NA, "unlock first", wait, wait, NA, wait))
  expect_equal(sw(), data.frame(scaffold = c(1L, 5L), check_switch = c(1L, 1L)))
  expect_equal(re_rev(con, 5, 1L), rev5)
  rev1 <- re_rev(con, 1); rev5 <- re_rev(con, 5)

  expect_equal(RiboPilot:::review_set_state(con, u, 0), c(NA, "unlock first", NA, NA, NA, NA))
  expect_equal(sw(), data.frame(scaffold = c(1L, 3:6), check_switch = rep(0L, 5)))
  expect_equal(c(re_rev(con, 1), re_rev(con, 5)), c(rev1, rev5))
  expect_true(all(DBI::dbGetQuery(con, "SELECT time_stamp FROM review_units")$time_stamp >= t0))
  expect_equal(re_rows(con, 2), 0L)

  # a locked contig's copy is left alone too
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold = 1")
  expect_equal(RiboPilot:::review_set_state(con, u[1, ], 1), "unlock first")
  expect_equal(re_unit(con, 1)$check_switch, 0L)
  expect_error(RiboPilot:::review_set_state(con, u, 2), "1 (Ready to run) or 0 (On hold)", fixed = TRUE)
})

test_that("each edit is one atomic write that also joins an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  one <- data.frame(ID = "S1", path = 1L, scaffold = 1L)
  edits <- list(
    function() RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 80),
    function() RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40),
    function() RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "start", 1),
    function() RiboPilot:::review_set_state(con, one, 1),
    function() RiboPilot:::review_set_state(con, one, 0))
  DBI::dbExecute(con, "CREATE TRIGGER re_boom BEFORE UPDATE ON review_units
    BEGIN SELECT RAISE(ABORT, 're_boom'); END")
  for (e in edits) {
    expect_error(e(), "re_boom")
    expect_equal(re_rows(con, 1), 0L)
  }
  DBI::dbExecute(con, "DROP TRIGGER re_boom")
  expect_false(RSQLite::sqliteIsTransacting(con))

  expect_error(DBI::dbWithTransaction(con, {
    for (e in edits) e()
    stop("outer failed")
  }), "outer failed")
  expect_equal(re_rows(con, 1), 0L)
  DBI::dbWithTransaction(con, {
    for (e in edits) e()
    RiboPilot:::review_reset_session(con, "S1", 1L, 1L, data.frame(element = "18S", start = 11L, stop = 40L,
      partial_start_user = NA, partial_stop_user = NA))
  })
  expect_equal(re_unit(con, 1)[, c("length", "check_switch")], data.frame(length = 60L, check_switch = 1L))
  expect_equal(re_feats(con, 1)$start[1:2], c(11L, 21L))
  expect_identical(re_originals(con), orig)
})
```

What the numbers pin. Plus strand, `review_trim(21, 80)` on the untouched 100 bp copy: the window is 21..80 (60 bp), every feature shifts by 20, 18S (11..40) loses its left end and becomes 1..20 with `partial_start_user = 1`, 28S (71..95) loses its right end and becomes 51..60 with `partial_stop_user = 1`; curation's `partial_stop = 1` on every feature is cleared by the invalidation, so the effective 3' marks are only the cut one. The second trim, 5..50 of the copy, is 25..70 of the original (`trim_start` 21 + 5 - 1 and 21 + 50 - 1); ITS2 (41..50 -> 37..46) ends exactly at the new length 46 and gets no mark, 28S (51..60 -> 47..56) falls outside and is deleted. Minus strand: `review_trim(11, 90)` cuts 28S (6..30) on the left, which is its 3' end (`partial_stop_user`), and 18S (61..90 -> 51..80) ends exactly at 80, so it is not marked; the next trim 1..70 cuts 18S on the right, its 5' end (`partial_start_user = 1`, replacing the 0 set before), while ITS1 keeps its user 0. The app fixture stores length 6120 for the sequence `ACGT`, so clamping with the stored length would give `trim_end` 6000 instead of 4. A swapped strand rule, a window computed from 1 instead of `trim_start`, `>=` for the right cut, a missing invalidation, a missing `edited`, a State that ignores the copy or the lock, an On hold that raises `rev`, a Reset that always invalidates or skips the user partials, a toggle that invalidates, and a trim outside one savepoint each fail at least one expectation here.

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-edits.R")'`
Expected: FAIL, `[ FAIL 8 | WARN 0 | SKIP 0 | PASS 0 ]`. Each of the 8 tests errors on the first call of a new function: `object 'review_trim' not found` (the three trim tests and the atomicity test), `object 'review_set_bounds' not found`, `object 'review_set_partial' not found`, `object 'review_reset_session' not found`, `object 'review_set_state' not found`.

- [ ] **Step 3: Add the edit functions to R/review.R**

Append the following to the end of `R/review.R`, after the last function already in the file (leave one blank line before it):

```r
#' Trim a contig's review copy to positions from..to of the copy (S5.2)
#'
#' The sequence is cut from the original by the new window. Features are
#' shifted, clamped, or dropped as trim_assembly did; a cut end gets
#' partial_*_user = 1 in the feature's own orientation. Makes the copy on first
#' touch and marks its checks not run.
#' @return invisible new length.
#' @noRd
review_trim <- function(con, id, path, scaffold, from, to) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  key <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  invisible(with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    r <- DBI::dbGetQuery(con,
      "SELECT r.sequence, r.trim_start, a.sequence AS orig FROM review_units r
       JOIN assemblies a ON a.ID = r.ID AND a.path = r.path AND a.scaffold = r.scaffold
       WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?", params = p)
    if (nrow(r) == 0 || is.na(r$sequence[1]) || !nzchar(r$sequence[1])) {
      stop(sprintf("No sequence for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
    }
    from <- max(1L, as.integer(from)); to <- min(nchar(r$sequence[1]), as.integer(to))
    if (!isTRUE(from <= to)) stop("Invalid trim range.", call. = FALSE)
    w <- as.integer(r$trim_start[1]) + c(from, to) - 1L
    sq <- substr(r$orig[1], w[1], w[2]); n <- nchar(sq)
    DBI::dbExecute(con, paste("UPDATE review_units SET sequence = ?, length = ?, trim_start = ?, trim_end = ?,
      edited = 1", key), params = c(list(sq, n, w[1], w[2]), p))
    ann <- DBI::dbGetQuery(con, paste("SELECT element, start, stop, strand FROM review_annotations", key),
      params = p)
    fkey <- paste(key, "AND element = ?")
    for (i in seq_len(nrow(ann))) {
      fp <- c(p, list(ann$element[i]))
      ns <- ann$start[i] - from + 1L; nt <- ann$stop[i] - from + 1L
      if (nt < 1 || ns > n) {
        DBI::dbExecute(con, paste("DELETE FROM review_annotations", fkey), params = fp)
        next
      }
      minus <- identical(ann$strand[i], "-")
      cut <- c(if (ns < 1) (if (minus) "partial_stop_user" else "partial_start_user"),
        if (nt > n) (if (minus) "partial_start_user" else "partial_stop_user"))
      for (col in cut) {
        DBI::dbExecute(con, sprintf("UPDATE review_annotations SET %s = 1 %s", col, fkey), params = fp)
      }
      ns <- max(1L, ns); nt <- min(n, nt)
      DBI::dbExecute(con, paste("UPDATE review_annotations SET start = ?, stop = ?, length = ?, edited = 1", fkey),
        params = c(list(ns, nt, nt - ns + 1L), fp))
    }
    review_invalidate(con, id, path, scaffold)
    n
  }))
}

#' Write one feature of a copy and mark feature and copy edited
#' @noRd
review_edit_feature <- function(con, id, path, scaffold, element, set, vals) {
  n <- DBI::dbExecute(con, sprintf("UPDATE review_annotations SET %s, edited = 1
    WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?", set),
    params = c(vals, list(id, path, scaffold, element)))
  if (n == 0) stop(sprintf("No feature %s on contig %s.%d.%d.", element, id, path, scaffold), call. = FALSE)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
}

#' Move a feature's boundaries on the copy (Edit mode nudge); checks go to not run
#' @return invisible new rev, as text.
#' @noRd
review_set_bounds <- function(con, id, path, scaffold, element, start, stop) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  start <- as.integer(start); stop <- as.integer(stop)
  invisible(with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    review_edit_feature(con, id, path, scaffold, element, "start = ?, stop = ?, length = ?",
      list(start, stop, stop - start + 1L))
    review_invalidate(con, id, path, scaffold)
  }))
}

#' Set the user's 5' (start) or 3' (stop) partial choice on a copy's feature
#'
#' No check reads partial ends, so the check state is kept.
#' @noRd
review_set_partial <- function(con, id, path, scaffold, element, end = c("start", "stop"), value) {
  end <- match.arg(end)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    review_edit_feature(con, id, path, scaffold, element, sprintf("partial_%s_user = ?", end),
      list(as.integer(value)))
  })
  invisible()
}

#' Write an Edit session's starting bounds and user partials back (Reset, Discard)
#'
#' backup: data.frame(element, start, stop, partial_start_user,
#' partial_stop_user) read from review_annotations_v when Edit was ticked.
#' Features the copy no longer has are skipped.
#' @return invisible TRUE when a boundary moved back and the checks were
#'   marked not run again; FALSE otherwise, including without a copy.
#' @noRd
review_reset_session <- function(con, id, path, scaffold, backup) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  invisible(with_savepoint(con, if (!review_has_copy(con, id, path, scaffold)) FALSE else {
    cur <- DBI::dbGetQuery(con,
      "SELECT element, start, stop FROM review_annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = p)
    b <- backup[backup$element %in% cur$element, , drop = FALSE]
    m <- match(b$element, cur$element)
    moved <- !identical(as.integer(c(cur$start[m], cur$stop[m])), as.integer(c(b$start, b$stop)))
    for (i in seq_len(nrow(b))) {
      DBI::dbExecute(con,
        "UPDATE review_annotations SET start = ?, stop = ?, length = ?, partial_start_user = ?,
           partial_stop_user = ?
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = c(list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, b$partial_start_user[i],
          b$partial_stop_user[i]), p, list(b$element[i])))
    }
    if (moved) review_invalidate(con, id, path, scaffold)
    moved
  }))
}

#' State on Annotation Review: 1 Ready to run (queue the recheck), 0 On hold
#'
#' Ready to run makes the copy if needed and marks its checks not run; On hold
#' sets check_switch 0 and keeps rev. Locked contigs are refused, and so is
#' Ready to run on a contig without a copy whose original the main pipeline
#' has queued or never annotated.
#' @param units data.frame(ID, path, scaffold).
#' @return character, NA where the state was set, else the reason.
#' @noRd
review_set_state <- function(con, units, state) {
  state <- as.integer(state)
  if (!isTRUE(state %in% 0:1)) stop("state must be 1 (Ready to run) or 0 (On hold).", call. = FALSE)
  vapply(seq_len(nrow(units)), function(i) {
    id <- units$ID[i]; path <- as.integer(units$path[i]); scaffold <- as.integer(units$scaffold[i])
    p <- list(id, path, scaffold)
    with_savepoint(con, {
      a <- DBI::dbGetQuery(con, "SELECT annotate_lock, annotate_switch, n_features FROM annotate
        WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
      if (isTRUE(a$annotate_lock[1] == 1)) {
        "unlock first"
      } else if (state == 1L && !review_has_copy(con, id, path, scaffold) &&
                 (a$annotate_switch[1] %in% c(1, 4) || is.na(a$n_features[1]))) {
        "press Update on the Assemble step first"
      } else {
        review_copy(con, id, path, scaffold)
        if (state == 1L) {
          review_invalidate(con, id, path, scaffold)
        } else {
          DBI::dbExecute(con, "UPDATE review_units SET check_switch = 0, time_stamp = ?
            WHERE ID = ? AND path = ? AND scaffold = ?", params = c(list(as.integer(Sys.time())), p))
        }
        NA_character_
      }
    })
  }, character(1))
}
```

Notes for the implementer:
- Never put `return()` inside a `with_savepoint()` block (it rolls the block back); each block above ends in its value. `review_reset_session()`'s no-copy case is the `if (...) FALSE else { ... }` expression for that reason.
- `review_trim()` clamps with `nchar()` of the copy's sequence, not `review_units.length` (hand-built fixtures store lengths that disagree with the sequence), and reads the original's sequence through the JOIN to cut the window. Its missing-unit error comes from `review_copy()` and matches `"No sequence"`, which `tests/testthat/test-batch-trim.R` ("one failing unit does not stop the rest") expects once Task 15 routes batch trims here.
- No new exported name, so NAMESPACE, man/, and `_pkgdown.yml` do not change. The user-facing strings say "contig" and "Assemble step".

Check the file parses: `Rscript -e 'invisible(parse("R/review.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-review-edits.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 80 ]`.

Then the other files that load `R/review.R` or scan the sources:
Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-copy", "review-window", "review-schema", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` (test-review-copy.R PASS 159, test-review-schema.R PASS 127, test-no-duplicate-definitions.R PASS 1, test-docs-site.R PASS 237; test-review-window.R as Task 3 left it).

- [ ] **Step 5: Commit**

```bash
git add R/review.R tests/testthat/test-review-edits.R
git commit -m "Review edit functions: trim, bounds, partials, reset, State"
```


### Task 5: ITSx states and stored origin

**Files:**
- Modify: `R/itsx_check.R` (end of `itsx_origin_scores`, origin block of `itsx_check_flags` :306-314, `itsx_check_units` :325-379, `itsx_store_warnings` :392-405, `itsx_check_one` :407-429; new `itsx_origin_flag`, `itsx_store_not_run`, `itsx_run_finished`)
- Modify: `R/pipeline_steps.R` (`reconcile_one` :102-130, `reset_unit_evidence` :243-260)
- Modify: `inst/nextflow/bin/rp_reconcile.R` (ITSx handler :14-17)
- Modify: `man/itsx_check_one.Rd` (regenerated, not edited by hand)
- Modify: `tests/testthat/helper-fake-tools.R` (fake ITSx :28-31)
- Test: `tests/testthat/test-itsx-check.R`, `tests/testthat/test-reconcile.R`, `tests/testthat/test-load-assemblies.R`

**Interfaces:**
- Consumes: `annotate.itsx_origin TEXT`, `annotate.itsx_scores TEXT`, and the `review_units` table, all created by `new_db()` since Task 1. Nothing else new.
- Produces:
  - `itsx_run_finished(out_prefix)`: TRUE when `<out_prefix>.summary.txt` contains "ITSx run finished" and `<out_prefix>.err` (absent counts as empty) has no line matching `ERROR` or `FATAL`.
  - `itsx_origin_flag(origin, scores, lineage)`: the `itsx_origin(<score comparison or group name>)` flag when the one-letter origin code is outside the groups the lineage expects, else `""`; `itsx_check_flags` calls it.
  - `itsx_check_units(con, units, out)`: data.frame(seqid, itsx_check, itsx_origin, itsx_scores, itsx_ITS1, itsx_5.8S, itsx_ITS2, itsx_notes). `itsx_origin` is the one-letter code and `itsx_scores` the verbatim group-score string, both `""` when ITSx reported none for a checked record and NA only on not_run rows.
  - `itsx_store_warnings(con, units, chk, review = FALSE)`: writes itsx_warnings, itsx_origin, itsx_scores to `annotate` (errors swallowed) or, with `review = TRUE`, to `review_units` (errors raised).
  - `itsx_store_not_run(con, id, path, scaffold, reason, review = FALSE)`: stores `not_run(<reason>)` with itsx_origin and itsx_scores NULL, through `itsx_store_warnings`.
  - `itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix)`: stores `not_run(ITSx did not finish)` (origin and scores NULL) unless `itsx_run_finished(out_prefix)`.
  - `reconcile_one()` clears itsx_warnings, itsx_origin, itsx_scores inside the annotations transaction; `reset_unit_evidence()` also clears chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores.
  - `rp_reconcile.R` stores `not_run(ITSx check failed: <message>)` when `itsx_check_one` throws.
  - The fake ITSx of `local_fake_pipeline_tools()` writes `<out>.summary.txt` containing "ITSx run finished".

A crashed stock ITSx run must read as `not_run`, never as "ITSx found nothing" (S6.3). ITSx exits 0 on fatal errors, so `itsx_check_one` now trusts a run only when its summary file says it finished and its captured stderr holds no `ERROR`/`FATAL` line; an exception inside the check is stored by `rp_reconcile.R` as `not_run(ITSx check failed: <message>)`, so NULL only ever means "never checked" (S6.3, S10). The check also stores ITSx's one-letter origin code and its group-score string in `annotate.itsx_origin`/`itsx_scores` (S4.1), and the origin flag moves into `itsx_origin_flag()` so Task 12 can recompute it after a taxonomy change. `reconcile_one` clears the three ITSx columns in the same transaction as the annotations rewrite, and `reset_unit_evidence` clears the six result columns of S5.3. `itsx_store_warnings` and `itsx_store_not_run` take `review = FALSE` now, because Task 10 stores a failed check on a copy with `itsx_store_not_run(review = TRUE)`; Task 8 only passes `review` through from `itsx_check_one`.

Tests that start `Rscript` use `skip_on_cran()`, as test-pipeline-drivers.R does, so the run commands below set `NOT_CRAN=true`; without it they are skipped, not run.

- [ ] **Step 1: Write the failing ITSx tests**

In `tests/testthat/test-itsx-check.R`, four edits to existing code, then one block appended at the end.

`mk_itsx_project` gains the two new columns (the store now writes three columns in one UPDATE, which fails whole, silently, on a table without them; the `itsx_check_one` and export-refresh tests read the stored value).

Find:

```r
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1)")
  dir
}
```

Replace with:

```r
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1)")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_origin TEXT")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_scores TEXT")
  dir
}
```

The env-gated real-ITSx test now gets the one-letter origin code (Task 18 ports this test later).

Find:

```r
  expect_equal(chk$itsx_origin, "Metazoa")
```

Replace with:

```r
  expect_equal(chk$itsx_origin, "M")
```

The `itsx_check_one` fixture needs the summary line a finished run writes, and a missing output now reads "ITSx did not finish".

Find:

```r
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, pre)
  expect_equal(out$itsx_check, "itsx_no_5.8S")
```

Replace with:

```r
  writeLines("ITSx run finished at Thu Oct  8 12:00:00 2026.", paste0(pre, ".summary.txt"))
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, pre)
  expect_equal(out$itsx_check, "itsx_no_5.8S")
```

Find:

```r
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, file.path(dir, "missing"))
  expect_match(out$itsx_check, "^not_run\\(")
```

Replace with:

```r
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, file.path(dir, "missing"))
  expect_equal(out$itsx_check, "not_run(ITSx did not finish)")
```

Append to the end of `tests/testthat/test-itsx-check.R`:

```r

# stock-ITSx output files under `pre`; finished = FALSE leaves the summary
# without its last line, as a killed run does
itsx_out <- function(pre, positions = character(0), results = character(0), finished = TRUE,
                     err = NULL) {
  writeLines(positions, paste0(pre, ".positions.txt"))
  writeLines(results, paste0(pre, ".extraction.results"))
  writeLines(c("Number of sequences in input file: 1",
    if (finished) "ITSx run finished at Thu Oct  8 12:00:00 2026."), paste0(pre, ".summary.txt"))
  if (!is.null(err)) writeLines(err, paste0(pre, ".err"))
  pre
}

test_that("an ITSx run counts as finished only with its summary line and no error on stderr", {
  pre <- file.path(withr::local_tempdir(), "unit.itsx_check")
  fin <- RiboPilot:::itsx_run_finished
  expect_false(fin(pre))
  itsx_out(pre, finished = FALSE)
  expect_false(fin(pre))
  itsx_out(pre)
  expect_true(fin(pre))
  writeLines("Thank you for using ITSx!", paste0(pre, ".err"))
  expect_true(fin(pre))
  for (e in c("FATAL ERROR :: Could not locate HMMER binaries!",
              "ERROR :: No input file was specified, and the --stdin option was set to F (false).",
              "FATAL ERROR! HMMER run seems to have been terminated before it was completed.")) {
    writeLines(c("Thank you for using ITSx!", e), paste0(pre, ".err"))
    expect_false(fin(pre), info = e)
  }
})

test_that("itsx_check_one stores the origin code and scores, and not_run when ITSx did not finish", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET lineage = 'Metazoa;Chordata' WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, itsx_warnings, itsx_origin, itsx_scores)
    VALUES ('S1', 1, 2, 'itsx_chimeric', 'F', 'F: 1 1e-10 40')")
  stored <- function() DBI::dbGetQuery(con, "SELECT itsx_warnings, itsx_origin, itsx_scores FROM annotate
    WHERE ID = 'S1' ORDER BY path, scaffold")
  d <- withr::local_tempdir()
  pos <- sub("^S1\t", "S1.1.1\t", pos_clean)
  res <- sub("^S1\t", "S1.1.1\t", sub("\tM\t", "\tF\t", res_clean, fixed = TRUE))
  sc <- "F: 4 1e-08 33.15, M: 4 2.9e-16 49.9"

  chk <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, itsx_out(file.path(d, "a"), pos, res))
  expect_equal(names(chk), c("seqid", "itsx_check", "itsx_origin", "itsx_scores", "itsx_ITS1",
    "itsx_5.8S", "itsx_ITS2", "itsx_notes"))
  expect_equal(unlist(stored()[1, ]), c(itsx_warnings = "itsx_origin(score: Fungi 33.15 vs Metazoa 49.9)",
    itsx_origin = "F", itsx_scores = sc))
  # checked, but ITSx gave no origin
  RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, itsx_out(file.path(d, "b"), pos))
  expect_equal(unlist(stored()[1, ]), c(itsx_warnings = "", itsx_origin = "", itsx_scores = ""))
  # a finished run that found nothing is a result, not a failure
  RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, itsx_out(file.path(d, "c")))
  expect_equal(stored()$itsx_warnings[1], "itsx_not_detected")
  # killed, failed with an error, or never ran: not_run, with no origin or scores
  runs <- list(itsx_out(file.path(d, "e"), pos, res, finished = FALSE),
    itsx_out(file.path(d, "f"), pos, res, err = "FATAL ERROR :: Could not locate HMMER binaries!"),
    file.path(d, "missing"))
  for (p in runs) {
    chk <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, p)
    expect_true(is.na(chk$itsx_origin) && is.na(chk$itsx_scores))
    got <- stored()
    expect_equal(got$itsx_warnings[1], "not_run(ITSx did not finish)", info = p)
    expect_true(is.na(got$itsx_origin[1]) && is.na(got$itsx_scores[1]), info = p)
    DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M', itsx_scores = 'M: 1 1e-10 40'
      WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  }
  expect_equal(unlist(stored()[2, ]), c(itsx_warnings = "itsx_chimeric", itsx_origin = "F",
    itsx_scores = "F: 1 1e-10 40"))
})

test_that("itsx_origin_flag builds the origin flag from a stored code and score string", {
  f <- RiboPilot:::itsx_origin_flag
  sc <- "F: 4 1e-08 33.15, M: 4 2.9e-16 49.9"
  expect_equal(f("F", sc, "Metazoa;Chordata"), "itsx_origin(score: Fungi 33.15 vs Metazoa 49.9)")
  expect_equal(f("F", "", "Metazoa;Chordata"), "itsx_origin(Fungi)")
  expect_equal(f("Z", NA_character_, "Metazoa;Chordata"), "itsx_origin(Z)")
  expect_equal(f("M", sc, "Metazoa;Chordata"), "")
  expect_equal(f("F", sc, "Bacteria;Proteobacteria"), "")
  expect_equal(f("F", sc, NA_character_), "")
  expect_equal(f("", sc, "Metazoa"), "")
  expect_equal(f(NA_character_, sc, "Metazoa"), "")
})

test_that("itsx_store_not_run writes the original or the review copy, never both", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = '', itsx_origin = 'M', itsx_scores = 'M: 1 1e-10 40'")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, itsx_warnings) VALUES ('S1', 1, 2, 'itsx_chimeric')")
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, rev, itsx_warnings, itsx_origin, itsx_scores)
    VALUES ('S1', 1, 1, 1, '', 'M', 'M: 1 1e-10 40')")
  get <- function(tbl) DBI::dbGetQuery(con, sprintf(
    "SELECT itsx_warnings, itsx_origin, itsx_scores FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = 1", tbl))

  RiboPilot:::itsx_store_not_run(con, "S1", 1L, 1L, "ITSx check failed: boom")
  an <- get("annotate")
  expect_equal(an$itsx_warnings, "not_run(ITSx check failed: boom)")
  expect_true(is.na(an$itsx_origin) && is.na(an$itsx_scores))
  expect_equal(get("review_units")$itsx_origin, "M")
  expect_equal(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate WHERE scaffold = 2")$itsx_warnings,
    "itsx_chimeric")

  RiboPilot:::itsx_store_not_run(con, "S1", 1L, 1L, "ITSx did not finish", review = TRUE)
  ru <- get("review_units")
  expect_equal(ru$itsx_warnings, "not_run(ITSx did not finish)")
  expect_true(is.na(ru$itsx_origin) && is.na(ru$itsx_scores))
  expect_equal(get("annotate")$itsx_warnings, "not_run(ITSx check failed: boom)")

  # a failed write raises in review mode and stays silent on the original
  for (tbl in c("annotate", "review_units")) {
    DBI::dbExecute(con, sprintf(
      "CREATE TRIGGER block_%s BEFORE UPDATE ON %s BEGIN SELECT RAISE(ABORT, 'blocked'); END", tbl, tbl))
  }
  expect_error(RiboPilot:::itsx_store_not_run(con, "S1", 1L, 1L, "x", review = TRUE), "blocked")
  expect_no_error(RiboPilot:::itsx_store_not_run(con, "S1", 1L, 1L, "x"))
})

test_that("the fake ITSx finishes like the real one, so pipeline tests store a real check", {
  skip_on_os("windows")
  tools <- local_fake_pipeline_tools()
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  unit <- file.path(d, "unit.fasta")
  writeLines(c(">S1.1.1", "ACGT"), unit)
  system2(file.path(tools$bin, "ITSx"), c("-i", unit, "-o", file.path(d, "unit.itsx_check")))
  expect_true(RiboPilot:::itsx_run_finished(file.path(d, "unit.itsx_check")))
  chk <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, file.path(d, "unit.itsx_check"))
  expect_equal(chk$itsx_check, "")
  expect_equal(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate WHERE ID = 'S1'")$itsx_warnings, "")
})
```

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'options(testthat.progress.max_fails = Inf); devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R")'` (the option keeps the reporter from quitting after 10 failures)
Expected: FAIL, `[ FAIL 16 | WARN 0 | SKIP 1 | PASS 74 ]`, from 6 tests; the skip is the real-ITSx test (no RIBOPILOT_TEST_PROJECT):
- "an ITSx run counts as finished only with its summary line and no error on stderr" and "the fake ITSx finishes like the real one, so pipeline tests store a real check": error `object 'itsx_run_finished' not found`.
- "itsx_origin_flag builds the origin flag from a stored code and score string": error `object 'itsx_origin_flag' not found`.
- "itsx_store_not_run writes the original or the review copy, never both": error `object 'itsx_store_not_run' not found`.
- "itsx_check_one stores the origin code and scores, and not_run when ITSx did not finish": failures (stored itsx_origin NA, expected "F"; a killed run stores the origin flag instead of `not_run(ITSx did not finish)`).
- "itsx_check_one stores stock-ITSx flags on the unit's annotate row": failure, `"not_run(no ITSx output)"` vs `"not_run(ITSx did not finish)"`.

- [ ] **Step 3: Implement the ITSx states, the stored origin and scores, and the fake summary**

In `R/itsx_check.R`, add `itsx_origin_flag()` after `itsx_origin_scores()`.

Find:

```r
    unname(ITSX_ORIGINS[best]) %|NA|% best, format(exp_sc[best]), tie)
}
```

Replace with:

```r
    unname(ITSX_ORIGINS[best]) %|NA|% best, format(exp_sc[best]), tie)
}

#' The itsx_origin flag for an ITSx origin code and the sample's lineage
#'
#' @param origin ITSx origin letter ("" or NA when ITSx gave none).
#' @param scores ITSx group-score string, verbatim.
#' @param lineage Sample lineage string or NA.
#' @return The itsx_origin flag when the origin is outside the groups the
#'   lineage expects, else "".
#' @noRd
itsx_origin_flag <- function(origin, scores, lineage) {
  expected <- itsx_expected_origins(lineage)
  if (is.null(expected) || is.na(origin) || !nzchar(origin) || origin %in% expected) return("")
  sc <- itsx_origin_scores(scores, origin, expected)
  sprintf("itsx_origin(%s)", if (nzchar(sc)) sc else unname(ITSX_ORIGINS[origin]) %|NA|% origin)
}
```

In `itsx_check_flags`, the origin block calls the helper.

Find:

```r
  expected <- itsx_expected_origins(lineage)
  if (nrow(my_res) > 0 && !is.null(expected)) {
    o <- my_res$origin[1]
    if (!is.na(o) && nzchar(o) && !(o %in% expected)) {
      sc <- itsx_origin_scores(my_res$scores[1], o, expected)
      flags <- c(flags, sprintf("itsx_origin(%s)", if (nzchar(sc)) sc else unname(ITSX_ORIGINS[o]) %|NA|% o))
    }
  }
```

Replace with:

```r
  if (nrow(my_res) > 0) {
    of <- itsx_origin_flag(my_res$origin[1], my_res$scores[1], lineage)
    if (nzchar(of)) flags <- c(flags, of)
  }
```

In `itsx_check_units`, three edits: the doc, the empty frame, and the per-record origin and scores.

Find:

```r
#' @return data.frame(seqid, itsx_check, itsx_origin, itsx_ITS1, itsx_5.8S,
#'   itsx_ITS2, itsx_notes). `itsx_check` is "" when clean, the flag list when
#'   not, or "not_run(<why>)" when ITSx could not be run.
```

Replace with:

```r
#' @return data.frame(seqid, itsx_check, itsx_origin, itsx_scores, itsx_ITS1,
#'   itsx_5.8S, itsx_ITS2, itsx_notes). `itsx_check` is "" when clean, the flag
#'   list when not, or "not_run(<why>)" when ITSx could not be run.
#'   `itsx_origin` is ITSx's one-letter origin code and `itsx_scores` its
#'   group-score string, verbatim; both are "" when ITSx reported none and NA
#'   only on not_run rows.
```

Find:

```r
    itsx_origin = rep(NA_character_, n), itsx_ITS1 = rep(NA_character_, n),
```

Replace with:

```r
    itsx_origin = rep(NA_character_, n), itsx_scores = rep(NA_character_, n),
    itsx_ITS1 = rep(NA_character_, n),
```

Find:

```r
    o <- res$origin[match(sid, res$seqid)]
    blank$itsx_origin[i] <- if (length(o) && !is.na(o)) unname(ITSX_ORIGINS[o]) %|NA|% o else NA_character_
```

Replace with:

```r
    j <- match(sid, res$seqid)
    blank$itsx_origin[i] <- res$origin[j] %|NA|% ""
    blank$itsx_scores[i] <- res$scores[j] %|NA|% ""
```

Replace `itsx_store_warnings` and add `itsx_store_not_run` and `itsx_run_finished` after it.

Find:

```r
#' Record each unit's flags in annotate.itsx_warnings
#'
#' Silent when the column is absent (a project made before it existed).
#' @noRd
itsx_store_warnings <- function(con, units, chk) {
  v <- chk$itsx_check[match(units$seqid, chk$seqid)]
  for (i in seq_len(nrow(units))) {
    tryCatch(DBI::dbExecute(con,
      "UPDATE annotate SET itsx_warnings = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(v[i], units$ID[i], units$path[i], units$scaffold[i])),
      error = function(e) NULL)
  }
  invisible(v)
}
```

Replace with:

```r
#' Record each unit's ITSx check: itsx_warnings, itsx_origin, itsx_scores
#'
#' Writes the annotate row, or the review copy (review_units) when `review` is
#' TRUE. On the annotate row a failed write is silent (a project made before
#' the columns existed); on a review copy it raises.
#' @noRd
itsx_store_warnings <- function(con, units, chk, review = FALSE) {
  m <- match(units$seqid, chk$seqid)
  sql <- sprintf("UPDATE %s SET itsx_warnings = ?, itsx_origin = ?, itsx_scores = ?
    WHERE ID = ? AND path = ? AND scaffold = ?", if (review) "review_units" else "annotate")
  for (i in seq_len(nrow(units))) {
    tryCatch(DBI::dbExecute(con, sql, params = list(chk$itsx_check[m[i]], chk$itsx_origin[m[i]],
      chk$itsx_scores[m[i]], units$ID[i], units$path[i], units$scaffold[i])),
      error = function(e) if (review) stop(e))
  }
  invisible(chk$itsx_check[m])
}

#' Store `not_run(<reason>)` as a unit's ITSx check, with no origin or scores
#' @noRd
itsx_store_not_run <- function(con, id, path, scaffold, reason, review = FALSE) {
  units <- data.frame(ID = id, path = as.integer(path), scaffold = as.integer(scaffold),
    seqid = paste(id, path, scaffold, sep = "."), stringsAsFactors = FALSE)
  itsx_store_warnings(con, units, itsx_check_units(con, units, list(ok = FALSE, message = reason)),
    review = review)
}

#' Did a stock ITSx run finish?
#'
#' ITSx exits 0 on fatal errors, so its exit status says nothing. A run counts
#' as finished when `<out_prefix>.summary.txt` contains "ITSx run finished" and
#' its captured stderr, `<out_prefix>.err` (absent counts as empty), has no
#' ERROR or FATAL line.
#' @noRd
itsx_run_finished <- function(out_prefix) {
  rd <- function(f) if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
  any(grepl("ITSx run finished", rd(paste0(out_prefix, ".summary.txt")), fixed = TRUE)) &&
    !any(grepl("ERROR|FATAL", rd(paste0(out_prefix, ".err"))))
}
```

In `itsx_check_one`, document the rule and use it.

Find:

```r
#' `<out_prefix>.*` files to reconcile, which calls this after the unit's
#' annotations are written.
```

Replace with:

```r
#' `<out_prefix>.*` files to reconcile, which calls this after the unit's
#' annotations are written. ITSx exits 0 on fatal errors, so a run counts only
#' when `<out_prefix>.summary.txt` says it finished and `<out_prefix>.err` has
#' no error line; otherwise the contig stores `not_run(ITSx did not finish)`.
```

Find:

```r
  out <- if (file.exists(paste0(out_prefix, ".positions.txt"))) {
    c(list(ok = TRUE), read_itsx_outputs(out_prefix))
  } else {
    list(ok = FALSE, message = "no ITSx output")
  }
```

Replace with:

```r
  out <- if (itsx_run_finished(out_prefix)) {
    c(list(ok = TRUE), read_itsx_outputs(out_prefix))
  } else {
    list(ok = FALSE, message = "ITSx did not finish")
  }
```

In `tests/testthat/helper-fake-tools.R`, the fake ITSx also writes the summary line a finished run ends with (without it, test-nf-main-run.R and test-pipeline-drivers.R would store `not_run(ITSx did not finish)`).

Find:

```r
        ' "5.8S: 2651-2807", "ITS2: 2808-3990", "LSU: 3991-6310", ""}\' "$in" > "$out.positions.txt"')),
```

Replace with:

```r
        ' "5.8S: 2651-2807", "ITS2: 2808-3990", "LSU: 3991-6310", ""}\' "$in" > "$out.positions.txt"'),
      'echo "ITSx run finished at $(date)." > "$out.summary.txt"'),
```

- [ ] **Step 4: Run the ITSx tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 113 ]` (29 tests; "a real ITSx run in the project's container over a genuine rDNA record is clean" skips without RIBOPILOT_TEST_PROJECT).

- [ ] **Step 5: Write the failing reconcile tests**

Append to the end of `tests/testthat/test-reconcile.R`:

```r

# barrnap GFF and tuned-ITSx positions for unit S1.1.1, as annotate_tools writes them
reconcile_inputs <- function(dir) {
  gff <- file.path(dir, "unit.barrnap.gff")
  writeLines(c(
    "S1.1.1\tbarrnap:0.9\trRNA\t351\t2150\t0\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA",
    "S1.1.1\tbarrnap:0.9\trRNA\t2651\t2807\t0\t+\t.\tName=5_8S_rRNA;product=5.8S ribosomal RNA",
    "S1.1.1\tbarrnap:0.9\trRNA\t3991\t6310\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA"), gff)
  pos <- file.path(dir, "unit.itsx.positions.txt")
  writeLines(paste0("S1.1.1\t6310 bp.\tSSU: 351-2150\tITS1: 2151-2650\t5.8S: 2651-2807\t",
    "ITS2: 2808-3990\tLSU: 3991-6310\tChimeric! "), pos)
  list(gff = gff, pos = pos)
}

test_that("reconcile_one writes chimeric and clears the stored ITSx check", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotate SET chimeric = 0, itsx_warnings = 'itsx_no_5.8S', itsx_origin = 'M',
    itsx_scores = 'M: 1 1e-10 40'")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, itsx_warnings, itsx_origin, itsx_scores)
    VALUES ('S1', 1, 2, '', 'M', 'M: 1 1e-10 40')")
  f <- reconcile_inputs(withr::local_tempdir())
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, gff = f$gff, positions = f$pos)
  an <- DBI::dbGetQuery(con, "SELECT annotate_switch, n_features, missing, chimeric, itsx_warnings,
    itsx_origin, itsx_scores FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")
  expect_equal(an$annotate_switch[1], 2L)
  expect_equal(an$n_features[1], 5L)
  expect_equal(an$missing[1], "")
  expect_equal(an$chimeric[1], 1L)
  expect_true(is.na(an$itsx_warnings[1]) && is.na(an$itsx_origin[1]) && is.na(an$itsx_scores[1]))
  expect_equal(an$itsx_origin[2], "M")
})

test_that("reconcile_one keeps the old features when the annotate update fails", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  f <- reconcile_inputs(withr::local_tempdir())
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, gff = f$gff, positions = f$pos)
  DBI::dbExecute(con,
    "CREATE TRIGGER block_annotate BEFORE UPDATE ON annotate BEGIN SELECT RAISE(ABORT, 'blocked'); END")
  expect_error(RiboPilot::reconcile_one(con, "S1", 1L, 1L), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations WHERE ID = 'S1'")$n, 5L)
})

test_that("rp_reconcile.R stores not_run when the ITSx check throws", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("pkgload")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!nzchar(src), "not a source checkout")
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  f <- reconcile_inputs(d)
  pre <- file.path(d, "unit.itsx_check")
  writeLines("ITSx run finished at Thu Oct  8 12:00:00 2026.", paste0(pre, ".summary.txt"))
  dir.create(paste0(pre, ".positions.txt"))   # unreadable: the check itself throws
  res <- processx::run("Rscript",
    c(file.path(system.file("nextflow", "bin", package = "RiboPilot"), "rp_reconcile.R"),
      DBI::dbGetInfo(con)$dbname, "S1", "1", "1", f$gff, f$pos, pre),
    env = c("current", RIBOPILOT_SRC = src), error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  an <- DBI::dbGetQuery(con, "SELECT n_features, itsx_warnings, itsx_origin FROM annotate
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  expect_equal(an$n_features, 5L)
  expect_match(an$itsx_warnings, "^not_run\\(ITSx check failed: .+\\)$")
  expect_true(is.na(an$itsx_origin))
})
```

- [ ] **Step 6: Run it to see it fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reconcile.R")'`
Expected: FAIL, `[ FAIL 3 | WARN 0 | SKIP 0 | PASS 49 ]`, one failure in each of 3 tests:
- "reconcile_one writes chimeric and clears the stored ITSx check": `is.na(an$itsx_warnings[1]) && is.na(an$itsx_origin[1]) && is.na(an$itsx_scores[1])` is FALSE.
- "reconcile_one keeps the old features when the annotate update fails": annotation count 0, expected 5 (the DELETE committed before the blocked UPDATE).
- "rp_reconcile.R stores not_run when the ITSx check throws": `an$itsx_warnings` is NA and does not match `^not_run\(ITSx check failed: .+\)$`. Without `NOT_CRAN=true` this test is skipped instead.

- [ ] **Step 7: Move the annotate update into the transaction and store the ITSx failure**

In `R/pipeline_steps.R` (`reconcile_one`), compute `missing` before the transaction.

Find:

```r
  ts <- as.integer(Sys.time())

  # one transaction: a failed append after a committed delete would leave the
```

Replace with:

```r
  ts <- as.integer(Sys.time())
  missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)

  # one transaction: a failed append after a committed delete would leave the
```

Then make the annotate UPDATE the last statement inside the transaction and clear the three ITSx columns.

Find:

```r
      DBI::dbAppendTable(con, "annotations", df)
    }
  })
  missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
  DBI::dbExecute(con,
    "UPDATE annotate SET n_features = ?, missing = ?, chimeric = ?, annotate_switch = 2, time_stamp = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(nrow(feat), paste(missing, collapse = ","), chimeric, ts, id, path, scaffold))
```

Replace with:

```r
      DBI::dbAppendTable(con, "annotations", df)
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET n_features = ?, missing = ?, chimeric = ?, annotate_switch = 2,
         itsx_warnings = NULL, itsx_origin = NULL, itsx_scores = NULL, time_stamp = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(nrow(feat), paste(missing, collapse = ","), chimeric, ts, id, path, scaffold))
  })
```

In `inst/nextflow/bin/rp_reconcile.R`, store the failure instead of only messaging. The store call is an unexported internal, so it is written `RiboPilot:::` (test-nf-bin-contract.R:119-130).

Find:

```r
    error = function(e) message("ITSx check skipped: ", conditionMessage(e)))
```

Replace with:

```r
    error = function(e) {
      message("ITSx check failed: ", conditionMessage(e))
      tryCatch(RiboPilot:::itsx_store_not_run(con, a[2], as.integer(a[3]), as.integer(a[4]),
        paste("ITSx check failed:", conditionMessage(e))),
        error = function(e2) message("ITSx check result not stored: ", conditionMessage(e2)))
    })
```

- [ ] **Step 8: Run the reconcile tests to see them pass**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-reconcile.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 52 ]` (15 tests).

- [ ] **Step 9: Write the failing re-assembly reset test**

In `tests/testthat/test-load-assemblies.R`, `mk_asm_db`'s annotate table gains the six columns `reset_unit_evidence` now clears (without them its UPDATE errors inside `load_assemblies`' transaction). No test in this file inserts into annotate by position.

Find:

```r
                       partial TEXT, problematic TEXT, export_group TEXT, n_features INT,
                       PRIMARY KEY (ID, path, scaffold))")
```

Replace with:

```r
                       partial TEXT, problematic TEXT, export_group TEXT, n_features INT,
                       chimeric INT, n_flags INT, curate_flags TEXT, itsx_warnings TEXT,
                       itsx_origin TEXT, itsx_scores TEXT, PRIMARY KEY (ID, path, scaffold))")
```

The reset test gives both units check results, then expects them cleared on the changed unit only.

Find:

```r
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, annotate_lock = 1, reviewed = 'yes',
    partial = 'yes', export_group = 'g1'")
```

Replace with:

```r
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, annotate_lock = 1, reviewed = 'yes',
    partial = 'yes', export_group = 'g1', chimeric = 1, n_flags = 2, curate_flags = 'short_18S',
    itsx_warnings = 'itsx_chimeric', itsx_origin = 'F', itsx_scores = 'F: 1 1e-10 40'")
```

Find:

```r
  expect_equal(an$export_group, c(NA, "g1"))
```

Replace with:

```r
  expect_equal(an$export_group, c(NA, "g1"))
  for (col in c("chimeric", "n_flags", "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores")) {
    expect_true(is.na(an[[col]][1]), info = col)
    expect_false(is.na(an[[col]][2]), info = col)
  }
```

- [ ] **Step 10: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-load-assemblies.R")'`
Expected: FAIL, `[ FAIL 6 | WARN 0 | SKIP 0 | PASS 125 ]`, all six in "a re-assembly that changes a unit's sequence clears its old evidence and lock": `is.na(an[[col]][1])` is FALSE for chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, and itsx_scores. The other 37 tests pass.

- [ ] **Step 11: Clear the check results on a re-assembled unit**

In `R/pipeline_steps.R` (`reset_unit_evidence`), two edits.

Find:

```r
#' Its gene coordinates, BLAST hits, screen verdict, chimera call, and depth then
#' describe another molecule, so they are deleted, and the unit's lock and
#' curator marks are cleared so it is annotated and reviewed again.
```

Replace with:

```r
#' Its gene coordinates, BLAST hits, screen verdict, chimera call, check
#' results, and depth then describe another molecule, so they are deleted, and
#' the unit's lock and curator marks are cleared so it is annotated and
#' reviewed again.
```

Find:

```r
       partial = 'no', problematic = 'no', export_group = NULL, n_features = NULL
     WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
```

Replace with:

```r
       partial = 'no', problematic = 'no', export_group = NULL, n_features = NULL,
       chimeric = NULL, n_flags = NULL, curate_flags = NULL, itsx_warnings = NULL,
       itsx_origin = NULL, itsx_scores = NULL
     WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
```

- [ ] **Step 12: Run the re-assembly tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-load-assemblies.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 131 ]` (38 tests).

- [ ] **Step 13: Regenerate the help page of itsx_check_one**

Its roxygen description changed in Step 3 and it is exported. Never run `devtools::document()` (NAMESPACE is hand-maintained).

Run: `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`
Then: `git status --short man/`
Expected: only ` M man/itsx_check_one.Rd`; its description now ends with the sentence about `\verb{<out_prefix>.summary.txt}` and `\verb{not_run(ITSx did not finish)}`. If any other man/ file changed, it is not this task's: restore it with `git checkout -- man/<that file>`.

- [ ] **Step 14: Run every test file this task touches or could break**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("itsx-check", "reconcile", "load-assemblies", "nf-bin-contract", "pipeline-drivers", "write-atomicity", "sample-state", "export-units", "export-headers", "nmnh", "export", "docs-site")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: PASS, every file reports `FAIL 0`. Reporter tallies at the time of writing, in this order: itsx-check `SKIP 1 | PASS 113` (the real-ITSx test skips without RIBOPILOT_TEST_PROJECT), reconcile `PASS 52`, load-assemblies `PASS 131`, nf-bin-contract `PASS 162`, pipeline-drivers `PASS 17`, write-atomicity `PASS 22`, sample-state `PASS 24`, export-units `PASS 53`, export-headers `PASS 60`, nmnh `PASS 116`, export `PASS 79`, docs-site `PASS 237`.

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-main-run.R")'`
Expected: PASS (4 tests). It runs the real main.nf with the fake tools, so the stock check in reconcile now depends on the fake summary file; it takes about 4 minutes.

- [ ] **Step 15: Commit**

```bash
git add R/itsx_check.R R/pipeline_steps.R inst/nextflow/bin/rp_reconcile.R man/itsx_check_one.Rd tests/testthat/helper-fake-tools.R tests/testthat/test-itsx-check.R tests/testthat/test-reconcile.R tests/testthat/test-load-assemblies.R
git commit -m "ITSx run states, stored origin and scores"
```


### Task 6: Review mode: curation, reference alignment, live alignment view

**Files:**
- Modify: `R/review.R` (append `review_tbl()` at the end of the file; Task 2 created it)
- Modify: `R/curate_rrna.R` (`curate_one()` and its roxygen block, lines 279-366)
- Modify: `R/blast_ref_align.R` (`compute_feature_msa()` roxygen and first two reads, lines 5-14; `compute_ref_alignment()` and its roxygen block, lines 71-134)
- Modify: `man/curate_one.Rd`, `man/compute_ref_alignment.Rd` (regenerated by roxygen, never edited by hand)
- Test: `tests/testthat/test-curate.R` (append)
- Test: `tests/testthat/test-ref-alignment.R` (append)

**Interfaces:**
- Consumes: Task 1: `REVIEW_MIRRORS` (named chr, original -> mirror for annotations, blast_hits, blast_screen, screen, chimera, chimera_windows, blast_ref_alignment); `new_db()` creates the 8 review tables and the 5 views, so `chimera_test_con()` has them. Task 2: `with_savepoint(con, code, name = "rp_review")`; `review_copy(con, id, path, scaffold)` (tests only).
- Produces:
  - `review_tbl(name, review = FALSE)` (R/review.R, `@noRd`): returns `name` when `review` is FALSE; in review mode `REVIEW_MIRRORS[[name]]` for the 7 mirrored tables, `"review_units"` for `"assemblies"` and `"annotate"`, otherwise `stop("No review table for <name>.")`. Callers route only the 7 mirrors, assemblies (sequence, length, blast_*), and annotate result columns (n_features, missing, chimeric, n_flags, curate_flags, itsx_*, time_stamp). samples, curate_opts, annotate review state, blast_searched, blast_ref_candidates, blast_ref_sequences, blast_ref_annotations, and coverage stay literal in both modes.
  - `curate_one(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE)` (exported). `ruleset` is the 5th positional argument, so callers always pass `review` by name.
  - `compute_ref_alignment(con, id, path = 1L, scaffold = 1L, review = FALSE)` (exported).
  - `compute_feature_msa(con, id, path = 1L, scaffold = 1L, element, max_hits = Inf)`: signature unchanged; reads `review_units_v` and `review_annotations_v`.

The recheck scores a review copy with the same functions the main pipeline uses (SPEC S7.3). This task adds `review_tbl()`, the one name map every review-mode check uses, and gives `curate_one()` and `compute_ref_alignment()` a `review = FALSE` argument that routes their per-contig reads and writes to the review tables, while the shared tables stay literal. In review mode curation replaces `partial_start`/`partial_stop` on every run instead of filling NULLs, because the user's choice now lives in `partial_start_user`/`partial_stop_user` (S4.4); original mode keeps COALESCE. `curate_one()`'s write block moves from `DBI::dbWithTransaction()` to `with_savepoint()`, so it also works inside an open transaction (the upgrade step, S11); its reads stay outside any transaction, as today. `compute_feature_msa()`, whose only caller is the Annotation Review details window (R/app_server.R:2684), reads the views, so the live alignment shows the copy's sequence and bounds (S8).

- [ ] **Step 1: Write the failing curation tests**

Append this at the end of `tests/testthat/test-curate.R` (one blank line before it). It reuses `good_feats()` and `FISH` from the top of that file and `chimera_test_con()` from helper-chimera.R.

```r
curate_seed <- function(con) {
  DBI::dbExecute(con, "UPDATE samples SET lineage = ?, lineage_taxids = '33208;7711;186623;;;;'",
    params = list(FISH))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 6000, ?)", params = list(strrep("A", 6000)))
  f <- good_feats()
  f$model_coverage[f$element == "28S"] <- 60
  DBI::dbAppendTable(con, "annotations", cbind(ID = "S1", path = 1L, scaffold = 1L, f))
  DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, lineage, lineage_taxids, is_nuclear_rrna)
    VALUES ('S1', 1, 1, 'Eukaryota; Sar; Alveolata; Dinophyceae; Symbiodinium', '2759;2698737;33630;2864;;;', 1)")
  DBI::dbExecute(con, "INSERT INTO blast_ref_alignment (ID, path, scaffold, accession, element, pid, aln_cov)
    VALUES ('S1', 1, 1, 'ACC1', '28S', 60, 0.9)")
}

test_that("review_tbl routes only the review-layer tables", {
  rt <- RiboPilot:::review_tbl
  expect_equal(rt("annotations"), "annotations")
  expect_equal(rt("annotate", review = FALSE), "annotate")
  expect_equal(unname(vapply(names(RiboPilot:::REVIEW_MIRRORS), rt, "", review = TRUE)),
    unname(RiboPilot:::REVIEW_MIRRORS))
  expect_equal(rt("blast_ref_alignment", TRUE), "review_ref_alignment")
  expect_equal(rt("assemblies", TRUE), "review_units")
  expect_equal(rt("annotate", TRUE), "review_units")
  for (t in c("samples", "curate_opts", "blast_searched", "blast_ref_candidates",
              "blast_ref_sequences", "blast_ref_annotations", "coverage")) {
    expect_error(rt(t, TRUE), t, fixed = TRUE)
  }
})

test_that("curate_one(review = TRUE) curates the copy and leaves the original alone", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  curate_seed(con)
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  orig <- function() list(
    DBI::dbGetQuery(con, "SELECT element, status, flags, partial_start, partial_stop
      FROM annotations ORDER BY element"),
    DBI::dbGetQuery(con, "SELECT n_flags, curate_flags, time_stamp FROM annotate"))
  before <- orig()
  expect_equal(before[[2]]$curate_flags, "contaminant")
  # the copy: 18S truncated, 28S complete, its own screen hit and reference identity
  DBI::dbExecute(con, "UPDATE review_annotations SET model_coverage = 70, partial_start_user = 0
    WHERE element = '18S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET model_coverage = 100 WHERE element = '28S'")
  DBI::dbExecute(con, "UPDATE review_blast_screen SET lineage = ?, lineage_taxids = NULL",
    params = list(FISH))
  DBI::dbExecute(con, "UPDATE review_ref_alignment SET element = '18S'")
  DBI::dbExecute(con, "UPDATE review_units SET time_stamp = 0")

  RiboPilot::curate_one(con, "S1", 1L, 1L, review = TRUE)
  expect_identical(orig(), before)
  ra <- DBI::dbGetQuery(con, "SELECT element, status, flags, partial_start, partial_stop,
    partial_start_user FROM review_annotations")
  rownames(ra) <- ra$element
  expect_equal(ra["18S", "flags"], "low_identity,truncated")
  expect_equal(ra["18S", "partial_start"], 1L)
  expect_true(is.na(ra["18S", "partial_stop"]))
  expect_equal(ra["18S", "partial_start_user"], 0L)
  expect_equal(ra["28S", "status"], "ok")
  expect_true(is.na(ra["28S", "partial_stop"]))  # replaced, not kept
  ru <- DBI::dbGetQuery(con, "SELECT n_flags, curate_flags, time_stamp FROM review_units")
  expect_equal(ru$n_flags, 0L)
  expect_true(is.na(ru$curate_flags))
  expect_gt(ru$time_stamp, 0)

  # original mode never touches the copy
  copy <- DBI::dbGetQuery(con, "SELECT element, status, flags, partial_start, partial_stop
    FROM review_annotations ORDER BY element")
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  expect_identical(DBI::dbGetQuery(con, "SELECT element, status, flags, partial_start, partial_stop
    FROM review_annotations ORDER BY element"), copy)
})

test_that("curate_one writes inside an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  curate_seed(con)
  DBI::dbBegin(con)
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT curate_flags FROM annotate")$curate_flags, "contaminant")
  DBI::dbRollback(con)
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT curate_flags FROM annotate")$curate_flags))
})
```

Why these values: the original has a Symbiodinium screen hit (contaminant), a low 28S identity, and a truncated 28S (partial_stop 1). The copy gets a fish screen hit, the low identity on 18S, a truncated 18S, and a complete 28S, so every read or write that hits the wrong table, and a COALESCE left in review mode (28S partial_stop would stay 1), fails an expectation. The sample and the original screen hit both carry taxids (the contaminant check compares taxids when both sides have them) and the copy's hit has none, so reading either the lineage or the taxids from the original screen row brings back the contaminant flag.

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-curate.R")'`
Expected: FAIL 3 (`[ FAIL 3 | WARN 0 | SKIP 0 | PASS 39 ]`): `object 'review_tbl' not found`; `unused argument (review = TRUE)`; `cannot start a transaction within a transaction` (from `DBI::dbWithTransaction` inside `curate_one`). The existing tests still pass.

- [ ] **Step 3: Add review_tbl()**

Append this at the end of `R/review.R` (one blank line before it):

```r
#' Table a check function reads or writes
#'
#' The name itself when review is FALSE. In review mode the 7 mirrored tables
#' map to their review_* mirror, and assemblies (sequence, length, blast_*) and
#' annotate (result columns only) map to review_units. Shared tables and
#' coverage are never routed, so asking for one is an error.
#' @noRd
review_tbl <- function(name, review = FALSE) {
  if (!review) return(name)
  if (name %in% names(REVIEW_MIRRORS)) return(REVIEW_MIRRORS[[name]])
  if (name %in% c("assemblies", "annotate")) return("review_units")
  stop("No review table for ", name, ".", call. = FALSE)
}
```

- [ ] **Step 4: Give curate_one() the review argument and a savepoint write**

In `R/curate_rrna.R`, Find:

```r
#' @param ruleset Optional ruleset override (else read from `curate_opts`).
#' @return (invisibly) the curate_rrna_core() result.
#' @export
curate_one <- function(con, id, path = 1L, scaffold = 1L, ruleset = NULL) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)

  feats <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop, length, strand, model_coverage FROM annotations
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
```

Replace with:

```r
#' @param ruleset Optional ruleset override (else read from `curate_opts`).
#' @param review TRUE curates the contig's review copy: it reads and writes the
#'   review tables and replaces the copy's partial ends on every run (the user's
#'   choice is kept in `partial_start_user` and `partial_stop_user`).
#' @return (invisibly) the curate_rrna_core() result.
#' @export
curate_one <- function(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  ann_tbl <- review_tbl("annotations", review)
  scr_tbl <- review_tbl("blast_screen", review)

  feats <- DBI::dbGetQuery(con, sprintf(
    "SELECT element, type, start, stop, length, strand, model_coverage FROM %s
     WHERE ID = ? AND path = ? AND scaffold = ?", ann_tbl),
    params = list(id, path, scaffold))
```

Find:

```r
  rp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT element, pid, aln_cov FROM blast_ref_alignment WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
```

Replace with:

```r
  rp <- tryCatch(DBI::dbGetQuery(con, sprintf(
    "SELECT element, pid, aln_cov FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("blast_ref_alignment", review)),
    params = list(id, path, scaffold)), error = function(e) data.frame())
```

Find:

```r
  scr <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage, is_nuclear_rrna FROM blast_screen WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  hit_lineage <- if (nrow(scr) > 0) scr$lineage[1] else NA_character_
  hit_is_nuc <- if (nrow(scr) > 0) as.integer(scr$is_nuclear_rrna[1]) else NA_integer_
  hit_taxids <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids t FROM blast_screen WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$t[1], error = function(e) NA_character_)
```

Replace with:

```r
  scr <- tryCatch(DBI::dbGetQuery(con, sprintf(
    "SELECT lineage, is_nuclear_rrna FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", scr_tbl),
    params = list(id, path, scaffold)), error = function(e) data.frame())
  hit_lineage <- if (nrow(scr) > 0) scr$lineage[1] else NA_character_
  hit_is_nuc <- if (nrow(scr) > 0) as.integer(scr$is_nuclear_rrna[1]) else NA_integer_
  hit_taxids <- tryCatch(DBI::dbGetQuery(con, sprintf(
    "SELECT lineage_taxids t FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", scr_tbl),
    params = list(id, path, scaffold))$t[1], error = function(e) NA_character_)
```

Find:

```r
  ts <- as.integer(Sys.time())
  DBI::dbWithTransaction(con, {
    if (nrow(res$features) > 0) {
      for (i in seq_len(nrow(res$features))) {
        DBI::dbExecute(con,
          "UPDATE annotations SET status = ?, flags = ?, partial_start = COALESCE(partial_start, ?),
             partial_stop = COALESCE(partial_stop, ?)
           WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
          params = list(res$features$status[i],
            if (nzchar(res$features$flags[i])) res$features$flags[i] else NA_character_,
            res$features$partial_start[i], res$features$partial_stop[i],
            id, path, scaffold, res$features$element[i]))
      }
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET n_flags = ?, curate_flags = ?, time_stamp = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(res$n_flags, res$notes, ts, id, path, scaffold))
  })
```

Replace with:

```r
  ts <- as.integer(Sys.time())
  partials <- if (review) "partial_start = ?, partial_stop = ?" else
    "partial_start = COALESCE(partial_start, ?), partial_stop = COALESCE(partial_stop, ?)"
  with_savepoint(con, {
    if (nrow(res$features) > 0) {
      for (i in seq_len(nrow(res$features))) {
        DBI::dbExecute(con, sprintf(
          "UPDATE %s SET status = ?, flags = ?, %s
           WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?", ann_tbl, partials),
          params = list(res$features$status[i],
            if (nzchar(res$features$flags[i])) res$features$flags[i] else NA_character_,
            res$features$partial_start[i], res$features$partial_stop[i],
            id, path, scaffold, res$features$element[i]))
      }
    }
    DBI::dbExecute(con, sprintf(
      "UPDATE %s SET n_flags = ?, curate_flags = ?, time_stamp = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      review_tbl("annotate", review)),
      params = list(res$n_flags, res$notes, ts, id, path, scaffold))
  })
```

The block stays write-first and holds only `DBI::dbExecute` calls, so nothing in it opens a transaction of its own.

- [ ] **Step 5: Run the curation tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-curate.R")'`
Expected: PASS (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 65 ]`). The original-mode test "curate_one fills partial ends once and keeps the curator's choice" still passes, because original mode keeps COALESCE.

- [ ] **Step 6: Write the failing alignment tests**

Append this at the end of `tests/testthat/test-ref-alignment.R` (one blank line before it). The existing tests above it keep using the hand-built `mk_aln_db()`, which has no review tables; that is fine because original mode reads the same tables as before.

```r
aln_seed <- function(con) {
  s <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 200, ?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length)
    VALUES ('S1', 1, 1, '18S', 1, 200, 200)")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates (ID, path, scaffold, rank, accession, species)
    VALUES ('S1', 1, 1, 1, 'ACC1', 'Danio rerio')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence) VALUES ('ACC1', ?)",
    params = list(s))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations (accession, element, start, stop, strand)
    VALUES ('ACC1', '18S', 1, 200, '+')")
  s
}

# The copy: one base changed every 10 bp over the first 100, and 18S cut to 1..100.
aln_edit_copy <- function(con, s) {
  b <- strsplit(s, "")[[1]]
  i <- seq(10, 100, by = 10)
  b[i] <- chartr("ACGT", "CATG", b[i])
  m <- paste(b, collapse = "")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?", params = list(m))
  DBI::dbExecute(con, "UPDATE review_annotations SET stop = 100, length = 100")
  m
}

test_that("compute_ref_alignment(review = TRUE) aligns the copy and writes only its mirror", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  s <- aln_seed(con)
  expect_equal(as.integer(RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)), 1L)
  aln_edit_copy(con, s)
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1, review = TRUE)
  expect_equal(as.integer(n), 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT element, pid, aln_len FROM review_ref_alignment"),
    data.frame(element = "18S", pid = 90, aln_len = 100L))
  expect_equal(DBI::dbGetQuery(con, "SELECT element, pid, aln_len FROM blast_ref_alignment"),
    data.frame(element = "18S", pid = 100, aln_len = 200L))
})

test_that("the live alignment reads the copy when there is one, else the original", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  s <- aln_seed(con)
  sample_row <- function() {
    m <- RiboPilot:::compute_feature_msa(con, "S1", 1L, 1L, "18S")
    gsub("-", "", as.character(m)[["sample|18S"]])
  }
  expect_equal(sample_row(), s)
  m <- aln_edit_copy(con, s)
  expect_equal(sample_row(), substr(m, 1, 100))
})
```

The copied `review_ref_alignment` row has the same key as the new one, so a review-mode DELETE that misses the mirror makes the append fail on the primary key, and one that hits `blast_ref_alignment` empties the original.

- [ ] **Step 7: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-ref-alignment.R")'`
Expected: FAIL 2 (`[ FAIL 2 | WARN 0 | SKIP 0 | PASS 11 ]`): `unused argument (review = TRUE)`; and `Expected sample_row() to equal substr(m, 1, 100)`, because `compute_feature_msa()` still reads the original 200 bp sequence.

- [ ] **Step 8: Route compute_ref_alignment() and point compute_feature_msa() at the views**

In `R/blast_ref_align.R`, Find:

```r
#' parsed feature for the element). Uses the sample's CURRENT boundaries, so the
#' view redraws on edits. Returns an aligned DNAStringSet (sample first) or NULL.
#' @noRd
compute_feature_msa <- function(con, id, path = 1L, scaffold = 1L, element, max_hits = Inf) {
  s_seq <- app_query_con(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    id, path, scaffold)$sequence
  sa <- app_query_con(con, "SELECT start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
    id, path, scaffold, element)
```

Replace with:

```r
#' parsed feature for the element). Reads the review layer (the contig's copy
#' when it has one), so the view redraws on edits. Returns an aligned
#' DNAStringSet (sample first) or NULL.
#' @noRd
compute_feature_msa <- function(con, id, path = 1L, scaffold = 1L, element, max_hits = Inf) {
  s_seq <- app_query_con(con, "SELECT sequence FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
    id, path, scaffold)$sequence
  sa <- app_query_con(con, "SELECT start, stop FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
    id, path, scaffold, element)
```

Find:

```r
#' @param path,scaffold Assembly unit within the sample.
#' @return (invisibly) number of element alignments written.
#' @export
compute_ref_alignment <- function(con, id, path = 1L, scaffold = 1L) {
  ts <- as.integer(Sys.time())
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  DBI::dbExecute(con, "DELETE FROM blast_ref_alignment WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
```

Replace with:

```r
#' @param path,scaffold Assembly unit within the sample.
#' @param review TRUE aligns the contig's review copy (its sequence and
#'   features) and writes `review_ref_alignment`; the reference tables are shared.
#' @return (invisibly) number of element alignments written.
#' @export
compute_ref_alignment <- function(con, id, path = 1L, scaffold = 1L, review = FALSE) {
  ts <- as.integer(Sys.time())
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  out_tbl <- review_tbl("blast_ref_alignment", review)
  DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", out_tbl),
    params = list(id, path, scaffold))
```

Find:

```r
  s_seq <- app_query_con(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    id, path, scaffold)$sequence
  r_seq <- app_query_con(con, "SELECT sequence FROM blast_ref_sequences WHERE accession = ?", acc)$sequence
```

Replace with:

```r
  s_seq <- app_query_con(con, sprintf("SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("assemblies", review)), id, path, scaffold)$sequence
  r_seq <- app_query_con(con, "SELECT sequence FROM blast_ref_sequences WHERE accession = ?", acc)$sequence
```

Find:

```r
  s_ann <- app_query_con(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    id, path, scaffold)
```

Replace with:

```r
  s_ann <- app_query_con(con, sprintf("SELECT element, start, stop FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("annotations", review)), id, path, scaffold)
```

Find:

```r
    DBI::dbAppendTable(con, "blast_ref_alignment", data.frame(ID = id, path = path, scaffold = scaffold,
```

Replace with:

```r
    DBI::dbAppendTable(con, out_tbl, data.frame(ID = id, path = path, scaffold = scaffold,
```

`DBI::dbAppendTable` names every column of the data frame in its INSERT, so it is not positional. The candidate, reference sequence, and reference annotation reads stay literal (shared tables), and the function keeps today's non-atomic shape (DELETE, then one append per element).

- [ ] **Step 9: Run the alignment tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-ref-alignment.R")'`
Expected: PASS (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 15 ]`).

- [ ] **Step 10: Regenerate the two help pages**

`curate_one` and `compute_ref_alignment` are exported and have man pages. NAMESPACE is hand-maintained: never run `devtools::document()`.

Run:

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
git status --short man/ NAMESPACE
```

Expected: exactly ` M man/compute_ref_alignment.Rd` and ` M man/curate_one.Rd`. Their usage lines now read `curate_one(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE)` and `compute_ref_alignment(con, id, path = 1L, scaffold = 1L, review = FALSE)`, each with a new `\item{review}{...}`. roxygen prints "Could not resolve link to topic" warnings for internal topics; they were there before this task. If any other man/ file is listed, put it back (`git checkout -- <file>` for a modified one, `rm <file>` for a new one) and mention it in your report; NAMESPACE must not be listed.

- [ ] **Step 11: Run the affected test files**

Run each:

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-curate.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-ref-alignment.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-no-duplicate-definitions.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-batch-trim.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-export.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-pipeline-drivers.R")'
```

Expected: PASS everywhere: test-curate.R 65, test-ref-alignment.R 15, and no failures in the rest. test-docs-site.R checks the regenerated Rd for backticks and `[fn()]` links. test-batch-trim.R and test-app-export.R run `curate_one()` through batch trim, the trim handler, and Save (original mode, now through `with_savepoint()`). `NOT_CRAN=true` makes test-pipeline-drivers.R run (not skip) its bin-script test, which calls `rp_curate.R` on a `new_db()` project; expect 0 skipped.

- [ ] **Step 12: Commit**

```bash
git add R/review.R R/curate_rrna.R R/blast_ref_align.R man/curate_one.Rd man/compute_ref_alignment.Rd tests/testthat/test-curate.R tests/testthat/test-ref-alignment.R
git commit -m "Review mode for curation and reference alignment"
```


### Task 7: Review mode: BLAST products and screen

**Files:**
- Modify: `R/screen.R` (`screen_one()` and its roxygen block, lines 308-376)
- Modify: `R/blast_genbank.R` (`parse_blast_candidates()` and its roxygen block, lines 40-135; `derive_blast_products()` and its roxygen block, lines 137-212; `write_blast_screen()` and its roxygen block, lines 266-292)
- Modify: `man/screen_one.Rd`, `man/parse_blast_candidates.Rd`, `man/derive_blast_products.Rd` (regenerated by roxygen, never edited by hand)
- Test: `tests/testthat/test-screen.R` (append)
- Test: `tests/testthat/test-blast-parse.R` (append)
- Test: `tests/testthat/test-write-atomicity.R` (append)

**Interfaces:**
- Consumes: Task 1: `new_db()` creates `review_units` and the 7 mirror tables (`review_blast_hits`, `review_blast_screen`, `review_screen` among them), so `chimera_test_con()` has them. Task 2: `review_copy(con, id, path, scaffold)` (tests only; it copies `blast_hits`, `blast_screen`, and `screen` rows into their mirrors and the `blast_*` summary into `review_units`). Task 6: `review_tbl(name, review = FALSE)` (R/review.R).
- Produces:
  - `screen_one(con, id, path = 1L, scaffold = 1L, force = FALSE, review = FALSE)` (exported). Review mode: the version skip reads `review_screen`; hits come from `review_blast_hits`; `samples`, `screen_params()`, and `blast_searched` stay literal; `variant_class`, `variant_group`, `variant_identity` are read from the original `screen` row before the transaction (NULL when it has none) and written into the new row; the write is a DELETE-first `DBI::dbWithTransaction()` on `review_screen`. Never calls `screen_variants()`. Returns the verdict invisibly.
  - `parse_blast_candidates(con, id, path = 1L, scaffold = 1L, hits_file, taxdump_dir = NULL, review = FALSE)` (exported). Review mode: hits (and the no-hits clear) go to `review_blast_hits`; `ensure_sample_lineage_taxids()` is never called; calls `derive_blast_products(..., review = TRUE)`.
  - `derive_blast_products(con, id, path = 1L, scaffold = 1L, ts = as.integer(Sys.time()), review = FALSE)` (exported). Review mode: reads `review_blast_hits`; never deletes or writes `blast_ref_candidates` (the candidate frame is still built in memory to pick the top hit and is returned); replaces the unit's `review_blast_screen` row; sets `review_units.blast_*` from the top candidate, or clears all six when there is no candidate or no `review_blast_hits` row; `blast_source` still comes from the original `blast_searched.source` (`local_rrna` for `local`, else `core_nt`). Note for Task 10: with zero `review_blast_hits` rows it clears the copied summary, so map mode skips this call when it wants to keep the summary (S7.3 step 2).
  - `write_blast_screen(con, id, path, scaffold, top, hits, ts = as.integer(Sys.time()), tbl = "blast_screen")` (`@noRd`): appends to `tbl`.
  - On a contig without a copy, review mode writes mirror rows that no view shows (no `review_units` row) and its `review_units` UPDATEs match nothing; no error. Callers (Tasks 10 and 12) only call it on copies.

The recheck scores a review copy with the same BLAST and screen code the main pipeline uses (SPEC S7.3 steps 2 and 3). This task gives `parse_blast_candidates()`, `derive_blast_products()`, and `screen_one()` a `review = FALSE` argument that routes their per-contig reads and writes through `review_tbl()` (Task 6): hits to `review_blast_hits`, the top screening hit to `review_blast_screen`, the `blast_*` summary to `review_units`, and the verdict to `review_screen`. Shared tables stay literal in both modes: `samples`, `blast_searched` (so the `local_rrna`/`core_nt` label and the "bundled rDNA database" hint still follow the original search), and `blast_ref_candidates`, which review mode never touches, because the recheck aligns against the original best reference (S4.2, S7.3 step 8). Review mode skips `ensure_sample_lineage_taxids()` (it writes `samples`), and review-mode `screen_one()` carries the original row's variant facts instead of running `screen_variants()`, because variant grouping compares every contig of a sample and belongs to the Assemble step (D1, S7.3 step 3); the " | variant: ..." clause `screen_variants()` appends to the original's `screen_reason` is not carried, since no review-layer reader shows `screen_reason`. Every `DBI::dbWithTransaction()` block keeps today's write-first shape, and original mode behaves as before.

- [ ] **Step 1: Write the failing tests**

Append this at the end of `tests/testthat/test-screen.R` (one blank line before it). It uses `FISH_T` and `FUNGI_T` from the top of that file and `chimera_test_con()` from helper-chimera.R.

```r
# The original holds three fish hits, an on_target verdict, and variant facts
# from the Assemble step; the copy holds three fungal hits and no screen row
# yet, as review_invalidate() leaves it.
review_screen_seed <- function(con) {
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = ?", params = list(FISH_T))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 8, 'ACGTACGT')")
  ins <- function(tbl, tax, lin) {
    for (i in 1:3) {
      DBI::dbExecute(con, sprintf(
        "INSERT INTO %s (ID, path, scaffold, hit_rank, hsp, bitscore, lineage_taxids, lineage, hit_class)
         VALUES ('S1', 1, 1, ?, 1, ?, ?, ?, 'rdna')", tbl),
        params = list(i, 900 - i, tax, lin))
    }
  }
  ins("blast_hits", FISH_T, "Eukaryota;Chordata;Actinopteri;o;f;g;s")
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "UPDATE screen SET variant_class = 'isomer', variant_group = 'S1.g1',
    variant_identity = 99.5, screen_reason = 'ORIG'")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "DELETE FROM review_screen")
  DBI::dbExecute(con, "DELETE FROM review_blast_hits")
  ins("review_blast_hits", FUNGI_T, "Eukaryota;Ascomycota;Sordariomycetes;o;f;g;s")
}

test_that("screen_one(review = TRUE) scores the copy's hits and carries the original's variant facts", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  review_screen_seed(con)
  orig <- function() DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason, screen_version,
    lca_n_hits, variant_class, variant_group, variant_identity, time_stamp FROM screen")
  before <- orig()
  expect_equal(before$screen_verdict, "on_target")
  copy_row <- function() DBI::dbGetQuery(con, "SELECT screen_verdict, screen_class, screen_reason,
    screen_version, variant_class, variant_group, variant_identity FROM review_screen")

  # the original row is current, so a version skip that read it would write nothing
  expect_equal(screen_one(con, "S1", 1, 1, review = TRUE), "off_target")
  r <- copy_row()
  expect_equal(nrow(r), 1L)
  expect_equal(r$screen_verdict, "off_target")
  expect_equal(r$screen_class, "fungal")
  expect_equal(r$screen_version, RiboPilot:::SCREEN_VERSION)
  expect_equal(r[, c("variant_class", "variant_group", "variant_identity")],
    data.frame(variant_class = "isomer", variant_group = "S1.g1", variant_identity = 99.5))
  expect_identical(orig(), before)

  # the version skip reads the copy's row
  DBI::dbExecute(con, "UPDATE review_screen SET screen_reason = 'SENTINEL'")
  screen_one(con, "S1", 1, 1, review = TRUE)
  expect_equal(copy_row()$screen_reason, "SENTINEL")
  screen_one(con, "S1", 1, 1, force = TRUE, review = TRUE)
  r <- copy_row()
  expect_equal(nrow(r), 1L)
  expect_false(identical(r$screen_reason, "SENTINEL"))
  expect_equal(r$variant_class, "isomer")

  # no hits on the copy: the search source still comes from blast_searched
  DBI::dbExecute(con, "DELETE FROM review_blast_hits")
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('S1', 1, 1, 'x', 'local')")
  expect_equal(screen_one(con, "S1", 1, 1, force = TRUE, review = TRUE), "no_reference")
  expect_match(copy_row()$screen_reason, "bundled")
  expect_identical(orig(), before)
})
```

Append this at the end of `tests/testthat/test-blast-parse.R` (one blank line before it). It uses `hit()` and `write_hits()` from the top of that file. The mock counts calls of `ensure_sample_lineage_taxids()`, which original mode makes because `chimera_test_con()` leaves `samples.lineage_taxids` NULL.

```r
test_that("parse_blast_candidates(review = TRUE) writes only the copy's tables", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  called <- 0L
  local_mocked_bindings(.package = "RiboPilot",
    ensure_sample_lineage_taxids = function(...) {
      called <<- called + 1L
      invisible(0L)
    })
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 8, 'ACGTACGT')")
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('S1', 1, 1, 'x', 'local')")
  parse_blast_candidates(con, "S1", 1, 1,
    write_hits(hit("ORIG1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600)))
  expect_equal(called, 1L)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  orig <- function() list(
    DBI::dbGetQuery(con, "SELECT hit_rank, hsp, accession, pident FROM blast_hits ORDER BY hit_rank, hsp"),
    DBI::dbGetQuery(con, "SELECT rank, accession, time_stamp FROM blast_ref_candidates ORDER BY rank"),
    DBI::dbGetQuery(con, "SELECT accession, is_nuclear_rrna, time_stamp FROM blast_screen"),
    DBI::dbGetQuery(con, "SELECT blast_accession, blast_species, blast_pident, blast_qcovs,
      blast_evalue, blast_source FROM assemblies"))
  before <- orig()
  expect_equal(before[[4]]$blast_accession, "ORIG1")
  summ <- function() DBI::dbGetQuery(con, "SELECT blast_accession, blast_species, blast_pident,
    blast_qcovs, blast_source FROM review_units")

  res <- parse_blast_candidates(con, "S1", 1, 1, write_hits(c(
    hit("CONTAM", "Saprolegnia parasitica strain X genome assembly", 99.0, 90, 5000, 6000, 1, 6000,
      staxid = 4769),
    hit("NEWREF", "Danio rerio 28S ribosomal RNA gene", 92.0, 70, 1200, 1400, 1, 1400))),
    review = TRUE)
  expect_equal(called, 1L)
  expect_identical(orig(), before)
  expect_equal(res$accession, "NEWREF")
  expect_equal(DBI::dbGetQuery(con, "SELECT accession FROM review_blast_hits ORDER BY hit_rank")$accession,
    c("CONTAM", "NEWREF"))
  expect_equal(DBI::dbGetQuery(con, "SELECT accession, is_nuclear_rrna FROM review_blast_screen"),
    data.frame(accession = "CONTAM", is_nuclear_rrna = 0L))
  expect_equal(summ(), data.frame(blast_accession = "NEWREF", blast_species = "Danio rerio",
    blast_pident = 92, blast_qcovs = 70, blast_source = "local_rrna"))

  # no candidate left: the summary is cleared, the screen hit stays
  DBI::dbExecute(con, "UPDATE review_blast_hits SET pident = 10, qcovs = 5, aln_len = 10")
  expect_equal(nrow(RiboPilot:::derive_blast_products(con, "S1", 1, 1, review = TRUE)), 0L)
  expect_true(all(vapply(summ(), is.na, logical(1))))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_blast_screen")$n, 1L)
  expect_identical(orig(), before)

  # no hits: the copy's hits, screen hit, and summary are cleared
  DBI::dbExecute(con, "UPDATE review_units SET blast_accession = 'X'")
  parse_blast_candidates(con, "S1", 1, 1, write_hits("NO_SIGNIFICANT_HITS"), review = TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_blast_hits")$n, 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_blast_screen")$n, 0L)
  expect_true(is.na(summ()$blast_accession))
  expect_identical(orig(), before)
})
```

Append this at the end of `tests/testthat/test-write-atomicity.R` (one blank line before it). It uses `block_inserts()` and `FISH_TAX` from the top of that file.

```r
test_that("a failed review_screen insert leaves the copy's stored verdict in place", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = ?", params = list(FISH_TAX))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 8, 'ACGTACGT')")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  for (i in 1:3) {
    DBI::dbExecute(con,
      "INSERT INTO review_blast_hits (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class)
       VALUES ('S1',1,1,?,1,?,?,'Eukaryota;Chordata;Actinopteri;o;f;g;s','rdna')",
      params = list(i, 900 - i, FISH_TAX))
  }
  screen_one(con, "S1", 1, 1, review = TRUE)
  before <- DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM review_screen")
  expect_equal(before$screen_verdict, "on_target")

  block_inserts(con, "review_screen")
  expect_error(screen_one(con, "S1", 1, 1, force = TRUE, review = TRUE), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM review_screen"), before)
})

test_that("a failed review_blast_hits insert leaves the copy's stored hits in place", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 8, 'ACGTACGT')")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  f <- withr::local_tempfile()
  writeLines(paste("AA1", "Danio rerio 18S ribosomal RNA gene", 99, 80, "1e-90", 7955,
    1500, 1600, 1, 1600, 1, 1600, sep = "\t"), f)
  parse_blast_candidates(con, "S1", 1, 1, f, review = TRUE)
  before <- DBI::dbGetQuery(con, "SELECT ID, accession FROM review_blast_hits")
  expect_equal(nrow(before), 1L)

  block_inserts(con, "review_blast_hits")
  expect_error(parse_blast_candidates(con, "S1", 1, 1, f, review = TRUE), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT ID, accession FROM review_blast_hits"), before)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_blast_screen")$n, 1L)
})
```

Why these values: the original and the copy disagree on everything a routed read or write can touch. The original's hits are fish (on_target) and the copy's fungal (off_target), so a review-mode hit read from `blast_hits` fails the verdict. The original `screen` row is current, so a version skip that reads it returns early and leaves `review_screen` empty. `review_screen` starts empty, so variant facts read from it instead of `screen` come back NULL. The BLAST copy gets a contaminant on top and a different reference, so every write that lands in `blast_hits`, `blast_ref_candidates`, `blast_screen`, or `assemblies` breaks `expect_identical(orig(), before)`, and every read from an original gives the wrong `review_blast_screen` or `review_units` values. The `blast_searched` row says `local`, so a hard-coded label or a routed `blast_searched` read fails `blast_source` and the "bundled" reason. A primary-key clash catches a DELETE aimed at the wrong table: the next append to the mirror fails.

- [ ] **Step 2: Run them to see them fail**

Run:

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-screen.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-blast-parse.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-write-atomicity.R")'
```

Expected: FAIL in each, every time with `unused argument (review = TRUE)`:
- test-screen.R `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 117 ]`
- test-blast-parse.R `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 70 ]` (`NOT_CRAN=true` makes its rp_searched.R test run instead of skip)
- test-write-atomicity.R `[ FAIL 2 | WARN 0 | SKIP 0 | PASS 22 ]`

- [ ] **Step 3: Give screen_one() the review argument**

In `R/screen.R`, Find:

```r
#' @param force Recompute even when the stored row is current.
#' @return (invisibly) the verdict written, or the stored one when skipped.
#' @export
screen_one <- function(con, id, path = 1L, scaffold = 1L, force = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  cur <- DBI::dbGetQuery(con,
    "SELECT screen_verdict, screen_version FROM screen
     WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
```

Replace with:

```r
#' @param force Recompute even when the stored row is current.
#' @param review TRUE scores the contig's review copy: hits from
#'   `review_blast_hits`, verdict to `review_screen`. The variant columns are
#'   carried over from the original `screen` row, because variant grouping
#'   compares every contig of a sample and is not rechecked.
#' @return (invisibly) the verdict written, or the stored one when skipped.
#' @export
screen_one <- function(con, id, path = 1L, scaffold = 1L, force = FALSE, review = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  out_tbl <- review_tbl("screen", review)
  cur <- DBI::dbGetQuery(con, sprintf(
    "SELECT screen_verdict, screen_version FROM %s
     WHERE ID = ? AND path = ? AND scaffold = ?", out_tbl), params = list(id, path, scaffold))
```

Find:

```r
  hits <- DBI::dbGetQuery(con,
    "SELECT hit_rank, hit_class, bitscore, lineage_taxids, lineage FROM blast_hits
     WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp", params = list(id, path, scaffold))
```

Replace with:

```r
  hits <- DBI::dbGetQuery(con, sprintf(
    "SELECT hit_rank, hit_class, bitscore, lineage_taxids, lineage FROM %s
     WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp", review_tbl("blast_hits", review)), params = list(id, path, scaffold))
```

Find:

```r
  # one transaction: a failed append after a committed delete would leave the
  # unit with no verdict at all
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con,
      "DELETE FROM screen WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    DBI::dbAppendTable(con, "screen", data.frame(
      ID = id, path = path, scaffold = scaffold,
      screen_verdict = r$verdict, screen_class = r$class,
      screen_confidence = r$confidence, screen_reason = r$reason,
      screen_version = SCREEN_VERSION, lca_taxid = r$lca_taxid,
      lca_taxon = r$lca_taxon, lca_rank = r$lca_rank, lca_n_hits = r$lca_n_hits,
      host_concordant = r$host_concordant, margin = r$margin,
      time_stamp = as.integer(Sys.time()), stringsAsFactors = FALSE))
  })
  invisible(r$verdict)
```

Replace with:

```r
  row <- data.frame(
    ID = id, path = path, scaffold = scaffold,
    screen_verdict = r$verdict, screen_class = r$class,
    screen_confidence = r$confidence, screen_reason = r$reason,
    screen_version = SCREEN_VERSION, lca_taxid = r$lca_taxid,
    lca_taxon = r$lca_taxon, lca_rank = r$lca_rank, lca_n_hits = r$lca_n_hits,
    host_concordant = r$host_concordant, margin = r$margin,
    time_stamp = as.integer(Sys.time()), stringsAsFactors = FALSE)
  if (review) {
    v <- DBI::dbGetQuery(con,
      "SELECT variant_class, variant_group, variant_identity FROM screen
       WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    row[names(v)] <- if (nrow(v) > 0) v[1, ] else NA
  }
  # one transaction: a failed append after a committed delete would leave the
  # unit with no verdict at all
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, sprintf(
      "DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", out_tbl),
      params = list(id, path, scaffold))
    DBI::dbAppendTable(con, out_tbl, row)
  })
  invisible(r$verdict)
```

The `samples` read, `screen_params()`, and the `blast_searched` read between these blocks stay literal. Original mode appends the same columns as before (the variant columns land NULL and `screen_sample()` fills them through `screen_variants()` afterwards); `DBI::dbAppendTable` names the data frame's columns in its INSERT, so it is not positional.

- [ ] **Step 4: Run the screen tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-screen.R")'`
Expected: PASS (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 131 ]`).

- [ ] **Step 5: Give parse_blast_candidates(), derive_blast_products(), and write_blast_screen() the review routing**

In `R/blast_genbank.R`, Find:

```r
#'   than treats as clean.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
parse_blast_candidates <- function(con, id, path = 1L, scaffold = 1L, hits_file, taxdump_dir = NULL) {
  ts <- as.integer(Sys.time())
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  clear_unit <- function() {
    DBI::dbExecute(con,
      "DELETE FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
  }

  h <- parse_blast_hits_file(hits_file)
  if (nrow(h) == 0) {
    DBI::dbWithTransaction(con, clear_unit())
    # blast_hits is now empty for this unit; let derive_blast_products clear
    # the derived tables (and any now-stale assemblies.blast_*) too.
    return(derive_blast_products(con, id, path, scaffold, ts = ts))
  }
```

Replace with:

```r
#'   than treats as clean.
#' @param review TRUE stores the hits of the contig's review copy in
#'   `review_blast_hits` and calls derive_blast_products() in review mode. The
#'   sample's lineage taxids are never filled in.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
parse_blast_candidates <- function(con, id, path = 1L, scaffold = 1L, hits_file, taxdump_dir = NULL,
                                   review = FALSE) {
  ts <- as.integer(Sys.time())
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  hits_tbl <- review_tbl("blast_hits", review)
  clear_unit <- function() {
    DBI::dbExecute(con,
      sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", hits_tbl),
      params = list(id, path, scaffold))
  }

  h <- parse_blast_hits_file(hits_file)
  if (nrow(h) == 0) {
    DBI::dbWithTransaction(con, clear_unit())
    # blast_hits is now empty for this unit; let derive_blast_products clear
    # the derived tables (and any now-stale assemblies.blast_*) too.
    return(derive_blast_products(con, id, path, scaffold, ts = ts, review = review))
  }
```

Find:

```r
  if (is.na(samp %||% NA) || !nzchar(samp %||% "")) {
```

Replace with:

```r
  if (!review && (is.na(samp %||% NA) || !nzchar(samp %||% ""))) {
```

Find:

```r
    DBI::dbAppendTable(con, "blast_hits", hits)
  })
  derive_blast_products(con, id, path, scaffold, ts = ts)
}
```

Replace with:

```r
    DBI::dbAppendTable(con, hits_tbl, hits)
  })
  derive_blast_products(con, id, path, scaffold, ts = ts, review = review)
}
```

Find:

```r
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
  clear_assembly_ref <- function() {
    DBI::dbExecute(con,
      "UPDATE assemblies SET blast_accession = NULL, blast_species = NULL, blast_pident = NULL,
         blast_qcovs = NULL, blast_evalue = NULL, blast_source = NULL
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
  }
  h <- DBI::dbGetQuery(con,
    "SELECT * FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp",
    params = list(id, path, scaffold))
```

Replace with:

```r
#' @param ts Timestamp to stamp the written rows with.
#' @param review TRUE derives the products of the contig's review copy: it
#'   reads `review_blast_hits`, replaces the row in `review_blast_screen`, and
#'   sets or clears the `blast_*` summary in `review_units`. The candidates
#'   only pick the top hit; `blast_ref_candidates` is never written, and the
#'   source label still comes from `blast_searched`.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
derive_blast_products <- function(con, id, path = 1L, scaffold = 1L,
                                  ts = as.integer(Sys.time()), review = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  screen_tbl <- review_tbl("blast_screen", review)
  unit_tbl <- review_tbl("assemblies", review)
  clear <- function() {
    for (tbl in c(if (!review) "blast_ref_candidates", screen_tbl)) {
      DBI::dbExecute(con,
        sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
        params = list(id, path, scaffold))
    }
  }
  clear_assembly_ref <- function() {
    DBI::dbExecute(con, sprintf(
      "UPDATE %s SET blast_accession = NULL, blast_species = NULL, blast_pident = NULL,
         blast_qcovs = NULL, blast_evalue = NULL, blast_source = NULL
       WHERE ID = ? AND path = ? AND scaffold = ?", unit_tbl),
      params = list(id, path, scaffold))
  }
  h <- DBI::dbGetQuery(con, sprintf(
    "SELECT hit_rank, hsp, accession, title, species, staxid, lineage, lineage_taxids,
       rank_distance, pident, qcovs, evalue, bitscore, aln_len, hit_class
     FROM %s WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp", review_tbl("blast_hits", review)),
    params = list(id, path, scaffold))
```

The hit read now names its columns in both modes (no `SELECT *` on a review table or its original): these 15 are every column the rest of the function, `ref_candidate_ok()`, and `write_blast_screen()` use. Every fixture that reaches this read builds `blast_hits` from `BLAST_HITS_DDL` or `new_db()`, which have all of them.

Find:

```r
    write_blast_screen(con, id, path, scaffold, subj[1, ], hits = h, ts = ts)
    if (!is.null(df)) {
      DBI::dbAppendTable(con, "blast_ref_candidates", df)
```

Replace with:

```r
    write_blast_screen(con, id, path, scaffold, subj[1, ], hits = h, ts = ts, tbl = screen_tbl)
    if (!is.null(df)) {
      if (!review) DBI::dbAppendTable(con, "blast_ref_candidates", df)
```

Find:

```r
      DBI::dbExecute(con,
        "UPDATE assemblies SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?,
           blast_evalue = ?, blast_source = ?
         WHERE ID = ? AND path = ? AND scaffold = ?",
```

Replace with:

```r
      DBI::dbExecute(con, sprintf(
        "UPDATE %s SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?,
           blast_evalue = ?, blast_source = ?
         WHERE ID = ? AND path = ? AND scaffold = ?", unit_tbl),
```

The `blast_searched` read between the last two blocks stays literal, so `blast_source` follows the original search in both modes.

Find:

```r
#' replaced.
#' @noRd
write_blast_screen <- function(con, id, path, scaffold, top, hits,
                               ts = as.integer(Sys.time())) {
```

Replace with:

```r
#' replaced.
#'
#' `tbl` is the table appended to: `blast_screen`, or `review_blast_screen`
#' for a review copy.
#' @noRd
write_blast_screen <- function(con, id, path, scaffold, top, hits,
                               ts = as.integer(Sys.time()), tbl = "blast_screen") {
```

Find:

```r
  DBI::dbAppendTable(con, "blast_screen", row)
```

Replace with:

```r
  DBI::dbAppendTable(con, tbl, row)
```

Both transactions in `derive_blast_products()` still open with a DELETE (`clear()`), so they stay write-first in both modes. The callers (`rp_blast.R`, `rp_searched.R`, `score_one()`, `rederive_blast_products()`) pass positional arguments only, so they are unchanged.

- [ ] **Step 6: Run the BLAST and atomicity tests to see them pass**

Run:

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-blast-parse.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-write-atomicity.R")'
```

Expected: PASS (`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 84 ]` and `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 29 ]`).

- [ ] **Step 7: Regenerate the three help pages**

`screen_one`, `parse_blast_candidates`, and `derive_blast_products` are exported and have man pages. NAMESPACE is hand-maintained: never run `devtools::document()`.

Run:

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
git status --short man/ NAMESPACE
```

Expected: exactly ` M man/derive_blast_products.Rd`, ` M man/parse_blast_candidates.Rd`, and ` M man/screen_one.Rd`. Their usage now ends in `review = FALSE`, each with a new `\item{review}{...}`. roxygen prints "Could not resolve link to topic" warnings for internal topics; they were there before this task. If any other man/ file is listed, put it back (`git checkout -- <file>` for a modified one, `rm <file>` for a new one) and mention it in your report; NAMESPACE must not be listed.

- [ ] **Step 8: Run the affected test files**

Run each:

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-screen.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-blast-parse.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-write-atomicity.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-screen-variants.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-taxonomy-assign.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-select-paths.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-scientific-defaults.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-units.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-no-duplicate-definitions.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-load-assemblies.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-bin-contract.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-db-dist.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-pipeline-drivers.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-main-run.R")'
```

Expected: 0 failures in every file. Known counts: test-screen.R 131, test-blast-parse.R 84, test-write-atomicity.R 29, test-screen-variants.R 132, test-taxonomy-assign.R 17, test-select-paths.R 30, test-scientific-defaults.R 34, test-export-units.R 53, test-db-dist.R 83, test-docs-site.R 237, test-no-duplicate-definitions.R 1. What the others cover: test-docs-site.R checks the regenerated Rd for raw backticks and `[fn()]` links; test-screen-variants.R and test-taxonomy-assign.R run original-mode `screen_one()` through `screen_sample()` and `rescreen_samples()`; test-load-assemblies.R and test-select-paths.R run `score_one()`, the in-package caller of `derive_blast_products()`; `NOT_CRAN=true` makes test-blast-parse.R and test-db-dist.R run `rp_searched.R`, test-pipeline-drivers.R run `rp_score.R`, and test-nf-main-run.R run the real main.nf with fake tools (a few minutes; it skips when Nextflow, sqlite3, or jq is missing).

- [ ] **Step 9: Commit**

```bash
git add R/screen.R R/blast_genbank.R man/screen_one.Rd man/parse_blast_candidates.Rd man/derive_blast_products.Rd tests/testthat/test-screen.R tests/testthat/test-blast-parse.R tests/testthat/test-write-atomicity.R
git commit -m "Review mode for BLAST products and screen"
```


### Task 8: Review mode: chimera and stock ITSx

**Files:**
- Modify: `R/chimera.R` (`chimera_unit()` and its roxygen block, lines 437-552)
- Modify: `R/itsx_check.R` (`itsx_check_units()`, `unit_sequences()`, `itsx_check_one()`, as Task 5 left them)
- Modify: `man/chimera_unit.Rd`, `man/itsx_check_one.Rd` (regenerated by roxygen, never edited by hand)
- Test: `tests/testthat/test-chimera.R` (append)
- Test: `tests/testthat/test-itsx-check.R` (append)

**Interfaces:**
- Consumes:
  - Task 1: `new_db()` creates `review_units`, the 7 mirror tables (this task uses `review_annotations`, `review_blast_hits`, `review_chimera`, `review_chimera_windows`), and the views, so `chimera_test_con()` has them.
  - Task 2: `review_copy(con, id, path, scaffold)` (tests only; copies the original's `chimera`, `chimera_windows`, and `annotations` rows into their mirrors).
  - Task 3: `review_coverage(con, id, path, scaffold)` -> `list(depths, clips)`: numeric tracks of the original coverage row cut to the copy's `[trim_start, trim_end]`, NULL for a missing row or track.
  - Task 5: `itsx_store_warnings(con, units, chk, review = FALSE)` (writes `review_units` when TRUE, raising errors); `itsx_run_finished(out_prefix)`; `itsx_check_one()` storing `not_run(ITSx did not finish)`; the test helper `itsx_out(pre, positions, results, finished = TRUE, err = NULL)` near the end of test-itsx-check.R.
  - Task 6: `review_tbl(name, review = FALSE)` (R/review.R).
- Produces:
  - `chimera_unit(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL, label_db = NULL, threads = 1L, workdir = tempfile("chimera_"), review = FALSE)` (exported). Review mode: sequence from `review_units`, features from `review_annotations`, hits from `review_blast_hits`, `chimeric` from `review_units`, depth and clip tracks from `review_coverage()`; writes `review_chimera` and `review_chimera_windows` (DELETE first, in one transaction). With no copy it returns NULL and writes nothing, as it does today for a unit with no sequence.
  - `unit_sequences(con, units, review = FALSE)`: `review_units.sequence` when TRUE.
  - `itsx_check_units(con, units, out, review = FALSE)`: features from `review_annotations` and sequences via `unit_sequences(review = TRUE)` when TRUE; output columns unchanged from Task 5.
  - `itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix, review = FALSE)` (exported): passes `review` to `itsx_check_units()` and `itsx_store_warnings()`, so the result lands on the copy's `review_units` row and the annotate row is untouched.

The recheck scores a review copy with the same chimera and stock-ITSx code the main pipeline uses (SPEC S7.3 steps 7 and 10). `chimera_unit()` gains `review = FALSE`: in review mode it reads the copy's sequence, features, hits, and ITSx chimera flag through `review_tbl()`, takes depth and clip tracks from `review_coverage()` (already numeric and cut to the copy's window, S4.5, so the text `track()` parser is bypassed), and writes the mirrors; samples stay literal. `itsx_check_one()`, `itsx_check_units()`, and `unit_sequences()` gain the same argument, so the stock check flags the copy's sequence and features and stores on the copy (Task 5 already gave `itsx_store_warnings()` its `review` argument). The write blocks keep today's shape: reads before the transaction, a DELETE first inside it. Intermediate writes are not guarded by `rev` (S5.4); `recheck_unit()` (Task 10) guards the final write.

- [ ] **Step 1: Write the failing chimera test**

Append this at the end of `tests/testthat/test-chimera.R`:

```r

test_that("chimera_unit(review = TRUE) scores the copy and writes only its mirrors", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  orig <- withr::with_seed(3, paste(sample(c("A", "C", "G", "T"), 4000, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence) VALUES ('S1', 1, 1, 4000, ?),
    ('S1', 1, 2, 4000, ?)", params = list(orig, orig))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 2)")
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = '33208;7711;186623;8030;8016;8017;8018' WHERE ID = 'S1'")
  DBI::dbExecute(con, "UPDATE annotate SET chimeric = 0")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop)
    VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1',1,1,?,?)",
    params = list(paste(rep(20, 4000), collapse = ","), paste(c(rep(0, 2199), 9, rep(0, 1800)), collapse = ",")))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, reason) VALUES ('S1',1,1,0,'clean','original')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident)
    VALUES ('S1',1,1,1,500,99), ('S1',1,1,3501,4000,99)")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  # the copy: first 1000 bp trimmed off (18S cut to 1..800), chimeric per its own ITSx run, its own hits
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?, length = 3000, trim_start = 1001, trim_end = 4000,
    chimeric = 1", params = list(substr(orig, 1001, 4000)))
  DBI::dbExecute(con, "UPDATE review_annotations SET start = max(start - 1000, 1), stop = stop - 1000")
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  DBI::dbAppendTable(con, "review_blast_hits", data.frame(ID = "S1", path = 1L, scaffold = 1L,
    hit_rank = 1:6, hsp = 1L, accession = paste0("A", 1:6), species = "Salmo salar",
    lineage = "Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar",
    lineage_taxids = fish, pident = 99, bitscore = 900,
    qstart = rep(c(1L, 1500L), each = 3), qend = rep(c(800L, 3000L), each = 3)))
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    writeLines(sprintf("NC_1\tSalmo salar 18S\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t1800\t%s", ids), out)
    0L
  }
  fake_rt <- function(fasta, outdir, ...) {
    writeLines(c("#idx target class strnd p/f unexpected_features", "1 u18S SSU.Eukarya plus PASS -",
      "2 u28S LSU.Eukarya plus PASS -"), file.path(outdir, "rt.ribotyper.short.out")); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = 8030L,
    lineage = "Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar", lineage_taxids = fish)
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, chimera_run_ribotyper = fake_rt,
    tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  originals <- function() list(DBI::dbGetQuery(con, "SELECT ID, path, scaffold, score, label, reason FROM chimera"),
    DBI::dbGetQuery(con, "SELECT ID, path, scaffold, start, stop FROM chimera_windows ORDER BY start"))
  before <- originals()

  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel", review = TRUE)
  expect_identical(originals(), before)
  ch <- jsonlite::fromJSON(r$channels)
  expect_equal(ch$reads$junctions, 1200L)
  expect_equal(ch$ribotyper$score, 0.2)
  expect_equal(jsonlite::fromJSON(r$gene_taxa)[, c("element", "n_hits")],
    data.frame(element = c("18S", "28S"), n_hits = c(3L, 3L)))
  rc <- DBI::dbGetQuery(con, "SELECT junctions, reason FROM review_chimera")
  expect_equal(rc$junctions, "1200")
  expect_equal(rc$reason, r$reason)
  w <- DBI::dbGetQuery(con, "SELECT start, stop, clip_frac FROM review_chimera_windows ORDER BY start")
  expect_equal(w$start, RiboPilot:::chimera_tile(3000)$start)
  expect_equal(max(w$stop), 3000L)
  expect_equal(w$start[w$clip_frac > 0], c(751L, 1001L))

  # no copy: nothing to score in review mode
  expect_null(RiboPilot:::chimera_unit(con, "S1", 1L, 2L, dbarg = "panel", review = TRUE))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_chimera WHERE scaffold = 2")$n, 0L)

  # original mode scores the original and leaves the copy alone
  RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel")
  expect_equal(DBI::dbGetQuery(con, "SELECT junctions, reason FROM review_chimera"), rc)
  expect_equal(DBI::dbGetQuery(con, "SELECT junctions FROM chimera")$junctions, "2200")
  expect_equal(DBI::dbGetQuery(con, "SELECT start FROM chimera_windows ORDER BY start")$start,
    RiboPilot:::chimera_tile(4000)$start)
  expect_equal(DBI::dbGetQuery(con, "SELECT start, stop, clip_frac FROM review_chimera_windows ORDER BY start"), w)
})
```

Why these values: the original has a clip pileup at 2200, so the copy (trimmed by 1000 bp on the left) must report it at 1200 from the cut coverage; reading the uncut track gives 2200. The copy's 18S is 1..800 and its own six hits sit on 1..800 and 1500..3000, so each element resolves from 3 hits only when both the copy's features and the copy's hits are read (the original has no hits, and its 18S and 28S are too long for those hits to cover half). `review_units.chimeric = 1` against `annotate.chimeric = 0` gives the ribotyper channel its 0.2 ITSx score only from the copy. Windows tile the 3000 bp copy, not the 4000 bp original. The copied `review_chimera` and `review_chimera_windows` rows share keys with the new ones, so a DELETE that misses a mirror fails the append on the primary key, and one that hits an original changes `originals()`. Contig 2 has an assemblies and an annotate row but no copy, so `review_units_v` would show it; review mode must read `review_units` and find nothing. After the original-mode run, the original windows tile 4000 bp and the mirror windows are unchanged. `chimera_run_ribotyper` is mocked so the test behaves the same with or without ribotyper on PATH.

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R")'`
Expected: FAIL, `[ FAIL 1 | WARN 0 | SKIP 2 | PASS 115 ]`: the new test errors with `unused argument (review = TRUE)` at its first `chimera_unit()` call. The two skips are the BLAST+ tests ("chimera_label_prefix builds a blast db for a labels fasta" and "synthetic fixtures get the expected label and a junction within tolerance"); on a machine with BLAST+ on PATH they run instead, and the PASS count is higher.

- [ ] **Step 3: Give chimera_unit() the review argument**

In `R/chimera.R`, Find:

```r
#' @param workdir Scratch directory for intermediate files.
#' @return (invisibly) the written `chimera` row, as a list.
#' @export
chimera_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                         label_db = NULL, threads = 1L, workdir = tempfile("chimera_")) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  label_db <- chimera_label_prefix(label_db, workdir)
  seq <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$sequence[1]
  if (is.na(seq %||% NA) || !nzchar(seq)) return(invisible(NULL))
  feats <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  hits <- tryCatch(DBI::dbGetQuery(con, "SELECT accession, bitscore, qstart, qend, lineage_taxids, lineage, species, pident FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  cov <- tryCatch(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  itsx <- tryCatch(DBI::dbGetQuery(con, "SELECT chimeric FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$chimeric[1], error = function(e) NA)
  track <- function(x) if (length(x) && !is.na(x) && nzchar(x)) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]) else NULL
```

Replace with:

```r
#' @param workdir Scratch directory for intermediate files.
#' @param review TRUE checks the contig's review copy: its sequence, features,
#'   hits, and ITSx chimera flag come from the review tables, coverage is cut
#'   to the copy's window, and the result goes to `review_chimera` and
#'   `review_chimera_windows`.
#' @return (invisibly) the written `chimera` row, as a list.
#' @export
chimera_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                         label_db = NULL, threads = 1L, workdir = tempfile("chimera_"), review = FALSE) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  label_db <- chimera_label_prefix(label_db, workdir)
  seq <- DBI::dbGetQuery(con, sprintf("SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("assemblies", review)), params = list(id, path, scaffold))$sequence[1]
  if (is.na(seq %||% NA) || !nzchar(seq)) return(invisible(NULL))
  feats <- DBI::dbGetQuery(con, sprintf("SELECT element, start, stop FROM %s WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    review_tbl("annotations", review)), params = list(id, path, scaffold))
  hits <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT accession, bitscore, qstart, qend, lineage_taxids, lineage, species, pident FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("blast_hits", review)), params = list(id, path, scaffold)), error = function(e) data.frame())
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  itsx <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT chimeric FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("annotate", review)), params = list(id, path, scaffold))$chimeric[1], error = function(e) NA)
  track <- function(x) if (length(x) && !is.na(x) && nzchar(x)) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]) else NULL
  cov <- if (review) review_coverage(con, id, path, scaffold) else {
    cv <- tryCatch(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold)), error = function(e) data.frame())
    list(depths = if (nrow(cv)) track(cv$depths[1]), clips = if (nrow(cv) && "clips" %in% names(cv)) track(cv$clips[1]))
  }
```

`cov` is now `list(depths, clips)` in both modes (`if` without `else` gives NULL); coverage stays a literal table read in original mode and is never routed through `review_tbl()`.

Find:

```r
  depths <- if (nrow(cov)) track(cov$depths[1]) else NULL
  clips <- if (nrow(cov) && "clips" %in% names(cov)) track(cov$clips[1]) else NULL
```

Replace with:

```r
  depths <- cov$depths
  clips <- cov$clips
```

Find:

```r
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbAppendTable(con, "chimera", row)
    if (!is.null(windows) && nrow(windows)) {
      windows$ID <- id; windows$path <- path; windows$scaffold <- scaffold
      DBI::dbAppendTable(con, "chimera_windows", windows[, c("ID", "path", "scaffold", "start", "stop", "accession",
```

Replace with:

```r
  ch_tbl <- review_tbl("chimera", review)
  win_tbl <- review_tbl("chimera_windows", review)
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", ch_tbl),
      params = list(id, path, scaffold))
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", win_tbl),
      params = list(id, path, scaffold))
    DBI::dbAppendTable(con, ch_tbl, row)
    if (!is.null(windows) && nrow(windows)) {
      windows$ID <- id; windows$path <- path; windows$scaffold <- scaffold
      DBI::dbAppendTable(con, win_tbl, windows[, c("ID", "path", "scaffold", "start", "stop", "accession",
```

`DBI::dbAppendTable` names every column of the data frame in its INSERT, so it is not positional. The block stays write-first (DELETE) with all reads and computation before it, and `chimera_unit()` is never called inside an open transaction (the recheck opens none), so `DBI::dbWithTransaction` stays.

- [ ] **Step 4: Run the chimera tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 2 | PASS 130 ]` (same two BLAST+ skips as Step 2). The three existing `chimera_unit` tests (mocked tools, label_db fallback, no_taxonomy) still pass in original mode.

- [ ] **Step 5: Write the failing stock-ITSx test**

Append this at the end of `tests/testthat/test-itsx-check.R` (after Task 5's tests; it uses Task 5's `itsx_out()` helper defined above it in the same file):

```r

test_that("itsx_check_one(review = TRUE) checks the copy's sequence and features and stores on the copy", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET lineage = 'Metazoa;Chordata' WHERE ID = 'S1'")
  s <- withr::with_seed(5, paste(sample(c("A", "C", "G", "T"), 1000, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence) VALUES ('S1', 1, 1, 1000, ?)",
    params = list(s))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES
    ('S1',1,1,'18S',1,300), ('S1',1,1,'ITS1',301,500), ('S1',1,1,'5.8S',501,660),
    ('S1',1,1,'ITS2',661,800), ('S1',1,1,'28S',801,1000)")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_chimeric', itsx_origin = 'F',
    itsx_scores = 'F: 1 1e-10 40'")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  # the copy: an N-run inside ITS2, and an ITS1 that ends 50 bp before the one ITSx finds
  m <- paste0(substr(s, 1, 700), strrep("N", 20), substr(s, 721, 1000))
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?", params = list(m))
  DBI::dbExecute(con, "UPDATE review_annotations SET stop = 450 WHERE element = 'ITS1'")
  units <- data.frame(ID = "S1", path = 1L, scaffold = c(1L, 2L))
  expect_equal(RiboPilot:::unit_sequences(con, units, review = TRUE), c(m, NA))
  expect_equal(RiboPilot:::unit_sequences(con, units), c(s, NA))

  get <- function(tbl) unlist(DBI::dbGetQuery(con, sprintf("SELECT itsx_warnings, itsx_origin, itsx_scores
    FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = 1", tbl)))
  sc <- "F: 4 1e-08 33.15, M: 4 2.9e-16 49.9"
  pre <- itsx_out(file.path(withr::local_tempdir(), "unit.itsx_check"),
    "S1.1.1\t1000 bp.\tSSU: 1-300\tITS1: 301-500\t5.8S: 501-660\tITS2: 661-800\tLSU: 801-1000\t",
    paste0("S1.1.1\t1000\tM\t0\t4\t2.9e-16\t49.9\t49.9\t1\t1000\t1_SSU\t4_LSU\t\t1_SSU,2_5.8,3_End,4_LSU\t", sc))
  flags <- "itsx_gap_in_ITS;itsx_boundary_shift(ITS1:50bp)"
  chk <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, pre, review = TRUE)
  expect_equal(chk$itsx_check, flags)
  expect_equal(get("review_units"), c(itsx_warnings = flags, itsx_origin = "M", itsx_scores = sc))
  expect_equal(get("annotate"), c(itsx_warnings = "itsx_chimeric", itsx_origin = "F", itsx_scores = "F: 1 1e-10 40"))

  # original mode checks the original and leaves the copy alone
  expect_equal(RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, pre)$itsx_check, "")
  expect_equal(get("annotate"), c(itsx_warnings = "", itsx_origin = "M", itsx_scores = sc))
  expect_equal(get("review_units")[["itsx_warnings"]], flags)

  # a run that did not finish is not_run on the copy only
  RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, file.path(dirname(pre), "missing"), review = TRUE)
  ru <- get("review_units")
  expect_equal(ru[["itsx_warnings"]], "not_run(ITSx did not finish)")
  expect_true(is.na(ru[["itsx_origin"]]) && is.na(ru[["itsx_scores"]]))
  expect_equal(get("annotate")[["itsx_warnings"]], "")
})
```

Why these values: ITSx reports ITS1 301-500, 5.8S 501-660, ITS2 661-800, which match the original's features exactly, so the original is clean. Only the copy's sequence has the 20-N run inside the ITS span (`itsx_gap_in_ITS`), and only the copy's ITS1 ends 50 bp early (`itsx_boundary_shift(ITS1:50bp)`), so reading the original sequence or the original features drops one flag. The origin M matches the Metazoa lineage, so no origin flag. The annotate and `review_units` rows start with different stored values, so a store on the wrong row shows.

- [ ] **Step 6: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R")'`
Expected: FAIL, `[ FAIL 1 | WARN 0 | SKIP 1 | PASS 113 ]`: the new test errors with `unused argument (review = TRUE)` in `RiboPilot:::unit_sequences(con, units, review = TRUE)`. The skip is the real-ITSx test (no RIBOPILOT_TEST_PROJECT).

- [ ] **Step 7: Give the stock-ITSx check the review argument**

In `R/itsx_check.R` (as Task 5 left it), Find:

```r
#'   and `message` when the run itself failed.
#' @return data.frame(seqid, itsx_check, itsx_origin, itsx_scores, itsx_ITS1,
```

Replace with:

```r
#'   and `message` when the run itself failed.
#' @param review TRUE reads each unit's sequence and features from its review copy.
#' @return data.frame(seqid, itsx_check, itsx_origin, itsx_scores, itsx_ITS1,
```

Find:

```r
itsx_check_units <- function(con, units, out) {
```

Replace with:

```r
itsx_check_units <- function(con, units, out, review = FALSE) {
```

Find:

```r
  seqs <- unit_sequences(con, units)
  its <- parse_itsx(out$positions)
```

Replace with:

```r
  seqs <- unit_sequences(con, units, review = review)
  its <- parse_itsx(out$positions)
```

Find:

```r
    ann <- DBI::dbGetQuery(con,
      "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
```

Replace with:

```r
    ann <- DBI::dbGetQuery(con, sprintf(
      "SELECT element, start, stop FROM %s WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      review_tbl("annotations", review)),
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
```

The sample lineage read just above it stays on `samples` in both modes.

Find:

```r
#' Scaffold sequence of each unit (NA when absent)
#' @noRd
unit_sequences <- function(con, units) {
  vapply(seq_len(nrow(units)), function(i) {
    r <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
```

Replace with:

```r
#' Sequence of each unit (NA when absent): the contig's, or its review copy's when `review` is TRUE
#' @noRd
unit_sequences <- function(con, units, review = FALSE) {
  vapply(seq_len(nrow(units)), function(i) {
    r <- DBI::dbGetQuery(con, sprintf(
      "SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
      review_tbl("assemblies", review)),
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
```

Find:

```r
#' @param out_prefix ITSx output prefix (files `<out_prefix>.positions.txt` etc.).
#' @return one-row data.frame as from `itsx_check_units()`, invisibly.
#' @export
itsx_check_one <- function(con, id, path = 1L, scaffold = 1L, out_prefix) {
```

Replace with:

```r
#' @param out_prefix ITSx output prefix (files `<out_prefix>.positions.txt` etc.).
#' @param review TRUE checks the contig's review copy (its sequence and
#'   features) and stores the result on the copy instead of the annotate row.
#' @return one-row data.frame as from `itsx_check_units()`, invisibly.
#' @export
itsx_check_one <- function(con, id, path = 1L, scaffold = 1L, out_prefix, review = FALSE) {
```

Find (Task 5's body of `itsx_check_one`):

```r
    list(ok = FALSE, message = "ITSx did not finish")
  }
  chk <- itsx_check_units(con, units, out)
  itsx_store_warnings(con, units, chk)
```

Replace with:

```r
    list(ok = FALSE, message = "ITSx did not finish")
  }
  chk <- itsx_check_units(con, units, out, review = review)
  itsx_store_warnings(con, units, chk, review = review)
```

`itsx_export_check()` further down keeps its original-mode calls (`unit_sequences(con, units)`, `itsx_check_units(con, units, out)`); Task 18 deletes it. `itsx_store_not_run()` needs no change: its not_run path returns before any sequence or feature read.

- [ ] **Step 8: Run the ITSx tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 124 ]` (the real-ITSx test still skips without RIBOPILOT_TEST_PROJECT).

- [ ] **Step 9: Regenerate the two help pages**

`chimera_unit` and `itsx_check_one` are exported and their signatures changed. NAMESPACE is hand-maintained: never run `devtools::document()`.

Run:

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
git status --short man/ NAMESPACE
```

Expected: exactly ` M man/chimera_unit.Rd` and ` M man/itsx_check_one.Rd`. Their usage now ends `workdir = tempfile("chimera_"),` / `review = FALSE` and `itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix, review = FALSE)`, each with a new `\item{review}{...}` (roxygen turns the backticks into `\code{...}`). roxygen prints "Could not resolve link to topic" warnings for internal topics; they were there before this task. If any other man/ file is listed, put it back (`git checkout -- <file>` for a modified one, `rm <file>` for a new one) and mention it in your report; NAMESPACE must not be listed.

- [ ] **Step 10: Run every test file this task touches or could break**

Run:

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("chimera", "itsx-check", "docs-site", "no-duplicate-definitions", "nf-bin-contract", "pipeline-drivers", "reconcile")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'
```

Expected: every file reports `FAIL 0`. Tallies at the time of writing (no BLAST+ on PATH): chimera `SKIP 2 | PASS 130`, itsx-check `SKIP 1 | PASS 124`, docs-site `PASS 237` (checks the regenerated Rd for raw backticks and `[fn()]` links), no-duplicate-definitions `PASS 1`, nf-bin-contract `PASS 162` (rp_chimera.R and rp_reconcile.R still call the functions with their old arguments), pipeline-drivers `PASS 17` (runs rp_reconcile.R, so `itsx_check_one` in original mode, through the bin script), reconcile `PASS 52`.

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-main-run.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 33 ]`. It runs the real main.nf with the fake tools, so rp_chimera.R scores every contig through the reworked original-mode coverage read; it takes about 4 minutes.

- [ ] **Step 11: Commit**

```bash
git add R/chimera.R R/itsx_check.R man/chimera_unit.Rd man/itsx_check_one.Rd tests/testthat/test-chimera.R tests/testthat/test-itsx-check.R
git commit -m "Review mode for chimera and stock ITSx check"
```


### Task 9: Shared rRNA tool script

**Files:**
- Create: `inst/nextflow/bin/rp_rrna_tools.sh` (git mode 100755)
- Modify: `inst/nextflow/modules/pipeline_r.nf` (process `annotate_tools`, lines 167-187)
- Modify: `inst/nextflow/modules/find_rrna.nf` (comment, line 32)
- Modify: `R/annotate.R` (roxygen of `rrna_kingdom_for`, lines 240-241)
- Modify: `tests/testthat/test-rrna-kingdom.R` (lines 157-170)
- Modify: `tests/testthat/test-nextflow.R` (lines 1583-1584)
- Test: `tests/testthat/test-rrna-tools.R` (new)

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `rp_rrna_tools.sh <fasta> <tool> <kingdom> <itsx_taxon> <itsx_opts> <rrna_opts> <cpus>`, called by name (Nextflow puts `inst/nextflow/bin` on the task PATH). Exactly 7 arguments, else a `usage:` line on stderr and exit 2 before any tool runs, so callers pass `"$ITSXO"` and `"$RRNAO"` quoted (either may be empty). Empty tool, kingdom, and taxon default to pybarrnap, euk, and all. Writes into the working directory: `unit.barrnap.gff`; `unit.itsx.positions.txt` (empty file when the tuned run wrote none); the stock run's `unit.itsx_check.*`, including `unit.itsx_check.positions.txt` (empty file when none) and `unit.itsx_check.err` (the stock run's stderr, always created); `unit.tools_status`, exactly three lines `barrnap=<exit>`, `itsx=<exit>`, `itsx_check=<exit>` (a tool is OK when its value is 0). Runs all three tools, then exits with the barrnap/pybarrnap exit status.
  - `annotate_tools` keeps its options read (pipeline_r.nf:155-165, from `read TOOL KINGDOM ITSXT ITSXO <<<` through the `RRNAO` query) byte for byte, for Task 11 to copy; its tool call is the single line `rp_rrna_tools.sh unit.fasta "\$TOOL" "\$KINGDOM" "\$ITSXT" "\$ITSXO" "\$RRNAO" ${task.cpus}`. Output tuple and `reconcile`'s input unchanged; `unit.tools_status` is not emitted.

`annotate_tools` runs barrnap/pybarrnap and both ITSx runs inline today (pipeline_r.nf:167-187). This task moves those commands, unchanged, into `inst/nextflow/bin/rp_rrna_tools.sh`, so the recheck (Task 11) runs exactly the commands the main pipeline runs (S7.1, S10). The script records each tool's exit status in `unit.tools_status` (read by `recheck_unit`, Task 10) and keeps the stock run's stderr in `unit.itsx_check.err`, which the existing `unit.itsx_check*` output glob already hands to `reconcile` (read by `itsx_run_finished`, Task 5). It exits with the barrnap/pybarrnap status and `annotate_tools` calls it without `||`, so a barrnap failure still fails the task and an out-of-memory exit still reaches the nextflow.config retry (Review Focus 3).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-rrna-tools.R`:

```r
# inst/nextflow/bin/rp_rrna_tools.sh runs barrnap or pybarrnap, the tuned ITSx
# run, and the stock ITSx check on one contig. The tools are stubs on PATH;
# each stub appends its argv to <name>.calls, one call per line, tab-separated.

rrna_sh <- function() system.file("nextflow", "bin", "rp_rrna_tools.sh", package = "RiboPilot")

mk_tool_stub <- function(dir, name, body = character(0)) {
  calls <- shQuote(file.path(dir, paste0(name, ".calls")))
  p <- file.path(dir, name)
  writeLines(c("#!/usr/bin/env bash",
    sprintf("printf '%%s\\t' \"$@\" >> %s; echo >> %s", calls, calls), body), p)
  Sys.chmod(p, "0755")
  p
}

stub_calls <- function(dir, name) {
  f <- file.path(dir, paste0(name, ".calls"))
  if (!file.exists(f)) return(list())
  lapply(readLines(f), function(l) strsplit(l, "\t", fixed = TRUE)[[1]])
}

run_rrna_sh <- function(dir, args) {
  processx::run("bash", c(rrna_sh(), args), wd = dir, error_on_status = FALSE,
    env = c("current", PATH = paste(dir, Sys.getenv("PATH"), sep = .Platform$path.sep)))
}

mk_unit_fasta <- function(dir) {
  writeLines(c(">S1.1.1", strrep("ACGT", 50)), file.path(dir, "unit.fasta"))
}

ITSX_OUT <- 'while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; esac; shift; done'

test_that("rp_rrna_tools.sh is executable and refuses a wrong argument count", {
  skip_on_os("windows")
  expect_equal(unname(file.access(rrna_sh(), 1)), 0L)
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap")
  mk_tool_stub(dir, "ITSx")
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "--a", "b", "", "2"))
  expect_equal(res$status, 2L)
  expect_match(res$stderr, "usage:", fixed = TRUE)
  expect_length(stub_calls(dir, "pybarrnap"), 0)
  expect_false(file.exists(file.path(dir, "unit.tools_status")))
})

test_that("rp_rrna_tools.sh runs pybarrnap and both ITSx runs with the pipeline's arguments", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap",
    "printf 'S1.1.1\\tbarrnap:0.9\\trRNA\\t1\\t180\\t0\\t+\\t.\\tName=18S_rRNA\\n'")
  mk_tool_stub(dir, "ITSx", c(ITSX_OUT,
    '[ "$out" = unit.itsx ] && echo "S1.1.1 tuned" > "$out.positions.txt"',
    'echo "$out note" >&2', "exit 0"))
  # empty tool, kingdom, and taxon take the defaults; both option strings word-split
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "--allow_single_domain 1e-9,0 -E 1e-3",
    "--evalue 1e-5 --lencutoff 0.5", "2"))
  expect_equal(res$status, 0L)
  expect_equal(stub_calls(dir, "pybarrnap"), list(c("-k", "euk", "-t", "2", "-q",
    "--evalue", "1e-5", "--lencutoff", "0.5", "unit.fasta")))
  expect_equal(stub_calls(dir, "ITSx"), list(
    c("-i", "unit.fasta", "-o", "unit.itsx", "--save_regions", "none", "--graphical", "F",
      "--cpu", "2", "-t", "all", "--preserve", "T", "--allow_single_domain", "1e-9,0", "-E", "1e-3"),
    c("-i", "unit.fasta", "-o", "unit.itsx_check", "--preserve", "T", "--detailed_results", "T",
      "--save_regions", "none", "--graphical", "F", "--cpu", "2")))
  expect_match(readLines(file.path(dir, "unit.barrnap.gff")), "Name=18S_rRNA", fixed = TRUE)
  expect_equal(readLines(file.path(dir, "unit.itsx.positions.txt")), "S1.1.1 tuned")
  expect_equal(readLines(file.path(dir, "unit.itsx_check.positions.txt")), character(0))
  # only the stock run's stderr is kept
  expect_equal(readLines(file.path(dir, "unit.itsx_check.err")), "unit.itsx_check note")
  expect_match(res$stderr, "unit.itsx note", fixed = TRUE)
  expect_equal(readLines(file.path(dir, "unit.tools_status")),
    c("barrnap=0", "itsx=0", "itsx_check=0"))
})

test_that("a barrnap failure sets the exit status after both ITSx runs", {
  skip_on_os("windows")
  for (case in list(list(tool = "barrnap", code = 1L), list(tool = "pybarrnap", code = 137L))) {
    dir <- withr::local_tempdir()
    mk_unit_fasta(dir)
    mk_tool_stub(dir, case$tool, sprintf("exit %d", case$code))
    mk_tool_stub(dir, "ITSx")
    res <- run_rrna_sh(dir, c("unit.fasta", case$tool, "euk", "", "", "--evalue 1e-5", "1"))
    expect_equal(res$status, case$code, info = case$tool)
    expect_equal(stub_calls(dir, case$tool)[[1]], if (case$tool == "barrnap") {
      c("--kingdom", "fun", "--quiet", "--evalue", "1e-5", "unit.fasta")
    } else {
      c("-k", "euk", "-t", "1", "-q", "--evalue", "1e-5", "unit.fasta")
    })
    expect_length(stub_calls(dir, "ITSx"), 2)
    expect_equal(readLines(file.path(dir, "unit.tools_status")),
      c(paste0("barrnap=", case$code), "itsx=0", "itsx_check=0"), info = case$tool)
    expect_true(file.exists(file.path(dir, "unit.itsx.positions.txt")))
    expect_true(file.exists(file.path(dir, "unit.itsx_check.positions.txt")))
  }
})

test_that("an ITSx failure is recorded but never fails the script", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap")
  mk_tool_stub(dir, "ITSx", c(ITSX_OUT,
    '[ "$out" = unit.itsx ] && exit 2',
    'echo "FATAL ERROR: hmmscan failed" >&2', "exit 3"))
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "", "", "1"))
  expect_equal(res$status, 0L)
  expect_equal(readLines(file.path(dir, "unit.tools_status")),
    c("barrnap=0", "itsx=2", "itsx_check=3"))
  expect_equal(readLines(file.path(dir, "unit.itsx.positions.txt")), character(0))
  expect_equal(readLines(file.path(dir, "unit.itsx_check.positions.txt")), character(0))
  expect_equal(readLines(file.path(dir, "unit.itsx_check.err")), "FATAL ERROR: hmmscan failed")
})

test_that("annotate_tools hands the tools to rp_rrna_tools.sh and keeps its outputs", {
  nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))
  start <- grep("^process annotate_tools \\{", nf)
  expect_length(start, 1)
  body <- nf[start:(start - 1 + which(nf[start:length(nf)] == "}")[1])]
  expect_equal(trimws(grep("^\\s*rp_rrna_tools\\.sh ", body, value = TRUE)),
    'rp_rrna_tools.sh unit.fasta "\\$TOOL" "\\$KINGDOM" "\\$ITSXT" "\\$ITSXO" "\\$RRNAO" ${task.cpus}')
  expect_false(any(grepl("barrnap --kingdom|pybarrnap -k|ITSx -i|tools_status", body)))
  expect_true(any(grepl('path("unit.itsx_check*")', body, fixed = TRUE)))
})

test_that("annotate_tools fails on a barrnap error and retries an out-of-memory exit", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  skip_if(!nzchar(Sys.which("sqlite3")), "sqlite3 not installed")
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  con <- chimera_test_con("S1")
  db <- con@dbname
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence)
    VALUES ('S1', 1, 1, 'ACGTACGTAC'), ('S1', 1, 2, 'ACGTACGTAC')")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold, annotate_opts)
    VALUES ('S1', 1, 1, 'default'), ('S1', 1, 2, 'default')")
  DBI::dbDisconnect(con)

  root <- withr::local_tempdir()
  seen <- shQuote(file.path(root, "seen"))
  # contig 1 is killed for memory on its first try; contig 2 always fails
  mk_tool_stub(root, "pybarrnap", c(
    'u=$(head -n 1 "${@: -1}" | cut -c2-)',
    sprintf('echo "$u" >> %s', seen),
    sprintf('n=$(grep -cx "$u" %s)', seen),
    'case "$u.$n" in S1.1.1.1) exit 137 ;; S1.1.2.*) exit 1 ;; esac'))
  mk_tool_stub(root, "ITSx")
  writeLines(c(
    sprintf("include { annotate_tools } from '%s'", file.path(nf_dir, "modules", "pipeline_r.nf")),
    "workflow {",
    "    annotate_tools(Channel.of(tuple('S1', 1, 1, 1, 1), tuple('S1', 1, 2, 1, 1)))",
    "    annotate_tools.out.view { \"EMITTED ${it[0]}.${it[1]}.${it[2]}\" }",
    "}"), file.path(root, "main.nf"))

  res <- processx::run("nextflow", c("run", file.path(root, "main.nf"),
    "-c", file.path(nf_dir, "nextflow.config"), "-work-dir", file.path(root, "work"),
    "-ansi-log", "false", "--db", db), wd = root,
    env = c(RiboPilot:::nf_env(),
      PATH = paste(root, file.path(nf_dir, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(sort(readLines(file.path(root, "seen"))), c("S1.1.1", "S1.1.1", "S1.1.2"))
  expect_match(res$stdout, "EMITTED S1.1.1", fixed = TRUE)
  expect_no_match(res$stdout, "EMITTED S1.1.2", fixed = TRUE)
})
```

In `tests/testthat/test-rrna-kingdom.R`, the test that read the kingdom `case` lines out of `annotate_tools` breaks on purpose (the lines move into the script); it now runs the script and compares with `rrna_kingdom_for()`.

Find:

```r
test_that("the barrnap branch of pipeline_r.nf translates the kingdom", {
  nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))
  start <- grep("^process annotate_tools \\{", nf)
  expect_length(start, 1)
  end <- start - 1 + which(nf[start:length(nf)] == "}")[1]
  body <- nf[start:end]

  bar <- grep("barrnap --kingdom", body)
  expect_length(bar, 1)
  expect_true(any(grepl("case \"\\\\\\$KINGDOM\" in euk\\|all\\) KINGDOM=fun", body[seq_len(bar)])))
  py <- grep("pybarrnap -k", body)
  expect_length(py, 1)
  expect_true(any(grepl("case \"\\\\\\$KINGDOM\" in fun\\) KINGDOM=euk", body[seq_len(py)])))
})
```

Replace with:

```r
test_that("rp_rrna_tools.sh translates the kingdom as rrna_kingdom_for does", {
  skip_on_os("windows")
  sh <- system.file("nextflow", "bin", "rp_rrna_tools.sh", package = "RiboPilot")
  kingdom_passed <- function(tool, kingdom) {
    dir <- withr::local_tempdir()
    mk_fake_tool(dir, tool, out = GFF_18S)
    mk_fake_tool(dir, "ITSx")
    withr::local_path(dir)
    mk_fasta(dir)
    status <- withr::with_dir(dir, system2("bash",
      shQuote(c(sh, "unit.fasta", tool, kingdom, "", "", "", "1")), stdout = FALSE, stderr = FALSE))
    expect_equal(status, 0)
    args <- strsplit(fake_args(dir, tool), " ", fixed = TRUE)[[1]]
    args[which(args %in% c("--kingdom", "-k")) + 1]
  }
  for (tool in c("barrnap", "pybarrnap")) {
    for (k in c("euk", "all", "fun", "bac", "arc", "")) {
      expect_equal(kingdom_passed(tool, k), RiboPilot:::rrna_kingdom_for(k, tool),
        info = paste(tool, k))
    }
  }
})
```

In `tests/testthat/test-nextflow.R` ("user rRNA caller flags reach both the per-unit and the scan commands"), the per-unit barrnap/pybarrnap command lines leave pipeline_r.nf on purpose; it now passes `"$RRNAO"` to the script, whose argv test-rrna-tools.R checks.

Find:

```r
  expect_true(any(grepl("--quiet \\$RRNAO unit.fasta", nf, fixed = TRUE)))
  expect_true(any(grepl("-q \\$RRNAO unit.fasta", nf, fixed = TRUE)))
```

Replace with:

```r
  expect_true(any(grepl('"\\$RRNAO" ${task.cpus}', nf, fixed = TRUE)))
```

- [ ] **Step 2: Run them to see them fail**

Run:
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rrna-tools.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rrna-kingdom.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nextflow.R")'`

Expected: FAIL. test-rrna-tools.R: `file.access(rrna_sh(), 1)` is -1, every script run exits 127 with `bash: : No such file or directory` (the script does not exist yet), and the annotate_tools test finds no `rp_rrna_tools.sh` line. test-rrna-kingdom.R: the new kingdom test fails (status 127, kingdom `character(0)`). test-nextflow.R: one failure, the `"\$RRNAO" ${task.cpus}` expectation. Without `NOT_CRAN=true` the Nextflow test in test-rrna-tools.R skips ("On CRAN"); with `NOT_CRAN=true` it passes already, because it pins today's annotate_tools failure and retry behavior, which this task must keep.

- [ ] **Step 3: Create the script**

Create `inst/nextflow/bin/rp_rrna_tools.sh`:

```bash
#!/bin/bash
# rRNA anchors (barrnap or pybarrnap), the tuned ITSx run, and the stock ITSx
# check on one contig. Writes unit.barrnap.gff, unit.itsx.*, unit.itsx_check.*
# (the stock run's stderr in unit.itsx_check.err), and unit.tools_status (each
# tool's exit status). Runs all three tools, then exits with the barrnap or
# pybarrnap status; ITSx never changes it.
set -euo pipefail
if [ $# -ne 7 ]; then
    echo "usage: rp_rrna_tools.sh <fasta> <tool> <kingdom> <itsx_taxon> <itsx_opts> <rrna_opts> <cpus>" >&2
    exit 2
fi
fa=$1; tool=${2:-pybarrnap}; kingdom=${3:-euk}; taxon=${4:-all}
# $itsxo and $rrnao stay unquoted below: each holds several flags
itsxo=$5; rrnao=$6; cpus=$7
b=0; i=0; c=0

# Same kingdom translation as R's rrna_kingdom_for(): the Perl barrnap fork
# calls its eukaryote Rfam models "fun" and has no euk/all; pybarrnap has no
# "fun". Passing the stored kingdom through untranslated makes the tool exit
# non-zero and the unit gets no anchors at all.
if [ "$tool" = "barrnap" ]; then
    case "$kingdom" in euk|all) kingdom=fun ;; esac
    barrnap --kingdom "$kingdom" --quiet $rrnao "$fa" > unit.barrnap.gff || b=$?
else
    case "$kingdom" in fun) kingdom=euk ;; esac
    pybarrnap -k "$kingdom" -t "$cpus" -q $rrnao "$fa" > unit.barrnap.gff || b=$?
fi

ITSx -i "$fa" -o unit.itsx --save_regions none --graphical F \
    --cpu "$cpus" -t "$taxon" --preserve T $itsxo || i=$?
[ -f unit.itsx.positions.txt ] || : > unit.itsx.positions.txt

# second, stock run (defaults, no organism group): what NCBI's Submission
# Portal will see. Its output becomes itsx_warnings.
ITSx -i "$fa" -o unit.itsx_check --preserve T --detailed_results T \
    --save_regions none --graphical F --cpu "$cpus" 2> unit.itsx_check.err || c=$?
[ -f unit.itsx_check.positions.txt ] || : > unit.itsx_check.positions.txt

printf 'barrnap=%s\nitsx=%s\nitsx_check=%s\n' "$b" "$i" "$c" > unit.tools_status
exit "$b"
```

Then make it executable (Nextflow calls bin scripts by name, and the Write tool creates 0644 files):

Run: `chmod +x inst/nextflow/bin/rp_rrna_tools.sh`

- [ ] **Step 4: Call the script from annotate_tools**

In `inst/nextflow/modules/pipeline_r.nf`, process `annotate_tools`. Leave everything above this block (the fasta write and the options read) as it is.

Find:

```
    # Same kingdom translation as R's rrna_kingdom_for(): the Perl barrnap fork
    # calls its eukaryote Rfam models "fun" and has no euk/all; pybarrnap has no
    # "fun". Passing the stored kingdom through untranslated makes the tool exit
    # non-zero and the unit gets no anchors at all.
    if [ "\$TOOL" = "barrnap" ]; then
        case "\$KINGDOM" in euk|all) KINGDOM=fun ;; esac
        barrnap --kingdom "\$KINGDOM" --quiet \$RRNAO unit.fasta > unit.barrnap.gff
    else
        case "\$KINGDOM" in fun) KINGDOM=euk ;; esac
        pybarrnap -k "\$KINGDOM" -t ${task.cpus} -q \$RRNAO unit.fasta > unit.barrnap.gff
    fi

    ITSx -i unit.fasta -o unit.itsx --save_regions none --graphical F \
        --cpu ${task.cpus} -t "\$ITSXT" --preserve T \$ITSXO || true
    [ -f unit.itsx.positions.txt ] || : > unit.itsx.positions.txt

    # second, stock run (defaults, no organism group): what NCBI's Submission
    # Portal will see. Reconcile turns its output into annotate.itsx_warnings.
    ITSx -i unit.fasta -o unit.itsx_check --preserve T --detailed_results T \
        --save_regions none --graphical F --cpu ${task.cpus} || true
    [ -f unit.itsx_check.positions.txt ] || : > unit.itsx_check.positions.txt
```

Replace with:

```
    # barrnap/pybarrnap, tuned ITSx, and the stock ITSx check. Called plainly so a
    # barrnap failure fails the task and an out-of-memory exit is retried.
    rp_rrna_tools.sh unit.fasta "\$TOOL" "\$KINGDOM" "\$ITSXT" "\$ITSXO" "\$RRNAO" ${task.cpus}
```

- [ ] **Step 5: Point the two kingdom comments at the script**

In `R/annotate.R` (roxygen of `rrna_kingdom_for`; `@noRd`, so no Rd change):

Find:

```r
#' place that mapping lives on the R side; the shell branch of
#' inst/nextflow/modules/pipeline_r.nf mirrors the same rule.
```

Replace with:

```r
#' place that mapping lives on the R side; inst/nextflow/bin/rp_rrna_tools.sh
#' mirrors the same rule.
```

In `inst/nextflow/modules/find_rrna.nf` (process `find_rrna_scan`):

Find:

```
    # Same kingdom translation as annotate_tools / R's rrna_kingdom_for(): the
```

Replace with:

```
    # Same kingdom translation as bin/rp_rrna_tools.sh / R's rrna_kingdom_for(): the
```

- [ ] **Step 6: Run the tests to see them pass**

Run:
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rrna-tools.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-rrna-kingdom.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nextflow.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-bin-contract.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-main-run.R")'`

Expected: PASS. test-rrna-tools.R 39 (the Nextflow test skips when Nextflow or sqlite3 is missing); test-rrna-kingdom.R 64; test-nextflow.R 298; test-nf-bin-contract.R 162; test-nf-main-run.R 33. test-nf-main-run.R runs the real main.nf with the fake tools (Nextflow puts `inst/nextflow/bin` on PATH), so it proves `annotate_tools` reaches the script and `reconcile` still gets its 6-tuple; it takes a few minutes.

- [ ] **Step 7: Commit**

```bash
git add --chmod=+x inst/nextflow/bin/rp_rrna_tools.sh
git add inst/nextflow/modules/pipeline_r.nf inst/nextflow/modules/find_rrna.nf R/annotate.R \
  tests/testthat/test-rrna-tools.R tests/testthat/test-rrna-kingdom.R tests/testthat/test-nextflow.R
git ls-files -s inst/nextflow/bin/rp_rrna_tools.sh
git commit -m "Move rRNA tool calls into rp_rrna_tools.sh"
```

Expected: the `git ls-files -s` line starts with `100755`, like every other file in `inst/nextflow/bin`.


### Task 10: recheck_unit

**Files:**
- Create: `R/recheck.R`
- Create: `tests/testthat/test-recheck.R`
- Test: `tests/testthat/test-recheck.R`, `tests/testthat/test-no-duplicate-definitions.R`, `tests/testthat/test-docs-site.R`

**Interfaces:**
- Consumes:
  - Task 1: `new_db()` creates `review_units` (with `rev INTEGER NOT NULL`, `check_switch`, `check_note`, `checked_at`, `n_features`, `missing`, `chimeric`, `itsx_*`, `blast_*`, `time_stamp`), the 7 mirror tables (`review_annotations` with `partial_start_user`/`partial_stop_user`, `review_blast_hits`, `review_blast_screen`, `review_screen`, `review_chimera`, `review_chimera_windows`, `review_ref_alignment`), and the 5 views.
  - Task 2: `with_savepoint(con, code, name = "rp_review")`; `review_counts(con, id, path, scaffold)` -> `list(n_features, missing)`; `review_drop(con, id, path, scaffold)` (tests).
  - Task 3: `review_map_hits(con, id, path, scaffold)` -> data.frame in `blast_hits` column order (original HSPs overlapping the window, clipped, shifted by `trim_start - 1`; clipped HSPs get `aln_len` = clipped span and `sstart`/`send` NA; `qcovs` = floor(100 * union of the subject's clipped spans / copy length + 0.5)); `review_coverage()` (through `chimera_unit`).
  - Task 4 (tests only): `review_trim(con, id, path, scaffold, from, to)` (copy coordinates; makes the copy, invalidates it); `review_set_state(con, units, state)` (1 = invalidate, 0 = `check_switch` 0 with `rev` unchanged).
  - Task 5: `itsx_store_not_run(con, id, path, scaffold, reason, review = FALSE)`; `itsx_check_one()` stores `not_run(ITSx did not finish)` unless `<out_prefix>.summary.txt` contains "ITSx run finished".
  - Task 6: `curate_one(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE)`; `compute_ref_alignment(con, id, path = 1L, scaffold = 1L, review = FALSE)`.
  - Task 7: `screen_one(con, id, path = 1L, scaffold = 1L, force = FALSE, review = FALSE)`; `parse_blast_candidates(con, id, path = 1L, scaffold = 1L, hits_file, taxdump_dir = NULL, review = FALSE)`; `derive_blast_products(con, id, path = 1L, scaffold = 1L, ts = as.integer(Sys.time()), review = FALSE)` (review mode clears `review_units.blast_*` when `review_blast_hits` is empty; `blast_source` from the original `blast_searched.source`).
  - Task 8: `chimera_unit(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL, label_db = NULL, threads = 1L, workdir = tempfile("chimera_"), review = FALSE)`; `itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix, review = FALSE)`.
  - Task 9 (file names only; the tests write these files by hand): `unit.tools_status` with lines `barrnap=<exit>`, `itsx=<exit>`, `itsx_check=<exit>`; `unit.barrnap.gff`; `unit.itsx.positions.txt`; `unit.itsx_check.*`.
- Produces (all in `R/recheck.R`, `@noRd`; bin scripts call them as `RiboPilot:::`):
  - `RECHECK_BLAST`: `list(outfmt = "6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid", task = "megablast", retry_task = "dc-megablast", max_hsps = 5L, retry_below_qcovs = 25, default_max_target_seqs = 50L)`.
  - `recheck_blast_local(seq, dbarg, max_target_seqs, extra_opts, workdir, threads = 1L, blastn = "blastn")` -> path of `<workdir>/hits.txt`, or NULL when megablast failed. Query header `>query`. Never remote.
  - `recheck_hits(con, id, path, scaffold, mode, dbarg = NULL, taxdump_dir = NULL, threads = 1L, workdir = tempfile("recheck_"))` -> `"local"` or `"map"`, whichever was used. `dbarg` NULL, `""`, or NA means no panel.
  - `recheck_model_coverage(con, id, path, scaffold, gff)` -> invisible numeric of the coverages written; no-op when `gff` is NULL.
  - `recheck_unit(con, id, path, scaffold, rev, mode, tools_dir, dbarg = NULL, taxdump_dir = NULL, label_db = NULL, threads = 1L, workdir = tempfile("recheck_"))` -> invisible `"done"`, `"stale"`, or `"failed"`. `rev` is the text `query_recheck` emitted. Workdirs used: `<workdir>/blast`, `<workdir>/chimera`.
  - `recheck_release_running(con)` -> integer count of copies moved from 4 to 1 (`check_note = 'recheck did not finish'`, `time_stamp` now).

This task adds the R side of the recheck (SPEC S7.3): `recheck_unit()` reruns every offline check on one contig's review copy in the spec's order (hits, screen, counts, barrnap model coverage, ITSx chimera flag, stock ITSx check, reference alignment, curation, chimera), through the shared check functions in review mode, so the main pipeline and the recheck score with the same code. Hits come from a fresh local megablast with blast.nf's settings and dc-megablast retry, or from the original hits mapped onto the copy's window; it never contacts NCBI (D1), and a missing panel or a failed search falls back to the mapped hits (Review Focus 4). A failed tool leaves its value as copied, and a crashed or throwing ITSx check stores `not_run(...)` instead of failing the contig (S6.3, Review Focus 3). The final write sets `check_switch = 2` only while the copy still has the run's `rev` and is queued or running, so an edit, On hold, or Restore during the run discards the result (S5.4); an error marks the copy failed under the same guard. `recheck_release_running()` returns leftover running copies to not run after an in-app run (S5.4); Task 13 calls it. Nothing calls `recheck_unit()` yet; Task 11's `rp_recheck.R` will.

Notes for this task:
- Never wrap the steps of `recheck_unit()` in a transaction or `with_savepoint()`: `parse_blast_candidates`, `derive_blast_products`, `screen_one`, `curate_one`, and `chimera_unit` each open their own `DBI::dbWithTransaction()`, which fails inside an open transaction.
- `rev` is a millisecond timestamp. Carry it as text and compare it only in SQL (`rev = CAST(? AS INTEGER)`).
- Always pass `review = TRUE` by name (`curate_one`'s 5th positional argument is `ruleset`).
- The tests mock check functions with `local_mocked_bindings(..., .package = "RiboPilot")`, which works under `devtools::load_all(".")`.
- Stage only the two new files. The user's uncommitted edits in `vignettes/Troubleshooting.Rmd` and `vignettes/Your-Own-Project.Rmd` stay unstaged.

- [ ] **Step 1: Write the failing tests for the BLAST helpers and model coverage**

Create `tests/testthat/test-recheck.R` with this content. The `rc_*` helpers at the top serve both test cycles of this task. The fixture contig is the 6310 bp `fixtures/chimera/clean.fa` scaffold (18S 351-2150, ITS1, 5.8S 2651-2807, ITS2, 28S 3991-6310), annotated through the real `reconcile_one()`. Trims cut only the 3' end, except in the mapping tests, which run no ITSx check.

```r
# recheck_unit and its helpers: the per-contig recheck of a review copy
# (SPEC S7.3, S5.4).

RC_SEQ <- paste(readLines(test_path("fixtures", "chimera", "clean.fa"))[-1], collapse = "")
RC_UNIT <- "ID = 'S1' AND path = 1 AND scaffold = 1"

# GFF lines as (py)barrnap writes them; cov < 100 marks a partial hit.
rc_gff <- function(file, element, start, stop, cov = 100, strand = "+") {
  nm <- c("18S" = "18S_rRNA", "5.8S" = "5_8S_rRNA", "28S" = "28S_rRNA")[element]
  note <- ifelse(cov < 100, sprintf(" (partial);Note=aligned only %s percent of the model", cov), "")
  writeLines(paste("S1.1.1", "barrnap:0.9", "rRNA", start, stop, "0", strand, ".",
    paste0("Name=", nm, ";product=", element, " ribosomal RNA", note), sep = "\t"), file)
}

# An ITSx positions line for a contig of `len` bp whose 28S ends at `lsu_end`.
rc_positions <- function(file, len, lsu_end, chimeric = FALSE) {
  writeLines(paste0("S1.1.1\t", len, " bp.\tSSU: 351-2150\tITS1: 2151-2650\t5.8S: 2651-2807\t",
    "ITS2: 2808-3990\tLSU: 3991-", lsu_end, "\t", if (chimeric) "Chimeric! " else ""), file)
}

# A 6310 bp contig S1.1.1 annotated from barrnap and ITSx output, done, named
# by the seed database, with flat coverage and, unless hits = FALSE, two
# original subjects: A1.1 (HSPs 1-3000 and 5801-6310) and B2.1 (6101-6310).
rc_project <- function(hits = TRUE, env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  db <- file.path(d, ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE, seeds_db = fa, labels_db = fa)
  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con), envir = env)
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore,
      edited, blast_accession, blast_species, blast_pident, blast_source, time_stamp)
    VALUES ('S1', 1, 1, 'linear', ?, ?, 0, 0, 'SEED1', 'Danio rerio', 99.4, 'seeds', 1)",
    params = list(nchar(RC_SEQ), RC_SEQ))
  rc_gff(file.path(d, "o.gff"), c("18S", "5.8S", "28S"), c(351, 2651, 3991), c(2150, 2807, 6310))
  rc_positions(file.path(d, "o.pos"), 6310, 6310)
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, file.path(d, "o.gff"), file.path(d, "o.pos"))
  DBI::dbExecute(con, paste("UPDATE annotate SET itsx_warnings = '' WHERE", RC_UNIT))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1', 1, 1, ?, ?)",
    params = list(paste(rep(30, 6310), collapse = ","), paste(rep(0, 6310), collapse = ",")))
  if (hits) {
    h <- data.frame(ID = "S1", path = 1L, scaffold = 1L, hit_rank = c(1L, 1L, 2L), hsp = c(1L, 2L, 1L),
      accession = c("A1.1", "A1.1", "B2.1"),
      title = c(rep("Danio rerio 18S ribosomal RNA gene, partial sequence", 2),
        "Danio aesculapii 28S ribosomal RNA gene, partial sequence"),
      species = c("Danio rerio", "Danio rerio", "Danio aesculapii"), staxid = c(7955L, 7955L, 1841699L),
      pident = c(99.5, 99.5, 98), qcovs = c(56, 56, 3), evalue = 0, bitscore = c(5000, 900, 380),
      aln_len = c(3000L, 510L, 210L), qstart = c(1L, 5801L, 6101L), qend = c(3000L, 6310L, 6310L),
      sstart = c(1L, 5801L, 1L), send = c(3000L, 6310L, 210L), slen = 7000L, is_nuclear_rrna = 1L,
      hit_class = "rdna", time_stamp = 1L)
    DBI::dbAppendTable(con, "blast_hits", h)
  }
  con
}

# Output of rp_rrna_tools.sh for a copy trimmed to 1..6000: barrnap finds the
# cut 28S at 70 percent of its model, and the tuned ITSx run calls it chimeric.
rc_tools <- function(status = c(barrnap = 0, itsx = 0, itsx_check = 0), finished = TRUE,
                     env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  if (!is.null(status)) writeLines(paste0(names(status), "=", status), file.path(d, "unit.tools_status"))
  rc_gff(file.path(d, "unit.barrnap.gff"), c("18S", "5.8S", "28S"), c(351, 2651, 3991),
    c(2150, 2807, 6000), cov = c(100, 100, 70))
  rc_positions(file.path(d, "unit.itsx.positions.txt"), 6000, 6000, chimeric = TRUE)
  rc_positions(file.path(d, "unit.itsx_check.positions.txt"), 6000, 6000)
  if (finished) writeLines("ITSx run finished at today.", file.path(d, "unit.itsx_check.summary.txt"))
  d
}

rc_rev <- function(con) {
  DBI::dbGetQuery(con, paste("SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE", RC_UNIT))$rev
}

rc_copy <- function(con) {
  DBI::dbGetQuery(con, paste("SELECT check_switch, check_note, checked_at, n_features, missing, chimeric,
    n_flags, curate_flags, itsx_warnings, blast_accession, blast_qcovs, blast_source, time_stamp
    FROM review_units WHERE", RC_UNIT))
}

rc_count <- function(con, tbl) {
  DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM", tbl, "WHERE", RC_UNIT))$n
}

rc_cov <- function(con, tbl = "review_annotations") {
  a <- DBI::dbGetQuery(con, paste("SELECT element, model_coverage FROM", tbl, "WHERE", RC_UNIT))
  stats::setNames(a$model_coverage, a$element)
}

# Every original table the recheck could touch, whole.
rc_originals <- function(con) {
  tb <- c("assemblies", "annotate", "annotations", "coverage", "screen", "chimera", "chimera_windows",
    "blast_hits", "blast_screen", "blast_ref_candidates", "blast_ref_alignment")
  stats::setNames(lapply(tb, function(t) DBI::dbReadTable(con, t)), tb)
}

# A blastn stand-in: logs each call's arguments to RC_BLAST_LOG (one per line,
# then "--"), prints the file RC_MEGA or RC_DC for -task megablast or
# dc-megablast, and exits 1 when that variable is "fail".
rc_blastn <- function(env = parent.frame()) {
  skip_on_os("windows")
  d <- withr::local_tempdir(.local_envir = env)
  f <- file.path(d, "blastn")
  writeLines(c("#!/usr/bin/env bash",
    'prev=""; task=""',
    'for a in "$@"; do printf "%s\\n" "$a" >> "$RC_BLAST_LOG"; [ "$prev" = "-task" ] && task=$a; prev=$a; done',
    'echo "--" >> "$RC_BLAST_LOG"',
    'if [ "$task" = "dc-megablast" ]; then out=$RC_DC; else out=$RC_MEGA; fi',
    '[ "$out" = "fail" ] && exit 1',
    '[ -n "$out" ] && cat "$out"',
    "exit 0"), f)
  Sys.chmod(f, "0755")
  withr::local_envvar(RC_BLAST_LOG = file.path(d, "calls.log"), RC_MEGA = "", RC_DC = "", .local_envir = env)
  f
}

rc_calls <- function() {
  log <- Sys.getenv("RC_BLAST_LOG")
  if (!file.exists(log)) return(list())
  x <- readLines(log)
  unname(split(x[x != "--"], cumsum(x == "--")[x != "--"]))
}

rc_hit_lines <- function(file, qcovs, accession = "P1.1") {
  writeLines(paste(accession, "Danio rerio 18S ribosomal RNA gene, partial sequence", 99.1, qcovs,
    "0.0", 7955, 4000, 2000, 1, 2000, 1, 2000, 7000, "query", sep = "\t"), file)
  file
}

test_that("the recheck's BLAST settings match the local search in blast.nf", {
  nf <- paste(readLines(system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")),
    collapse = "\n")
  b <- RiboPilot:::RECHECK_BLAST
  for (s in c(sprintf('OUTFMT="%s"', b$outfmt), sprintf("-max_hsps %d -task %s", b$max_hsps, b$task),
              paste("-task", b$retry_task), sprintf("-lt %d ]", b$retry_below_qcovs),
              sprintf("COALESCE(b.max_target_seqs,%d)", b$default_max_target_seqs))) {
    expect_true(grepl(s, nf, fixed = TRUE), info = s)
  }
})

test_that("recheck_blast_local runs megablast with the sample's options and keeps a good result", {
  blastn <- rc_blastn()
  wd <- withr::local_tempdir()
  Sys.setenv(RC_MEGA = rc_hit_lines(file.path(wd, "mega.txt"), 80))
  out <- RiboPilot:::recheck_blast_local("ACGTACGT", "/db/panel /db/decoys", 7L,
    " -evalue 1e-20  -perc_identity 90 ", file.path(wd, "run"), threads = 2L, blastn = blastn)
  expect_equal(out, file.path(wd, "run", "hits.txt"))
  expect_equal(readLines(out), readLines(file.path(wd, "mega.txt")))
  expect_equal(readLines(file.path(wd, "run", "query.fa")), c(">query", "ACGTACGT"))
  calls <- rc_calls()
  expect_length(calls, 1)
  expect_equal(calls[[1]], c("-db", "/db/panel /db/decoys", "-query", file.path(wd, "run", "query.fa"),
    "-outfmt", RiboPilot:::RECHECK_BLAST$outfmt, "-max_target_seqs", "7", "-max_hsps", "5",
    "-task", "megablast", "-num_threads", "2", "-evalue", "1e-20", "-perc_identity", "90"))
})

test_that("recheck_blast_local retries with dc-megablast as blast.nf does", {
  blastn <- rc_blastn()
  wd <- withr::local_tempdir()
  mega <- rc_hit_lines(file.path(wd, "mega.txt"), 24, "MEGA.1")
  dc <- rc_hit_lines(file.path(wd, "dc.txt"), 90, "DC.1")
  empty <- file.path(wd, "empty.txt")
  file.create(empty)
  run <- function(m, d, sub) {
    Sys.setenv(RC_MEGA = m, RC_DC = d)
    unlink(Sys.getenv("RC_BLAST_LOG"))
    out <- RiboPilot:::recheck_blast_local("ACGT", "panel", 50L, "", file.path(wd, sub), blastn = blastn)
    list(lines = readLines(out), tasks = vapply(rc_calls(), function(a) a[which(a == "-task") + 1], ""))
  }
  # best qcovs under 25: the retry's hits replace megablast's
  r <- run(mega, dc, "a")
  expect_equal(r$tasks, c("megablast", "dc-megablast"))
  expect_equal(r$lines, readLines(dc))
  # no hits at all: the retry runs too
  r <- run(empty, dc, "b")
  expect_equal(r$tasks, c("megablast", "dc-megablast"))
  expect_equal(r$lines, readLines(dc))
  # a retry with no output keeps the sentinel, a failed retry keeps megablast's hits
  expect_equal(run(empty, empty, "c")$lines, "NO_SIGNIFICANT_HITS")
  expect_equal(run(mega, "fail", "d")$lines, readLines(mega))
  # 25 or more: no retry
  r <- run(rc_hit_lines(file.path(wd, "m25.txt"), 25), dc, "e")
  expect_equal(r$tasks, "megablast")
  # megablast itself failed: no hits file
  Sys.setenv(RC_MEGA = "fail")
  expect_null(RiboPilot:::recheck_blast_local("ACGT", "panel", 50L, "", file.path(wd, "f"), blastn = blastn))
})

test_that("map mode moves the original hits onto the copy's window", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 101L, 6000L)
  before <- rc_originals(con)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  h <- DBI::dbGetQuery(con, paste("SELECT accession, qstart, qend, aln_len, sstart, qcovs
    FROM review_blast_hits WHERE", RC_UNIT, "ORDER BY hit_rank, hsp"))
  expect_equal(h$accession, c("A1.1", "A1.1"))
  expect_equal(h$qstart, c(1L, 5701L))
  expect_equal(h$qend, c(2900L, 5900L))
  expect_equal(h$aln_len, c(2900L, 200L))
  expect_equal(h$sstart, c(NA_integer_, NA_integer_))
  expect_equal(h$qcovs, c(53, 53))
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "A1.1")
  expect_equal(u$blast_qcovs, 53)
  expect_equal(rc_count(con, "review_blast_screen"), 1)
  expect_equal(rc_originals(con), before)
})

test_that("map mode keeps the copied summary only when the original has no hits", {
  con <- rc_project(hits = FALSE)
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  # rows a discarded recheck left behind
  DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank, hsp, accession)
    VALUES ('S1', 1, 1, 1, 1, 'OLD.1')")
  DBI::dbExecute(con, "INSERT INTO review_blast_screen (ID, path, scaffold, accession) VALUES ('S1', 1, 1, 'OLD.1')")
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  expect_equal(rc_count(con, "review_blast_hits"), 0)
  expect_equal(rc_count(con, "review_blast_screen"), 0)
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "SEED1")
  expect_equal(u$blast_source, "seeds")

  # the original has hits, but none inside the window: the summary goes
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 3001L, 5800L)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  expect_equal(rc_count(con, "review_blast_hits"), 0)
  expect_true(is.na(rc_copy(con)$blast_accession))
})

test_that("local mode searches the panel with the sample's options and never NCBI", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  DBI::dbExecute(con, "UPDATE blast_opts SET extra_opts = '-evalue 1e-5', max_target_seqs = 9
    WHERE blast_opts = 'default'")
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('S1', 1, 1, 'x', 'local')")
  blastn <- rc_blastn()
  withr::local_envvar(PATH = paste(dirname(blastn), Sys.getenv("PATH"), sep = .Platform$path.sep))
  wd <- withr::local_tempdir()
  Sys.setenv(RC_MEGA = rc_hit_lines(file.path(wd, "mega.txt"), 95))
  before <- rc_originals(con)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = "panel"), "local")
  a <- rc_calls()[[1]]
  expect_false("-remote" %in% a)
  expect_equal(a[which(a == "-max_target_seqs") + 1], "9")
  expect_equal(utils::tail(a, 2), c("-evalue", "1e-5"))
  expect_equal(DBI::dbGetQuery(con, paste("SELECT accession FROM review_blast_hits WHERE", RC_UNIT))$accession,
    "P1.1")
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "P1.1")
  expect_equal(u$blast_source, "local_rrna")
  expect_equal(rc_originals(con), before)

  # a failed search falls back to the mapped original hits
  Sys.setenv(RC_MEGA = "fail")
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = "panel"), "map")
  expect_equal(rc_copy(con)$blast_accession, "A1.1")
})

test_that("local mode without a panel maps the hits and never calls blastn", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  blastn <- rc_blastn()
  withr::local_envvar(PATH = paste(dirname(blastn), Sys.getenv("PATH"), sep = .Platform$path.sep))
  for (dbarg in list(NULL, "", NA_character_)) {
    expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = dbarg), "map")
  }
  expect_length(rc_calls(), 0)
  expect_equal(rc_copy(con)$blast_accession, "A1.1")
})

test_that("recheck_model_coverage refreshes only barrnap features, by element, strand, and overlap", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  DBI::dbExecute(con, paste("UPDATE review_annotations SET tool = 'barrnap+ITSx+ITSx', model_coverage = 55
    WHERE element = '5.8S' AND", RC_UNIT))
  DBI::dbExecute(con, paste("UPDATE review_annotations SET model_coverage = 42 WHERE element = 'ITS1' AND", RC_UNIT))
  gff <- withr::local_tempfile()
  writeLines(c(
    paste("S1.1.1", "barrnap:0.9", "rRNA", 1, 400, "0", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA (partial);Note=aligned only 20 percent of the model", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 351, 2150, "0", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA (partial);Note=aligned only 91 percent of the model", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 2651, 2807, "0", "+", ".",
      "Name=5_8S_rRNA;product=5.8S ribosomal RNA", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 3991, 6000, "0", "-", ".",
      "Name=28S_rRNA;product=28S ribosomal RNA", sep = "\t")), gff)
  orig <- rc_cov(con, "annotations")
  RiboPilot:::recheck_model_coverage(con, "S1", 1L, 1L, NULL)
  expect_equal(rc_cov(con)[c("18S", "28S")], orig[c("18S", "28S")])
  RiboPilot:::recheck_model_coverage(con, "S1", 1L, 1L, gff)
  cv <- rc_cov(con)
  expect_equal(unname(cv[c("18S", "5.8S", "28S", "ITS1")]), c(91, 100, 0, 42))
  expect_true(is.na(cv[["ITS2"]]))
  expect_equal(rc_cov(con, "annotations"), orig)
})
```

What the numbers mean: `review_trim(..., 101L, 6000L)` keeps original bases 101-6000 (a 5900 bp copy). A1.1's HSPs become 1-2900 and 5701-5900, both clipped (so `sstart` is NA), and its `qcovs` is floor(100 * 3100 / 5900 + 0.5) = 53; B2.1 (6101-6310) falls outside. In the model coverage test the 18S takes the hit that overlaps it most (91, not 20), the 5.8S (tool `barrnap+ITSx+ITSx`) matches a full hit (100), the 28S hit is on the other strand so the feature gets 0, and the ITSx-sourced ITS1 and ITS2 keep their values.

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-recheck.R")'`
Expected: FAIL, `[ FAIL 8 | WARN 0 | SKIP 0 | PASS 0 ]`, every test erroring with `object 'RECHECK_BLAST' not found`, `object 'recheck_blast_local' not found`, `object 'recheck_hits' not found`, or `object 'recheck_model_coverage' not found`.

- [ ] **Step 3: Create R/recheck.R with the BLAST helpers and model coverage**

Create `R/recheck.R`:

```r
# The recheck: reruns every offline check on one contig's review copy
# (dev/specs/2026-10-08-contig-recheck-design.md, section 7.3).

#' blastn settings of the local search in inst/nextflow/modules/blast.nf
#'
#' The recheck repeats that search in R. test-recheck.R checks these values
#' against blast.nf, so the two cannot drift.
#' @noRd
RECHECK_BLAST <- list(
  outfmt = "6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid",
  task = "megablast",
  retry_task = "dc-megablast",
  max_hsps = 5L,
  retry_below_qcovs = 25,
  default_max_target_seqs = 50L
)

#' Local megablast of one sequence, with blast.nf's dc-megablast retry
#'
#' Retries with dc-megablast when megablast finds nothing or its best query
#' coverage is under `RECHECK_BLAST$retry_below_qcovs`, and keeps the retry only
#' when it exits 0 with output. Never searches NCBI.
#' @param seq Query sequence.
#' @param dbarg blastn `-db` value (rp_dbarg.R's panel list), passed as one argument.
#' @param max_target_seqs,extra_opts The sample's BLAST options; `extra_opts`
#'   is split on whitespace, as blast.nf's unquoted `${XO}` is.
#' @param workdir Folder for query.fa and the hits files.
#' @return Path of the hits file, or NULL when megablast failed.
#' @noRd
recheck_blast_local <- function(seq, dbarg, max_target_seqs, extra_opts, workdir,
                                threads = 1L, blastn = "blastn") {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  query <- file.path(workdir, "query.fa")
  writeLines(c(">query", seq), query)
  xo <- strsplit(trimws(extra_opts %|NA|% ""), "[[:space:]]+")[[1]]
  run <- function(task, out) {
    st <- suppressWarnings(system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query),
      "-outfmt", shQuote(RECHECK_BLAST$outfmt), "-max_target_seqs", max_target_seqs,
      "-max_hsps", RECHECK_BLAST$max_hsps, "-task", task, "-num_threads", threads,
      shQuote(xo[nzchar(xo)])), stdout = out, stderr = FALSE))
    identical(as.integer(st), 0L)
  }
  hits <- file.path(workdir, "hits.txt")
  if (!run(RECHECK_BLAST$task, hits)) return(NULL)
  if (file.size(hits) == 0) writeLines("NO_SIGNIFICANT_HITS", hits)
  h <- parse_blast_hits_file(hits)
  if (nrow(h) == 0 || max(h$qcovs, 0, na.rm = TRUE) < RECHECK_BLAST$retry_below_qcovs) {
    dc <- file.path(workdir, "hits_dc.txt")
    if (run(RECHECK_BLAST$retry_task, dc) && isTRUE(file.size(dc) > 0)) file.rename(dc, hits)
  }
  hits
}

#' Fresh BLAST hits for a review copy (SPEC S7.3 step 2)
#'
#' `local` searches the panel again with the sample's options. Without a panel
#' (`dbarg` NULL or "") or when that search fails, it maps the original hits
#' instead, as `map` does. A copy whose original has no hits at all keeps the
#' BLAST summary it copied (a seed-database name, for example). Writes
#' review_blast_hits, review_blast_screen, and review_units.blast_*.
#' @return "local" or "map", whichever was used.
#' @noRd
recheck_hits <- function(con, id, path, scaffold, mode, dbarg = NULL, taxdump_dir = NULL,
                         threads = 1L, workdir = tempfile("recheck_")) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  unit <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  p <- list(id, path, scaffold)
  if (identical(mode, "local") && isTRUE(nzchar(dbarg, keepNA = TRUE))) {
    o <- DBI::dbGetQuery(con, "SELECT b.extra_opts, b.max_target_seqs FROM assemble a
      JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = ?", params = list(id))
    seq <- DBI::dbGetQuery(con, paste("SELECT sequence FROM review_units", unit), params = p)$sequence
    hits <- recheck_blast_local(seq, dbarg,
      o$max_target_seqs[1] %|NA|% RECHECK_BLAST$default_max_target_seqs,
      o$extra_opts[1] %|NA|% "", workdir, threads)
    if (!is.null(hits)) {
      parse_blast_candidates(con, id, path, scaffold, hits, taxdump_dir = taxdump_dir, review = TRUE)
      return("local")
    }
  }
  if (DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM blast_hits", unit), params = p)$n == 0) {
    with_savepoint(con, for (t in c("review_blast_hits", "review_blast_screen")) {
      DBI::dbExecute(con, paste("DELETE FROM", t, unit), params = p)
    })
    return("map")
  }
  m <- review_map_hits(con, id, path, scaffold)
  with_savepoint(con, {
    DBI::dbExecute(con, paste("DELETE FROM review_blast_hits", unit), params = p)
    DBI::dbAppendTable(con, "review_blast_hits", m)
  })
  derive_blast_products(con, id, path, scaffold, review = TRUE)
  "map"
}

#' Refresh the model coverage of a copy's barrnap features (SPEC S7.3 step 5)
#'
#' Each barrnap-sourced feature (tool "barrnap", "barrnap+ITSx", ...) takes the
#' coverage of the barrnap hit with its element and strand that overlaps it
#' most; a feature no hit overlaps gets 0. Other features are left alone.
#' @param gff barrnap GFF of the copy's sequence, or NULL to change nothing.
#' @noRd
recheck_model_coverage <- function(con, id, path, scaffold, gff) {
  if (is.null(gff)) return(invisible())
  p <- list(id, as.integer(path), as.integer(scaffold))
  bar <- if (file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
  f <- DBI::dbGetQuery(con, "SELECT element, start, stop, strand, tool FROM review_annotations
    WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  f <- f[which(startsWith(f$tool, "barrnap")), , drop = FALSE]
  cov <- vapply(seq_len(nrow(f)), function(i) {
    b <- if (nrow(bar)) bar[bar$element == f$element[i] & bar$strand == f$strand[i], , drop = FALSE] else bar
    ov <- if (nrow(b)) pmin(b$stop, f$stop[i]) - pmax(b$start, f$start[i]) + 1 else 0
    if (max(ov) < 1) 0 else b$model_cov[which.max(ov)]
  }, numeric(1))
  with_savepoint(con, for (i in seq_len(nrow(f))) {
    DBI::dbExecute(con, "UPDATE review_annotations SET model_coverage = ?
      WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(list(cov[i]), p, list(f$element[i])))
  })
  invisible(cov)
}
```

`extra_opts` is split in R and each word quoted, instead of handing the raw string to the shell: blast.nf's unquoted `${XO}` word-splits the same way, and a `;` in a stored option cannot run a command. Check it parses: `Rscript -e 'invisible(parse("R/recheck.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-recheck.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 56 ]`.

- [ ] **Step 5: Write the failing recheck_unit tests**

Append this to the end of `tests/testthat/test-recheck.R` (one blank line before it). Mid-run edits are injected by mocking the last check (`chimera_unit`), so the final guard, not the start check, has to discard the result. The happy path sets the copy to running with no counts (as `query_recheck` leaves it), and gives the original a best reference whose 18S differs at every tenth base, so curation's `low_identity` flag shows the reference alignment ran on the copy before curation; the 28S `truncated` flag shows the same for the new model coverage.

```r
test_that("recheck_unit reruns every check on the copy and marks it done", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  # running, as query_recheck leaves it, with no counts yet
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 4, n_features = NULL, missing = NULL
    WHERE", RC_UNIT))
  # the original's best reference: its 18S differs at every tenth base
  ref <- strsplit(RC_SEQ, "")[[1]]
  i <- seq(360L, 2150L, by = 10L)
  ref[i] <- chartr("ACGT", "CGTA", ref[i])
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates (ID, path, scaffold, rank, accession)
    VALUES ('S1', 1, 1, 1, 'REF1.1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence, length) VALUES ('REF1.1', ?, ?)",
    params = list(paste(ref, collapse = ""), length(ref)))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations (accession, element, start, stop, strand)
    VALUES ('REF1.1', '18S', 351, 2150, '+')")
  before <- rc_originals(con)
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "done")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 2L)
  expect_true(is.na(u$check_note))
  expect_gte(u$checked_at, t0)
  expect_gte(u$time_stamp, t0)
  expect_equal(u$n_features, 5L)
  expect_equal(u$missing, "")
  expect_equal(u$chimeric, 1L)
  expect_false(is.na(u$n_flags))
  expect_equal(u$itsx_warnings, "")
  expect_equal(u$blast_accession, "A1.1")
  expect_equal(rc_cov(con)[c("18S", "28S")], c(`18S` = 100, `28S` = 70))
  ra <- DBI::dbGetQuery(con, paste("SELECT element, pid FROM review_ref_alignment WHERE", RC_UNIT))
  expect_equal(ra$element, "18S")
  expect_lt(ra$pid, 95)
  # curation saw the new model coverage and reference alignment, and the
  # chimera check the new ITSx flag
  fl <- DBI::dbGetQuery(con, paste("SELECT element, flags FROM review_annotations WHERE", RC_UNIT))
  expect_equal(fl$flags[fl$element == "28S"], "truncated")
  expect_match(fl$flags[fl$element == "18S"], "low_identity", fixed = TRUE)
  expect_match(DBI::dbGetQuery(con, paste("SELECT reason FROM review_chimera WHERE", RC_UNIT))$reason,
    "ITSx chimeric", fixed = TRUE)
  for (t in c("review_blast_hits", "review_blast_screen", "review_screen", "review_chimera")) {
    expect_gt(rc_count(con, t), 0, label = t)
  }
  expect_equal(rc_rev(con), rev)
  expect_equal(rc_originals(con), before)
})

test_that("failed tools leave their values as copied and the contig still finishes", {
  for (status in list(c(barrnap = 1, itsx = 1, itsx_check = 0), NULL)) {
    con <- rc_project()
    RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
    tools <- rc_tools(status = status, finished = FALSE)
    expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", tools), "done")
    u <- rc_copy(con)
    expect_equal(u$check_switch, 2L)
    expect_equal(u$chimeric, 0L)
    expect_equal(u$itsx_warnings, "not_run(ITSx did not finish)")
    expect_equal(rc_cov(con)[c("18S", "28S")], c(`18S` = 100, `28S` = 100))
  }
})

test_that("an ITSx check that throws is stored as not_run and the contig still finishes", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  local_mocked_bindings(itsx_check_one = function(...) stop("bad ITSx table"), .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", rc_tools()), "done")
  u <- DBI::dbGetQuery(con, paste("SELECT check_switch, itsx_warnings, itsx_origin, itsx_scores
    FROM review_units WHERE", RC_UNIT))
  expect_equal(u$check_switch, 2L)
  expect_equal(u$itsx_warnings, "not_run(ITSx check failed: bad ITSx table)")
  expect_true(is.na(u$itsx_origin))
  expect_true(is.na(u$itsx_scores))
})

test_that("an error in a step marks the copy failed, under the rev guard", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 4, time_stamp = 0 WHERE", RC_UNIT))
  local_mocked_bindings(screen_one = function(...) stop("screen broke"), .package = "RiboPilot")
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "failed")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 3L)
  expect_equal(u$check_note, "screen broke")
  expect_gte(u$time_stamp, t0)
  expect_equal(rc_count(con, "review_chimera"), 0)

  # an edit before the error: the failure is not recorded
  RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 1)
  rev <- rc_rev(con)
  local_mocked_bindings(screen_one = function(con, id, path, scaffold, ...) {
    RiboPilot:::review_trim(con, id, path, scaffold, 1L, 5900L)
    stop("screen broke")
  }, .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "failed")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 1L)
  expect_true(is.na(u$check_note))
})

test_that("a copy that changed before the run starts is not rechecked", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5900L)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_screen"), 0)
  expect_equal(rc_copy(con)$check_switch, 1L)

  rev <- rc_rev(con)
  RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 0)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_screen"), 0)
  RiboPilot:::review_drop(con, "S1", 1L, 1L)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_units"), 0)
  expect_equal(rc_count(con, "review_screen"), 0)
})

test_that("the rev guard discards a result after an edit, On hold, or Restore then a new trim", {
  edits <- list(
    edit = function(con) RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5900L),
    on_hold = function(con) RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 0),
    restore_and_trim = function(con) {
      RiboPilot:::review_drop(con, "S1", 1L, 1L)
      Sys.sleep(0.01)
      RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5800L)
    })
  for (nm in names(edits)) {
    con <- rc_project()
    RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
    rev <- rc_rev(con)
    # the edit lands while the last check runs
    local_mocked_bindings(chimera_unit = function(con, ...) edits[[nm]](con), .package = "RiboPilot")
    expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale", info = nm)
    u <- rc_copy(con)
    expect_equal(u$check_switch, if (nm == "on_hold") 0L else 1L, info = nm)
    expect_true(is.na(u$checked_at), info = nm)
    if (nm != "on_hold") expect_false(identical(rc_rev(con), rev), info = nm)
  }
})

test_that("recheck_release_running puts copies left running back to not run", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, check_switch,
        check_note, rev, time_stamp) VALUES ('S1', 1, ?, 'ACGT', 4, ?, ?, CAST(? AS INTEGER), 5)",
      params = list(s, c(4L, 4L, 1L, 2L)[s], c(NA, "old", NA, NA)[s], "1760000000000"))
  }
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_release_running(con), 2L)
  r <- DBI::dbGetQuery(con, "SELECT check_switch, check_note, time_stamp FROM review_units ORDER BY scaffold")
  expect_equal(r$check_switch, c(1L, 1L, 1L, 2L))
  expect_equal(r$check_note, c("recheck did not finish", "recheck did not finish", NA, NA))
  expect_true(all(r$time_stamp[1:2] >= t0))
  expect_equal(r$time_stamp[3:4], c(5L, 5L))
  expect_equal(RiboPilot:::recheck_release_running(con), 0L)
})
```

The `Sys.sleep(0.01)` keeps the recreated copy's `rev` (a millisecond timestamp plus one) from ever equaling the discarded one.

- [ ] **Step 6: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-recheck.R")'`
Expected: FAIL, `[ FAIL 7 | WARN 0 | SKIP 0 | PASS 56 ]`: six tests error with `object 'recheck_unit' not found` and one with `object 'recheck_release_running' not found`.

- [ ] **Step 7: Add recheck_unit() and recheck_release_running()**

Append this to the end of `R/recheck.R` (one blank line before it):

```r
#' Recheck one contig's review copy (SPEC S7.3)
#'
#' Reruns every offline check on the copy in the order of the spec, then marks
#' it done, but only while the copy still has the `rev` the run was queued with
#' and is still queued or running (1 or 4); an edit, On hold, or Restore during
#' the run discards the result. An error in a step marks the copy failed under
#' the same guard. Opens no transaction around the steps: each check opens its
#' own.
#' @param rev The copy's `rev` as text, as query_recheck emitted it.
#' @param mode "local" (search the panel again) or "map" (map the original hits).
#' @param tools_dir Folder holding rp_rrna_tools.sh's output for the copy.
#' @return "done", "stale", or "failed", invisibly.
#' @noRd
recheck_unit <- function(con, id, path, scaffold, rev, mode, tools_dir, dbarg = NULL,
                         taxdump_dir = NULL, label_db = NULL, threads = 1L,
                         workdir = tempfile("recheck_")) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  unit <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  guard <- paste(unit, "AND rev = CAST(? AS INTEGER) AND check_switch IN (1, 4)")
  p <- list(id, path, scaffold)
  pg <- c(p, list(as.character(rev)))
  if (DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM review_units", guard), params = pg)$n == 0) {
    return(invisible("stale"))
  }
  st <- file.path(tools_dir, "unit.tools_status")
  st <- if (file.exists(st)) trimws(readLines(st, warn = FALSE)) else character(0)
  ok <- function(tool) paste0(tool, "=0") %in% st
  err <- tryCatch({
    recheck_hits(con, id, path, scaffold, mode, dbarg = dbarg, taxdump_dir = taxdump_dir,
      threads = threads, workdir = file.path(workdir, "blast"))
    screen_one(con, id, path, scaffold, force = TRUE, review = TRUE)
    cnt <- review_counts(con, id, path, scaffold)
    DBI::dbExecute(con, paste("UPDATE review_units SET n_features = ?, missing = ?", unit),
      params = c(list(cnt$n_features, cnt$missing), p))
    recheck_model_coverage(con, id, path, scaffold,
      if (ok("barrnap")) file.path(tools_dir, "unit.barrnap.gff"))
    if (ok("itsx")) {
      pos <- readLines(file.path(tools_dir, "unit.itsx.positions.txt"), warn = FALSE)
      DBI::dbExecute(con, paste("UPDATE review_units SET chimeric = ?", unit),
        params = c(list(as.integer(itsx_is_chimeric(pos))), p))
    }
    tryCatch(itsx_check_one(con, id, path, scaffold,
      out_prefix = file.path(tools_dir, "unit.itsx_check"), review = TRUE),
      error = function(e) itsx_store_not_run(con, id, path, scaffold,
        paste("ITSx check failed:", conditionMessage(e)), review = TRUE))
    compute_ref_alignment(con, id, path, scaffold, review = TRUE)
    curate_one(con, id, path, scaffold, review = TRUE)
    chimera_unit(con, id, path, scaffold, dbarg = dbarg, taxdump_dir = taxdump_dir,
      label_db = label_db, threads = threads, workdir = file.path(workdir, "chimera"), review = TRUE)
    NULL
  }, error = function(e) conditionMessage(e))
  ts <- as.integer(Sys.time())
  if (is.null(err)) {
    n <- DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 2, checked_at = ?,
      check_note = NULL, time_stamp = ?", guard), params = c(list(ts, ts), pg))
    return(invisible(if (n == 1) "done" else "stale"))
  }
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 3, check_note = ?, time_stamp = ?",
    guard), params = c(list(err, ts), pg))
  invisible("failed")
}

#' After an in-app recheck run: copies left running go back to not run
#'
#' @return The number of copies released.
#' @noRd
recheck_release_running <- function(con) {
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 1, check_note = 'recheck did not finish',
    time_stamp = ? WHERE check_switch = 4", params = list(as.integer(Sys.time())))
}
```

The start check uses the same guard as the final write (`rev` and `check_switch IN (1, 4)`), so a copy put On hold before its task started is skipped without running any check. The stock ITSx check always runs (its own summary file decides finished or not, S6.3); only barrnap and the tuned ITSx run are gated on `unit.tools_status`. Check it parses: `Rscript -e 'invisible(parse("R/recheck.R"))'`.

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-recheck.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 127 ]`.

Then the files that scan every R file (duplicate definitions, banned strings):
Run: `Rscript -e 'devtools::load_all("."); for (f in c("no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: both `FAIL 0`.

Then confirm the new files are ASCII only:
Run: `grep -nP '[^\x00-\x7F]' R/recheck.R tests/testthat/test-recheck.R`
Expected: no output.

- [ ] **Step 9: Commit**

```bash
git add R/recheck.R tests/testthat/test-recheck.R
git commit -m "Add recheck_unit for review copies"
```


### Task 11: Recheck workflow and command

**Files:**
- Create: `inst/nextflow/main_recheck.nf`
- Create: `inst/nextflow/modules/recheck.nf`
- Create: `inst/nextflow/bin/rp_recheck.R` (git mode 100755)
- Modify: `R/recheck.R` (append `RECHECK_QUEUE_SQL` after `recheck_release_running`, which Task 10 created)
- Modify: `R/nf_command.R` (header comment, lines 1-7; append `nf_recheck_command` and `nf_recheck_pending` after `nf_pending_count`)
- Modify: `inst/nextflow/nextflow.config` (OOM selector, line 59)
- Modify: `inst/config.NMNH_Hydra` (memory selector, line 53)
- Modify: `tests/testthat/test-nextflow.R` (lines 124-128, line 1576, new tests appended at the end)
- Modify: `tests/testthat/test-nf-pending.R` (comments at lines 1-3 and 58-59, new tests appended at the end)
- Modify: `tests/testthat/test-utils.R` (lines 70-75)
- Test: `tests/testthat/test-nf-recheck-run.R` (new)

**Interfaces:**
- Consumes:
  - Task 9: `rp_rrna_tools.sh <fasta> <tool> <kingdom> <itsx_taxon> <itsx_opts> <rrna_opts> <cpus>` (exactly 7 arguments; writes `unit.barrnap.gff`, `unit.itsx.positions.txt`, `unit.itsx_check.*`, `unit.tools_status`; exits with the barrnap/pybarrnap status), and `annotate_tools`' options read (pipeline_r.nf, from `read TOOL KINGDOM ITSXT ITSXO <<<` through the `RRNAO` query), unchanged by Task 9.
  - Task 10: `recheck_unit(con, id, path, scaffold, rev, mode, tools_dir, dbarg = NULL, taxdump_dir = NULL, label_db = NULL, threads = 1L, workdir = tempfile("recheck_"))` returning invisible `"done"`, `"stale"`, or `"failed"`; `recheck_release_running(con)` (the last function of `R/recheck.R`).
  - Task 1: the `review_units` table in `new_db()`. Task 4 (test only): `review_trim(con, id, path, scaffold, from, to)`, `review_set_state(con, units, state)`.
- Produces:
  - `RECHECK_QUEUE_SQL` (`R/recheck.R`, `@noRd`): one SELECT of `ID, path, scaffold, rev (text), cpus, mem, mode`; no trailing semicolon, no `$`, `"`, `\`, or `;`.
  - `nf_recheck_command(path)` (`R/nf_command.R`, `@noRd`): `c("run", <main_recheck.nf>, "-c", <project>/.config, "-work-dir", <project>/work_recheck, "--db", <project>/.sqlite, "--projectDir", <project>, ["--rawDir", <rawDir>], nf_pkg_args(), "--taxdump_dir", <project>/.taxonomy, ["--blastdb", <panel>], <container bind options as nf_command>)`. No `-resume`, `--minDepth`, `--ref_dir`, or `--ncbi_api_key`; no stop on a blank rawDir; no side effects; no empty values. Launch it from `<project>/.recheck`.
  - `nf_recheck_pending(path)` (`R/nf_command.R`, `@noRd`): integer, the row count of `RECHECK_QUEUE_SQL`, one `open_project_db` call.
  - `inst/nextflow/main_recheck.nf`: `query_recheck` -> `recheck_tools` -> `recheck`, every process `cache false`.
  - `inst/nextflow/bin/rp_recheck.R <db> <id> <path> <scaffold> <rev> <mode> <tools_dir> [dbarg] [taxdump_dir] [label_db]`: exits 0 for done, stale, and failed and prints `recheck: <id> <path> <scaffold> - <result>`.

This task adds the second Nextflow workflow that reruns the offline checks on review copies (S7.1, S7.2). `query_recheck` takes the queued copies and marks them running (4) in one sqlite3 transaction, emitting each copy's `rev` and a BLAST mode (S5.4, S7.2). The selection is `RECHECK_QUEUE_SQL`, kept word for word in `modules/recheck.nf` and counted by `nf_recheck_pending` (S7.4); a test runs the query through the machine's sqlite3 CLI (3.37 here, older than RSQLite's SQLite). `recheck_tools` writes the copy's sequence only while its `rev` is unchanged and calls Task 9's `rp_rrna_tools.sh` with the options `annotate_tools` reads; a tool failure is recorded for `recheck_unit`, and only out-of-memory exits fail the task so the nextflow.config retry can give it more memory (Review Focus 3). `recheck` resolves the local panel and labels as the main pipeline's `chimera` process does and calls `rp_recheck.R`, which calls Task 10's `recheck_unit()`. `nf_recheck_command()` builds the launch command (S7.4); the real-run test launches it from `<project>/.recheck` with `RIBOPILOT_NO_LOCAL_BLAST=1`, so a contig whose original search was local falls back to mapped hits and still finishes without calling NCBI (Review Focus 4).

- [ ] **Step 1: Write the failing tests**

In `tests/testthat/test-nextflow.R`, the list of workflow files allows `main_recheck.nf` (a second workflow on purpose, S7.1).

Find:

```r
  # main_reference.nf is the reference-database workflow, not a second WF1
  wf <- setdiff(list.files(nf_dir, pattern = "\\.nf$"),
    c("nextflow.config", "main_reference.nf"))
```

Replace with:

```r
  # main_reference.nf builds reference databases and main_recheck.nf rechecks
  # review copies; neither is a second WF1
  wf <- setdiff(list.files(nf_dir, pattern = "\\.nf$"),
    c("nextflow.config", "main_reference.nf", "main_recheck.nf"))
```

The OOM selector gains `recheck_tools` on purpose (S7.2: same retry as `annotate_tools`).

Find:

```r
  expect_match(cfg, "withName: 'assemble|fallback_assemble|preprocess|annotate_tools|find_rrna_scan'",
```

Replace with:

```r
  expect_match(cfg, "withName: 'assemble|fallback_assemble|preprocess|annotate_tools|find_rrna_scan|recheck_tools'",
```

Append the recheck tests at the end of the file.

Find:

```r
  main <- readLines(system.file("nextflow", "main.nf", package = "RiboPilot"))
  expect_true(any(grepl("r.rrna_kingdom, r.rrna_opts)", main, fixed = TRUE)))
})
```

Replace with:

```r
  main <- readLines(system.file("nextflow", "main.nf", package = "RiboPilot"))
  expect_true(any(grepl("r.rrna_kingdom, r.rrna_opts)", main, fixed = TRUE)))
})

# ---- the recheck workflow: main_recheck.nf and modules/recheck.nf ----

recheck_nf <- readLines(system.file("nextflow", "modules", "recheck.nf", package = "RiboPilot"))
main_recheck_nf <- readLines(system.file("nextflow", "main_recheck.nf", package = "RiboPilot"))

recheck_body <- function(lines, proc) {
  start <- grep(paste0("^process ", proc, " \\{"), lines)
  expect_length(start, 1)
  lines[start:(start - 1 + which(lines[start:length(lines)] == "}")[1])]
}

test_that("main_recheck.nf feeds query_recheck through recheck_tools into recheck", {
  expect_true(any(main_recheck_nf ==
    "include { query_recheck; recheck_tools; recheck } from './modules/recheck.nf'"))
  code <- paste(main_recheck_nf[!grepl("^\\s*//", main_recheck_nf)], collapse = "\n")
  expect_match(code, "query_recheck().splitCsv(header: true, sep: '\\t')", fixed = TRUE)
  expect_match(code, "tuple(r.ID, r.path, r.scaffold, r.rev, r.cpus, r.mem, r.mode)", fixed = TRUE)
  expect_match(code, "recheck_tools(units)", fixed = TRUE)
  expect_match(code, "recheck(recheck_tools.out)", fixed = TRUE)
  # every field the workflow reads is a column query_recheck writes
  hdr <- grep("printf 'ID", recheck_body(recheck_nf, "query_recheck"), value = TRUE)
  expect_length(hdr, 1)
  cols <- strsplit(sub(".*printf '", "", sub("\\\\n' > .*", "", hdr)), "\\\\t")[[1]]
  expect_equal(cols, c("ID", "path", "scaffold", "rev", "cpus", "mem", "mode"))
  used <- sub("^r\\.", "", regmatches(code, gregexpr("\\br\\.[A-Za-z_]+", code))[[1]])
  expect_true(all(used %in% cols))
  # uncached, and the per-contig processes are tagged by contig
  for (p in c("query_recheck", "recheck_tools", "recheck")) {
    expect_true(any(grepl("^\\s*cache false$", recheck_body(recheck_nf, p))), info = p)
  }
  for (p in c("recheck_tools", "recheck")) {
    body <- recheck_body(recheck_nf, p)
    expect_true(any(grepl('tag "${id}.${upath}.${uscaf}"', body, fixed = TRUE)), info = p)
    expect_true(any(grepl("errorStrategy 'ignore'", body, fixed = TRUE)), info = p)
  }
  expect_true(any(grepl("executor 'local'", recheck_body(recheck_nf, "query_recheck"), fixed = TRUE)))
  expect_true(all(c("params.blastdb = ''", "params.taxdump_dir = ''") %in% recheck_nf))
})

test_that("query_recheck selects RECHECK_QUEUE_SQL word for word and marks it running", {
  body <- paste(recheck_body(recheck_nf, "query_recheck"), collapse = "\n")
  sql <- regmatches(body, regexpr("CREATE TEMP TABLE q AS\\s*\\K[^;]*", body, perl = TRUE))
  norm <- function(x) gsub("\\s+", " ", trimws(x))
  expect_equal(norm(sql), norm(RiboPilot:::RECHECK_QUEUE_SQL))
  # it travels inside a Nextflow string and a double-quoted shell argument
  expect_false(any(strsplit(RiboPilot:::RECHECK_QUEUE_SQL, "")[[1]] %in% c("$", "\"", "\\", ";")))
  expect_match(body, "BEGIN IMMEDIATE;", fixed = TRUE)
  expect_match(body, "UPDATE review_units SET check_switch = 4", fixed = TRUE)
  expect_match(body, "COMMIT;\" >> recheck_units.tsv", fixed = TRUE)
})

test_that("recheck_tools reads the tool options exactly as annotate_tools does", {
  opts <- function(lines, proc) {
    body <- recheck_body(lines, proc)
    from <- grep("read TOOL KINGDOM ITSXT ITSXO <<<", body, fixed = TRUE)
    to <- grep("RRNAO=\"", body, fixed = TRUE)
    expect_length(from, 1)
    expect_length(to, 1)
    body[from:(to - 1 + grep("\")\"$", body[to:length(body)])[1])]
  }
  got <- opts(recheck_nf, "recheck_tools")
  expect_gt(length(got), 5)
  expect_identical(got, opts(pipeline_r_nf, "annotate_tools"))
})

test_that("recheck_tools runs the tools only on the queued rev and fails only on out-of-memory exits", {
  body <- recheck_body(recheck_nf, "recheck_tools")
  txt <- paste(body, collapse = "\n")
  expect_match(txt, "'>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM review_units", fixed = TRUE)
  expect_match(txt, "AND rev=${rev}", fixed = TRUE)
  expect_match(txt, "check_switch IN (1,4) AND sequence IS NOT NULL", fixed = TRUE)
  expect_match(txt, "if [ -s unit.fasta ]; then", fixed = TRUE)
  call <- grep("rp_rrna_tools.sh ", body, value = TRUE, fixed = TRUE)
  expect_length(call, 1)
  expect_equal(trimws(call), paste('rp_rrna_tools.sh unit.fasta "\\$TOOL" "\\$KINGDOM" "\\$ITSXT"',
    '"\\$ITSXO" "\\$RRNAO" ${task.cpus} || rc=\\$?'))
  codes <- grep('case "\\$rc" in', body, value = TRUE, fixed = TRUE)
  expect_length(codes, 1)
  codes <- as.integer(strsplit(sub(".* in ([0-9|]+)\\).*", "\\1", codes), "|", fixed = TRUE)[[1]])
  expect_setequal(codes, RiboPilot:::NF_RETRY_CODES)
  expect_match(txt, "path('unit.itsx_check*'), optional: true", fixed = TRUE)
})

test_that("the recheck process resolves the panel and labels as the chimera process does", {
  block <- function(lines, proc) {
    body <- recheck_body(lines, proc)
    from <- grep('DBARG=""', body, fixed = TRUE)
    to <- grep("2>/dev/null || echo \"\")", body, fixed = TRUE)
    expect_length(from, 1)
    body[from:to[length(to)]]
  }
  expect_identical(block(recheck_nf, "recheck"), block(pipeline_r_nf, "chimera"))
  call <- gsub("\\s+", " ", paste(recheck_body(recheck_nf, "recheck"), collapse = " "))
  expect_match(call, paste("RP_THREADS=${task.cpus} RIBOPILOT_SRC='${params.ribopilot_src}' rp_recheck.R \\",
    "'${params.db}' '${id}' ${upath} ${uscaf} ${rev} ${mode} . \"\\$DBARG\" '${params.taxdump_dir}' \"\\$LABELS\""),
    fixed = TRUE)
})

test_that("nf_recheck_command runs main_recheck.nf in work_recheck, without -resume", {
  skip_if_not_installed("withr")
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  mk_ref_config_project(proj, refs, blastdb = "/opt/rp/ribopilot_rrna")
  writeLines(c("params {", "  rawDir = '/data/raw'", "  minDepth = 5",
    "  ncbi_api_key = 'abc123'", "}"), file.path(proj, ".config"))
  a <- RiboPilot:::nf_recheck_command(proj)
  val <- function(flag) a[which(a == flag) + 1]
  expect_equal(a[1:2], c("run", system.file("nextflow", "main_recheck.nf", package = "RiboPilot")))
  expect_true(file.exists(a[2]))
  expect_equal(val("-c"), file.path(proj, ".config"))
  expect_equal(val("-work-dir"), file.path(proj, "work_recheck"))
  expect_equal(val("--db"), file.path(proj, ".sqlite"))
  expect_equal(val("--projectDir"), proj)
  expect_equal(val("--rawDir"), "/data/raw")
  expect_equal(val("--taxdump_dir"), file.path(proj, ".taxonomy"))
  expect_equal(val("--blastdb"), "/opt/rp/ribopilot_rrna")
  expect_equal(val("--ribopilot_version"), as.character(utils::packageVersion("RiboPilot")))
  for (f in c("-resume", "--minDepth", "--ref_dir", "--ncbi_api_key")) expect_false(f %in% a, info = f)
  expect_true(all(nzchar(a)))
  # the main pipeline's command for the same project does carry them
  expect_true(all(c("--minDepth", "--ref_dir", "--ncbi_api_key") %in%
    RiboPilot:::nf_command(proj, resume = FALSE)))
})

test_that("nf_recheck_command leaves out a blank rawDir and blastdb instead of refusing", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  # reads but no rawDir: the main pipeline refuses, the recheck reads no reads
  mk_rawdir_project(proj, raw = "NA", r1 = "a.fq")
  expect_error(RiboPilot:::nf_command(proj, resume = FALSE), "No rawDir set")
  a <- RiboPilot:::nf_recheck_command(proj)
  expect_false("--rawDir" %in% a)
  expect_false("--blastdb" %in% a)
  expect_true(all(nzchar(a)))
})

test_that("nf_recheck_command mounts the reference folders as nf_command does", {
  skip_if_not_installed("withr")
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  dir.create(file.path(refs, "seeds")); dir.create(file.path(refs, "labels"))
  s <- file.path(refs, "seeds", "x_seeds.fasta"); l <- file.path(refs, "labels", "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l)
  a <- RiboPilot:::nf_recheck_command(proj)
  b <- RiboPilot:::nf_command(proj, resume = FALSE)
  opt <- "-process.containerOptions"
  expect_length(which(a == opt), 1L)
  expect_equal(a[which(a == opt) + 1], b[which(b == opt) + 1])
  expect_true(grepl(file.path(refs, "seeds"), a[which(a == opt) + 1], fixed = TRUE))
})

test_that("building the recheck command records no panel provenance", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(d, ".config"))
  panel_dir <- withr::local_tempdir()
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(panel_dir, "manifest.json"))
  prefix <- file.path(panel_dir, "ribopilot_rrna")
  file.create(paste0(prefix, ".nin"))
  mk_ref_config_project(d, NA_character_, blastdb = prefix)
  expect_true("--blastdb" %in% RiboPilot:::nf_recheck_command(d))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::blastdb_provenance(con), "")
})
```

(`pipeline_r_nf`, `mk_ref_config_project`, `mk_rawdir_project`, and `mk_mount_project` are already defined in test-nextflow.R, at lines 41, 356, 1441, and 238; `NF_RETRY_CODES` is in R/nf_run.R:16.)

In `tests/testthat/test-nf-pending.R`, the header comment and one test comment say there is one workflow; they now name both Update buttons.

Find:

```r
# One workflow, one Update button, so one number: samples queued for assembly
# plus units queued for annotation only. It must agree with what the two
# Nextflow queries select, or the modal promises work the run will not do.
```

Replace with:

```r
# Each Update button shows one number, and it must agree with what its
# workflow's queries select, or the run window promises work the run will not
# do. The Assemble step's Update counts samples queued for assembly plus units
# queued for annotation only (main.nf); Annotation Review's counts the queued
# review copies (main_recheck.nf).
```

Find:

```r
  # There is one workflow; a caller passing "annotate" would silently get a
  # count for something that no longer exists.
```

Replace with:

```r
  # The main pipeline is one workflow; a caller passing "annotate" would
  # silently get a count for something that no longer exists.
```

Append the queue tests at the end of the file. The regex in "every nf_command() caller matches its signature" (`nf_command\\([^)]*\\)`) does not match `nf_recheck_command(`, so the last new test checks the recheck builder's callers itself (`app_server` and Task 13's `build_recheck_submit_script`; it skips a caller that does not exist yet).

Find:

```r
        expect_true(a %in% ok, info = paste(nm, cl, a))
      }
    }
  }
})
```

Replace with:

```r
        expect_true(a %in% ok, info = paste(nm, cl, a))
      }
    }
  }
})

# ---- Annotation Review's Update: the recheck queue ----

recheck_queue_ids <- c("QUEUED", "RUNNING", "LOCAL", "NO_BLAST", "DONE", "HELD", "FAILED",
  "LOCKED", "DROPPED", "ASM_OPEN", "TAX_OPEN", "NO_COPY")

# one contig per sample; only the first four are queued for a recheck
fill_recheck_queue <- function(con) {
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = CASE WHEN ID = 'TAX_OPEN' THEN 0 ELSE 1 END")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = CASE WHEN ID = 'ASM_OPEN' THEN 0 ELSE 1 END")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2,
    annotate_lock = CASE WHEN ID = 'LOCKED' THEN 1 ELSE 0 END")
  DBI::dbExecute(con, "UPDATE annotate_opts SET cpus = 3, memory = 5 WHERE annotate_opts = 'default'")
  DBI::dbExecute(con, "INSERT INTO blast_opts (blast_opts, run_blast) VALUES ('off', 0)")
  DBI::dbExecute(con, "UPDATE assemble SET blast_opts = 'off' WHERE ID = 'NO_BLAST'")
  for (id in recheck_queue_ids) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length, ignore)
      VALUES (?, 1, 1, 'ACGTACGT', 8, ?)", params = list(id, as.integer(id == "DROPPED")))
  }
  cs <- c(QUEUED = 1L, RUNNING = 4L, LOCAL = 1L, NO_BLAST = 1L, DONE = 2L, HELD = 0L,
    FAILED = 3L, LOCKED = 1L, DROPPED = 1L, ASM_OPEN = 1L, TAX_OPEN = 1L)
  for (i in seq_along(cs)) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length,
      trim_start, trim_end, edited, check_switch, rev, time_stamp)
      VALUES (?, 1, 1, 'ACGTACGT', 8, 1, 8, 1, ?, CAST(? AS INTEGER), 0)",
      params = list(names(cs)[i], cs[[i]], format(1791234567890 + i, scientific = FALSE)))
  }
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('QUEUED', 1, 1, 'x', 'remote'), ('LOCAL', 1, 1, 'x', 'local'),
           ('NO_BLAST', 1, 1, 'x', 'local')")
}

test_that("the recheck queue takes queued and running copies of locked, kept contigs", {
  con <- chimera_test_con(recheck_queue_ids)
  on.exit(DBI::dbDisconnect(con))
  fill_recheck_queue(con)
  got <- DBI::dbGetQuery(con, RiboPilot:::RECHECK_QUEUE_SQL)
  expect_named(got, c("ID", "path", "scaffold", "rev", "cpus", "mem", "mode"))
  got <- got[order(got$ID), ]
  expect_equal(got$ID, c("LOCAL", "NO_BLAST", "QUEUED", "RUNNING"))
  # rev as exact text: the workflow hands it back to the guard unchanged
  expect_identical(got$rev, c("1791234567893", "1791234567894", "1791234567891", "1791234567892"))
  # local only when the sample searches and the contig's original search was local
  expect_equal(got$mode, c("local", "map", "map", "map"))
  expect_equal(unique(got$cpus), 3L)
  expect_equal(unique(got$mem), 5L)
  expect_identical(RiboPilot:::nf_recheck_pending(dirname(DBI::dbGetInfo(con)$dbname)), 4L)
})

test_that("query_recheck runs on the sqlite3 CLI and marks what it takes as running", {
  skip_if(!nzchar(Sys.which("sqlite3")), "sqlite3 not installed")
  nf <- readLines(system.file("nextflow", "modules", "recheck.nf", package = "RiboPilot"))
  start <- grep("^process query_recheck \\{", nf)
  body <- nf[start:(start - 1 + which(nf[start:length(nf)] == "}")[1])]
  from <- grep("sqlite3 -cmd '.timeout 10000' -separator", body, fixed = TRUE)
  to <- grep("COMMIT;\" >> recheck_units.tsv", body, fixed = TRUE)
  expect_length(from, 1)
  expect_length(to, 1)
  sql <- sub("\" >> recheck_units.tsv$", "", paste(body[(from + 1):to], collapse = "\n"))

  con <- chimera_test_con(recheck_queue_ids)
  on.exit(DBI::dbDisconnect(con))
  fill_recheck_queue(con)
  db <- DBI::dbGetInfo(con)$dbname
  want <- DBI::dbGetQuery(con, RiboPilot:::RECHECK_QUEUE_SQL)
  want[] <- lapply(want, as.character)
  t0 <- as.integer(Sys.time())
  # the machine's sqlite3 CLI, which can be older than the SQLite inside RSQLite
  res <- processx::run("sqlite3", c("-cmd", ".timeout 10000", "-separator", "\t", db, sql))
  out <- utils::read.delim(text = res$stdout, header = FALSE, col.names = names(want),
    colClasses = "character")
  expect_equal(out[order(out$ID), ], want[order(want$ID), ], ignore_attr = TRUE)

  st <- DBI::dbGetQuery(con, "SELECT ID, check_switch, time_stamp FROM review_units ORDER BY ID")
  taken <- st$ID %in% want$ID
  expect_equal(st$check_switch[taken], rep(4L, 4))
  expect_true(all(st$time_stamp[taken] >= t0))
  expect_equal(stats::setNames(st$check_switch[!taken], st$ID[!taken]),
    c(ASM_OPEN = 1L, DONE = 2L, DROPPED = 1L, FAILED = 3L, HELD = 0L, LOCKED = 1L, TAX_OPEN = 1L))
  expect_true(all(st$time_stamp[!taken] == 0))
  # running copies stay in the queue, so the next Update takes them again
  expect_identical(RiboPilot:::nf_recheck_pending(dirname(db)), 4L)
})

test_that("nf_recheck_command and nf_recheck_pending take only the project", {
  expect_named(formals(RiboPilot:::nf_recheck_command), "path")
  expect_named(formals(RiboPilot:::nf_recheck_pending), "path")
  # the callers the app launch adds pass the project and nothing else
  ns <- asNamespace("RiboPilot")
  for (nm in c("app_server", "build_recheck_submit_script")) {
    fn <- get0(nm, envir = ns, inherits = FALSE)
    if (is.null(fn)) next
    txt <- paste(deparse(fn), collapse = " ")
    for (cl in regmatches(txt, gregexpr("nf_recheck_command\\([^)]*\\)", txt))[[1]]) {
      args <- trimws(strsplit(sub("^nf_recheck_command\\(", "", sub("\\)$", "", cl)), ",")[[1]])
      expect_length(args, 1)
      expect_false(grepl("=", sub("^path\\s*=", "", args), fixed = TRUE), label = paste(nm, cl))
    }
  }
})
```

In `tests/testthat/test-utils.R`, `nf_recheck_pending` opens the project with the busy timeout too, so the count goes from 9 to 10.

Find:

```r
  RiboPilot:::nf_pending_count(dir)
  RiboPilot:::project_has_reads(dir)
```

Replace with:

```r
  RiboPilot:::nf_pending_count(dir)
  RiboPilot:::nf_recheck_pending(dir)
  RiboPilot:::project_has_reads(dir)
```

Find:

```r
  expect_equal(opened, 9L)
```

Replace with:

```r
  expect_equal(opened, 10L)
```

Create `tests/testthat/test-nf-recheck-run.R`. The first test runs `recheck_tools` alone with stub tools (Review Focus 3): a pybarrnap exit 1 is recorded and the outputs still go on, an exit 137 is retried with more memory, an exit 137 on every try drops the contig, and a copy queued with an older `rev` never reaches the tools. The second runs the real `main_recheck.nf` with the app's own command from `<project>/.recheck` on a project the real `main.nf` has run, with `RIBOPILOT_NO_LOCAL_BLAST=1` (Review Focus 4). The fake tools answer from one fixed 6310 bp contig (helper-fake-tools.R), so the trim keeps the 5' end (from = 1).

```r
# The recheck workflow with its tools faked: recheck_tools alone, for its exit
# rule and its rev guard; then the real main_recheck.nf, launched with the
# app's own command from <project>/.recheck, on a project the real main.nf has
# run (helper-fake-tools.R).

test_that("recheck_tools records tool failures, retries out-of-memory exits, and skips stale copies", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  skip_if(!nzchar(Sys.which("sqlite3")), "sqlite3 not installed")
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  con <- chimera_test_con("S1")
  db <- DBI::dbGetInfo(con)$dbname
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length)
      VALUES ('S1', 1, ?, 'ACGTACGTAC', 10)", params = list(s))
    DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold, annotate_opts)
      VALUES ('S1', 1, ?, 'default')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length,
      trim_start, trim_end, edited, check_switch, rev, time_stamp)
      VALUES ('S1', 1, ?, 'ACGTACGT', 8, 1, 8, 1, 4, CAST(? AS INTEGER), 0)",
      params = list(s, "1791234567890"))
  }
  DBI::dbDisconnect(con)

  root <- withr::local_tempdir()
  seen <- shQuote(file.path(root, "seen"))
  # contig 1: pybarrnap fails; contig 2: killed for memory once;
  # contig 4: killed for memory every time
  writeLines(c("#!/usr/bin/env bash",
    'u=$(head -n 1 "${@: -1}" | cut -c2-)',
    sprintf('echo "$u" >> %s', seen),
    sprintf('n=$(grep -cx "$u" %s)', seen),
    'case "$u.$n" in S1.1.1.*) exit 1 ;; S1.1.2.1) exit 137 ;; S1.1.4.*) exit 137 ;; esac'),
    file.path(root, "pybarrnap"))
  writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(root, "ITSx"))
  Sys.chmod(file.path(root, c("pybarrnap", "ITSx")), "0755")
  # contig 3 is queued with an older rev than the copy holds now
  writeLines(c(
    sprintf("include { recheck_tools } from '%s'", file.path(nf_dir, "modules", "recheck.nf")),
    "workflow {",
    "    recheck_tools(Channel.of(",
    "        tuple('S1', 1, 1, '1791234567890', 1, 1, 'map'),",
    "        tuple('S1', 1, 2, '1791234567890', 1, 1, 'local'),",
    "        tuple('S1', 1, 3, '1791234567889', 1, 1, 'map'),",
    "        tuple('S1', 1, 4, '1791234567890', 1, 1, 'map')))",
    "    recheck_tools.out.view { \"EMITTED ${it[0]}.${it[1]}.${it[2]} ${it[3]} ${it[4]} ${it[5].text.readLines().join(',')}\" }",
    "}"), file.path(root, "main.nf"))

  res <- processx::run("nextflow", c("run", file.path(root, "main.nf"),
    "-c", file.path(nf_dir, "nextflow.config"), "-work-dir", file.path(root, "work"),
    "-ansi-log", "false", "--db", db), wd = root,
    env = c(RiboPilot:::nf_env(),
      PATH = paste(root, file.path(nf_dir, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)
  # the fasta header is ID.path.scaffold; the stale copy never reached the tools
  expect_equal(sort(readLines(file.path(root, "seen"))),
    c("S1.1.1", "S1.1.2", "S1.1.2", rep("S1.1.4", 4)))
  expect_match(res$stdout, "EMITTED S1.1.1 1791234567890 map barrnap=1,itsx=0,itsx_check=0", fixed = TRUE)
  expect_match(res$stdout, "EMITTED S1.1.2 1791234567890 local barrnap=0,itsx=0,itsx_check=0", fixed = TRUE)
  expect_no_match(res$stdout, "EMITTED S1.1.3", fixed = TRUE)
  expect_no_match(res$stdout, "EMITTED S1.1.4", fixed = TRUE)
})

test_that("main_recheck.nf rechecks queued copies from .recheck, mapped and with no panel", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
  tools <- local_fake_pipeline_tools()
  p <- file.path(withr::local_tempdir(), "proj")
  suppressMessages(suppressWarnings(RiboPilot::new_test_project(p, n = 0, executor = "conda",
    custom_seeds_db = tools$scaffold, custom_labels_db = tools$scaffold, Rproj = FALSE, fetch_metadata = FALSE)))
  db <- file.path(p, ".sqlite")
  con <- RiboPilot:::open_project_db(db)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  DBI::dbDisconnect(con)
  no_panels <- withr::local_tempdir()
  path_env <- paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels, PATH = path_env),
    error_on_status = FALSE, timeout = 600, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)

  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  u <- q("SELECT ID, path, scaffold FROM annotate ORDER BY ID")
  expect_equal(u$ID, c("CHIM00", "CHIM01"))
  # CHIM00: 3' end trimmed, so its hits are mapped onto the shorter copy.
  # CHIM01: set to Ready to run untouched; its original search was local, and
  # the recheck has no panel, so it falls back to mapping too.
  RiboPilot:::review_trim(con, "CHIM00", u$path[1], u$scaffold[1], 1L, 6000L)
  expect_true(is.na(RiboPilot:::review_set_state(con, u[2, ], 1L)))
  DBI::dbExecute(con, "UPDATE blast_searched SET source = 'local' WHERE ID = 'CHIM01'")
  before <- q("SELECT ID, CAST(rev AS TEXT) AS rev, check_switch FROM review_units ORDER BY ID")
  expect_equal(before$check_switch, c(1L, 1L))
  expect_identical(RiboPilot:::nf_recheck_pending(p), 2L)
  orig <- q("SELECT ID, length, sequence FROM assemblies ORDER BY ID")

  # every blastn call is logged, then answered by the fake
  spy <- withr::local_tempdir()
  writeLines(c("#!/usr/bin/env bash",
    sprintf('echo "$*" >> %s', shQuote(file.path(spy, "calls"))),
    sprintf('exec %s "$@"', shQuote(file.path(tools$bin, "blastn")))), file.path(spy, "blastn"))
  Sys.chmod(file.path(spy, "blastn"), "0755")
  rc_dir <- file.path(p, ".recheck")
  dir.create(rc_dir)
  res2 <- processx::run("nextflow", RiboPilot:::nf_recheck_command(p), wd = rc_dir,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels, RIBOPILOT_NO_LOCAL_BLAST = "1",
      PATH = paste(spy, path_env, sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 600, cleanup_tree = TRUE)
  expect_equal(res2$status, 0L, info = res2$stderr)

  after <- q("SELECT ID, CAST(rev AS TEXT) AS rev, check_switch, check_note, checked_at,
    length, n_features, itsx_warnings, blast_accession FROM review_units ORDER BY ID")
  expect_equal(after$check_switch, c(2L, 2L), info = paste(after$check_note, collapse = "; "))
  expect_identical(after$rev, before$rev)
  expect_false(anyNA(after$checked_at))
  expect_equal(after$length, c(6000L, 6310L))
  expect_equal(after$n_features, c(5L, 5L))
  # the stock ITSx output reached the check (a run it cannot read is not_run)
  expect_false(any(is.na(after$itsx_warnings) | startsWith(after$itsx_warnings, "not_run(")))
  expect_equal(after$blast_accession, c("FAKE0001.1", "FAKE0001.1"))
  expect_equal(q("SELECT qend FROM review_blast_hits ORDER BY ID")$qend, c(6000L, 6310L))
  for (tbl in c("review_screen", "review_chimera")) {
    expect_equal(q(paste("SELECT ID FROM", tbl, "ORDER BY ID"))$ID, c("CHIM00", "CHIM01"), info = tbl)
  }
  calls <- if (file.exists(file.path(spy, "calls"))) readLines(file.path(spy, "calls")) else character(0)
  expect_false(any(grepl("-remote", calls, fixed = TRUE)))
  expect_identical(q("SELECT ID, length, sequence FROM assemblies ORDER BY ID"), orig)
  # its own launch folder and work folder: the main pipeline's history still
  # holds its one run
  expect_true(dir.exists(file.path(p, "work_recheck")))
  expect_true(file.exists(file.path(rc_dir, ".nextflow", "history")))
  expect_length(readLines(file.path(p, ".nextflow", "history")), 1L)
  expect_identical(RiboPilot:::nf_recheck_pending(p), 0L)
})
```

- [ ] **Step 2: Run them to see them fail**

Run:
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-pending.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-utils.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nextflow.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-recheck-run.R")'`

Expected: FAIL.
- test-nf-pending.R: 3 errors: `object 'RECHECK_QUEUE_SQL' not found`; `argument of length 0` (with a `file("")` warning: recheck.nf does not exist); `object 'nf_recheck_command' not found`.
- test-utils.R: 1 error in "grouping, export, and the run count open the project with the busy timeout": `object 'nf_recheck_pending' not found`.
- test-nextflow.R: the OOM selector expectation fails (no `|recheck_tools`), and every new recheck test fails or errors (`recheck_nf` and `main_recheck_nf` are `character(0)` with a `file("")` warning; `RECHECK_QUEUE_SQL` and `nf_recheck_command` not found). "there is ONE workflow" still passes (it only allows the new file).
- test-nf-recheck-run.R: the first test fails with `res$status` 1 (Nextflow cannot include `modules/recheck.nf`) and an error reading `seen`; the second errors with `object 'nf_recheck_pending' not found`. Without `NOT_CRAN=true` both skip ("On CRAN").

- [ ] **Step 3: Add the queue SQL and the launch helpers**

In `R/recheck.R`, append `RECHECK_QUEUE_SQL` after Task 10's last function.

Find:

```r
recheck_release_running <- function(con) {
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 1, check_note = 'recheck did not finish',
    time_stamp = ? WHERE check_switch = 4", params = list(as.integer(Sys.time())))
}
```

Replace with:

```r
recheck_release_running <- function(con) {
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 1, check_note = 'recheck did not finish',
    time_stamp = ? WHERE check_switch = 4", params = list(as.integer(Sys.time())))
}

#' The copies the recheck workflow takes
#'
#' Queued or running copies (check_switch 1 or 4) of unlocked, kept contigs on
#' samples whose assembly and taxonomy are locked, with the compute settings of
#' their annotate options and a BLAST mode: "local" when the sample searches
#' and the contig's original search was local, "map" otherwise (SPEC S7.2).
#' query_recheck in inst/nextflow/modules/recheck.nf runs this text through the
#' sqlite3 CLI, and nf_recheck_pending() counts its rows; test-nextflow.R keeps
#' the two identical. One SELECT with no trailing semicolon and no dollar sign,
#' double quote, or backslash, so it passes through the Nextflow script string
#' and the shell unchanged.
#' @noRd
RECHECK_QUEUE_SQL <- "SELECT r.ID, r.path, r.scaffold, CAST(r.rev AS TEXT) AS rev,
    COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem,
    CASE WHEN COALESCE(bo.run_blast,1) = 1 AND EXISTS (SELECT 1 FROM blast_searched bs
      WHERE bs.ID = r.ID AND bs.path = r.path AND bs.scaffold = r.scaffold AND bs.source = 'local')
    THEN 'local' ELSE 'map' END AS mode
  FROM review_units r
  JOIN annotate an ON an.ID = r.ID AND an.path = r.path AND an.scaffold = r.scaffold
  JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
  JOIN assemblies asm ON asm.ID = r.ID AND asm.path = r.path AND asm.scaffold = r.scaffold
  JOIN assemble a ON a.ID = r.ID
  JOIN samples s ON s.ID = r.ID
  LEFT JOIN blast_opts bo ON bo.blast_opts = a.blast_opts
  WHERE r.check_switch IN (1,4) AND COALESCE(an.annotate_lock,0) = 0
    AND COALESCE(asm.ignore,0) = 0 AND a.assemble_lock = 1 AND s.taxonomy_lock = 1"
```

In `R/nf_command.R`, the header comment says there is one workflow; it now names both builders.

Find:

```r
# Shared Nextflow command builder used by run_assembly() (CLI) and the app's
# Update modal (processx). Returns the argument vector that follows the
# `nextflow` binary, so both a `bash -c "cd .. && nextflow <args>"` call and
# `processx::process$new("nextflow", args, wd = path)` can use it.
#
# There is ONE workflow. Re-annotation used to be WF2 (main_annotate.nf) and is
# now an arm of WF1: dev/design/15_WF2_DECOMMISSION.md.
```

Replace with:

```r
# Nextflow command builders: nf_command() for the main pipeline, used by
# run_assembly() (CLI) and the Assemble step's Update (processx), and
# nf_recheck_command() for the recheck behind the Update on Annotation Review.
# Each returns the argument vector that follows the `nextflow` binary, so both a
# `bash -c "cd .. && nextflow <args>"` call and
# `processx::process$new("nextflow", args, wd = path)` can use it.
#
# The main pipeline is ONE workflow: re-annotation used to be WF2
# (main_annotate.nf) and is now an arm of WF1 (dev/design/15_WF2_DECOMMISSION.md).
# The recheck (main_recheck.nf) is a second workflow on purpose: it reads and
# writes only the review copies, never the originals
# (dev/specs/2026-10-08-contig-recheck-design.md, section 7).
```

Append the two helpers after `nf_pending_count`.

Find:

```r
        AND NOT (a.assemble_switch IN (1,4) AND COALESCE(a.assemble_lock,0) = 0)")$n
  as.integer(n_asm + n_ann)
}
```

Replace with:

```r
        AND NOT (a.assemble_switch IN (1,4) AND COALESCE(a.assemble_lock,0) = 0)")$n
  as.integer(n_asm + n_ann)
}

#' Build the Nextflow argument vector for the recheck workflow
#'
#' main_recheck.nf reruns the offline checks on the review copies queued on
#' Annotation Review. It is launched from `<project>/.recheck` with its own work
#' folder, and every recheck process is uncached, so there is no `-resume` and
#' the main pipeline's run history is never touched. Never stops on a blank
#' rawDir (the recheck reads no reads), never passes an empty parameter value
#' (Nextflow turns `--x ''` into the string "true"), and has no side effects: it
#' also backs the command preview.
#' @param path Project directory (contains `.sqlite`, `.config`).
#' @return character vector of args (starting with "run").
#' @noRd
nf_recheck_command <- function(path) {
  cfg <- file.path(path, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg) else character(0)
  raw_dir <- config_param(cfg_lines, "rawDir") %|NA|% ""
  if (raw_dir %in% c("NA", "NULL")) raw_dir <- ""
  blastdb <- nf_blastdb(path)
  c("run", app_sys("nextflow", "main_recheck.nf"), "-c", cfg,
    "-work-dir", file.path(path, "work_recheck"),
    "--db", file.path(path, ".sqlite"),
    "--projectDir", path,
    if (nzchar(raw_dir)) c("--rawDir", raw_dir),
    nf_pkg_args(),
    "--taxdump_dir", file.path(path, ".taxonomy"),
    if (nzchar(blastdb)) c("--blastdb", blastdb),
    container_bind_opts(cfg_lines,
      reference_mount_dirs(path,
        c(path, if (nzchar(raw_dir)) raw_dir, config_mounted_dirs(cfg_lines)))))
}

#' Count the contigs the recheck workflow would take (Annotation Review's Update)
#'
#' The rows of RECHECK_QUEUE_SQL, the selection query_recheck makes.
#' @noRd
nf_recheck_pending <- function(path) {
  con <- open_project_db(file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  as.integer(nrow(DBI::dbGetQuery(con, RECHECK_QUEUE_SQL)))
}
```

`config_param()` returns NA for an absent key (R/run_assembly.R:30-38), so `%|NA|%` gives "" there; the `"NA"`/`"NULL"` check matches nf_command's `cfg_val`. `nrow()` rather than `SELECT COUNT(*) FROM (...)`, so the constant never needs to be wrapped. Check both files parse: `Rscript -e 'invisible(parse("R/recheck.R")); invisible(parse("R/nf_command.R"))'`.

- [ ] **Step 4: Create the workflow, its module, and the bin script**

Create `inst/nextflow/main_recheck.nf`:

```groovy
// RiboPilot recheck: reruns every offline check on the review copies queued on
// Annotation Review (dev/specs/2026-10-08-contig-recheck-design.md, section 7).
// Each copy goes through the same rRNA tool script as the main pipeline's
// annotate_tools (bin/rp_rrna_tools.sh), then rp_recheck.R scores it with the
// main pipeline's check functions in review mode. No reads, no NCBI.
//
// Launched from <project>/.recheck with its own work folder and no -resume
// (nf_recheck_command() in R/nf_command.R), so the main pipeline's run history
// is never touched.

include { query_recheck; recheck_tools; recheck } from './modules/recheck.nf'

params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.blastdb = ''

workflow {
    units = query_recheck().splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.ID, r.path, r.scaffold, r.rev, r.cpus, r.mem, r.mode) }
    recheck_tools(units)
    recheck(recheck_tools.out)
}
```

Create `inst/nextflow/modules/recheck.nf`. The options read in `recheck_tools` (from `read TOOL KINGDOM ITSXT ITSXO <<<` through the `RRNAO` query) is a byte-for-byte copy of `annotate_tools`' lines in pipeline_r.nf, and the `DBARG`/`LABELS` block of `recheck` is a byte-for-byte copy of the `chimera` process's; test-nextflow.R compares both. The `rev` reaches the shell as the TSV text and is never turned into a number.

```groovy
// Processes of the recheck workflow (main_recheck.nf): rerun every offline
// check on the review copies queued on Annotation Review. All are uncached and,
// after the query, run per contig. dev/specs/2026-10-08-contig-recheck-design.md,
// sections 5.4 and 7.

// Declared here, not only in main_recheck.nf: includes are resolved before the
// including script's body runs, so a main_recheck.nf default is still null when
// this module reads it, and a null interpolates into the shell as 'null'.
params.blastdb = ''
params.taxdump_dir = ''

// Take the queued copies and mark them running (4) in one transaction, with
// the rev each copy had at that moment. An edit after this raises the rev, so
// recheck_tools and the recheck discard the result. The selection is
// RECHECK_QUEUE_SQL in R/recheck.R, word for word (test-nextflow.R checks it);
// nf_recheck_pending() counts the same rows.
process query_recheck {
    cache false
    executor 'local'
    // A read that loses a lock race must not terminate the run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3
    output:
    path 'recheck_units.tsv'
    script:
    """
    printf 'ID\tpath\tscaffold\trev\tcpus\tmem\tmode\n' > recheck_units.tsv
    sqlite3 -cmd '.timeout 10000' -separator '\t' '${params.db}' "
      BEGIN IMMEDIATE;
      CREATE TEMP TABLE q AS
      SELECT r.ID, r.path, r.scaffold, CAST(r.rev AS TEXT) AS rev,
          COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem,
          CASE WHEN COALESCE(bo.run_blast,1) = 1 AND EXISTS (SELECT 1 FROM blast_searched bs
            WHERE bs.ID = r.ID AND bs.path = r.path AND bs.scaffold = r.scaffold AND bs.source = 'local')
          THEN 'local' ELSE 'map' END AS mode
        FROM review_units r
        JOIN annotate an ON an.ID = r.ID AND an.path = r.path AND an.scaffold = r.scaffold
        JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
        JOIN assemblies asm ON asm.ID = r.ID AND asm.path = r.path AND asm.scaffold = r.scaffold
        JOIN assemble a ON a.ID = r.ID
        JOIN samples s ON s.ID = r.ID
        LEFT JOIN blast_opts bo ON bo.blast_opts = a.blast_opts
        WHERE r.check_switch IN (1,4) AND COALESCE(an.annotate_lock,0) = 0
          AND COALESCE(asm.ignore,0) = 0 AND a.assemble_lock = 1 AND s.taxonomy_lock = 1;
      SELECT * FROM q;
      UPDATE review_units SET check_switch = 4, time_stamp = CAST(strftime('%s','now') AS INTEGER)
        WHERE EXISTS (SELECT 1 FROM q WHERE q.ID = review_units.ID AND q.path = review_units.path
          AND q.scaffold = review_units.scaffold);
      COMMIT;" >> recheck_units.tsv
    """
}

// barrnap/pybarrnap, tuned ITSx, and the stock ITSx check on one copy, through
// the script annotate_tools uses, with the same options read. The fasta is
// written only while the copy still has the queued rev; otherwise nothing runs
// and nothing is emitted. A tool failure is recorded in unit.tools_status for
// the recheck; only an out-of-memory exit fails the task, so nextflow.config
// retries it with more memory.
process recheck_tools {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    cpus { nthreads as int }
    // scales on retry; see the OOM retry policy in nextflow.config
    memory { "${(mem as int) * task.attempt} GB" }
    input:
    tuple val(id), val(upath), val(uscaf), val(rev), val(nthreads), val(mem), val(mode)
    output:
    tuple val(id), val(upath), val(uscaf), val(rev), val(mode), path('unit.tools_status'),
          path('unit.barrnap.gff'), path('unit.itsx.positions.txt'), path('unit.itsx_check*'), optional: true
    script:
    """
    sqlite3 -cmd '.timeout 10000' '${params.db}' \
      "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM review_units
       WHERE ID='${id}' AND path=${upath} AND scaffold=${uscaf} AND rev=${rev}
         AND check_switch IN (1,4) AND sequence IS NOT NULL;" > unit.fasta

    read TOOL KINGDOM ITSXT ITSXO <<< "\$(sqlite3 -cmd '.timeout 10000' -separator ' ' '${params.db}' "
      SELECT COALESCE(ao.rrna_tool,'pybarrnap'), COALESCE(ao.rrna_kingdom,'euk'),
             COALESCE(ao.itsx_taxon,'all'), COALESCE(ao.itsx_opts,'')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"
    TOOL="\${TOOL:-pybarrnap}"; KINGDOM="\${KINGDOM:-euk}"; ITSXT="\${ITSXT:-all}"
    # own query: an empty itsx_opts would otherwise shift the columns under read
    RRNAO="\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "
      SELECT COALESCE(CASE WHEN ao.rrna_tool = 'barrnap' THEN ao.barrnap_opts ELSE ao.pybarrnap_opts END, '')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"

    rc=0
    if [ -s unit.fasta ]; then
        rp_rrna_tools.sh unit.fasta "\$TOOL" "\$KINGDOM" "\$ITSXT" "\$ITSXO" "\$RRNAO" ${task.cpus} || rc=\$?
    fi
    case "\$rc" in 104|134|137|139|140|143|247) exit "\$rc" ;; esac
    """
}

// Score the copy: recheck_unit() reruns the checks on the review tables and
// marks the copy done only if its rev is unchanged. The local panel and the
// labels database are resolved as the main pipeline's chimera process does.
// rp_recheck.R exits 0 for done, stale, and failed; the database holds which.
process recheck {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    input:
    tuple val(id), val(upath), val(uscaf), val(rev), val(mode), path(status), path(gff), path(pos), path(chk)
    script:
    """
    DBARG=""
    if [ -z "\${RIBOPILOT_NO_LOCAL_BLAST:-}" ]; then
        DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}' 2>/dev/null || echo "")
    fi
    LABELS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "
      SELECT COALESCE(ao.labels_db,'') FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      WHERE a.ID='${id}' LIMIT 1;" 2>/dev/null || echo "")
    RP_THREADS=${task.cpus} RIBOPILOT_SRC='${params.ribopilot_src}' rp_recheck.R \
      '${params.db}' '${id}' ${upath} ${uscaf} ${rev} ${mode} . "\$DBARG" '${params.taxdump_dir}' "\$LABELS"
    """
}
```

Notes on this file: a stale `rev` leaves `unit.fasta` empty, so the tools never run, `unit.tools_status` is never written, and the `optional: true` output emits nothing (no "Missing output file(s)" failure). When barrnap fails, `rp_rrna_tools.sh` still writes all outputs and the status file, so the copy reaches `recheck` and `recheck_unit` reads the failure from `unit.tools_status`. The input of `recheck_tools` is named `mem` because the OOM closure in nextflow.config and the Hydra `clusterOptions` read it. Empty `DBARG` and `LABELS` reach rp_recheck.R as empty positional arguments, which its `opt()` turns into NULL (never `--x ''` on the Nextflow command line).

Create `inst/nextflow/bin/rp_recheck.R` (it follows test-nf-bin-contract.R: the `rp_file` guard, `rp_args(a, 7, ...)` at the start of a line before `local({`, one indented `on.exit(DBI::dbDisconnect(con))`, and the internal `recheck_unit` only through `RiboPilot:::`):

```r
#!/usr/bin/env Rscript
# Args: <db> <id> <path> <scaffold> <rev> <mode> <tools_dir> [dbarg] [taxdump_dir] [label_db]
# Exits 0 when the recheck ran, whether the copy ended done, stale, or failed:
# the database holds the outcome.
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 7, "rp_recheck.R <db> <id> <path> <scaffold> <rev> <mode> <tools_dir> [dbarg] [taxdump_dir] [label_db]")
opt <- function(i) if (length(a) >= i && nzchar(a[i])) a[i] else NULL
r <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  RiboPilot:::recheck_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), rev = a[5], mode = a[6],
    tools_dir = a[7], dbarg = opt(8), taxdump_dir = opt(9), label_db = opt(10),
    threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "recheck_work"))
})
cat("recheck:", a[2], a[3], a[4], "-", r, "\n")
```

Make it executable (Nextflow calls bin scripts by name, and the Write tool creates 0644 files):

Run: `chmod +x inst/nextflow/bin/rp_recheck.R`

- [ ] **Step 5: Give recheck_tools the out-of-memory retry**

In `inst/nextflow/nextflow.config`:

Find:

```groovy
    withName: 'assemble|fallback_assemble|preprocess|annotate_tools|find_rrna_scan' {
```

Replace with:

```groovy
    withName: 'assemble|fallback_assemble|preprocess|annotate_tools|find_rrna_scan|recheck_tools' {
```

In `inst/config.NMNH_Hydra` (the cluster memory request; R/generate_config.R is unchanged):

Find:

```groovy
    withName: 'preprocess|assemble|fallback_assemble|annotate_tools|find_rrna_scan' {
```

Replace with:

```groovy
    withName: 'preprocess|assemble|fallback_assemble|annotate_tools|find_rrna_scan|recheck_tools' {
```

The config selector overrides the process's own `errorStrategy 'ignore'`, so exit 137 is retried with `memory` scaled by `task.attempt` and any other failure is still ignored. test-nextflow.R's "an out-of-memory kill retries with more memory" test reads this selector and checks that `recheck_tools` declares `memory { ...mem... task.attempt }`.

- [ ] **Step 6: Run the tests to see them pass**

Run:
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-pending.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-utils.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nextflow.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-bin-contract.R")'`
`NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-recheck-run.R")'`
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-no-duplicate-definitions.R")'`
`ls -l inst/nextflow/bin/rp_recheck.R`
`grep -nP '[^\x00-\x7F]' R/recheck.R R/nf_command.R inst/nextflow/main_recheck.nf inst/nextflow/modules/recheck.nf inst/nextflow/bin/rp_recheck.R tests/testthat/test-nf-recheck-run.R`

Expected: PASS. test-nf-pending.R 33; test-nextflow.R 389; test-nf-recheck-run.R 28 (it runs the real main.nf and then main_recheck.nf with the fake tools; a few minutes); test-utils.R, test-nf-bin-contract.R (now also checks rp_recheck.R), and test-no-duplicate-definitions.R pass with no failures. `rp_recheck.R` is `-rwxrwxr-x`. The ASCII grep prints nothing.

- [ ] **Step 7: Commit**

```bash
git add --chmod=+x inst/nextflow/bin/rp_recheck.R
git add inst/nextflow/main_recheck.nf inst/nextflow/modules/recheck.nf R/recheck.R R/nf_command.R \
  inst/nextflow/nextflow.config inst/config.NMNH_Hydra \
  tests/testthat/test-nextflow.R tests/testthat/test-nf-pending.R tests/testthat/test-utils.R \
  tests/testthat/test-nf-recheck-run.R
git ls-files -s inst/nextflow/bin/rp_recheck.R
git commit -m "Add recheck workflow and launch command"
```

Expected: the `git ls-files -s` line starts with `100755`, like every other file in `inst/nextflow/bin`.


### Task 12: Main pipeline follows copies; taxonomy refresh

**Files:**
- Modify: `R/review.R` (append `review_follow_original()` at the end of the file; Task 2 created it)
- Modify: `R/screen.R` (`screen_sample()` and its `@return` line, lines 383-394)
- Modify: `R/tax_search.R` (`rescreen_samples()` and its roxygen, lines 103-112)
- Modify: `R/itsx_check.R` (new `itsx_refresh_origin()` right after `itsx_origin_flag()`, which Task 5 added)
- Modify: `R/app_opts.R` (`opts_apply()`, the per-sample loop at lines 617-623)
- Modify: `inst/nextflow/bin/rp_reconcile.R`, `inst/nextflow/bin/rp_blast.R`, `inst/nextflow/bin/rp_searched.R`, `inst/nextflow/bin/rp_curate.R`, `inst/nextflow/bin/rp_chimera.R`, `inst/nextflow/bin/rp_screen.R`
- Modify: `man/screen_sample.Rd` (regenerated by roxygen, never edited by hand)
- Create: `tests/testthat/test-review-follow.R`
- Test: `tests/testthat/test-screen.R` (lines 392-397, one test appended), `tests/testthat/test-taxonomy-assign.R` (`mk_tax_db()` lines 12-13, one test appended), `tests/testthat/test-app-opts.R` (one test appended), `tests/testthat/test-blast-parse.R` (fixture at lines 376-378), `tests/testthat/test-db-dist.R` (`mk_searched_script_db()` lines 352-369)

**Interfaces:**
- Consumes:
  - Task 1: `new_db()` builds the 8 review tables and 5 views (so `chimera_test_con()` has them); `REVIEW_UNITS_DDL`, `REVIEW_TABLES`; `mk_app_project()` carries the review schema; `fixture_mark_done(con)` (tests).
  - Task 2: `with_savepoint(con, code, name = "rp_review")`, `review_has_copy(con, id, path, scaffold)`, `review_drop(con, id, path, scaffold)`, `review_invalidate(con, id, path, scaffold)`, `review_copy(con, id, path, scaffold)` (tests).
  - Task 5: `itsx_origin_flag(origin, scores, lineage)`; `annotate.itsx_origin` / `itsx_scores` stored by the stock check; `rp_reconcile.R`'s `itsx_store_not_run` handler (quoted below).
  - Task 6: `curate_one(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE)`.
  - Task 7: `screen_one(con, id, path = 1L, scaffold = 1L, force = FALSE, review = FALSE)`.
- Produces:
  - `review_follow_original(con, id, path, scaffold)` (R/review.R, `@noRd`): reads `review_units` first and returns `invisible("none")` without touching `annotate` when the contig has no copy; else one `with_savepoint()` that re-reads the copy's `edited` and the contig's `annotate_lock`: locked -> `"locked"` (nothing written); `edited = 0` -> `review_drop()`, `"dropped"`; else `review_invalidate()`, `"invalidated"`. Returns that string invisibly.
  - `screen_sample(con, id, force = FALSE, variants = TRUE)`: still returns the unit count, now with `attr(, "written")` = `data.frame(path, scaffold)` of the units whose screen row `screen_one()` rewrote (not the ones it skipped as current).
  - `rescreen_samples(con, ids)`: per sample, DELETE its screen rows and `screen_sample(con, id)` (variants on), then `tryCatch(curate_one(con, id, p, s))` for every unit with `annotate.n_features IS NOT NULL`, then `screen_one(con, id, p, s, force = TRUE, review = TRUE)` and `curate_one(con, id, p, s, review = TRUE)` for every copy at `check_switch = 2` (locked included), then `UPDATE review_units SET rev = rev + 1, time_stamp = ? WHERE ID = ? AND check_switch IN (1, 4)`; finally `itsx_refresh_origin(con, ids)`. Never calls `review_follow_original`.
  - `itsx_refresh_origin(con, ids)` (R/itsx_check.R, `@noRd`): for `annotate` rows and `review_units` rows at `check_switch = 2` of samples `ids` with a non-empty `itsx_origin` and an `itsx_warnings` that is not NULL, not `itsx_not_detected`, and not `not_run(...)`: drops the `itsx_origin(...)` token, appends `itsx_origin_flag(itsx_origin, itsx_scores, samples.lineage)` when non-empty, joins with `;`, and UPDATEs the row when the text changed. Returns the number of rows rewritten, invisibly.
  - The six per-contig bin scripts call `RiboPilot:::review_follow_original(con, a[2], <path>, <scaffold>)` once, after the contig's own writes; `rp_chimera.R` only when `chimera_unit()` returned a result; `rp_screen.R` only for the units in `attr(n, "written")`.
  - `opts_apply()` follows every copy of each target sample when the option set queues `annotate` (annotation and curation options).

The original tables keep the main pipeline's results and edits live on review copies, so when the main pipeline rewrites a contig's originals its copy has to follow (S5.3): a copy made only by State is dropped so the contig reads the new originals, an edited copy keeps its edits and goes back to "not run", and a locked copy is left alone because the recheck skips locked contigs. `review_follow_original()` does this, and the six per-contig bin scripts (`rp_reconcile.R`, `rp_blast.R`, `rp_searched.R`, `rp_curate.R`, `rp_chimera.R`, `rp_screen.R`) and `opts_apply()` call it (S5.3 table, S10). `rp_screen.R` screens a whole sample and `screen_one()` skips rows that are already current, so `screen_sample()` now reports which units it really rewrote and only those are followed; it learns that from SQLite's `total_changes()` around each `screen_one()` call, so `screen_one()`'s verdict return (Task 7's function) is untouched. A taxonomy change (S5.3 last rows, S6.3) refreshes the originals (screen, curation, ITSx origin flag) and the done copies in place, and raises `rev` on queued or running copies so an in-flight recheck result is discarded.

Tests that start `Rscript` use `skip_on_cran()`, so the run commands below set `NOT_CRAN=true`; without it they are skipped, not run.

- [ ] **Step 1: Write the failing tests for following copies and for the screen's written units**

Create `tests/testthat/test-review-follow.R`:

```r
# The main pipeline's reruns and option changes follow a contig's review copy (SPEC S5.3).

# Five annotated, done 100 bp contigs of S1, each with one feature and a
# current screen row.
rf_units <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:5) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence)
      VALUES ('S1', 1, ?, 'linear', 100, ?)", params = list(s, strrep("ACGT", 25)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        annotate_opts, curate_opts, n_features, missing, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, 'default', 'default', 1, '', '')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand, status)
      VALUES ('S1', 1, ?, '18S', 1, 90, '+', 'ok')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version)
      VALUES ('S1', 1, ?, 'on_target', ?)", params = list(s, RiboPilot:::SCREEN_VERSION))
  }
}

test_that("review_follow_original drops a State copy, invalidates an edited one, and spares a locked one", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rf_units(con)
  for (s in 2:5) RiboPilot:::review_copy(con, "S1", 1L, s)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE scaffold IN (3, 4)")
  DBI::dbExecute(con, "UPDATE review_annotations SET start = 11 WHERE scaffold IN (3, 4)")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold IN (4, 5)")
  follow <- function(s) RiboPilot:::review_follow_original(con, "S1", 1L, s)
  n_rows <- function(s) vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = ?", t),
    params = list(s))$n, integer(1))
  snap <- function(s) list(n_rows(s),
    DBI::dbGetQuery(con, "SELECT edited, check_switch, CAST(rev AS TEXT) AS rev, time_stamp
      FROM review_units WHERE scaffold = ?", params = list(s)),
    DBI::dbGetQuery(con, "SELECT start, status FROM review_annotations WHERE scaffold = ?",
      params = list(s)))
  rev3 <- DBI::dbGetQuery(con, "SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE scaffold = 3")$rev
  before <- lapply(4:5, snap)
  expect_equal(unname(n_rows(3)), c(1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L))

  expect_invisible(follow(1))
  expect_equal(follow(1), "none")
  expect_equal(sum(n_rows(1)), 0L)

  expect_equal(follow(2), "dropped")
  expect_equal(sum(n_rows(2)), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy FROM review_units_v WHERE scaffold = 2")$has_copy, 0L)

  expect_equal(follow(3), "invalidated")
  expect_equal(DBI::dbGetQuery(con, "SELECT edited, check_switch, rev = CAST(? AS INTEGER) + 1 AS bumped
    FROM review_units WHERE scaffold = 3", params = list(rev3)),
    data.frame(edited = 1L, check_switch = 1L, bumped = 1L))
  a3 <- DBI::dbGetQuery(con, "SELECT start, status FROM review_annotations WHERE scaffold = 3")
  expect_equal(a3$start, 11L)
  expect_true(is.na(a3$status))
  expect_equal(n_rows(3)[["review_screen"]], 0L)

  expect_equal(follow(4), "locked")
  expect_equal(follow(5), "locked")
  expect_equal(lapply(4:5, snap), before)

  # inside an open transaction it nests as a savepoint
  DBI::dbBegin(con)
  expect_equal(follow(3), "invalidated")
  DBI::dbRollback(con)
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("a contig without a copy is answered from review_units alone", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  expect_equal(RiboPilot:::review_follow_original(con, "S1", 1L, 1L), "none")
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("each per-contig pipeline script follows the copy of the contig it wrote", {
  bin <- system.file("nextflow", "bin", package = "RiboPilot")
  for (f in c("rp_reconcile.R", "rp_blast.R", "rp_searched.R", "rp_curate.R", "rp_chimera.R",
              "rp_screen.R")) {
    txt <- readLines(file.path(bin, f), warn = FALSE)
    expect_equal(sum(grepl("RiboPilot:::review_follow_original(con, a[2], ", txt, fixed = TRUE)), 1L,
      info = f)
  }
  chim <- readLines(file.path(bin, "rp_chimera.R"), warn = FALSE)
  expect_true(any(grepl("if (!is.null(r)) RiboPilot:::review_follow_original", chim, fixed = TRUE)))
})

test_that("rp_screen.R follows only the contigs whose screen it rewrote, and rp_curate.R its contig", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("pkgload")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!nzchar(src), "not a source checkout")
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rf_units(con)
  for (s in 1:3) RiboPilot:::review_copy(con, "S1", 1L, s)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE scaffold = 3")
  # only contigs 2 and 3 lack a current screen row, so the screen rewrites only theirs
  DBI::dbExecute(con, "UPDATE screen SET screen_version = NULL WHERE scaffold IN (2, 3)")
  run <- function(script, args = character(0)) processx::run("Rscript",
    c(file.path(system.file("nextflow", "bin", package = "RiboPilot"), script),
      DBI::dbGetInfo(con)$dbname, "S1", args),
    env = c("current", RIBOPILOT_SRC = src), error_on_status = FALSE)
  copies <- function() DBI::dbGetQuery(con,
    "SELECT scaffold, check_switch FROM review_units ORDER BY scaffold")

  res <- run("rp_screen.R")
  expect_equal(res$status, 0L, info = res$stderr)
  expect_match(res$stdout, "screen: S1 - 5 unit(s)", fixed = TRUE)
  expect_equal(copies(), data.frame(scaffold = c(1L, 3L), check_switch = c(2L, 1L)))

  res <- run("rp_curate.R", c("1", "1"))
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(copies(), data.frame(scaffold = 3L, check_switch = 1L))
})
```

In `tests/testthat/test-screen.R`, the count now carries the `written` attribute, so compare the bare number (an attribute makes `expect_equal()` fail):

Find:

```r
  n <- screen_sample(con, "S1")
  expect_equal(n, 2L)
```

Replace with:

```r
  n <- screen_sample(con, "S1")
  expect_equal(as.integer(n), 2L)
```

Append to the end of `tests/testthat/test-screen.R` (one blank line before it; it reuses `mk_screen_db()`, `add_hits()`, and `FISH_T` from that file):

```r

test_that("screen_sample names the contigs whose screen row it rewrote", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T); add_hits(con, 2, FISH_T)
  written <- function(force = FALSE) attr(screen_sample(con, "S1", force = force), "written")
  expect_equal(written(), data.frame(path = c(1L, 1L), scaffold = c(1L, 2L)))
  # both rows are current now, so neither is rewritten
  expect_equal(nrow(written()), 0L)
  DBI::dbExecute(con, "UPDATE screen SET screen_version = 0 WHERE scaffold = 2")
  expect_equal(written(), data.frame(path = 1L, scaffold = 2L))
  expect_equal(written(force = TRUE), data.frame(path = c(1L, 1L), scaffold = c(1L, 2L)))
})
```

Two existing tests run the real `rp_searched.R` against hand-built databases that have no `review_units`; once the script asks for the contig's copy they would stop with `no such table: review_units`. Give both fixtures the table (no copy, so the script answers "none" and never reads `annotate`).

In `tests/testthat/test-blast-parse.R` ("a locally searched unit is labelled local_rrna, not core_nt"):

Find:

```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,?)", params = list(FISH_T))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
```

Replace with:

```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,?)", params = list(FISH_T))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
```

In `tests/testthat/test-db-dist.R` (`mk_searched_script_db()`), two edits:

Find:

```r
# rp_searched.R's DB fixture: everything derive_blast_products() and
# mark_unit_searched() touch, plus reference_config for the provenance write.
```

Replace with:

```r
# rp_searched.R's DB fixture: everything derive_blast_products() and
# mark_unit_searched() touch, plus reference_config for the provenance write
# and review_units for the review-copy check.
```

Find:

```r
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  DBI::dbDisconnect(con)
}
```

Replace with:

```r
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbDisconnect(con)
}
```

- [ ] **Step 2: Run them to see them fail**

Run: `NOT_CRAN=true Rscript -e 'options(testthat.progress.max_fails = Inf); devtools::load_all("."); for (f in c("review-follow", "screen")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: FAIL.
- test-review-follow.R `[ FAIL 11 | WARN 0 | SKIP 0 | PASS 4 ]`: the first two tests error with `object 'review_follow_original' not found`; the static test fails 7 times (`sum(...)` is 0, expected 1, once per script, and the `rp_chimera.R` guard is missing); the script test fails twice (all three copies still at `check_switch` 2 after `rp_screen.R`, and copy 1 still there after `rp_curate.R`).
- test-screen.R `FAIL 4`: in "screen_sample names the contigs whose screen row it rewrote", every `written()` is NULL. The other tests (Task 7's included) pass.

test-blast-parse.R and test-db-dist.R pass both before and after this task's code; their fixture edits are what keep them green after Step 3.

- [ ] **Step 3: Add review_follow_original(), the written units, and the calls in the bin scripts**

Append to the end of `R/review.R` (one blank line before it):

```r

#' Follow a main-pipeline rerun on a contig's review copy (S5.3)
#'
#' The pipeline's per-contig scripts call this after they rewrite a contig's
#' originals, and option changes call it at once. An unlocked copy made only
#' by State (edited 0) is deleted, so the contig reads its new originals; an
#' unlocked edited copy keeps its edits and its checks go back to not run; a
#' locked copy is left alone, since the recheck skips locked contigs.
#' @return invisible "none" (no copy), "locked", "dropped", or "invalidated".
#' @noRd
review_follow_original <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  if (!review_has_copy(con, id, path, scaffold)) return(invisible("none"))
  invisible(with_savepoint(con, {
    u <- DBI::dbGetQuery(con,
      "SELECT r.edited, COALESCE(an.annotate_lock, 0) AS locked FROM review_units r
       LEFT JOIN annotate an ON an.ID = r.ID AND an.path = r.path AND an.scaffold = r.scaffold
       WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?", params = list(id, path, scaffold))
    if (nrow(u) == 0) {
      "none"
    } else if (u$locked[1] == 1) {
      "locked"
    } else if (u$edited[1] == 0) {
      review_drop(con, id, path, scaffold)
      "dropped"
    } else {
      review_invalidate(con, id, path, scaffold)
      "invalidated"
    }
  }))
}
```

The copy check outside the savepoint is a plain read, so a contig without a copy (nearly every contig the pipeline writes) never takes the write lock. Inside the savepoint the state is read again, because another process may have restored or locked the contig in between; `review_drop()` and `review_invalidate()` each open their own savepoint inside this one.

In `R/screen.R`, `screen_sample()` records which units `screen_one()` rewrote. Find:

```r
#' @return (invisibly) the number of units screened.
#' @export
screen_sample <- function(con, id, force = FALSE, variants = TRUE) {
  u <- DBI::dbGetQuery(con,
    "SELECT path, scaffold FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
    params = list(id))
  for (i in seq_len(nrow(u))) screen_one(con, id, u$path[i], u$scaffold[i], force = force)
  # after the loop, never before: screen_one() rewrites a row with DELETE plus
  # append, so a variant column written first would be silently wiped
  if (isTRUE(variants)) screen_variants(con, id)
  invisible(nrow(u))
}
```

Replace with:

```r
#' @return (invisibly) the number of units screened, with attribute "written":
#'   data.frame(path, scaffold) of the units whose screen row was rewritten
#'   rather than skipped as current.
#' @export
screen_sample <- function(con, id, force = FALSE, variants = TRUE) {
  u <- DBI::dbGetQuery(con,
    "SELECT path, scaffold FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
    params = list(id))
  # screen_one() writes nothing for a row it skips
  n_changes <- function() DBI::dbGetQuery(con, "SELECT total_changes() AS n")$n
  wrote <- logical(nrow(u))
  for (i in seq_len(nrow(u))) {
    n0 <- n_changes()
    screen_one(con, id, u$path[i], u$scaffold[i], force = force)
    wrote[i] <- n_changes() > n0
  }
  # after the loop, never before: screen_one() rewrites a row with DELETE plus
  # append, so a variant column written first would be silently wiped
  if (isTRUE(variants)) screen_variants(con, id)
  written <- data.frame(path = u$path[wrote], scaffold = u$scaffold[wrote])
  invisible(structure(nrow(u), written = written))
}
```

`total_changes()` counts the rows this connection has written; `screen_one()` either returns before any write or rewrites the row (DELETE plus append), so the count moves exactly when it wrote.

In `inst/nextflow/bin/rp_reconcile.R` (this is Task 5's handler), follow after the annotation and the ITSx check. Find:

```r
        error = function(e2) message("ITSx check result not stored: ", conditionMessage(e2)))
    })
})
```

Replace with:

```r
        error = function(e2) message("ITSx check result not stored: ", conditionMessage(e2)))
    })
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
})
```

In `inst/nextflow/bin/rp_blast.R`, find:

```r
    hits_file = a[5], taxdump_dir = taxdump_dir)
})
```

Replace with:

```r
    hits_file = a[5], taxdump_dir = taxdump_dir)
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
})
```

In `inst/nextflow/bin/rp_searched.R`, find:

```r
  derive_blast_products(con, a[2], as.integer(a[3]), as.integer(a[4]))
})
```

Replace with:

```r
  derive_blast_products(con, a[2], as.integer(a[3]), as.integer(a[4]))
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
})
```

In `inst/nextflow/bin/rp_curate.R` (the last line prints `res$n_flags`, so the block still returns curate_one's result), find:

```r
  curate_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
})
```

Replace with:

```r
  res <- curate_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
  res
})
```

In `inst/nextflow/bin/rp_chimera.R` (`chimera_unit()` returns NULL when it skipped the contig and wrote nothing), find:

```r
  chimera_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), dbarg = opt(5),
    taxdump_dir = opt(6), label_db = opt(7), threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "chimera_work"))
})
```

Replace with:

```r
  r <- chimera_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), dbarg = opt(5),
    taxdump_dir = opt(6), label_db = opt(7), threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "chimera_work"))
  if (!is.null(r)) RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
  r
})
```

In `inst/nextflow/bin/rp_screen.R` (keep the literal `variants = FALSE`; test-screen-variants.R:452 greps for it), find:

```r
  screen_sample(con, a[2], variants = FALSE)
})
```

Replace with:

```r
  n <- screen_sample(con, a[2], variants = FALSE)
  w <- attr(n, "written")
  for (i in seq_len(nrow(w))) RiboPilot:::review_follow_original(con, a[2], w$path[i], w$scaffold[i])
  n
})
```

`review_follow_original` is an internal, so every bin script writes it `RiboPilot:::` (test-nf-bin-contract.R:119-130 fails a bare call to an unexported function). The calls are not wrapped in tryCatch: if one fails, the task fails and Nextflow retries it, rather than leaving a done copy that no longer matches its original.

Check the files parse: `for f in R/review.R R/screen.R inst/nextflow/bin/rp_reconcile.R inst/nextflow/bin/rp_blast.R inst/nextflow/bin/rp_searched.R inst/nextflow/bin/rp_curate.R inst/nextflow/bin/rp_chimera.R inst/nextflow/bin/rp_screen.R; do Rscript -e "invisible(parse('$f'))"; done`

- [ ] **Step 4: Run the tests to see them pass**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("review-follow", "screen", "blast-parse", "db-dist", "nf-bin-contract", "screen-variants")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: every file `FAIL 0 | WARN 0 | SKIP 0`. test-review-follow.R `PASS 31`; test-db-dist.R `PASS 83`; test-nf-bin-contract.R `PASS 162`; test-screen-variants.R `PASS 132`; test-screen.R and test-blast-parse.R include Task 7's tests in their counts (without them: 120 and 68).

- [ ] **Step 5: Write the failing taxonomy-refresh test**

In `tests/testthat/test-taxonomy-assign.R`, `rescreen_samples()` now reads `annotate` and `review_units`, which `mk_tax_db()` lacks; the existing re-screen test (:49-59) would stop with `no such table`. Give the fixture both (empty `annotate`, so nothing is curated there).

Find:

```r
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon, taxonomy_switch, taxonomy_lock)
```

Replace with:

```r
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  # the re-screen also reads annotate and the review copies
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, n_features INT,
    itsx_warnings TEXT, itsx_origin TEXT, itsx_scores TEXT)")
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon, taxonomy_switch, taxonomy_lock)
```

Append to the end of `tests/testthat/test-taxonomy-assign.R` (one blank line before it). The lineage strings are the ones the test taxdump gives taxids 3497 (the plant Morus) and 37577 (the bird Morus).

```r

PLANT_LIN <- "Viridiplantae;Streptophyta;Magnoliopsida;Rosales;Moraceae;Morus;Morus alba"
PLANT_T <- "33090;35493;3398;3744;3487;3497;3498"
BIRD_LIN <- "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;unclassified Morus species"
BIRD_T <- "33208;7711;8782;30410;30446;37577;"
ITSX_SC <- "M: 4 2.9e-16 49.9, T: 4 1e-08 33.15"

# S1 has just moved from the plant Morus to the bird Morus. Its four contigs
# and S2's one hold plant hits and were checked under the plant lineage;
# S1.1.4 was never annotated. Copies: S1.1.1 done and locked, holding bird
# hits; S1.1.2 not run; S1.1.3 running; S1.1.4 on hold; S2.1.1 not run.
rescreen_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE samples SET lineage = ?, lineage_taxids = ? WHERE ID = 'S1'",
    params = list(BIRD_LIN, BIRD_T))
  units <- data.frame(ID = c("S1", "S1", "S1", "S1", "S2"), scaffold = c(1:4, 1L))
  for (i in seq_len(nrow(units))) {
    p <- list(units$ID[i], units$scaffold[i])
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence)
      VALUES (?, 1, ?, 'linear', 100, ?)", params = c(p, strrep("ACGT", 25)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        annotate_opts, curate_opts, n_features, n_flags, curate_flags, itsx_warnings, itsx_origin,
        itsx_scores)
      VALUES (?, 1, ?, 2, 0, 'default', 'default', 1, 99, 'stale', '', 'T', ?)", params = c(p, ITSX_SC))
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length, strand)
      VALUES (?, 1, ?, '18S', 1, 90, 90, '+')", params = p)
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, lineage, lineage_taxids,
      is_nuclear_rrna) VALUES (?, 1, ?, ?, ?, 1)", params = c(p, PLANT_LIN, PLANT_T))
    for (k in 1:3) DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp,
        lineage, lineage_taxids, bitscore, hit_class) VALUES (?, 1, ?, ?, 1, ?, ?, ?, 'rdna')",
      params = c(p, k, PLANT_LIN, PLANT_T, 1000 - k))
  }
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M',
    itsx_warnings = 'itsx_chimeric;itsx_origin(score: Metazoa 49.9 vs Tracheophyta 33.15)'
    WHERE ID = 'S1' AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'not_run(ITSx did not finish)'
    WHERE ID = 'S1' AND scaffold = 3")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL, itsx_origin = NULL,
    itsx_warnings = 'itsx_origin(Tracheophyta)' WHERE ID = 'S1' AND scaffold = 4")

  for (i in seq_len(nrow(units))) RiboPilot:::review_copy(con, units$ID[i], 1L, units$scaffold[i])
  DBI::dbExecute(con, "UPDATE review_units SET rev = 1000, time_stamp = 1, n_flags = 99,
    curate_flags = 'stale', itsx_warnings = '', itsx_origin = 'T', itsx_scores = ?, check_switch = 1",
    params = list(ITSX_SC))
  for (s in 1:4) {
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = ? WHERE ID = 'S1' AND scaffold = ?",
      params = list(c(2L, 1L, 4L, 0L)[s], s))
  }
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE ID = 'S1' AND scaffold = 1")
  DBI::dbExecute(con, "DELETE FROM review_blast_hits WHERE ID = 'S1' AND scaffold = 1")
  for (k in 1:3) DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank,
      hsp, lineage, lineage_taxids, bitscore, hit_class) VALUES ('S1', 1, 1, ?, 1, ?, ?, ?, 'rdna')",
    params = list(k, "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;Morus bassanus",
      paste0(BIRD_T, "37578"), 1000 - k))
}

test_that("a taxon change re-curates the contigs, refreshes done copies, and requeues running ones", {
  con <- chimera_test_con(c("S1", "S2"))
  on.exit(DBI::dbDisconnect(con))
  rescreen_seed(con)
  RiboPilot:::rescreen_samples(con, "S1")
  q <- function(sql) DBI::dbGetQuery(con, sql)

  # originals: screened against the new lineage, curated when annotated
  expect_equal(q("SELECT screen_verdict FROM screen ORDER BY ID, scaffold")$screen_verdict,
    rep("off_target", 4))
  an <- q("SELECT ID, scaffold, curate_flags, itsx_warnings FROM annotate ORDER BY ID, scaffold")
  expect_equal(an$curate_flags[4:5], c("stale", "stale"))
  for (i in 1:3) {
    expect_true("contaminant" %in% strsplit(an$curate_flags[i], "; ")[[1]], info = an$curate_flags[i])
  }
  expect_equal(an$itsx_warnings, c("itsx_chimeric",
    "itsx_origin(score: Tracheophyta 33.15 vs Metazoa 49.9)", "not_run(ITSx did not finish)",
    "itsx_origin(Tracheophyta)", ""))

  # copies: the done one re-scored from its own hits; queued and running ones requeued
  ru <- q("SELECT ID, scaffold, check_switch, CAST(rev AS TEXT) AS rev, time_stamp, curate_flags,
    itsx_warnings FROM review_units ORDER BY ID, scaffold")
  expect_equal(ru$check_switch, c(2L, 1L, 4L, 0L, 1L))
  expect_equal(ru$rev, c("1000", "1001", "1001", "1000", "1000"))
  expect_equal(ru$time_stamp > 1, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_false(identical(ru$curate_flags[1], "stale"))
  expect_equal(ru$curate_flags[2:5], rep("stale", 4))
  expect_equal(ru$itsx_warnings, c("itsx_origin(score: Tracheophyta 33.15 vs Metazoa 49.9)",
    "", "", "", ""))
  expect_equal(q("SELECT ID, scaffold, screen_verdict FROM review_screen"),
    data.frame(ID = "S1", scaffold = 1L, screen_verdict = "on_target"))
})
```

- [ ] **Step 6: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-taxonomy-assign.R")'`
Expected: FAIL, `[ FAIL 9 | WARN 0 | SKIP 0 | PASS 21 ]`, all in "a taxon change re-curates the contigs, refreshes done copies, and requeues running ones": `"contaminant" %in% ...` is FALSE three times (curate_flags still `stale`); the originals' `itsx_warnings` keep their plant-lineage origin flags; `ru$rev` is all `"1000"`; no copy's `time_stamp` moved; the done copy's `curate_flags` is still `stale` and its `itsx_warnings` is `""`; `review_screen` has no row. The originals' screen verdicts already pass, and the other 6 tests pass.

- [ ] **Step 7: Refresh the originals, the done copies, and the ITSx origin on a taxonomy change**

In `R/itsx_check.R`, add `itsx_refresh_origin()` right after `itsx_origin_flag()` (Task 5's code). Find:

```r
  sprintf("itsx_origin(%s)", if (nzchar(sc)) sc else unname(ITSX_ORIGINS[origin]) %|NA|% origin)
}
```

Replace with:

```r
  sprintf("itsx_origin(%s)", if (nzchar(sc)) sc else unname(ITSX_ORIGINS[origin]) %|NA|% origin)
}

#' Judge the stored ITSx origin again after a sample's lineage changed (S6.3)
#'
#' For the originals and the done review copies (locked ones too) of samples
#' `ids`, replaces the itsx_origin flag in itsx_warnings with the one the
#' stored origin code and scores give under the sample's current lineage.
#' Rows without a stored origin (checked before it was stored), never
#' checked, not run, or with no ITS detected are left alone.
#' @return (invisibly) the number of rows rewritten.
#' @noRd
itsx_refresh_origin <- function(con, ids) {
  n <- 0L
  for (tbl in c("annotate", "review_units")) {
    for (id in ids) {
      d <- DBI::dbGetQuery(con, sprintf(
        "SELECT x.path, x.scaffold, x.itsx_warnings, x.itsx_origin, x.itsx_scores, s.lineage
         FROM %s x JOIN samples s ON s.ID = x.ID
         WHERE x.ID = ? AND x.itsx_origin IS NOT NULL AND x.itsx_origin != ''
           AND x.itsx_warnings IS NOT NULL AND x.itsx_warnings != 'itsx_not_detected'
           AND x.itsx_warnings NOT LIKE 'not_run(%%'%s", tbl,
        if (tbl == "review_units") " AND x.check_switch = 2" else ""), params = list(id))
      for (i in seq_len(nrow(d))) {
        tok <- strsplit(d$itsx_warnings[i], ";", fixed = TRUE)[[1]]
        tok <- c(tok[nzchar(tok) & !startsWith(tok, "itsx_origin(")],
          itsx_origin_flag(d$itsx_origin[i], d$itsx_scores[i], d$lineage[i]))
        new <- paste(tok[nzchar(tok)], collapse = ";")
        if (identical(new, d$itsx_warnings[i])) next
        DBI::dbExecute(con, sprintf(
          "UPDATE %s SET itsx_warnings = ? WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
          params = list(new, id, d$path[i], d$scaffold[i]))
        n <- n + 1L
      }
    }
  }
  invisible(n)
}
```

`itsx_check_flags()` joins its flags with `;` and puts the origin flag last, so dropping the old token and appending the new one gives the text a fresh check would store. `%%` is a literal `%` inside `sprintf()`.

In `R/tax_search.R`, find:

```r
#' Recompute the contamination screen of samples whose taxon changed
#' @noRd
rescreen_samples <- function(con, ids) {
  ids <- intersect(ids, DBI::dbGetQuery(con, "SELECT DISTINCT ID FROM assemblies")$ID)
  for (id in ids) {
    DBI::dbExecute(con, "DELETE FROM screen WHERE ID = ?", params = list(id))
    screen_sample(con, id)
  }
  invisible(length(ids))
}
```

Replace with:

```r
#' Recompute the lineage-dependent checks of samples whose taxon changed
#'
#' Per sample: the screen of every contig (then the variant facts), curation of
#' every annotated contig, then the screen and curation of each review copy
#' whose checks are done (locked ones too). Copies queued or running get rev
#' raised, so a recheck already under way is discarded and the contig reruns
#' with the new lineage. Last, the stored ITSx origin is judged again.
#' @noRd
rescreen_samples <- function(con, ids) {
  ids <- intersect(ids, DBI::dbGetQuery(con, "SELECT DISTINCT ID FROM assemblies")$ID)
  for (id in ids) {
    DBI::dbExecute(con, "DELETE FROM screen WHERE ID = ?", params = list(id))
    screen_sample(con, id)
    u <- DBI::dbGetQuery(con,
      "SELECT path, scaffold FROM annotate WHERE ID = ? AND n_features IS NOT NULL", params = list(id))
    for (i in seq_len(nrow(u))) {
      tryCatch(curate_one(con, id, u$path[i], u$scaffold[i]), error = function(e) NULL)
    }
    cp <- DBI::dbGetQuery(con,
      "SELECT path, scaffold FROM review_units WHERE ID = ? AND check_switch = 2", params = list(id))
    for (i in seq_len(nrow(cp))) {
      screen_one(con, id, cp$path[i], cp$scaffold[i], force = TRUE, review = TRUE)
      curate_one(con, id, cp$path[i], cp$scaffold[i], review = TRUE)
    }
    DBI::dbExecute(con,
      "UPDATE review_units SET rev = rev + 1, time_stamp = ? WHERE ID = ? AND check_switch IN (1, 4)",
      params = list(as.integer(Sys.time()), id))
  }
  itsx_refresh_origin(con, ids)
  invisible(length(ids))
}
```

The originals are refreshed first, because the copy's screen carries the original screen row's variant facts, which `screen_sample()` has just rewritten. A never-annotated contig (`n_features` NULL) is not curated: curation would give it a full set of "missing" flags. `rev` is only ever changed in SQL (it does not fit an R integer). Its callers (`resolve_samples_taxonomy`, `assign_samples_taxid`, and the app at R/app_server.R:1151) run it outside any transaction, and every function it calls opens its own.

Check the files parse: `for f in R/itsx_check.R R/tax_search.R; do Rscript -e "invisible(parse('$f'))"; done`

- [ ] **Step 8: Run the test to see it pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-taxonomy-assign.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 30 ]` (7 tests).

- [ ] **Step 9: Write the failing option-change test**

Append to the end of `tests/testthat/test-app-opts.R` (one blank line before it):

```r

test_that("saving annotation options follows the sample's review copies at once", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  fixture_mark_done(con)
  for (s in 1:2) RiboPilot:::review_copy(con, "S1", 1L, s)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE scaffold = 2")
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$flushReact()
    session$setInputs(set_annotate_opts = 1)
    session$setInputs(no_name = "fresh", no_edit = TRUE, no_cpus = 2, no_mem = 4,
      no_tool = "pybarrnap", no_kingdom = "euk", no_itsxtaxon = "all", no_itsxopts = "")
    session$setInputs(update_no = 1)
    # the State copy is gone; the edited copy keeps its edits and waits for Update
    expect_equal(RiboPilot:::app_query(
      "SELECT path, scaffold, edited, check_switch FROM review_units ORDER BY path, scaffold"),
      data.frame(path = 1L, scaffold = 2L, edited = 1L, check_switch = 1L))
  })
})
```

- [ ] **Step 10: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-opts.R")'`
Expected: FAIL, `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 63 ]`: in "saving annotation options follows the sample's review copies at once" the query still returns both copies, each at `check_switch` 2.

- [ ] **Step 11: Follow the copies in opts_apply()**

In `R/app_opts.R` (`opts_apply()`), only the annotation and curation option sets queue `annotate`; their samples' copies follow at once. The calls sit inside the same `app_write()` block (a tryCatch, not a transaction; `rows_upsert()` above has already committed), so a failure shows the usual "Not saved" alert. Find:

```r
      DBI::dbExecute(con, sprintf("UPDATE %s SET %s = 1 WHERE ID = ?", sp$queue_tbl, sp$queue),
        params = list(id))
    }
  })
```

Replace with:

```r
      DBI::dbExecute(con, sprintf("UPDATE %s SET %s = 1 WHERE ID = ?", sp$queue_tbl, sp$queue),
        params = list(id))
      if (sp$queue_tbl == "annotate") {
        cp <- DBI::dbGetQuery(con, "SELECT path, scaffold FROM review_units WHERE ID = ?",
          params = list(id))
        for (i in seq_len(nrow(cp))) review_follow_original(con, id, cp$path[i], cp$scaffold[i])
      }
    }
  })
```

Check: `Rscript -e 'invisible(parse("R/app_opts.R"))'`

- [ ] **Step 12: Run it to see it pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-opts.R")'`
Expected: PASS, `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 64 ]`.

- [ ] **Step 13: Regenerate the help page of screen_sample**

Its `@return` changed and it is exported. Never run `devtools::document()` (NAMESPACE is hand-maintained).

Run: `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`
Then: `git status --short man/`
Expected: only ` M man/screen_sample.Rd`, whose `\value{}` now reads `(invisibly) the number of units screened, with attribute "written": data.frame(path, scaffold) of the units whose screen row was rewritten rather than skipped as current.` (over three lines). If any other man/ file changed, it is not this task's: restore it with `git checkout -- man/<that file>`.

- [ ] **Step 14: Run every test file this task touches or could break**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("review-follow", "screen", "screen-variants", "taxonomy-assign", "app-opts", "blast-parse", "db-dist", "nf-bin-contract", "pipeline-drivers", "reconcile", "itsx-check", "curate", "write-atomicity", "docs-site", "no-duplicate-definitions")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: every file reports `FAIL 0`. The only skip is test-itsx-check.R's real-ITSx test (no RIBOPILOT_TEST_PROJECT). test-pipeline-drivers.R and test-reconcile.R run the real `rp_reconcile.R` and `rp_curate.R` on `new_db()` projects without copies, so the follow answers "none" there.

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nf-main-run.R")'`
Expected: PASS (4 tests). It runs the real main.nf with the fake tools, so all six changed bin scripts run on a fresh project; it takes about 4 minutes.

- [ ] **Step 15: Commit**

```bash
git add R/review.R R/screen.R R/tax_search.R R/itsx_check.R R/app_opts.R man/screen_sample.Rd inst/nextflow/bin/rp_reconcile.R inst/nextflow/bin/rp_blast.R inst/nextflow/bin/rp_searched.R inst/nextflow/bin/rp_curate.R inst/nextflow/bin/rp_chimera.R inst/nextflow/bin/rp_screen.R tests/testthat/test-review-follow.R tests/testthat/test-screen.R tests/testthat/test-taxonomy-assign.R tests/testthat/test-app-opts.R tests/testthat/test-blast-parse.R tests/testthat/test-db-dist.R
git commit -m "Main pipeline follows review copies; taxonomy refresh"
```


### Task 13: Recheck launch from Annotation Review

**Files:**
- Modify: `R/run_reports.R` (new `run_failed_rechecks()` after `run_failed_samples()` at :203-211; `nf_run_report()` at :249-250, :261-262, :267, :289)
- Modify: `R/app_run_reports.R` (Workflow label, :58)
- Modify: `R/app_modals.R` (`run_modal()`, :514-568)
- Modify: `R/app_run.R` (new `build_recheck_submit_script()` after `build_reference_submit_script()` at :159-169; `nf_run_window()` doc comment, :228 and :236)
- Modify: `R/nf_command.R` (new `nf_recheck_waiting()` at the end of the file)
- Modify: `R/app_ui.R` (`annot_ctrls`, :89-98; `rr_open` tooltip, :110)
- Modify: `R/app_server.R` (recheck block inserted after the pipeline `nf_run_window()` call at :1170-1179, before `register_app_lifecycle(session)` at :1183)
- Test: `tests/testthat/test-nf-run.R`, `tests/testthat/test-run-reports.R`, `tests/testthat/test-app-headless.R`, `tests/testthat/test-nf-command-quoting.R`, `tests/testthat/test-hydra-setup.R`, `tests/testthat/test-app-run.R`, `tests/testthat/test-app-opts.R` (:129-135)

**Interfaces:**
- Consumes: `nf_recheck_command(path)` (Task 11; argument vector after `nextflow`, with `-work-dir <project>/work_recheck`, `--projectDir <project>`, no `-resume`, no side effects) and `nf_recheck_pending(path)` (Task 11; integer, rows of `RECHECK_QUEUE_SQL`), both in `R/nf_command.R`; `recheck_release_running(con)` (Task 10, `R/recheck.R`; copies at 4 go to 1 with `check_note = 'recheck did not finish'` and a new `time_stamp`); `REVIEW_UNITS_DDL` and the `review_units` table that `new_db()` creates (Task 1).
- Produces:
  - `nf_recheck_waiting(path)` -> integer (`R/nf_command.R`)
  - `run_modal(n, headless = FALSE, sched = NULL, script = "", pfx = "", title = NULL, resume = TRUE)`; spinner id `paste0(pfx, "rp_gears")`
  - `build_recheck_submit_script(project, sched = detect_scheduler(), base = job_basename("recheck"))`
  - `run_failed_rechecks(con, since, until = NA)` -> `data.frame(ID = "<ID>.<path>.<scaffold>", note)`
  - App ids: input `annot_run_btn`; recheck window inputs `rc_nf_start`, `rc_nf_stop`, `rc_nf_close`, `rc_nf_close_confirm`, `rc_submit_job`, `rc_save_script`, `rc_submit_script_text`; outputs `rc_nf_cmd_text`, `rc_nf_prog_all`; panel `rc_nf_progress`; spinner `rc_rp_gears`
  - Run workflow name `recheck`: files `.runs/{jobs,logs,nextflow,reports}/recheck_<timestamp>*`, Run Reports label "Recheck"

Annotation Review gets its own Update button that launches the recheck workflow the way the Assemble step's Update launches the main pipeline (S7.4): run lock check, "Nothing to update" with the number of contigs waiting for the Assemble step's Update, a run window with `rc_` ids, no Resume box, Run from App from `<project>/.recheck`, or a job script on a cluster host. The run report lists the contigs the run left failed (3) or still running (4) and says to press Update on Annotation Review when stopped (S7.4). After an in-app run ends, `on_done` puts contigs still at 4 back in the queue with the note "recheck did not finish" (S5.4); `nf_run_window()` writes the report before `on_done`, so they are reported as failures first. The Assemble step's Update count, observer, and window stay exactly as they are; the recheck observers are separate.

- [ ] **Step 1: Write the failing run report tests**

Append to the end of `tests/testthat/test-nf-run.R`:

```r

test_that("a recheck report lists the contigs the run left failed or running", {
  d <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, rev, check_switch,
    check_note, time_stamp) VALUES
    ('S1', 1, 10, 1, 4, NULL, 200), ('S1', 1, 2, 1, 3, 'chimera check: no hits', 200),
    ('OLD', 1, 1, 1, 3, 'an earlier run', 50), ('DONE', 1, 1, 1, 2, NULL, 200),
    ('QUEUED', 1, 1, 1, 1, NULL, 200)")
  f <- RiboPilot:::run_failed_rechecks(con, 100)
  expect_equal(f$ID, c("S1.1.2", "S1.1.10"))
  expect_equal(f$note, c("chimera check: no hits", "recheck did not finish"))
  expect_length(RiboPilot:::run_failed_rechecks(con, 100, until = 150)$ID, 0)
  DBI::dbDisconnect(con)
  base <- RiboPilot:::job_basename("recheck", as.POSIXct(100, origin = "1970-01-01"))
  log <- file.path(d, paste0(base, ".nextflow.log"))
  r <- RiboPilot:::nf_run_report(d, log, 0L)
  expect_equal(r$type, "warning")
  expect_match(r$lead, "Run finished, but 2 contig(s) failed.", fixed = TRUE)
  expect_match(r$text, "Failed contigs (contig | note):", fixed = TRUE)
  expect_match(r$text, "S1.1.10 | recheck did not finish", fixed = TRUE)
  expect_false(grepl("Failed samples", r$text, fixed = TRUE))
  expect_false(grepl("assembler", r$text, fixed = TRUE))
  stopped <- RiboPilot:::nf_run_report(d, log, 1L, how = "stopped")
  expect_equal(stopped$lead,
    "Run stopped. Nextflow stopped its tasks. Press Update on Annotation Review.")
})
```

Append to the end of `tests/testthat/test-run-reports.R`:

```r

test_that("Run Reports names a recheck run as such", {
  d <- withr::local_tempdir()
  put_log(d, "recheck_2026-09-08_12-51-00")
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    session$setInputs(rr_open = 1)
    expect_equal(h$runs()$workflow, "recheck")
    expect_match(output$rr_tbl, "Recheck", fixed = TRUE)
    expect_false(grepl("Pipeline", output$rr_tbl, fixed = TRUE))
  })
})
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-nf-run.R", "test-run-reports.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: FAIL. test-nf-run.R: `object 'run_failed_rechecks' not found`. test-run-reports.R: `Expected output$rr_tbl to match string "Recheck"` (the table says "Pipeline"), and the `grepl("Pipeline", ...)` expectation is TRUE.

- [ ] **Step 3: Report recheck failures by contig**

In `R/run_reports.R`, find:

```r
    params = list(floor(as.numeric(since)), until)), error = function(e) empty)
}
```

Replace with:

```r
    params = list(floor(as.numeric(since)), until)), error = function(e) empty)
}

#' Contigs a recheck run left failed or still running, with their notes
#'
#' A contig still running has no note yet: the app writes "recheck did not
#' finish" only after the report.
#' @param since,until Epoch seconds (or POSIXct) bounding the run.
#' @return data.frame(ID, note); ID is `<ID>.<path>.<scaffold>`.
#' @noRd
run_failed_rechecks <- function(con, since, until = NA) {
  empty <- data.frame(ID = character(0), note = character(0), stringsAsFactors = FALSE)
  if (is.null(con) || length(since) == 0 || is.na(since)) return(empty)
  until <- if (is.na(until)) 4e9 else as.numeric(until) + 1
  tryCatch(DBI::dbGetQuery(con,
    "SELECT r.ID || '.' || r.path || '.' || r.scaffold AS ID,
            COALESCE(r.check_note, 'recheck did not finish') AS note
       FROM review_units r
      WHERE r.check_switch IN (3, 4) AND r.time_stamp >= ? AND r.time_stamp <= ?
      ORDER BY r.ID, r.path, r.scaffold",
    params = list(floor(as.numeric(since)), until)), error = function(e) empty)
}
```

Find:

```r
  bad <- if (workflow != "reference") run_failed_samples(con, since, f$ended) else
    run_failed_samples(NULL, NA)
```

Replace with:

```r
  recheck <- workflow == "recheck"
  bad <- if (recheck) run_failed_rechecks(con, since, f$ended) else
    if (workflow != "reference") run_failed_samples(con, since, f$ended) else
    run_failed_samples(NULL, NA)
```

Find:

```r
        lead = paste("Run stopped. Nextflow stopped its tasks.", if (workflow != "reference")
          "Press Update with Resume ticked to continue where it left off." else "Press Build again to restart it.")),
```

Replace with:

```r
        lead = paste("Run stopped. Nextflow stopped its tasks.", switch(workflow,
          reference = "Press Build again to restart it.",
          recheck = "Press Update on Annotation Review.",
          "Press Update with Resume ticked to continue where it left off."))),
```

Find:

```r
        paste(c(if (n_t) sprintf("%d task(s)", n_t), if (n_s) sprintf("%d sample(s)", n_s)),
```

Replace with:

```r
        paste(c(if (n_t) sprintf("%d task(s)", n_t),
          if (n_s) sprintf(if (recheck) "%d contig(s)" else "%d sample(s)", n_s)),
```

Find:

```r
  samples <- if (n_s == 0) character(0) else c("",
```

Replace with:

```r
  samples <- if (n_s == 0) character(0) else if (recheck) c("",
    "Failed contigs (contig | note):", sprintf("%s | %s", bad$ID, bad$note)) else c("",
```

In `R/app_run_reports.R`, find:

```r
      Workflow = ifelse(df$workflow == "reference", "Reference build", "Pipeline"),
```

Replace with:

```r
      Workflow = ifelse(df$workflow == "reference", "Reference build",
        ifelse(df$workflow == "recheck", "Recheck", "Pipeline")),
```

- [ ] **Step 4: Run them to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-nf-run.R", "test-run-reports.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: PASS, 0 failures (test-nf-run.R 14 tests, two of them skipped by `skip_on_cran()` under `test_file`; test-run-reports.R 10 tests). The pipeline tests at test-nf-run.R:83-127 (stopped lead says "Resume", sample failures carry the "Too few reads" hint) still pass.

- [ ] **Step 5: Write the failing run window test**

Append to the end of `tests/testthat/test-app-headless.R`:

```r

test_that("the recheck window prefixes every id and has no Resume box", {
  has_id <- function(html, id) grepl(sprintf('id="%s"', id), html, fixed = TRUE)
  local_ids <- c("nf_cmd_text", "nf_progress", "nf_prog_all", "nf_start", "nf_stop",
    "nf_close", "rp_gears")
  html <- as.character(RiboPilot:::run_modal(2, pfx = "rc_",
    title = "RiboPilot recheck - 2 contig(s) to update", resume = FALSE))
  for (id in local_ids) {
    expect_true(has_id(html, paste0("rc_", id)), info = id)
    expect_false(has_id(html, id), info = id)
  }
  expect_true(grepl("RiboPilot recheck - 2 contig(s) to update", html, fixed = TRUE))
  expect_false(grepl("nf_resume", html, fixed = TRUE))
  hl <- as.character(RiboPilot:::run_modal(2, headless = TRUE,
    sched = list(cmd = "qsub", name = "SGE/PBS"), script = "x", pfx = "rc_", resume = FALSE))
  for (id in c("submit_script_text", "submit_job", "save_script", "nf_close", "nf_cmd_text")) {
    expect_true(has_id(hl, paste0("rc_", id)), info = id)
    expect_false(has_id(hl, id), info = id)
  }
  # the Assemble step's window keeps its ids, title, and Resume box
  def <- as.character(RiboPilot:::run_modal(3))
  for (id in c(local_ids, "nf_resume")) expect_true(has_id(def, id), info = id)
  expect_true(grepl("RiboPilot pipeline - 3 sample(s) or contig(s) to update", def, fixed = TRUE))
})
```

Ids are matched as `id="rc_nf_start"` because a bare `grepl("nf_start")` also matches `rc_nf_start`.

- [ ] **Step 6: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-headless.R")'`
Expected: FAIL, `unused arguments (pfx = "rc_", title = "RiboPilot recheck - 2 contig(s) to update", resume = FALSE)`.

- [ ] **Step 7: Give run_modal a prefix, a title, and an optional Resume box**

In `R/app_modals.R`, find:

```r
#' cluster job script to Submit as job / Save instead.
#' @noRd
run_modal <- function(n, headless = FALSE, sched = NULL, script = "") {
  wf <- "RiboPilot pipeline"
```

Replace with:

```r
#' cluster job script to Submit as job / Save instead.
#' @param pfx Prefix on every input and output id: "" for the Assemble step's
#'   Update, "rc_" for the recheck on Annotation Review.
#' @param title Window title; NULL gives the pipeline's.
#' @param resume Offer the Resume box; the recheck has no -resume.
#' @noRd
run_modal <- function(n, headless = FALSE, sched = NULL, script = "", pfx = "",
                      title = NULL, resume = TRUE) {
  id <- function(x) paste0(pfx, x)
  title <- title %||% sprintf("RiboPilot pipeline - %d sample(s) or contig(s) to update", n)
```

Find:

```r
      tags$textarea(id = "submit_script_text", rows = 11,
```

Replace with:

```r
      tags$textarea(id = id("submit_script_text"), rows = 11,
```

Find:

```r
    shinyjs::hidden(div(id = "nf_progress",
```

Replace with:

```r
    shinyjs::hidden(div(id = id("nf_progress"),
```

Find:

```r
          textOutput("nf_prog_all", container = tags$code)))))
  }
  footer <- if (headless) {
    mp_footer(
      primary = if (!is.null(sched)) actionButton("submit_job", "Submit as job"),
      dismiss = NULL,
      extra = tagList(
        actionButton("save_script", "Save script"),
        actionButton("nf_close", "Close")))
  } else {
    mp_footer(
      primary = actionButton("nf_start", "Run from App"), dismiss = NULL,
      extra = tagList(
        shinyjs::hidden(actionButton("nf_stop", "Stop")),
        actionButton("nf_close", "Close")))
  }
  modalDialog(
    title = run_modal_title(sprintf("%s - %d sample(s) or contig(s) to update", wf, n), "rp_gears"),
    size = "l",
    shinyWidgets::prettyCheckbox("nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("nf_cmd_text")),
```

Replace with:

```r
          textOutput(id("nf_prog_all"), container = tags$code)))))
  }
  footer <- if (headless) {
    mp_footer(
      primary = if (!is.null(sched)) actionButton(id("submit_job"), "Submit as job"),
      dismiss = NULL,
      extra = tagList(
        actionButton(id("save_script"), "Save script"),
        actionButton(id("nf_close"), "Close")))
  } else {
    mp_footer(
      primary = actionButton(id("nf_start"), "Run from App"), dismiss = NULL,
      extra = tagList(
        shinyjs::hidden(actionButton(id("nf_stop"), "Stop")),
        actionButton(id("nf_close"), "Close")))
  }
  modalDialog(
    title = run_modal_title(title, id("rp_gears")),
    size = "l",
    if (resume) shinyWidgets::prettyCheckbox(id("nf_resume"), "Resume previous run?",
      value = TRUE, status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput(id("nf_cmd_text"))),
```

The pipeline's call in `R/app_server.R` (`run_modal(n, headless = hl, sched = sched, script = script)`) stays as it is and renders exactly as before.

- [ ] **Step 8: Run it to see it pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-headless.R")'`
Expected: PASS, 8 tests, 0 failures (the test at :27-37 still finds the unprefixed ids of the pipeline window).

- [ ] **Step 9: Write the failing submit script tests**

Append to the end of `tests/testthat/test-nf-command-quoting.R`:

```r

test_that("the recheck submit script makes and enters .recheck, with no -resume", {
  skip_if_not_installed("withr")
  proj <- file.path(withr::local_tempdir(), "my project")
  dir.create(proj)
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "}"))
  script <- RiboPilot:::build_recheck_submit_script(proj, sched = NULL)
  words <- shell_words(script)
  expect_equal(words[which(words == "--projectDir") + 1], proj)
  expect_equal(words[which(words == "-work-dir") + 1], file.path(proj, "work_recheck"))
  expect_false("-resume" %in% words)
  expect_equal(dirname(words[which(words == "-log") + 1]), file.path(proj, ".runs", "nextflow"))
  expect_equal(utils::tail(words, 2), c("-ansi-log", "false"))
  lines <- strsplit(script, "\n", fixed = TRUE)[[1]]
  rc <- file.path(proj, ".recheck")
  expect_lt(match(paste("mkdir -p", shQuote(rc)), lines), match(paste("cd", shQuote(rc)), lines))
  expect_lt(match(paste("cd", shQuote(rc)), lines), grep("^'nextflow'", lines))
})
```

Append to the end of `tests/testthat/test-hydra-setup.R`:

```r

test_that("a recheck submission's script and logs share one base name", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    nf_recheck_command = function(path) c("run", "main_recheck.nf"),
    .package = "RiboPilot")
  base <- RiboPilot:::job_basename("recheck", as.POSIXct("2026-10-08 09:45:03", tz = "UTC"))
  expect_equal(base, "recheck_2026-10-08_09-45-03")
  txt <- RiboPilot:::build_recheck_submit_script(dir, sched = list(cmd = "qsub"), base = base)
  expect_match(txt, file.path(dir, ".runs", "logs", paste0(base, ".$JOB_ID.log")), fixed = TRUE)
  expect_match(txt, file.path(dir, ".runs", "nextflow", paste0(base, ".nextflow.log")), fixed = TRUE)
  expect_match(txt, paste0("#$ -N ", base), fixed = TRUE)
  expect_true(dir.exists(file.path(dir, ".runs", "logs")))
  expect_true(dir.exists(file.path(dir, ".runs", "nextflow")))
})
```

- [ ] **Step 10: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-nf-command-quoting.R", "test-hydra-setup.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: FAIL, `object 'build_recheck_submit_script' not found` in both new tests.

- [ ] **Step 11: Add the recheck job script builder**

In `R/app_run.R`, find:

```r
  cmd <- c(paste("cd", shQuote(file.path(project, ".reference"))),
    paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "))
  job_script(project, cmd, sched, base)
}
```

Replace with:

```r
  cmd <- c(paste("cd", shQuote(file.path(project, ".reference"))),
    paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "))
  job_script(project, cmd, sched, base)
}

#' Cluster job script for a recheck run
#'
#' Launched from its own folder, so the pipeline's -resume history in the
#' project folder is never touched; the job creates that folder itself.
#' @param base Timestamped name for this submission's script and logs.
#' @noRd
build_recheck_submit_script <- function(project, sched = detect_scheduler(),
                                        base = job_basename("recheck")) {
  rc <- file.path(project, ".recheck")
  args <- c(nf_log_opt(project, base), nf_recheck_command(project), "-ansi-log", "false")
  cmd <- c(paste("mkdir -p", shQuote(rc)), paste("cd", shQuote(rc)),
    paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "))
  job_script(project, cmd, sched, base)
}
```

Find:

```r
#' Shared by the Update window (`pfx = ""`) and the Reference build window
```

Replace with:

```r
#' Shared by the Assemble step's Update window (`pfx = ""`), the recheck window
#' on Annotation Review (`pfx = "rc_"`), and the Reference build window
```

Find:

```r
#' @param workflow Workflow named in the run's files ("pipeline", "reference").
```

Replace with:

```r
#' @param workflow Workflow named in the run's files ("pipeline", "recheck",
#'   "reference").
```

- [ ] **Step 12: Run them to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-nf-command-quoting.R", "test-hydra-setup.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: PASS, 0 failures (test-nf-command-quoting.R 11 tests; test-hydra-setup.R 15 tests, one skipped by `skip_on_cran()` under `test_file`).

- [ ] **Step 13: Write the failing app tests**

In `tests/testthat/test-app-run.R`, `mock_run()` also stands in for the recheck command (with the real one-argument signature, so an extra argument at a call site fails the launch) and records the launch folder. Find:

```r
    nf_command = function(...) c("run", "main.nf"),
    nf_spawn = function(cmd, args, wd, env) {
      seen$args <- args
```

Replace with:

```r
    nf_command = function(...) c("run", "main.nf"),
    nf_recheck_command = function(path) c("run", "main_recheck.nf"),
    nf_spawn = function(cmd, args, wd, env) {
      seen$args <- args
      seen$wd <- wd
```

Append to the end of `tests/testthat/test-app-run.R`:

```r

# ---- Update on Annotation Review: the recheck window ----

# S1 with one kept, annotated contig whose review copy is queued for a recheck
mk_recheck_project <- function(env = parent.frame()) {
  d <- mk_run_project(env)
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, n_features = 3")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
    VALUES ('S1', 1, 1, 4, 'ACGT', 0)")
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length,
    trim_start, trim_end, edited, check_switch, rev, time_stamp)
    VALUES ('S1', 1, 1, 'ACGT', 4, 1, 4, 1, 1, 1, 100)")
  DBI::dbDisconnect(con)
  d
}

# what query_recheck does when the run picks the contig up
start_copy <- function(d) {
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 4, time_stamp = ?",
    params = list(as.integer(Sys.time())))
  DBI::dbDisconnect(con)
}

copy_state <- function(d) {
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "SELECT check_switch, check_note FROM review_units")
}

test_that("contigs waiting for the Assemble step are unlocked, kept, uncopied, and not done", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "DELETE FROM annotate")
  u <- data.frame(scaffold = 1:7, annotate_switch = c(1L, 2L, 2L, 4L, 1L, 1L, 1L),
    annotate_lock = c(0L, 0L, 0L, NA, 1L, 0L, 0L), n_features = c(NA, 3L, NA, 3L, NA, NA, NA),
    ignore = c(0L, 0L, 0L, 0L, 0L, 1L, 0L))
  for (i in seq_len(nrow(u))) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES ('S1', 1, ?, 4, 'ACGT', ?)", params = list(u$scaffold[i], u$ignore[i]))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
      n_features, annotate_opts, curate_opts) VALUES ('S1', 1, ?, ?, ?, ?, 'default', 'default')",
      params = list(u$scaffold[i], u$annotate_switch[i], u$annotate_lock[i], u$n_features[i]))
  }
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, check_switch, rev, time_stamp)
    VALUES ('S1', 1, 7, 1, 1, 100)")
  DBI::dbDisconnect(con)
  # 1 queued, 3 never annotated, 4 running with a NULL lock; not 2 (done), 5 (locked),
  # 6 (dropped), or 7 (has a copy)
  expect_identical(RiboPilot:::nf_recheck_waiting(d), 3L)
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  expect_identical(RiboPilot:::nf_recheck_waiting(d), 0L)
})

test_that("with nothing queued, the recheck says how many contigs wait for the Assemble step", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
    VALUES ('S1', 1, 1, 4, 'ACGT', 0)")
  DBI::dbDisconnect(con)
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
  })
  expect_equal(seen$alert$title, "Nothing to update")
  expect_match(seen$alert$msg, "Waiting for the Assemble step's Update: 1 contig.", fixed = TRUE)
  expect_null(seen$p)
})

test_that("the recheck will not open while another run holds the project", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  d <- mk_recheck_project()
  withr::local_dir(d)
  m <- RiboPilot:::run_lock_take(d, "pipeline")
  other <- processx::process$new("sleep", "30")
  Sys.unsetenv(m)
  on.exit(other$kill(), add = TRUE)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
  })
  expect_match(seen$alert$msg, "already going")
  expect_null(seen$p)
})

test_that("Update on Annotation Review opens the recheck window, not the pipeline's", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  testthat::local_mocked_bindings(run_modal = function(n, ...) {
    seen$modal <- list(n = n, ...)
    shiny::modalDialog()
  }, .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    expect_equal(output$rc_nf_cmd_text, "nextflow run main_recheck.nf")
  })
  expect_equal(seen$modal$n, 1L)
  expect_equal(seen$modal$pfx, "rc_")
  expect_false(seen$modal$resume)
  expect_equal(seen$modal$title, "RiboPilot recheck - 1 contig(s) to update")
})

test_that("a recheck runs from .recheck, reports contigs it left running, then queues them again", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_nf_start = 1)
    expect_true(file.exists(RiboPilot:::run_lock_path(d)))
    start_copy(d)
    seen$p$finish()
    session$elapse(300)
  })
  expect_equal(seen$wd, file.path(d, ".recheck"))
  expect_true(dir.exists(file.path(d, ".recheck")))
  expect_true("main_recheck.nf" %in% seen$args)
  lg <- seen$args[which(seen$args == "-log") + 1]
  expect_equal(dirname(lg), file.path(d, ".runs", "nextflow"))
  expect_match(basename(lg), "^recheck_")
  # the report is written before the release, so the contig counts as a failure
  expect_equal(seen$alert$title, "Run finished with failures")
  expect_match(seen$alert$msg, "S1.1.1 | recheck did not finish", fixed = TRUE)
  expect_length(list.files(file.path(d, ".runs", "reports"), "^recheck_.*_report\\.txt$"), 1)
  expect_equal(copy_state(d)$check_switch, 1L)
  expect_equal(copy_state(d)$check_note, "recheck did not finish")
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("a stopped recheck points back to Annotation Review and queues its contigs again", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_nf_start = 1)
    start_copy(d)
    session$setInputs(rc_nf_stop = 1)
    session$elapse(300)
  })
  expect_true(isTRUE(seen$interrupted))
  expect_equal(seen$alert$title, "Run stopped")
  expect_match(seen$alert$msg, "Press Update on Annotation Review.", fixed = TRUE)
  expect_equal(copy_state(d)$check_switch, 1L)
})

test_that("a headless recheck saves its job script under the name of its logs", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  withr::local_options(list(RiboPilot.headless = TRUE))
  seen <- new.env(); mock_run(seen)
  testthat::local_mocked_bindings(detect_scheduler = function() NULL, .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_submit_script_text = "   ")
    session$setInputs(rc_save_script = 1)
  })
  f <- Sys.glob(file.path(d, ".runs", "jobs", "recheck_*.sh"))
  expect_length(f, 1)
  expect_length(Sys.glob(file.path(d, ".runs", "jobs", "pipeline_*.sh")), 0)
  s <- readLines(f)
  base <- sub("\\.sh$", "", basename(f))
  expect_true(any(grepl(paste0(base, ".nextflow.log"), s, fixed = TRUE)))
  expect_true(paste("cd", shQuote(file.path(d, ".recheck"))) %in% s)
  expect_true(any(grepl("'main_recheck.nf'", s, fixed = TRUE)))
  expect_match(seen$alert$msg, "Saved job script to", fixed = TRUE)
})

test_that("Submit as job hands the recheck script to the scheduler once", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  d <- mk_recheck_project()
  withr::local_dir(d)
  withr::local_options(list(RiboPilot.headless = TRUE))
  bin <- withr::local_tempdir()
  calls <- file.path(bin, "calls.txt")
  writeLines(c("#!/bin/sh", sprintf('echo "$@" >> %s', shQuote(calls)),
    'echo "Your job 12345 has been submitted"'), file.path(bin, "qsub"))
  Sys.chmod(file.path(bin, "qsub"), "0755")
  withr::local_path(bin, action = "prefix")
  skip_if(!identical(RiboPilot:::detect_scheduler()$cmd, "qsub"), "another scheduler on PATH")
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_submit_job = 1)
    while (!later::loop_empty()) later::run_now()
    expect_match(seen$alert$msg, "Submitted via qsub", fixed = TRUE)
    session$setInputs(rc_submit_job = 2)
    while (!later::loop_empty()) later::run_now()
    expect_match(seen$alert$msg, "already submitted", fixed = TRUE)
  })
  expect_length(readLines(calls), 1L)
  expect_match(readLines(calls), "/\\.runs/jobs/recheck_[0-9_-]+\\.sh")
})
```

In `tests/testthat/test-app-opts.R`, the old test pinned "no Update on Annotate" and would keep passing after the button lands (the button sits past its 400-character window), so it is replaced on purpose by a positive check. Find:

```r
test_that("the Annotate tab has no Update button", {
  ui <- paste(deparse(RiboPilot:::app_ui), collapse = "\n")
  i <- regexpr("input.mode == 'Annotate'", ui, fixed = TRUE)
  expect_true(i > 0)
  block <- substr(ui, i, i + 400)
  expect_false(grepl("run_btn", block, fixed = TRUE))
})
```

Replace with:

```r
test_that("the Annotation Review toolbar has its own Update button", {
  ui <- paste(deparse(RiboPilot:::app_ui), collapse = "\n")
  i <- regexpr("input.mode == 'Annotate'", ui, fixed = TRUE)
  j <- regexpr("input.mode == 'Export'", ui, fixed = TRUE)
  expect_true(i > 0 && j > i)
  block <- substr(ui, i, j)
  expect_true(grepl("annot_run_btn", block, fixed = TRUE))
  expect_false(grepl("asmb_run_btn", block, fixed = TRUE))
  html <- as.character(RiboPilot:::app_ui(NULL))
  expect_true(grepl('id="annot_run_btn"', html, fixed = TRUE))
})
```

- [ ] **Step 14: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-app-run.R", "test-app-opts.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: FAIL. test-app-run.R: the 12 existing tests pass; the 8 new ones fail, in order, with `object 'nf_recheck_waiting' not found`; `seen$alert$title` NULL instead of "Nothing to update"; `seen$alert$msg` NULL (no lock message); `The test referenced an output that hasn't been defined yet: output$rc_nf_cmd_text`; `attempt to apply non-function` at `seen$p$finish()` (nothing launched); `isTRUE(seen$interrupted)` FALSE; no `recheck_*.sh` saved; and `seen$alert$msg` NULL where "Submitted via qsub" is expected. test-app-opts.R: `grepl("annot_run_btn", block, fixed = TRUE)` is FALSE, and so is the `id="annot_run_btn"` check.

- [ ] **Step 15: Add the waiting count, the Update button, and the recheck window**

Append to the end of `R/nf_command.R`:

```r

#' Count the contigs that wait for the Assemble step's Update
#'
#' Unlocked, kept contigs of assembly-locked samples that have no review copy
#' and whose original is queued for the main pipeline or was never annotated.
#' The recheck cannot take them; "Nothing to update" says how many there are.
#' @noRd
nf_recheck_waiting <- function(path) {
  con <- open_project_db(file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  as.integer(DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       JOIN assemble a ON a.ID = an.ID
      WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0
        AND COALESCE(an.annotate_lock, 0) = 0
        AND NOT EXISTS (SELECT 1 FROM review_units r
                         WHERE r.ID = an.ID AND r.path = an.path AND r.scaffold = an.scaffold)
        AND (an.annotate_switch IN (1,4) OR an.n_features IS NULL)")$n)
}
```

In `R/app_ui.R`, find:

```r
            # No Update button here: one workflow, launched from Assemble.
            shiny::conditionalPanel("input.mode == 'Annotate'", id = "annot_ctrls", class = "mp-toolbar",
```

Replace with:

```r
            shiny::conditionalPanel("input.mode == 'Annotate'", id = "annot_ctrls", class = "mp-toolbar",
```

Find:

```r
                title = MP_TRIM_TIP()[["restore"]], needs_selection = TRUE)),
```

Replace with:

```r
                title = MP_TRIM_TIP()[["restore"]], needs_selection = TRUE),
              btn("annot_run_btn", "Update", icon = mp_icon("circle-play"), emphasis = "primary",
                title = "Review and launch the recheck of the contigs that are Ready to run")),
```

Find:

```r
                title = "Read the reports of past pipeline runs and reference builds"),
```

Replace with:

```r
                title = "Read the reports of past pipeline runs, rechecks, and reference builds"),
```

In `R/app_server.R`, insert the recheck block between the pipeline's run window and `register_app_lifecycle(session)` (the lifecycle grace period must start after every run window is wired). Leave the `tax_run_btn`/`asmb_run_btn` observer, `nf_pending_count`, and the pipeline's save and submit observers unchanged. Find:

```r
    on_done = function(how) refresh(isolate(refresh()) + 1))

  # Stop the app when the last browser session closes (survives refresh). After
```

Replace with:

```r
    on_done = function(how) refresh(isolate(refresh()) + 1))

  # ---- Recheck: Update on Annotation Review reruns the checks on review copies ----
  rc_submit_base <- reactiveVal(NULL)   # timestamp shared by one submission's files

  output$rc_nf_cmd_text <- renderText({
    input$annot_run_btn
    paste(c("nextflow", nf_recheck_command(app_project())), collapse = " ")
  })

  observeEvent(input$annot_run_btn, ignoreInit = TRUE, {
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    n <- nf_recheck_pending(app_project())
    if (n == 0) {
      w <- nf_recheck_waiting(app_project())
      rp_alert(paste(c("No contigs are Ready to run. A trim or a moved gene boundary queues a",
        "contig, and so does setting its State to Ready to run. Locked contigs are skipped.",
        if (w > 0) paste0("Waiting for the Assemble step's Update: ", mp_n(w, "contig"), ".")),
        collapse = " "), type = "warning", title = "Nothing to update")
      return()
    }
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    if (hl) rc_submit_base(job_basename("recheck"))
    script <- if (hl) build_recheck_submit_script(app_project(), sched = sched,
      base = rc_submit_base()) else ""
    showModal(run_modal(n, headless = hl, sched = sched, script = script, pfx = "rc_",
      title = sprintf("RiboPilot recheck - %d contig(s) to update", n), resume = FALSE))
  })

  write_recheck_script <- function() {
    base <- rc_submit_base() %||% job_basename("recheck")
    typed <- input$rc_submit_script_text %||% ""
    script <- if (nzchar(trimws(typed))) typed else
      build_recheck_submit_script(app_project(), base = base)
    f <- job_path(app_project(), "script", paste0(base, ".sh"))
    writeLines(script, f)
    try(Sys.chmod(f, "0755"), silent = TRUE)
    f
  }
  observeEvent(input$rc_save_script, {
    f <- write_recheck_script()
    rp_alert(paste0("Saved job script to ", f), type = "message")
  })
  rc_job_sent <- reactiveVal(FALSE)
  rc_resubmit_ok <- reactiveVal(FALSE)
  observeEvent(input$rc_submit_job, {
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    if (isTRUE(rc_job_sent()) && !isTRUE(rc_resubmit_ok())) {
      rc_resubmit_ok(TRUE)
      rp_alert(paste("A recheck job for this project was already submitted from this window and",
        "may still be queued or running. Press Submit as job again to submit another."),
        type = "warning")
      return()
    }
    rc_resubmit_ok(FALSE)
    sched <- detect_scheduler()
    if (is.null(sched)) {
      rp_alert("No scheduler (sbatch/qsub/bsub) found on PATH.", type = "warning"); return()
    }
    f <- write_recheck_script()
    shinyWidgets::sendSweetAlert(title = "Submitting...",
      text = paste0("Handing the job script to ", sched$cmd,
        ". Hold tight, a busy cluster can take a few seconds to answer."),
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        out <- tryCatch(submit_job(sched$cmd, f), error = function(e) conditionMessage(e))
        rc_job_sent(TRUE)
        rp_alert(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " ")),
          type = "message")
        removeModal()
      })
    }, delay = 0.1)
  })

  # Launched from .recheck so the pipeline's -resume history is never touched.
  # nf_run_window writes the report before on_done, so contigs the run left
  # running are reported as failures before they go back to the queue.
  nf_run_window(input, output, session, pfx = "rc_", gears = "rc_rp_gears",
    project = app_project, kind = "recheck", workflow = "recheck",
    launch = function(nflog) {
      wd <- file.path(app_project(), ".recheck")
      dir.create(wd, recursive = TRUE, showWarnings = FALSE)
      nf_spawn("nextflow", c("-log", nflog, nf_recheck_command(app_project()), "-ansi-log"),
        wd, nf_env())
    },
    on_done = function(how) {
      try(local({
        con <- app_con()
        on.exit(DBI::dbDisconnect(con))
        recheck_release_running(con)
      }), silent = TRUE)
      refresh(isolate(refresh()) + 1)
    })

  # Stop the app when the last browser session closes (survives refresh). After
```

Notes for this step: `output$rc_nf_cmd_text` reads `input$annot_run_btn` so the preview is recomputed on each press (testServer computes it on every flush, so `nf_recheck_command()` must stay free of side effects, as Task 11 defines it). `launch` creates `.recheck` itself because processx refuses a missing working directory. The idle `rc_` window's session-end handler returns at once, so a session that never ran a recheck still writes exactly one report (test-app-run.R:134-148).

- [ ] **Step 16: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("test-app-run.R", "test-app-opts.R", "test-app-headless.R", "test-run-reports.R", "test-nf-run.R", "test-nf-command-quoting.R", "test-hydra-setup.R", "test-nf-pending.R", "test-app-scaffold-uid.R", "test-reference-modal.R", "test-app-table.R", "test-app-export.R", "test-docs-site.R", "test-no-duplicate-definitions.R", "test-batch-trim.R", "test-blast-provenance.R", "test-notes-partial.R", "test-specimen-app.R", "test-reference-app.R", "test-ui-helpers.R")) testthat::test_file(file.path("tests/testthat", f))'`
Expected: PASS, 0 failures in every file. test-app-run.R has 20 tests (12 existing, 8 new). The `skip_on_cran()` tests in test-nf-run.R (2) and test-hydra-setup.R (1) are skipped under `test_file`, as before. test-app-scaffold-uid.R:571-601 still saves `pipeline_*.sh` through the unprefixed `submit_script_text`/`save_script`, and test-docs-site.R finds no "<Step> tab" text in the new strings.

- [ ] **Step 17: Commit**

```bash
git add R/run_reports.R R/app_run_reports.R R/app_modals.R R/app_run.R R/nf_command.R R/app_ui.R R/app_server.R tests/testthat/test-nf-run.R tests/testthat/test-run-reports.R tests/testthat/test-app-headless.R tests/testthat/test-nf-command-quoting.R tests/testthat/test-hydra-setup.R tests/testthat/test-app-run.R tests/testthat/test-app-opts.R
git commit -m "Update button on Annotation Review launches the recheck"
```


### Task 14: Annotation Review table, State column, Checks

**Files:**
- Modify: `R/review.R` (append `review_unit_state()` and `review_state_text()` at the end of the file)
- Modify: `R/export.R` (`unit_checks()`, :247-298)
- Modify: `R/app_server.R` (Assemble branch ITSx roll-up :184-190; Annotate branch of `all_rows_base` :219-245; `with_checks()` :253-255; `output$tbl` colDefs :320-342 and the State column :383-384)
- Modify: `R/utils_reactable.R` (`rt_dynamicIcon()` :50-91; `rt_itsx_badge()` :227-242; new `rt_gated()` after it)
- Modify: `R/constants.R` (`MP_STATE_META` :208-209, `MP_STATE_CODES` :221-227)
- Modify: `R/app_ui_helpers.R` (`mp_state_icon_classes()` :165-170)
- Modify: `R/app_ui.R` (`filter_pickers()` :5-15)
- Create: `tests/testthat/test-review-state.R`
- Modify: `tests/testthat/test-export.R` (:286-325), `tests/testthat/test-app-export.R` (`mk_export_project()` :37; Task 1's view test at the end), `tests/testthat/test-app-opts.R` (:122), `tests/testthat/test-app-scaffold-uid.R` (:758-774), `tests/testthat/test-ui-helpers.R` (:33-34, one test appended), `tests/testthat/test-ui-reactable-helpers.R` (tests appended), `tests/testthat/test-app-table.R` (one test appended)
- Test: `tests/testthat/test-review-state.R`, `tests/testthat/test-export.R`, `tests/testthat/test-app-table.R`, `tests/testthat/test-app-scaffold-uid.R`, `tests/testthat/test-ui-helpers.R`, `tests/testthat/test-ui-reactable-helpers.R`

**Interfaces:**
- Consumes: `review_units_v`, `review_screen_v`, `review_chimera_v`, `review_annotations_v` (Task 1; `review_units_v` columns `ID, path, scaffold, has_copy, ..., check_switch, check_note, ..., n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, ..., blast_species, blast_accession, blast_pident, blast_qcovs, blast_source, ...`); `fixture_review_schema()` via `fixture_fill_schema()` in `mk_app_project()` / `mk_export_project()` (Task 1); `R/review.R` exists (Task 2); `nf_recheck_waiting(path)` (Task 13, R/nf_command.R), used only by a test; `unit_key(id, path, scaffold)` (R/export.R, unchanged).
- Produces:
  - `review_unit_state(con, units = NULL)` (R/review.R): data.frame(ID, path, scaffold, has_copy, check_switch, annotate_switch, annotate_lock, n_features, check_note, state). `units = NULL`: every kept contig (`COALESCE(assemblies.ignore, 0) = 0`) of samples with `assemble_lock = 1`, ordered by ID, path, scaffold. `units` = data.frame(ID, path, scaffold): one row per row of `units`, in its order, matched in R (no IN list); a unit with no `review_units_v` + `annotate` row gets NA columns and state "not run". `annotate_switch`, `annotate_lock`, `n_features` come from `annotate` (the original), `has_copy`, `check_switch`, `check_note` from `review_units_v`. `state`: "waiting" when `has_copy = 0` and (`annotate_switch` IN (1, 4) or `n_features` IS NULL), the condition `nf_recheck_waiting()` counts; otherwise from `check_switch`: 2 "done", 1 or NULL "not run", 4 "running", 3 "failed", 0 "on hold". This is the only definition of waiting (Tasks 15, 16, 17, 21 call it). `rev` is never selected.
  - `review_state_text(state)` (R/review.R): done -> "done", not run -> "not run", running -> "running", failed -> "failed", on hold -> "on hold", waiting -> "waiting for Assemble Update"; NA stays NA.
  - `unit_checks(con, by = c("unit", "sample"), review = FALSE)` (R/export.R): `review = TRUE` reads `review_units_v`, `review_screen_v`, `review_chimera_v`, `review_annotations_v` instead of `annotate`, `screen`, `chimera`, `annotations`; a group with any contig that is not done gets `checks` = its `review_state_text()` values (unique, ", "-joined) and `checks_tone = "info"`. Both modes: an `itsx_warnings` that is NULL or starts with `not_run(` adds the amber item "ITSx not checked" (after "ITSx"), so the cell is never "ok". Rows come ordered by ID, path, scaffold.
  - `rt_gated(cell = NULL)` (R/utils_reactable.R): wraps a cell renderer; a row whose `check_state` is set renders `<span class="mp-pill mp-pill-info" title="<check_note>">check_state</span>`, other rows (and tables without a `check_state` column) render through `cell`, or the plain value when `cell` is NULL. Task 20 reuses it on the Export step.
  - `rt_dynamicIcon(icons = NULL, labels = NULL, inputId = NULL, note_col = NULL)`: with `note_col`, a row's non-empty note is appended to the glyph's title as "<label> Note: <note>".
  - `rt_itsx_badge()`: NULL -> grey "not checked"; `''` -> blank; `not_run(...)` -> grey "not run"; ticks only text matching `/^0 of [0-9]+ contigs$/`; anything else an amber pill.
  - `MP_STATE_META[["5"]]` = list(label = "Waiting for Assemble Update", icon = "fa fa-clock", tip naming the Assemble step's Update); `MP_STATE_CODES$annotate = c("5", "0", "1", "4", "2", "3")`; `MP_STATE_SETTABLE` unchanged; `mp_state_icon_classes()` tone for 5 is `mp-fg-neutral`; the Annotation Review State filter offers `mp_state_choices("annotate")`.
  - Annotation Review table (`all_rows_base`, mode "Annotate"): `length, n_features, missing, blast_species, blast_accession, blast_pident, blast_qcovs, blast_source, n_flags, curate_flags, itsx_warnings` from `review_units_v`, `screen_verdict` from `review_screen_v`, `chimera_label, chimera_reason` from `review_chimera_v`; marks, notes, group, exported, time_stamp from `annotate`. Column `annotate_switch` = `COALESCE(review_units_v.check_switch, 1)`, set to 5 where `review_unit_state()` says waiting. Hidden columns `check_state` (NA when done, else `review_state_text(state)`) and `check_note`. `with_checks(d, by, review = FALSE)`; the Annotate branch passes `review = TRUE`, the Assemble branch keeps `with_checks("sample")` on originals.
  - Assemble step ITSx roll-up (column `itsx_warnings` of a sample row): over kept contigs only; NULL when none has a stored value; else `"<flagged> of <kept> contigs"`, followed by `", <k> not checked"` when k kept contigs are NULL or `not_run(...)`.

Annotation Review now reads contig data and check results through the review layer (S8): a contig with a review copy shows the copy, any other contig its original. The State column shows the checks' state (`review_units_v.check_switch`, S6.1), and a contig without a copy whose original is still queued for, or was never finished by, the main pipeline shows the new display-only code 5 "Waiting for Assemble Update" (S8, contract section 1). While a contig's checks are not done, every check cell and the Checks column show that state in an info pill instead of a value, which also hides rows a discarded recheck left behind (S5.4, S8). The Checks column and the Assemble step's ITSx roll-up apply the ITSx rule of S6.3: a NULL or `not_run` result is an amber "ITSx not checked" item and never reads as clean, and the roll-up counts only kept contigs. The State button still writes `annotate.annotate_switch` until Task 16 rewrites it on copies; the details window, Trim, and Lock keep reading originals until Tasks 15-17 and 21.

- [ ] **Step 1: Write the failing tests for the state and the Checks column**

Create `tests/testthat/test-review-state.R`:

```r
# Check state on Annotation Review (SPEC S6.1, S6.3, S8).

# Contigs 1-10 of S1 (assembly locked) cover every state, contig 11 is dropped,
# and S2's assembly is not locked. Copies hold only the columns read here.
rst_units <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon) VALUES ('S2', 'Danio rerio')")
  DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_lock) VALUES ('S2', 0)")
  orig <- data.frame(
    ID = c(rep("S1", 11), "S2"), scaffold = c(1:11, 1L),
    annotate_switch = c(1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 4L, 2L, 1L),
    annotate_lock = c(rep(0L, 9), 1L, 0L, 0L),
    n_features = c(2L, NA, rep(2L, 10)),
    itsx_warnings = c(NA, NA, NA, "", NA, "", "", "", "", NA, "", NA),
    ignore = c(rep(0L, 10), 1L, 0L))
  for (i in seq_len(nrow(orig))) {
    o <- orig[i, ]
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES (?, 1, ?, 4, 'ACGT', ?)", params = list(o$ID, o$scaffold, o$ignore))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
      n_features, itsx_warnings) VALUES (?, 1, ?, ?, ?, ?, ?)",
      params = list(o$ID, o$scaffold, o$annotate_switch, o$annotate_lock, o$n_features, o$itsx_warnings))
  }
  # copies: 5 failed (its original is queued again; the copy wins), 6 running,
  # 7 on hold, 8 with no state stored, 9 done
  cp <- data.frame(scaffold = 5:9, check_switch = c(3L, 4L, 0L, NA, 2L),
    check_note = c("boom", NA, NA, NA, NA))
  for (i in seq_len(nrow(cp))) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, trim_start,
      trim_end, edited, check_switch, check_note, rev, n_features, itsx_warnings)
      VALUES ('S1', 1, ?, 'ACGT', 4, 1, 4, 1, ?, ?, CAST(? AS INTEGER), 2, '')",
      params = list(cp$scaffold[i], cp$check_switch[i], cp$check_note[i], "1760000000000"))
  }
}

test_that("review_unit_state names every state, and waiting matches the recheck's waiting count", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rst_units(con)
  st <- RiboPilot:::review_unit_state(con)
  expect_equal(names(st), c("ID", "path", "scaffold", "has_copy", "check_switch", "annotate_switch",
    "annotate_lock", "n_features", "check_note", "state"))
  expect_equal(st$ID, rep("S1", 10))
  expect_equal(st$scaffold, 1:10)
  expect_equal(st$has_copy, c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L))
  expect_equal(st$state, c("waiting", "waiting", "not run", "done", "failed", "running", "on hold",
    "not run", "done", "waiting"))
  expect_equal(st$check_note, c(NA, NA, NA, NA, "boom", NA, NA, NA, NA, NA))
  expect_equal(RiboPilot:::review_state_text(st$state), c("waiting for Assemble Update",
    "waiting for Assemble Update", "not run", "done", "failed", "running", "on hold", "not run",
    "done", "waiting for Assemble Update"))
  # the Update alert counts the same contigs, less the locked one
  expect_equal(RiboPilot:::nf_recheck_waiting(dirname(con@dbname)), 2L)

  # given units: their order, any lock or drop state, and a unit that does not exist
  u <- data.frame(ID = c("S1", "S2", "S1", "SX"), path = 1L, scaffold = c(9L, 1L, 11L, 1L))
  st <- RiboPilot:::review_unit_state(con, u)
  expect_equal(st[c("ID", "path", "scaffold")], u)
  expect_equal(st$state, c("done", "waiting", "done", "not run"))
  expect_equal(nrow(RiboPilot:::review_unit_state(con, u[0, ])), 0L)
})

test_that("unit_checks on the review layer reads the copies and shows the state of unfinished checks", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
      VALUES ('S1', 1, ?, 4, 'ACGT')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, n_features, itsx_warnings)
      VALUES ('S1', 1, ?, ?, 2, '')", params = list(s, if (s == 4) 1L else 2L))
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, ?, 'on_target')",
      params = list(s))
  }
  copy <- function(s, sw, flags, itsx) DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold,
      sequence, length, trim_start, trim_end, edited, check_switch, rev, n_features, curate_flags,
      itsx_warnings) VALUES ('S1', 1, ?, 'ACG', 3, 1, 3, 1, ?, CAST(? AS INTEGER), 2, ?, ?)",
    params = list(s, sw, "1760000000000", flags, itsx))
  # 2: a done copy the screen calls off target, with a contaminant flag, a
  # flagged feature, and an ITSx check that could not decide
  copy(2L, 2L, "contaminant", "not_run(ITSx did not finish)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 2, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, status)
    VALUES ('S1', 1, 2, '28S', 'flag')")
  # 3: an edited copy not rechecked yet, under the severe chimera a discarded recheck left
  copy(3L, 1L, "", NA)
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label) VALUES ('S1', 1, 3, 'severe')")

  u <- RiboPilot:::unit_checks(con, "unit", review = TRUE)
  expect_equal(u$scaffold, 1:4)
  expect_equal(u$checks, c("ok", "off target, contaminant, ITSx not checked, feature flags", "not run",
    "waiting for Assemble Update"))
  expect_equal(u$checks_tone, c("success", "danger", "info", "info"))

  # the Assemble step reads the originals, which no copy changed
  o <- RiboPilot:::unit_checks(con, "unit")
  expect_equal(o$checks, rep("ok", 4))
  expect_equal(o$checks_tone, rep("success", 4))
})
```

In `tests/testthat/test-export.R`, the existing `unit_checks` test stores no ITSx value, which now reads as "ITSx not checked"; give its units a clean ITSx result and check the new item at the end.

Find:

```r
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1',1,?)", params = list(s))
```

Replace with:

```r
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, itsx_warnings) VALUES ('S1',1,?,'')",
      params = list(s))
```

Find:

```r
  expect_equal(RiboPilot:::unit_checks(con, "unit")[1, c("checks", "checks_tone")],
    data.frame(checks = "no reference", checks_tone = "neutral"))
})
```

Replace with:

```r
  expect_equal(RiboPilot:::unit_checks(con, "unit")[1, c("checks", "checks_tone")],
    data.frame(checks = "no reference", checks_tone = "neutral"))

  # never checked, or checked without a verdict: an amber item, never ok
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = NULL WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'not_run(ITSx did not finish)' WHERE scaffold = 2")
  u <- RiboPilot:::unit_checks(con, "unit")
  expect_equal(u$checks, c("ITSx not checked", "minor chimera, ITSx not checked, feature flags",
    "off target, contaminant"))
  expect_equal(u$checks_tone, c("warning", "warning", "danger"))
})
```

In `tests/testthat/test-app-export.R`, `mk_export_project()`'s locked contig gets a clean stored ITSx result, so it reads as done on Annotation Review (S6.1) and its Assemble Checks cell stays "contaminant" (:414). Find:

```r
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,2,1,'default','default','yes','no','no',3,'',0,NULL,NULL,NULL,1)")
```

Replace with:

```r
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,2,1,'default','default','yes','no','no',3,'',0,NULL,'',NULL,1)")
```

Task 1's view test at the end of the same file pinned the old fixture state. Find:

```r
  # locked and annotated, but no stored ITSx value: not done until fixture_mark_done()
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v"),
    data.frame(has_copy = 0L, length = 600L, check_switch = 1L))
```

Replace with:

```r
  # locked, annotated, and checked by ITSx: done
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v"),
    data.frame(has_copy = 0L, length = 600L, check_switch = 2L))
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-state", "export", "app-export")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-state.R `[ FAIL 2 | WARN 0 | SKIP 0 | PASS 0 ]`: `object 'review_unit_state' not found` and `unused argument (review = TRUE)`. test-export.R `[ FAIL 2 | WARN 0 | SKIP 0 | PASS 79 ]`: the new `u$checks` expectation gets `"no reference"` and `"minor chimera, feature flags"` (no "ITSx not checked" item yet), and `u$checks_tone` starts with `"neutral"`. test-app-export.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 142 ]` (the fixture change and the updated view test agree before and after this task).

- [ ] **Step 3: Add the state functions and the review mode of unit_checks**

Append to the end of `R/review.R`:

```r

#' Check state of contigs, as Annotation Review shows it
#'
#' One row per kept contig of the assembly-locked samples, ordered by ID, path,
#' scaffold; or, with `units` (data.frame(ID, path, scaffold)), one row per
#' unit in its order, matched in R. "waiting": no copy, and the original is
#' queued for the main pipeline or was never annotated (the contigs
#' nf_recheck_waiting() counts). Otherwise check_switch decides: 2 done,
#' 1 or NULL not run, 4 running, 3 failed, 0 on hold.
#' @noRd
review_unit_state <- function(con, units = NULL) {
  d <- DBI::dbGetQuery(con, paste(
    "SELECT v.ID, v.path, v.scaffold, v.has_copy, v.check_switch, an.annotate_switch,
       an.annotate_lock, an.n_features, v.check_note
     FROM review_units_v v
     JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold",
    if (is.null(units)) {
      "JOIN assemblies asm ON asm.ID = v.ID AND asm.path = v.path AND asm.scaffold = v.scaffold
     JOIN assemble a ON a.ID = v.ID
     WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0"
    },
    "ORDER BY v.ID, v.path, v.scaffold"))
  if (!is.null(units)) {
    d <- d[match(unit_key(units$ID, units$path, units$scaffold), unit_key(d$ID, d$path, d$scaffold)), ,
      drop = FALSE]
    d$ID <- as.character(units$ID)
    d$path <- as.integer(units$path)
    d$scaffold <- as.integer(units$scaffold)
    rownames(d) <- NULL
  }
  waiting <- d$has_copy %in% 0 & (d$annotate_switch %in% c(1, 4) | is.na(d$n_features))
  st <- unname(c(`0` = "on hold", `1` = "not run", `2` = "done", `3` = "failed",
    `4` = "running")[as.character(d$check_switch)])
  d$state <- as.character(ifelse(waiting, "waiting", ifelse(is.na(st), "not run", st)))
  d
}

#' Words for review_unit_state()'s state, as the table cells show them
#' @noRd
review_state_text <- function(state) {
  unname(c(done = "done", `not run` = "not run", running = "running", failed = "failed",
    `on hold` = "on hold", waiting = "waiting for Assemble Update")[state])
}
```

In `R/export.R`, find:

```r
#' Findings are listed worst first. The red ones (off target, contaminant,
#' severe chimera, ITSx chimeric) are the tags `has_screen_risk()` reads, so
#' Lock and Export ask about them. "ok" needs every unit screened on_target;
#' with nothing flagged and no usable screen the cell says "no reference".
#' Only kept scaffolds count.
#' @param con Open project DB connection.
#' @param by "unit" (one row per ID, path, scaffold) or "sample" (one per ID).
#' @return data.frame(ID, \[path, scaffold,\] checks, checks_tone)
#' @noRd
unit_checks <- function(con, by = c("unit", "sample")) {
  by <- match.arg(by)
  u <- DBI::dbGetQuery(con,
    "SELECT an.ID, an.path, an.scaffold, sc.screen_verdict, ch.label AS chimera_label,
       an.chimeric, an.curate_flags, an.itsx_warnings,
       (SELECT COUNT(*) FROM annotations f WHERE f.ID = an.ID AND f.path = an.path
          AND f.scaffold = an.scaffold AND f.status = 'flag') AS n_feat_flags
     FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
       LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
     WHERE COALESCE(asm.ignore, 0) = 0")
```

Replace with:

```r
#' Findings are listed worst first. The red ones (off target, contaminant,
#' severe chimera, ITSx chimeric) are the tags `has_screen_risk()` reads, so
#' Lock and Export ask about them. "ok" needs every unit screened on_target;
#' with nothing flagged and no usable screen the cell says "no reference". An
#' ITSx result that is missing or not_run is the amber "ITSx not checked".
#' Only kept scaffolds count.
#' @param con Open project DB connection.
#' @param by "unit" (one row per ID, path, scaffold) or "sample" (one per ID).
#' @param review TRUE reads the review layer (Annotation Review): a contig
#'   whose checks are not done shows its state instead of findings.
#' @return data.frame(ID, \[path, scaffold,\] checks, checks_tone)
#' @noRd
unit_checks <- function(con, by = c("unit", "sample"), review = FALSE) {
  by <- match.arg(by)
  tb <- if (review) {
    c("review_units_v", "review_screen_v", "review_chimera_v", "review_annotations_v")
  } else {
    c("annotate", "screen", "chimera", "annotations")
  }
  u <- DBI::dbGetQuery(con, sprintf(
    "SELECT an.ID, an.path, an.scaffold, sc.screen_verdict, ch.label AS chimera_label,
       an.chimeric, an.curate_flags, an.itsx_warnings,
       (SELECT COUNT(*) FROM %4$s f WHERE f.ID = an.ID AND f.path = an.path
          AND f.scaffold = an.scaffold AND f.status = 'flag') AS n_feat_flags
     FROM %1$s an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       LEFT JOIN %2$s sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
       LEFT JOIN %3$s ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
     WHERE COALESCE(asm.ignore, 0) = 0
     ORDER BY an.ID, an.path, an.scaffold", tb[1], tb[2], tb[3], tb[4]))
```

Find:

```r
    "ITSx" = nzchar(itsx) & !startsWith(itsx, "not_run("),
```

Replace with:

```r
    "ITSx" = nzchar(itsx) & !startsWith(itsx, "not_run("),
    "ITSx not checked" = is.na(u$itsx_warnings) | startsWith(itsx, "not_run("),
```

Find:

```r
  rows <- split(seq_len(nrow(u)), factor(key, levels = unique(key)))
  cell <- vapply(rows, function(i) {
    h <- colSums(hit[i, , drop = FALSE]) > 0
```

Replace with:

```r
  gate <- rep(NA_character_, nrow(u))
  if (review) {
    st <- review_unit_state(con, u[c("ID", "path", "scaffold")])$state
    gate <- ifelse(st %in% "done", NA_character_, review_state_text(st))
  }
  rows <- split(seq_len(nrow(u)), factor(key, levels = unique(key)))
  cell <- vapply(rows, function(i) {
    g <- unique(stats::na.omit(gate[i]))
    if (length(g)) return(c(paste(g, collapse = ", "), "info"))
    h <- colSums(hit[i, , drop = FALSE]) > 0
```

Check the files parse: `Rscript -e 'invisible(parse("R/review.R")); invisible(parse("R/export.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-state", "export", "app-export")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-state.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`, test-export.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 81 ]`, test-app-export.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 142 ]`.

- [ ] **Step 5: Write the failing tests for the table, the State column, and the ITSx display**

In `tests/testthat/test-ui-helpers.R` (Annotation Review's State column gains code 5), find:

```r
  expect_equal(names(mp_state_choices("annotate")),
               c("On hold", "Ready to run", "In progress", "Success", "Failed"))
```

Replace with:

```r
  expect_equal(names(mp_state_choices("annotate")),
               c("Waiting for Assemble Update", "On hold", "Ready to run", "In progress", "Success", "Failed"))
```

Append to the end of `tests/testthat/test-ui-helpers.R`:

```r

test_that("only Annotation Review offers Waiting for Assemble Update, and nobody can set it", {
  expect_match(as.character(filter_pickers("Annotate")), 'value="5"', fixed = TRUE)
  expect_false(grepl('value="5"', as.character(filter_pickers("Assemble")), fixed = TRUE))
  expect_equal(unname(mp_state_icon_classes("annotate")[["5"]]), "fa fa-clock mp-fg-neutral")
  expect_match(mp_state_titles("annotate")[["5"]], "Update on the Assemble step", fixed = TRUE)
  expect_equal(MP_STATE_CODES$assemble, c("0", "1", "4", "2", "3"))
  expect_equal(MP_STATE_SETTABLE, c("0", "1", "4", "2", "3"))
})
```

Append to the end of `tests/testthat/test-ui-reactable-helpers.R`:

```r

test_that("rt_gated shows the state of unfinished checks, else the wrapped cell", {
  js <- as.character(rt_gated(rt_flag_badge()))
  expect_match(js, "row['check_state']", fixed = TRUE)
  expect_match(js, "row['check_note']", fixed = TRUE)
  expect_match(js, "mp-pill mp-pill-info", fixed = TRUE)
  expect_match(js, as.character(rt_flag_badge()), fixed = TRUE)
  expect_match(as.character(rt_gated()), "String(v)", fixed = TRUE)
})

test_that("the ITSx badge shows a missing result as not checked and ticks only a clean roll-up", {
  js <- as.character(rt_itsx_badge())
  expect_match(js, "not checked", fixed = TRUE)
  re <- regmatches(js, regexpr("\\^0 of [^/]+\\$", js))
  expect_length(re, 1)
  expect_true(grepl(re, "0 of 4 contigs"))
  expect_false(any(grepl(re, c("0 of 4 contigs, 1 not checked", "1 of 4 contigs", "10 of 40 contigs"))))
})

test_that("the State glyph names the recheck's note only when asked to", {
  js <- as.character(rt_dynamicIcon(c(`3` = "fa fa-x"), labels = c(`3` = "Failed"), note_col = "check_note"))
  expect_match(js, "(cellInfo.row || {})['check_note']", fixed = TRUE)
  expect_false(grepl("note", as.character(rt_lock_icon("lock_row")), fixed = TRUE))
})
```

In `tests/testthat/test-app-scaffold-uid.R`, the roll-up test checked a copy of the roll-up SQL that no longer matches the app; read the real one through the Assemble table instead. Find:

```r
test_that("ITSx warnings reach the annotate rows, the scaffold view, and the sample roll-up", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_chimeric' WHERE path = 1 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = '' WHERE path = 1 AND scaffold = 2")
  DBI::dbDisconnect(con)
  local_mocked_bindings(app_project = function() dir, .package = "RiboPilot")
  ev <- RiboPilot:::unit_evidence("S1")
  expect_equal(ev$itsx_warnings[ev$path == 1 & ev$scaffold == 1], "itsx_chimeric")
  roll <- RiboPilot:::app_query(
    "SELECT (SELECT CASE WHEN COUNT(x.itsx_warnings) = 0 THEN NULL ELSE
        SUM(CASE WHEN COALESCE(x.itsx_warnings,'') <> ''
                 AND x.itsx_warnings NOT LIKE 'not_run(%' THEN 1 ELSE 0 END)
        || ' of ' || COUNT(x.itsx_warnings) || ' units' END
      FROM annotate x WHERE x.ID = s.ID) AS itsx_warnings FROM samples s")
  expect_equal(roll$itsx_warnings, "1 of 2 units")   # two checked, two not yet
})
```

Replace with:

```r
test_that("ITSx warnings reach the annotate rows, the scaffold view, and the sample roll-up", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_chimeric' WHERE path = 1 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = '' WHERE path = 1 AND scaffold = 2")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  ev <- RiboPilot:::unit_evidence("S1")
  expect_equal(ev$itsx_warnings[ev$path == 1 & ev$scaffold == 1], "itsx_chimeric")
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }
  js <- as.character(RiboPilot:::rt_itsx_badge())
  clean <- regmatches(js, regexpr("\\^0 of [^/]+\\$", js))
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    roll <- function(k) { session$setInputs(refresh = k); all_rows()$itsx_warnings }
    expect_equal(roll(1), "1 of 4 contigs, 2 not checked")
    set("UPDATE annotate SET itsx_warnings = 'not_run(ITSx did not finish)' WHERE path = 2 AND scaffold = 1")
    set("UPDATE annotate SET itsx_warnings = '' WHERE path = 2 AND scaffold = 2")
    expect_equal(roll(2), "1 of 4 contigs, 1 not checked")
    set("UPDATE annotate SET itsx_warnings = '' WHERE path = 1 AND scaffold = 1")
    expect_equal(roll(3), "0 of 4 contigs, 1 not checked")
    expect_false(grepl(clean, roll(4)))
    # a dropped contig does not count
    set("UPDATE assemblies SET ignore = 1 WHERE path = 2 AND scaffold = 1")
    expect_equal(roll(5), "0 of 3 contigs")
    expect_true(grepl(clean, roll(6)))
    set("UPDATE annotate SET itsx_warnings = NULL")
    expect_true(is.na(roll(7)))
  })
})
```

Append to the end of `tests/testthat/test-app-table.R`:

```r

test_that("Annotation Review reads the review layer and shows the state of unfinished checks", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  # p1.s1 done with an ITSx flag; p1.s2 still queued for the main pipeline;
  # p2.s1 a copy whose recheck failed; p2.s2 a rechecked copy
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, itsx_warnings = 'itsx_chimeric'
    WHERE path = 1 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, itsx_warnings = '' WHERE path = 2")
  copy <- function(s, len, sw, note, species) DBI::dbExecute(con, "INSERT INTO review_units (ID, path,
      scaffold, sequence, length, trim_start, trim_end, edited, check_switch, check_note, rev,
      n_features, itsx_warnings, blast_species)
    VALUES ('S1', 2, ?, 'ACGT', ?, 1, ?, 1, ?, ?, CAST(? AS INTEGER), 3, '', ?)",
    params = list(s, len, len, sw, note, "1760000000000", species))
  copy(1L, 500L, 3L, "tool crashed", "Copy sp")
  copy(2L, 800L, 2L, NA, "Danio rerio")
  # p2.s2's original screen and chimera results differ from its copy's
  DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 2, 2, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1', 2, 2, 'severe')")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 2, 2, 'on_target')")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    d <- all_rows()
    d <- d[order(d$path, d$scaffold), ]
    expect_equal(d$length, c(6120L, 1442L, 500L, 800L))
    expect_equal(d$blast_species, c("Danio rerio", "Danio rerio", "Copy sp", "Danio rerio"))
    expect_equal(d$screen_verdict, c(NA, NA, NA, "on_target"))
    expect_true(all(is.na(d$chimera_label)))
    expect_equal(d$annotate_switch, c(2L, 5L, 3L, 2L))
    expect_equal(d$check_state, c(NA, "waiting for Assemble Update", "failed", NA))
    expect_equal(d$check_note, c(NA, NA, "tool crashed", NA))
    expect_equal(d$checks, c("ITSx", "waiting for Assemble Update", "failed", "ok"))
    expect_equal(d$checks_tone, c("warning", "info", "info", "success"))

    session$setInputs(Annotate_state_filter = "5")
    expect_equal(unit_key(samples_data()$ID, samples_data()$path, samples_data()$scaffold),
      unit_key("S1", 1L, 2L))
    j <- jsonlite::fromJSON(as.character(output$tbl), simplifyVector = FALSE)
    cols <- j$x$tag$attribs$columns
    names(cols) <- vapply(cols, function(c) c$id %||% "", "")
    expect_false(cols$check_state$show)
    expect_false(cols$check_note$show)
    for (k in c("n_flags", "itsx_warnings", "screen_verdict", "chimera_label", "blast_species",
                "blast_accession", "blast_pident", "blast_qcovs", "blast_source")) {
      expect_match(cols[[k]]$cell, "row['check_state']", fixed = TRUE, info = k)
      expect_true(isTRUE(cols[[k]]$html), info = k)
    }
    expect_match(cols$annotate_switch$cell, "(cellInfo.row || {})['check_note']", fixed = TRUE)
    expect_match(cols$annotate_switch$cell, "Waiting for Assemble Update", fixed = TRUE)
  })
})
```

In `tests/testthat/test-app-opts.R`, the anchor on the Annotate query moves with its first line (the State now comes from `review_units_v`); the 700 characters after it must still carry no option sets. Find:

```r
  i <- regexpr("an.annotate_lock, an.annotate_switch", srv, fixed = TRUE)
```

Replace with:

```r
  i <- regexpr("an.annotate_lock, COALESCE(v.check_switch, 1) AS annotate_switch", srv, fixed = TRUE)
```

- [ ] **Step 6: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("ui-helpers", "ui-reactable-helpers", "app-scaffold-uid", "app-table", "app-opts")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-ui-helpers.R `[ FAIL 3 | WARN 0 | SKIP 0 | PASS 109 ]`: the `mp_state_choices("annotate")` names lack "Waiting for Assemble Update"; the new test finds no `value="5"` and errors with `subscript out of bounds` on `mp_state_icon_classes("annotate")[["5"]]`. test-ui-reactable-helpers.R `[ FAIL 5 | WARN 0 | SKIP 0 | PASS 35 ]`: `could not find function "rt_gated"`, `expect_length(re, 1)` fails and `grepl()` errors with `invalid 'pattern' argument` (the badge still tests `/^0 of /`), and `unused argument (note_col = "check_note")`. test-app-scaffold-uid.R `[ FAIL 4 | WARN 0 | SKIP 0 | PASS 122 ]`: the roll-up reads `"1 of 2 contigs"`, not `"1 of 4 contigs, 2 not checked"`. test-app-table.R: the new test stops after 10 failures, starting with `d$length` (`5980`, `890` from the originals), `d$blast_species`, `d$annotate_switch` (`2 1 2 2`), and `d$check_state` (NULL). test-app-opts.R `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 62 ]`: `expect_true(i > 0)` fails (the new anchor is not in the code yet).

- [ ] **Step 7: Add state 5, the gated cells, and the review-layer query**

In `R/constants.R`, find:

```r
MP_STATE_META <- list(
  `0` = list(label = "On hold",      icon = "fa fa-hourglass",
```

Replace with:

```r
MP_STATE_META <- list(
  `5` = list(label = "Waiting for Assemble Update", icon = "fa fa-clock",
             tip = "The main pipeline has not finished this contig. Press Update on the Assemble step first."),
  `0` = list(label = "On hold",      icon = "fa fa-hourglass",
```

Find:

```r
#' Codes each module can display.
#' @noRd
MP_STATE_CODES <- list(
  taxonomy = c("0", "1", "4", "2", "3"),
  assemble = c("0", "1", "4", "2", "3"),
  annotate = c("0", "1", "4", "2", "3")
)
```

Replace with:

```r
#' Codes each module can display. Annotation Review's 5 is display-only.
#' @noRd
MP_STATE_CODES <- list(
  taxonomy = c("0", "1", "4", "2", "3"),
  assemble = c("0", "1", "4", "2", "3"),
  annotate = c("5", "0", "1", "4", "2", "3")
)
```

In `R/app_ui_helpers.R`, find:

```r
  tone <- c(`0` = "mp-fg-neutral", `1` = "mp-fg-neutral", `4` = "mp-fg-neutral",
            `2` = "mp-fg-success", `3` = "mp-fg-danger")
```

Replace with:

```r
  tone <- c(`5` = "mp-fg-neutral", `0` = "mp-fg-neutral", `1` = "mp-fg-neutral", `4` = "mp-fg-neutral",
            `2` = "mp-fg-success", `3` = "mp-fg-danger")
```

In `R/app_ui.R`, find:

```r
#' assemble/annotate lock (0/1) and switch (0/1/2/3/4) codes. Each mode gets its
```

Replace with:

```r
#' assemble/annotate lock (0/1) and switch (0/1/2/3/4) codes, plus 5 (waiting for
#' the Assemble step's Update) on Annotation Review. Each mode gets its
```

Find:

```r
filter_pickers <- function(mode) {
  shiny::tagList(
    mp_filter_picker(paste0(mode, "_lock_filter"), "Lock:", LOCK_CHOICES, width = "140px"),
    mp_filter_picker(paste0(mode, "_state_filter"), "State:", STATE_CHOICES, width = "140px")
  )
}
```

Replace with:

```r
filter_pickers <- function(mode) {
  states <- if (mode == "Annotate") mp_state_choices("annotate") else STATE_CHOICES
  shiny::tagList(
    mp_filter_picker(paste0(mode, "_lock_filter"), "Lock:", LOCK_CHOICES, width = "140px"),
    mp_filter_picker(paste0(mode, "_state_filter"), "State:", states, width = "140px")
  )
}
```

In `R/utils_reactable.R`, find:

```r
#' @param inputId optional Shiny input id. When supplied the icon is a button
#'   that sends its 1-based row index to this input.
#'
#' @noRd
rt_dynamicIcon <- function(icons = NULL, labels = NULL, inputId = NULL) {
```

Replace with:

```r
#' @param inputId optional Shiny input id. When supplied the icon is a button
#'   that sends its 1-based row index to this input.
#' @param note_col optional column of the row; its text, when the row has one,
#'   is added to the glyph's name.
#'
#' @noRd
rt_dynamicIcon <- function(icons = NULL, labels = NULL, inputId = NULL, note_col = NULL) {
  note <- if (is.null(note_col)) "" else sprintf(
    "var note = (cellInfo.row || {})['%s'];
      if (label && note) label += ' Note: ' + String(note).replace(/&/g, '&amp;')
        .replace(/</g, '&lt;').replace(/'/g, '&#39;');", note_col)
```

Find:

```r
      var label = labels[value];
      var name = label ?
```

Replace with:

```r
      var label = labels[value];
      %s
      var name = label ?
```

Find:

```r
    mp_js_obj(icons), mp_js_obj(labels), inputId %||% ""
```

Replace with:

```r
    mp_js_obj(icons), mp_js_obj(labels), note, inputId %||% ""
```

Find:

```r
#' Render the stock-ITSx check: tick when clean, grey when not run, amber pill
#' with the flag list otherwise (or the 'n of m contigs' roll-up on a sample row)
#' @noRd
rt_itsx_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined || v === '') return ''
      v = String(v)
      var esc = v.replace(/\"/g, '&quot;')
      if (v.indexOf('not_run(') === 0) return '<span title=\"' + esc + '\" style=\"color:#999;\">not run</span>'
      if (/^0 of /.test(v)) return '<span class=\"mp-fg-success\">&#10003;</span>'
      var txt = v.split(';').join('; ')
      return '<span title=\"' + esc + '\" class=\"mp-pill mp-pill-warning\">' + txt + '</span>'
    }")
}
```

Replace with:

```r
#' Render the stock-ITSx check: grey when never checked or not run, blank when
#' clean, amber pill with the flag list otherwise. A sample row's roll-up ticks
#' only when it reads exactly '0 of <m> contigs'.
#' @noRd
rt_itsx_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value
      if (v === null || v === undefined) return '<span style=\"color:#999;\">not checked</span>'
      v = String(v)
      if (v === '') return ''
      var esc = v.replace(/\"/g, '&quot;')
      if (v.indexOf('not_run(') === 0) return '<span title=\"' + esc + '\" style=\"color:#999;\">not run</span>'
      if (/^0 of [0-9]+ contigs$/.test(v)) return '<span class=\"mp-fg-success\">&#10003;</span>'
      var txt = v.split(';').join('; ')
      return '<span title=\"' + esc + '\" class=\"mp-pill mp-pill-warning\">' + txt + '</span>'
    }")
}

#' Show a check's state instead of its result while the checks are not done
#'
#' A row with a `check_state` ("not run", "running", "failed", "on hold", or
#' "waiting for Assemble Update") gets an info pill, titled with the row's
#' `check_note` when it has one. Other rows render through `cell`, or show the
#' plain value when `cell` is NULL.
#' @noRd
rt_gated <- function(cell = NULL) {
  inner <- if (is.null(cell)) {
    "function(cellInfo) { var v = cellInfo.value; return (v === null || v === undefined) ? '' : String(v); }"
  } else {
    as.character(cell)
  }
  htmlwidgets::JS(sprintf(
    "function(cellInfo) {
      var row = cellInfo.row || {};
      var st = row['check_state'];
      if (st) {
        var note = row['check_note'];
        var tip = note ? ' title=\"' + String(note).replace(/&/g, '&amp;').replace(/\"/g, '&quot;')
          .replace(/</g, '&lt;') + '\"' : '';
        return '<span class=\"mp-pill mp-pill-info\"' + tip + '>' + st + '</span>';
      }
      return (%s)(cellInfo);
    }", inner))
}
```

In `R/app_server.R`, the Assemble roll-up. Find:

```r
           -- per-sample roll-up of the stock-ITSx check: NULL until any unit
           -- has been checked, else 'n of m contigs' flagged
           (SELECT CASE WHEN COUNT(x.itsx_warnings) = 0 THEN NULL ELSE
              SUM(CASE WHEN COALESCE(x.itsx_warnings,'') <> ''
                       AND x.itsx_warnings NOT LIKE 'not_run(%' THEN 1 ELSE 0 END)
              || ' of ' || COUNT(x.itsx_warnings) || ' contigs' END
            FROM annotate x WHERE x.ID = s.ID) AS itsx_warnings,
```

Replace with:

```r
           -- per-sample roll-up of the stock-ITSx check over kept contigs: NULL
           -- until any is checked, else 'n of m contigs' flagged, plus how many
           -- have no result or a not_run one, so it never reads as clean then
           (SELECT CASE WHEN COUNT(x.itsx_warnings) = 0 THEN NULL ELSE
              SUM(CASE WHEN COALESCE(x.itsx_warnings,'') <> ''
                       AND x.itsx_warnings NOT LIKE 'not_run(%' THEN 1 ELSE 0 END)
              || ' of ' || COUNT(*) || ' contigs' ||
              CASE WHEN SUM(x.itsx_warnings IS NULL OR x.itsx_warnings LIKE 'not_run(%') = 0 THEN ''
                ELSE ', ' || SUM(x.itsx_warnings IS NULL OR x.itsx_warnings LIKE 'not_run(%') ||
                  ' not checked' END END
            FROM annotate x
              JOIN assemblies y ON y.ID = x.ID AND y.path = x.path AND y.scaffold = x.scaffold
            WHERE x.ID = s.ID AND COALESCE(y.ignore,0) = 0) AS itsx_warnings,
```

The Annotate branch. Find:

```r
      # annotate shows one row per non-ignored assembly UNIT of an assembly-locked
      # sample (one (ID, path, scaffold) each, annotated independently)
      app_query(paste(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ",
           s.Taxon, s.sci_name, s.taxid,
           asm.length, an.n_features, an.missing,
           -- blast_source says WHICH search named the species: the local seed
           -- BLAST, the local panel BLAST, or the remote core_nt one. Without
           -- it a seed or panel hit reads as a GenBank identification.
           -- blast_qcovs is written by the remote and local-panel paths only
           -- and its column was already defined but never selected.
           asm.blast_species, asm.blast_accession, asm.blast_pident, asm.blast_qcovs,
           asm.blast_source,
           an.n_flags, an.curate_flags, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason, an.reviewed,
           an.problematic, an.partial, an.annotate_notes,
           an.export_group,", UNIT_EXPORTED_SQL, ", an.time_stamp,
           'details' AS details, 'output' AS output
         FROM annotate an
           JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0")) |> with_checks("unit")
    }
```

Replace with:

```r
      # annotate shows one row per non-ignored assembly UNIT of an assembly-locked
      # sample (one (ID, path, scaffold) each, annotated independently). Contig
      # data and check results come from the review layer (the contig's review
      # copy, else its original); the marks, notes, and group from annotate.
      d <- app_query(paste(
        "SELECT an.annotate_lock, COALESCE(v.check_switch, 1) AS annotate_switch, s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ",
           s.Taxon, s.sci_name, s.taxid,
           v.length, v.n_features, v.missing,
           -- blast_source says WHICH search named the species: the local seed
           -- BLAST, the local panel BLAST, or the remote core_nt one. Without
           -- it a seed or panel hit reads as a GenBank identification.
           -- blast_qcovs is written by the remote and local-panel paths only
           -- and its column was already defined but never selected.
           v.blast_species, v.blast_accession, v.blast_pident, v.blast_qcovs,
           v.blast_source,
           v.n_flags, v.curate_flags, v.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason, an.reviewed,
           an.problematic, an.partial, an.annotate_notes,
           an.export_group,", UNIT_EXPORTED_SQL, ", an.time_stamp,
           'details' AS details, 'output' AS output
         FROM annotate an
           JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
           LEFT JOIN review_screen_v sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN review_chimera_v ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0")) |> with_checks("unit", review = TRUE)
      # The State column shows the checks' state, 5 for a contig waiting for the
      # Assemble step's Update; check_state gates the check cells (rt_gated).
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      st <- review_unit_state(con, d[c("ID", "path", "scaffold")])
      d$annotate_switch[st$state %in% "waiting"] <- 5L
      d$check_state <- ifelse(st$state %in% "done", NA_character_, review_state_text(st$state))
      d$check_note <- st$check_note
      d
    }
```

Find:

```r
  with_checks <- function(d, by) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ck <- unit_checks(con, by)
```

Replace with:

```r
  with_checks <- function(d, by, review = FALSE) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ck <- unit_checks(con, by, review = review)
```

The column definitions in `output$tbl`. Find:

```r
      checks_tone = reactable::colDef(show = FALSE),
```

Replace with:

```r
      checks_tone = reactable::colDef(show = FALSE),
      check_state = reactable::colDef(show = FALSE),
      check_note = reactable::colDef(show = FALSE),
```

Find:

```r
      n_flags = cd("n_flags", width = 100, html = TRUE, filterable = FALSE, align = "center",
        cell = rt_flag_badge()),
```

Replace with:

```r
      n_flags = cd("n_flags", width = 100, html = TRUE, filterable = FALSE, align = "center",
        cell = rt_gated(rt_flag_badge())),
```

Find:

```r
      itsx_warnings = cd("itsx_warnings", minWidth = 130, html = TRUE, cell = rt_itsx_badge()),
      screen_verdict = cd("screen_verdict", width = 110, html = TRUE, align = "center",
        cell = rt_screen_badge()),
      chimera_label = cd("chimera_label", width = 100, html = TRUE, align = "center",
        cell = rt_chimera_badge()),
```

Replace with:

```r
      itsx_warnings = cd("itsx_warnings", minWidth = 130, html = TRUE, cell = rt_gated(rt_itsx_badge())),
      screen_verdict = cd("screen_verdict", width = 110, html = TRUE, align = "center",
        cell = rt_gated(rt_screen_badge())),
      chimera_label = cd("chimera_label", width = 100, html = TRUE, align = "center",
        cell = rt_gated(rt_chimera_badge())),
```

Find:

```r
      blast_species = cd("blast_species", minWidth = 160, html = TRUE, cell = rt_longtext()),
      blast_accession = cd("blast_accession", minWidth = 130, html = TRUE, cell = rt_ncbi_link()),
      blast_pident = cd("blast_pident", width = 130),
      blast_qcovs = cd("blast_qcovs", width = 140),
      blast_source = cd("blast_source", width = 100, html = TRUE, align = "center",
        cell = rt_blast_source_badge()),
```

Replace with:

```r
      blast_species = cd("blast_species", minWidth = 160, html = TRUE, cell = rt_gated(rt_longtext())),
      blast_accession = cd("blast_accession", minWidth = 130, html = TRUE, cell = rt_gated(rt_ncbi_link())),
      blast_pident = cd("blast_pident", width = 130, html = TRUE, cell = rt_gated()),
      blast_qcovs = cd("blast_qcovs", width = 140, html = TRUE, cell = rt_gated()),
      blast_source = cd("blast_source", width = 100, html = TRUE, align = "center",
        cell = rt_gated(rt_blast_source_badge())),
```

Find:

```r
        cell = rt_dynamicIcon(icons = mp_state_icon_classes(mm$pfx), labels = mp_state_titles(mm$pfx),
          inputId = "state_row"))
```

Replace with:

```r
        cell = rt_dynamicIcon(icons = mp_state_icon_classes(mm$pfx), labels = mp_state_titles(mm$pfx),
          inputId = "state_row", note_col = "check_note"))
```

Rows without a `check_state` (the Assemble and Export steps) render exactly as before. Check the files parse: `Rscript -e 'for (f in c("R/constants.R", "R/app_ui_helpers.R", "R/app_ui.R", "R/utils_reactable.R", "R/app_server.R")) invisible(parse(f))'`.

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-state", "export", "ui-helpers", "ui-reactable-helpers", "app-scaffold-uid", "app-table", "app-opts", "app-export")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS, every file `FAIL 0` (PASS 16, 81, 115, 46, 129, 79, 63, 142).

Then every other file that renders or reads the app tables, the views, or the state registry:
Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "notes-partial", "specimen-app", "batch-trim", "blast-provenance", "seqview", "app-run", "app-utils", "app-headless", "docs-site", "no-duplicate-definitions")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0` (PASS 127, 9, 99, 51, 23, 17, 31, 26, 14, 237, 1). On a heavily loaded machine the two `open_path` timing tests in test-app-utils.R can time out; rerun that file alone.

- [ ] **Step 9: Commit**

```bash
git add R/review.R R/export.R R/app_server.R R/utils_reactable.R R/constants.R R/app_ui_helpers.R R/app_ui.R tests/testthat/test-review-state.R tests/testthat/test-export.R tests/testthat/test-app-export.R tests/testthat/test-app-opts.R tests/testthat/test-app-scaffold-uid.R tests/testthat/test-ui-helpers.R tests/testthat/test-ui-reactable-helpers.R tests/testthat/test-app-table.R
git commit -m "Annotation Review table reads the review layer; check state column"
```


### Task 15: Trim and Restore on copies

**Files:**
- Modify: `R/export.R` (`unit_trimmed()`, repo :225-238)
- Modify: `R/pipeline_steps.R` (`batch_trim_plan()` and `batch_trim_apply()`, repo :775-827; delete `trim_assembly()` :562-641, `ensure_assembly_backup()` :661-665, `restore_assembly()` :719-764; keep `write_coverage_summary()` and `trim_bounds()`)
- Modify: `R/chimera.R` (delete `chimera_snapshot()`, repo :588-598; roxygen title of `chimera_restore()` :600)
- Modify: `R/app_server.R` (batch Trim / Restore block, repo :905-948; trim and restore handlers and the trim window, repo :2873-2997)
- Modify: `R/app_modals.R` (`trim_modal()` help text :656-658; `batch_trim_modal()` title :673; `batch_trim_summary()` :689-699; `plot_coverage()` :908-934)
- Modify: `R/constants.R` (`MP_TRIM_TIP()`, :402-407)
- Modify: `tests/testthat/test-export.R` (:57-78), `tests/testthat/test-batch-trim.R`, `tests/testthat/test-app-export.R` (trim parts, :562-632), `tests/testthat/test-chimera.R` (delete :433-456)
- Rewrite: `tests/testthat/test-trim-coverage.R`
- Test: `tests/testthat/test-batch-trim.R`, `tests/testthat/test-export.R`, `tests/testthat/test-app-export.R`, `tests/testthat/test-trim-coverage.R`

Line numbers are those of the repo before this plan; earlier tasks shift some of them by a few lines, so find each block by its text.

**Interfaces:**
- Consumes:
  - Task 1: views `review_units_v` (columns used here: `has_copy, orig_length, sequence, length, trim_start, trim_end, check_switch`), `review_annotations_v`, `review_chimera_v`; table `review_units`; test helpers `fixture_review_schema(con)`, `mk_app_project()`, `mk_export_project()` carrying the review schema; `chimera_test_con()` builds a real `new_db()`.
  - Task 2: `review_has_copy(con, id, path, scaffold)` (logical(1)), `review_copy(con, id, path, scaffold)`, `review_drop(con, id, path, scaffold)` (one savepoint deleting the unit's rows in all 8 review tables), `REVIEW_EFF_PARTIALS`.
  - Task 3: `review_coverage(con, id, path, scaffold)` -> `list(depths, clips)`: the original coverage track, cut to the copy's window when there is a copy; `NULL` for a missing row, track, or empty cut.
  - Task 4: `review_trim(con, id, path, scaffold, from, to)` -> invisible new length. `from`/`to` are copy coordinates; it makes the copy on first touch, cuts the sequence from the original by the window, shifts and clamps `review_annotations` (cut ends get `partial_*_user = 1`), sets `edited = 1`, and invalidates (check_switch 1). Errors: `"No sequence for contig <ID>.<path>.<scaffold>."` (from `review_copy()` when the unit has no `review_units_v` row, or from `review_trim()` when the copy has no sequence) and `"Invalid trim range."`.
- Produces:
  - `batch_trim_plan(con, u, flank = 0L, restore = FALSE)`: the `running` argument is gone. Reads `review_units_v` / `review_annotations_v`. Reasons: `"annotation locked"`; restore: `"no edits to restore"` for a unit without a copy (every unit with a copy qualifies, trimmed or not); trim: `"no sequence"`, `"missing <elements>"`, `"nothing to trim"`. `from`/`to` are copy coordinates.
  - `batch_trim_apply(con, p, restore = FALSE)`: `review_drop()` (restore) or `review_trim()` (trim) per row, no `curate_one()`; error message per unit, `NA` where it worked.
  - `unit_trimmed(con, id, path, scaffold)`: `c(isTRUE(trim_start > 1), isTRUE(trim_end < orig_length))` from `review_units_v`; `c(FALSE, FALSE)` when the unit has no row or the database has no review layer.
  - `plot_coverage(id, path = 1L, scaffold = 1L)` and `plot_trim(id, path = 1L, scaffold = 1L, from = NA, to = NA)` (signatures unchanged) draw the review layer: track from `review_coverage()`, features from `review_annotations_v`, in copy coordinates.
  - In `app_server()`: `end_edit_session()` (session-local helper) clears `edit_active`, `edit_backup`, `edit_dirty` and unticks Edit; called after every trim and restore.
  - Removed: `trim_assembly()`, `restore_assembly()`, `ensure_assembly_backup()`, `chimera_snapshot()`. Kept for Task 22: `write_coverage_summary()`, `chimera_shift()`, `chimera_restore()` (direct tests in test-trim-coverage.R).

Trim (single and batch) and Restore now write only the contig's review copy (SPEC S5.2 trim and restore rows; S4.1: originals stay as the pipeline left them for the Assemble step). A trim makes the copy on first touch through `review_trim()`, which also marks the copy's checks not run; Restore calls `review_drop()`, which deletes the copy with its trims, boundary edits, and 5'/3' choices while the marks, notes, group, and lock in `annotate` stay (D5). Batch Restore offers every contig with a copy, not only trimmed ones (S5.2), and the batch "in the active run" skip goes (S5.4). The inline `curate_one()` calls after a trim or restore go (D6); the window text says the checks show not run until Update. `unit_trimmed()` now reads the kept window instead of searching the backup (S9), and the trim window plans, plots, and validates from the copy. The old in-place trim functions lose their last callers and are deleted; the coverage and chimera helpers the Task 22 upgrade still needs keep direct tests. Edit mode, the 5'/3' toggles, Save, Reset, and the feature table move to the review layer in Task 16, so this task's tests observe only the trim window, the views, and the original tables.

- [ ] **Step 1: Write the failing batch and window tests**

In `tests/testthat/test-export.R`, the trimmed-end test built its trimmed state from `assembly_backup`; it now checks the window on the review layer.

Find:

```r
test_that("a scaffold end marks a gene partial only on an untrimmed unit", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT,
    start INT, stop INT, strand TEXT, partial_start INT, partial_stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, RiboPilot:::ASSEMBLY_BACKUP_DDL)
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 300)", params = list(strrep("A", 300)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',1,300,'+',NULL,NULL)")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t>300\trRNA")
  # trimmed: the snapshot holds a longer sequence, so the ends are cut, not assembly ends
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence) VALUES ('S1',1,1,?)",
    params = list(paste0(strrep("C", 50), strrep("A", 300), strrep("G", 50))))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n1\t300\trRNA")
  # only the 3' tail cut: base 1 is still an assembly end
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = ?",
    params = list(paste0(strrep("A", 300), strrep("G", 50))))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t300\trRNA")
  # restored: snapshot and sequence agree again
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = ?", params = list(strrep("A", 300)))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t>300\trRNA")
})
```

Replace with:

```r
test_that("unit_trimmed reads the kept window of the review copy", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence) VALUES ('S1', 1, 1, 400, ?)",
    params = list(strrep("A", 400)))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, strand)
    VALUES ('S1', 1, 1, '18S', 'rRNA', 1, 400, '+')")
  tr <- function() RiboPilot:::unit_trimmed(con, "S1", 1L, 1L)
  # no copy: both ends are assembly ends
  expect_equal(tr(), c(FALSE, FALSE))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t>400\trRNA")
  expect_equal(RiboPilot:::unit_trimmed(con, "S9", 1L, 1L), c(FALSE, FALSE))

  win <- function(a, b) {
    DBI::dbExecute(con, "DELETE FROM review_units")
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, trim_start, trim_end, rev)
      VALUES ('S1', 1, 1, ?, ?, ?, ?, 1)", params = list(strrep("A", b - a + 1L), b - a + 1L, a, b))
    tr()
  }
  expect_equal(win(51L, 350L), c(TRUE, TRUE))
  expect_equal(win(1L, 350L), c(FALSE, TRUE))
  expect_equal(win(51L, 400L), c(TRUE, FALSE))
  expect_equal(win(1L, 400L), c(FALSE, FALSE))

  # a database without the review layer reads as untrimmed
  mem <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(mem), add = TRUE)
  expect_equal(RiboPilot:::unit_trimmed(mem, "S1", 1L, 1L), c(FALSE, FALSE))
})
```

In `tests/testthat/test-batch-trim.R`, every reader moves to the review layer, Restore reasons change, the run skip goes, and the tests that drove `trim_assembly()`/`restore_assembly()` move to the copy functions.

Find:

```r
len <- function(con) DBI::dbGetQuery(con, "SELECT length FROM assemblies ORDER BY path, scaffold")$length
```

Replace with:

```r
len <- function(con) DBI::dbGetQuery(con, "SELECT length FROM review_units_v ORDER BY path, scaffold")$length
```

Find:

```r
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)$reason,
    c("not trimmed", "not trimmed", "annotation locked", "not trimmed"))
})

test_that("a run skips queued units", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2")
  expect_equal(RiboPilot:::batch_trim_plan(con, units4)$reason, c(NA, "missing ITS2, 28S", "annotation locked", NA))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE path = 1 AND scaffold = 1")
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, running = TRUE)$reason,
    c("in the active run", "missing ITS2, 28S", "annotation locked", NA))
  expect_true(is.na(RiboPilot:::batch_trim_plan(con, units4)$reason[1]))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 4")
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, running = TRUE)$reason[4], "in the active run")
})

test_that("batch trim cuts to gene ends or a flank, sets no partial ends, and restores", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ok <- function(p) p[is.na(p$reason), , drop = FALSE]

  err <- RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(err, c(NA_character_, NA_character_))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assembly_backup")$n, 2L)
  expect_false(any(DBI::dbGetQuery(con,
    "SELECT COALESCE(partial_start, 0) + COALESCE(partial_stop, 0) n FROM annotations")$n > 0))

  # trimmed again with a 20 bp flank: the unit already sits at its gene ends,
  # so the flank is clipped to the scaffold and nothing changes
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, flank = 20)$reason[1], "nothing to trim")

  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c(NA, "not trimmed", "annotation locked", NA))
  RiboPilot:::batch_trim_apply(con, ok(p), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))

  # flank mode, then gene ends: restore still returns the original
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, flank = 20)))
  expect_equal(len(con), c(540L, 600L, 600L, 540L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT start, stop FROM annotations WHERE path = 1 AND scaffold = 1 AND element = '18S'"),
    data.frame(start = 51L, stop = 200L))
})
```

Replace with:

```r
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)$reason,
    c("no edits to restore", "no edits to restore", "annotation locked", "no edits to restore"))
})

test_that("contigs queued for a run are trimmed like any other", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 4")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 4")
  expect_equal(RiboPilot:::batch_trim_plan(con, units4)$reason, c(NA, "missing ITS2, 28S", "annotation locked", NA))
  expect_false("running" %in% names(formals(RiboPilot:::batch_trim_plan)))
})

test_that("batch trim cuts review copies to gene ends or a flank, and restore deletes them", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ok <- function(p) p[is.na(p$reason), , drop = FALSE]
  orig <- function() list(
    DBI::dbGetQuery(con, "SELECT path, scaffold, length, sequence FROM assemblies ORDER BY path, scaffold"),
    DBI::dbGetQuery(con, "SELECT path, scaffold, element, start, stop, partial_start, partial_stop
      FROM annotations ORDER BY path, scaffold, start"))
  before <- orig()

  err <- RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(err, c(NA_character_, NA_character_))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT path, scaffold, trim_start, trim_end, edited, check_switch
    FROM review_units ORDER BY path, scaffold"),
    data.frame(path = c(1L, 2L), scaffold = c(1L, 2L), trim_start = 51L, trim_end = 550L,
      edited = 1L, check_switch = 1L))
  # a cut exactly at the gene ends leaves no partial end
  expect_false(any(DBI::dbGetQuery(con,
    "SELECT COALESCE(partial_start_user, 0) + COALESCE(partial_stop_user, 0) n FROM review_annotations")$n > 0))
  expect_equal(orig(), before)

  # trimmed again with a 20 bp flank: the copy already sits at its gene ends,
  # so the flank is clipped to the contig and nothing changes
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, flank = 20)$reason[1], "nothing to trim")

  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c(NA, "no edits to restore", "annotation locked", NA))
  expect_equal(RiboPilot:::batch_trim_apply(con, ok(p), restore = TRUE), c(NA_character_, NA_character_))
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_units")$n, 0L)

  # flank mode, then gene ends: the second trim plans from the copy, in copy coordinates
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, flank = 20)))
  expect_equal(len(con), c(540L, 600L, 600L, 540L))
  p <- RiboPilot:::batch_trim_plan(con, units4)
  expect_equal(c(p$from[1], p$to[1]), c(21L, 520L))
  RiboPilot:::batch_trim_apply(con, ok(p))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT trim_start, trim_end FROM review_units WHERE path = 1 AND scaffold = 1"),
    data.frame(trim_start = 51L, trim_end = 550L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT start, stop FROM review_annotations_v WHERE path = 1 AND scaffold = 1 AND element = '18S'"),
    data.frame(start = 51L, stop = 200L))
  expect_equal(orig(), before)
})
```

Find:

```r
  lens <- function() RiboPilot:::app_query("SELECT length FROM assemblies ORDER BY path, scaffold")$length

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = seq_len(nrow(samples_data())), annot_trim_btn = 1)
    html <- output$batch_trim_summary$html
    expect_match(html, "2 contigs are trimmed", fixed = TRUE)
```

Replace with:

```r
  lens <- function() RiboPilot:::app_query("SELECT length FROM review_units_v ORDER BY path, scaffold")$length
  tip <- RiboPilot:::MP_TRIM_TIP()
  expect_false(grepl("running run", tip[["trim"]], fixed = TRUE))
  expect_match(tip[["restore"]], "boundary edits", fixed = TRUE)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = seq_len(nrow(samples_data())), annot_trim_btn = 1)
    html <- output$batch_trim_summary$html
    expect_match(html, "2 contigs are trimmed; their checks show not run until you press Update.", fixed = TRUE)
    expect_false(grepl("re-curated", html, fixed = TRUE))
```

Find:

```r
    session$setInputs(annot_restore_btn = 1)
    expect_match(output$batch_trim_summary$html, "S1 p1 s2: not trimmed", fixed = TRUE)
```

Replace with:

```r
    session$setInputs(annot_restore_btn = 1)
    expect_match(output$batch_trim_summary$html, "S1 p1 s2: no edits to restore", fixed = TRUE)
    expect_match(output$batch_trim_summary$html, "boundary edits, and 5'/3' choices are deleted", fixed = TRUE)
```

Find:

```r
test_that("only the trimmed end of a unit stops being an assembly end", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotations SET start = 1 WHERE path = 1 AND scaffold = 1 AND element = '18S'")
  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, 1, 550)
  tr <- RiboPilot:::unit_trimmed(con, "S1", 1L, 1L)
  expect_equal(tr, c(FALSE, TRUE))
  d <- DBI::dbGetQuery(con, "SELECT * FROM annotations WHERE path = 1 AND scaffold = 1 ORDER BY start")
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p5[1])
  expect_false(e$p3[5])
  # on the minus strand the untouched left end is the 3' end
  d$strand <- "-"
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p3[1])
  expect_false(e$p5[5])
  RiboPilot:::restore_assembly(con, "S1", 1L, 1L)
  expect_equal(RiboPilot:::unit_trimmed(con, "S1", 1L, 1L), c(FALSE, FALSE))
})

test_that("a trim after a restore snapshots the edits made in between", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  el18 <- function() DBI::dbGetQuery(con,
    "SELECT start FROM annotations WHERE path = 1 AND scaffold = 1 AND element = '18S'")$start
  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, 51, 550)
  RiboPilot:::restore_assembly(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE annotations SET start = 41 WHERE path = 1 AND scaffold = 1 AND element = '18S'")
  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, 41, 550)
  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, 2, 500)
  RiboPilot:::restore_assembly(con, "S1", 1L, 1L)
  expect_equal(el18(), 41L)
  expect_equal(len(con)[1], 600L)
})
```

Replace with:

```r
test_that("only the trimmed end of a copy stops being an assembly end", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotations SET start = 1 WHERE path = 1 AND scaffold = 1 AND element = '18S'")
  expect_equal(RiboPilot:::unit_trimmed(con, "S1", 1L, 1L), c(FALSE, FALSE))
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1, 550)
  tr <- RiboPilot:::unit_trimmed(con, "S1", 1L, 1L)
  expect_equal(tr, c(FALSE, TRUE))
  d <- DBI::dbGetQuery(con, paste("SELECT element, start, stop, strand,", RiboPilot:::REVIEW_EFF_PARTIALS,
    "FROM review_annotations_v WHERE path = 1 AND scaffold = 1 ORDER BY start"))
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p5[1])
  expect_false(e$p3[5])
  # on the minus strand the untouched left end is the 3' end
  d$strand <- "-"
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p3[1])
  expect_false(e$p5[5])
  RiboPilot:::batch_trim_apply(con, units4[1, ], restore = TRUE)
  expect_equal(RiboPilot:::unit_trimmed(con, "S1", 1L, 1L), c(FALSE, FALSE))
})

test_that("batch Restore offers every contig with a copy, not only trimmed ones", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  # an untrimmed copy, as State makes one
  RiboPilot:::review_copy(con, "S1", 2L, 2L)
  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c("no edits to restore", "no edits to restore", "annotation locked", NA))
  expect_equal(RiboPilot:::batch_trim_apply(con, p[4, ], restore = TRUE), NA_character_)
  expect_false(RiboPilot:::review_has_copy(con, "S1", 2L, 2L))
})
```

Find:

```r
  db <- function() RiboPilot:::app_query(
    "SELECT element, start, stop FROM annotations WHERE path = 2 AND scaffold = 2 ORDER BY start")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- which(d$path == 2 & d$scaffold == 2)
    session$setInputs(ann_details = row)
    expect_equal(feature_rows()$start[1], 51L)
    session$setInputs(annotate_close = 1)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1)

    session$setInputs(batch_trim_mode = "flank", batch_trim_flank = NA)
    expect_match(output$batch_trim_summary$html, "Enter a flank of at least 1 bp.", fixed = TRUE)
    session$setInputs(batch_trim_go = 1)
    expect_equal(db()$start[1], 51L)

    session$setInputs(batch_trim_mode = "annotated", batch_trim_go = 2)
    session$setInputs(ann_details = row)
    f <- feature_rows()
    expect_equal(f$start, db()$start)
    expect_equal(f$stop, db()$stop)
    expect_equal(unname(trim_span()), c(1, 500, 500))
  })
})
```

Replace with:

```r
  db <- function() RiboPilot:::app_query(
    "SELECT element, start, stop FROM review_annotations_v WHERE path = 2 AND scaffold = 2 ORDER BY start")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- which(d$path == 2 & d$scaffold == 2)
    session$setInputs(ann_details = row)
    expect_equal(feature_rows()$start[1], 51L)
    expect_equal(unname(trim_span()), c(51, 550, 600))
    # an Edit session left open when the window closes
    session$setInputs(feat_edit = TRUE)
    expect_false(is.null(edit_backup()))
    session$setInputs(annotate_close = 1)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1)

    session$setInputs(batch_trim_mode = "flank", batch_trim_flank = NA)
    expect_match(output$batch_trim_summary$html, "Enter a flank of at least 1 bp.", fixed = TRUE)
    session$setInputs(batch_trim_go = 1)
    expect_equal(db()$start[1], 51L)
    expect_false(is.null(edit_backup()))

    # the trim ends that session, so a later Reset cannot write pre-trim coordinates
    session$setInputs(batch_trim_mode = "annotated", batch_trim_go = 2)
    expect_equal(db()$start[1], 1L)
    expect_null(edit_backup())
    expect_false(edit_active())
    session$setInputs(ann_details = row)
    expect_equal(unname(trim_span()), c(1, 500, 500))
  })
})
```

The feature table (`feature_rows()`) still reads the original tables until Task 16, so this test no longer compares it after the trim; `trim_span()` reads the review layer and is re-run by the same `edit_tick()` the feature table uses.

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("export", "batch-trim")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL.
- test-export.R: "unit_trimmed reads the kept window of the review copy" fails at `win(51L, 350L)`: `c(FALSE, FALSE)` instead of `c(TRUE, TRUE)` (today's `unit_trimmed()` reads `assembly_backup`).
- test-batch-trim.R: restore reasons are `"not trimmed"`; `"running" %in% names(formals(...))` is TRUE; the batch trim writes `assemblies` (`orig()` changes, `review_units` is empty); `unit_trimmed()` gives `c(FALSE, FALSE)` after `review_trim()`; Restore skips the untrimmed copy ("not trimmed"); the summary has no "their checks show not run" text and still says "re-curated"; `edit_backup()` is not NULL after the batch trim; `MP_TRIM_TIP()` still says "running run".

- [ ] **Step 3: Move batch Trim / Restore and unit_trimmed to the review layer**

In `R/export.R`:

Find:

```r
#' Is each end of a unit trimmed now? (compared with its pre-trim snapshot)
#' @return c(left, right); any() says whether the unit is trimmed at all.
#' @noRd
unit_trimmed <- function(con, id, path, scaffold) {
  b <- tryCatch(DBI::dbGetQuery(con,
    "SELECT b.sequence bak, a.sequence cur FROM assembly_backup b JOIN assemblies a
       ON a.ID = b.ID AND a.path = b.path AND a.scaffold = b.scaffold
     WHERE b.ID = ? AND b.path = ? AND b.scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold))), error = function(e) data.frame())
  if (nrow(b) == 0 || is.na(b$bak[1])) return(c(FALSE, FALSE))
  if (is.na(b$cur[1])) return(c(TRUE, TRUE))
  off <- as.integer(regexpr(b$cur[1], b$bak[1], fixed = TRUE))
  if (off < 1) return(c(TRUE, TRUE))
  c(off > 1, off + nchar(b$cur[1]) - 1 < nchar(b$bak[1]))
}
```

Replace with:

```r
#' Is each end of a unit trimmed? Read from the kept window of its review copy
#'
#' A unit without a copy, or a database without the review layer, is untrimmed.
#' Readers of the original tables pass c(FALSE, FALSE) instead of calling this.
#' @return c(left, right); any() says whether the unit is trimmed at all.
#' @noRd
unit_trimmed <- function(con, id, path, scaffold) {
  w <- tryCatch(DBI::dbGetQuery(con,
    "SELECT trim_start, trim_end, orig_length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold))), error = function(e) data.frame())
  if (nrow(w) == 0) return(c(FALSE, FALSE))
  c(isTRUE(w$trim_start[1] > 1), isTRUE(w$trim_end[1] < w$orig_length[1]))
}
```

The `tryCatch` stays: many export tests build in-memory databases without the review layer.

In `R/pipeline_steps.R`:

Find:

```r
#' Which units a batch Trim (or Restore) can act on
#'
#' @param con Open DBI connection.
#' @param u data.frame(ID, path, scaffold).
#' @param flank bp kept past the gene ends (0 = gene ends).
#' @param restore TRUE plans a Restore instead: only trimmed units qualify.
#' @param running TRUE when a run holds the project: queued units are skipped.
#' @return `u` with from, to and reason (NA where the unit qualifies).
#' @noRd
batch_trim_plan <- function(con, u, flank = 0L, restore = FALSE, running = FALSE) {
  u$from <- rep(NA_integer_, nrow(u)); u$to <- u$from; u$reason <- rep(NA_character_, nrow(u))
  key <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  for (i in seq_len(nrow(u))) {
    q <- function(sql) DBI::dbGetQuery(con, paste(sql, key), params = list(u$ID[i], u$path[i], u$scaffold[i]))
    if (isTRUE(q("SELECT annotate_lock FROM annotate")$annotate_lock[1] == 1)) {
      u$reason[i] <- "annotation locked"; next
    }
    st <- DBI::dbGetQuery(con, "SELECT assemble_switch FROM assemble WHERE ID = ?",
      params = list(u$ID[i]))
    if (running && (isTRUE(st$assemble_switch[1] %in% c(1, 4)) ||
        isTRUE(q("SELECT annotate_switch FROM annotate")$annotate_switch[1] %in% c(1, 4)))) {
      u$reason[i] <- "in the active run"; next
    }
    if (restore) {
      if (!any(unit_trimmed(con, u$ID[i], u$path[i], u$scaffold[i]))) u$reason[i] <- "not trimmed"
      next
    }
    n <- q("SELECT LENGTH(sequence) n FROM assemblies")$n
    if (!length(n) || is.na(n[1]) || n[1] == 0) { u$reason[i] <- "no sequence"; next }
    a <- q("SELECT element, start, stop FROM annotations")
```

Replace with:

```r
#' Which units a batch Trim (or Restore) can act on
#'
#' Plans from the review layer: a contig with a review copy is trimmed from the
#' copy, in copy coordinates.
#' @param con Open DBI connection.
#' @param u data.frame(ID, path, scaffold).
#' @param flank bp kept past the gene ends (0 = gene ends).
#' @param restore TRUE plans a Restore instead: every unit with a review copy qualifies.
#' @return `u` with from, to and reason (NA where the unit qualifies).
#' @noRd
batch_trim_plan <- function(con, u, flank = 0L, restore = FALSE) {
  u$from <- rep(NA_integer_, nrow(u)); u$to <- u$from; u$reason <- rep(NA_character_, nrow(u))
  key <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  for (i in seq_len(nrow(u))) {
    q <- function(sql) DBI::dbGetQuery(con, paste(sql, key), params = list(u$ID[i], u$path[i], u$scaffold[i]))
    if (isTRUE(q("SELECT annotate_lock FROM annotate")$annotate_lock[1] == 1)) {
      u$reason[i] <- "annotation locked"; next
    }
    if (restore) {
      if (!review_has_copy(con, u$ID[i], u$path[i], u$scaffold[i])) u$reason[i] <- "no edits to restore"
      next
    }
    n <- q("SELECT LENGTH(sequence) n FROM review_units_v")$n
    if (!length(n) || is.na(n[1]) || n[1] == 0) { u$reason[i] <- "no sequence"; next }
    a <- q("SELECT element, start, stop FROM review_annotations_v")
```

Find:

```r
#' Trim (or restore) each planned unit and re-curate it; one failure does not
#' stop the rest.
#' @param con Open DBI connection. @param p Eligible rows of [batch_trim_plan()].
#' @param restore TRUE restores instead of trimming.
#' @return Error message per unit, NA where it worked.
#' @noRd
batch_trim_apply <- function(con, p, restore = FALSE) {
  vapply(seq_len(nrow(p)), function(i) tryCatch({
    if (restore) restore_assembly(con, p$ID[i], p$path[i], p$scaffold[i])
    else trim_assembly(con, p$ID[i], p$path[i], p$scaffold[i], p$from[i], p$to[i])
    tryCatch(curate_one(con, p$ID[i], p$path[i], p$scaffold[i]), error = function(e) NULL)
    NA_character_
  }, error = function(e) conditionMessage(e)), character(1))
}
```

Replace with:

```r
#' Trim (or restore) each planned unit's review copy; one failure does not stop
#' the rest. A trim makes the copy on first touch; a restore deletes it.
#' @param con Open DBI connection. @param p Eligible rows of [batch_trim_plan()].
#' @param restore TRUE restores instead of trimming.
#' @return Error message per unit, NA where it worked.
#' @noRd
batch_trim_apply <- function(con, p, restore = FALSE) {
  vapply(seq_len(nrow(p)), function(i) tryCatch({
    if (restore) review_drop(con, p$ID[i], p$path[i], p$scaffold[i])
    else review_trim(con, p$ID[i], p$path[i], p$scaffold[i], p$from[i], p$to[i])
    NA_character_
  }, error = function(e) conditionMessage(e)), character(1))
}
```

In `R/app_server.R` (batch block):

Find:

```r
  # ---- batch Trim / Restore on ticked units ----
  batch_trim <- reactiveVal(NULL)
```

Replace with:

```r
  # ---- batch Trim / Restore on ticked units ----
  # A trim or restore rewrites the unit's review copy, so it ends any
  # boundary-edit session left open in the details window.
  end_edit_session <- function() {
    edit_active(FALSE); edit_backup(NULL); edit_dirty(FALSE)
    shinyjs::hide("feat_edit_ctrls")
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
  }
  batch_trim <- reactiveVal(NULL)
```

Find:

```r
    batch_trim_plan(con, b$u, batch_flank(), b$restore,
      running = !is.null(run_lock_holder(app_project())))
```

Replace with:

```r
    batch_trim_plan(con, b$u, batch_flank(), b$restore)
```

Find:

```r
    err <- batch_trim_apply(con, ok, b$restore)
    removeModal(); batch_trim(NULL); edit_tick(edit_tick() + 1); refresh(refresh() + 1)
```

Replace with:

```r
    err <- batch_trim_apply(con, ok, b$restore)
    end_edit_session()
    removeModal(); batch_trim(NULL); edit_tick(edit_tick() + 1); refresh(refresh() + 1)
```

`edit_active`, `edit_backup`, and `edit_dirty` are defined further down in `app_server()`; the helper only runs inside observers, after the whole server function has run.

The trim window's span (`trim_span()`, which the batch test reads after a trim) plans from the review layer too:

Find:

```r
    sp <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    n <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
```

Replace with:

```r
    sp <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    n <- app_query("SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
```

In `R/app_modals.R`:

Find:

```r
    title = if (restore) "Restore trimmed contigs" else "Trim contigs", easyClose = TRUE,
```

Replace with:

```r
    title = if (restore) "Restore contigs" else "Trim contigs", easyClose = TRUE,
```

Find:

```r
    tags$p(paste0(mp_n(ok, "contig"), if (ok == 1) " is" else " are",
      if (restore) " restored to full length and re-curated." else " trimmed and re-curated.")),
```

Replace with:

```r
    tags$p(paste0(mp_n(ok, "contig"), if (ok == 1) " is" else " are",
      if (restore) " restored to the original: trims, boundary edits, and 5'/3' choices are deleted."
      else " trimmed; their checks show not run until you press Update.")),
```

In `R/constants.R`:

Find:

```r
      ") are trimmed; locked contigs and contigs queued in a running run are skipped."),
    restore = paste0("Undo the trim on the selected contigs. Untrimmed and locked contigs are skipped."))
```

Replace with:

```r
      ") are trimmed; locked contigs are skipped. Trimmed contigs show their checks as not run ",
      "until you press Update."),
    restore = paste0("Restore the selected contigs to their originals: trims, boundary edits, and ",
      "5'/3' choices are deleted. Contigs with no edits and locked contigs are skipped."))
```

Check: `Rscript -e 'for (f in c("R/export.R", "R/pipeline_steps.R", "R/app_server.R", "R/app_modals.R", "R/constants.R")) invisible(parse(f))'`

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("export", "batch-trim")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS, both files `FAIL 0`.

Do not run test-trim-coverage.R yet: `trim_assembly()` (deleted in Step 9) used `unit_trimmed()` to decide when to snapshot, so its old tests fail from here until Step 9 replaces them.

- [ ] **Step 5: Write the failing trim window tests**

In `tests/testthat/test-app-export.R`:

Find:

```r
# trim_assembly() writes these; fixture_review_schema() adds any the fixture lacks
```

Replace with:

```r
# the review layer reads these; fixture_review_schema() adds any the fixture lacks
```

The locked unit must not get a review copy from the refused edits (S5.1):

Find:

```r
    "SELECT asm.length, an.problematic, an.annotate_notes,
       (SELECT stop FROM annotations WHERE element = '18S') AS stop18
     FROM assemblies asm JOIN annotate an ON an.ID = asm.ID")
```

Replace with:

```r
    "SELECT asm.length, an.problematic, an.annotate_notes,
       (SELECT stop FROM annotations WHERE element = '18S') AS stop18,
       (SELECT COUNT(*) FROM review_units) AS n_copies
     FROM assemblies asm JOIN annotate an ON an.ID = asm.ID")
```

Find:

```r
test_that("a second trim offers the unit's current span, not the first trim's", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  add_trim_cols(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(trim_mode = "annotated")
    session$flushReact()
    expect_equal(unname(trim_span()), c(10, 420, 600))
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    session$flushReact()
    expect_equal(unname(trim_span()), c(1, 411, 411))

    # a 5 bp flank around the trimmed unit keeps every base of 18S
    session$setInputs(trim_mode = "flank", trim_flank = 5)
    session$flushReact()
    expect_equal(trim_preset(), c(1, 411))
    session$setInputs(trim_from = 1, trim_to = 411, trim_apply = 2)
    session$flushReact()
    expect_equal(RiboPilot:::app_query("SELECT start, stop FROM annotations WHERE element = '18S'"),
      data.frame(start = 1L, stop = 191L))
  })
})
```

Replace with:

```r
test_that("a trim writes the review copy, and a second trim offers the copy's span", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  add_trim_cols(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0")
  DBI::dbExecute(con, "UPDATE coverage SET depths = ?", params = list(paste(1:600, collapse = ",")))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1', 1, 1, 'severe', '300')")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  alerts <- character(0)
  local_mocked_bindings(rp_alert = function(msg, ...) alerts <<- c(alerts, msg))

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(trim_mode = "annotated")
    session$flushReact()
    expect_equal(unname(trim_span()), c(10, 420, 600))
    # the junction hint is a check result, shown only once the contig's checks are done
    expect_error(output$trim_chimera_hint)
    c2 <- RiboPilot:::app_con(); fixture_mark_done(c2); DBI::dbDisconnect(c2)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_equal(output$trim_chimera_hint, "suggested junctions: 300")

    active_feature("18S")
    session$setInputs(feat_edit = TRUE)
    expect_false(is.null(edit_backup()))
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    session$flushReact()
    expect_equal(utils::tail(alerts, 1), "Contig trimmed. Its checks show not run until you press Update.")
    # the trim ends the Edit session, so a later Reset cannot write pre-trim coordinates
    expect_null(edit_backup())
    expect_false(edit_active())
    expect_equal(unname(trim_span()), c(1, 411, 411))
    expect_equal(RiboPilot:::app_query("SELECT trim_start, trim_end, check_switch FROM review_units"),
      data.frame(trim_start = 10L, trim_end = 420L, check_switch = 1L))

    # not run on the copy: no hint; once done, the copy's own junctions
    expect_error(output$trim_chimera_hint)
    c2 <- RiboPilot:::app_con()
    DBI::dbExecute(c2, "UPDATE review_units SET check_switch = 2")
    DBI::dbExecute(c2, "INSERT INTO review_chimera (ID, path, scaffold, label, junctions) VALUES ('S1', 1, 1, 'severe', '291')")
    DBI::dbDisconnect(c2)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_equal(output$trim_chimera_hint, "suggested junctions: 291")

    # the trim window plots the copy: the coverage cut to the window, genes in copy coordinates
    g <- RiboPilot:::plot_trim("S1", 1L, 1L, 1, 411)
    expect_equal(range(g$data$pos), c(1L, 411L))
    expect_equal(g$data$depth[1], 10)
    f <- g$layers[[1]]$data
    expect_equal(c(f$start[f$element == "18S"], f$stop[f$element == "18S"]), c(1L, 191L))

    # a 5 bp flank around the trimmed copy keeps every base of 18S
    session$setInputs(trim_mode = "flank", trim_flank = 5)
    session$flushReact()
    expect_equal(trim_preset(), c(1, 411))
    session$setInputs(trim_from = 1, trim_to = 411, trim_apply = 2)
    session$flushReact()
    expect_equal(RiboPilot:::app_query("SELECT start, stop FROM review_annotations_v WHERE element = '18S'"),
      data.frame(start = 1L, stop = 191L))
    # the original stays as the pipeline left it, for the Assemble step
    expect_equal(RiboPilot:::app_query("SELECT start, stop FROM annotations WHERE element = '18S'"),
      data.frame(start = 10L, stop = 200L))
    expect_equal(RiboPilot:::app_query("SELECT length FROM assemblies")$length, 600L)
  })
})

test_that("Restore deletes the copy and keeps the marks; a contig without a copy has nothing to restore", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, annotate_notes = 'keep me', export_group = 'g1'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  alerts <- character(0)
  local_mocked_bindings(rp_alert = function(msg, ...) alerts <<- c(alerts, msg))
  asked <- character(0)
  local_mocked_bindings(mp_confirm = function(id, title, text, ...) asked <<- c(asked, text))
  copies <- function() RiboPilot:::app_query("SELECT COUNT(*) AS n FROM review_units")$n
  marks <- function() RiboPilot:::app_query(
    "SELECT reviewed, problematic, partial, annotate_notes, export_group, annotate_lock FROM annotate")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(restore_btn = 1)
    expect_equal(utils::tail(alerts, 1), "This contig has no edits to restore.")
    expect_length(asked, 0)

    # marks write annotate only and never make a copy
    session$setInputs(mark_partial = 1)
    expect_equal(copies(), 0)
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    expect_equal(copies(), 1)
    before <- marks()
    expect_equal(before$partial, "yes")

    active_feature("18S")
    session$setInputs(feat_edit = TRUE)
    session$setInputs(restore_btn = 2)
    expect_match(asked[1], "411 bp now and 600 bp originally", fixed = TRUE)
    expect_match(asked[1], "Marks, notes, group, and lock stay.", fixed = TRUE)
    session$setInputs(restore_confirm = TRUE)
    expect_equal(copies(), 0)
    expect_equal(marks(), before)
    expect_null(edit_backup())
    expect_equal(utils::tail(alerts, 1), "Contig restored to its original (600 bp).")
    expect_equal(unname(trim_span()), c(10, 420, 600))
  })
})
```

- [ ] **Step 6: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-export.R")'`
Expected: FAIL in the two new tests: the trim calls `trim_assembly()`, so `review_units` stays empty (`data.frame` with 0 rows), `trim_span()` after the trim still gives `c(1, 411, 411)` only because the original was rewritten, and the original 18S check gets `start = 1` instead of 10; the alert still reads "Contig trimmed; feature coordinates + coverage shifted."; `edit_backup()` is not NULL; `output$trim_chimera_hint` shows the original junctions before the contig is done; `restore_btn` on a contig without a copy asks for confirmation instead of alerting. Every other test in the file passes.

- [ ] **Step 7: Move the trim window and the single-contig Trim and Restore to the review layer**

In `R/app_server.R` (trim and restore handlers; `end_edit_session()` is the helper added in Step 3):

Find:

```r
    s <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(s) == 0 || is.na(s$length[1])) {
      rp_alert("No contig sequence to trim.", type = "warning"); return()
    }
    span <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
```

Replace with:

```r
    s <- app_query("SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(s) == 0 || is.na(s$length[1])) {
      rp_alert("No contig sequence to trim.", type = "warning"); return()
    }
    span <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
```

The junction hint is a check result: it reads the copy's chimera rows, in copy coordinates, and only when the contig's checks are done. `edit_tick()` redraws it after a trim.

Find:

```r
  output$trim_chimera_hint <- renderText({
    u <- sel_unit(); req(u)
    ch <- app_query("SELECT junctions FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
```

Replace with:

```r
  output$trim_chimera_hint <- renderText({
    edit_tick(); u <- sel_unit(); req(u)
    ch <- app_query("SELECT c.junctions FROM review_chimera_v c JOIN review_units_v v
      ON v.ID = c.ID AND v.path = c.path AND v.scaffold = c.scaffold
      WHERE c.ID = ? AND c.path = ? AND c.scaffold = ? AND v.check_switch = 2", params = uargs())
```

Find:

```r
    b <- input$trim_brush; u <- sel_unit(); req(b, u)
    n <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
```

Replace with:

```r
    b <- input$trim_brush; u <- sel_unit(); req(b, u)
    n <- app_query("SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
```

Find:

```r
    n <- DBI::dbGetQuery(con, "SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length[1]
    if (is.na(from) || is.na(to) || from < 1 || from > to || (!is.na(n) && to > n)) {
      rp_alert(sprintf("Invalid trim range. From must be at least 1 and no more than To, and To no more than the contig's current length (%s bp).",
        n %|NA|% "?"), type = "error"); return()
    }
    ok <- tryCatch({
      trim_assembly(con, u$ID, u$path, u$scaffold, from, to)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # re-curate after the shift
      TRUE
    }, error = function(e) { rp_alert(conditionMessage(e), type = "error"); FALSE })
    if (isTRUE(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()  # back to the details modal, now showing trimmed coords
      rp_alert("Contig trimmed; feature coordinates + coverage shifted.", type = "message")
    }
  })
```

Replace with:

```r
    n <- DBI::dbGetQuery(con, "SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length[1]
    if (is.na(from) || is.na(to) || from < 1 || from > to || (!is.na(n) && to > n)) {
      rp_alert(sprintf("Invalid trim range. From must be at least 1 and no more than To, and To no more than the contig's current length (%s bp).",
        n %|NA|% "?"), type = "error"); return()
    }
    if (!app_write(review_trim(con, u$ID, u$path, u$scaffold, from, to))) return()
    end_edit_session()
    edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    show_annotate_unit()  # back to the details modal, now showing trimmed coords
    rp_alert("Contig trimmed. Its checks show not run until you press Update.", type = "message")
  })
```

Find:

```r
  # ---- restore the active unit to its original (pre-trim) length ----
  observeEvent(input$restore_btn, {
    u <- sel_unit(); req(u)
    if (edits_block() || lock_blocked()) return()
    len <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length[1]
    mp_confirm("restore_confirm", title = "Restore the full assembly?",
      text = sprintf("The contig is %s bp now. Restoring undoes the trim and re-runs curation.", len %||% "?"),
      action_label = "Restore", danger = TRUE)
  })
  observeEvent(input$restore_confirm, {
    req(isTRUE(input$restore_confirm))
    u <- sel_unit(); req(u)
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ok <- tryCatch({
      len <- restore_assembly(con, u$ID, u$path, u$scaffold)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # recompute flags for the full unit
      len
    }, error = function(e) { rp_alert(conditionMessage(e), type = "warning"); NULL })
    if (!is.null(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()
      rp_alert(paste0("Contig restored to original length (", ok, " bp)."), type = "message")
    }
  })
```

Replace with:

```r
  # ---- restore the active unit: delete its review copy ----
  observeEvent(input$restore_btn, {
    u <- sel_unit(); req(u)
    if (edits_block() || lock_blocked()) return()
    v <- app_query("SELECT has_copy, length, orig_length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    if (!isTRUE(v$has_copy[1] == 1)) {
      rp_alert("This contig has no edits to restore.", type = "warning"); return()
    }
    mp_confirm("restore_confirm", title = "Restore the original contig?",
      text = sprintf(paste("The contig is %s bp now and %s bp originally. Restore deletes its trims,",
        "boundary edits, 5'/3' choices, and the checks run on them. Marks, notes, group, and lock stay."),
        v$length[1] %|NA|% "?", v$orig_length[1] %|NA|% "?"),
      action_label = "Restore", danger = TRUE)
  })
  observeEvent(input$restore_confirm, {
    req(isTRUE(input$restore_confirm))
    u <- sel_unit(); req(u)
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(review_drop(con, u$ID, u$path, u$scaffold))) return()
    end_edit_session()
    edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    show_annotate_unit()
    len <- DBI::dbGetQuery(con, "SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length[1]
    rp_alert(paste0("Contig restored to its original (", len %|NA|% "?", " bp)."), type = "message")
  })
```

In `R/app_modals.R`:

Find:

```r
    opts_help("Trimming cuts the flanks off the contig and shifts every ",
      "feature coordinate + the coverage track to match. Features that fall entirely ",
      "outside the kept range are removed. Use \"Restore assembly\" to undo."),
```

Replace with:

```r
    opts_help("Trimming cuts the flanks off this contig's review copy and shifts its ",
      "feature coordinates and coverage track to match. Features that fall entirely ",
      "outside the kept range are removed. The original stays as it is for the Assemble step. ",
      "The contig's checks show not run until you press Update. Use \"Restore assembly\" to ",
      "undo every edit."),
```

Find:

```r
#' Read-coverage track for one assembly unit, with feature regions shaded
#' @noRd
plot_coverage <- function(id, path = 1L, scaffold = 1L) {
  if (is.null(id)) {
    return(NULL)
  }
  cov <- app_query("SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(cov) == 0 || is.na(cov$depths[1]) || !nzchar(cov$depths[1])) {
```

Replace with:

```r
#' Read-coverage track for one unit on the review layer, with feature regions shaded
#'
#' A review copy shows the original track cut to its kept window, and its genes,
#' in copy coordinates.
#' @noRd
plot_coverage <- function(id, path = 1L, scaffold = 1L) {
  if (is.null(id)) {
    return(NULL)
  }
  con <- app_con(); on.exit(DBI::dbDisconnect(con))
  dep <- review_coverage(con, id, path, scaffold)$depths
  if (!length(dep)) {
```

Find:

```r
  dep <- parse_depths(cov$depths)
  df <- data.frame(pos = seq_along(dep), depth = dep)
  feats <- app_query("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
```

Replace with:

```r
  df <- data.frame(pos = seq_along(dep), depth = dep)
  feats <- app_query("SELECT element, start, stop FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
```

`plot_trim()` calls `plot_coverage()` and needs no change; the trim window is the only caller of both.

Check: `Rscript -e 'for (f in c("R/app_server.R", "R/app_modals.R")) invisible(parse(f))'`

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("app-export", "batch-trim")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS, both files `FAIL 0`.

- [ ] **Step 9: Delete the in-place trim functions and port their helper tests**

Nothing calls `trim_assembly()`, `restore_assembly()`, `ensure_assembly_backup()`, or `chimera_snapshot()` any more. `write_coverage_summary()`, `chimera_shift()`, and `chimera_restore()` stay: the Task 22 upgrade writes a trimmed contig's original back with them.

In `R/pipeline_steps.R`:

Find:

```r
#' Trim one assembly unit (ID,path,scaffold) to \[from, to\], shifting its
#' annotations and coverage to match. Features that fall entirely outside the
#' kept range are dropped; overlapping features are clamped. Used by the Annotate
#' tab's "Trim assembly" control.
#'
#' @param con Open DBI connection.
#' @param id,path,scaffold Assembly unit.
#' @param from,to 1-based inclusive kept range in the unit's coordinates.
#' @return (invisibly) the new unit length.
#' @noRd
trim_assembly <- function(con, id, path = 1L, scaffold = 1L, from, to) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  s <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(s) == 0 || is.na(s$sequence[1]) || !nzchar(s$sequence[1])) {
    stop("No sequence for contig ", id, ".", path, ".", scaffold, call. = FALSE)
  }
  n <- nchar(s$sequence[1])
  from <- max(1L, as.integer(from)); to <- min(n, as.integer(to))
  if (from > to) stop("Invalid trim range.", call. = FALSE)

  # Snapshot the untrimmed unit state, so [restore_assembly()] can undo this.
  # A trim of an already trimmed unit keeps the snapshot, so restore always
  # returns the full unit.
  ensure_assembly_backup(con)
  if (!any(unit_trimmed(con, id, path, scaffold))) {
    ann0 <- DBI::dbGetQuery(con, "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    cov0 <- DBI::dbGetQuery(con, "SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    DBI::dbExecute(con,
      "INSERT OR REPLACE INTO assembly_backup (ID, path, scaffold, sequence, length, annotations, coverage, chimera, time_stamp)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(id, path, scaffold, s$sequence[1], n,
        as.character(jsonlite::toJSON(ann0, dataframe = "rows", na = "null")),
        if (nrow(cov0)) cov0$depths[1] else NA_character_,
        chimera_snapshot(con, id, path, scaffold), as.integer(Sys.time())))
  }
  new_seq <- substr(s$sequence[1], from, to); new_len <- nchar(new_seq)

  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(new_seq, new_len, id, path, scaffold))

  ann <- DBI::dbGetQuery(con, "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  for (i in seq_len(nrow(ann))) {
    ns <- ann$start[i] - from + 1; nt <- ann$stop[i] - from + 1
    if (nt < 1 || ns > new_len) {
      DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = list(id, path, scaffold, ann$element[i]))
      next
    }
    # a cut through the gene makes that end partial
    minus <- identical(ann$strand[i], "-")
    cut <- c(if (ns < 1) (if (minus) "partial_stop" else "partial_start"),
      if (nt > new_len) (if (minus) "partial_start" else "partial_stop"))
    for (col in cut) {
      DBI::dbExecute(con, sprintf(
        "UPDATE annotations SET %s = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?", col),
        params = list(id, path, scaffold, ann$element[i]))
    }
    ns <- max(1L, ns); nt <- min(new_len, nt)
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = list(ns, nt, nt - ns + 1, id, path, scaffold, ann$element[i]))
  }

  cov <- DBI::dbGetQuery(con, "SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  dep <- if (nrow(cov) > 0) parse_depths(cov$depths) else numeric(0)
  # `from` is clamped to the sequence, which can outrun a shorter depth track
  # (truncated `samtools depth -a`, or a re-assembly coverage never re-ran on)
  if (length(dep) > 0 && from <= length(dep)) {
    write_coverage_summary(con, id, path, scaffold, dep[seq(from, min(length(dep), to))])
  }
  DBI::dbExecute(con, "UPDATE coverage SET clips = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  chimera_shift(con, id, path, scaffold, from - 1L, new_len)
  invisible(new_len)
}

#' Rewrite one unit's coverage row from a depth track
```

Replace with:

```r
#' Rewrite one unit's coverage row from a depth track
```

Find:

```r
#' Ensure the assembly_backup table exists (idempotent).
#' @noRd
ensure_assembly_backup <- function(con) {
  DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)
}

#' Reconcile a unit's annotate row after its assemblies.ignore flag changed
```

Replace with:

```r
#' Reconcile a unit's annotate row after its assemblies.ignore flag changed
```

Find:

```r
#' Restore a sample's assembly to its original (pre-trim) state
#'
#' Reverses one or more `trim_assembly()` calls by writing back the snapshot
#' taken at the first trim: the full-length sequence, the original annotation
#' rows, and the untrimmed coverage track. The snapshot is left in place so a
#' fresh trim + restore returns to the same original.
#' @param con Open project DB connection. @param id,path,scaffold Assembly unit.
#' @return (invisibly) the restored length.
#' @noRd
restore_assembly <- function(con, id, path = 1L, scaffold = 1L) {
  ensure_assembly_backup(con)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  bak <- DBI::dbGetQuery(con,
    "SELECT sequence, length, annotations, coverage, chimera FROM assembly_backup WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(bak) == 0 || is.na(bak$sequence[1])) {
    stop("No original assembly on record for this contig (nothing to restore).", call. = FALSE)
  }
  cur <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$sequence[1]
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(bak$sequence[1], bak$length[1], id, path, scaffold))

  # one transaction: a failed restore after a committed delete would leave the
  # unit with neither the trimmed annotations nor the snapshot ones
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    if (!is.na(bak$annotations[1]) && nzchar(bak$annotations[1])) {
      ann <- jsonlite::fromJSON(bak$annotations[1])
      if (is.data.frame(ann) && nrow(ann) > 0) {
        ann$ID <- id; ann$path <- path; ann$scaffold <- scaffold
        DBI::dbAppendTable(con, "annotations", ann)
      }
    }
  })

  write_coverage_summary(con, id, path, scaffold, parse_depths(bak$coverage))
  if (!is.na(bak$chimera[1])) {
    chimera_restore(con, id, path, scaffold, bak$chimera[1])
  } else {
    off <- if (is.na(cur)) -1L else as.integer(regexpr(cur, bak$sequence[1], fixed = TRUE))
    if (off > 1L) chimera_shift(con, id, path, scaffold, 1L - off, bak$length[1])
  }
  invisible(bak$length[1])
}

#' Kept window for a trim preset: gene ends widened by `flank`, clipped to the unit
```

Replace with:

```r
#' Kept window for a trim preset: gene ends widened by `flank`, clipped to the unit
```

In `R/chimera.R`:

Find:

```r
#' A unit's chimera result and windows as JSON, kept with the trim backup
#' @noRd
chimera_snapshot <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  tryCatch({
    ch <- DBI::dbGetQuery(con, "SELECT * FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    if (!nrow(ch)) return(NA_character_)
    w <- DBI::dbGetQuery(con, "SELECT * FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    as.character(jsonlite::toJSON(list(chimera = ch, windows = w), dataframe = "rows", na = "null", digits = NA))
  }, error = function(e) NA_character_)
}

#' Put back the chimera result and windows saved by `chimera_snapshot()`
```

Replace with:

```r
#' Put back a chimera result and windows stored as JSON (assembly_backup.chimera)
```

In `tests/testthat/test-chimera.R`, delete the test that trims with `trim_assembly()` (Task 17 adds the review-layer chimera plot test):

Find:

```r
test_that("plot_chimera keeps the original check after a trim and marks the cut", {
  con <- chimera_test_con()
  db <- DBI::dbGetInfo(con)$dbname
  sq <- withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 4500, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,?,4500)", params = list(sq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand) VALUES ('S1',1,1,'18S',1001,2800,'+')")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, junctions) VALUES ('S1',1,1,0.9,'severe','3500')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident) VALUES ('S1',1,1,1,500,99), ('S1',1,1,3001,3500,97)")
  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, 501, 3000)
  DBI::dbExecute(con, "UPDATE annotations SET stop = stop - 100")
  DBI::dbDisconnect(con)
  g <- RiboPilot:::plot_chimera("S1", 1L, 1L, db = db)
  b <- ggplot2::ggplot_build(g)
  expect_equal(b$layout$panel_params[[1]]$x.range, c(0, 4500))
  geoms <- vapply(g$layers, function(l) class(l$geom)[1], "")
  gene <- b$data[[which(geoms == "GeomRect")[1]]]
  expect_equal(c(gene$xmin, gene$xmax), c(1001, 2700))
  expect_equal(nrow(b$data[[which(geoms == "GeomSegment")[1]]]), 2L)
  expect_true(3500 %in% unlist(lapply(b$data[geoms == "GeomVline"], `[[`, "xintercept")))
  cut <- unlist(lapply(b$data[geoms == "GeomVline"], `[[`, "xintercept"))
  expect_true(all(c(500, 3000) %in% cut))
})

test_that("plot_chimera with NULL junctions/segments has no junction line", {
```

Replace with:

```r
test_that("plot_chimera with NULL junctions/segments has no junction line", {
```

Overwrite `tests/testthat/test-trim-coverage.R` with the following. The trim and restore tests become direct tests of the helpers Task 22 keeps; the short-track cases now belong to `review_coverage()` (Task 3) and the cut-end partials to `review_trim()` (Task 4), which test them.

```r
# Helpers the schema step 4 upgrade uses to write a trimmed contig's original
# back: the coverage summary (median_depth is the column the contig table and
# triage read) and the chimera rows, moved by a shift or put back from the JSON
# an older trim stored in assembly_backup.chimera.

mk_trim_db <- function(dir, seq, depths) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT,
                       length INT, sequence TEXT, edited INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
                       element TEXT, start INT, stop INT, length INT, edited INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT,
                       mean_depth REAL, median_depth REAL, min_depth INT, max_depth INT,
                       depths TEXT, clips TEXT, time_stamp INT)")
  DBI::dbExecute(con, RiboPilot:::CHIMERA_DDL)
  DBI::dbExecute(con, RiboPilot:::CHIMERA_WINDOWS_DDL)
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?,?,0)",
    params = list(nchar(seq), seq))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, mean_depth, median_depth, min_depth,
                       max_depth, depths, time_stamp) VALUES('S1',1,1,?,?,?,?,?,0)",
    params = list(mean(depths), stats::median(depths), min(depths), max(depths),
      paste(depths, collapse = ",")))
  con
}

cov_row <- function(con) {
  DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, min_depth, max_depth, depths
                          FROM coverage WHERE ID='S1' AND path=1 AND scaffold=1")
}

# the unit's chimera rows as JSON, in the form an older trim stored them
chimera_json <- function(con) {
  ch <- DBI::dbGetQuery(con, "SELECT * FROM chimera")
  w <- DBI::dbGetQuery(con, "SELECT * FROM chimera_windows")
  as.character(jsonlite::toJSON(list(chimera = ch, windows = w), dataframe = "rows", na = "null", digits = NA))
}

test_that("write_coverage_summary rewrites median_depth with the track and skips an empty one", {
  dir <- withr::local_tempdir()
  dep <- c(rep(2, 50), rep(20, 50))
  con <- mk_trim_db(dir, strrep("ACGT", 25), dep)
  on.exit(DBI::dbDisconnect(con))

  expect_true(RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep[51:100]))
  cv <- cov_row(con)
  expect_equal(cv$median_depth, 20)
  expect_equal(cv$mean_depth, 20)
  expect_equal(cv$min_depth, 20L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 50L)

  RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep)
  cv <- cov_row(con)
  expect_equal(cv$median_depth, stats::median(dep))
  expect_equal(cv$mean_depth, mean(dep))
  expect_equal(cv$min_depth, 2L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 100L)

  expect_false(RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, numeric(0)))
  expect_equal(cov_row(con), cv)
})

test_that("chimera_shift moves junctions and windows both ways, and chimera_restore puts the JSON back", {
  dir <- withr::local_tempdir()
  seq <- paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")
  con <- mk_trim_db(dir, seq, rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1',1,1,'severe','30,70')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  ch <- function() DBI::dbGetQuery(con, "SELECT junctions, segments FROM chimera")
  win <- function() DBI::dbGetQuery(con, "SELECT start, stop FROM chimera_windows ORDER BY start")
  snap <- chimera_json(con)

  # keep 41..100: what was cut goes, the rest moves 40 bp left
  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 40L, 60L)
  expect_equal(ch()$junctions, "30")
  expect_equal(win()$start, 11L)
  expect_equal(jsonlite::fromJSON(ch()$segments)$start, 11L)

  # the reverse shift puts the kept rows back in full-length coordinates
  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, -40L, 100L)
  expect_equal(ch()$junctions, "70")
  expect_equal(win()$start, 51L)

  RiboPilot:::chimera_restore(con, "S1", 1L, 1L, snap)
  expect_equal(ch()$junctions, "30,70")
  expect_equal(win()$start, c(1L, 51L))
})

test_that("shifting away every junction clears the chimera label, and chimera_restore brings it back", {
  dir <- withr::local_tempdir()
  seq <- paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")
  con <- mk_trim_db(dir, seq, rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, junctions, reason)
    VALUES ('S1',1,1,0.9,'severe','30','1-40 C1; 41-90 C2; junction near 30')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  ch <- function() DBI::dbGetQuery(con, "SELECT label, score, junctions, reason FROM chimera")
  snap <- chimera_json(con)

  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 44L, 56L)
  expect_equal(ch()$label, "clean")
  expect_equal(ch()$score, 0)
  expect_match(ch()$reason, "junctions trimmed away", fixed = TRUE)
  expect_false(grepl("41-90", ch()$reason, fixed = TRUE))

  RiboPilot:::chimera_restore(con, "S1", 1L, 1L, snap)
  expect_equal(ch()$label, "severe")
  expect_equal(ch()$junctions, "30")
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM chimera_windows")), 2L)
})
```

Check that nothing refers to the deleted functions:
Run: `grep -rn "trim_assembly\|restore_assembly\|ensure_assembly_backup\|chimera_snapshot" R/ inst/ tests/testthat/`
Expected: no output.

Check: `Rscript -e 'for (f in c("R/pipeline_steps.R", "R/chimera.R", "tests/testthat/test-trim-coverage.R", "tests/testthat/test-chimera.R")) invisible(parse(f))'`

- [ ] **Step 10: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("batch-trim", "export", "app-export", "trim-coverage", "chimera")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS, every file `FAIL 0`.

Then the files that reach `unit_trimmed()`, the batch or trim code, or the app fixtures:
Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-edits", "review-copy", "review-window", "seqview", "export-opts", "export-headers", "export-units", "nmnh", "app-table", "app-opts", "notes-partial", "app-scaffold-uid", "load-assemblies", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0`.

Then the whole suite once: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: `FAIL 0`.

- [ ] **Step 11: Commit**

```bash
git add R/export.R R/pipeline_steps.R R/chimera.R R/app_server.R R/app_modals.R R/constants.R tests/testthat/test-export.R tests/testthat/test-batch-trim.R tests/testthat/test-app-export.R tests/testthat/test-chimera.R tests/testthat/test-trim-coverage.R
git commit -m "Trim and Restore write review copies"
```


### Task 16: Edit mode and State on copies

**Files:**
- Create: `tests/testthat/test-review-handlers.R`
- Modify: `R/app_server.R` (State handler `observeEvent(input$update_state, ...)` :845-853; details window: `edit_backup` declaration :2446, `sel_unit()` observer :2456, `unit_partial_rows()`/`partial_state()` :2464-2477, `nudge_ctx` :2651-2655, `feature_rows` :2715-2716, `restore_edit_backup()`/`save_edits()` :2801-2821, the Edit backup query :2831, the notes reload flag :3000-3043, `nudge_boundary()` :3063-3066, `toggle_partial()` :3075-3080)
- Modify: `R/app_ui_helpers.R` (`mp_state_modal()` :193-202)
- Modify: `tests/testthat/test-app-export.R` (the partial toggles test, :713-762)
- Modify: `tests/testthat/test-ui-helpers.R` (one test appended)
- Test: `tests/testthat/test-review-handlers.R`, `tests/testthat/test-app-export.R`, `tests/testthat/test-ui-helpers.R`

**Interfaces:**
- Consumes (Task 2, `R/review.R`): `REVIEW_EFF_PARTIALS`, `review_copy(con, id, path, scaffold)` (tests only). (Task 4, `R/review.R`): `review_set_bounds(con, id, path, scaffold, element, start, stop)`, `review_set_partial(con, id, path, scaffold, element, end = c("start", "stop"), value)`, `review_reset_session(con, id, path, scaffold, backup)` with `backup = data.frame(element, start, stop, partial_start_user, partial_stop_user)`, `review_set_state(con, units, state)` returning NA or a reason (`"unlock first"`, `"press Update on the Assemble step first"`) per unit. (Task 15, `R/export.R`): `unit_trimmed(con, id, path, scaffold)` reading `review_units_v`. (Task 1): `review_units_v`, `review_annotations_v`, test helpers `mk_app_project()`, `fixture_mark_done(con)`.
- Produces: in `app_server()`, `reload_on_close` (reactiveVal; replaces `notes_dirty`: set it TRUE and the table reloads when the details window closes) and `set_check_state(u, st)` (State on Annotation Review). `mp_state_modal(ids, current, mode = "annotate")` offers only codes `c("1", "0")`.

Annotation Review's Edit mode now works on the contig's review copy (S5.1, S5.2): a boundary nudge calls `review_set_bounds()` (copy made on first touch, checks invalidated), a 5'/3' toggle calls `review_set_partial()` (checks unchanged), and Reset or Discard call `review_reset_session()` with the session's starting bounds and user partials, which invalidates again only when a coordinate changed. Save ends the session and no longer runs `curate_one()` (D6); flags show "not run" until Update. The details window's feature readers (first feature, feature table, partial state, nudge context) read the sequence and length from `review_units_v` and the features from `review_annotations_v` with `REVIEW_EFF_PARTIALS` in an explicit column list, never `SELECT *` (a `SELECT *, COALESCE(...) AS partial_start` returns the raw column under `$partial_start`). Because an edit changes the contig's State, nudges, toggles, and resets set the existing notes flag (renamed `reload_on_close`) so the table reloads when the window closes, not on every click. State on Annotation Review offers only Ready to run and On hold (D7) and goes through `review_set_state()`; locked and waiting contigs are left unchanged and listed with the reason. The marks, notes, and group still write `annotate` only and never make a copy (S5.1).

- [ ] **Step 1: Write the failing Edit mode tests**

Create `tests/testthat/test-review-handlers.R`:

```r
# Edit mode and State on Annotation Review write the contig's review copy
# (SPEC S5.1, S5.2, D6, D7); the original tables never change.

# mk_app_project() with every contig 600 bp; p1.s1 holds 18S, ITS1, and 5.8S on
# the + strand. Every contig is done (annotated, ITSx stored) and unlocked.
mk_handler_project <- function(env = parent.frame()) {
  dir <- mk_app_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 600", params = list(strrep("ACGT", 150)))
  for (f in list(c("18S", 10, 200), c("ITS1", 201, 260), c("5.8S", 261, 420))) {
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, strand, status, flags)
      VALUES ('S1', 1, 1, ?, ?, ?, ?, '+', 'ok', '')",
      params = list(f[1], RiboPilot:::rrna_element_type(f[1]), as.integer(f[2]), as.integer(f[3])))
  }
  fixture_mark_done(con)
  dir
}

hq <- function(sql, params = NULL) RiboPilot:::app_query(sql, params)
hx <- function(sql) {
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, sql)
}

test_that("a nudge edits the review copy and queues its checks; marks and Save make no copy and run no curation", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  curated <- 0L
  local_mocked_bindings(curate_one = function(...) curated <<- curated + 1L)
  copies <- function() hq("SELECT COUNT(*) AS n FROM review_units")$n

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(mark_reviewed = 1, mark_partial = 1, annotate_notes = "check the 18S start")
    session$setInputs(tbl__reactable__selected = 1:4, reviewed_btn = 1, partial_btn = 1)
    session$setInputs(modal_hidden = 1)
    expect_equal(copies(), 0)

    expect_equal(active_feature(), "18S")
    session$setInputs(feat_edit = TRUE)
    r0 <- refresh()
    session$setInputs(nt_5_out = 1)
    expect_equal(hq("SELECT start FROM annotations WHERE element = '18S'")$start, 10L)
    expect_equal(hq("SELECT start, stop, length, edited, status FROM review_annotations WHERE element = '18S'"),
      data.frame(start = 9L, stop = 200L, length = 192L, edited = 1L, status = NA_character_))
    expect_equal(hq("SELECT path, scaffold, edited, check_switch FROM review_units"),
      data.frame(path = 1L, scaffold = 1L, edited = 1L, check_switch = 1L))
    expect_equal(feature_rows()$start, c(9L, 201L, 261L))
    expect_equal(nudge_ctx()$start, 9L)

    # the table reloads once the window closes, not on every nudge
    expect_equal(refresh(), r0)
    session$setInputs(modal_hidden = 2)
    expect_equal(refresh(), r0 + 1)

    session$setInputs(feat_save = 1)
    expect_equal(curated, 0L)
    expect_equal(hq("SELECT check_switch FROM review_units")$check_switch, 1L)
  })
})

test_that("Reset and Discard put back the session's bounds and 5'/3' choices", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  # curation's verdict, which Reset must not turn into a user choice
  hx("UPDATE annotations SET partial_stop = 1 WHERE element = '18S'")
  feat <- function() hq("SELECT start, stop, partial_start_user, partial_stop_user FROM review_annotations WHERE element = '18S'")
  unit <- function() hq("SELECT check_switch, CAST(rev AS TEXT) AS rev FROM review_units")
  raised <- function(rev) hq("SELECT rev - CAST(? AS INTEGER) AS d FROM review_units", list(rev))$d

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(feat_edit = TRUE)
    session$setInputs(toggle_partial_start = 1)
    expect_equal(feat(), data.frame(start = 10L, stop = 200L, partial_start_user = 1L, partial_stop_user = NA_integer_))
    u0 <- unit()
    # a 5'/3' choice leaves the checks as they were
    expect_equal(u0$check_switch, 2L)
    session$setInputs(feat_reset = 1)
    expect_equal(feat(), data.frame(start = 10L, stop = 200L, partial_start_user = NA_integer_,
      partial_stop_user = NA_integer_))
    expect_equal(unit(), u0)

    session$setInputs(nt_5_out = 1)
    expect_equal(feat()$start, 9L)
    expect_equal(raised(u0$rev), 1)
    hx("UPDATE review_units SET check_switch = 2")
    session$setInputs(feat_reset = 2)
    expect_equal(feat()$start, 10L)
    expect_equal(unit()$check_switch, 1L)
    expect_equal(raised(u0$rev), 2)

    # leaving Edit with unsaved nudges and choosing Discard does the same
    session$setInputs(feat_edit = FALSE)
    session$setInputs(feat_edit = TRUE)
    session$setInputs(nt_3_in = 1)
    expect_equal(feat()$stop, 199L)
    session$setInputs(feat_edit = FALSE)
    session$setInputs(unsaved_confirm = FALSE)
    expect_equal(feat()$stop, 200L)
    expect_equal(raised(u0$rev), 4)
  })
  expect_equal(hq("SELECT start, stop, partial_start, partial_stop FROM annotations WHERE element = '18S'"),
    data.frame(start = 10L, stop = 200L, partial_start = NA_integer_, partial_stop = 1L))
})

test_that("the details window reads the copy: features, effective 5'/3' ends, sequence, and length", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  # the copy keeps 181..600 of the original: 18S is gone and 5.8S runs to the copy's end
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?, length = 420, trim_start = 181, trim_end = 600",
    params = list(strrep("T", 420)))
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE element = '18S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET start = start - 180, stop = stop - 180")
  DBI::dbExecute(con, "UPDATE review_annotations SET stop = 420 WHERE element = '5.8S'")
  # curation said 5' partial; the user said no
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start = 1, partial_start_user = 0 WHERE element = 'ITS1'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    expect_equal(active_feature(), "ITS1")
    f <- feature_rows()
    expect_equal(f$element, c("ITS1", "5.8S"))
    expect_equal(f$start, c(21L, 81L))
    expect_equal(f$partial_start, c(0L, 0L))
    expect_equal(f$partial_stop, c(0L, 1L))
    expect_equal(output$feat_edit_active, "Active feature: ITS1")
    active_feature("5.8S")
    session$flushReact()
    expect_equal(output$feat_edit_active, "Active feature: 5.8S (partial 3')")
    cx <- nudge_ctx()
    expect_equal(c(cx$start, cx$stop), c(81L, 420L))
    expect_equal(cx$seq, strrep("T", 420))
  })
})
```

In `tests/testthat/test-app-export.R`, the partial toggles test now reads the user columns of the review layer (`as.integer()`: a contig without a copy reads the view's `NULL AS partial_start_user`, which RSQLite returns as logical NA), and its scaffold-end case edits the copy the first toggle made (the original features stay as they were):

Find:

```r
  p <- function() unname(unlist(RiboPilot:::app_query(
    "SELECT partial_start, partial_stop FROM annotations WHERE element = '18S'")))
```

Replace with:

```r
  p <- function() as.integer(unlist(RiboPilot:::app_query(
    "SELECT partial_start_user, partial_stop_user FROM review_annotations_v WHERE element = '18S'")))
```

Find:

```r
    c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE annotations SET start = 1, partial_start = NULL WHERE element = '18S'")
```

Replace with:

```r
    c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE review_annotations SET start = 1, partial_start_user = NULL WHERE element = '18S'")
```

Find:

```r
    session$setInputs(toggle_partial_stop = 2)
    expect_equal(p(), c(0L, 1L))
  })
})
```

Replace with:

```r
    session$setInputs(toggle_partial_stop = 2)
    expect_equal(p(), c(0L, 1L))
  })
  # the choices live on the review copy; the original feature never changes
  expect_equal(RiboPilot:::app_query("SELECT start, partial_start, partial_stop FROM annotations WHERE element = '18S'"),
    data.frame(start = 10L, partial_start = NA_integer_, partial_stop = NA_integer_))
})
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-handlers", "app-export")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-handlers.R: in the first test the nudge writes the original (`annotations` 18S start is `9`, not `10`), `review_units` and `review_annotations` stay empty (0-row data frames), the closing window does not reload the table, and `curated` is `1` after Save; in the second test `feat()` has 0 rows and `u0$check_switch` is empty; in the third test `feature_rows()$element` is `c("18S", "ITS1", "5.8S")` (the originals) and `active_feature()` is `"18S"`. test-app-export.R, 6 failures, all in the partial toggles test: `p()` is `NA NA` after each toggle (the toggles wrote `annotations`), `output$feat_edit_active` lacks `5'`, and the new last expectation finds the original 18S `partial_start`/`partial_stop` set to `1`. Every other test in both files passes.

- [ ] **Step 3: Point Edit mode at the review copy (R/app_server.R)**

The edit session's backup now holds the user partial columns, which is what `review_reset_session()` writes back.

Find:

```r
  edit_backup <- reactiveVal(NULL)    # annotations snapshot for Reset
```

Replace with:

```r
  edit_backup <- reactiveVal(NULL)    # bounds and 5'/3' choices at Edit, for Reset
```

Find:

```r
      els <- app_query("SELECT element FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
```

Replace with:

```r
      els <- app_query("SELECT element FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
```

Find:

```r
  # effective 5'/3' partial state of the unit's features (curator's flag, else automatic)
  unit_partial_rows <- function(d) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    n <- DBI::dbGetQuery(con, "SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length
```

Replace with:

```r
  # effective 5'/3' partial state of the unit's features (user's choice, else
  # curation's, else automatic); d comes from review_annotations_v with REVIEW_EFF_PARTIALS
  unit_partial_rows <- function(d) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    n <- DBI::dbGetQuery(con, "SELECT length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length
```

Find:

```r
    unit_partial_rows(app_query("SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(uargs(), list(el))))
```

Replace with:

```r
    unit_partial_rows(app_query(paste0("SELECT element, start, stop, strand, ", REVIEW_EFF_PARTIALS,
      " FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?"),
      params = c(uargs(), list(el))))
```

Find:

```r
    a <- app_query("SELECT start, stop, length, type FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(uargs(), list(el)))
    if (nrow(a) == 0) return(NULL)
    seq <- app_query("SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
```

Replace with:

```r
    a <- app_query("SELECT start, stop, length, type FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(uargs(), list(el)))
    if (nrow(a) == 0) return(NULL)
    seq <- app_query("SELECT sequence FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
```

Find:

```r
    unit_partial_rows(app_query("SELECT element, type, start, stop, length, strand, partial_start, partial_stop, status, flags
       FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", params = uargs()))
```

Replace with:

```r
    unit_partial_rows(app_query(paste0("SELECT element, type, start, stop, length, strand, ", REVIEW_EFF_PARTIALS,
      ", status, flags FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start"),
      params = uargs()))
```

Find:

```r
  restore_edit_backup <- function() {
    b <- edit_backup()
    if (is.null(b) || !nrow(b)) return(TRUE)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    app_write(DBI::dbWithTransaction(con, for (i in seq_len(nrow(b))) {
      DBI::dbExecute(con,
        "UPDATE annotations SET start = ?, stop = ?, length = ?, partial_start = ?, partial_stop = ?
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, b$partial_start[i], b$partial_stop[i],
          u$ID, u$path, u$scaffold, b$element[i]))
    }))
  }
  save_edits <- function() {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # refresh flags after edits
    edit_backup(NULL); edit_dirty(FALSE)
    edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    rp_toast("Boundary edits saved.", type = "message")
  }
```

Replace with:

```r
  restore_edit_backup <- function() {
    b <- edit_backup()
    if (is.null(b) || !nrow(b)) return(TRUE)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    ok <- app_write(review_reset_session(con, u$ID, u$path, u$scaffold, b))
    if (ok) reload_on_close(TRUE)
    ok
  }
  save_edits <- function() {
    edit_backup(NULL); edit_dirty(FALSE)
    edit_tick(edit_tick() + 1)
    rp_toast("Boundary edits saved.", type = "message")
  }
```

Find:

```r
      edit_backup(app_query("SELECT element, start, stop, partial_start, partial_stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
```

Replace with:

```r
      edit_backup(app_query("SELECT element, start, stop, partial_start_user, partial_stop_user FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ?",
```

The notes flag becomes the window's reload flag (five places):

Find:

```r
  # the table reloads once the window closes, not on every keystroke
  notes_dirty <- reactiveVal(FALSE)
  observeEvent(input$modal_hidden, {
    if (!notes_dirty()) return()
    notes_dirty(FALSE); refresh(refresh() + 1)
  })
```

Replace with:

```r
  # notes and edits reach the table once the window closes, not on every keystroke or nudge
  reload_on_close <- reactiveVal(FALSE)
  observeEvent(input$modal_hidden, {
    if (!reload_on_close()) return()
    reload_on_close(FALSE); refresh(refresh() + 1)
  })
```

Find:

```r
      params = c(list(new), uargs())))) notes_dirty(TRUE)
```

Replace with:

```r
      params = c(list(new), uargs())))) reload_on_close(TRUE)
```

Find:

```r
    notes_dirty(TRUE)
    updateTextAreaInput(session, "annotate_notes", value = "")
```

Replace with:

```r
    reload_on_close(TRUE)
    updateTextAreaInput(session, "annotate_notes", value = "")
```

Find:

```r
    if (app_write(DBI::dbExecute(con, "UPDATE assemble SET user_notes = ? WHERE ID = ?", params = list(new, id)))) notes_dirty(TRUE)
```

Replace with:

```r
    if (app_write(DBI::dbExecute(con, "UPDATE assemble SET user_notes = ? WHERE ID = ?", params = list(new, id)))) reload_on_close(TRUE)
```

Find:

```r
    notes_dirty(TRUE)
    updateTextAreaInput(session, "assemble_user_notes", value = "")
```

Replace with:

```r
    reload_on_close(TRUE)
    updateTextAreaInput(session, "assemble_user_notes", value = "")
```

The nudge and the toggles write the copy:

Find:

```r
    if (!app_write(DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = list(start, stop, stop - start + 1, u$ID, u$path, u$scaffold, cx$el)))) return()
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
```

Replace with:

```r
    if (!app_write(review_set_bounds(con, u$ID, u$path, u$scaffold, cx$el, start, stop))) return()
    edit_dirty(TRUE); reload_on_close(TRUE); edit_tick(edit_tick() + 1)
```

Find:

```r
    to <- if (isTRUE(partial_state(el)[[col]][1] == 1)) 0L else 1L
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con, sprintf(
      "UPDATE annotations SET %s = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      col), params = c(list(to), uargs(), list(el))))) return()
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
```

Replace with:

```r
    to <- if (isTRUE(partial_state(el)[[col]][1] == 1)) 0L else 1L
    u <- sel_unit()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(review_set_partial(con, u$ID, u$path, u$scaffold, el, sub("partial_", "", col), to))) return()
    edit_dirty(TRUE); reload_on_close(TRUE); edit_tick(edit_tick() + 1)
```

Check: `Rscript -e 'invisible(parse("R/app_server.R"))'` parses, and `grep -n "notes_dirty\|UPDATE annotations" R/app_server.R` prints nothing. Other `FROM annotations` and `FROM assemblies` reads in the file (trim window, curation panel, chimera panel, Assemble step) belong to Tasks 15 and 17; leave them.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-handlers", "app-export", "notes-partial", "batch-trim")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-handlers.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 33 ]`; test-app-export.R, test-notes-partial.R, and test-batch-trim.R `FAIL 0` (counts as Tasks 14 and 15 left them; test-notes-partial.R covers the renamed notes flag, test-batch-trim.R the Edit session a trim ends).

- [ ] **Step 5: Write the failing State tests**

Append to `tests/testthat/test-review-handlers.R`:

```r

test_that("State on Annotation Review queues or holds the recheck on copies and refuses locked and waiting contigs", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE path = 2 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE path = 2 AND scaffold = 2")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  said <- character(0)
  local_mocked_bindings(mp_alert = function(title, text = NULL, ...) said <<- c(said, paste(title, text)))
  copies <- function() hq("SELECT path, scaffold, check_switch FROM review_units ORDER BY path, scaffold")
  rev_of <- function(p, s) hq("SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE path = ? AND scaffold = ?",
    list(p, s))$rev

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    row <- function(p, s) { d <- samples_data(); which(d$path == p & d$scaffold == s) }
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(2, 1), row(2, 2)))
    session$setInputs(annot_state_btn = 1)
    session$setInputs(new_state = "1", update_state = 1)
    expect_equal(copies(), data.frame(path = 1L, scaffold = 1L, check_switch = 1L))
    expect_equal(said, paste("Left 2 contigs unchanged",
      "S1 p2.s1 (unlock first), S1 p2.s2 (press Update on the Assemble step first)."))

    # On hold keeps rev; Ready to run again raises it
    r1 <- rev_of(1, 1)
    session$setInputs(state_row = as.character(row(1, 1)))
    session$setInputs(new_state = "0", update_state = 2)
    expect_equal(copies()$check_switch, 0L)
    expect_equal(rev_of(1, 1), r1)
    session$setInputs(state_row = as.character(row(1, 1)))
    session$setInputs(new_state = "1", update_state = 3)
    expect_equal(copies()$check_switch, 1L)
    expect_equal(hq("SELECT rev - CAST(? AS INTEGER) AS d FROM review_units", list(r1))$d, 1)

    # On hold on a contig without a copy makes one, held
    session$setInputs(state_row = as.character(row(1, 2)))
    session$setInputs(new_state = "0", update_state = 4)
    expect_equal(copies(), data.frame(path = c(1L, 1L), scaffold = c(1L, 2L), check_switch = c(1L, 0L)))
  })
  # the main pipeline's queue is not touched
  expect_equal(hq("SELECT annotate_switch FROM annotate ORDER BY path, scaffold")$annotate_switch, c(2L, 2L, 2L, 1L))
})
```

Append to `tests/testthat/test-ui-helpers.R`:

```r

test_that("the State window on Annotation Review offers only Ready to run and On hold", {
  html <- as.character(RiboPilot:::mp_state_modal(c("S1 p1.s1", "S1 p1.s2"), mode = "annotate"))
  expect_match(html, "Set state for 2 contigs")
  expect_match(html, "Ready to run.*On hold")
  expect_false(grepl("In progress|Success|Failed", html))
  expect_match(html, "Checked again by the next Update on Annotation Review", fixed = TRUE)
})
```

- [ ] **Step 6: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-handlers", "ui-helpers")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-handlers.R `[ FAIL 7 | WARN 0 | SKIP 0 | PASS 34 ]`, all in the new State test: State still writes `annotate.annotate_switch` for every selected contig (the last expectation reads `1 0 1 1`, not `2 2 2 1`), so `copies()` has 0 rows, `said` is empty (no "Left 2 contigs unchanged" alert), and the rev checks get a 0-length result. test-ui-helpers.R, 3 failures in the new test: `"Ready to run.*On hold"` does not match (MP_STATE_SETTABLE lists On hold first), the window still lists In progress, Success, and Failed, and the Annotation Review wording is missing.

- [ ] **Step 7: State goes through review_set_state()**

In `R/app_ui_helpers.R`:

Find:

```r
mp_state_modal <- function(ids, current = character(0), mode = "assemble") {
  codes <- MP_STATE_SETTABLE
  effect <- c(
    `0` = "Skipped by the next Update.",
    `1` = "Processed by the next Update; results this module already stored are replaced.",
    `4` = "Picked up by the next run; normally set by the pipeline itself.",
    `2` = "Treated as finished; the next Update skips it.",
    `3` = "Treated as finished with an error; the next Update skips it."
  )
```

Replace with:

```r
mp_state_modal <- function(ids, current = character(0), mode = "assemble") {
  annot <- identical(mode, "annotate")
  codes <- if (annot) c("1", "0") else MP_STATE_SETTABLE
  effect <- c(
    `0` = "Skipped by the next Update.",
    `1` = "Processed by the next Update; results this module already stored are replaced.",
    `4` = "Picked up by the next run; normally set by the pipeline itself.",
    `2` = "Treated as finished; the next Update skips it.",
    `3` = "Treated as finished with an error; the next Update skips it."
  )
  if (annot) effect[c("1", "0")] <- c(
    "Checked again by the next Update on Annotation Review; its check results are cleared now.",
    "Skipped by the Update on Annotation Review until set to Ready to run.")
```

In `R/app_server.R`:

Find:

```r
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state)
    set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st), u = state_units())
    removeModal()
  })
```

Replace with:

```r
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state)
    if (identical(mm$pfx, "annotate")) set_check_state(state_units(), st)
    else set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st), u = state_units())
    removeModal()
  })
  # State on Annotation Review queues (1) or holds (0) the recheck of each
  # contig's review copy; locked and waiting contigs are listed, not changed.
  set_check_state <- function(u, st) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ok <- app_write(why <- review_set_state(con, u, st))
    refresh(refresh() + 1)
    if (!ok) return()
    done <- is.na(why)
    if (any(done)) mp_toast(paste(mp_n(sum(done), "contig"), "updated."))
    if (!all(done)) mp_alert(paste("Left", mp_n(sum(!done), "contig"), "unchanged"),
      text = paste0(mp_id_list(sprintf("%s p%s.s%s (%s)", u$ID[!done], u$path[!done],
        u$scaffold[!done], why[!done])), "."), type = "warning")
  }
```

`refresh()` runs even after a failed write because `review_set_state()` commits unit by unit, so units before the failing one have changed. Check: `Rscript -e 'invisible(parse("R/app_server.R")); invisible(parse("R/app_ui_helpers.R"))'`.

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-handlers", "ui-helpers", "app-table", "app-export", "notes-partial", "batch-trim", "app-scaffold-uid", "seqview", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-handlers.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 41 ]`; test-ui-helpers.R gains 4 passing expectations over what Task 14 left; every other file `FAIL 0`. test-app-table.R's row-icon State test runs on the Assemble step and still writes `assemble.assemble_switch`; test-docs-site.R checks the new modal text says "contig" and "step" only.

- [ ] **Step 9: Commit**

```bash
git add R/app_server.R R/app_ui_helpers.R tests/testthat/test-review-handlers.R tests/testthat/test-app-export.R tests/testthat/test-ui-helpers.R
git commit -m "Edit mode and State write review copies"
```


### Task 17: Details window readers and the Assemble chimera panel

**Files:**
- Modify: `R/seqview.R` (`seqview_send()`, :71-90)
- Modify: `R/app_modals.R` (delete `chimera_original()`, :1016-1044; `plot_chimera()`, :1046-1123)
- Modify: `R/app_server.R` (Annotation Review viewer observer, :2499-2502; chimera panel outputs from `hit_unit` to `output$chimera_windows`, :2506-2603; curation panel `flag_rows` and `output$curate_notes`, :2766-2780; line numbers as of HEAD 36564f2, Tasks 13-16 move them)
- Test: `tests/testthat/test-chimera.R`, `tests/testthat/test-seqview.R`

**Interfaces:**
- Consumes:
  - Task 1: the views `review_units_v`, `review_annotations_v`, `review_chimera_v`, `review_chimera_windows_v`; test helper `fixture_mark_done(con)`.
  - Task 2: `REVIEW_EFF_PARTIALS` (`"COALESCE(partial_start_user, partial_start) AS partial_start, COALESCE(partial_stop_user, partial_stop) AS partial_stop"`); `review_copy(con, id, path, scaffold)` and `review_invalidate(con, id, path, scaffold)` (tests only).
  - Task 3: `review_coverage(con, id, path, scaffold)` -> `list(depths, clips)`, cut to the copy's window.
  - Task 4: `review_trim(con, id, path, scaffold, from, to)` (tests only).
  - Task 14: `review_unit_state(con, units = NULL)` -> data.frame with a `state` column (`"done"`, `"not run"`, `"running"`, `"failed"`, `"on hold"`, `"waiting"`); `review_state_text(state)`.
  - Task 15: `unit_trimmed(con, id, path, scaffold)` reads the window from `review_units_v`; the old trim test in test-chimera.R is gone.
- Produces:
  - `seqview_send(session, canvas, input, id, path, scaffold, selected = NULL, review = FALSE)`: `review = TRUE` reads the sequence from `review_units_v`, features from `review_annotations_v` with `REVIEW_EFF_PARTIALS`, tracks from `review_coverage()`, and passes `unit_trimmed()` to `partial_ends()`; `FALSE` reads the originals with `trimmed = c(FALSE, FALSE)`. The Annotation Review viewer passes `TRUE`; the Assemble viewer keeps the default.
  - `plot_chimera(id, path = 1L, scaffold = 1L, db = NULL, review = FALSE)`: `review = TRUE` reads `review_chimera_v`, `review_chimera_windows_v`, `review_annotations_v`, and `review_units_v.length`, and returns an empty plot labelled `"Not run yet. Press Update."` unless `review_unit_state()` says `"done"`; `FALSE` reads the originals and draws no trim overlay. `chimera_original()` is deleted.
  - Inside `app_server()`: reactive `sel_check_state()` (the `state` of the contig open in Annotation Review; depends on `edit_tick()` and `refresh()`), reactive `chim_unit()`, and helper `chim_query(cols, tbl, tail = "")`, which reads `review_<tbl>_v` on Annotation Review (only once the contig is done) and `<tbl>` on the Assemble step.

The details window of Annotation Review now shows the review copy, and the Assemble window keeps showing the originals (S8, details window and Assemble step). The sequence viewer reads the copy's sequence, its features with effective partials, and the original coverage cut to the copy's window; the Assemble viewer passes no trimmed ends, because originals are never trimmed (S4.4, S4.5). The chimera panel shared by both windows picks its layer from `input$mode`: on Annotation Review it reads the review views and shows "Not run yet. Press Update." until the contig's checks are done, which also hides rows a discarded recheck left behind (S5.4); on the Assemble step it shows the original check, and the pre-trim juggling in `chimera_original()` goes. The curation warnings panel reads `review_units_v` and `review_annotations_v` and is gated on the same state, so a nudge or trim that clears the copy's checks no longer shows "No warnings" (D6).

- [ ] **Step 1: Write the failing tests for the two readers**

Task 15 deleted the old test "plot_chimera keeps the original check after a trim and marks the cut" (it called `trim_assembly()`); this task replaces it. Check: `grep -c "keeps the original check after a trim" tests/testthat/test-chimera.R` prints `0`.

Append to the end of `tests/testthat/test-chimera.R`:

```r

test_that("plot_chimera shows the copy on Annotation Review once its checks are done, and the original on the Assemble step", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  db <- DBI::dbGetInfo(con)$dbname
  sq <- withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 4500, TRUE), collapse = ""))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,?,4500)", params = list(sq))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, n_features = 1, itsx_warnings = NULL")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand) VALUES ('S1',1,1,'18S',1001,2800,'+')")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, junctions) VALUES ('S1',1,1,0.9,'severe','3500')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident) VALUES ('S1',1,1,1,500,99), ('S1',1,1,3001,3500,97)")
  draw <- function(review) {
    g <- RiboPilot:::plot_chimera("S1", 1L, 1L, db = db, review = review)
    list(b = ggplot2::ggplot_build(g), geoms = vapply(g$layers, function(l) class(l$geom)[1], ""))
  }
  layer <- function(d, geom) d$b$data[d$geoms == geom]
  note <- function(d) unlist(lapply(layer(d, "GeomText"), `[[`, "label"))
  shows <- function(d, xmax, gene, junc, nwin) {
    expect_equal(d$b$layout$panel_params[[1]]$x.range, c(0, xmax))
    r <- layer(d, "GeomRect")[[1]]
    expect_equal(c(r$xmin, r$xmax), gene)
    expect_equal(unique(unlist(lapply(layer(d, "GeomVline"), `[[`, "xintercept"))), junc)
    expect_equal(nrow(layer(d, "GeomSegment")[[1]]), nwin)
  }
  waiting <- "Not run yet. Press Update."

  # no copy: the original counts as done once its ITSx value is stored
  expect_equal(note(draw(TRUE)), waiting)
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = ''")
  shows(draw(TRUE), 4500, c(1001, 2800), 3500, 2L)

  # a trimmed copy waits for its recheck, even when a discarded run left rows behind
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 501, 3000)
  expect_equal(note(draw(TRUE)), waiting)
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, score, label, junctions) VALUES ('S1',1,1,0.4,'minor','1200')")
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start, stop, pident) VALUES ('S1',1,1,1,500,98)")
  expect_equal(note(draw(TRUE)), waiting)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2")
  shows(draw(TRUE), 2500, c(501, 2300), 1200, 1L)

  # the Assemble step keeps the original check, with no trim marks
  shows(draw(FALSE), 4500, c(1001, 2800), 3500, 2L)
  expect_false("trim" %in% note(draw(FALSE)))
  expect_false(exists("chimera_original", envir = asNamespace("RiboPilot"), inherits = FALSE))
})
```

Append to the end of `tests/testthat/test-seqview.R`:

```r

test_that("the Annotation Review viewer shows the copy and the Assemble viewer the original", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  sq <- strrep("ACGTACGTAC", 10)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 100 WHERE path = 1 AND scaffold = 1",
    params = list(sq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,40,'+'), ('S1',1,1,'28S',61,90,'+')")
  DBI::dbExecute(con, "UPDATE coverage SET depths = ?, clips = ? WHERE path = 1 AND scaffold = 1",
    params = list(paste(1:100 * 10, collapse = ","), paste(c(rep(0, 20), 105, rep(0, 79)), collapse = ",")))
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 90)
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  feat <- function(el) Filter(function(f) f$row == el, sent$features)[[1]]

  RiboPilot:::seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", "S1", 1L, 1L, review = TRUE)
  expect_equal(sent$len, 70L)
  expect_equal(sent$seq, substr(sq, 21, 90))
  expect_equal(c(feat("18S")$pos1, feat("18S")$pos2), c(1L, 20L))
  expect_true(feat("18S")$partial5)    # cut by the trim
  expect_equal(c(feat("28S")$pos1, feat("28S")$pos2), c(41L, 70L))
  expect_false(feat("28S")$partial3)   # trimmed exactly at its end: not an assembly end
  expect_equal(sent$depth, as.integer(21:90 * 10))
  expect_equal(sent$clip[1], 0.5)
  expect_equal(sent$input, "seqv_ann_pick")

  RiboPilot:::seqview_send(session, "seqv_asm-canvas", NULL, "S1", 1L, 1L)
  expect_equal(sent$len, 100L)
  expect_equal(sent$seq, sq)
  expect_equal(c(feat("18S")$pos1, feat("18S")$pos2), c(1L, 40L))
  expect_true(feat("18S")$partial5)    # starts at the contig's first base
  expect_false(feat("28S")$partial3)
  expect_equal(sent$depth, as.integer(1:100 * 10))
  expect_equal(sent$clip[21], 0.5)
})
```

What the expectations pin: `review_trim(21, 90)` cuts 18S (1..40) on the left, so the copy's 18S is 1..20 with `partial_start_user = 1`, and leaves 28S (61..90 -> 41..70) ending exactly at the copy's right end with no user partial. Reading the raw `partial_start` (the `SELECT *` trap) loses the 18S mark; passing `c(FALSE, FALSE)` on the copy marks 28S 3' as an assembly end; calling `unit_trimmed()` for the originals (trimmed on both ends now) drops the original 18S 5' mark; an uncut coverage track shifts depth and clips by 20.

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R"); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: FAIL. test-chimera.R: 1 failure, `Error in RiboPilot:::plot_chimera("S1", 1L, 1L, db = db, review = review): unused argument (review = review)`. test-seqview.R: 1 failure, `Error in RiboPilot:::seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", "S1", 1L, 1L, review = TRUE): unused argument (review = TRUE)`. Every other test in both files passes or skips as before.

- [ ] **Step 3: Make seqview_send read either layer (R/seqview.R)**

Find:

```r
#' @param canvas Canvas id. @param input Input id for gene clicks, or NULL.
#' @noRd
seqview_send <- function(session, canvas, input, id, path, scaffold, selected = NULL) {
  p <- list(id, path, scaffold)
  s <- app_query("SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  if (!nrow(s) || is.na(s$sequence[1])) return(invisible())
  f <- app_query("SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", params = p)
  if (nrow(f)) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    e <- partial_ends(f, nchar(s$sequence[1]), unit_trimmed(con, id, path, scaffold))
    f$partial_start <- as.integer(e$p5); f$partial_stop <- as.integer(e$p3)
  }
  cov <- app_query("SELECT * FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  dep <- if (nrow(cov)) parse_depths(cov$depths[1]) else NULL
  clp <- if (nrow(cov) && !is.null(cov$clips)) parse_depths(cov$clips[1]) else NULL
  pl <- seqview_payload(s$sequence[1], f, dep, clp, selected)
```

Replace with:

```r
#' @param canvas Canvas id. @param input Input id for gene clicks, or NULL.
#' @param review TRUE (Annotation Review): the review layer, with the copy's
#'   coverage window and trimmed ends. FALSE (Assemble step): the originals.
#' @noRd
seqview_send <- function(session, canvas, input, id, path, scaffold, selected = NULL, review = FALSE) {
  con <- app_con(); on.exit(DBI::dbDisconnect(con))
  q <- function(sql, tail = "") DBI::dbGetQuery(con,
    paste(sql, "WHERE ID = ? AND path = ? AND scaffold = ?", tail), params = list(id, path, scaffold))
  s <- q(paste("SELECT sequence FROM", if (review) "review_units_v" else "assemblies"))
  if (!nrow(s) || is.na(s$sequence[1])) return(invisible())
  f <- q(paste("SELECT element, start, stop, strand,", if (review) REVIEW_EFF_PARTIALS else "partial_start, partial_stop",
    "FROM", if (review) "review_annotations_v" else "annotations"), "ORDER BY start")
  if (nrow(f)) {
    e <- partial_ends(f, nchar(s$sequence[1]), if (review) unit_trimmed(con, id, path, scaffold) else c(FALSE, FALSE))
    f$partial_start <- as.integer(e$p5); f$partial_stop <- as.integer(e$p3)
  }
  tr <- if (review) review_coverage(con, id, path, scaffold) else {
    cov <- q("SELECT * FROM coverage")
    list(depths = if (nrow(cov)) parse_depths(cov$depths[1]),
      clips = if (nrow(cov) && !is.null(cov$clips)) parse_depths(cov$clips[1]))
  }
  pl <- seqview_payload(s$sequence[1], f, tr$depths, tr$clips, selected)
```

`coverage` stays `SELECT *` because hand-built fixtures lack `clips`. Check: `Rscript -e 'invisible(parse("R/seqview.R"))'`.

- [ ] **Step 4: Make plot_chimera read either layer; delete chimera_original (R/app_modals.R)**

Find:

```r
#' A unit's chimera check as first run: before its first trim, in untrimmed
#' coordinates, with the current genes moved to match and the kept span
#' @param q Query function (sql, params).
#' @noRd
chimera_original <- function(q, id, path, scaffold) {
  ch <- q("SELECT label, junctions, segments, channels FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  win <- q("SELECT start, stop, pident, lineage FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", list(id, path, scaffold))
  feats <- q("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  n <- q("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))$length
  # a trimmed unit shows the check from before its first trim, in untrimmed coordinates
  kept <- NULL
  bak <- tryCatch(q("SELECT b.sequence bak, b.length, b.chimera, a.sequence cur FROM assembly_backup b
    JOIN assemblies a ON a.ID = b.ID AND a.path = b.path AND a.scaffold = b.scaffold
    WHERE b.ID = ? AND b.path = ? AND b.scaffold = ?", list(id, path, scaffold)), error = function(e) data.frame())
  if (nrow(bak) && !is.na(bak$bak[1]) && !is.na(bak$cur[1])) {
    off <- as.integer(regexpr(bak$cur[1], bak$bak[1], fixed = TRUE)) - 1L
    if (off >= 0L && nchar(bak$cur[1]) < bak$length[1]) {
      kept <- c(off, off + nchar(bak$cur[1]))
      feats$start <- feats$start + off; feats$stop <- feats$stop + off
      n <- bak$length[1]
      if (!is.na(bak$chimera[1])) {
        s <- jsonlite::fromJSON(bak$chimera[1])
        ch <- s$chimera
        win <- if (is.data.frame(s$windows)) s$windows[order(s$windows$start), ] else win[0, ]
      }
    }
  }
  list(ch = ch, win = win, feats = feats, n = n, kept = kept)
}

plot_chimera <- function(id, path = 1L, scaffold = 1L, db = NULL) {
  q <- if (is.null(db)) app_query else {
    con <- open_project_db(db); on.exit(DBI::dbDisconnect(con))
    function(sql, params = NULL) DBI::dbGetQuery(con, sql, params = params)
  }
  empty <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = msg) + ggplot2::theme_void()
  o <- chimera_original(q, id, path, scaffold)
  ch <- o$ch; win <- o$win; feats <- o$feats; n <- o$n; kept <- o$kept
```

Replace with:

```r
#' Chimera check plot for one unit
#'
#' review = TRUE (Annotation Review) reads the review layer and shows a note
#' instead until the contig's checks are done. FALSE (Assemble step) reads the
#' originals.
#' @noRd
plot_chimera <- function(id, path = 1L, scaffold = 1L, db = NULL, review = FALSE) {
  con <- if (is.null(db)) app_con() else open_project_db(db)
  on.exit(DBI::dbDisconnect(con))
  q <- function(sql, params = list(id, path, scaffold)) DBI::dbGetQuery(con, sql, params = params)
  empty <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = msg) + ggplot2::theme_void()
  if (review && !identical(review_unit_state(con, data.frame(ID = id, path = as.integer(path),
      scaffold = as.integer(scaffold)))$state, "done")) {
    return(empty("Not run yet. Press Update."))
  }
  tbl <- function(x) if (review) paste0("review_", x, "_v") else x
  unit <- " WHERE ID = ? AND path = ? AND scaffold = ?"
  ch <- q(paste0("SELECT label, junctions, segments FROM ", tbl("chimera"), unit))
  win <- q(paste0("SELECT start, stop, pident FROM ", tbl("chimera_windows"), unit, " ORDER BY start"))
  feats <- q(paste0("SELECT element, start, stop FROM ", tbl("annotations"), unit))
  n <- q(paste0("SELECT length FROM ", if (review) "review_units_v" else "assemblies", unit))$length
```

Find:

```r
  if (!is.null(kept)) {
    cuts <- c(if (kept[1] > 0) kept[1], if (kept[2] < n) kept[2])
    g <- g +
      ggplot2::annotate("rect", xmin = c(0, kept[2]), xmax = c(kept[1], n), ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::geom_vline(xintercept = cuts, linetype = "dashed", colour = "#111111", linewidth = 0.7) +
      ggplot2::geom_text(data = fr(data.frame(x = cuts, y = 12), "genes"), ggplot2::aes(x = x, y = y),
        label = "trim", vjust = -0.3, size = 3.2, colour = "#111111")
  }
  g + ggplot2::facet_grid(row ~ ., scales = "free_y", space = "free_y", switch = "y",
```

Replace with:

```r
  g + ggplot2::facet_grid(row ~ ., scales = "free_y", space = "free_y", switch = "y",
```

The samples query further down (`q("SELECT lineage_taxids FROM samples WHERE ID = ?", list(id))`) keeps working with the new `q`. Check: `Rscript -e 'invisible(parse("R/app_modals.R"))'` and `grep -rn "chimera_original" R/` lists only `R/app_server.R` (Step 8 removes it).

- [ ] **Step 5: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R"); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: PASS, `FAIL 0` in both files. The older plot_chimera tests ("builds from stored windows and segments", "with NULL junctions/segments has no junction line") pass unchanged with the default `review = FALSE`.

- [ ] **Step 6: Write the failing tests for the details windows**

Append to the end of `tests/testthat/test-chimera.R`:

```r

test_that("the chimera panel reads the review layer on Annotation Review and the originals on the Assemble step", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, reason, junctions, db_source)
    VALUES ('S1',1,1,'severe','original reason','3500','/x/orig.fa')")
  fixture_mark_done(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label, reason, junctions, db_source)
    VALUES ('S1',1,1,'minor','copy reason','1200','/x/copy.fa')")
  drew <- NULL
  local_mocked_bindings(plot_chimera = function(id, path = 1L, scaffold = 1L, db = NULL, review = FALSE) {
    drew <<- review
    ggplot2::ggplot()
  })

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    # not run: the row a discarded recheck left stays hidden
    expect_error(output$chimera_result)
    expect_error(output$chimera_db_source)

    DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2")
    refresh(refresh() + 1)
    session$flushReact()
    expect_match(output$chimera_result$html, "minor", fixed = TRUE)
    expect_match(output$chimera_result$html, "copy reason", fixed = TRUE)
    expect_equal(output$chimera_db_source, "windows searched against: copy.fa")
    expect_match(output$chimera_junctions$html, '"position":[1200]', fixed = TRUE)
    output$chimera_plot
    expect_true(drew)

    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()
    expect_match(output$chimera_result$html, "severe", fixed = TRUE)
    expect_equal(output$chimera_db_source, "windows searched against: orig.fa")
    expect_match(output$chimera_junctions$html, '"position":[3500]', fixed = TRUE)
    output$chimera_plot
    expect_false(drew)
  })
})
```

Append to the end of `tests/testthat/test-seqview.R`:

```r

test_that("Annotation Review shows curation warnings only once the contig's checks are done", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand, status, flags)
    VALUES ('S1',1,1,'18S',1,40,'+','flag','original feature flag'), ('S1',1,1,'28S',61,90,'+','ok',NULL)")
  DBI::dbExecute(con, "UPDATE annotate SET n_flags = 1, curate_flags = 'original contig flag'
    WHERE path = 1 AND scaffold = 1")
  fixture_mark_done(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    expect_match(output$curate_notes$html, "1 contig warning(s): ", fixed = TRUE)
    expect_match(output$curate_notes$html, "original contig flag", fixed = TRUE)
    expect_equal(flag_rows()$flags, "original feature flag")

    # an edit makes a copy and clears its checks
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_match(output$curate_notes$html, "Checks not run. Curation warnings appear once they are done.", fixed = TRUE)
    expect_error(flag_rows())

    # the recheck finished
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, n_flags = 2, curate_flags = 'copy contig flag'")
    DBI::dbExecute(con, "UPDATE review_annotations SET status = 'flag', flags = 'copy feature flag' WHERE element = '28S'")
    refresh(refresh() + 1)
    session$flushReact()
    expect_match(output$curate_notes$html, "2 contig warning(s): ", fixed = TRUE)
    expect_match(output$curate_notes$html, "copy contig flag", fixed = TRUE)
    expect_equal(flag_rows()$element, "28S")
    expect_equal(flag_rows()$flags, "copy feature flag")
  })
})

test_that("each details window sends its own layer to the sequence viewer", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  sent <- list()
  local_mocked_bindings(seqview_send = function(session, canvas, input, id, path, scaffold, selected = NULL,
    review = FALSE) sent[[canvas]] <<- review)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()
  })
  expect_true(sent[["seqv_ann-canvas"]])
  expect_false(sent[["seqv_asm-canvas"]])
})
```

- [ ] **Step 7: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R"); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: FAIL. test-chimera.R `FAIL 6`, all in "the chimera panel reads the review layer on Annotation Review and the originals on the Assemble step": `Expected output$chimera_result to throw a error.` (and the same for `chimera_db_source`), then `minor`, `copy reason`, and `copy.fa` mismatches (the outputs still read the original `chimera` row), and `output$chimera_junctions` errors with `could not find function "chimera_original"`. test-seqview.R `FAIL 7`: in "Annotation Review shows curation warnings only once the contig's checks are done" the first three expectations pass, then the `Checks not run.` note does not match, `Expected flag_rows() to throw a error.`, and the `2 contig warning(s)`, `copy contig flag`, `28S`, and `copy feature flag` expectations fail (the panel still reads `annotate` and `annotations`); in "each details window sends its own layer to the sequence viewer" `Expected sent[["seqv_ann-canvas"]] to be TRUE.`

- [ ] **Step 8: Wire the details windows (R/app_server.R)**

The Annotation Review viewer reads the copy. Find:

```r
    try(seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", u$ID, u$path, u$scaffold, active_feature()), silent = TRUE)
```

Replace with:

```r
    try(seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", u$ID, u$path, u$scaffold, active_feature(),
      review = TRUE), silent = TRUE)
```

The Assemble viewer call (`seqview_send(session, "seqv_asm-canvas", NULL, id, h$path, h$scaffold)`) stays as it is: the default `review = FALSE` reads the originals.

The open contig's check state, and the chimera panel's layer. Find:

```r
  # Chimera panels follow the annotate modal's unit, or the assemble modal's
  # highlighted scaffold.
  hit_unit <- reactive({
```

Replace with:

```r
  # Check state of the contig open in Annotation Review; its check results show
  # only when it is done.
  sel_check_state <- reactive({
    edit_tick(); refresh(); u <- sel_unit(); req(u)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    st <- review_unit_state(con, data.frame(ID = u$ID, path = as.integer(u$path),
      scaffold = as.integer(u$scaffold)))$state
    req(length(st) == 1)
    st
  })

  # Chimera panels follow the annotate modal's unit, or the assemble modal's
  # highlighted scaffold. Annotation Review reads the review layer once the
  # contig's checks are done; the Assemble step reads the originals.
  chim_unit <- reactive({
    edit_tick(); refresh(); u <- hit_unit()
    review <- identical(input$mode, "Annotate")
    if (review) req(identical(sel_check_state(), "done"))
    c(u, review = review)
  })
  chim_query <- function(cols, tbl, tail = "") {
    u <- chim_unit()
    app_query(sprintf("SELECT %s FROM %s WHERE ID = ? AND path = ? AND scaffold = ?%s", cols,
      if (u$review) paste0("review_", tbl, "_v") else tbl, tail), list(u$ID, u$path, u$scaffold))
  }
  hit_unit <- reactive({
```

The plot (it shows its own "Not run yet" note, so it does not go through `chim_unit()`). Find:

```r
    edit_tick(); refresh(); u <- hit_unit(); grid::grid.draw(fixed_gutter(plot_chimera(u$ID, u$path, u$scaffold)))
```

Replace with:

```r
    edit_tick(); refresh(); u <- hit_unit()
    grid::grid.draw(fixed_gutter(plot_chimera(u$ID, u$path, u$scaffold, review = identical(input$mode, "Annotate"))))
```

Find:

```r
  output$chimera_result <- renderUI({
    edit_tick(); refresh(); u <- hit_unit()
    ch <- app_query("SELECT label, reason FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
```

Replace with:

```r
  output$chimera_result <- renderUI({
    ch <- chim_query("label, reason", "chimera")
```

Find:

```r
    edit_tick(); refresh(); u <- hit_unit()
    ch <- app_query("SELECT db_source FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
```

Replace with:

```r
    ch <- chim_query("db_source", "chimera")
```

Find:

```r
    edit_tick(); refresh()
    u <- hit_unit()
    ch <- app_query("SELECT segments FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$segments[1]))
```

Replace with:

```r
    u <- hit_unit()
    ch <- chim_query("segments", "chimera")
    req(nrow(ch) == 1, !is.na(ch$segments[1]))
```

Junctions are in the coordinates of the layer shown, so the "now trimmed away" note goes. Find:

```r
    edit_tick(); refresh(); u <- hit_unit()
    o <- chimera_original(app_query, u$ID, u$path, u$scaffold)
    req(NROW(o$ch) == 1)
    j <- chimera_junction_table(o$ch$junctions[1], o$ch$channels[1])
    if (!nrow(j)) return(tags$p(class = "mp-help-text", "No possible junction was found."))
    if (!is.null(o$kept)) j$why <- ifelse(j$position <= o$kept[1] | j$position > o$kept[2],
      paste(j$why, "This part is now trimmed away."), j$why)
```

Replace with:

```r
    ch <- chim_query("junctions, channels", "chimera")
    req(nrow(ch) == 1)
    j <- chimera_junction_table(ch$junctions[1], ch$channels[1])
    if (!nrow(j)) return(tags$p(class = "mp-help-text", "No possible junction was found."))
```

Find:

```r
    edit_tick(); refresh()
    u <- hit_unit()
    w <- app_query("SELECT start, stop, species, lineage, pident, bitscore, rank_distance, depth_median, clip_frac FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      list(u$ID, u$path, u$scaffold))
```

Replace with:

```r
    u <- hit_unit()
    w <- chim_query("start, stop, species, lineage, pident, bitscore, rank_distance, depth_median, clip_frac",
      "chimera_windows", " ORDER BY start")
```

Find:

```r
    ch <- app_query("SELECT segments FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    samp <- app_query("SELECT lineage_taxids FROM samples WHERE ID = ?", list(u$ID))
    cols <- if (nrow(ch) == 1 && !is.na(ch$segments[1])) {
```

Replace with:

```r
    ch <- chim_query("segments", "chimera")
    samp <- app_query("SELECT lineage_taxids FROM samples WHERE ID = ?", list(u$ID))
    cols <- if (nrow(ch) == 1 && !is.na(ch$segments[1])) {
```

The curation warnings panel reads the review layer and waits for done. Find:

```r
  flag_rows <- reactive({
    edit_tick(); u <- sel_unit(); req(u)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND status = 'flag' ORDER BY start",
      params = uargs())
  })
  output$curate_notes <- renderUI({
    u <- sel_unit(); req(u)
    a <- app_query("SELECT n_flags, curate_flags FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
```

Replace with:

```r
  flag_rows <- reactive({
    req(identical(sel_check_state(), "done"))
    app_query("SELECT element, flags FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? AND status = 'flag' ORDER BY start",
      params = uargs())
  })
  output$curate_notes <- renderUI({
    st <- sel_check_state()
    if (st != "done") {
      return(span(class = "mp-help-text",
        sprintf("Checks %s. Curation warnings appear once they are done.", review_state_text(st))))
    }
    a <- app_query("SELECT n_flags, curate_flags FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
```

When the contig is not done, `flag_tbl` and the `flag_jump` observer stop at `flag_rows()`'s `req()`, so the table stays empty under the note. Checks: `Rscript -e 'invisible(parse("R/app_server.R"))'`; `grep -rn "chimera_original" R/` prints nothing; `grep -n "FROM chimera WHERE\|FROM chimera_windows WHERE" R/app_server.R` prints nothing (Task 15 already moved `output$trim_chimera_hint` to `review_chimera_v`).

- [ ] **Step 9: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-chimera.R"); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: PASS, `FAIL 0` in both files.

Then every file that runs `app_server` or reads the touched functions (the batch-trim test mocks `plot_chimera()` and still sees a redraw after a trim, because `output$chimera_plot` keeps `edit_tick()` and `refresh()`; the Assemble examine tests read originals as before):
Run: `Rscript -e 'devtools::load_all("."); for (f in c("batch-trim", "app-scaffold-uid", "app-export", "app-table", "app-opts", "notes-partial", "blast-provenance", "specimen-app", "app-run", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0`. If `tests/testthat/test-review-handlers.R` exists (Task 16), run it too: `FAIL 0`.

- [ ] **Step 10: Commit**

```bash
git add R/seqview.R R/app_modals.R R/app_server.R tests/testthat/test-chimera.R tests/testthat/test-seqview.R
git commit -m "Details window reads the review layer; Assemble chimera panel shows originals"
```


### Task 18: Export without the ITSx run

**Files:**
- Modify: `R/itsx_check.R` (header comment :1-11; delete `itsx_check_command` :34-77 and `run_itsx_check` :79-116; `itsx_check_units` roxygen `@param out` :323-324; delete `itsx_export_check` :431-462)
- Modify: `R/export.R` (`submission_readme` :299-380, `export_summary_df` :382-463, `export_bundle` roxygen and signature :672-691, its ITSx block :834-852)
- Modify: `man/export_bundle.Rd` (regenerated, not edited by hand)
- Modify: `R/app_modals.R` (`export_modal` :809-812, `export_check_summary_ui` :819-857)
- Modify: `R/app_server.R` (`take_export_snap` :1195-1199, `run_export` :1902-1916)
- Test: `tests/testthat/test-itsx-check.R`, `tests/testthat/test-export.R` (:326-336), `tests/testthat/test-export-units.R` (:35-36, :115, :169-173, :184, :293), `tests/testthat/test-export-headers.R` (:187), `tests/testthat/test-nmnh.R` (:218), `tests/testthat/test-app-export.R` (:907)

**Interfaces:**
- Consumes:
  - Task 1: `annotate.itsx_origin TEXT` and `annotate.itsx_scores TEXT` in `new_db()`.
  - Task 5: `itsx_run_finished(out_prefix)`; `itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix)` storing `itsx_warnings`, the one-letter `itsx_origin`, and `itsx_scores` on the annotate row; Task 5's edits to `tests/testthat/test-itsx-check.R` (the two `ALTER TABLE annotate ADD COLUMN itsx_origin/itsx_scores` lines at the end of `mk_itsx_project`, `expect_equal(chk$itsx_origin, "M")` in the real-ITSx test, the summary line and `"not_run(ITSx did not finish)"` in the `itsx_check_one` test, and the block Task 5 appended at the end of the file, which starts with the comment line quoted in Step 3 above its `itsx_out` helper).
  - Task 9: `inst/nextflow/bin/rp_rrna_tools.sh <fasta> <tool> <kingdom> <itsx_taxon> <itsx_opts> <rrna_opts> <cpus>` (exactly 7 arguments; empty tool, kingdom, taxon default to pybarrnap, euk, all), writing `unit.itsx_check.*` (with `unit.itsx_check.summary.txt` and `unit.itsx_check.err`) and `unit.tools_status` (`itsx_check=<exit>`) into the working directory.
  - Existing: `config_container_engine(cfg_lines)`, `config_container_image(cfg_lines)`, `singularity_image_ref(img, project, cfg_lines)`, `chimera_test_con(id = "S1")`, `mp_id_list(ids, max)`, `mp_n(n, noun)`.
- Produces:
  - `export_bundle(path = ".", ids = NULL, units = NULL, group = NULL, out_dir = file.path(path, "out"), label = NULL, markers = TRUE, report = TRUE, template = NULL, run_table2asn = FALSE, zip = TRUE, export_opts = "default", fasta_header = NULL, fasta_header_marker = NULL)`: no `itsx_check` or `itsx_cpu` argument; never runs ITSx, never writes `itsx_warnings`, writes no `itsx_check/` folder and no `<label>_itsx_check.csv`; still warns `ITSx check flags on N of M contig(s) (still exported; see itsx_check column): ...` when stored flags exist.
  - `export_summary_df(con, units)`: last column `itsx_check` = `annotate.itsx_warnings` of the unit, NA written as `not_run(not checked)` (Task 19 adds `review`).
  - `submission_readme(label)`: no `itsx_check` argument; an `ITSx CHECK` section (flag sources, flag legend, `not_run` legend with the advice, batch-size caveat) replaces `FINAL ITSx CHECK`.
  - `export_check_summary_ui(smry, not_done = 0L)`: always one ITSx line `ITSx check: <f> flagged, <c> clean, <n> not checked.` (`not_run(...)` counts as not checked), followed by ` Flagged: <seqid (flags)>, .... See the itsx_check column in the summary CSV.` (bold amber) when any is flagged, else by ` README.txt says how to check them again.` (amber) when any is not checked, else green with a tick unless a submission-risk flag is shown; when `not_done > 0` a further line `<mp_n(not_done, "contig")> left out: checks not finished. Press Update on Annotation Review, then export again.` Task 20 passes `not_done`.
  - Gone: `itsx_check_command`, `run_itsx_check`, `itsx_export_check`, the app input `export_itsx_check`.

Export stops running its own stock ITSx pass and reports the ITSx check already stored for each contig (D4, S9). The summary CSV's `itsx_check` column now always exists and carries `annotate.itsx_warnings`, with NA written as `not_run(not checked)`; Task 19 switches the read to `review_units_v`. The app's "Final ITSx check" box, its progress text, `<label>_itsx_check.csv`, and the `itsx_check/` folder go; the Export done dialog counts flagged, clean, and not checked contigs and names the flagged ones, and takes the not-done count Task 20 will pass (S6.2). The README explains where the flags come from, each `not_run` reason (including texts left by the removed export run in older projects) and what to do, and the batch-size caveat (S9). The real-ITSx test is ported onto `itsx_check_one`: it runs `rp_rrna_tools.sh` in the project's container through a test-only command builder and stays gated on `RIBOPILOT_TEST_PROJECT` (S12).

- [ ] **Step 1: Write the failing tests**

In `tests/testthat/test-itsx-check.R`, replace the two export-wiring tests and the real-ITSx test with the new export tests and the ported real-ITSx test. The region below is the file after Task 5 (its real-ITSx test already expects `"M"`).

Find:

```r
test_that("export_bundle writes a not_run itsx_check column when ITSx is unavailable", {
  skip_if_not_installed("withr")
  dir <- mk_itsx_project()
  withr::local_envvar(PATH = tempdir())   # no ITSx, docker, or singularity here
  g <- suppressWarnings(suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = FALSE, report = FALSE, zip = FALSE)))
  smry <- utils::read.csv(file.path(g, "export_summary.csv"), stringsAsFactors = FALSE)
  expect_true("itsx_check" %in% names(smry))
  expect_match(smry$itsx_check, "^not_run\\(")
  expect_match(paste(readLines(file.path(g, "README.txt")), collapse = "\n"), "FINAL ITSx CHECK")
  expect_equal(attr(g, "summary")$itsx_check, smry$itsx_check)
})

test_that("itsx_check = FALSE leaves the column and README section out", {
  skip_if_not_installed("withr")
  dir <- mk_itsx_project()
  g <- suppressWarnings(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = FALSE, report = FALSE, zip = FALSE, itsx_check = FALSE))
  smry <- utils::read.csv(file.path(g, "export_summary.csv"), stringsAsFactors = FALSE)
  expect_false("itsx_check" %in% names(smry))
  expect_false(file.exists(file.path(g, "export_itsx_check.csv")))
  expect_false(grepl("FINAL ITSx CHECK", paste(readLines(file.path(g, "README.txt")), collapse = "\n")))
})

test_that("a real ITSx run in the project's container over a genuine rDNA record is clean", {
  skip_if_not_installed("withr")
  proj <- Sys.getenv("RIBOPILOT_TEST_PROJECT", "")
  skip_if(!nzchar(proj), "RIBOPILOT_TEST_PROJECT not set")
  src <- file.path(proj, ".sqlite")
  skip_if(!file.exists(src), "test project database not present")
  skip_if(!file.exists(file.path(proj, ".config")), "test project has no .config")
  probe <- RiboPilot:::itsx_check_command(proj, tempdir())
  skip_if(!nzchar(probe$bin), probe$message)
  scon <- DBI::dbConnect(RSQLite::SQLite(), src)
  on.exit(DBI::dbDisconnect(scon), add = TRUE)
  u <- DBI::dbGetQuery(scon, "SELECT an.ID, an.path, an.scaffold FROM annotate an
    JOIN annotations a ON a.ID = an.ID AND a.path = an.path AND a.scaffold = an.scaffold
    WHERE a.element = 'ITS1' LIMIT 1")
  skip_if(nrow(u) == 0, "no ITS1-annotated unit in the local test project")
  seq <- DBI::dbGetQuery(scon, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = unname(as.list(u[1, ])))$sequence[1]
  ann <- DBI::dbGetQuery(scon, "SELECT element, type, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    params = unname(as.list(u[1, ])))

  dir <- withr::local_tempdir()
  file.copy(file.path(proj, ".config"), file.path(dir, ".config"))   # same container
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT, problematic TEXT,
    partial TEXT, reviewed TEXT, chimeric INT, n_features INT, missing TEXT, export_group TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('REAL','Fish sp.','Metazoa;Chordata')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('REAL',1,1,?,?,'linear',0)", params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('REAL',1,1,1,'no','no','yes',0,5,'',NULL)")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('REAL',1)")
  for (i in seq_len(nrow(ann))) {
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('REAL',1,1,?,?,?,?)",
      params = list(ann$element[i], ann$type[i], ann$start[i], ann$stop[i]))
  }
  DBI::dbDisconnect(con)

  g <- suppressWarnings(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = FALSE, report = FALSE, zip = FALSE, itsx_cpu = 4L))
  smry <- utils::read.csv(file.path(g, "export_summary.csv"), stringsAsFactors = FALSE)
  expect_true(smry$itsx_check %in% c("", NA))   # read.csv turns "" into NA
  chk <- utils::read.csv(file.path(g, "export_itsx_check.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(chk$itsx_origin, "M")
  expect_match(chk$itsx_ITS1, "^[0-9]+-[0-9]+$")
  expect_true(file.exists(file.path(g, "itsx_check", "itsx.positions.txt")))
})
```

Replace with:

```r
test_that("export_summary_df carries the stored ITSx check, and NA as not checked", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  vals <- c("", "itsx_no_5.8S;itsx_origin(Fungi)", "not_run(ITSx did not finish)", NA)
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1', 1, ?, 'ACGT', 4)",
      params = list(s))
    DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1', 1, ?)", params = list(s))
    DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = ? WHERE ID = 'S1' AND path = 1 AND scaffold = ?",
      params = list(vals[s], s))
  }
  units <- data.frame(ID = "S1", path = 1L, scaffold = 4:1, seqid = paste0("S1_p1_s", 4:1))
  smry <- RiboPilot:::export_summary_df(con, units)
  expect_equal(smry$itsx_check, c("not_run(not checked)", "not_run(ITSx did not finish)",
    "itsx_no_5.8S;itsx_origin(Fungi)", ""))
  expect_equal(utils::tail(names(smry), 1), "itsx_check")
})

test_that("export reports the stored ITSx check and never runs ITSx", {
  skip_on_os("windows")
  con <- chimera_test_con()
  dir <- dirname(DBI::dbGetInfo(con)$dbname)
  seq <- paste0(strrep("A", 200), strrep("C", 200))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length, topology, ignore)
    VALUES ('S1', 1, 1, ?, 400, 'linear', 0)", params = list(seq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop)
    VALUES ('S1', 1, 1, '18S', 'rRNA', 1, 200), ('S1', 1, 1, '28S', 'rRNA', 201, 400)")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, annotate_lock = 1, problematic = 'no',
    n_features = 2, missing = '', itsx_warnings = 'itsx_chimeric'")
  DBI::dbDisconnect(con)
  # an ITSx first on PATH that leaves a mark when called
  bin <- withr::local_tempdir()
  writeLines(c("#!/bin/sh", sprintf("touch '%s'", file.path(bin, "called")), "exit 1"), file.path(bin, "ITSx"))
  Sys.chmod(file.path(bin, "ITSx"), "0755")
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep))

  expect_warning(g <- suppressMessages(export_bundle(dir, label = "g1", markers = FALSE,
    report = FALSE, zip = FALSE)), "ITSx check flags on 1 of 1 contig(s)", fixed = TRUE)
  expect_false(file.exists(file.path(bin, "called")))
  expect_false(dir.exists(file.path(g, "itsx_check")))
  expect_false(file.exists(file.path(g, "g1_itsx_check.csv")))
  smry <- utils::read.csv(file.path(g, "g1_summary.csv"), stringsAsFactors = FALSE)
  expect_equal(smry$itsx_check, "itsx_chimeric")
  expect_equal(attr(g, "summary")$itsx_check, "itsx_chimeric")
  readme <- paste(readLines(file.path(g, "README.txt")), collapse = "\n")
  expect_false(grepl("FINAL ITSx CHECK", readme, fixed = TRUE))
  for (x in c("ITSx CHECK", "Export does not run ITSx", "not_run(ITSx did not finish)",
              "not_run(ITSx check failed: <why>)", "not_run(not checked before upgrade)",
              "not_run(not checked)", "not_run(no ITSx output)", "not_run(docker not on PATH)",
              "batch size")) {
    expect_match(readme, x, fixed = TRUE)
  }
  con <- RiboPilot:::open_project_db(file.path(dir, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate")$itsx_warnings, "itsx_chimeric")
  for (fn in c("itsx_export_check", "run_itsx_check", "itsx_check_command")) {
    expect_false(exists(fn, envir = asNamespace("RiboPilot"), inherits = FALSE), info = fn)
  }
})

# Run `cmd` in the project's container the way Nextflow runs the pipeline's
# tools: docker or singularity per .config, straight from PATH without one.
# `workdir` is mounted at its own path and is the working directory.
container_run <- function(project, workdir, cmd) {
  cfg <- readLines(file.path(project, ".config"), warn = FALSE)
  engine <- RiboPilot:::config_container_engine(cfg)
  img <- RiboPilot:::config_container_image(cfg)
  wd <- normalizePath(workdir)
  if (is.na(engine)) {
    skip_if(!nzchar(Sys.which("ITSx")), "ITSx not on PATH (no container enabled in .config)")
    return(processx::run(cmd[1], cmd[-1], wd = wd, error_on_status = FALSE))
  }
  skip_if(!nzchar(img), "no process.container image in .config")
  if (identical(engine, "docker")) {
    skip_if(!nzchar(Sys.which("docker")), "docker not on PATH")
    user <- paste(system2("id", "-u", stdout = TRUE), system2("id", "-g", stdout = TRUE), sep = ":")
    return(processx::run("docker", c("run", "--rm", "-u", user, "-v", paste0(wd, ":", wd), "-w", wd,
      img, cmd), wd = wd, error_on_status = FALSE))
  }
  bin <- Sys.which(c("singularity", "apptainer"))
  skip_if(!any(nzchar(bin)), "singularity/apptainer not on PATH")
  processx::run(unname(bin[nzchar(bin)][1]), c("exec", "--bind", wd, "--pwd", wd,
    RiboPilot:::singularity_image_ref(img, project, cfg), cmd), wd = wd, error_on_status = FALSE)
}

test_that("a real stock ITSx run in the project's container over a genuine rDNA record is clean", {
  skip_if_not_installed("withr")
  proj <- Sys.getenv("RIBOPILOT_TEST_PROJECT", "")
  skip_if(!nzchar(proj), "RIBOPILOT_TEST_PROJECT not set")
  src <- file.path(proj, ".sqlite")
  skip_if(!file.exists(src), "test project database not present")
  skip_if(!file.exists(file.path(proj, ".config")), "test project has no .config")
  # read a copy: nothing is written next to the project's database
  db <- file.path(withr::local_tempdir(), "src.sqlite")
  file.copy(src, db)
  scon <- DBI::dbConnect(RSQLite::SQLite(), db)
  withr::defer(DBI::dbDisconnect(scon))
  u <- DBI::dbGetQuery(scon, "SELECT an.ID, an.path, an.scaffold FROM annotate an
    JOIN annotations a ON a.ID = an.ID AND a.path = an.path AND a.scaffold = an.scaffold
    WHERE a.element = 'ITS1' ORDER BY an.ID, an.path, an.scaffold LIMIT 1")
  skip_if(nrow(u) == 0, "no ITS1-annotated unit in the local test project")
  key <- unname(as.list(u[1, ]))
  seq <- DBI::dbGetQuery(scon, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = key)$sequence[1]
  ann <- DBI::dbGetQuery(scon, "SELECT element, type, start, stop FROM annotations
    WHERE ID = ? AND path = ? AND scaffold = ?", params = key)

  wd <- withr::local_tempdir()
  file.copy(system.file("nextflow", "bin", "rp_rrna_tools.sh", package = "RiboPilot"), wd)
  writeLines(c(">REAL.1.1", seq), file.path(wd, "unit.fasta"))
  res <- container_run(proj, wd, c("bash", "rp_rrna_tools.sh", "unit.fasta", "", "", "", "", "", "4"))
  st <- file.path(wd, "unit.tools_status")
  expect_true("itsx_check=0" %in% (if (file.exists(st)) readLines(st)), info = res$stderr)
  pre <- file.path(wd, "unit.itsx_check")
  expect_true(RiboPilot:::itsx_run_finished(pre))

  con <- chimera_test_con("REAL")
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET lineage = 'Metazoa;Chordata'")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('REAL', 1, 1, ?, ?)",
    params = list(seq, nchar(seq)))
  DBI::dbAppendTable(con, "annotations", cbind(data.frame(ID = "REAL", path = 1L, scaffold = 1L), ann))
  chk <- RiboPilot:::itsx_check_one(con, "REAL", 1L, 1L, pre)
  expect_equal(chk$itsx_check, "")
  expect_equal(chk$itsx_origin, "M")
  expect_match(chk$itsx_scores, "M: ", fixed = TRUE)
  expect_match(chk$itsx_ITS1, "^[0-9]+-[0-9]+$")
  expect_equal(unlist(DBI::dbGetQuery(con, "SELECT itsx_warnings, itsx_origin FROM annotate
    WHERE ID = 'REAL' AND path = 1 AND scaffold = 1")), c(itsx_warnings = "", itsx_origin = "M"))
})
```

In `tests/testthat/test-export.R`, replace the done-dialog test.

Find:

```r
test_that("the export done dialog names flagged records and is green only when nothing is flagged", {
  smry <- data.frame(seqid = c("A", "B_p1_s1"), itsx_check = c("", ""),
    submission_risk = c("", "off_target(Actiniaria);contaminant"))
  html <- as.character(RiboPilot:::export_check_summary_ui(smry))
  expect_match(html, "B_p1_s1 (BLAST hits agree on Actiniaria; the best BLAST hit is from another phylum)",
    fixed = TRUE)
  expect_false(grepl("mp-fg-success", html, fixed = TRUE))
  expect_lt(regexpr("Submission-risk", html), regexpr("Final ITSx", html))
  smry$submission_risk <- ""
  expect_match(as.character(RiboPilot:::export_check_summary_ui(smry)), "mp-fg-success", fixed = TRUE)
})
```

Replace with:

```r
test_that("the export done dialog counts the stored ITSx check and is green only when all is clean", {
  smry <- data.frame(seqid = c("A", "B_p1_s1", "C", "D"),
    itsx_check = c("", "itsx_chimeric", "not_run(ITSx did not finish)", "not_run(not checked)"),
    submission_risk = c("", "off_target(Actiniaria);contaminant", "", ""))
  html <- as.character(RiboPilot:::export_check_summary_ui(smry))
  expect_match(html, "B_p1_s1 (BLAST hits agree on Actiniaria; the best BLAST hit is from another phylum)",
    fixed = TRUE)
  expect_match(html, "ITSx check: 1 flagged, 1 clean, 2 not checked. Flagged: B_p1_s1 (itsx_chimeric).",
    fixed = TRUE)
  expect_lt(regexpr("Submission-risk", html), regexpr("ITSx check:", html))
  expect_false(grepl("mp-fg-success", html, fixed = TRUE))
  expect_false(grepl("left out", html, fixed = TRUE))

  # nothing flagged, but a contig not checked: never green
  smry$submission_risk <- ""
  smry$itsx_check <- c("", "", "", "not_run(not checked)")
  html <- as.character(RiboPilot:::export_check_summary_ui(smry, not_done = 2L))
  expect_match(html, "ITSx check: 0 flagged, 3 clean, 1 not checked.", fixed = TRUE)
  expect_false(grepl("mp-fg-success", html, fixed = TRUE))
  expect_match(html, "2 contigs left out: checks not finished. Press Update on Annotation Review", fixed = TRUE)

  smry$itsx_check <- ""
  smry$submission_risk <- c("", "off_target(Actiniaria);contaminant", "", "")
  expect_false(grepl("mp-fg-success", as.character(RiboPilot:::export_check_summary_ui(smry)), fixed = TRUE))
  smry$submission_risk <- ""
  expect_match(as.character(RiboPilot:::export_check_summary_ui(smry)), "mp-fg-success", fixed = TRUE)
  # the Export window no longer offers its own ITSx run
  expect_false(grepl("export_itsx_check", as.character(RiboPilot:::export_modal("g1")), fixed = TRUE))
})
```

In `tests/testthat/test-export-units.R`, the exact summary CSV gains the `itsx_check` column; this fixture's annotate table has no `itsx_warnings` column, so the value is `not_run(not checked)`.

Find:

```r
        "\"variant_group\",\"variant_identity\",\"n_variants\",\"submission_risk\",\"chimera_warning\""),
      "\"S2\",\"S2\",1,1,\"Genus speciesB\",400,2,\"\",\"yes\",\"no\",\"no\",NA,NA,NA,NA,\"\",\"\""))
```

Replace with:

```r
        "\"variant_group\",\"variant_identity\",\"n_variants\",\"submission_risk\",\"chimera_warning\",",
        "\"itsx_check\""),
      paste0("\"S2\",\"S2\",1,1,\"Genus speciesB\",400,2,\"\",\"yes\",\"no\",\"no\",NA,NA,NA,NA,\"\",\"\",",
        "\"not_run(not checked)\"")))
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'options(testthat.progress.max_fails = Inf); devtools::load_all("."); for (f in c("itsx-check", "export", "export-units")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: FAIL, in these four tests only (every other test in the three files passes, and the real-ITSx test skips with "RIBOPILOT_TEST_PROJECT not set"):
- test-itsx-check.R "export_summary_df carries the stored ITSx check, and NA as not checked": `smry$itsx_check` is NULL, and the last column is `chimera_warning` (2 failures).
- test-itsx-check.R "export reports the stored ITSx check and never runs ITSx" (17 failures): no `ITSx check flags` warning; the fake ITSx was called (`called` exists); the `itsx_check/` folder exists; the CSV, the `summary` attribute, and the stored `itsx_warnings` read `not_run(ITSx failed: )` instead of `itsx_chimeric`; the README still has `FINAL ITSx CHECK` and none of the new legend; the three functions still exist.
- test-export.R "the export done dialog counts the stored ITSx check and is green only when all is clean": no `ITSx check:` line (2 failures, the second because `regexpr` gives -1), then the test stops with the error `unused argument (not_done = 2L)`.
- test-export-units.R "bundle file contents are exactly as specified": the CSV header and row have no `itsx_check` (1 failure).

- [ ] **Step 3: Remove the export ITSx run and report the stored check**

In `R/itsx_check.R`, four edits. The header comment describes the pipeline-side check:

Find:

```r
# Final pre-export ITSx check, run the way NCBI's Submission Portal runs it:
# stock ITSx, default parameters, no organism-group restriction. The pipeline's
# own ITSx pass (annotation step) uses the project's itsx_taxon / itsx_opts and
# happens before curation, so it cannot stand in for what the Portal will see.
```

Replace with:

```r
# Stock ITSx check, run the way NCBI's Submission Portal runs it: stock ITSx,
# default parameters, no organism-group restriction. rp_rrna_tools.sh runs it
# beside the tuned ITSx pass, which uses the project's itsx_taxon / itsx_opts and
# so cannot stand in for what the Portal will see, in the annotation step and in
# the recheck. itsx_check_one() stores the result; export only reports it.
```

Delete `itsx_check_command` and `run_itsx_check` (everything from the first roxygen line of `itsx_check_command` down to, not including, the roxygen of `read_itsx_outputs`):

Find:

```r
#' How to invoke ITSx for a project: through its container, or from PATH
#'
#' The pipeline's tools live in the project's container (docker / singularity
#' per `.config`), so the check runs ITSx there. Only a project on the conda
#' executor (no container enabled) falls back to an ITSx on PATH.
#' @param project Project directory (holds `.config`).
#' @param workdir Directory the container must see: the check's input FASTA
#'   and output land here. Bind-mounted at the same path.
#' @param which Lookup for executables (injectable for tests).
#' @return list(bin, prefix, message): `bin` is the executable to call and
#'   `prefix` the arguments that go before ITSx's own; `bin` is "" and
#'   `message` says why when the check cannot run.
#' @noRd
itsx_check_command <- function(project, workdir, which = Sys.which) {
  cfg <- file.path(project, ".config")
  lines <- if (file.exists(cfg)) readLines(cfg, warn = FALSE) else character(0)
  engine <- config_container_engine(lines)
  none <- function(msg) list(bin = "", prefix = character(0), message = msg)
  if (is.na(engine)) {
    bin <- unname(which("ITSx"))
    if (!nzchar(bin)) return(none("ITSx not on PATH (no container enabled in .config)"))
    return(list(bin = bin, prefix = character(0), message = ""))
  }
  img <- config_container_image(lines)
  if (!nzchar(img)) return(none("no process.container image in .config"))
  wd <- normalizePath(workdir, mustWork = FALSE)
  if (identical(engine, "docker")) {
    bin <- unname(which("docker"))
    if (!nzchar(bin)) return(none("docker not on PATH"))
    ids <- vapply(c("-u", "-g"), function(f) {
      tryCatch(trimws(system2("id", f, stdout = TRUE, stderr = FALSE))[1], error = function(e) "")
    }, character(1))
    user <- if (all(nzchar(ids))) c("-u", paste(ids, collapse = ":")) else character(0)
    return(list(bin = bin,
      prefix = c("run", "--rm", user, "-v", shQuote(paste0(wd, ":", wd)), "-w", shQuote(wd), img, "ITSx"),
      message = ""))
  }
  w <- which(c("singularity", "apptainer"))
  w <- w[nzchar(w)]
  if (length(w) < 1) return(none("singularity/apptainer not on PATH"))
  list(bin = unname(w[1]),
    prefix = c("exec", "--bind", shQuote(wd), "--pwd", shQuote(wd), singularity_image_ref(img, project, lines), "ITSx"),
    message = "")
}

#' Run stock ITSx on a FASTA
#'
#' Only output-shaping options are passed (`--preserve`, `--detailed_results`,
#' `--save_regions`, `--graphical`, `--cpu`). Detection parameters are left at
#' their defaults on purpose: that is what the Portal runs.
#' @param fasta Input FASTA; must sit inside the directory `cmd` was built for.
#' @param out_prefix Output prefix (ITSx appends its own suffixes); same rule.
#' @param cpu Threads.
#' @param cmd list(bin, prefix, message) from `itsx_check_command()`.
#' @return list(ok, positions, results, problematic, no_detections, chimeric,
#'   message); `ok = FALSE` and a message when ITSx could not run.
#' @noRd
run_itsx_check <- function(fasta, out_prefix, cpu = 2L, cmd) {
  empty <- list(ok = FALSE, positions = character(0), results = character(0),
    problematic = character(0), no_detections = character(0), chimeric = character(0),
    message = "")
  if (!nzchar(cmd$bin %||% "")) {
    empty$message <- cmd$message %||% "ITSx unavailable"
    return(empty)
  }
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  err <- tempfile("rp_itsx_check_err_")
  on.exit(unlink(err), add = TRUE)
  fasta <- normalizePath(fasta, mustWork = FALSE)
  out_prefix <- file.path(normalizePath(dirname(out_prefix), mustWork = FALSE), basename(out_prefix))
  st <- suppressWarnings(system2(cmd$bin, c(cmd$prefix,
    "-i", shQuote(fasta), "-o", shQuote(out_prefix),
    "--preserve", "T", "--detailed_results", "T", "--save_regions", "none",
    "--graphical", "F", "--cpu", as.integer(cpu)),
    stdout = FALSE, stderr = err))
  if (!identical(as.integer(st), 0L)) {
    empty$message <- paste("ITSx failed:",
      paste(utils::tail(tryCatch(readLines(err, warn = FALSE), error = function(e) ""), 3L),
        collapse = " | "))
    return(empty)
  }
  c(list(ok = TRUE), read_itsx_outputs(out_prefix), list(message = ""))
}

#' Read the files ITSx wrote under `out_prefix`
```

Replace with:

```r
#' Read the files ITSx wrote under `out_prefix`
```

The `out` argument of `itsx_check_units` no longer comes from `run_itsx_check()`:

Find:

```r
#' @param out list from `run_itsx_check()` / `read_itsx_outputs()`, plus `ok`
#'   and `message` when the run itself failed.
```

Replace with:

```r
#' @param out list from `read_itsx_outputs()` plus `ok = TRUE`, or
#'   `list(ok = FALSE, message = <why>)` when ITSx did not finish.
```

Delete `itsx_export_check`, the last function of the file today:

Find:

```r
#' Run the final ITSx check over an export's records
#'
#' Writes the ITSx output under `<gdir>/itsx_check/` and a per-record CSV
#' `<gdir>/<label>_itsx_check.csv`, refreshes `annotate.itsx_warnings`, and
#' returns one row per unit.
#' @param con Project DB connection.
#' @param units data.frame(ID, path, scaffold, seqid) as built by [export_bundle()].
#' @param gdir Bundle directory.
#' @param label Bundle label.
#' @param cpu Threads for ITSx.
#' @param project Project directory; its `.config` decides the container.
#' @param cmd Override for `itsx_check_command()` (tests).
#' @return data.frame as from `itsx_check_units()`.
#' @noRd
itsx_export_check <- function(con, units, gdir, label, cpu = 2L, project = ".", cmd = NULL) {
  if (nrow(units) == 0) return(itsx_check_units(con, units, list(ok = TRUE)))
  seqs <- unit_sequences(con, units)
  cdir <- file.path(gdir, "itsx_check")
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  fa <- file.path(cdir, "records.fasta")
  ok <- !is.na(seqs) & nzchar(seqs)
  writeLines(as.vector(rbind(paste0(">", units$seqid[ok]), seqs[ok])), fa)

  if (is.null(cmd)) cmd <- itsx_check_command(project, cdir)
  out <- run_itsx_check(fa, file.path(cdir, "itsx"), cpu = cpu, cmd = cmd)
  chk <- itsx_check_units(con, units, out)
  itsx_store_warnings(con, units, chk)
  if (isTRUE(out$ok)) {
    utils::write.csv(chk, file.path(gdir, paste0(label, "_itsx_check.csv")), row.names = FALSE)
  }
  chk
}
```

Replace with nothing: delete the whole block, plus the blank line above it, so the function before it (`itsx_check_one`) keeps a single blank line after its closing brace (or before whatever function a later task added below it).

In `R/export.R`, the README loses its argument and its "FINAL ITSx CHECK" section gives way to one that explains the stored check.

Find:

```r
submission_readme <- function(label, itsx_check = TRUE) {
```

Replace with:

```r
submission_readme <- function(label) {
```

Find:

```r
    if (itsx_check) c(
    "FINAL ITSx CHECK",
    "----------------",
    "Stock ITSx (default parameters, no organism group, the Portal's setup) was",
    paste0("run over ", label, ".fasta; see the itsx_check column in ", label, "_summary.csv"),
    paste0("and the per-record detail in ", label, "_itsx_check.csv (raw ITSx output in"),
    "itsx_check/). It ran inside the project's container, and the same flags",
    "were written back to the app's ITSx Warnings column. Empty means ITSx saw",
    "nothing to flag; not_run(...) means ITSx was unavailable. Flags, and the",
    "Portal problem each maps to:",
```

Replace with:

```r
    "ITSx CHECK",
    "----------",
    paste0("The itsx_check column in ", label, "_summary.csv is the stock ITSx check stored"),
    "for each contig, the same as the app's ITSx Warnings column. RiboPilot runs",
    "stock ITSx (default parameters, no organism group, the Portal's setup) on each",
    "contig in the annotation step, and again when Update on Annotation Review",
    "rechecks an edited contig. Export does not run ITSx. Empty means ITSx saw",
    "nothing to flag. Flags, and the Portal problem each maps to:",
```

Find:

```r
    "  itsx_not_detected         ITSx found no ITS; the Portal cannot annotate it",
    "")), collapse = "\n")
}
```

Replace with:

```r
    "  itsx_not_detected         ITSx found no ITS; the Portal cannot annotate it",
    "not_run(...) means the check could not decide; the contig is exported as it",
    "is. To check it again, unlock it on Annotation Review, set its State to Ready",
    "to run, and press Update. The reasons:",
    "  not_run(ITSx did not finish)         ITSx crashed or was stopped, for",
    "                                       example when it ran out of memory",
    "  not_run(ITSx check failed: <why>)    RiboPilot could not read the ITSx output",
    "  not_run(not checked before upgrade)  the contig was locked when the project",
    "                                       was upgraded and had no stored check",
    "  not_run(not checked)                 no check is stored for this contig",
    "  any other not_run(...)               written by an older RiboPilot, by its",
    "                                       export ITSx run or its annotation step,",
    "                                       for example not_run(no ITSx output) or",
    "                                       not_run(docker not on PATH)",
    "ITSx calls can shift with batch size. RiboPilot checks one contig at a time,",
    "while the Portal checks the whole upload together, so the Portal may call a",
    "borderline contig differently.",
    ""), collapse = "\n")
}
```

`export_summary_df` gains the `itsx_check` column. Four edits:

Find:

```r
#' `submission_risk` column (`risk_flags`) so the CSV and HTML report both
#' surface portal-parity risk flags.
```

Replace with:

```r
#' `submission_risk` column (`risk_flags`) so the CSV and HTML report both
#' surface portal-parity risk flags, and an `itsx_check` column: the stored
#' stock ITSx check (annotate.itsx_warnings), `not_run(not checked)` when none
#' is stored.
```

Find:

```r
  nts <- fetch("SELECT ID, path, scaffold, curate_flags FROM annotate WHERE ID IN (%s)")
```

Replace with:

```r
  nts <- fetch("SELECT ID, path, scaffold, curate_flags FROM annotate WHERE ID IN (%s)")
  itx <- fetch("SELECT ID, path, scaffold, itsx_warnings FROM annotate WHERE ID IN (%s)")
```

Find:

```r
  nti <- if (nrow(nts)) match(ukey, unit_key(nts$ID, nts$path, nts$scaffold)) else rep(NA_integer_, n)
```

Replace with:

```r
  nti <- if (nrow(nts)) match(ukey, unit_key(nts$ID, nts$path, nts$scaffold)) else rep(NA_integer_, n)
  itsx <- if (nrow(itx)) itx$itsx_warnings[match(ukey, unit_key(itx$ID, itx$path, itx$scaffold))] else rep(NA_character_, n)
  itsx[is.na(itsx)] <- "not_run(not checked)"
```

Find:

```r
    chimera_warning = ifelse(chm_label %in% "minor", chm_reason, ""),
    stringsAsFactors = FALSE)
```

Replace with:

```r
    chimera_warning = ifelse(chm_label %in% "minor", chm_reason, ""),
    itsx_check = itsx,
    stringsAsFactors = FALSE)
```

`export_bundle` loses `itsx_check` and `itsx_cpu` (roxygen and signature), and reports the stored check instead of running ITSx.

Find:

```r
#' @param itsx_check Run stock ITSx (default parameters, no organism group) over
#'   the exported records, the way NCBI's Submission Portal does, and report
#'   what it would flag in the `itsx_check` column of the summary CSV and in
#'   `<label>_itsx_check.csv`. Warn-only. ITSx runs inside the project's
#'   container (docker / singularity per `.config`); a conda-executor project
#'   uses the ITSx on PATH. Otherwise the column reads `not_run(...)`.
#' @param itsx_cpu Threads for the ITSx check.
#' @return (invisibly) the group output directory, with the summary data.frame
#'   attached as attribute `summary`.
#' @details A unit is exported when its sample's assembly and the unit's
#'   annotation are both locked, its contig is kept, and it is not excluded
#'   from export.
```

Replace with:

```r
#' @return (invisibly) the group output directory, with the summary data.frame
#'   attached as attribute `summary`.
#' @details A unit is exported when its sample's assembly and the unit's
#'   annotation are both locked, its contig is kept, and it is not excluded
#'   from export. The `itsx_check` column of the summary CSV reports the stock
#'   ITSx check stored for each contig (`not_run(not checked)` when none is
#'   stored); export does not run ITSx.
```

Find:

```r
                          fasta_header_marker = NULL, itsx_check = TRUE,
                          itsx_cpu = 2L) {
```

Replace with:

```r
                          fasta_header_marker = NULL) {
```

Find:

```r
  # final ITSx check: stock run, as the Portal does it (warn-only)
  if (isTRUE(itsx_check)) {
    chk <- itsx_export_check(con, units, gdir, label, cpu = itsx_cpu, project = path)
    smry$itsx_check <- chk$itsx_check[match(smry$seqid, chk$seqid)]
    not_run <- startsWith(smry$itsx_check, "not_run(")
    if (any(not_run)) {
      message("ITSx check skipped: ", sub("^not_run\\((.*)\\)$", "\\1", smry$itsx_check[not_run][1]))
    }
    flagged <- which(nzchar(smry$itsx_check) & !not_run)
    if (length(flagged) > 0) {
      warning("ITSx check flags on ", length(flagged), " of ", nrow(smry),
        " contig(s) (still exported; see itsx_check column): ",
        paste(sprintf("%s[%s]", smry$seqid[flagged], smry$itsx_check[flagged]), collapse = ", "),
        call. = FALSE)
    }
  }
  smry <- .meta_join(smry, con)
  utils::write.csv(smry, file.path(gdir, paste0(label, "_summary.csv")), row.names = FALSE)
  writeLines(submission_readme(label, itsx_check = isTRUE(itsx_check)), file.path(gdir, "README.txt"))
```

Replace with:

```r
  # stored stock ITSx check, what the Portal will see (warn-only)
  flagged <- which(nzchar(smry$itsx_check) & !startsWith(smry$itsx_check, "not_run("))
  if (length(flagged) > 0) {
    warning("ITSx check flags on ", length(flagged), " of ", nrow(smry),
      " contig(s) (still exported; see itsx_check column): ",
      paste(sprintf("%s[%s]", smry$seqid[flagged], smry$itsx_check[flagged]), collapse = ", "),
      call. = FALSE)
  }
  smry <- .meta_join(smry, con)
  utils::write.csv(smry, file.path(gdir, paste0(label, "_summary.csv")), row.names = FALSE)
  writeLines(submission_readme(label), file.path(gdir, "README.txt"))
```

In `R/app_modals.R`, the Export window drops the "Final ITSx check" box.

Find:

```r
    opts_help("One HTML page with per-element alignments across the group."),
    mp_checkbox("export_itsx_check", "Final ITSx check", value = flag("export_itsx_check", TRUE)),
    opts_help("Runs stock ITSx over the exported records the way the NCBI Submission ",
      "Portal does and lists what it would flag in the summary CSV. Can take a few minutes.")),
```

Replace with:

```r
    opts_help("One HTML page with per-element alignments across the group.")),
```

Replace `export_check_summary_ui` (roxygen and body):

Find:

```r
#' Reads the summary attached to [export_bundle()]'s return value; names the
#' records the submission-risk gate flagged, and why, and says how many the
#' ITSx check flagged, so a curator sees it before opening the CSV.
#' @param smry data.frame from `export_summary_df()` (possibly with `itsx_check`).
#' @noRd
export_check_summary_ui <- function(smry) {
  if (is.null(smry) || nrow(smry) == 0) return(NULL)
  n <- nrow(smry)
  risky <- nzchar(smry$submission_risk %||% "")
  lines <- list()
  if (any(risky)) {
    lines[[1]] <- tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
      sprintf("Submission-risk flags on %d of %d record(s): %s. See the submission_risk column in the summary CSV.",
        sum(risky), n, mp_id_list(sprintf("%s (%s)", smry$seqid[risky],
          risk_reasons(smry$submission_risk[risky])), max = 10)))
  }
  if (!is.null(smry$itsx_check)) {
    not_run <- startsWith(smry$itsx_check, "not_run(")
    if (all(not_run)) {
      lines[[length(lines) + 1]] <- tags$p(class = "mp-fg-warning",
        "ITSx check not run: ", sub("^not_run\\((.*)\\)$", "\\1", smry$itsx_check[1]), ".")
    } else {
      flagged <- sum(nzchar(smry$itsx_check) & !not_run)
      lines[[length(lines) + 1]] <- if (flagged == 0) {
        tags$p(class = if (!any(risky)) "mp-fg-success", if (!any(risky)) HTML("&#10003; "),
          sprintf("Final ITSx check: nothing flagged on %d record(s).", n))
      } else {
        tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
          sprintf("Final ITSx check flagged %d of %d record(s): %s. See the itsx_check column in the summary CSV.",
            flagged, n, paste(smry$seqid[nzchar(smry$itsx_check) & !not_run], collapse = ", ")))
      }
    }
  }
  if (length(lines) == 0) return(NULL)
  do.call(tagList, lines)
}
```

Replace with:

```r
#' Reads the summary attached to [export_bundle()]'s return value; names the
#' records the submission-risk gate flagged, and why, and counts the stored
#' ITSx check (flagged, clean, not checked), naming the flagged records, so a
#' curator sees both before opening the CSV.
#' @param smry data.frame from `export_summary_df()`.
#' @param not_done Contigs of the group left out because their checks are not
#'   finished.
#' @noRd
export_check_summary_ui <- function(smry, not_done = 0L) {
  if (is.null(smry) || nrow(smry) == 0) return(NULL)
  n <- nrow(smry)
  risky <- nzchar(smry$submission_risk %||% "")
  lines <- list()
  if (any(risky)) {
    lines[[1]] <- tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
      sprintf("Submission-risk flags on %d of %d record(s): %s. See the submission_risk column in the summary CSV.",
        sum(risky), n, mp_id_list(sprintf("%s (%s)", smry$seqid[risky],
          risk_reasons(smry$submission_risk[risky])), max = 10)))
  }
  itsx <- smry$itsx_check
  unchecked <- startsWith(itsx, "not_run(")
  flagged <- nzchar(itsx) & !unchecked
  counts <- sprintf("ITSx check: %d flagged, %d clean, %d not checked.",
    sum(flagged), sum(!flagged & !unchecked), sum(unchecked))
  lines[[length(lines) + 1]] <- if (any(flagged)) {
    tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
      sprintf("%s Flagged: %s. See the itsx_check column in the summary CSV.", counts,
        mp_id_list(sprintf("%s (%s)", smry$seqid[flagged], itsx[flagged]), max = 10)))
  } else if (any(unchecked)) {
    tags$p(class = "mp-fg-warning", sprintf("%s README.txt says how to check them again.", counts))
  } else {
    tags$p(class = if (!any(risky)) "mp-fg-success", if (!any(risky)) HTML("&#10003; "), counts)
  }
  if (not_done > 0) {
    lines[[length(lines) + 1]] <- tags$p(class = "mp-fg-warning",
      sprintf("%s left out: checks not finished. Press Update on Annotation Review, then export again.",
        mp_n(not_done, "contig")))
  }
  do.call(tagList, lines)
}
```

In `R/app_server.R`, the Export snapshot and `run_export` stop naming the removed input.

Find:

```r
    export_snap(lapply(stats::setNames(nm = c("export_group_sel", "export_template",
      "fasta_header", "fasta_header_marker", "export_markers", "export_report",
      "export_itsx_check")), function(id) input[[id]]))
```

Replace with:

```r
    export_snap(lapply(stats::setNames(nm = c("export_group_sel", "export_template",
      "fasta_header", "fasta_header_marker", "export_markers", "export_report")),
      function(id) input[[id]]))
```

Find:

```r
  run_export <- function(u, grp) {
    detail <- if (isTRUE(input$export_itsx_check)) {
      "Writing files and running the final ITSx check. Hold tight, this can take a few minutes."
    } else {
      "Writing files. Hold tight."
    }
    gdir <- withProgress(message = "Exporting...", detail = detail, value = 0.5, tryCatch(
```

Replace with:

```r
  run_export <- function(u, grp) {
    gdir <- withProgress(message = "Exporting...", detail = "Writing files. Hold tight.", value = 0.5, tryCatch(
```

Find:

```r
        report = isTRUE(input$export_report),
        itsx_check = isTRUE(input$export_itsx_check)),
```

Replace with:

```r
        report = isTRUE(input$export_report)),
```

Existing tests that pinned the removed pieces, updated in this step:

`tests/testthat/test-itsx-check.R`: the header comment, then three deletions. These tests call `run_itsx_check` and `itsx_check_command`, which no longer exist; the `singularity_image_ref` test (which those functions used) stays.

Find:

```r
# Final pre-export ITSx check: parsers over stock-ITSx output, the pure flag
# function, and the export_bundle wiring.
```

Replace with:

```r
# Stock ITSx check: parsers over stock-ITSx output, the pure flag function,
# storing the result, and what export reports from it.
```

Find:

```r
test_that("run_itsx_check reports not-run when there is no way to call ITSx", {
  out <- RiboPilot:::run_itsx_check("x.fa", tempfile(),
    cmd = list(bin = "", prefix = character(0), message = "docker not on PATH"))
  expect_false(out$ok)
  expect_equal(out$message, "docker not on PATH")
})

# fake executable lookup: everything named in `have` exists
fake_which <- function(have) function(x) stats::setNames(ifelse(x %in% have, paste0("/bin/", x), ""), x)

write_cfg <- function(dir, ...) writeLines(c(...), file.path(dir, ".config"))

test_that("itsx_check_command runs ITSx through docker when .config enables it", {
  dir <- withr::local_tempdir()
  write_cfg(dir, "docker {", "  enabled = true", "}", "process {", "  container = 'ribopilot:test'", "}")
  wd <- file.path(dir, "out", "itsx_check")
  cmd <- RiboPilot:::itsx_check_command(dir, wd, which = fake_which("docker"))
  expect_equal(cmd$bin, "/bin/docker")
  expect_equal(cmd$prefix[1:2], c("run", "--rm"))
  expect_true(shQuote(paste0(wd, ":", wd)) %in% cmd$prefix)   # workdir mounted at itself
  expect_equal(utils::tail(cmd$prefix, 2), c("ribopilot:test", "ITSx"))
  # engine enabled but nothing to run it with
  expect_equal(RiboPilot:::itsx_check_command(dir, wd, which = fake_which("ITSx"))$bin, "")
})

test_that("itsx_check_command uses singularity or apptainer for a singularity config", {
  dir <- withr::local_tempdir()
  sif <- file.path(dir, "rp.sif"); file.create(sif)
  write_cfg(dir, "singularity {", "  enabled = true", "}", "process {", paste0("  container = '", sif, "'"), "}")
  wd <- file.path(dir, "wd")
  cmd <- RiboPilot:::itsx_check_command(dir, wd, which = fake_which("apptainer"))
  expect_equal(cmd$bin, "/bin/apptainer")
  expect_equal(cmd$prefix, c("exec", "--bind", shQuote(wd), "--pwd", shQuote(wd), sif, "ITSx"))
  expect_match(RiboPilot:::itsx_check_command(dir, wd, which = fake_which("ITSx"))$message, "singularity")
})

test_that("itsx_check_command falls back to PATH only for a container-less config", {
  dir <- withr::local_tempdir()
  write_cfg(dir, "process {", "  executor = 'local'", "}")
  cmd <- RiboPilot:::itsx_check_command(dir, dir, which = fake_which("ITSx"))
  expect_equal(cmd$bin, "/bin/ITSx")
  expect_equal(cmd$prefix, character(0))
  expect_match(RiboPilot:::itsx_check_command(dir, dir, which = fake_which("docker"))$message, "not on PATH")
  # no .config at all behaves the same
  expect_equal(RiboPilot:::itsx_check_command(tempfile(), dir, which = fake_which("ITSx"))$bin, "/bin/ITSx")
})

# Minimal project for the export wiring (same shape as test-export-units.R)
```

Replace with:

```r
# Minimal project for itsx_check_one
```

Find:

```r
test_that("itsx_check_command hands singularity a docker:// URI for a Docker Hub image", {
  dir <- withr::local_tempdir()
  write_cfg(dir, "singularity {", "  enabled = true", "}", "process {", "  container = 'macguigand/ribopilot:0.1'", "}")
  cmd <- RiboPilot:::itsx_check_command(dir, dir, which = fake_which("singularity"))
  expect_equal(cmd$prefix[6], "docker://macguigand/ribopilot:0.1")
})

test_that("itsx_check_one stores stock-ITSx flags on the unit's annotate row", {
```

Replace with:

```r
test_that("itsx_check_one stores stock-ITSx flags on the unit's annotate row", {
```

The next Find ends on the first comment line of the block Task 5 appended, which stays.

Find:

```r
test_that("export refreshes itsx_warnings on the annotate row", {
  skip_if_not_installed("withr")
  dir <- mk_itsx_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_warnings TEXT")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_chimeric'")
  DBI::dbDisconnect(con)
  withr::local_envvar(PATH = tempdir())
  suppressWarnings(suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = FALSE, report = FALSE, zip = FALSE)))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  expect_match(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate")$itsx_warnings, "^not_run\\(")
})

test_that("the ITSx container call keeps a folder with a space as one argument", {
  skip_on_os("windows")
  root <- file.path(withr::local_tempdir(), "Batch 1")
  wd <- file.path(root, "itsx_check")
  dir.create(wd, recursive = TRUE)
  write_cfg(root, "docker {", "  enabled = true", "}", "process {", "  container = 'ribopilot:test'", "}")
  seen <- file.path(root, "args.txt")
  fake <- file.path(root, "docker")
  writeLines(c("#!/bin/sh", sprintf("for a in \"$@\"; do printf '%%s\\n' \"$a\"; done > '%s'", seen),
    "exit 1"), fake)
  Sys.chmod(fake, "0755")
  cmd <- RiboPilot:::itsx_check_command(root, wd, which = function(x) stats::setNames(fake, x))
  fa <- file.path(wd, "records.fasta")
  writeLines(c(">a", "ACGT"), fa)
  RiboPilot:::run_itsx_check(fa, file.path(wd, "itsx"), cmd = cmd)
  a <- readLines(seen)
  wdn <- normalizePath(wd)
  expect_equal(a[which(a == "-v") + 1], paste0(wdn, ":", wdn))
  expect_equal(a[which(a == "-w") + 1], wdn)
})

# stock-ITSx output files under `pre`; finished = FALSE leaves the summary
```

Replace with:

```r
# stock-ITSx output files under `pre`; finished = FALSE leaves the summary
```

`tests/testthat/test-export-units.R`: four callers pass the removed `itsx_check = FALSE` (an unused argument now).

Find:

```r
    markers = FALSE, report = FALSE, zip = FALSE, itsx_check = FALSE, ...))
```

Replace with:

```r
    markers = FALSE, report = FALSE, zip = FALSE, ...))
```

Find:

```r
      markers = FALSE, report = FALSE, zip = FALSE, itsx_check = FALSE, ids = "S1")),
```

Replace with:

```r
      markers = FALSE, report = FALSE, zip = FALSE, ids = "S1")),
```

Find:

```r
    markers = TRUE, report = FALSE, zip = FALSE, itsx_check = FALSE, label = "notax")))
```

Replace with:

```r
    markers = TRUE, report = FALSE, zip = FALSE, label = "notax")))
```

Find:

```r
    markers = FALSE, report = FALSE, itsx_check = FALSE, zip = TRUE, ...)))
```

Replace with:

```r
    markers = FALSE, report = FALSE, zip = TRUE, ...)))
```

`tests/testthat/test-export-headers.R` (same reason):

Find:

```r
    zip = FALSE, itsx_check = FALSE, label = "vou", fasta_header = "{seqid} [specimen-voucher={voucher}]")),
```

Replace with:

```r
    zip = FALSE, label = "vou", fasta_header = "{seqid} [specimen-voucher={voucher}]")),
```

`tests/testthat/test-nmnh.R` (same reason):

Find:

```r
    report = FALSE, zip = FALSE, itsx_check = FALSE, label = "g1",
```

Replace with:

```r
    report = FALSE, zip = FALSE, label = "g1",
```

`tests/testthat/test-app-export.R` sets the input the Export window no longer has:

Find:

```r
      export_markers = FALSE, export_report = FALSE, export_itsx_check = FALSE)
```

Replace with:

```r
      export_markers = FALSE, export_report = FALSE)
```

Check that nothing still names the removed pieces:

Run: `grep -rn "itsx_check_command\|run_itsx_check\|itsx_export_check\|export_itsx_check\|itsx_cpu\|itsx_check = \(TRUE\|FALSE\)\|FINAL ITSx\|Final ITSx\|_itsx_check.csv" R/ inst/ tests/testthat/`
Expected: exactly four lines, all in the new tests: in tests/testthat/test-itsx-check.R the `g1_itsx_check.csv` line, the `FINAL ITSx CHECK` line, and the `for (fn in c("itsx_export_check", "run_itsx_check", "itsx_check_command"))` line; in tests/testthat/test-export.R the `export_modal("g1")` line. (vignettes/GenBank-Submission.Rmd still describes the old run; Task 23 rewrites it.)

- [ ] **Step 4: Regenerate the help page of export_bundle**

Its arguments and details changed in Step 3 and it is exported. Never run `devtools::document()` (NAMESPACE is hand-maintained).

Run: `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`
Then: `git status --short man/`
Expected: only ` M man/export_bundle.Rd`; its usage ends `fasta_header_marker = NULL\n)`, the `itsx_check` and `itsx_cpu` items are gone, and the details end with `\verb{not_run(not checked)} when none is stored); export does not run ITSx.` If any other man/ file changed, it is not this task's: restore it with `git checkout -- man/<that file>`.

- [ ] **Step 5: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("itsx-check", "export", "export-units", "export-headers", "nmnh", "export-report", "export-opts", "utils", "specimen-app", "docs-site", "app-export")) testthat::test_file(file.path("tests/testthat", paste0("test-", f, ".R")))'`
Expected: PASS, every file `FAIL 0`. test-itsx-check.R reports one skip, the real-ITSx test ("RIBOPILOT_TEST_PROJECT not set"), and has 7 fewer tests than before this task: 10 removed (the old real-ITSx test among them) and 3 added (the two export tests and the ported real-ITSx test). test-export-report.R keeps its skip_on_cran skip. Tallies when this was written on top of Task 5 only (earlier tasks add tests to some of these files, so counts can be higher): itsx-check `SKIP 1 | PASS 110` (22 tests), export `PASS 86`, export-units `PASS 53`, export-headers `PASS 60`, nmnh `PASS 116`, export-report `SKIP 1 | PASS 18`, export-opts `PASS 43`, utils `PASS 31`, specimen-app `PASS 99`, docs-site `PASS 237`, app-export `PASS 135`.

When a project with ITS1-annotated contigs and a working container config is at hand (for example `$HOME/Documents/DJMCompBio/RiboPilot_test/2026_10_08`, image `ribopilot:devel`; the test reads a copy of its database and writes nothing into the project), also run the ported real-ITSx test:

Run: `RIBOPILOT_TEST_PROJECT=$HOME/Documents/DJMCompBio/RiboPilot_test/2026_10_08 Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R", reporter = "progress")'`
Expected: PASS with no skip; the container run takes about a minute. It pins Task 5's and Task 9's behaviour on real stock ITSx output (finished-run detection, one-letter origin `M`, a scores string), so it also passes before this task's Step 3.

- [ ] **Step 6: Commit**

```bash
git add R/itsx_check.R R/export.R R/app_modals.R R/app_server.R man/export_bundle.Rd tests/testthat/test-itsx-check.R tests/testthat/test-export.R tests/testthat/test-export-units.R tests/testthat/test-export-headers.R tests/testthat/test-nmnh.R tests/testthat/test-app-export.R
git commit -m "Export reports the stored ITSx check; no export ITSx run"
```


### Task 19: Export reads the review layer

**Files:**
- Create: `tests/testthat/test-export-review.R`, `tests/testthat/test-export-views.R`
- Modify: `R/export.R` (`feature_table_text_con`, `export_summary_df`, `exportable_units`, `export_left_out`, `feature_map_svg`, `export_report_html`, `export_bundle`, `export_project`; as Tasks 14, 15, and 18 left them)
- Modify: `R/export_opts.R` (`export_header_data`, repo :133-185)
- Modify: `R/export_report.R` (`report_elements`, `gather_feature_seqs`)
- Modify: `inst/export_report.Rmd` (:23)
- Modify: `R/app_server.R` (`screen_flagged` and its two Export-flow calls, repo :1839, :1882-1887, :1893)
- Modify: `man/export_bundle.Rd`, `man/export_project.Rd` (regenerated)
- Modify (fixtures): `tests/testthat/test-export.R` (:12-29, :31-55, :87-104, :106-130, :258-263), `tests/testthat/test-export-units.R` (`mk_project`, the exact-CSV test), `tests/testthat/test-export-headers.R` (`mk_hdr_project`), `tests/testthat/test-export-opts.R` (`mk_opts_db`), `tests/testthat/test-export-report.R` (`mk_report_project`), `tests/testthat/test-nmnh.R` (the NMNH export test)
- Test: `tests/testthat/test-export-review.R`, `tests/testthat/test-export-views.R`

**Interfaces:**
- Consumes: Task 1: views `review_units_v`, `review_annotations_v`, `review_screen_v`, `review_chimera_v`; test helpers `fixture_review_schema(con)`, `fixture_mark_done(con)`, `mk_app_project()`. Task 2: `REVIEW_EFF_PARTIALS`, `review_copy(con, id, path, scaffold)`. Task 4: `review_trim(con, id, path, scaffold, from, to)`. Task 15: `unit_trimmed(con, id, path, scaffold)` reading `review_units_v`. Task 18: `export_bundle()` without `itsx_check`/`itsx_cpu`, `export_summary_df()` with its `itsx_check` column (`itx` fetch).
- Produces:
  - `export_summary_df(con, units, review = FALSE)`: `review = TRUE` reads length, n_features, missing, chimeric, curate_flags, itsx_warnings, and sequence from `review_units_v`, features from `review_annotations_v`, the screen verdict and LCA from `review_screen_v`, the chimera label and reason from `review_chimera_v`; reviewed, partial, problematic from `annotate`; the variant facts (variant_class, variant_group, variant_identity, n_variants) from the ORIGINAL `screen` in both modes. `review = FALSE` is unchanged (originals).
  - `exportable_units(con)`: data.frame(ID, path, scaffold, export_group, done); `done` is 1 when `review_units_v.check_switch = 2`, else 0. Rows are every exportable unit (both locks, kept, not excluded), done or not.
  - `export_bundle(...)`: same signature as after Task 18. Drops units with `done == 0` after the group / ids / units filters; when none remain and some were dropped, `stop("No contigs to export: the checks of <N> contig(s) are not finished. Press Update on Annotation Review, then export again.")`; keeps the old "No contigs to export: a contig needs its assembly ..." stop when nothing was dropped. SeqIDs still count every `exportable_units()` row. The console message adds `"\n  <N> contig(s) left out: checks not finished"` when N > 0. `attr(gdir, "not_done")` = data.frame(ID, path, scaffold) of the dropped units (0 rows when none). Calls `export_summary_df(con, units, review = TRUE)`. Still one `open_project_db` per call.
  - `export_left_out(con, group)`: data.frame(ID, path, scaffold, why); also lists exportable units that are not done, with why `'checks not finished'`; a dropped contig reads `'contig dropped'`.
  - `export_project(path, out_dir)`: a group with no done unit is skipped with the message `"Skipped export group <g>: the checks of its contigs are not finished. Press Update on Annotation Review, then export again."`; returns (invisibly) the groups written.
  - `feature_table_text_con`, `export_header_data`, `report_elements`, `gather_feature_seqs`, `feature_map_svg`, `export_report_html` (default `smry`), and `inst/export_report.Rmd` read per-contig data only through the views (effective partials via `REVIEW_EFF_PARTIALS`); marks from `annotate`, sample columns from `samples`. `unit_fasta_text` stays on `assemblies`.
  - App (`R/app_server.R`): `screen_flagged(u, review = FALSE)`; `screen_then_proceed` and the `export_screen_skip` observer pass `review = TRUE`. `lock_flagged` still calls `screen_flagged(u)` (Task 21 adds its review argument).
  - `tests/testthat/test-export-views.R`: `EXPORT_VIEW_FILES`, `EXPORT_ORIG_TABLES`, `EXPORT_CHECK_COLS`, `EXPORT_VIEW_ALLOW`, `export_strings(f)`, `sql_tables(txt)`; Task 20 extends the allowlist and file list.

Export now writes what Annotation Review shows: every per-contig value comes from the review views, so a trimmed or edited contig ships its copy's sequence, features, effective partial ends, and check results (S9, S4.4). The original tables stay the Assemble step's: `export_summary_df()` gains `review` for the lock check (S8), and its variant facts stay on the original `screen` (an Assemble fact, D1). Contigs whose checks are not done are left out of the bundle as a backstop (S6.2): `exportable_units()` gains `done`, `export_bundle()` drops not-done units after its filters but names records from every exportable unit, `export_left_out()` reports them, and `export_project()` skips a group with none done. A static test pins the rule that export code names original per-contig tables only in the allowed places (S12). The hand-built export fixtures gain the review tables and views (`fixture_review_schema`) and, where they export, done checks (`fixture_mark_done`).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-export-review.R`:

```r
# Export reads the review layer (S9) and leaves out contigs whose checks are
# not done (S6.2).

# A new_db project with one minus-strand contig S1.1.1, both locks on, its
# original checks done, in export group g1. The 28S sits at the left end and
# the 18S at the right end, so both are cut by a trim of both ends.
mk_review_export <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"),
    file.path(d, "mapping.csv"), row.names = FALSE)
  writeLines(c(">18S", "ACGT"), file.path(d, "x.fasta"))
  RiboPilot::new_db(db_path = file.path(d, ".sqlite"), mapping_fn = file.path(d, "mapping.csv"),
    check_reads = FALSE, seeds_db = file.path(d, "x.fasta"), labels_db = file.path(d, "x.fasta"))
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con), envir = env)
  set.seed(7)
  seq <- paste(sample(c("A", "C", "G", "T"), 600, replace = TRUE), collapse = "")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 1, 'linear', 600, ?, 0)", params = list(seq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length, strand)
    VALUES ('S1', 1, 1, '28S', 'rRNA', 1, 200, 200, '-'), ('S1', 1, 1, '18S', 'rRNA', 401, 600, 200, '-')")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, itsx_warnings = '', annotate_lock = 1,
    problematic = 'no', reviewed = 'yes', n_features = 2, missing = '', chimeric = 0, export_group = 'g1'")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  list(con = con, dir = d, seq = seq)
}

review_bundle <- function(dir, ...) {
  suppressWarnings(suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
    report = FALSE, zip = FALSE, ...)))
}

test_that("a trimmed minus-strand copy exports its own sequence, features, and partial ends", {
  p <- mk_review_export()
  con <- p$con
  # copy window 51..487 of the original; 28S cut on the left, 18S on the right
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 51L, 487L)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, itsx_warnings = 'itsx_no_ITS1', chimeric = 1")
  cp <- substr(p$seq, 51, 487)

  g <- review_bundle(p$dir, group = "g1", markers = TRUE)
  fa <- readLines(file.path(g, "samples", "S1.fasta"))
  expect_equal(fa[2], cp)
  expect_match(fa[1], "28S ribosomal RNA gene and 18S ribosomal RNA gene, partial sequence.", fixed = TRUE)
  # minus strand: written stop first; the cut ends are partial in the feature's own orientation
  expect_equal(readLines(file.path(g, "samples", "S1.tbl")),
    c(">Feature S1", "150\t>1\trRNA", "\t\t\tproduct\t28S ribosomal RNA",
      "<437\t351\trRNA", "\t\t\tproduct\t18S ribosomal RNA"))
  expect_equal(readLines(file.path(g, "markers", "18S", "S1_18S.fasta"))[2], substr(cp, 351, 437))
  smry <- utils::read.csv(file.path(g, "g1_summary.csv"), stringsAsFactors = FALSE)
  expect_equal(smry$length, 437L)
  expect_equal(smry$itsx_check, "itsx_no_ITS1")
  expect_match(smry$submission_risk, "chimeric(ITSx)", fixed = TRUE)

  # the originals stay as the main pipeline left them
  expect_equal(DBI::dbGetQuery(con, "SELECT start, stop FROM annotations ORDER BY start")$stop, c(200L, 600L))
  o <- RiboPilot:::export_summary_df(con, data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1"))
  expect_equal(o$length, 600L)
  expect_equal(o$itsx_check, "")
  expect_false(grepl("chimeric", o$submission_risk))
  h <- RiboPilot:::export_header_data(con, data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1"))
  expect_equal(h$length, 437L)
})

test_that("an untrimmed contig without a copy marks its assembly ends partial", {
  p <- mk_review_export()
  g <- review_bundle(p$dir, group = "g1")
  expect_equal(readLines(file.path(g, "samples", "S1.tbl")),
    c(">Feature S1", "200\t>1\trRNA", "\t\t\tproduct\t28S ribosomal RNA",
      "<600\t401\trRNA", "\t\t\tproduct\t18S ribosomal RNA"))
  expect_equal(readLines(file.path(g, "samples", "S1.fasta"))[2], p$seq)
})

test_that("contigs whose checks are not done are left out, and record names do not change", {
  p <- mk_review_export()
  con <- p$con
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 2, 'linear', 600, ?, 0)", params = list(p$seq))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, itsx_warnings,
    annotate_lock, problematic, export_group) VALUES ('S1', 1, 2, 2, NULL, 1, 'no', 'g1')")
  ex <- RiboPilot:::exportable_units(con)
  expect_equal(ex$done, c(1L, 0L))
  expect_equal(RiboPilot:::export_left_out(con, "g1")$why, "checks not finished")

  msg <- testthat::capture_messages(g <- suppressWarnings(export_bundle(p$dir,
    out_dir = file.path(p$dir, "out"), group = "g1", report = FALSE, zip = FALSE, markers = FALSE)))
  expect_true(any(grepl("1 contig(s) left out: checks not finished", msg, fixed = TRUE)))
  # S1 still has two exportable contigs, so the one written keeps its suffix
  expect_equal(basename(Sys.glob(file.path(g, "samples", "*.fasta"))), "S1_p1_s1.fasta")
  expect_equal(attr(g, "not_done"), data.frame(ID = "S1", path = 1L, scaffold = 2L))
  expect_equal(nrow(attr(g, "summary")), 1L)

  # a trimmed copy is not done until its recheck finishes
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 51L, 487L)
  expect_error(review_bundle(p$dir, group = "g1"),
    "the checks of 2 contig(s) are not finished. Press Update on Annotation Review", fixed = TRUE)
  expect_equal(RiboPilot:::export_left_out(con, "g1")$why, rep("checks not finished", 2))

  # nothing locked: the old message, not the checks one
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0")
  expect_error(review_bundle(p$dir, group = "g1"), "needs its assembly and its annotation locked")
  expect_equal(RiboPilot:::export_left_out(con, "g1")$why, rep("annotation not locked", 2))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1")
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE scaffold = 2")
  expect_equal(RiboPilot:::export_left_out(con, "g1")$why, c("checks not finished", "contig dropped"))
})

test_that("export_project skips a group whose contigs are not done", {
  p <- mk_review_export()
  con <- p$con
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 2, 'linear', 600, ?, 0)", params = list(p$seq))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, itsx_warnings,
    annotate_lock, problematic, export_group) VALUES ('S1', 1, 2, 1, NULL, 1, 'no', 'g2')")
  seen <- character(0)
  testthat::local_mocked_bindings(export_bundle = function(path, group = NULL, ...) {
    seen <<- c(seen, group); invisible(path)
  }, .package = "RiboPilot")
  msg <- testthat::capture_messages(res <- RiboPilot:::export_project(p$dir))
  expect_equal(res, "g1")
  expect_equal(seen, "g1")
  expect_true(any(grepl("Skipped export group g2: the checks of its contigs are not finished", msg, fixed = TRUE)))
})

test_that("the report and its helpers read the copy", {
  skip_if_not_installed("Biostrings")
  p <- mk_review_export()
  con <- p$con
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 51L, 487L)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, itsx_warnings = ''")
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE element = '28S'")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::report_elements(con, units), "18S")
  expect_equal(as.character(RiboPilot:::gather_feature_seqs(con, units, "18S")[[1]]),
    substr(p$seq, 401, 487))
  svg <- RiboPilot:::feature_map_svg(con, "S1")
  expect_equal(lengths(regmatches(svg, gregexpr("<polygon", svg, fixed = TRUE))), 1L)

  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  out <- file.path(p$dir, "report.html")
  suppressMessages(RiboPilot:::export_report(con, units, out, "g1", alignments = FALSE))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # written by the Rmd, not by the fallback report
  expect_match(html, "Per-sample feature maps", fixed = TRUE)
  expect_match(html, "<td[^>]*>\\s*437\\s*</td>")
  expect_match(html, ">420</text>", fixed = TRUE)
  expect_false(grepl(">600</text>", html, fixed = TRUE))
})

test_that("the Export flow asks about an off-target verdict the copy's screen found", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1, export_group = 'g1', reviewed = 'no'
    WHERE path = 1 AND scaffold = 1")
  fixture_mark_done(con)
  DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 1, 'on_target')")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_screen SET screen_verdict = 'off_target', lca_taxon = 'Actiniaria'")
  DBI::dbDisconnect(con)
  u <- data.frame(ID = "S1", path = 1L, scaffold = 1L)
  fasta <- file.path(dir, "out", "export", "g1", "g1.fasta")

  shiny::testServer(RiboPilot:::app_server, {
    expect_equal(nrow(screen_flagged(u)), 0L)
    expect_equal(nrow(screen_flagged(u, review = TRUE)), 1L)
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    session$flushReact(); session$elapse(600); session$flushReact()
    session$setInputs(export_go = 1)
    session$flushReact(); session$elapse(600); session$flushReact()
    expect_false(file.exists(fasta))
    session$setInputs(export_screen_skip = 1)
    session$flushReact(); session$elapse(600); session$flushReact()
    expect_false(file.exists(fasta))
    session$setInputs(export_go = 2)
    session$flushReact(); session$elapse(600); session$flushReact()
    expect_warning(session$setInputs(export_screen_keep = 1), "off_target(Actiniaria)", fixed = TRUE)
    session$flushReact(); session$elapse(600); session$flushReact()
    expect_true(file.exists(fasta))
  })
})
```

Create `tests/testthat/test-export-views.R`:

```r
# Export reads per-contig data only through the review views (S9). A static
# check over the SQL strings of the export code: the original per-contig tables
# and the check columns of annotate are named only where the plan allows it.

export_view_root <- function() {
  root <- testthat::test_path("../..")
  if (!file.exists(file.path(root, "R", "export.R"))) skip("package sources not available")
  root
}

EXPORT_VIEW_FILES <- c("R/export.R", "R/export_opts.R", "R/export_report.R")
EXPORT_ORIG_TABLES <- c("assemblies", "annotations", "screen", "chimera", "chimera_windows",
  "coverage", "blast_hits", "blast_screen", "blast_ref_alignment")
EXPORT_CHECK_COLS <- c("n_features", "missing", "chimeric", "n_flags", "curate_flags",
  "itsx_warnings", "itsx_origin", "itsx_scores")
# Top-level object -> original tables its SQL may name. unit_fasta_text is Copy
# FASTA on the Assemble step; the others read assemblies only for ignore; the
# review = FALSE branches of unit_checks and export_summary_df serve the
# Assemble step.
EXPORT_VIEW_ALLOW <- list(
  unit_fasta_text = "assemblies", UNIT_MULTI_SQL = "assemblies",
  exportable_units = "assemblies", export_left_out = "assemblies",
  unit_checks = c(EXPORT_ORIG_TABLES, "annotate"),
  export_summary_df = c(EXPORT_ORIG_TABLES, "annotate"))

# Every string constant of a file, with the top-level object it sits in
export_strings <- function(f) {
  exprs <- parse(f, keep.source = TRUE)
  pd <- utils::getParseData(exprs)
  out <- lapply(seq_along(exprs), function(i) {
    e <- exprs[[i]]
    nm <- if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=")) as.character(e[[2]]) else ""
    sr <- attr(exprs, "srcref")[[i]]
    s <- pd[pd$token == "STR_CONST" & pd$line1 >= sr[1] & pd$line2 <= sr[3], ]
    if (nrow(s) == 0) return(NULL)
    data.frame(fn = nm, file = basename(f),
      txt = vapply(s$text, function(x) eval(parse(text = x)), character(1)),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

sql_tables <- function(txt) {
  m <- regmatches(txt, gregexpr("\\b(FROM|JOIN)\\s+[A-Za-z_]+", txt, perl = TRUE))[[1]]
  unique(sub("^(FROM|JOIN)\\s+", "", m))
}

test_that("export code names the original per-contig tables only where allowed", {
  root <- export_view_root()
  s <- do.call(rbind, lapply(file.path(root, EXPORT_VIEW_FILES), export_strings))
  expect_gt(nrow(s), 50)
  bad <- character(0)
  for (i in seq_len(nrow(s))) {
    tb <- sql_tables(s$txt[i])
    allow <- EXPORT_VIEW_ALLOW[[s$fn[i]]] %||% character(0)
    orig <- setdiff(intersect(tb, EXPORT_ORIG_TABLES), allow)
    if (length(orig)) bad <- c(bad, sprintf("%s %s reads %s", s$file[i], s$fn[i], paste(orig, collapse = ",")))
    if (s$fn[i] %in% c("unit_checks", "export_summary_df")) next
    cols <- paste(EXPORT_CHECK_COLS, collapse = "|")
    if (grepl(sprintf("\\ban\\.(%s)\\b", cols), s$txt[i], perl = TRUE) ||
        ("annotate" %in% tb && grepl(sprintf("(?<![.\\w])(%s)\\b", cols), s$txt[i], perl = TRUE))) {
      bad <- c(bad, sprintf("%s %s reads a check column of annotate", s$file[i], s$fn[i]))
    }
  }
  expect_equal(bad, character(0))
  # the functions allowed to read originals also read the views
  for (fn in c("unit_checks", "export_summary_df")) {
    txt <- paste(s$txt[s$fn == fn], collapse = "\n")
    for (v in c("review_units_v", "review_annotations_v", "review_screen_v", "review_chimera_v")) {
      expect_match(txt, v, fixed = TRUE, info = paste(fn, v))
    }
  }
})

test_that("every export summary in the export code reads the review layer", {
  root <- export_view_root()
  n <- 0
  for (f in file.path(root, EXPORT_VIEW_FILES)) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    fc <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "export_summary_df"]
    for (id in pd$parent[match(fc, pd$id)]) {
      expect_match(utils::getParseText(pd, id), "review = TRUE", fixed = TRUE, info = basename(f))
      n <- n + 1
    }
  }
  expect_gte(n, 2)
  rmd <- readLines(file.path(root, "inst", "export_report.Rmd"), warn = FALSE)
  calls <- grep("export_summary_df(", rmd, fixed = TRUE, value = TRUE)
  expect_length(calls, 1)
  expect_match(calls, "review = TRUE", fixed = TRUE)
  expect_false(any(grepl("dbGetQuery|dbExecute|dbReadTable|dbSendQuery", rmd)))
})
```

- [ ] **Step 2: Run them to see them fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-review.R")'`
Expected: FAIL, about 26 failures and 1 error in 5 of the 6 tests (the untrimmed test passes): the trimmed copy exports the original (sequence of 600 bases, defline `complete sequence`, `.tbl` lines `200\t>1` and `<600\t401`, summary length 600, `itsx_check` NA, `submission_risk` `misassembly` without `chimeric(ITSx)`, header length 600); `ex$done` is NULL and `export_left_out()` returns 0 rows; the not-done contig is exported; `export_project()` calls `export_bundle` for g2 too; `report_elements()` still lists `28S` and the report shows `>600</text>`; in the app test `screen_flagged(u, review = TRUE)` errors with `unused argument (review = TRUE)`.

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-views.R")'`
Expected: FAIL, 8 failures: `bad` lists `export.R feature_table_text_con reads annotations`, `... reads assemblies`, `export.R feature_map_svg reads annotations`, `export.R export_bundle reads assemblies`, `export.R export_bundle reads annotations`, `export_opts.R export_header_data reads assemblies`, `export_opts.R export_header_data reads a check column of annotate`, `export_opts.R export_header_data reads annotations`, `export_report.R report_elements reads annotations`, `export_report.R gather_feature_seqs reads annotations`, `export_report.R gather_feature_seqs reads assemblies`; export_summary_df names none of the four views; the two `export_summary_df(con, units)` calls in R/export.R and the Rmd call lack `review = TRUE`.

- [ ] **Step 3: Switch R/export.R to the review layer**

In `R/export.R`, `feature_table_text_con()`:

Find:
```r
  d <- DBI::dbGetQuery(con,
    "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  hdr <- paste0(">Feature ", seqid)
```

Replace with:
```r
  d <- DBI::dbGetQuery(con, paste("SELECT element, type, start, stop, strand,", REVIEW_EFF_PARTIALS,
    "FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start"),
    params = list(id, path, scaffold))
  hdr <- paste0(">Feature ", seqid)
```

Find:
```r
  st <- DBI::dbGetQuery(con, "SELECT sequence, length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
```

Replace with:
```r
  st <- DBI::dbGetQuery(con, "SELECT sequence, length FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
```

`export_summary_df()` (as Task 18 left it):

Find:
```r
#' stock ITSx check (annotate.itsx_warnings), `not_run(not checked)` when none
#' is stored.
#' @noRd
export_summary_df <- function(con, units) {
```

Replace with:
```r
#' stock ITSx check (annotate.itsx_warnings), `not_run(not checked)` when none
#' is stored.
#' @param review TRUE reads the per-contig values from the review layer (export,
#'   Annotation Review); FALSE from the originals (Assemble step). The marks
#'   come from annotate and the variant facts from the original screen in both.
#' @noRd
export_summary_df <- function(con, units, review = FALSE) {
```

Find:
```r
  m <- fetch(
    "SELECT s.ID, an.path, an.scaffold, asm.length, an.n_features,
       an.missing, an.reviewed, an.partial, an.problematic
     FROM samples s
     JOIN annotate an ON an.ID = s.ID
     JOIN assemblies asm ON asm.ID = s.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE s.ID IN (%s)")
  chi <- fetch("SELECT ID, path, scaffold, chimeric FROM annotate WHERE ID IN (%s)")
  chm <- fetch("SELECT ID, path, scaffold, label, reason FROM chimera WHERE ID IN (%s)")
  ann <- fetch("SELECT ID, path, scaffold, element, start, stop FROM annotations
                WHERE ID IN (%s) ORDER BY ID, path, scaffold, start")
  asm <- fetch("SELECT ID, path, scaffold, sequence FROM assemblies WHERE ID IN (%s)")
  scr <- fetch("SELECT ID, path, scaffold, variant_class, variant_group, variant_identity
                FROM screen WHERE ID IN (%s)")
  sv <- fetch("SELECT ID, path, scaffold, screen_verdict, lca_taxon FROM screen WHERE ID IN (%s)")
  nts <- fetch("SELECT ID, path, scaffold, curate_flags FROM annotate WHERE ID IN (%s)")
  itx <- fetch("SELECT ID, path, scaffold, itsx_warnings FROM annotate WHERE ID IN (%s)")
```

Replace with:
```r
  tb <- if (review) {
    c(res = "review_units_v", seq = "review_units_v", ann = "review_annotations_v",
      scr = "review_screen_v", chm = "review_chimera_v")
  } else {
    c(res = "annotate", seq = "assemblies", ann = "annotations", scr = "screen", chm = "chimera")
  }
  m <- fetch(if (review) {
    "SELECT s.ID, an.path, an.scaffold, v.length, v.n_features,
       v.missing, an.reviewed, an.partial, an.problematic
     FROM samples s
     JOIN annotate an ON an.ID = s.ID
     JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
     WHERE s.ID IN (%s)"
  } else {
    "SELECT s.ID, an.path, an.scaffold, asm.length, an.n_features,
       an.missing, an.reviewed, an.partial, an.problematic
     FROM samples s
     JOIN annotate an ON an.ID = s.ID
     JOIN assemblies asm ON asm.ID = s.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE s.ID IN (%s)"
  })
  chi <- fetch(paste("SELECT ID, path, scaffold, chimeric FROM", tb[["res"]], "WHERE ID IN (%s)"))
  chm <- fetch(paste("SELECT ID, path, scaffold, label, reason FROM", tb[["chm"]], "WHERE ID IN (%s)"))
  ann <- fetch(paste("SELECT ID, path, scaffold, element, start, stop FROM", tb[["ann"]],
    "WHERE ID IN (%s) ORDER BY ID, path, scaffold, start"))
  asm <- fetch(paste("SELECT ID, path, scaffold, sequence FROM", tb[["seq"]], "WHERE ID IN (%s)"))
  scr <- fetch("SELECT ID, path, scaffold, variant_class, variant_group, variant_identity
                FROM screen WHERE ID IN (%s)")
  sv <- fetch(paste("SELECT ID, path, scaffold, screen_verdict, lca_taxon FROM", tb[["scr"]], "WHERE ID IN (%s)"))
  nts <- fetch(paste("SELECT ID, path, scaffold, curate_flags FROM", tb[["res"]], "WHERE ID IN (%s)"))
  itx <- fetch(paste("SELECT ID, path, scaffold, itsx_warnings FROM", tb[["res"]], "WHERE ID IN (%s)"))
```

`exportable_units()`:

Find:
```r
#' The units an export writes, with their export group
#' @noRd
exportable_units <- function(con) {
  DBI::dbGetQuery(con, paste(
    "SELECT an.ID, an.path, an.scaffold, an.export_group FROM annotate an
     JOIN assemble a ON a.ID = an.ID
     JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE", EXPORTABLE_WHERE, "ORDER BY an.ID, an.path, an.scaffold"))
}
```

Replace with:
```r
#' The units an export may write, with their export group
#'
#' `done` is 1 when the contig's checks are done (review_units_v.check_switch
#' = 2). Only done units are written; every row counts toward the sibling count
#' that names records, so record names do not change when checks finish.
#' @noRd
exportable_units <- function(con) {
  DBI::dbGetQuery(con, paste(
    "SELECT an.ID, an.path, an.scaffold, an.export_group,
       CASE WHEN v.check_switch = 2 THEN 1 ELSE 0 END AS done
     FROM annotate an
     JOIN assemble a ON a.ID = an.ID
     JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     LEFT JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
     WHERE", EXPORTABLE_WHERE, "ORDER BY an.ID, an.path, an.scaffold"))
}
```

`export_left_out()`:

Find:
```r
            WHEN COALESCE(an.annotate_lock,0) != 1 THEN 'annotation not locked'
            ELSE 'contig dropped' END AS why
     FROM annotate an
     LEFT JOIN assemble a ON a.ID = an.ID
     LEFT JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE an.export_group = ? AND NOT (", EXPORTABLE_WHERE, ")
     ORDER BY an.ID, an.path, an.scaffold"), params = list(group))
```

Replace with:
```r
            WHEN COALESCE(an.annotate_lock,0) != 1 THEN 'annotation not locked'
            WHEN COALESCE(asm.ignore,0) != 0 THEN 'contig dropped'
            ELSE 'checks not finished' END AS why
     FROM annotate an
     LEFT JOIN assemble a ON a.ID = an.ID
     LEFT JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     LEFT JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
     WHERE an.export_group = ? AND (NOT (", EXPORTABLE_WHERE, ")
       OR (asm.ID IS NOT NULL AND v.check_switch IS NOT 2))
     ORDER BY an.ID, an.path, an.scaffold"), params = list(group))
```

`feature_map_svg()`:

Find:
```r
  d <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  if (nrow(d) == 0) {
    return("<em>no annotations</em>")
```

Replace with:
```r
  d <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  if (nrow(d) == 0) {
    return("<em>no annotations</em>")
```

`export_report_html()`:

Find:
```r
export_report_html <- function(con, units, out_html, label, smry = export_summary_df(con, units)) {
```

Replace with:
```r
export_report_html <- function(con, units, out_html, label, smry = export_summary_df(con, units, review = TRUE)) {
```

`export_bundle()` roxygen (Task 18's text):

Find:
```r
#' @return (invisibly) the group output directory, with the summary data.frame
#'   attached as attribute `summary`.
#' @details A unit is exported when its sample's assembly and the unit's
#'   annotation are both locked, its contig is kept, and it is not excluded
#'   from export. The `itsx_check` column of the summary CSV reports the stock
```

Replace with:
```r
#' @return (invisibly) the group output directory, with the summary data.frame
#'   attached as attribute `summary` and the contigs left out because their
#'   checks are not finished as attribute `not_done` (data.frame(ID, path,
#'   scaffold)).
#' @details A unit is exported when its sample's assembly and the unit's
#'   annotation are both locked, its contig is kept, it is not excluded from
#'   export, and its checks are done (press Update on Annotation Review after
#'   an edit). Every value written comes from the reviewed copy of the contig
#'   when it has one. The `itsx_check` column of the summary CSV reports the stock
```

`export_bundle()` body:

Find:
```r
    sel <- sel[key(sel) %in% key(units), , drop = FALSE]
  }
  if (nrow(sel) == 0) {
```

Replace with:
```r
    sel <- sel[key(sel) %in% key(units), , drop = FALSE]
  }
  not_done <- sel[sel$done != 1, c("ID", "path", "scaffold"), drop = FALSE]
  rownames(not_done) <- NULL
  sel <- sel[sel$done == 1, , drop = FALSE]
  if (nrow(sel) == 0 && nrow(not_done) > 0) {
    stop("No contigs to export: the checks of ", nrow(not_done), " contig(s) are not finished. ",
      "Press Update on Annotation Review, then export again.", call. = FALSE)
  }
  if (nrow(sel) == 0) {
```

Find:
```r
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, upath, uscaf))$sequence[1]
```

Replace with:
```r
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, upath, uscaf))$sequence[1]
```

Find:
```r
      ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
        params = list(id, upath, uscaf))
```

Replace with:
```r
      ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
        params = list(id, upath, uscaf))
```

Find:
```r
  smry <- export_summary_df(con, units)
  risky <- which(nzchar(smry$submission_risk))
```

Replace with:
```r
  smry <- export_summary_df(con, units, review = TRUE)
  risky <- which(nzchar(smry$submission_risk))
```

Find:
```r
  message("Exported ", nrow(units), " contig(s) to ", gdir,
    "\n  GenBank Portal takes FASTA only: upload ", label,
```

Replace with:
```r
  message("Exported ", nrow(units), " contig(s) to ", gdir,
    if (nrow(not_done) > 0) paste0("\n  ", nrow(not_done), " contig(s) left out: checks not finished"),
    "\n  GenBank Portal takes FASTA only: upload ", label,
```

Find:
```r
  attr(gdir, "summary") <- smry
  invisible(gdir)
```

Replace with:
```r
  attr(gdir, "summary") <- smry
  attr(gdir, "not_done") <- not_done
  invisible(gdir)
```

`export_project()`:

Find:
```r
#' both its assembly and its annotation locked into one bundle named "export".
#' @param path Project directory.
#' @param out_dir Project "out" directory; bundles land in `<out_dir>/export`.
#' @return (invisibly) the bundle names written.
#' @export
export_project <- function(path = ".", out_dir = file.path(path, "out")) {
  con <- open_project_db(file.path(path, ".sqlite"))
  g <- tryCatch(exportable_units(con)$export_group, finally = DBI::dbDisconnect(con))
  groups <- sort(unique(g[nzchar(g)]))
  if (length(groups) == 0) {
    export_bundle(path, out_dir = out_dir, zip = TRUE)
    return(invisible("export"))
  }
  for (g in groups) export_bundle(path, group = g, out_dir = out_dir, zip = TRUE)
  invisible(groups)
}
```

Replace with:
```r
#' both its assembly and its annotation locked into one bundle named "export".
#' Only contigs whose checks are done are written; a group with none is
#' skipped with a message.
#' @param path Project directory.
#' @param out_dir Project "out" directory; bundles land in `<out_dir>/export`.
#' @return (invisibly) the bundle names written.
#' @export
export_project <- function(path = ".", out_dir = file.path(path, "out")) {
  con <- open_project_db(file.path(path, ".sqlite"))
  ex <- tryCatch(exportable_units(con), finally = DBI::dbDisconnect(con))
  g <- ex$export_group
  groups <- sort(unique(g[nzchar(g)]))
  if (length(groups) == 0) {
    export_bundle(path, out_dir = out_dir, zip = TRUE)
    return(invisible("export"))
  }
  done <- groups[groups %in% ex$export_group[ex$done == 1]]
  for (g in setdiff(groups, done)) {
    message("Skipped export group ", g, ": the checks of its contigs are not finished. ",
      "Press Update on Annotation Review, then export again.")
  }
  for (g in done) export_bundle(path, group = g, out_dir = out_dir, zip = TRUE)
  invisible(done)
}
```

`unit_fasta_text()` (Copy FASTA on the Assemble step) keeps reading `assemblies`; `unit_checks()` is Task 14's and already reads the views in review mode.

- [ ] **Step 4: Switch the header data, the report helpers, and the report template**

In `R/export_opts.R`, `export_header_data()` (sequence was fetched but never used):

Find:
```r
    asm <- DBI::dbGetQuery(con,
      "SELECT length, topology, sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    ann <- DBI::dbGetQuery(con,
      "SELECT n_features, missing, partial, reviewed, export_group
         FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    feat <- DBI::dbGetQuery(con,
      "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(id, p, s))
```

Replace with:
```r
    asm <- DBI::dbGetQuery(con,
      "SELECT length, topology, n_features, missing FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    ann <- DBI::dbGetQuery(con,
      "SELECT partial, reviewed, export_group
         FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    feat <- DBI::dbGetQuery(con, paste("SELECT element, start, stop, strand,", REVIEW_EFF_PARTIALS,
      "FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start"),
      params = list(id, p, s))
```

Find:
```r
      n_features = if (nrow(ann)) ann$n_features[1] else NA_integer_,
      missing = if (nrow(ann)) ann$missing[1] else NA_character_,
```

Replace with:
```r
      n_features = if (nrow(asm)) asm$n_features[1] else NA_integer_,
      missing = if (nrow(asm)) asm$missing[1] else NA_character_,
```

In `R/export_report.R`, `report_elements()`:

Find:
```r
      "SELECT DISTINCT element FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
```

Replace with:
```r
      "SELECT DISTINCT element FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ?",
```

`gather_feature_seqs()`:

Find:
```r
      "SELECT start, stop FROM annotations
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ? ORDER BY start",
```

Replace with:
```r
      "SELECT start, stop FROM review_annotations_v
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ? ORDER BY start",
```

Find:
```r
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))$sequence
    if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) next
```

Replace with:
```r
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))$sequence
    if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) next
```

In `inst/export_report.Rmd` (its feature maps, elements, and alignments already go through `feature_map_svg`, `report_elements`, and `gather_feature_seqs`):

Find:
```r
smry <- rp$export_summary_df(con, units)
```

Replace with:
```r
smry <- rp$export_summary_df(con, units, review = TRUE)
```

- [ ] **Step 5: The Export flow's contamination prompt reads the review layer**

In `R/app_server.R` (Task 21 passes `review` from `lock_flagged`; leave that call as it is):

Find:
```r
  screen_then_proceed <- function(u, grp) {
    flagged <- screen_flagged(u)
```

Replace with:
```r
  screen_then_proceed <- function(u, grp) {
    flagged <- screen_flagged(u, review = TRUE)
```

Find:
```r
  screen_flagged <- function(u) {
    if (NROW(u) == 0) return(data.frame())
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    s <- export_summary_df(con, cbind(u, seqid = u$ID))
```

Replace with:
```r
  screen_flagged <- function(u, review = FALSE) {
    if (NROW(u) == 0) return(data.frame())
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    s <- export_summary_df(con, cbind(u, seqid = u$ID), review = review)
```

Find:
```r
    u <- export_pending()
    f <- screen_flagged(u)
```

Replace with:
```r
    u <- export_pending()
    f <- screen_flagged(u, review = TRUE)
```

Check: `Rscript -e 'for (f in c("R/export.R", "R/export_opts.R", "R/export_report.R", "R/app_server.R")) invisible(parse(f))'` prints nothing.

- [ ] **Step 6: Give the hand-built export fixtures the review layer**

The views read columns and tables these fixtures lack (a broken view fails every SELECT from it), so each fixture runs `fixture_review_schema(con)` after its positional INSERTs; fixtures that export also run `fixture_mark_done(con)`, since export now leaves out contigs whose checks are not done. Fixtures with no `annotate` table get an `annotate` row so `review_units_v` (whose fallback joins `annotate`) returns the contig's sequence and length.

In `tests/testthat/test-export.R` (feature table test, :12-29):

Find:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',201,260)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "<1\t100\trRNA")            # 5' partial
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',201,260)")
  fixture_review_schema(con)
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 1)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "<1\t100\trRNA")            # 5' partial
```

(partial ends test, :31-55)

Find:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',300,460,'+',0,1)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',300,460,'+',0,1)")
  fixture_review_schema(con)
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 1)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
```

(unit marked partial, :87-104)

Find:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',100,260,0,0)")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n100\t260\trRNA")
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',100,260,0,0)")
  fixture_review_schema(con)
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n100\t260\trRNA")
```

(partial end exported partial, :106-130)

Find:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',155,1488,'length,truncated',0,1)")
  expect_match
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',155,1488,'length,truncated',0,1)")
  fixture_review_schema(con)
  expect_match
```

(feature map with no annotations, :258-263)

Find:
```r
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT, stop INT)")
  expect_equal(RiboPilot:::feature_map_svg(con, "S1"), "<em>no annotations</em>")
```

Replace with:
```r
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT, stop INT)")
  fixture_review_schema(con)
  expect_equal(RiboPilot:::feature_map_svg(con, "S1"), "<em>no annotations</em>")
```

The summary test at :132-176 keeps reading originals (`export_summary_df()` defaults to `review = FALSE`) and needs no change.

In `tests/testthat/test-export-units.R`, `mk_project()` (the `ALTER TABLE annotate ADD COLUMN export_time_stamp` above stays, since it runs first):

Find:
```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL)")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1),('S2',1)")
  dir
}
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL)")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1),('S2',1)")
  fixture_review_schema(con)
  fixture_mark_done(con)
  dir
}
```

The exact-CSV test: a done contig carries the stored ITSx value (`''` from `fixture_mark_done`), no longer the not-checked text:

Find:
```r
      paste0("\"S2\",\"S2\",1,1,\"Genus speciesB\",400,2,\"\",\"yes\",\"no\",\"no\",NA,NA,NA,NA,\"\",\"\",",
        "\"not_run(not checked)\"")))
```

Replace with:
```r
      paste0("\"S2\",\"S2\",1,1,\"Genus speciesB\",400,2,\"\",\"yes\",\"no\",\"no\",NA,NA,NA,NA,\"\",\"\",",
        "\"\"")))
```

In `tests/testthat/test-export-headers.R`, `mk_hdr_project()`:

Find:
```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL,'Reef B')")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1),('S2',1)")
  dir
}
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL,'Reef B')")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1),('S2',1)")
  fixture_review_schema(con)
  fixture_mark_done(con)
  dir
}
```

In `tests/testthat/test-export-opts.R`, `mk_opts_db()`:

Find:
```r
      params = list(f[1], as.integer(f[2]), as.integer(f[3])))
  }
  con
}
```

Replace with:
```r
      params = list(f[1], as.integer(f[2]), as.integer(f[3])))
  }
  fixture_review_schema(con)
  con
}
```

In `tests/testthat/test-export-report.R`, `mk_report_project()`:

Find:
```r
      params = list(id, paste("Genus species", k)))
  }
  dir
}
```

Replace with:
```r
      params = list(id, paste("Genus species", k)))
  }
  fixture_review_schema(con)
  dir
}
```

In `tests/testthat/test-nmnh.R` ("export writes NMNH values and drops the empty modifier"):

Find:
```r
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1)")
  DBI::dbDisconnect(con)
  local_mocked_bindings(.gbif_get = function(path) stop("offline"))
```

Replace with:
```r
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1)")
  fixture_review_schema(con)
  fixture_mark_done(con)
  DBI::dbDisconnect(con)
  local_mocked_bindings(.gbif_get = function(path) stop("offline"))
```

No change is needed in `tests/testthat/test-app-export.R` (`mk_export_project()` already stores annotate_switch 2 and itsx_warnings `''`, and runs `fixture_fill_schema()`), `tests/testthat/test-itsx-check.R` (Task 18 removed its export through `mk_itsx_project`; its summary test marks its contig done), or `tests/testthat/test-utils.R` (export_bundle still opens the database once; a `new_db()` project with nothing exportable still stops with "No contigs to export").

- [ ] **Step 7: Regenerate the help pages of export_bundle and export_project**

Run: `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`

Then `git status --short man/ DESCRIPTION` must list only `man/export_bundle.Rd` and `man/export_project.Rd`; restore any other file roxygen rewrote (for example a link target that differs between R versions, or a `RoxygenNote` line added to DESCRIPTION) with `git checkout -- <file>`. Expected diff: in `export_bundle.Rd` the `\value` gains "and the contigs left out because their checks are not finished as attribute \code{not_done} (data.frame(ID, path, scaffold))" and `\details` gains "and its checks are done (press Update on Annotation Review after an edit). Every value written comes from the reviewed copy of the contig when it has one."; in `export_project.Rd` the description gains "Only contigs whose checks are done are written; a group with none is skipped with a message." Never run `devtools::document()` (NAMESPACE is hand-maintained).

- [ ] **Step 8: Run the tests to see them pass**

Run each:
```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-review.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-views.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-units.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-headers.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-opts.R")'
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-report.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nmnh.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-itsx-check.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-utils.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-export.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-specimen-app.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'
```
Expected: FAIL 0 in every file. test-export-review.R PASS 41 (with pandoc; without it the report part skips), test-export-views.R PASS 16, test-export.R PASS 92, test-export-headers.R PASS 60, test-export-opts.R PASS 43, test-export-report.R PASS 22, test-nmnh.R PASS 116, test-utils.R PASS 31; the other files keep the counts the earlier tasks left.

- [ ] **Step 9: Commit**

```bash
git add R/export.R R/export_opts.R R/export_report.R inst/export_report.Rmd R/app_server.R man/export_bundle.Rd man/export_project.Rd tests/testthat/test-export-review.R tests/testthat/test-export-views.R tests/testthat/test-export.R tests/testthat/test-export-units.R tests/testthat/test-export-headers.R tests/testthat/test-export-opts.R tests/testthat/test-export-report.R tests/testthat/test-nmnh.R
git commit -m "Export reads the review layer; leaves out unfinished checks"
```


### Task 20: Export step table

**Files:**
- Modify: `R/export.R` (new constant `EXPORT_STEP_SQL` after `UNIT_EXPORTED_SQL`)
- Modify: `R/app_server.R` (`all_rows_base` Export branch; `output$tbl` Export colDefs; `export_units_df`; `with_marker_fields`; `run_export`)
- Test: `tests/testthat/test-export-views.R`, `tests/testthat/test-app-export.R`

**Interfaces:**
- Consumes: `exportable_units(con)` with the `done` column (T19); `export_left_out(con, group)` with why `'checks not finished'` (T19); `export_check_summary_ui(smry, not_done = 0L)` (T18); `rt_gated(cell = NULL)` and the hidden `check_state` colDef (T14); `review_copy(con, id, path, scaffold)` (T02); the views `review_units_v`, `review_screen_v`, `review_chimera_v`, `review_annotations_v` (T01); `export_strings(f)`, `sql_tables(txt)`, `EXPORT_VIEW_ALLOW`, `EXPORT_ORIG_TABLES` in test-export-views.R (T19).
- Produces: `EXPORT_STEP_SQL` (chr, R/export.R, @noRd): the Export step table query. Columns: ID, path, scaffold, path_n, scaffold_n, Taxon, sci_name, taxid, length, n_features, missing, n_flags, curate_flags, annotate_notes, itsx_warnings, screen_verdict, chimera_label, chimera_reason, reviewed, partial, export_group, exported, check_state, time_stamp. `exported` and `check_state` are `'checks not finished, left out'` when `COALESCE(review_units_v.check_switch, 1) != 2`; else `exported` is today's text and `check_state` NULL. In app_server, `export_units_df(group = NULL)` returns only `done == 1` rows.

The Export step table reads the review layer (S8 Export step): length and checks from `review_units_v`, screen from `review_screen_v`, chimera from `review_chimera_v`, the marks from `annotate`. Its query moves into R/export.R as `EXPORT_STEP_SQL`, so the static check in test-export-views.R covers it. Contigs whose checks are not done stay listed but show "checks not finished, left out" in the exported column (an Export-only expression; `UNIT_EXPORTED_SQL` stays as it is for Annotation Review) and their check cells show the same info pill through `rt_gated()`, which hides check values left from a discarded recheck. `export_units_df()` keeps only done contigs, so the selected count, the group count, the Export button, the header previews, the NMNH lookup, and the contamination prompt (S6.2) all leave them out, while `app_export_seqids()` keeps counting every exportable row so record names never change (S6.2). The marker header preview reads its sample element from `review_annotations_v`, and the Export done dialog counts the group's contigs left out because their checks are not finished (S9).

- [ ] **Step 1: Write the failing tests**

In `tests/testthat/test-export-views.R`, `getParseData()` stores a long string constant as a placeholder like `[1043 chars quoted with '"']`, which `parse()` cannot read; `EXPORT_STEP_SQL` holds such a string. Read the text with `getParseText()` instead.

Find:

```r
      txt = vapply(s$text, function(x) eval(parse(text = x)), character(1)),
```

Replace with:

```r
      txt = vapply(s$id, function(i) eval(parse(text = utils::getParseText(pd, i))), character(1)),
```

`EXPORT_STEP_SQL` joins `assemblies` for the `ignore` gate of `EXPORTABLE_WHERE`, like `exportable_units`.

Find:

```r
  unit_fasta_text = "assemblies", UNIT_MULTI_SQL = "assemblies",
```

Replace with:

```r
  unit_fasta_text = "assemblies", UNIT_MULTI_SQL = "assemblies", EXPORT_STEP_SQL = "assemblies",
```

Append to the end of `tests/testthat/test-export-views.R`:

```r

test_that("the Export step table and its marker fields read the review layer", {
  root <- export_view_root()
  s <- export_strings(file.path(root, "R", "export.R"))
  txt <- paste(s$txt[s$fn == "EXPORT_STEP_SQL"], collapse = "\n")
  for (v in c("review_units_v", "review_screen_v", "review_chimera_v")) {
    expect_match(txt, v, fixed = TRUE, info = v)
  }
  f <- file.path(root, "R", "app_server.R")
  expect_match(paste(readLines(f, warn = FALSE), collapse = "\n"), "app_query(EXPORT_STEP_SQL)", fixed = TRUE)
  pd <- utils::getParseData(parse(f, keep.source = TRUE))
  sym <- pd$id[pd$token == "SYMBOL" & pd$text == "with_marker_fields"]
  defs <- vapply(sym, function(i) utils::getParseText(pd, pd$parent[pd$id == pd$parent[pd$id == i]]), "")
  def <- grep("^with_marker_fields <- function", defs, value = TRUE)
  expect_length(def, 1)
  expect_match(def, "FROM review_annotations_v", fixed = TRUE)
  expect_equal(intersect(sql_tables(def), EXPORT_ORIG_TABLES), character(0))
})
```

Append to the end of `tests/testthat/test-app-export.R`:

```r

# S1 p1.s1 done, with a review copy (length 500, 4 flags); S1 p1.s2 locked in
# g1 but not checked by ITSx yet, so its checks are not done.
mk_export_not_done <- function(env = parent.frame()) {
  dir <- mk_export_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 2, 'linear', 400, ?, 0)", params = list(strrep("ACGT", 100)))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
    reviewed, problematic, partial, n_features, n_flags, itsx_warnings, time_stamp)
    VALUES ('S1', 1, 2, 2, 1, 'no', 'no', 'no', 0, 0, NULL, 1)")
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_units SET length = 500, n_flags = 4
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  dir
}

test_that("the Export step reads the review layer and leaves out contigs whose checks are not done", {
  skip_if_not_installed("withr")
  dir <- mk_export_not_done()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "DELETE FROM annotations")
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE element != '18S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET element = '28S'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export", export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = TRUE, export_report = FALSE)
    settle(session)
    d <- samples_data()
    d <- d[order(d$scaffold), ]
    expect_equal(d$scaffold, c(1L, 2L))
    expect_equal(d$length, c(500L, 400L))
    expect_equal(d$n_flags, c(4L, 0L))
    expect_equal(d$exported, c("", "checks not finished, left out"))
    expect_equal(d$check_state, c(NA, "checks not finished, left out"))
    j <- jsonlite::fromJSON(as.character(output$tbl), simplifyVector = FALSE)
    cd <- Filter(function(c) identical(c$id, "exported"), j$x$tag$attribs$columns)[[1]]
    expect_match(paste(unlist(cd$cell), collapse = " "), "check_state", fixed = TRUE)

    u <- export_units_df("g1")
    expect_equal(u$scaffold, 1L)
    expect_equal(nrow(export_units_df()), 1L)
    expect_equal(output$export_sel_n, "1 contig locked & ready")
    con <- RiboPilot:::app_con()
    expect_equal(app_export_seqids(con, u), "S1_p1_s1")
    DBI::dbDisconnect(con)
    expect_equal(export_hdr_data()$seqid, "S1_p1_s1")
    expect_equal(export_hdr_data_marker()$element, "28S")
    expect_match(output$export_summary$html, "1 contig in group", fixed = TRUE)
    expect_match(output$export_summary$html, "S1 p1.s2 (checks not finished)", fixed = TRUE)
  })
})

test_that("Export writes only the done contig and its done dialog counts the one left out", {
  skip_if_not_installed("withr")
  dir <- mk_export_not_done()
  withr::local_dir(dir)
  shown <- character(0)
  local_mocked_bindings(showModal = function(ui, ...) shown <<- c(shown, as.character(ui)))

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export", export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    session$setInputs(export_go = 1)
    settle(session)
  })
  fa <- readLines(file.path(dir, "out", "export", "g1", "g1.fasta"))
  expect_equal(grep("^>", fa, value = TRUE), grep("^>S1_p1_s1", fa, value = TRUE))
  expect_length(grep("^>", fa), 1)
  done <- grep("Export complete", shown, value = TRUE)
  expect_length(done, 1)
  expect_match(done, "1 contig left out: checks not finished.", fixed = TRUE)
})
```

- [ ] **Step 2: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-views.R")'`
Expected: FAIL 6 in "the Export step table and its marker fields read the review layer" (no `EXPORT_STEP_SQL`, so its three view matches fail; `app_query(EXPORT_STEP_SQL)` not found; `with_marker_fields` reads `FROM annotations`). The two older tests in the file pass.

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-export.R")'`
Expected: FAIL in the two new tests only: `d$length` is `600 400` and `d$n_flags` `0 0` (originals), `d$exported` is `"" ""`, `d$check_state` is NULL, `export_units_df("g1")` has two rows, `export_hdr_data_marker()$element` is `"18S"`, and the done dialog lacks "1 contig left out: checks not finished.". Every older test passes.

- [ ] **Step 3: Add EXPORT_STEP_SQL to R/export.R**

Find:

```r
UNIT_EXPORTED_SQL <- "CASE WHEN an.export_time_stamp IS NULL THEN ''
   ELSE 'exported ' || date(an.export_time_stamp, 'unixepoch', 'localtime') END AS exported"
```

Replace with:

```r
UNIT_EXPORTED_SQL <- "CASE WHEN an.export_time_stamp IS NULL THEN ''
   ELSE 'exported ' || date(an.export_time_stamp, 'unixepoch', 'localtime') END AS exported"

#' The Export step table: the units an export may write (EXPORTABLE_WHERE), one
#' row per (ID, path, scaffold), with contig data and checks from the review
#' layer and the marks from annotate. A contig whose checks are not done shows
#' "checks not finished, left out" in exported and in check_state, which gates
#' the check cells (rt_gated).
#' @noRd
EXPORT_STEP_SQL <- paste(
  "SELECT s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ", s.Taxon, s.sci_name, s.taxid, v.length,
     v.n_features, v.missing,
     v.n_flags, v.curate_flags, an.annotate_notes, v.itsx_warnings, sc.screen_verdict,
     ch.label AS chimera_label, ch.reason AS chimera_reason,
     an.reviewed, an.partial, an.export_group,
     CASE WHEN COALESCE(v.check_switch, 1) != 2 THEN 'checks not finished, left out'
          WHEN an.export_time_stamp IS NULL THEN ''
          ELSE 'exported ' || date(an.export_time_stamp, 'unixepoch', 'localtime') END AS exported,
     CASE WHEN COALESCE(v.check_switch, 1) != 2 THEN 'checks not finished, left out' END AS check_state,
     an.time_stamp
   FROM annotate an JOIN samples s ON s.ID = an.ID
     JOIN assemble a ON a.ID = an.ID
     JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     LEFT JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
     LEFT JOIN review_screen_v sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
     LEFT JOIN review_chimera_v ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
   WHERE", EXPORTABLE_WHERE)
```

- [ ] **Step 4: Point the Export step in R/app_server.R at the review layer**

The table query (in `all_rows_base`).

Find:

```r
    } else if (mode == "Export") {
      # the units an export writes (EXPORTABLE_WHERE), one row per (ID, path, scaffold)
      app_query(paste(
        "SELECT s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ", s.Taxon, s.sci_name, s.taxid, asm.length,
           an.n_features, an.missing,
           an.n_flags, an.curate_flags, an.annotate_notes, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason,
           an.reviewed, an.partial,
           an.export_group,", UNIT_EXPORTED_SQL, ", an.time_stamp
         FROM annotate an JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE", EXPORTABLE_WHERE))
```

Replace with:

```r
    } else if (mode == "Export") {
      app_query(EXPORT_STEP_SQL)
```

The exported column on the Export step shows the left-out pill (Annotation Review keeps the shared green pill).

Find:

```r
    if (mode == "Export") {
      cols[["ID"]] <- cd("ID", sticky = "left", width = mp_fit_width(d$ID))
    }
```

Replace with:

```r
    if (mode == "Export") {
      cols[["ID"]] <- cd("ID", sticky = "left", width = mp_fit_width(d$ID))
      cols[["exported"]] <- cd("exported", width = 210, html = TRUE, cell = rt_gated(rt_text_pill("success")))
    }
```

Only done contigs are written, counted, or previewed. `app_export_seqids()` stays as it is: it reads `exportable_units(con)` directly and counts every row.

Find:

```r
  # Exportable units (exportable_units()), narrowed to one export group when given.
  export_units_df <- function(group = NULL) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    locked <- exportable_units(con)
```

Replace with:

```r
  # Exportable units whose checks are done, narrowed to one export group when given.
  export_units_df <- function(group = NULL) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    locked <- exportable_units(con)
    locked <- locked[locked$done == 1, , drop = FALSE]
```

Find:

```r
    el <- app_query("SELECT element FROM annotations LIMIT 1")$element
```

Replace with:

```r
    el <- app_query("SELECT element FROM review_annotations_v LIMIT 1")$element
```

The done dialog counts the group's contigs left out because their checks are not finished. `run_export` gets only done contigs, so `attr(gdir, "not_done")` is empty here; the count comes from `export_left_out()`.

Find:

```r
    if (is.null(gdir)) return(invisible())
    export_done_dir(gdir)
```

Replace with:

```r
    if (is.null(gdir)) return(invisible())
    export_done_dir(gdir)
    con <- app_con()
    out <- tryCatch(export_left_out(con, grp), finally = DBI::dbDisconnect(con))
```

Find:

```r
      export_check_summary_ui(attr(gdir, "summary")),
```

Replace with:

```r
      export_check_summary_ui(attr(gdir, "summary"), not_done = sum(out$why == "checks not finished")),
```

- [ ] **Step 5: Run the tests to see them pass**

Run:

```bash
for f in export-views app-export export-review specimen-app nmnh chimera app-table; do
  Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-$f.R')"
done
```

Expected: PASS, FAIL 0 in every file (in validation: export-views 23, app-export 160, export-review 37 with 1 skip, specimen-app 99, nmnh 116, chimera 115 with 2 skips, app-table 47).

- [ ] **Step 6: Commit**

```bash
git add R/export.R R/app_server.R tests/testthat/test-export-views.R tests/testthat/test-app-export.R
git commit -m "Export step table reads the review layer; leaves out unfinished contigs"
```


### Task 21: Lock gate

**Files:**
- Modify: `R/review.R` (append `review_lock_gate()` at the end of the file)
- Modify: `R/app_server.R` (`lock_flagged` repo :765-768, the lock check in `lock_rows` repo :817, `lock_unit` repo :3133-3134)
- Modify: `tests/testthat/test-app-scaffold-uid.R` (two lock tests, repo :475-528 and :530-570)
- Test: `tests/testthat/test-lock-gate.R` (new)

**Interfaces:**
- Consumes: `review_unit_state(con, units)` (Task 14; columns include `has_copy` and `state` in "waiting", "not run", "running", "failed", "on hold", "done"); `unit_key(id, path, scaffold)` (R/export.R); `screen_flagged(u, review = FALSE)` (Task 19); `review_copy()` (Task 2); `fixture_mark_done(con)`, `mk_app_project()` (Task 1); `chimera_test_con()` (helper-chimera.R).
- Produces: `review_lock_gate(con, units)` in R/review.R: character vector, one element per row of `units` (data.frame ID, path, scaffold), NA where Lock may go ahead (done, or `annotate.problematic = 'yes'`), else the next-step message of SPEC S6.2. In R/app_server.R: `lock_flagged(u, review = FALSE)` and the local guard `need_checks_done(u)` (TRUE when every unit passes the gate; otherwise a warning alert titled "Checks not finished" and FALSE).

Lock on Annotation Review (toolbar `annot_lock_btn`, the row lock icon `lock_row`, and Lock / Lock & Close in the details window) now refuses contigs whose checks are not done, naming the next step per state (SPEC S6.2, D2). One refused contig refuses the whole click, as the one-path rule does; nothing is locked. Exclude & Lock (`problematic_btn`, `mark_problematic`) stays ungated, and a plain Lock on an excluded contig passes too: an excluded contig is never exported (contract plan decision). The flagged-contig question before a lock reads the review layer on Annotation Review and the originals on the Assemble step (SPEC S8, "The lock check reuses the export risk path"). The Assemble and Taxonomy locks are unchanged.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-lock-gate.R`:

```r
# Lock on Annotation Review waits for finished checks (S6.2); Exclude & Lock does not.

lock_gate_units <- function(con) {
  seq <- paste(rep("ACGT", 50), collapse = "")
  DBI::dbExecute(con, "DELETE FROM annotate")
  # scaffold: annotate_switch, n_features, itsx_warnings, problematic
  rows <- list(c(1, 2, 3, "", "no"), c(2, 1, 3, NA, "no"), c(3, 2, 3, NA, "no"), c(4, 2, 3, "", "no"),
    c(5, 2, 3, "", "no"), c(6, 2, 3, "", "no"), c(7, 2, 3, "", "no"), c(8, 2, 3, "", "no"),
    c(9, 1, 3, NA, "yes"), c(10, 2, NA, NA, "no"))
  for (r in rows) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', 200, ?, 0)", params = list(as.integer(r[1]), seq))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features,
      itsx_warnings, problematic) VALUES ('S1', 1, ?, ?, 0, ?, ?, ?)",
      params = list(as.integer(r[1]), as.integer(r[2]), as.integer(r[3]), r[4], r[5]))
  }
  for (s in 4:8) RiboPilot:::review_copy(con, "S1", 1L, s)
  for (x in list(c(4, 1), c(5, 4), c(6, 3), c(7, 0), c(8, 2))) {
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = ? WHERE ID = 'S1' AND path = 1 AND scaffold = ?",
      params = list(as.integer(x[2]), as.integer(x[1])))
  }
}

test_that("review_lock_gate names the next step for each state and passes done and excluded contigs", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  lock_gate_units(con)
  u <- data.frame(ID = "S1", path = 1L, scaffold = c(8L, 1:7, 9:10))
  expect_identical(RiboPilot:::review_lock_gate(con, u), c(
    NA,                                                  # copy, done
    NA,                                                  # no copy, done
    "Press Update on the Assemble step first.",          # no copy, queued for the main pipeline
    "Set State to Ready to run, then press Update.",     # no copy, ITSx never checked
    "Press Update first.",                               # copy, not run
    "Wait for the running Update to finish.",            # copy, running
    "Set State to Ready to run, then press Update.",     # copy, failed
    "Set State to Ready to run, then press Update.",     # copy, on hold
    NA,                                                  # excluded, though queued
    "Press Update on the Assemble step first."))         # no copy, never annotated
  expect_identical(RiboPilot:::review_lock_gate(con, u[0, ]), character(0))
})

test_that("Lock on Annotation Review refuses contigs whose checks are not done", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  said <- character(0)
  local_mocked_bindings(mp_alert = function(title, text = NULL, ...) said <<- c(said, paste(title, text)))
  lk <- function() RiboPilot:::app_query(
    "SELECT COALESCE(annotate_lock, 0) AS v FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")$v
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- function(p, s) which(d$path == p & d$scaffold == s)[1]

    # every contig is waiting for the main pipeline: toolbar, row icon, and details window refuse
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(1, 2)), annot_lock_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_match(said[length(said)], "Checks not finished", fixed = TRUE)
    expect_match(said[length(said)], "S1 p1.s1, S1 p1.s2: Press Update on the Assemble step first.", fixed = TRUE)
    session$setInputs(lock_row = as.character(row(1, 1)))
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    n <- length(said)
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(lock_and_close = 1)
    session$flushReact()
    session$setInputs(modal_lock = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_equal(length(said), n + 2L)

    # one done contig with one not done: nothing is locked
    set("UPDATE annotate SET annotate_switch = 2, itsx_warnings = '' WHERE path = 1 AND scaffold = 1")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(1, 2)), annot_lock_btn = 2)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_match(said[length(said)], "S1 p1.s2: Press Update on the Assemble step first.", fixed = TRUE)
    expect_false(grepl("p1.s1", said[length(said)], fixed = TRUE))

    # the done contig alone locks, from the toolbar and from the details window
    session$setInputs(tbl__reactable__selected = row(1, 1), annot_lock_btn = 3)
    session$flushReact()
    expect_equal(lk(), c(1L, 0L, 0L, 0L))
    set("UPDATE annotate SET annotate_lock = 0")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(modal_lock = 2)
    session$flushReact()
    expect_equal(lk(), c(1L, 0L, 0L, 0L))

    # Exclude & Lock is exempt, and so is a plain Lock on an excluded contig
    set("UPDATE annotate SET annotate_lock = 0")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(tbl__reactable__selected = row(1, 2), problematic_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 1L, 0L, 0L))
    set("UPDATE annotate SET problematic = 'yes' WHERE path = 2 AND scaffold = 1")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(lock_row = as.character(row(2, 1)))
    session$flushReact()
    expect_equal(lk(), c(0L, 1L, 1L, 0L))
  })
})

test_that("Lock on Annotation Review asks about flags found on the review copy", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  fixture_mark_done(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, itsx_warnings = '', curate_flags = 'contaminant'
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  local_mocked_bindings(mp_alert = function(...) NULL)
  lk <- function() RiboPilot:::app_query(
    "SELECT COALESCE(annotate_lock, 0) AS v FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")$v
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- function(p, s) which(d$path == p & d$scaffold == s)[1]
    session$setInputs(tbl__reactable__selected = row(1, 1), annot_lock_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    session$setInputs(lock_confirm = TRUE)
    expect_equal(lk(), c(1L, 0L, 0L, 0L))

    # the details window asks too
    set("UPDATE annotate SET annotate_lock = 0")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(modal_lock = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
  })
})
```

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lock-gate.R")'`
Expected: FAIL (about 12 failures and 1 error in all 3 tests): the first test errors with `object 'review_lock_gate' not found`; in the second, the waiting contigs get locked (`lk()` is not all 0) and no "Checks not finished" alert is recorded; in the third, Lock locks at once because the flag check reads the original `annotate.curate_flags` (NULL), not the copy's `contaminant`.

- [ ] **Step 3: Add review_lock_gate()**

Append to the end of `R/review.R`:

```r
#' Why each unit cannot be locked on Annotation Review (S6.2)
#'
#' NA where Lock may go ahead: the checks are done, or the contig is excluded
#' (an excluded contig is never exported). Otherwise the next step to take.
#' @noRd
review_lock_gate <- function(con, units) {
  if (!NROW(units)) return(character(0))
  st <- review_unit_state(con, units)
  ex <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate WHERE problematic = 'yes'")
  again <- "Set State to Ready to run, then press Update."
  msg <- unname(c(waiting = "Press Update on the Assemble step first.", `not run` = again,
    running = "Wait for the running Update to finish.", failed = again, `on hold` = again,
    done = NA_character_)[st$state])
  msg[st$state == "not run" & st$has_copy %in% 1] <- "Press Update first."
  msg[unit_key(units$ID, units$path, units$scaffold) %in% unit_key(ex$ID, ex$path, ex$scaffold)] <- NA_character_
  msg
}
```

- [ ] **Step 4: Gate the Annotation Review locks**

In `R/app_server.R`, Find:

```r
  lock_flagged <- function(u) {
    f <- screen_flagged(u)
    if (nrow(f)) f[!(f$problematic %in% "yes"), , drop = FALSE] else f
  }
```

Replace with:

```r
  lock_flagged <- function(u, review = FALSE) {
    f <- screen_flagged(u, review = review)
    if (nrow(f)) f[!(f$problematic %in% "yes"), , drop = FALSE] else f
  }
  # Annotation Review locks only contigs whose checks are done (S6.2).
  need_checks_done <- function(u) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    g <- review_lock_gate(con, u)
    bad <- !is.na(g)
    if (!any(bad)) return(TRUE)
    by <- split(sprintf("%s p%s.s%s", u$ID, u$path, u$scaffold)[bad], g[bad])
    rp_alert(paste0("Lock needs finished checks, so nothing was locked. ",
      paste0(vapply(by, mp_id_list, ""), ": ", names(by), collapse = " ")),
      type = "warning", title = "Checks not finished")
    FALSE
  }
```

Find (in `lock_rows`):

```r
    f <- lock_flagged(units)
    if (nrow(f)) return(confirm_flagged_lock(f, do_lock))
```

Replace with:

```r
    ann <- identical(mm$pfx, "annotate")
    if (ann && !need_checks_done(units)) return()
    f <- lock_flagged(units, review = ann)
    if (nrow(f)) return(confirm_flagged_lock(f, do_lock))
```

(`units` there is the unlocked part of the selection; the gate runs only on a lock, never on an unlock, which returns earlier.)

Find (in `lock_unit`):

```r
    if (!need_one_path(as.data.frame(u))) return()
    f <- lock_flagged(as.data.frame(u))
```

Replace with:

```r
    if (!need_one_path(as.data.frame(u)) || !need_checks_done(as.data.frame(u))) return()
    f <- lock_flagged(as.data.frame(u), review = TRUE)
```

Leave `problematic_btn` and `mark_problematic` as they are (Exclude & Lock stays ungated).

- [ ] **Step 5: Update the two existing Annotation Review lock tests**

`mk_app_project()` seeds every contig as waiting for the main pipeline (`annotate_switch = 1`, `itsx_warnings` NULL), so the gate now refuses the locks these tests expect; mark the contigs done first.

In `tests/testthat/test-app-scaffold-uid.R`, Find:

```r
test_that("a second assembly path is refused from the toolbar and from Lock & Close", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
```

Replace with:

```r
test_that("a second assembly path is refused from the toolbar and from Lock & Close", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  fixture_mark_done(con)
  DBI::dbDisconnect(con)
```

Find:

```r
  DBI::dbExecute(con, "UPDATE annotate SET problematic = 'yes' WHERE path = 1 AND scaffold = 2")
```

Replace with:

```r
  fixture_mark_done(con)
  DBI::dbExecute(con, "UPDATE annotate SET problematic = 'yes' WHERE path = 1 AND scaffold = 2")
```

The lock tests in `tests/testthat/test-app-export.R` use `mk_export_project()`, whose contig is done (`annotate_switch = 2`, `itsx_warnings = ''`), and `tests/testthat/test-app-table.R` locks only on the Assemble step; neither needs a change.

- [ ] **Step 6: Run the tests to see them pass**

Run:

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lock-gate.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-scaffold-uid.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-export.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-app-table.R")'
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-export-review.R")'
```

Expected: PASS, 0 failures in each (test-lock-gate.R 18 expectations; test-app-scaffold-uid.R 31 tests). test-app-export.R may show one error in "export_project skips a group with no unit that has both locks" only where no `zip` binary is installed; it is unrelated to this task.

- [ ] **Step 7: Commit**

```bash
git add R/review.R R/app_server.R tests/testthat/test-lock-gate.R tests/testthat/test-app-scaffold-uid.R
git commit -m "Lock on Annotation Review waits for finished checks"
```


### Task 22: Upgrade step 4 conversion; assembly_backup removed

**Files:**
- Create: `R/review_upgrade.R`
- Create: `tests/testthat/test-review-upgrade.R`
- Modify: `R/chimera.R` (`chimera_shift()` write block, repo :575; `chimera_restore()` roxygen title and write block, repo :600-605)
- Modify: `R/db_migrate.R` (the step 4 `fn` Task 1 appended to `db_migrations()`)
- Modify: `R/db_schema.R` (delete `ASSEMBLY_BACKUP_DDL`, repo :145-158)
- Modify: `R/init_db.R` (`new_db()`, repo :481-482)
- Modify: `R/pipeline_steps.R` (`UNIT_EVIDENCE_TABLES`, repo :239-241 as Task 2 left it; `load_assemblies()`, repo :547-549)
- Modify: `tests/testthat/helper-app-fixture.R` (`fixture_fill_schema()`, :5-7)
- Modify: `tests/testthat/test-load-assemblies.R` (`mk_asm_db()`, :20)
- Modify: `tests/testthat/test-backwards-compat.R` (parity test, as Task 1 left it)
- Modify: `tests/testthat/test-trim-coverage.R` (as Task 15 rewrote it: header comment; one test appended)
- Test: `tests/testthat/test-review-upgrade.R`, `tests/testthat/test-trim-coverage.R`, `tests/testthat/test-backwards-compat.R`

Line numbers are those of the repo before this plan; earlier tasks shift some of them, so find each block by its text.

**Interfaces:**
- Consumes:
  - Task 1: `db_migrations()[[4]]` (note `"review copies for edited contigs; ITSx origin and scores"`; fn = the two `ALTER TABLE annotate ADD COLUMN` statements, then `review_apply_schema(con)`), the 8 review tables, and `review_units_v` (columns used here: `ID, path, scaffold, has_copy, orig_length, check_switch`; fallback `check_switch` = `annotate_switch`, except 2 with NULL `itsx_warnings` reads 1). `tests/testthat/fixtures/schema_v1.sql` holds `assembly_backup`; the parity test in test-backwards-compat.R ends with `expect_setequal(names(bc_views(old)), names(RiboPilot:::REVIEW_VIEWS_DDL))`.
  - Task 2 (R/review.R): `with_savepoint(con, code, name = "rp_review")` (BEGIN IMMEDIATE outside a transaction, SAVEPOINT inside one; an error rolls back that level and propagates), `review_copy(con, id, path, scaffold)` (copies the unit's current originals; `trim_start = 1`, `trim_end = length`, `edited = 0`, check state as the view showed it), `review_counts(con, id, path, scaffold)` -> `list(n_features, missing)`, `review_invalidate(con, id, path, scaffold)` (check_switch 1, rev + 1, check results cleared, the six check mirrors emptied, `status`, `flags`, `partial_start`, `partial_stop` of the copy's features cleared). `UNIT_EVIDENCE_TABLES` ends with `"assembly_backup", REVIEW_TABLES)`.
  - Task 6: `curate_one(con, id, path, scaffold)` writes through `with_savepoint()`, so it works inside the upgrade transaction.
  - Task 15: `chimera_shift(con, id, path, scaffold, offset, len)` (new position = old - offset), `chimera_restore(con, id, path, scaffold, snap)`, `write_coverage_summary(con, id, path, scaffold, dep)`, all in their repo form apart from `chimera_restore()`'s roxygen title; test-trim-coverage.R with `mk_trim_db()` and `chimera_json()`. Task 15 deleted `trim_assembly()`, `restore_assembly()`, `ensure_assembly_backup()`, and `chimera_snapshot()` and moved `unit_trimmed()` to the review layer; Task 17 removed `chimera_original()`. After them only the lines this task deletes still name `assembly_backup`.
  - Existing: `parse_depths()` (R/coverage.R), `record_schema_event()`, `project_schema_version()`.
- Produces:
  - `db_migrations()[[4]]`'s fn ends with `review_upgrade(con)` and `DBI::dbExecute(con, "DROP TABLE assembly_backup")`. A version 4 project upgraded to 5 has review copies for its old trims and edits and no `assembly_backup`; `new_db()` no longer creates the table.
  - `review_upgrade(con)` (R/review_upgrade.R, `@noRd`): the S11.3-11.4 conversion. Each contig converts in its own `with_savepoint(name = "rp_upgrade")`; a contig that errors is rolled back and listed, never stopping the step. One `message()` (only when there is something to say) lists: the count converted, `Converted with problems: <ID.path.scaffold> (<problem>); ...`, `Not converted, left as they were: <ID.path.scaffold> (<error>); ...`, `Locked contigs marked done without a recheck: <ID.path.scaffold>, ...` (S11.4), and `Trimmed contigs left as they are (dropped or never annotated): ...`.
  - Internal helpers in the same file: `review_upgrade_unit(con, u)` (one contig; returns its problems), `review_upgrade_original(con, u, off)` (writes a trimmed original back from its backup row), `review_upgrade_done(con, id, path, scaffold)` (check_switch 2, counts from `review_counts()`, `COALESCE(itsx_warnings, 'not_run(not checked before upgrade)')`).
  - `chimera_shift()` and `chimera_restore()` write through `with_savepoint()`: they work alone and inside an open transaction; signatures unchanged.
  - Removed: `ASSEMBLY_BACKUP_DDL`, the `assembly_backup` entry of `UNIT_EVIDENCE_TABLES`, and `load_assemblies()`'s backup delete.

Older projects hold trims and edits in the original tables, with the pre-trim state in `assembly_backup` (S2). Schema step 4 now converts them (S11.3): every contig whose sequence differs from its backup, or that has a feature with `annotations.edited = 1`, gets a review copy of its current state with `edited = 1`, the window found by locating the current sequence in the backup sequence, and the partial ends of edited features moved to the user columns; a trimmed contig's original is then written back from the backup (sequence, length, features, depth track, chimera snapshot or a reverse shift), curated again, and its `itsx_warnings` cleared. Unlocked copies are queued (S5.2 invalidate); locked copies are marked done, with `not_run(not checked before upgrade)` for a NULL ITSx value, and so is every locked contig without a copy that does not read as done (S11.4). Because a failed upgrade locks the user out of the project (Review Focus 1), each contig converts in its own savepoint and is listed instead of failing the step, units are selected by joining `assemblies` and `annotate` (orphan backup rows are only listed), a backup without a sequence or one that does not contain the current sequence keeps the copy at full length without a write-back, and a chimera snapshot that cannot be put back falls back to the reverse shift. The step then drops `assembly_backup` (S11.5, S4.1), and the table leaves `new_db()`, the fixtures, and the sweeps. `chimera_shift()` and `chimera_restore()` move to `with_savepoint()` first, because the conversion calls them inside the upgrade transaction.

- [ ] **Step 1: Write the failing test for the chimera helpers inside a transaction**

In `tests/testthat/test-trim-coverage.R`, the header comment names the table this task drops:

Find:

```r
# triage read) and the chimera rows, moved by a shift or put back from the JSON
# an older trim stored in assembly_backup.chimera.
```

Replace with:

```r
# triage read) and the chimera rows, moved by a shift or put back from the JSON
# an older trim saved.
```

Append to the end of `tests/testthat/test-trim-coverage.R`:

```r

test_that("chimera_shift and chimera_restore work inside an open transaction", {
  dir <- withr::local_tempdir()
  con <- mk_trim_db(dir, strrep("ACGT", 25), rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1',1,1,'severe','30,70')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  win <- function() DBI::dbGetQuery(con, "SELECT start FROM chimera_windows ORDER BY start")$start
  snap <- chimera_json(con)

  DBI::dbWithTransaction(con, {
    RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 40L, 60L)
    expect_equal(win(), 11L)
    RiboPilot:::chimera_restore(con, "S1", 1L, 1L, snap)
  })
  expect_equal(win(), c(1L, 51L))

  # a snapshot that no longer fits the table leaves the rows as they were
  bad <- sub('"label":', '"retired":1,"label":', snap, fixed = TRUE)
  DBI::dbWithTransaction(con, {
    expect_error(RiboPilot:::chimera_restore(con, "S1", 1L, 1L, bad))
    expect_equal(win(), c(1L, 51L))
    RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 40L, 60L)
  })
  expect_equal(win(), 11L)
  expect_false(RSQLite::sqliteIsTransacting(con))
})
```

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-trim-coverage.R")'`
Expected: FAIL, `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 27 ]`: "chimera_shift and chimera_restore work inside an open transaction" errors with `cannot start a transaction within a transaction` (chimera_shift's `DBI::dbWithTransaction()` inside the test's transaction).

- [ ] **Step 3: Write the chimera rows through with_savepoint() (R/chimera.R)**

Find (in `chimera_shift()`):

```r
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    if (nrow(w)) DBI::dbAppendTable(con, "chimera_windows", w)
```

Replace with:

```r
  with_savepoint(con, {
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    if (nrow(w)) DBI::dbAppendTable(con, "chimera_windows", w)
```

Find:

```r
#' Put back a chimera result and windows stored as JSON (assembly_backup.chimera)
#' @noRd
chimera_restore <- function(con, id, path, scaffold, snap) {
  s <- jsonlite::fromJSON(snap)
  p <- list(id, as.integer(path), as.integer(scaffold))
  DBI::dbWithTransaction(con, {
```

Replace with:

```r
#' Put back a chimera result and windows saved as JSON by an older trim
#' @noRd
chimera_restore <- function(con, id, path, scaffold, snap) {
  s <- jsonlite::fromJSON(snap)
  p <- list(id, as.integer(path), as.integer(scaffold))
  with_savepoint(con, {
```

Both blocks keep their reads before the write and hold only `DBI::dbExecute()` and `DBI::dbAppendTable()` calls, so nothing inside opens a transaction. Check: `Rscript -e 'invisible(parse("R/chimera.R"))'`.

- [ ] **Step 4: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("trim-coverage", "chimera")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-trim-coverage.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 33 ]`; test-chimera.R `FAIL 0` (its BLAST tests skip without blastn).

- [ ] **Step 5: Write the failing upgrade tests**

Create `tests/testthat/test-review-upgrade.R`:

```r
# Schema step 4 turns the trims and edits an older RiboPilot wrote into the
# original tables into review copies, then drops assembly_backup (SPEC S11.3-11.5).

RU_SEQ <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
RU_OTHER <- withr::with_seed(8, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
RU_DEP <- c(rep(5, 100), rep(15, 100))

# A version 4 database: the 0.1.0 schema plus steps 1-3, with sample S1
ru_v4 <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  for (st in RiboPilot:::db_migrations()[1:3]) st$fn(con)
  RiboPilot:::record_schema_event(con, "created", 4L)
  DBI::dbExecute(con, "INSERT INTO samples (ID) VALUES ('S1')")
  con
}

# Contig S1.1.<s>: 200 bp, three + strand features, a depth track, and a
# chimera result over three windows
ru_unit <- function(con, s, lock = 0L, sw = 2L, itsx = "itsx_x") {
  x <- function(sql, ...) DBI::dbExecute(con, sql, params = list(...))
  x("INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore, edited)
     VALUES ('S1', 1, ?, 'linear', 200, ?, 0, 0)", s, RU_SEQ)
  x("INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features, missing,
       n_flags, curate_flags, itsx_warnings) VALUES ('S1', 1, ?, ?, ?, 99, 'x', 7, 'stale', ?)",
    s, sw, lock, itsx)
  for (f in list(list("18S", 21L, 80L), list("ITS1", 81L, 120L), list("28S", 121L, 190L))) {
    x("INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length, strand,
         model_coverage, tool, status, edited)
       VALUES ('S1', 1, ?, ?, 'rRNA', ?, ?, ?, '+', 99, 'barrnap', 'stale', 0)",
      s, f[[1]], f[[2]], f[[3]], f[[3]] - f[[2]] + 1L)
  }
  x("INSERT INTO coverage (ID, path, scaffold, mean_depth, median_depth, min_depth, max_depth, depths)
     VALUES ('S1', 1, ?, 10, 10, 5, 15, ?)", s, paste(RU_DEP, collapse = ","))
  x("INSERT INTO chimera (ID, path, scaffold, score, label, junctions)
     VALUES ('S1', 1, ?, 0.9, 'severe', '30,150')", s)
  for (w in c(1L, 51L, 101L)) {
    x("INSERT INTO chimera_windows (ID, path, scaffold, start, stop, lineage, lineage_taxids, pident)
       VALUES ('S1', 1, ?, ?, ?, 'K;P;C1;o;f;g;s', '1;2;3;;;;', 99)", s, w, w + 49L)
  }
}

# What the old trim_assembly() left after cutting S1.1.<s> to 41..170: the
# backup row (snapshots as JSON), the cut sequence, features shifted with cut
# ends partial and edited = 1, the depth track cut, and the chimera rows shifted
ru_trim <- function(con, s, snap = TRUE) {
  rows <- function(t) DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s WHERE scaffold = ?",
    paste(DBI::dbListFields(con, t), collapse = ", "), t), params = list(s))
  chim <- if (snap) {
    as.character(jsonlite::toJSON(list(chimera = rows("chimera"), windows = rows("chimera_windows")),
      dataframe = "rows", na = "null", digits = NA))
  } else NA_character_
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length, annotations,
      coverage, chimera, time_stamp) VALUES ('S1', 1, ?, ?, 200, ?, ?, ?, 0)",
    params = list(s, RU_SEQ, as.character(jsonlite::toJSON(rows("annotations"), dataframe = "rows",
      na = "null")), paste(RU_DEP, collapse = ","), chim))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 130, edited = 1 WHERE scaffold = ?",
    params = list(substr(RU_SEQ, 41, 170), s))
  DBI::dbExecute(con, "UPDATE annotations SET start = MAX(start - 40, 1), stop = MIN(stop - 40, 130),
      edited = 1, partial_start = CASE WHEN start < 41 THEN 1 ELSE partial_start END,
      partial_stop = CASE WHEN stop > 170 THEN 1 ELSE partial_stop END
    WHERE scaffold = ?", params = list(s))
  DBI::dbExecute(con, "UPDATE annotations SET length = stop - start + 1 WHERE scaffold = ?", params = list(s))
  RiboPilot:::write_coverage_summary(con, "S1", 1L, s, RU_DEP[41:170])
  RiboPilot:::chimera_shift(con, "S1", 1L, s, 40L, 130L)
}

test_that("step 4 turns old trims and edits into review copies and lists what it could not convert", {
  dir <- withr::local_tempdir()
  con <- ru_v4(dir)
  for (s in 1:12) {
    ru_unit(con, s, lock = as.integer(s %in% c(4, 5, 12)), sw = if (s == 4) 1L else 2L,
      itsx = if (s %in% c(4, 5)) NA_character_ else if (s == 12) "" else "itsx_x")
  }
  # 1 trimmed; 2 trimmed, then restored; 3 nudged, its 5' end set by hand;
  # 4 locked, State set by hand; 5 trimmed and locked, no ITSx value, no chimera
  # snapshot; 6 its sequence is not in the backup; 7 backup without a sequence;
  # 8 snapshots holding a column the tables no longer have; 9 its features
  # cannot be written back; 10 trimmed, then lost its annotate row; 11 fails
  # part way; 12 locked and done
  for (s in c(1, 5:11)) ru_trim(con, s, snap = s != 5)
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length)
    VALUES ('S1', 1, 2, ?, 200)", params = list(RU_SEQ))
  DBI::dbExecute(con, "UPDATE annotations SET start = 25, length = 56, edited = 1, partial_start = 1
    WHERE scaffold = 3 AND element = '18S'")
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = ? WHERE scaffold = 6", params = list(RU_OTHER))
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = NULL WHERE scaffold = 7")
  DBI::dbExecute(con, "UPDATE assembly_backup SET annotations = replace(annotations, ?, ?),
      chimera = replace(chimera, ?, ?) WHERE scaffold = 8",
    params = list('"element":', '"retired":1,"element":', '"label":', '"retired":1,"label":'))
  DBI::dbExecute(con, "CREATE TRIGGER ru_boom BEFORE INSERT ON annotations WHEN NEW.scaffold = 9
    BEGIN SELECT RAISE(ABORT, 'ru_boom'); END")
  DBI::dbExecute(con, "DELETE FROM annotate WHERE scaffold = 10")
  DBI::dbDisconnect(con)

  real <- RiboPilot:::review_invalidate
  local_mocked_bindings(review_invalidate = function(con, id, path, scaffold) {
    if (scaffold == 11) stop("disk full")
    real(con, id, path, scaffold)
  })
  msgs <- testthat::capture_messages(backwards_compatibility(dir))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::project_schema_version(file.path(dir, ".sqlite")), 5L)
  expect_false("assembly_backup" %in% DBI::dbListTables(con))

  # the copies hold the trimmed or edited state; unlocked ones wait for Update
  nr <- "not_run(not checked before upgrade)"
  expect_equal(DBI::dbGetQuery(con, "SELECT scaffold, length, trim_start, trim_end, edited, check_switch,
      n_features, missing, n_flags, itsx_warnings FROM review_units ORDER BY scaffold"),
    data.frame(scaffold = c(1L, 3:9), length = c(130L, 200L, 200L, rep(130L, 5)),
      trim_start = c(41L, 1L, 1L, 41L, 1L, 1L, 41L, 1L), trim_end = c(170L, 200L, 200L, 170L, 130L, 130L, 170L, 130L),
      edited = c(1L, 1L, 0L, rep(1L, 5)), check_switch = c(1L, 1L, 2L, 2L, rep(1L, 4)),
      n_features = 3L, missing = "5.8S,ITS2", n_flags = c(NA, NA, 7L, 7L, rep(NA, 4)),
      itsx_warnings = c(NA, NA, nr, nr, rep(NA, 4))))
  expect_equal(DBI::dbGetQuery(con, "SELECT sequence FROM review_units WHERE scaffold IN (1, 5, 8)")$sequence,
    rep(substr(RU_SEQ, 41, 170), 3))
  rf <- function(s) DBI::dbGetQuery(con, "SELECT element, start, stop, partial_start, partial_stop,
      partial_start_user, partial_stop_user, status FROM review_annotations WHERE scaffold = ? ORDER BY start",
    params = list(s))
  feats <- data.frame(element = c("18S", "ITS1", "28S"), start = c(1L, 41L, 81L), stop = c(40L, 80L, 130L),
    partial_start = c(1L, NA, NA), partial_stop = c(NA, NA, 1L), partial_start_user = c(1L, NA, NA),
    partial_stop_user = c(NA, NA, 1L), status = "stale")
  # locked: the copy keeps its check results
  expect_equal(rf(5), feats)
  # unlocked: the cut ends stay as the user's choice, curation's values are cleared
  feats[c("partial_start", "partial_stop", "status")] <- list(NA_integer_, NA_integer_, NA_character_)
  expect_equal(rf(1), feats)
  expect_equal(rf(3)[c("element", "start", "partial_start_user", "partial_start")],
    data.frame(element = c("18S", "ITS1", "28S"), start = c(25L, 81L, 121L),
      partial_start_user = c(1L, NA, NA), partial_start = NA_integer_))
  expect_equal(DBI::dbGetQuery(con, "SELECT scaffold, start FROM review_chimera_windows ORDER BY scaffold, start"),
    data.frame(scaffold = c(4L, 4L, 4L, 5L, 5L), start = c(1L, 51L, 101L, 11L, 61L)))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_annotations WHERE scaffold = 11")$n, 0L)

  # trimmed originals are put back and curated again; the others keep their state
  expect_equal(DBI::dbGetQuery(con, "SELECT length, sequence = ? AS full FROM assemblies ORDER BY scaffold",
      params = list(RU_SEQ)),
    data.frame(length = c(rep(200L, 5), 130L, 130L, 200L, 130L, 130L, 130L, 200L),
      full = c(rep(1L, 5), 0L, 0L, 1L, 0L, 0L, 0L, 1L)))
  expect_equal(DBI::dbGetQuery(con, "SELECT start FROM annotations WHERE element = '18S' ORDER BY scaffold")$start,
    c(21L, 21L, 25L, 21L, 21L, 1L, 1L, 21L, 1L, 1L, 1L, 21L))
  expect_equal(DBI::dbGetQuery(con, "SELECT element, start, stop, edited FROM annotations WHERE scaffold = 8
      ORDER BY start"),
    data.frame(element = c("18S", "ITS1", "28S"), start = c(21L, 81L, 121L), stop = c(80L, 120L, 190L),
      edited = 0L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations
      WHERE scaffold IN (1, 5, 8) AND status = 'stale'")$n, 0L)
  an <- DBI::dbGetQuery(con, "SELECT scaffold, itsx_warnings, curate_flags FROM annotate ORDER BY scaffold")
  back <- an$scaffold %in% c(1, 5, 8)
  expect_true(all(is.na(an$itsx_warnings[back])))
  expect_false(any(an$curate_flags[back] %in% "stale"))
  expect_equal(an$itsx_warnings[an$scaffold %in% c(3, 6, 7, 9, 11)], rep("itsx_x", 5))
  expect_equal(an$curate_flags[!back], rep("stale", 8))
  cv <- DBI::dbGetQuery(con, "SELECT median_depth, depths FROM coverage WHERE scaffold IN (1, 9) ORDER BY scaffold")
  expect_equal(cv$median_depth, c(10, 15))
  expect_equal(lengths(strsplit(cv$depths, ",")), c(200L, 130L))
  ch <- function(s) list(
    DBI::dbGetQuery(con, "SELECT label, junctions FROM chimera WHERE scaffold = ?", params = list(s)),
    DBI::dbGetQuery(con, "SELECT start FROM chimera_windows WHERE scaffold = ? ORDER BY start",
      params = list(s))$start)
  expect_equal(ch(1), list(data.frame(label = "severe", junctions = "30,150"), c(1L, 51L, 101L)))
  # no snapshot, or one that no longer fits: the rows are shifted back
  expect_equal(ch(5), list(data.frame(label = "severe", junctions = "150"), c(51L, 101L)))
  expect_equal(ch(8), ch(5))
  expect_equal(ch(9)[[2]], c(11L, 61L))

  v <- DBI::dbGetQuery(con, "SELECT scaffold, has_copy, orig_length, check_switch FROM review_units_v
    ORDER BY scaffold")
  expect_equal(v, data.frame(scaffold = c(1:9, 11:12), has_copy = c(1L, 0L, rep(1L, 7), 0L, 0L),
    orig_length = c(rep(200L, 5), 130L, 130L, 200L, 130L, 130L, 200L),
    check_switch = c(1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L)))

  m <- grep("S1.1.", msgs, fixed = TRUE, value = TRUE)
  expect_length(m, 1)
  expect_match(m, "Moved the edits of 7 contig(s)", fixed = TRUE)
  expect_match(m, "S1.1.6 (original not restored: the trimmed sequence is not in the saved original)",
    fixed = TRUE)
  expect_match(m, "S1.1.7 (original not restored: the saved original has no sequence)", fixed = TRUE)
  expect_match(m, "S1.1.8 (chimera result shifted back instead of restored: ", fixed = TRUE)
  expect_match(m, "S1.1.9 (original not restored: ru_boom", fixed = TRUE)
  expect_match(m, "Not converted, left as they were: S1.1.11 (disk full).", fixed = TRUE)
  expect_match(m, "Locked contigs marked done without a recheck: S1.1.4.", fixed = TRUE)
  expect_match(m, "Trimmed contigs left as they are (dropped or never annotated): S1.1.10.", fixed = TRUE)
  expect_no_match(m, "S1\\.1\\.(1|2|3|5|12)\\b", perl = TRUE)
})
```

The fixture rebuilds what the deleted `trim_assembly()` stored (git show a25fd25:R/pipeline_steps.R, lines 572-641): the annotations snapshot is every annotations row as JSON, the chimera snapshot `list(chimera, windows)` as JSON, `coverage` the depth CSV. Contig 8's snapshots carry a `retired` key the tables do not have, so the features still go back (by column name) while the chimera result falls back to the reverse shift; contig 9's trigger makes the write-back fail half way; the mock makes contig 11 fail after its write-back, so the whole contig must roll back. Contig 10 has no annotate row, so it is only listed.

In `tests/testthat/test-backwards-compat.R`, the parity test also checks that the table is gone (new_db() must match):

Find:

```r
  expect_setequal(names(bc_views(old)), names(RiboPilot:::REVIEW_VIEWS_DDL))
```

Replace with:

```r
  expect_setequal(names(bc_views(old)), names(RiboPilot:::REVIEW_VIEWS_DDL))
  expect_false("assembly_backup" %in% names(bc_schema(old)))
```

- [ ] **Step 6: Run them to see them fail**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-upgrade", "backwards-compat")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: FAIL. test-review-upgrade.R `[ FAIL 28 | WARN 0 | SKIP 0 | PASS 15 ]`, failing from its first check after the upgrade: `"assembly_backup" %in% DBI::dbListTables(con)` is TRUE, `review_units` is empty, the originals still hold the trimmed sequences (lengths 130 for contigs 1, 5, 8), and no message names a contig (`expect_length(m, 1)` gets 0). test-backwards-compat.R `[ FAIL 1 | WARN 0 | SKIP 0 | PASS 38 ]`: `"assembly_backup" %in% names(bc_schema(old))` is TRUE.

- [ ] **Step 7: Convert old edits in step 4 and remove assembly_backup**

Create `R/review_upgrade.R`:

```r
# Schema step 4 (version 4 to 5): the trims and edits an older RiboPilot wrote
# into the original tables become review copies, and trimmed originals are put
# back from assembly_backup (SPEC S11.3-11.4).

#' Convert old trims and edits into review copies (schema step 4)
#'
#' Runs inside the upgrade transaction. Each contig converts in its own
#' savepoint, so a contig that fails is rolled back and listed instead of
#' stopping the upgrade. One message lists the converted contigs, the ones with
#' a problem or left as they were, and the locked ones marked done.
#' @noRd
review_upgrade <- function(con) {
  uid <- function(d) paste(d$ID, d$path, d$scaffold, sep = ".")
  u <- DBI::dbGetQuery(con,
    "SELECT a.ID, a.path, a.scaffold, a.sequence AS cur, b.ID IS NOT NULL AS has_bak,
       b.sequence AS bak, b.annotations AS bak_ann, b.coverage AS bak_cov, b.chimera AS bak_chim,
       COALESCE(an.annotate_lock, 0) AS locked
     FROM assemblies a
     JOIN annotate an ON an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold
     LEFT JOIN assembly_backup b ON b.ID = a.ID AND b.path = a.path AND b.scaffold = a.scaffold
     WHERE (b.sequence IS NOT NULL AND a.sequence IS NOT b.sequence)
        OR EXISTS (SELECT 1 FROM annotations f WHERE f.ID = a.ID AND f.path = a.path
                     AND f.scaffold = a.scaffold AND f.edited = 1)
     ORDER BY a.ID, a.path, a.scaffold")
  n <- 0L
  bad <- skip <- done <- character()
  for (i in seq_len(nrow(u))) {
    r <- tryCatch(with_savepoint(con, review_upgrade_unit(con, u[i, ]), name = "rp_upgrade"),
      error = function(e) e)
    if (inherits(r, "error")) {
      skip <- c(skip, sprintf("%s (%s)", uid(u[i, ]), conditionMessage(r)))
    } else {
      n <- n + 1L
      if (length(r)) bad <- c(bad, sprintf("%s (%s)", uid(u[i, ]), paste(r, collapse = "; ")))
    }
  }

  lk <- DBI::dbGetQuery(con,
    "SELECT v.ID, v.path, v.scaffold FROM review_units_v v
     JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold
     WHERE v.has_copy = 0 AND an.annotate_lock = 1 AND v.check_switch IS NOT 2
     ORDER BY v.ID, v.path, v.scaffold")
  for (i in seq_len(nrow(lk))) {
    x <- lk[i, ]
    r <- tryCatch(with_savepoint(con, {
      review_copy(con, x$ID, x$path, x$scaffold)
      review_upgrade_done(con, x$ID, x$path, x$scaffold)
    }, name = "rp_upgrade"), error = function(e) e)
    if (inherits(r, "error")) {
      skip <- c(skip, sprintf("%s (%s)", uid(x), conditionMessage(r)))
    } else {
      done <- c(done, uid(x))
    }
  }

  left <- uid(DBI::dbGetQuery(con,
    "SELECT b.ID, b.path, b.scaffold FROM assembly_backup b
     JOIN assemblies a ON a.ID = b.ID AND a.path = b.path AND a.scaffold = b.scaffold
     WHERE b.sequence IS NOT NULL AND a.sequence IS NOT b.sequence
       AND NOT EXISTS (SELECT 1 FROM annotate an
         WHERE an.ID = b.ID AND an.path = b.path AND an.scaffold = b.scaffold)
     ORDER BY b.ID, b.path, b.scaffold"))

  txt <- c(
    if (n) sprintf(paste("Moved the edits of %d contig(s) to review copies on Annotation Review; the",
      "Assemble step shows the original contigs. Press Update on Annotation Review to recheck the",
      "unlocked ones."), n),
    if (length(bad)) paste0("Converted with problems: ", paste(bad, collapse = "; "), "."),
    if (length(skip)) paste0("Not converted, left as they were: ", paste(skip, collapse = "; "), "."),
    if (length(done)) paste0("Locked contigs marked done without a recheck: ", paste(done, collapse = ", "), "."),
    if (length(left)) paste0("Trimmed contigs left as they are (dropped or never annotated): ",
      paste(left, collapse = ", "), "."))
  if (length(txt)) message(paste(txt, collapse = "\n"))
  invisible()
}

#' Convert one edited contig of an older project; returns its problems
#'
#' The copy takes the contig's current (trimmed or edited) state, with the
#' partial ends of edited features kept as the user's choice. A trimmed
#' contig's original is written back from its backup row and curated again;
#' when that is not possible the original stays as it is and the copy keeps
#' its full length.
#' @return character, empty when the contig converted cleanly.
#' @noRd
review_upgrade_unit <- function(con, u) {
  id <- u$ID; path <- as.integer(u$path); scaffold <- as.integer(u$scaffold)
  key <- list(id, path, scaffold)
  review_copy(con, id, path, scaffold)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = key)
  DBI::dbExecute(con,
    "UPDATE review_annotations SET partial_start_user = partial_start, partial_stop_user = partial_stop
     WHERE ID = ? AND path = ? AND scaffold = ? AND edited = 1", params = key)
  probs <- character()
  if (u$has_bak == 1 && !identical(u$cur, u$bak)) {
    off <- if (anyNA(c(u$cur, u$bak)) || !nzchar(u$cur)) -1L else
      as.integer(regexpr(u$cur, u$bak, fixed = TRUE))
    if (is.na(u$bak)) {
      probs <- "original not restored: the saved original has no sequence"
    } else if (off < 1L) {
      probs <- "original not restored: the trimmed sequence is not in the saved original"
    } else {
      r <- tryCatch(with_savepoint(con, review_upgrade_original(con, u, off), name = "rp_upgrade_orig"),
        error = function(e) e)
      probs <- if (inherits(r, "error")) paste("original not restored:", conditionMessage(r)) else c(r,
        tryCatch({
          curate_one(con, id, path, scaffold)
          NULL
        }, error = function(e) paste("curation not rerun on the original:", conditionMessage(e))))
    }
  }
  if (u$locked == 1) review_upgrade_done(con, id, path, scaffold) else review_invalidate(con, id, path, scaffold)
  probs
}

#' Write a trimmed contig's original back from its backup row (schema step 4)
#'
#' off is where the trimmed sequence starts in the saved original; it also sets
#' the copy's window. The stored ITSx value may describe the trimmed sequence,
#' so it is cleared.
#' @return a problem when the saved chimera result could not be put back and
#'   the rows were shifted back instead, else character(0).
#' @noRd
review_upgrade_original <- function(con, u, off) {
  id <- u$ID; path <- as.integer(u$path); scaffold <- as.integer(u$scaffold)
  key <- list(id, path, scaffold)
  ann <- if (!is.na(u$bak_ann)) jsonlite::fromJSON(u$bak_ann)
  if (!is.data.frame(ann) && !identical(ann, list())) stop("the saved features cannot be read", call. = FALSE)
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ? WHERE ID = ? AND path = ? AND scaffold = ?",
    params = c(list(u$bak, nchar(u$bak)), key))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = key)
  if (is.data.frame(ann) && nrow(ann)) {
    ann$ID <- id; ann$path <- path; ann$scaffold <- scaffold
    DBI::dbAppendTable(con, "annotations", ann[intersect(names(ann), DBI::dbListFields(con, "annotations"))])
  }
  write_coverage_summary(con, id, path, scaffold, parse_depths(u$bak_cov))
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
    params = key)
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = ?, trim_end = ? WHERE ID = ? AND path = ? AND scaffold = ?",
    params = c(list(off, off + nchar(u$cur) - 1L), key))
  shift_back <- function() if (off > 1L) chimera_shift(con, id, path, scaffold, 1L - off, nchar(u$bak))
  if (is.na(u$bak_chim)) {
    shift_back()
    return(character())
  }
  tryCatch({
    chimera_restore(con, id, path, scaffold, u$bak_chim)
    character()
  }, error = function(e) {
    shift_back()
    paste("chimera result shifted back instead of restored:", conditionMessage(e))
  })
}

#' Mark a converted copy of a locked contig done without a recheck (schema step 4)
#' @noRd
review_upgrade_done <- function(con, id, path, scaffold) {
  cnt <- review_counts(con, id, path, scaffold)
  DBI::dbExecute(con,
    "UPDATE review_units SET check_switch = 2, n_features = ?, missing = ?,
       itsx_warnings = COALESCE(itsx_warnings, 'not_run(not checked before upgrade)'), time_stamp = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(cnt$n_features, cnt$missing, as.integer(Sys.time()), id, as.integer(path),
      as.integer(scaffold)))
}
```

How it fits together:
- `review_copy()` runs first, so the copy captures the trimmed or edited state and the stored check values (including `itsx_warnings`) before the original is written back. The copy keeps `trim_start = 1` and `trim_end = length` unless the write-back succeeds, which is what "untrimmed full window, no write-back" means for a backup without a sequence or one that does not contain the current sequence.
- Every block passed to `with_savepoint()` ends on its value and never calls `return()` itself (a `return()` inside a called function, as in `review_upgrade_original()`, is fine). Nothing inside calls `DBI::dbWithTransaction()`: `review_copy()`, `review_invalidate()`, `curate_one()` (Task 6), `chimera_shift()`, and `chimera_restore()` (Step 3) all use `with_savepoint()`.
- The write-back sits in its own savepoint, so a failure there (contig 9) rolls back only the write-back and the contig still converts with the copy at full length; a failure after it (contig 11) rolls back the whole contig.
- Features go back by column name (`intersect()` with the table's fields), so a snapshot key the table no longer has is ignored. A chimera snapshot that cannot be appended rolls back inside `chimera_restore()` and the rows are shifted back instead, as `restore_assembly()` did when there was no snapshot.
- S11.4 uses `check_switch IS NOT 2`, since a NULL `annotate_switch` derives to NULL. It runs after the conversion, so converted contigs (which have a copy) are not counted twice.
- A trimmed copy's `review_blast_hits` stay in original coordinates (the old trim never shifted `blast_hits`); `review_invalidate()` deletes them on unlocked copies, and nothing reads hit coordinates of a locked copy (only `chimera_unit()` does, and the recheck skips locked contigs).
Check: `Rscript -e 'invisible(parse("R/review_upgrade.R"))'`.

In `R/db_migrate.R` (step 4, as Task 1 wrote it):

Find:

```r
        review_apply_schema(con)
      })
  )
```

Replace with:

```r
        review_apply_schema(con)
        review_upgrade(con)
        DBI::dbExecute(con, "DROP TABLE assembly_backup")
      })
  )
```

Step 4 is not released yet (Task 1 added it on this branch), so extending it is allowed; steps 1-3 stay as they are.

In `R/db_schema.R`:

Find:

```r
#' DDL for the assembly_backup table (original pre-trim state, for undo/restore).
#' @noRd
ASSEMBLY_BACKUP_DDL <- "CREATE TABLE IF NOT EXISTS assembly_backup (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  sequence TEXT,
  length INTEGER,
  annotations TEXT,
  coverage TEXT,
  chimera TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for the reference table (one row per clade reference-database group).
```

Replace with:

```r
#' DDL for the reference table (one row per clade reference-database group).
```

In `R/init_db.R` (`new_db()`; the annotations lines are Task 1's):

Find:

```r
  # assembly_backup (original pre-trim state so a trim can be undone) ----
  DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)

  # annotations (nrDNA feature model per assembly unit; see ANNOTATIONS_DDL) ----
```

Replace with:

```r
  # annotations (nrDNA feature model per assembly unit; see ANNOTATIONS_DDL) ----
```

In `R/pipeline_steps.R` (`UNIT_EVIDENCE_TABLES` as Task 2 left it):

Find:

```r
  "blast_ref_alignment", "assembly_backup", REVIEW_TABLES)
```

Replace with:

```r
  "blast_ref_alignment", REVIEW_TABLES)
```

Find (in `load_assemblies()`):

```r
  # A (re)assembly invalidates any pre-trim undo snapshots for this sample.
  tryCatch(DBI::dbExecute(con, "DELETE FROM assembly_backup WHERE ID = ?", params = list(ID)),
    error = function(e) NULL)
  # triage closes the sample, and it never runs for one with no unit to annotate
```

Replace with:

```r
  # triage closes the sample, and it never runs for one with no unit to annotate
```

In `tests/testthat/helper-app-fixture.R` (`fixture_fill_schema()`):

Find:

```r
                   RiboPilot:::SCREEN_DDL, RiboPilot:::ASSEMBLY_BACKUP_DDL,
                   RiboPilot:::REFERENCE_DDL, RiboPilot:::REFERENCE_CONFIG_DDL,
```

Replace with:

```r
                   RiboPilot:::SCREEN_DDL, RiboPilot:::REFERENCE_DDL, RiboPilot:::REFERENCE_CONFIG_DDL,
```

In `tests/testthat/test-load-assemblies.R` (`mk_asm_db()`; nothing there reads the table any more):

Find:

```r
  DBI::dbExecute(con, "CREATE TABLE assembly_backup(ID TEXT, path INT, scaffold INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
```

Replace with:

```r
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
```

Check: `Rscript -e 'for (f in c("R/db_migrate.R", "R/db_schema.R", "R/init_db.R", "R/pipeline_steps.R", "tests/testthat/helper-app-fixture.R", "tests/testthat/test-load-assemblies.R")) invisible(parse(f))'`.

Check that only the conversion still names the table:
Run: `grep -rn "assembly_backup\|ASSEMBLY_BACKUP_DDL\|ensure_assembly_backup\|chimera_snapshot" R/ inst/ tests/testthat/*.R`
Expected: only `R/review_upgrade.R` (its header comment and the two queries), `R/db_migrate.R` (the `DROP TABLE` line), `tests/testthat/test-review-upgrade.R` (it builds old projects with the table), and the new `expect_false("assembly_backup" ...)` line in `tests/testthat/test-backwards-compat.R`. `tests/testthat/fixtures/schema_v1.sql` keeps its table (not matched by `tests/testthat/*.R`). Any other hit (for example `R/app_modals.R` or `tests/testthat/test-export.R`) means Task 15 or Task 17 was not applied; stop and report it.

- [ ] **Step 8: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); for (f in c("review-upgrade", "backwards-compat", "trim-coverage")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: PASS. test-review-upgrade.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 34 ]`, test-backwards-compat.R `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 39 ]`, test-trim-coverage.R `PASS 33`.

Then the files that upgrade a version 1 database through step 4, build `new_db()`, run `load_assemblies()` or `curate_one()`, or use the app fixtures (the fixtures no longer create `assembly_backup`, and nothing reads it):
Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); for (f in c("review-schema", "review-copy", "load-assemblies", "init-assemblies", "select-paths", "curate", "chimera", "pipeline-drivers", "write-atomicity", "batch-trim", "export", "app-export", "app-table", "app-scaffold-uid", "reference-schema", "nextflow", "no-duplicate-definitions", "docs-site")) testthat::test_file(sprintf("tests/testthat/test-%s.R", f))'`
Expected: every file `FAIL 0`.

Then the whole suite once: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'`
Expected: `FAIL 0`. On a heavily loaded machine the two `open_path` tests in test-app-utils.R can time out waiting for their marker file; rerun that file alone, where they pass.

- [ ] **Step 9: Commit**

```bash
git add R/review_upgrade.R R/chimera.R R/db_migrate.R R/db_schema.R R/init_db.R R/pipeline_steps.R tests/testthat/test-review-upgrade.R tests/testthat/test-trim-coverage.R tests/testthat/test-backwards-compat.R tests/testthat/helper-app-fixture.R tests/testthat/test-load-assemblies.R
git commit -m "Schema step 4 converts old trims to review copies; drop assembly_backup"
```

From here on `backwards_compatibility()` may run on real projects (Global Constraints). Run it first on a copied `.sqlite` of a real version 4 project and read its message.


### Task 23: Docs and stale comments

**Files:**
- Modify: `vignettes/Test-Project-Annotate.Rmd` (lines 81, 83, 84)
- Modify: `vignettes/GenBank-Submission.Rmd` (lines 150-154, 183-187)
- Modify: `vignettes/Test-Project-Export.Rmd` (lines 42-48)
- Modify: `vignettes/Test-Project-Assemble.Rmd` (line 83)
- Modify: `vignettes/Curation-Flags.Rmd` (lines 45-47, 68, 70)
- Modify: `vignettes/Citing.Rmd` (line 51)
- Modify: `vignettes/Caches-and-Offline.Rmd` (line 42)
- Modify: `vignettes/RiboPilot.Rmd` (after the state table, line 131)
- Modify: `vignettes/Troubleshooting.Rmd` (lines 49-50, 59-60, 81, 102-103 of the committed file), only the edits the user approves in Step 8
- Modify: `R/constants.R` (`SCAFFOLD_COL_HELP$itsx_warnings`, lines 150-151)
- Modify: `inst/nextflow/main.nf` (comment, lines 20-24)
- Modify: `inst/nextflow/modules/pipeline_r.nf` (comment above `query_annotate_units`, lines 47-51)
- Modify: `NEWS.md` (lines 64 and 116-117; new section before `## Re-running samples` at line 59)
- Modify: `dev/design/01_ARCHITECTURE.md` (lines 106-107, 114-121)
- Modify: `dev/design/15_WF2_DECOMMISSION.md` (after line 3)
- Test: `tests/testthat/test-docs-site.R`

**Interfaces:**
- Consumes: no code. The text describes what earlier tasks built: the Update button on Annotation Review and its run window (Task 13); the State column's "Ready to run", "On hold", and "Waiting for Assemble Update" and the Checks item "ITSx not checked" (Task 14); Restore deleting the copy (Task 15); State offering Ready to run and On hold (Task 16); the chimera panel's "Not run yet. Press Update." (Task 17); the removed export ITSx run and the Export done window's flagged / clean / not checked counts (Task 18); the lock refusal "Press Update first." (Task 21); the not-run strings `not_run(ITSx did not finish)`, `not_run(ITSx check failed: <message>)` (Task 5), `not_run(not checked before upgrade)` (Task 22), and `not_run(not checked)` (Task 18); the processes `query_recheck`, `recheck_tools`, `recheck`, the `recheck_tools` out-of-memory retry, and the `.recheck/` and `work_recheck/` folders (Task 11).
- Produces: the list of screenshots to retake (Step 12), which Task 24 repeats.

The user docs still describe the old design: trims recalculate flags at once, Restore only undoes a trim, and export runs a "final ITSx check" that writes `<group>_itsx_check.csv` and `itsx_check/`. This task brings the pages SPEC S13 names in line with S5.2 (edits mark checks not run; Restore deletes the copy), S6.2 (lock gate, export backstop), S6.3 (ITSx states), S7 (Update on Annotation Review), and S9 (no ITSx run at export, stored flags, batch-size caveat, `not_run` legend), drops "updated again by export" from the ITSx Warnings tooltip, and fixes the stale comments S10 assigns here (`inst/nextflow/main.nf:20-24` and `pipeline_r.nf:47-51` claim a trim queues re-annotation; no code does). It adds the NEWS entry and the two design-doc notes of S13. `vignettes/Troubleshooting.Rmd` holds uncommitted edits of the user's own, so its four small edits happen only if the user agrees (Step 8), and only those lines are committed. `vignettes/Your-Own-Project.Rmd` (also holding user edits) needs no change and is never touched. Two additions outside the S13 list: the project-folder table in Caches-and-Offline.Rmd gains the recheck's `work_recheck/` and `.recheck/` folders, and the getting-started page (RiboPilot.Rmd), whose state table lists five states, gains a note on what State offers on Annotation Review and on the new "Waiting for Assemble Update" state.

- [ ] **Step 1: Write the failing test**

In `tests/testthat/test-docs-site.R`, append a test after the last one.

Find:

```r
  expect_match(format(cit, style = "text"), "9.9.9", fixed = TRUE)
})
```

Replace with:

```r
  expect_match(format(cit, style = "text"), "9.9.9", fixed = TRUE)
})

test_that("docs describe the stored ITSx check and Update on Annotation Review", {
  root <- testthat::test_path("../..")
  vig <- file.path(root, "vignettes")
  skip_if_not(dir.exists(vig), "package sources not available")
  flat <- function(f) gsub("\\s+", " ", paste(readLines(f, warn = FALSE), collapse = " "))
  pages <- c(list.files(vig, "[.]Rmd$", full.names = TRUE), file.path(root, "README.Rmd"))
  docs <- vapply(pages, flat, character(1))
  stale <- c("final ITSx check", "_itsx_check.csv", "itsx_check/",
    "keep the values from the last pipeline run", "quality flags are recalculated right away")
  bad <- as.character(unlist(lapply(stale, function(b) {
    hit <- grepl(b, docs, fixed = TRUE)
    if (any(hit)) paste0(b, " (", basename(pages[hit]), ")")
  })))
  expect_equal(bad, character(0))
  gb <- flat(file.path(vig, "GenBank-Submission.Rmd"))
  legend <- c("`not_run(ITSx did not finish)`", "`not_run(ITSx check failed: <message>)`",
    "`not_run(not checked before upgrade)`", "`not_run(not checked)`")
  expect_equal(legend[!vapply(legend, grepl, logical(1), x = gb, fixed = TRUE)], character(0))
  ann <- flat(file.path(vig, "Test-Project-Annotate.Rmd"))
  expect_true(grepl("press **Update** in the Annotation Review toolbar", ann, fixed = TRUE))
  expect_true(grepl("\"Press Update first.\"", ann, fixed = TRUE))
  expect_true(grepl("it undoes every trim, gene boundary edit, and 5'/3' choice", ann, fixed = TRUE))
  tip <- RiboPilot:::SCAFFOLD_COL_HELP$itsx_warnings
  expect_false(grepl("export", tip, fixed = TRUE))
  expect_true(grepl("'not checked'", tip, fixed = TRUE))
})
```

- [ ] **Step 2: Run it to see it fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'`
Expected: FAIL 7, all in the new test (PASS 237 on the repo as of this plan; earlier tasks may shift that count). The new test fails on: `bad` listing "final ITSx check (Citing.Rmd)", "(Curation-Flags.Rmd)", "(GenBank-Submission.Rmd)", "(Test-Project-Assemble.Rmd)", "(Test-Project-Export.Rmd)", "_itsx_check.csv (GenBank-Submission.Rmd)", "itsx_check/ (GenBank-Submission.Rmd)", "keep the values from the last pipeline run (Test-Project-Annotate.Rmd)", and "quality flags are recalculated right away (Test-Project-Annotate.Rmd)"; all four legend strings missing; all three `grepl(...)` checks on Test-Project-Annotate.Rmd FALSE; `grepl("export", tip, ...)` TRUE; `grepl("'not checked'", tip, ...)` FALSE. The five existing tests still pass.

- [ ] **Step 3: Update the vignettes**

Every Find below occurs exactly once in its file. Lines 80 and 83 of Test-Project-Annotate.Rmd end with a space; the Finds there stop before it, so leave it.

In `vignettes/Test-Project-Annotate.Rmd`:

Find:

```markdown
to undo a trim.
```

Replace with:

```markdown
to undo a trim. Restore returns the contig to its original sequence, features, and
check results: it undoes every trim, gene boundary edit, and 5'/3' choice on the
contig. Your Reviewed, Exclude, and Partial marks, notes, and export group stay
as they are.
```

Find:

```markdown
the contig's quality flags are recalculated right away.
```

Replace with:

```markdown
the contig's checks are marked "not run".
```

Find:

```markdown
BLAST, chimera checks, ITSx, and coverage results keep the values from the last pipeline run.
```

Replace with:

```markdown
Trims and edits change a working copy of the contig; the original, which the
Assemble step shows, keeps the pipeline's results. Until the checks run again,
the contig's State shows **Ready to run**, and its check columns and **Checks**
say "not run". To run them, press **Update** in the Annotation Review toolbar,
then **Run from App** (on a cluster, the window offers a job script, as in the
Assemble step). Update reruns every check that works offline on each queued
contig: feature counts, how complete each gene is, the ITSx checks, the BLAST
hits (a new search of the bundled rDNA database when the contig was found
there, otherwise its original hits cut to the trimmed contig), the BLAST
screen, the reference alignment, the quality flags, and the chimera check. It
does not re-annotate the contig, and it never contacts NCBI.

**State** in Annotation Review offers **Ready to run**, which queues a contig
for the next Update, and **On hold**, which keeps it out. A contig the
pipeline has not finished annotating shows **Waiting for Assemble Update**:
press **Update** in the Assemble step first.

**Lock** refuses a contig whose checks are not finished, and says what to do
next, for example "Press Update first." **Exclude & Lock** always works.
Export leaves out any contig whose checks are not finished.

If a contig you edited is annotated again in the Assemble step, for example
after you change its annotation options, it keeps your edits and its checks are
queued again. **Restore** then gives you the new annotation.
```

In `vignettes/GenBank-Submission.Rmd`:

Find:

```markdown
**The final ITSx check.** NCBI runs a standard ITSx check on every rRNA and ITS
submission. It may suggest that you remove records judged chimeric, gapped, or misassembled by ITSx.
RiboPilot runs the same check on your records after writing them. The result is printed at export and
recorded in the `itsx_check` column of `<group>_summary.csv` and in
`<group>_itsx_check.csv`. It only warns: nothing is removed.
```

Replace with:

```markdown
**The ITSx check.** NCBI runs a standard ITSx check on every rRNA and ITS
submission. It may suggest that you remove records judged chimeric, gapped, or misassembled by ITSx.
RiboPilot runs the same stock ITSx check on every contig when it annotates it, and again when
**Update** on Annotation Review rechecks a contig you trimmed or edited. Export does not run
ITSx: it reports the stored result in the `itsx_check` column of `<group>_summary.csv`, and the
Export done window counts the flagged, clean, and not checked records and names the flagged
ones. The **ITSx Warnings** column in the app shows the same flags. It only warns: nothing is
removed.

RiboPilot checks one contig at a time, while NCBI runs ITSx on your whole upload at once.
ITSx's calls can shift a little with the set of sequences it is given, so on borderline
records NCBI's check may differ from RiboPilot's.

An `itsx_check` value that starts with `not_run` means ITSx gave no answer for that record.
The record still counts as checked, so it can be locked and exported, and NCBI checks it
again when you upload. To check it again in RiboPilot, unlock the contig in Annotation
Review, set its State to Ready to run, press **Update**, and lock it again once its checks
are finished.

| Value | What it means |
|---|---|
| `not_run(ITSx did not finish)` | ITSx crashed or was stopped before it finished this contig. |
| `not_run(ITSx check failed: <message>)` | Reading the ITSx result failed; the message says why. |
| `not_run(not checked before upgrade)` | The contig was locked in an older RiboPilot version that stored no ITSx check for it. |
| `not_run(not checked)` | No ITSx check is stored for the record. |
| any other `not_run(...)` | Stored by an older RiboPilot version, for example `not_run(no ITSx output)`. |
```

Find:

```text
    <group>_itsx_check.csv  what stock ITSx says about each record
    README.txt              how to submit
    samples/                one .fasta and .tbl per record
    markers/                one folder per marker (18S, ITS1, and so on)
    itsx_check/             raw ITSx output behind the check
```

Replace with:

```text
    README.txt              how to submit
    samples/                one .fasta and .tbl per record
    markers/                one folder per marker (18S, ITS1, and so on)
```

In `vignettes/Test-Project-Export.Rmd`:

Find:

```markdown
The export ends with a summary. If any exported record still carries a
submission-risk flag, the summary names it and why, above the result of the
final ITSx check; it is green only when nothing was flagged. In the test
project CHIM01 was left out, so the summary only shows the final ITSx check,
which flags one record:

![Export summary, with the final ITSx check flagging one record](figures/get-started/export-done.png)
```

Replace with:

```markdown
The export ends with a summary. If any exported record still carries a
submission-risk flag, the summary names it and why, above the ITSx counts: how
many records the stored ITSx check flagged, how many are clean, and how many
were not checked, with the flagged records named. It is green only when nothing
was flagged. Export does not run ITSx itself; see
[the ITSx check](GenBank-Submission.html#check-before-you-submit). In the test
project CHIM01 was left out, so the summary only shows the ITSx counts:

![Export summary, with the ITSx counts](figures/get-started/export-done.png)
```

In `vignettes/Test-Project-Assemble.Rmd`:

Find:

```markdown
  [the final ITSx check](GenBank-Submission.html#check-before-you-submit).
```

Replace with:

```markdown
  [the ITSx check](GenBank-Submission.html#check-before-you-submit).
```

In `vignettes/Curation-Flags.Rmd`:

Find:

```markdown
check finds for the contig (see [the final ITSx check](GenBank-Submission.html#check-before-you-submit)).
For example
`itsx_chimeric`, `itsx_no_5.8S`, or `itsx_boundary_shift(...)`. The **Chimera** column reports the
```

Replace with:

```markdown
check finds for the contig (see [the ITSx check](GenBank-Submission.html#check-before-you-submit)).
For example
`itsx_chimeric`, `itsx_no_5.8S`, or `itsx_boundary_shift(...)`. A grey "not run" means ITSx gave
no answer for the contig, and a grey "not checked" means no ITSx check is stored yet; the Checks
column lists both as "ITSx not checked". The **Chimera** column reports the
```

Find:

```markdown
| minor chimera, uncertain screen, ITSx, warnings, feature flags | amber | Worth a look, but nothing asks when you Lock or Export. |
```

Replace with:

```markdown
| minor chimera, uncertain screen, ITSx, ITSx not checked, warnings, feature flags | amber | Worth a look, but nothing asks when you Lock or Export. |
```

Find:

```markdown
| no reference | grey | Nothing is flagged, but the BLAST screen had no usable match to judge by. |
```

Replace with:

```markdown
| no reference | grey | Nothing is flagged, but the BLAST screen had no usable match to judge by. |

On Annotation Review, a contig whose checks are not finished shows its state in the
Checks column and in every check column (Warnings, ITSx Warnings, Screen, Chimera,
and the BLAST columns) instead, for example "not run" after a trim, until
**Update** on Annotation Review reruns its checks.
```

In `vignettes/Citing.Rmd`:

Find:

```markdown
Finds ITS1 and ITS2, and runs the final ITSx check at export.
```

Replace with:

```markdown
Finds ITS1 and ITS2, and runs the stock ITSx check on each contig.
```

In `vignettes/Caches-and-Offline.Rmd`:

Find:

```markdown
| `work_ref/`, `.reference/` | Working files of reference builds. | Safe to delete once every group is Ready. |
```

Replace with:

```markdown
| `work_ref/`, `.reference/` | Working files of reference builds. | Safe to delete once every group is Ready. |
| `work_recheck/`, `.recheck/` | Working files of the recheck that **Update** on Annotation Review runs. | Safe to delete when no recheck is running. |
```

In `vignettes/RiboPilot.Rmd`:

Find:

```markdown
</table>
```

Replace with:

```markdown
</table>

On Annotation Review, **State** offers only Ready to run, which queues a contig
for that step's own **Update**, and On hold, which keeps it out. A contig the
pipeline has not finished annotating shows a sixth state there, **Waiting for
Assemble Update**; see [Test project: Annotation Review](Test-Project-Annotate.html).
```

- [ ] **Step 4: Update the ITSx Warnings tooltip**

In `R/constants.R` (`SCAFFOLD_COL_HELP$itsx_warnings`):

Find:

```r
    "the annotation step and updated again by export. Empty means clean; 'not run'",
    "means ITSx produced no output."),
```

Replace with:

```r
    "the annotation step, and by Update on Annotation Review after a trim or an edit.",
    "Empty means clean; 'not run' means ITSx gave no answer (hover for why);",
    "'not checked' means no ITSx check is stored yet."),
```

- [ ] **Step 5: Run the tests to see them pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'`
Expected: PASS, FAIL 0 (the Step 2 count plus the 7 new expectations: 244 on the repo as of this plan).

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-ui-helpers.R")'`
Expected: PASS, FAIL 0 (MP_COL_TIPS reuses the tooltip; it must stay non-empty).

Run: `grep -nP '[^\x00-\x7F]' vignettes/Test-Project-Annotate.Rmd vignettes/GenBank-Submission.Rmd vignettes/Test-Project-Export.Rmd vignettes/Test-Project-Assemble.Rmd vignettes/Curation-Flags.Rmd vignettes/Citing.Rmd vignettes/Caches-and-Offline.Rmd vignettes/RiboPilot.Rmd R/constants.R`
Expected: no output (ASCII only).

- [ ] **Step 6: Fix the stale comments in the Nextflow files**

In `inst/nextflow/main.nf`:

Find:

```groovy
// Re-annotation lives here too. There is no second workflow: units queued after
// review (a trim, an options change, a restored scaffold, an edited sequence)
// enter through query_annotate_units and skip straight to the annotation chain,
// paying for neither assembly nor a second remote BLAST. See
// dev/design/15_WF2_DECOMMISSION.md.
```

Replace with:

```groovy
// Re-annotation lives here too: units queued for annotation alone (an
// annotation or curation options change, a scaffold put back) enter through
// query_annotate_units and skip straight to the annotation chain, paying for
// neither assembly nor a second remote BLAST. See
// dev/design/15_WF2_DECOMMISSION.md. A trim or a boundary edit on Annotation
// Review queues nothing here: it changes the contig's review copy, and
// main_recheck.nf reruns the checks on that copy without re-annotating it.
```

In `inst/nextflow/modules/pipeline_r.nf`:

Find:

```groovy
// WF2: assembled samples queued for (re-)annotation, with their annotate_opts
// compute settings.
// Units queued for annotation ONLY: a trim, an annotate_opts or curate_opts
// change, a restored scaffold, an edited sequence. WF2 used to own this; there
// is one workflow now (dev/design/15_WF2_DECOMMISSION.md).
```

Replace with:

```groovy
// Assembled units queued for (re-)annotation, with their annotate_opts
// compute settings.
// Units queued for annotation ONLY: an annotate_opts or curate_opts change, or
// a scaffold put back. WF2 used to own this; it is an arm of main.nf now
// (dev/design/15_WF2_DECOMMISSION.md). Trims and boundary edits change review
// copies, which main_recheck.nf rechecks; they never queue units here.
```

Run: `grep -n -E 'a trim, an|an edited sequence|restored scaffold' inst/nextflow/main.nf inst/nextflow/modules/pipeline_r.nf`
Expected: no output.

- [ ] **Step 7: NEWS entry and design-doc notes**

In `NEWS.md`:

Find:

```markdown
* Changing annotation options re-runs only annotation, and keeps trims.
```

Replace with:

```markdown
* Changing annotation options re-runs only annotation. A contig you trimmed or
  edited on Annotation Review keeps its edits, and its checks are queued for
  Update on Annotation Review; Restore takes the new annotation.
```

Find:

```markdown
* The chimera check plot, segments, and windows redraw after a trim, a Restore,
  or a feature boundary edit, instead of only when another unit is picked.
```

Replace with:

```markdown
* The chimera check panel follows the contig's edits: after a trim or a feature
  boundary edit it says "Not run yet. Press Update." until the recheck has run,
  and after a Restore it shows the original check again.
```

Find:

```markdown
## Re-running samples
```

Replace with:

```markdown
## Annotation Review edits and Update

* Trims, gene boundary edits, and 5'/3' choices on Annotation Review change a
  working copy of the contig. The original contig and its check results stay as
  the pipeline made them, and the Assemble step shows only the original.
* A trim or a boundary edit marks the contig's checks "not run". The new
  **Update** button on Annotation Review reruns every offline check on those
  contigs, one contig at a time: feature counts, gene completeness, ITSx, the
  BLAST hits, the BLAST screen, the reference alignment, the quality flags, and
  the chimera check. It never re-annotates and never contacts NCBI. Quality
  flags are no longer recalculated in the app right after an edit.
* State on Annotation Review offers only Ready to run and On hold, and queues
  contigs for that Update. A contig whose original still waits for the
  Assemble step's Update says so.
* Lock on Annotation Review refuses a contig whose checks are not finished and
  says what to do next. Exclude & Lock still works. Export leaves such contigs
  out and counts them.
* Restore deletes the working copy, so the contig shows its original sequence,
  features, and checks again. It undoes trims, boundary edits, and 5'/3'
  choices; marks, notes, and export group stay. Restore on the toolbar works on
  any contig with a working copy, not only trimmed ones.
* Export no longer runs ITSx. The `itsx_check` column of `<group>_summary.csv`
  and the Export done window report the ITSx check stored for each contig, and
  `<group>_itsx_check.csv` and the `itsx_check/` folder are gone. An ITSx run
  that crashed reads `not_run(ITSx did not finish)`, never clean, and the
  Checks column shows it as "ITSx not checked".
* Existing projects need `backwards_compatibility()` once. It turns earlier
  trims and edits into working copies and puts the original contigs back.

## Re-running samples
```

In `dev/design/01_ARCHITECTURE.md`:

Find:

```markdown
state as it is; `query_annotate_units` takes those units without
re-assembling, so trims and edits survive.
```

Replace with:

```markdown
state as it is; `query_annotate_units` takes those units without
re-assembling. Re-annotation rewrites only the original tables. Trims and
edits live in review copies (the `review_*` tables), which it never rewrites:
an unlocked edited copy keeps its edits and is queued for the recheck, an
unedited copy is dropped, and a locked copy is left alone (see
`dev/specs/2026-10-08-contig-recheck-design.md`, section 5.3).
```

Find:

```markdown
## Nextflow workflow

ONE workflow, `inst/nextflow/main.nf`. It takes two kinds of work: samples
queued for assembly, and units queued for annotation only (a trim, an options
change, a restored scaffold, an edited sequence). The second kind skips
preprocess, assembly, scoring, the remote BLAST and coverage. There is no WF2:
`main_annotate.nf` and `run_annotate()` were deleted, see
`15_WF2_DECOMMISSION.md`. Modules mirror MitoPilot's stage split:
```

Replace with:

```markdown
## Nextflow workflows

Two workflows work on a project's samples; a third, `main_reference.nf`,
builds reference databases.

`inst/nextflow/main.nf`, launched by Update on the Assemble step, takes two
kinds of work: samples queued for assembly, and units queued for annotation
only (an annotation or curation options change, a scaffold put back). The
second kind skips preprocess, assembly, scoring, the remote BLAST and
coverage. There is no WF2: `main_annotate.nf` and `run_annotate()` were
deleted, see `15_WF2_DECOMMISSION.md`.

`inst/nextflow/main_recheck.nf`, launched by Update on Annotation Review,
reruns the offline checks on the review copies of contigs the user trimmed or
edited, one contig at a time (`query_recheck -> recheck_tools -> recheck`). It
never re-annotates and never contacts NCBI. It calls the same R check
functions (with `review = TRUE`) and the same `bin/rp_rrna_tools.sh` as
main.nf, so there is no second copy of the pipeline to drift. See
`dev/specs/2026-10-08-contig-recheck-design.md`.

main.nf's modules mirror MitoPilot's stage split:
```

In `dev/design/15_WF2_DECOMMISSION.md`:

Find:

```markdown
Written 2026-07-31. Branch `wf2-decommission`, cut from `main`.
```

Replace with:

```markdown
Written 2026-07-31. Branch `wf2-decommission`, cut from `main`.

## Status: a second workflow is back, with a different job

The contig recheck (`dev/specs/2026-10-08-contig-recheck-design.md`) adds
`inst/nextflow/main_recheck.nf` and an Update button on Annotation Review, so
decision 4 below ("One Update button") no longer holds. Re-annotation still
has one entry point, the annotate-only arm of main.nf. The recheck never
re-annotates: it reruns the offline checks on the review copies of contigs the
user trimmed or edited, through the same R check functions and the same
`bin/rp_rrna_tools.sh` as main.nf, so there is no second copy to drift. A trim
never queued re-annotation, although comments in main.nf and pipeline_r.nf
said so; they now point to main_recheck.nf.
```

Run: `grep -nP '[^\x00-\x7F]' NEWS.md dev/design/01_ARCHITECTURE.md dev/design/15_WF2_DECOMMISSION.md inst/nextflow/main.nf inst/nextflow/modules/pipeline_r.nf`
Expected: no output.

- [ ] **Step 8: Ask the user about Troubleshooting.Rmd**

`vignettes/Troubleshooting.Rmd` holds uncommitted edits of the user's own. Do not edit, stage, or revert it before the user answers. Stop and ask the user (if you are a subagent, report back and have the controller ask), word for word:

> vignettes/Troubleshooting.Rmd has uncommitted edits of yours. The recheck adds three Nextflow steps and a new kind of run report that the page does not mention. May I make these small edits and commit only my lines? Your own edits stay uncommitted, as they are now.
> 1. Process table: a row for `query_recheck`, `recheck_tools`, and `recheck`.
> 2. Out-of-memory sentence: `recheck_tools` joins the steps that retry with more memory.
> 3. Report description: the Workflow line can name a recheck.
> 4. Report description: after a recheck the report lists failed contigs instead of failed samples.
> Answer yes, no, or the numbers you want.

If the answer is no, skip Step 9. Otherwise make in Step 9 only the edits the user approved.

- [ ] **Step 9: Edit Troubleshooting.Rmd and stage only these lines**

The approved edits go into two files: a copy of the committed file, which is what gets staged, and the user's working copy, so the user sees them. Before any edit, save the user's own changes and the committed file:

```bash
mkdir -p /tmp/claude-1000/rp-t23
git diff -U0 vignettes/Troubleshooting.Rmd | grep -E '^[+-]' | grep -v -E '^(\+\+\+|---) ' > /tmp/claude-1000/rp-t23/user_before.txt
git show HEAD:vignettes/Troubleshooting.Rmd > /tmp/claude-1000/rp-t23/Troubleshooting.Rmd
```

Then make each approved edit below, with the Edit tool, in both `/tmp/claude-1000/rp-t23/Troubleshooting.Rmd` and `vignettes/Troubleshooting.Rmd`. Each Find occurs once in both. Edit 1:

Find:

```markdown
| `triage` | Proposes which contigs to keep. |
```

Replace with:

```markdown
| `triage` | Proposes which contigs to keep. |
| `query_recheck`, `recheck_tools`, `recheck` | The recheck that **Update** on Annotation Review runs: pick the queued contigs, run barrnap and ITSx on each, then rerun its checks. |
```

Edit 2:

Find:

```markdown
`find_rrna_scan`, and `annotate_tools`) retry by themselves after codes 104, 134, 137, 139, 140, 143,
```

Replace with:

```markdown
`find_rrna_scan`, `annotate_tools`, and `recheck_tools`) retry by themselves after codes 104, 134, 137, 139, 140, 143,
```

Edit 3:

Find:

```markdown
- **Workflow** through **Log**: which workflow ran (the pipeline, or a
  reference build), when, how it ended and with what exit status, whether it
```

Replace with:

```markdown
- **Workflow** through **Log**: which workflow ran (the pipeline, a recheck, or
  a reference build), when, how it ended and with what exit status, whether it
```

Edit 4:

Find:

```markdown
- **Failed samples**: samples that finished without a usable assembly, with the
  pipeline's note and the likely cause.
```

Replace with:

```markdown
- **Failed samples**: samples that finished without a usable assembly, with the
  pipeline's note and the likely cause. After a recheck (**Update** on
  Annotation Review), **Failed contigs** lists instead each contig whose
  recheck failed or did not finish, with its note.
```

Then stage the edited committed copy, never the working copy, and check that the user's own changes are untouched:

```bash
git update-index --cacheinfo 100644,$(git hash-object -w /tmp/claude-1000/rp-t23/Troubleshooting.Rmd),vignettes/Troubleshooting.Rmd
git diff --cached vignettes/Troubleshooting.Rmd
git diff -U0 vignettes/Troubleshooting.Rmd | grep -E '^[+-]' | grep -v -E '^(\+\+\+|---) ' > /tmp/claude-1000/rp-t23/user_after.txt
diff /tmp/claude-1000/rp-t23/user_before.txt /tmp/claude-1000/rp-t23/user_after.txt && echo "user edits intact"
```

Expected: `git diff --cached` shows only the approved lines (none of the user's own changes, such as "What each part means:" becoming "What each section means:"); the last command prints `user edits intact`. If it does not, run `git restore --staged vignettes/Troubleshooting.Rmd` (this touches only the index, never the user's file) and stop and report.

- [ ] **Step 10: Run the tests and render the changed pages**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-docs-site.R")'`
Expected: PASS, FAIL 0 (the test reads the working copy of Troubleshooting.Rmd, user edits included; it holds none of the stale phrases).

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-ui-helpers.R")'`
Expected: PASS, FAIL 0.

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-nextflow.R")'`
Expected: PASS, FAIL 0 (comment edits only; the file parses main.nf and pipeline_r.nf).

Run:

```bash
Rscript -e 'd <- tempfile(); dir.create(d); for (f in c("Test-Project-Annotate", "GenBank-Submission", "Test-Project-Export", "Test-Project-Assemble", "Curation-Flags", "Citing", "Caches-and-Offline", "RiboPilot", "Troubleshooting")) rmarkdown::render(file.path("vignettes", paste0(f, ".Rmd")), output_dir = d, intermediates_dir = d, quiet = TRUE); h <- paste(readLines(file.path(d, "GenBank-Submission.html"), warn = FALSE), collapse = " "); stopifnot(grepl("<td><code>not_run(not checked before upgrade)</code></td>", h, fixed = TRUE)); cat("rendered\n")'
git status --short vignettes
```

Expected: `rendered` (the legend renders as a table), and `git status --short vignettes` shows no `??` entries (render writes only to the temp folder).

- [ ] **Step 11: Commit**

Never `git add` `vignettes/Troubleshooting.Rmd` or `vignettes/Your-Own-Project.Rmd`, and never `git add -A` or `git add .`: that would stage the user's own edits.

```bash
git add tests/testthat/test-docs-site.R vignettes/Test-Project-Annotate.Rmd vignettes/GenBank-Submission.Rmd vignettes/Test-Project-Export.Rmd vignettes/Test-Project-Assemble.Rmd vignettes/Curation-Flags.Rmd vignettes/Citing.Rmd vignettes/Caches-and-Offline.Rmd vignettes/RiboPilot.Rmd R/constants.R inst/nextflow/main.nf inst/nextflow/modules/pipeline_r.nf NEWS.md dev/design/01_ARCHITECTURE.md dev/design/15_WF2_DECOMMISSION.md
git diff --cached --name-only
git commit -m "Docs for the contig recheck; fix stale comments"
git status --short
```

Expected: `git diff --cached --name-only` lists those 15 files, plus `vignettes/Troubleshooting.Rmd` only when Step 9 ran. After the commit, `git status --short` still shows ` M vignettes/Troubleshooting.Rmd` and ` M vignettes/Your-Own-Project.Rmd` (the user's own edits) and nothing else from this task.

- [ ] **Step 12: Report the screenshots to retake**

The figures show the old app; this task does not retake them. End your report with this list for the user:

- `vignettes/figures/get-started/annotation-table.png`: the Annotation Review toolbar now has **Update**, and the State column can show "Waiting for Assemble Update".
- `vignettes/figures/get-started/export-window.png`: the "Final ITSx check" box is gone.
- `vignettes/figures/get-started/export-done.png`: the summary shows flagged, clean, and not checked counts instead of the final ITSx check.


### Task 24: Full verification

**Files:**
- Create: nothing in the repository. Scratch only: helper scripts, logs, and screenshots under `/tmp/rp-t24/`; a fresh test project at `~/Documents/DJMCompBio/RiboPilot_test/recheck_t24`; a local Docker image `ribopilot:recheck`.
- Modify: nothing.
- Test: the whole suite (`tests/testthat/`, run on an archive copy of HEAD), mutation runs on a second copy, upgrades of copies of the local test projects, and a scripted drive of the app.

**Interfaces:**
- Consumes (all from Tasks 1-23, by the names and texts the contract fixes):
  - `current_schema_version()` = 5; `REVIEW_VIEWS_DDL` (5 views); tables `review_units`, `review_annotations`, `review_screen`; `review_units_v` columns `has_copy`, `check_switch`, `itsx_warnings`, `length` (Task 1).
  - Copy semantics: a trim makes a copy with `edited = 1`, `check_switch = 1`, the window in `trim_start`/`trim_end`, and the cut end in `partial_start_user`/`partial_stop_user` (Tasks 2, 4, 15).
  - Bin and R internals in the image: `recheck_unit`, `review_follow_original`, `itsx_run_finished`, `itsx_store_not_run`, `review_copy` (Tasks 2, 5, 10, 12).
  - App ids and texts: `annot_run_btn`, `rc_nf_start`, `rc_nf_close`, run window title `RiboPilot recheck - %d contig(s) to update` without a Resume box (Task 13); the Annotation Review State column keeps the field name `annotate_switch` and holds the copy's check state (Task 14); Lock refusal text `Press Update first.` (Task 21); the Export window has no "Final ITSx check" box, export writes no `itsx_check/` folder and no `<group>_itsx_check.csv`, and the summary's `itsx_check` column holds the stored value with NA written as `not_run(not checked)` (Task 18); export reads the copy (Task 19).
  - Existing app ids this task clicks: `mode`, `tax_run_btn`, `tax_lock_btn`, `ref_suggest_btn`, `ref_diff_apply`, `ref_view`, `ref_row_action`, `ref_act_build`, `asmb_run_btn`, `nf_start`, `nf_close`, `asmb_lock_btn`, `ann_details`, `trim_btn`, `trim_mode`, `trim_from`, `trim_to`, `trim_apply`, `annotate_close`, `annot_lock_btn`, `assign_group_btn`, `group_name`, `assign_group_go`, `group_keep_one`, `export_btn`, `export_group_sel`, `export_go`, `export_screen_skip`, `export_overwrite_go`.
- Produces: nothing in the code. A verification report for the user (final message) and screenshots in `/tmp/rp-t24/shots/`.

This task changes no code. It is the end-to-end check of SPEC S12: static checks on the finished branch, the whole suite, mutation checks on the key new tests (S12: "break the code on purpose and confirm the test fails"), upgrades of copies of the real version 2, 3, and 4 test projects (Review Focus 1 and 5, S11), a rebuilt image (S12: pipeline R steps run the package installed in the container), and the app drive S12 names: trim, Update on Annotation Review, Lock, Export, on a fresh test project. It ends with the screenshots to retake (S13). Never touch the user's real projects: every upgrade and app check below uses a copy. Do not fix code here. When a check fails, stop, keep its output, name the task that owns the code (contract section 4), and report to the user; a fix becomes a follow-up commit after the user agrees. Budget: suite 5-15 min, mutations 20-40 min, image 10-30 min, app drive 45-60 min.

- [ ] **Step 1: Confirm the branch state**

Run:

```bash
cd /home/dmacguig/Documents/GitHub/RiboPilot
git branch --show-current
git status --short
git log --oneline a25fd25..HEAD
git ls-files dev/plans/contig-recheck-wip | wc -l
git log --oneline a25fd25..HEAD -- vignettes/Troubleshooting.Rmd vignettes/Your-Own-Project.Rmd
```

Expected: branch `feat/contig-recheck`; `git status --short` shows only ` M vignettes/Troubleshooting.Rmd` and ` M vignettes/Your-Own-Project.Rmd` (the user's own edits; never stage, revert, or overwrite them); the log shows the spec commit and one commit per task 1-23; `0` files of the WIP folder are tracked. The last command lists a commit only if Task 23 touched Troubleshooting.Rmd after asking the user; say which in the report. If any task is missing a commit, stop: this task runs only on the finished branch.

- [ ] **Step 2: Static checks**

Run each command from the repository root and compare with its expectation:

```bash
cd /home/dmacguig/Documents/GitHub/RiboPilot
grep -rn "assembly_backup\|ASSEMBLY_BACKUP_DDL\|ensure_assembly_backup\|chimera_snapshot" R/ inst/ tests/testthat/*.R
```

Expected: every hit is in the step-4 conversion (`R/db_migrate.R` or `R/review_upgrade.R`) or in a test that builds a version-4 database for it (`tests/testthat/test-review-upgrade.R`); `ASSEMBLY_BACKUP_DDL` and `ensure_assembly_backup` give no hit at all (contract section 3, Upgrade).

```bash
grep -rnw "trim_assembly\|restore_assembly\|ensure_assembly_backup\|chimera_original\|itsx_export_check\|run_itsx_check\|itsx_check_command" R/ inst/ man/ NAMESPACE _pkgdown.yml
grep -rnw "trim_assembly\|restore_assembly\|ensure_assembly_backup\|chimera_original\|itsx_export_check\|run_itsx_check" tests/testthat/
grep -rn "There is no second workflow\|There is ONE workflow\|updated again by export\|No Update button here\|Final ITSx check\|final ITSx check\|re-runs curation\|recalculated right away\|in the active run\|trimmed and re-curated\|restored to full length and re-curated" R/ inst/ README.Rmd vignettes/*.Rmd
```

Expected: the first two print nothing (a hit in a comment is a stale comment: report it). The third prints nothing from `R/` or `inst/`; a hit in a vignette is acceptable only when the sentence says export no longer runs ITSx (read each one and report any that still describes the old behavior).

```bash
grep -rPzo 'SELECT\s+\*[^;"]*?FROM\s+review_\w+' R/ inst/ | tr '\0' '\n'
grep -rPzo 'INSERT\s+(OR\s+\w+\s+)?INTO\s+review_\w+\s+(VALUES|SELECT)' R/ inst/ | tr '\0' '\n'
git diff --stat a25fd25 HEAD -- NAMESPACE
git diff --name-only a25fd25 HEAD -- man/ | grep -v -E "^man/(curate_one|compute_ref_alignment|screen_one|parse_blast_candidates|derive_blast_products|chimera_unit|itsx_check_one|export_bundle)\.Rd$"
git ls-files -s inst/nextflow/bin | awk '$1 != "100755"'
git diff --name-only --diff-filter=d a25fd25 HEAD | xargs -r grep -lP '[^\x00-\x7F]'
(cd vignettes && grep -oh "figures/[^)]*\.png" *.Rmd | sort -u | while read f; do [ -f "$f" ] || echo "missing vignettes/$f"; done)
```

Expected: no output from any of them. In order: no `SELECT *` from a review table or view, no positional INSERT into a review table (contract section 1), NAMESPACE unchanged, only the Rd files of the eight changed exported functions changed, every bin file mode 100755, every changed file ASCII, and every figure a vignette names exists. A `missing` line is a screenshot Task 23 referenced but nobody has taken yet: add it to the screenshot list in Step 16 instead of failing.

- [ ] **Step 3: Copy HEAD and upgrade copies of the local test projects**

Make the scratch copy of the committed branch that Steps 3-4 use:

```bash
rm -rf /tmp/rp-t24 && mkdir -p /tmp/rp-t24/head
git -C /home/dmacguig/Documents/GitHub/RiboPilot archive HEAD | tar -x -C /tmp/rp-t24/head
```

Create `/tmp/rp-t24/upgrade.R` with the Write tool:

```r
devtools::load_all("/tmp/rp-t24/head", quiet = TRUE)
src <- path.expand("~/Documents/DJMCompBio/RiboPilot_test")
pre_sql <- "SELECT
  (SELECT COUNT(*) FROM assemblies a JOIN annotate an
     ON an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold) AS units,
  (SELECT COUNT(*) FROM assembly_backup) AS backups,
  (SELECT COUNT(*) FROM annotations WHERE edited = 1) AS edited_features,
  (SELECT COUNT(*) FROM annotate WHERE annotate_lock = 1
     AND NOT (annotate_switch = 2 AND itsx_warnings IS NOT NULL)) AS locked_not_done"
for (p in c("2026_09_30", "2026_10_01", "2026_10_08")) {
  d <- file.path("/tmp/rp-t24/upgrade", p)
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE)
  stopifnot(system2("sqlite3", c("-readonly", shQuote(file.path(src, p, ".sqlite")),
    shQuote(paste0(".backup ", file.path(d, ".sqlite"))))) == 0)
  file.copy(file.path(src, p, ".config"), d)
  con <- open_project_db(file.path(d, ".sqlite"))
  pre <- DBI::dbGetQuery(con, pre_sql)
  DBI::dbDisconnect(con)
  message("== ", p)
  print(pre)
  backwards_compatibility(d)
  con <- open_project_db(file.path(d, ".sqlite"))
  v <- DBI::dbGetQuery(con, "SELECT schema_version FROM schema_history ORDER BY rowid DESC LIMIT 1")$schema_version
  n <- vapply(names(REVIEW_VIEWS_DDL), function(x)
    DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", x))$n, numeric(1))
  print(n)
  ru <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold, edited, check_switch, itsx_warnings FROM review_units")
  print(ru)
  stopifnot(v == 5, n[["review_units_v"]] == pre$units, !DBI::dbExistsTable(con, "assembly_backup"),
    DBI::dbGetQuery(con, "PRAGMA integrity_check")[[1]] == "ok")
  # no edits and no locked unit short of done: S11 gives no copies
  if (pre$backups == 0 && pre$edited_features == 0 && pre$locked_not_done == 0) stopifnot(nrow(ru) == 0)
  DBI::dbDisconnect(con)
}
```

Run: `Rscript /tmp/rp-t24/upgrade.R`
Expected: no error. 2026_09_30 (version 2) upgrades through versions 3, 4, and 5; 2026_10_01 (version 3) through 4 and 5; 2026_10_08 (version 4) to 5, each printing the step-4 note "review copies for edited contigs; ITSx origin and scores". On 2026-10-09 the pre-counts were `units` 0, 13, and 18 with no backups, no edited features, and no locked unit short of done, so every `review_units` is empty and `review_units_v` has 0, 13, and 18 rows. If a project's pre-counts are non-zero now, check its copies against S11.3-11.4 by hand and report them. The originals under `~/Documents/DJMCompBio/RiboPilot_test/` are only read (`sqlite3 -readonly`). The upgraded `2026_10_08` copy is the `RIBOPILOT_TEST_PROJECT` of Step 4, and `2026_10_01` is opened in the app in Step 14.

- [ ] **Step 4: Run the whole suite**

Create `/tmp/rp-t24/suite.R` with the Write tool:

```r
devtools::load_all("/tmp/rp-t24/head", quiet = TRUE)
r <- as.data.frame(testthat::test_dir("/tmp/rp-t24/head/tests/testthat", reporter = "summary",
  stop_on_failure = FALSE))
saveRDS(r, "/tmp/rp-t24/suite.rds")
bad <- r[r$failed > 0 | r$error, c("file", "test", "failed", "error")]
sk <- r[r$skipped, c("file", "test")]
cat(sprintf("\n%d tests: %d failing, %d skipped, %.0f s\n", nrow(r), nrow(bad), nrow(sk), sum(r$real)))
print(bad, row.names = FALSE)
print(sk, row.names = FALSE)
```

Run it in the background (Bash `run_in_background: true`; it takes 5-15 minutes) from a shell whose PATH has no conda environment with barrnap or pybarrnap (check `which barrnap pybarrnap` prints nothing):

```bash
cd /tmp/rp-t24/head && NOT_CRAN=true RIBOPILOT_TEST_PROJECT=/tmp/rp-t24/upgrade/2026_10_08 Rscript /tmp/rp-t24/suite.R > /tmp/rp-t24/suite.log 2>&1; tail -80 /tmp/rp-t24/suite.log
```

Expected: `0 failing`. The suite runs on the archive copy so the user's uncommitted vignette edits are not tested as part of the branch; `NOT_CRAN=true` runs the real-Nextflow tests (test-nf-main-run.R, test-nf-recheck-run.R) and the env-gated real stock-ITSx test on the upgraded copy (it needs Docker and the `ribopilot:devel` image that copy's `.config` names). Acceptable skips: the chimera fixture test (`blast not installed`, `no local panel`, or `no local taxdump`), the GBIF and GEOME API tests when offline, Windows-only and pandoc skips. A skip for `RIBOPILOT_TEST_PROJECT not set`, `Nextflow not installed`, or `On CRAN` means the environment above was not applied: fix it and rerun. Known non-branch failures: test-rrna-kingdom.R fails 4 expectations when a real barrnap is on PATH (rerun that file with a clean PATH); a GBIF/GEOME API test can fail on a network hiccup. For any failure you believe is not this branch's, rerun that one file on the base commit and report both results. The function below takes the test file name; the example reruns test-gbif-api.R:

```bash
base_rerun() {
  [ -d /tmp/rp-t24/base/R ] || { mkdir -p /tmp/rp-t24/base && git -C /home/dmacguig/Documents/GitHub/RiboPilot archive a25fd25 | tar -x -C /tmp/rp-t24/base; }
  (cd /tmp/rp-t24/base && NOT_CRAN=true Rscript -e "devtools::load_all('.'); testthat::test_file('tests/testthat/$1')")
}
base_rerun test-gbif-api.R
```

- [ ] **Step 5: Mutation checks that apply themselves**

Create `/tmp/rp-t24/mut.R` with the Write tool:

```r
# Mutation checks on a scratch copy of the committed branch. Source this file,
# then call baseline(), mutate(), check(), and restore().
REPO <- Sys.getenv("RP_REPO", "/home/dmacguig/Documents/GitHub/RiboPilot")
WORK <- Sys.getenv("RP_MUT", "/tmp/rp-t24/mut")
if (!dir.exists(file.path(WORK, "R"))) {
  dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
  system(sprintf("git -C %s archive HEAD | tar -x -C %s", shQuote(REPO), shQuote(WORK)))
}

# Failing expectations plus errors over the test files; NA when no listed file
# exists or the run crashed (the output is in last_run.log).
run_tests <- function(tests, not_cran = FALSE) {
  have <- file.exists(file.path(WORK, "tests/testthat", tests))
  if (any(!have)) message("  no such test file: ", paste(tests[!have], collapse = ", "))
  tests <- tests[have]
  if (!length(tests)) return(NA_integer_)
  out <- tempfile(fileext = ".rds")
  code <- sprintf(paste0("devtools::load_all(%s, quiet = TRUE); ",
    "r <- do.call(rbind, lapply(%s, function(f) as.data.frame(testthat::test_file(",
    "file.path(%s, 'tests/testthat', f), reporter = 'silent', stop_on_failure = FALSE)))); ",
    "saveRDS(r, %s)"),
    deparse(WORK), paste(deparse(tests), collapse = ""), deparse(WORK), deparse(out))
  log <- file.path(WORK, "last_run.log")
  system2("Rscript", c("-e", shQuote(code)), stdout = log, stderr = log,
    env = if (not_cran) "NOT_CRAN=true" else character(0))
  if (!file.exists(out)) return(NA_integer_)
  r <- readRDS(out)
  as.integer(sum(r$failed) + sum(r$error))
}

baseline <- function(tests, not_cran = FALSE) {
  bad <- run_tests(tests, not_cran)
  message(sprintf("baseline %s: %s failing", paste(tests, collapse = ", "), bad))
  invisible(bad)
}

verdict <- function(name, bad, n = NA) {
  res <- if (is.na(bad)) "CRASHED (see last_run.log)" else if (bad > 0) "CAUGHT" else "MISSED"
  message(sprintf("%-20s %s%s", name, res, if (is.na(n)) "" else sprintf(" (%d edit(s))", n)))
  invisible(res)
}

# Replace every fixed-string `from` with `to` in one file, run the tests, put
# the file back.
mutate <- function(name, file, from, to, tests, not_cran = FALSE) {
  f <- file.path(WORK, file)
  if (!file.exists(f)) {
    message(sprintf("%-20s NOT FOUND (no file %s)", name, file))
    return(invisible("NOT FOUND"))
  }
  orig <- readChar(f, file.size(f), useBytes = TRUE)
  n <- lengths(regmatches(orig, gregexpr(from, orig, fixed = TRUE)))
  if (n == 0) {
    message(sprintf("%-20s NOT FOUND in %s", name, file))
    return(invisible("NOT FOUND"))
  }
  on.exit(writeChar(orig, f, eos = NULL, useBytes = TRUE))
  writeChar(gsub(from, to, orig, fixed = TRUE), f, eos = NULL, useBytes = TRUE)
  verdict(name, run_tests(tests, not_cran), n)
}

# After a hand edit of a file under WORK: run the tests, then call restore().
check <- function(name, tests, not_cran = FALSE) verdict(name, run_tests(tests, not_cran))

restore <- function(file) {
  system(sprintf("git -C %s show HEAD:%s > %s", shQuote(REPO), shQuote(file),
    shQuote(file.path(WORK, file))))
  invisible(TRUE)
}
```

Create `/tmp/rp-t24/mutations.R` with the Write tool. Each entry breaks one rule of the SPEC with the exact text the contract or an earlier task fixes, and names the test file that must catch it:

```r
source(file.path(Sys.getenv("RP_T24", "/tmp/rp-t24"), "mut.R"))
mt <- function(name, file, from, to, tests, not_cran = FALSE)
  list(name = name, file = file, from = from, to = to, tests = tests, not_cran = not_cran)

M <- list(
  mt("views-fallback", "R/db_schema.R", 'WHERE NOT ", REVIEW_COPY_SQL)', 'WHERE ", REVIEW_COPY_SQL)',
    "test-review-schema.R"),
  mt("done-rule", "R/db_schema.R",
    "CASE WHEN an.annotate_switch = 2 AND an.itsx_warnings IS NULL THEN 1 ELSE an.annotate_switch END",
    "an.annotate_switch", "test-review-schema.R"),
  mt("invalidate-rev", "R/review.R", "rev = rev + 1", "rev = rev",
    c("test-review-copy.R", "test-recheck.R")),
  mt("guard-rev", "R/recheck.R", "AND rev = CAST(? AS INTEGER)",
    "AND CAST(? AS INTEGER) IS NOT NULL", "test-recheck.R"),
  mt("guard-state", "R/recheck.R", "AND check_switch IN (1, 4)", "AND check_switch IN (0, 1, 2, 3, 4)",
    "test-recheck.R"),
  mt("blast-constants", "R/recheck.R", "max_hsps = 5L", "max_hsps = 4L", "test-recheck.R"),
  mt("queue-copy", "inst/nextflow/modules/recheck.nf", "COALESCE(an.annotate_lock,0) = 0",
    "COALESCE(an.annotate_lock,0) >= 0", c("test-nextflow.R", "test-nf-pending.R")),
  mt("recheck-run", "inst/nextflow/bin/rp_recheck.R", "RiboPilot:::recheck_unit(",
    "RiboPilot:::recheck_unit_gone(", "test-nf-recheck-run.R", not_cran = TRUE),
  mt("tools-exit", "inst/nextflow/bin/rp_rrna_tools.sh", 'exit "$b"', "exit 0", "test-rrna-tools.R"),
  mt("itsx-error-line", "R/itsx_check.R", 'grepl("ERROR|FATAL"', 'grepl("NO_SUCH_LINE"',
    "test-itsx-check.R"),
  mt("curate-replace", "R/curate_rrna.R",
    'partials <- if (review) "partial_start = ?, partial_stop = ?" else',
    'partials <- if (FALSE) "partial_start = ?, partial_stop = ?" else', "test-curate.R"),
  mt("release-running", "R/app_server.R", "recheck_release_running(con)", "invisible(con)",
    "test-app-run.R"),
  mt("waiting-count", "R/nf_command.R", "AND NOT EXISTS (SELECT 1 FROM review_units r",
    "AND EXISTS (SELECT 1 FROM review_units r", "test-app-run.R"),
  mt("lineage-rev", "R/tax_search.R", "rev = rev + 1", "rev = rev",
    c("test-review-follow.R", "test-taxonomy-assign.R")),
  mt("export-static", "R/export.R", "review_annotations_v", "annotations",
    c("test-export-views.R", "test-export-review.R")),
  mt("export-not-run", "R/export.R", "not_run(not checked)", "not_run(nope)",
    c("test-export.R", "test-export-units.R", "test-export-review.R", "test-itsx-check.R"))
)

for (nc in c(FALSE, TRUE)) {
  fs <- unique(unlist(lapply(Filter(function(m) m$not_cran == nc, M), `[[`, "tests")))
  for (f in fs) baseline(f, not_cran = nc)
}
res <- vapply(M, function(m) mutate(m$name, m$file, m$from, m$to, m$tests, m$not_cran), character(1))
print(data.frame(mutation = vapply(M, `[[`, "", "name"), result = res), row.names = FALSE)
```

What each one breaks: views-fallback, a contig with a copy also returns its original rows (S4.3); done-rule, a done original without a stored ITSx value reads as done (S6.1); invalidate-rev and guard-rev, an edit during a run no longer discards the result (S5.4); guard-state, On hold no longer discards it (S5.4); blast-constants, the recheck BLAST drifts from blast.nf (S7.3); queue-copy, query_recheck and `nf_recheck_pending` select differently (S7.4); recheck-run, the real `main_recheck.nf` run (S12); tools-exit, a barrnap failure no longer fails `annotate_tools` (Review Focus 3); itsx-error-line, a crashed ITSx reads as clean (S6.3); curate-replace, curation keeps a stale partial on a copy (S4.4); release-running, contigs left at 4 stay running after an in-app run (S5.4); waiting-count, "Nothing to update" counts the wrong contigs (S7.4); lineage-rev, a taxonomy change no longer requeues queued copies (S5.3); export-static, export reads an original table (S9, S12 static check); export-not-run, the summary's NA legend (S9).

Run it in the background (20-40 minutes): `Rscript /tmp/rp-t24/mutations.R > /tmp/rp-t24/mutations.log 2>&1; tail -40 /tmp/rp-t24/mutations.log`
Expected: every `baseline` line says `0 failing`, and every mutation says `CAUGHT`. A baseline that is not 0 makes that file's verdicts meaningless: report it (Step 4 should already have shown the failure). For the other results:
- `NOT FOUND`: the code differs from the text quoted here. Open that file in `/tmp/rp-t24/mut`, find the line that implements the rule named above, break it the same way by hand with the Edit tool (only under `/tmp/rp-t24/mut`), then run the check, restore, and compare command the way Step 6 does, filled in with that entry's name, test files, and file (add `not_cran = TRUE` for recheck-run). For example, for guard-rev: `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("guard-rev", "test-recheck.R"); restore("R/recheck.R")' && cmp /tmp/rp-t24/mut/R/recheck.R /home/dmacguig/Documents/GitHub/RiboPilot/R/recheck.R`. If no such line exists, report the rule as untested.
- `no such test file`: the owning task put the test elsewhere. Find it by grepping `/tmp/rp-t24/mut/tests/testthat` for the function the rule lives in (for example `grep -rln "recheck_unit" /tmp/rp-t24/mut/tests/testthat`) and rerun the mutation by hand with that file.
- `MISSED`: no test pins the rule. Do not add tests here; report the rule, the file, and the test file it should live in.
- `CRASHED`: read `/tmp/rp-t24/mut/last_run.log`; a mutation must leave the code loadable, so a crash means the replacement text broke the parse: redo it by hand.

- [ ] **Step 6: Mutation checks made by hand**

These rules have no fixed text to replace. For each one: make the break with the Edit tool in the named file under `/tmp/rp-t24/mut` (never in the repository), run its command (which checks, restores the file, and compares it with the repository; `cmp` must print nothing), and record the verdict before the next one.

1. trim-orientation (Review Focus 2). `R/review.R`, `review_trim()`: make a left cut on a `-` strand feature set `partial_start_user` instead of `partial_stop_user` (the `+` rule on both strands).
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("trim-orientation", "test-review-edits.R"); restore("R/review.R")' && cmp /tmp/rp-t24/mut/R/review.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review.R`
2. copy-leftovers (S5.1, S5.4). `R/review.R`, `review_copy()`: delete the statement that removes the unit's rows from the seven mirror tables before copying.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("copy-leftovers", "test-review-copy.R"); restore("R/review.R")' && cmp /tmp/rp-t24/mut/R/review.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review.R`
3. state-waiting (S5.2). `R/review.R`, `review_set_state()`: drop the "press Update on the Assemble step first" refusal.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("state-waiting", "test-review-edits.R"); restore("R/review.R")' && cmp /tmp/rp-t24/mut/R/review.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review.R`
4. follow-locked (S5.3). `R/review.R`, `review_follow_original()`: treat a locked copy as unlocked (skip the `annotate_lock` branch).
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("follow-locked", "test-review-follow.R"); restore("R/review.R")' && cmp /tmp/rp-t24/mut/R/review.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review.R`
5. lock-gate (D2, S6.2). `R/review.R`, `review_lock_gate()`: return NA (lockable) for a copy whose state is not run.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("lock-gate", "test-lock-gate.R"); restore("R/review.R")' && cmp /tmp/rp-t24/mut/R/review.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review.R`
6. map-no-hits (S7.3 step 2). `R/recheck.R`, `recheck_hits()`: in map mode with no original hits, also set the copy's `blast_accession`, `blast_species`, `blast_pident`, `blast_qcovs`, `blast_evalue`, and `blast_source` to NULL.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("map-no-hits", "test-recheck.R"); restore("R/recheck.R")' && cmp /tmp/rp-t24/mut/R/recheck.R /home/dmacguig/Documents/GitHub/RiboPilot/R/recheck.R`
7. evidence-tables (S5.3, re-assembly row). `R/pipeline_steps.R`, `UNIT_EVIDENCE_TABLES`: remove the eight review tables from the vector.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("evidence-tables", c("test-review-copy.R", "test-load-assemblies.R")); restore("R/pipeline_steps.R")' && cmp /tmp/rp-t24/mut/R/pipeline_steps.R /home/dmacguig/Documents/GitHub/RiboPilot/R/pipeline_steps.R`
8. export-backstop (D2, S6.2). `R/export.R`, `export_bundle()`: remove the filter that drops units with `done == 0`.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("export-backstop", "test-export-review.R"); restore("R/export.R")' && cmp /tmp/rp-t24/mut/R/export.R /home/dmacguig/Documents/GitHub/RiboPilot/R/export.R`
9. upgrade-unit (Review Focus 1). The step-4 conversion (`R/review_upgrade.R`, or `R/db_migrate.R` if Task 22 kept it there): in the handler that records one unit's error for the step's message, rethrow the error instead (`stop(e)`). With the conversion in `R/db_migrate.R`, use that path in both places of the command.
   `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("upgrade-unit", "test-review-upgrade.R"); restore("R/review_upgrade.R")' && cmp /tmp/rp-t24/mut/R/review_upgrade.R /home/dmacguig/Documents/GitHub/RiboPilot/R/review_upgrade.R`
10. recheck-tools-exit (Review Focus 3). `inst/nextflow/modules/recheck.nf`, process `recheck_tools`: exit with `rc` for every non-zero `rc` (ignore the out-of-memory code list).
    `Rscript -e 'source("/tmp/rp-t24/mut.R"); check("recheck-tools-exit", c("test-nextflow.R", "test-nf-recheck-run.R"), not_cran = TRUE); restore("inst/nextflow/modules/recheck.nf")' && cmp /tmp/rp-t24/mut/inst/nextflow/modules/recheck.nf /home/dmacguig/Documents/GitHub/RiboPilot/inst/nextflow/modules/recheck.nf`

Expected: `CAUGHT` for each. Handle `MISSED`, a missing test file, or a rule you cannot find as in Step 5. Afterwards confirm that every mutation was undone; this prints nothing:

```bash
cd /tmp/rp-t24/mut && git -C /home/dmacguig/Documents/GitHub/RiboPilot ls-files R inst | while read f; do cmp -s "$f" "/home/dmacguig/Documents/GitHub/RiboPilot/$f" || echo "differs: $f"; done
```

- [ ] **Step 7: Rebuild the image**

Pipeline R steps in Docker load the RiboPilot installed in the image, while the bin scripts come from the host; the branch's bin scripts call new internals (`RiboPilot:::review_follow_original`, `RiboPilot:::recheck_unit`, `RiboPilot:::itsx_store_not_run`), so no older image can run this branch. Build a separate local tag from the repository root, in the background (10-30 minutes; it reuses the panel tarballs already staged in `docker/` and replaces only the ignored `docker/RiboPilot_*.tar.gz`):

```bash
cd /home/dmacguig/Documents/GitHub/RiboPilot && docker/deploy-local.sh recheck > /tmp/rp-t24/image.log 2>&1; tail -3 /tmp/rp-t24/image.log
```

Expected: `Local image build successful: ribopilot:recheck (rDNA panel 2026.08.18 baked in)`. Then:

```bash
docker run --rm ribopilot:recheck Rscript -e 'ns <- asNamespace("RiboPilot"); f <- c("recheck_unit", "review_follow_original", "itsx_run_finished", "itsx_store_not_run", "review_copy"); stopifnot(all(f %in% ls(ns, all.names = TRUE))); cat("image ok\n")'
git -C /home/dmacguig/Documents/GitHub/RiboPilot status --short
```

Expected: `image ok`; git status unchanged from Step 1.

- [ ] **Step 8: Create a fresh test project**

`~/Documents/DJMCompBio/RiboPilot_test/recheck_t24` must not exist yet; if it does (an earlier attempt), ask the user before deleting it. Create `/tmp/rp-t24/d0_project.R` with the Write tool:

```r
devtools::load_all("/home/dmacguig/Documents/GitHub/RiboPilot", quiet = TRUE)
proj <- path.expand("~/Documents/DJMCompBio/RiboPilot_test/recheck_t24")
stopifnot(!file.exists(file.path(proj, ".sqlite")))
dir.create(file.path(proj, "data"), recursive = TRUE, showWarnings = FALSE)
# the filtered reads of an earlier test project spare the release download
src <- path.expand(file.path("~/Documents/DJMCompBio/RiboPilot_test/2026_10_08/data",
  c("SRR22396758_R1.fastq.gz", "SRR22396758_R2.fastq.gz")))
if (all(file.exists(src))) file.copy(src, file.path(proj, "data"))
new_test_project(proj, n = 1, container = "ribopilot:recheck",
  ref_dir = path.expand("~/ribopilot_ref_dbs"), Rproj = FALSE, fetch_metadata = FALSE)
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
s <- DBI::dbGetQuery(con, "SELECT ID FROM samples ORDER BY ID")$ID
v <- DBI::dbGetQuery(con, "SELECT schema_version FROM schema_history ORDER BY rowid DESC LIMIT 1")$schema_version
DBI::dbDisconnect(con)
print(s)
stopifnot(identical(s, c("CHIM00", "CHIM01", "SRR22396758")), v == 5)
cat(grep("container", readLines(file.path(proj, ".config")), value = TRUE), sep = "\n")
```

Run: `Rscript /tmp/rp-t24/d0_project.R`
Expected: the three samples (one fish genome skim, SRR22396758 *Upeneus parvus*, plus the bundled CHIM00 and CHIM01, all one reference group), schema version 5, and a `.config` container line naming `ribopilot:recheck`. `ref_dir` points at `~/ribopilot_ref_dbs`, where the Syngnathiformes database of the earlier test projects already sits, so the Reference step wires it instead of downloading. If the reads were not copied, the download needs `GITHUB_PAT` (see the error text).

- [ ] **Step 9: Start the app and a headless browser**

Check that ports 8791 and 9222 are free: `ss -ltn | grep -E ':(8791|9222)\b'` prints nothing. (If one is taken, pick others, and prefix every `Rscript /tmp/rp-t24/d*.R` call below with `RP_APP_URL=http://127.0.0.1:<app port> RP_CDP_PORT=<chrome port>`.)

Create `/tmp/rp-t24/app.R` with the Write tool:

```r
a <- commandArgs(TRUE)
devtools::load_all(a[1], quiet = TRUE)
shiny::runApp(RiboPilot(a[2], port = as.integer(a[3]), launch.browser = FALSE, headless = FALSE))
```

Create `/tmp/rp-t24/cdp.R` with the Write tool. It attaches every call to the one browser tab, so the app keeps a single Shiny session (an in-app run stops when its session ends):

```r
# Drive the app tab of the headless Chrome listening on port 9222. Every call
# attaches to that one tab, so the app keeps a single Shiny session.
CDP_PORT <- as.integer(Sys.getenv("RP_CDP_PORT", "9222"))
APP_URL <- Sys.getenv("RP_APP_URL", "http://127.0.0.1:8791")
SHOTS <- Sys.getenv("RP_SHOTS", "/tmp/rp-t24/shots")
PROJ <- path.expand(Sys.getenv("RP_PROJ", "~/Documents/DJMCompBio/RiboPilot_test/recheck_t24"))

tab <- local({
  b <- NULL
  function() {
    if (!is.null(b)) return(b)
    cm <- chromote::Chromote$new(browser = chromote::ChromeRemote$new(host = "127.0.0.1", port = CDP_PORT))
    tg <- cm$Target$getTargets()$targetInfos
    pg <- Filter(function(t) identical(t$type, "page") && startsWith(t$url, APP_URL), tg)
    if (!length(pg)) stop("No Chrome tab on ", APP_URL)
    b <<- chromote::ChromoteSession$new(parent = cm, targetId = pg[[1]]$targetId)
    b
  }
})

js <- function(expr) {
  r <- tab()$Runtime$evaluate(expr, returnByValue = TRUE, awaitPromise = TRUE)
  if (!is.null(r$exceptionDetails)) stop("JS error: ", r$exceptionDetails$text)
  r$result$value
}

# Wait until Shiny is idle.
idle <- function(timeout = 60) {
  t0 <- Sys.time()
  Sys.sleep(0.5)
  while (isTRUE(js("document.documentElement.classList.contains('shiny-busy')"))) {
    if (difftime(Sys.time(), t0, units = "secs") > timeout) stop("app still busy after ", timeout, " s")
    Sys.sleep(0.5)
  }
  Sys.sleep(0.5)
  invisible(TRUE)
}

# Click a visible, enabled element by id.
click <- function(id, timeout = 60) {
  ok <- js(sprintf("(() => { const e = document.getElementById('%s');
    if (!e || e.disabled || e.offsetParent === null) return false; e.click(); return true; })()", id))
  if (!isTRUE(ok)) stop("No visible, enabled element #", id)
  idle(timeout)
}

# Click the last visible, enabled button whose text is exactly `label`.
click_text <- function(label, timeout = 60) {
  ok <- js(sprintf("(() => { const want = %s;
    const els = [...document.querySelectorAll('button, a.btn, .swal2-confirm, .swal2-cancel')]
      .filter(e => e.offsetParent !== null && !e.disabled && e.textContent.trim() === want);
    if (!els.length) return false; els[els.length - 1].click(); return true; })()",
    jsonlite::toJSON(label, auto_unbox = TRUE)))
  if (!isTRUE(ok)) stop("No visible, enabled button labelled '", label, "'")
  idle(timeout)
}

# Set an input's value on the server (text, numeric, radio, picker, row index).
set_input <- function(id, value) {
  js(sprintf("Shiny.setInputValue('%s', %s, {priority: 'event'})", id,
    jsonlite::toJSON(value, auto_unbox = TRUE)))
  idle()
}

# Click one option of a radio group (the step bar, the Reference view switch).
pick <- function(name, value) {
  ok <- js(sprintf("(() => { const e = document.querySelector('input[name=\"%s\"][value=\"%s\"]');
    if (!e) return false; e.click(); return true; })()", name, value))
  if (!isTRUE(ok)) stop("No option ", value, " in ", name)
  idle()
}
set_mode <- function(mode) pick("mode", mode)

# The main table's data (every column, hidden ones too), as character.
rows <- function(table = "tbl") {
  d <- js(sprintf("Reactable.getState('%s').data", table))
  if (!length(d)) return(data.frame())
  keys <- unique(unlist(lapply(d, names)))
  as.data.frame(lapply(stats::setNames(keys, keys), function(k)
    vapply(d, function(r) { v <- r[[k]]; if (is.null(v) || length(v) != 1) NA_character_ else as.character(v) },
      character(1))), stringsAsFactors = FALSE)
}

# Tick exactly these rows (1-based indices into rows()).
select_rows <- function(i, table = "tbl") {
  js(sprintf("Reactable.getInstance('%s').setRowsSelected(%s)", table,
    jsonlite::toJSON(as.character(i - 1L))))
  idle()
}

shot <- function(name) {
  dir.create(SHOTS, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(SHOTS, paste0(name, ".png"))
  tab()$screenshot(f, show = FALSE)
  message("saved ", f)
  invisible(f)
}

modal_text <- function() {
  js("(() => [...document.querySelectorAll('.modal.show .modal-content, .modal.in .modal-content, .swal2-popup')]
    .filter(e => e.offsetParent !== null).map(e => e.innerText).join('\\n----\\n'))()")
}

page_text <- function() js("document.body.innerText")

# Poll the page until its text matches `pattern`; FALSE on timeout.
wait_text <- function(pattern, seconds = 300) {
  t0 <- Sys.time()
  repeat {
    if (grepl(pattern, page_text())) return(invisible(TRUE))
    if (difftime(Sys.time(), t0, units = "secs") > seconds) return(invisible(FALSE))
    Sys.sleep(2)
  }
}

# Press the confirm button of an open alert, if there is one.
close_alert <- function() {
  js("(() => { const e = document.querySelector('.swal2-confirm');
    if (!e || e.offsetParent === null) return false; e.click(); return true; })()")
  idle()
}

# Read-only query on the project database.
q <- function(sql, ...) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(PROJ, ".sqlite"), flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, sql, params = if (...length()) list(...) else NULL)
}

reports <- function(workflow) {
  list.files(file.path(PROJ, ".runs", "reports"), paste0("^", workflow, "_.*_report\\.txt$"))
}

# Press a run window's Run from App, remembering which reports existed before.
start_run <- function(id, workflow) {
  saveRDS(reports(workflow), file.path(dirname(SHOTS), paste0(workflow, "_before.rds")))
  click(id)
}

# Print the report of the run start_run() began; NULL when it is not done
# within `minutes` (call again: the run goes on in the app).
wait_run <- function(workflow, minutes = 9) {
  before <- readRDS(file.path(dirname(SHOTS), paste0(workflow, "_before.rds")))
  t0 <- Sys.time()
  repeat {
    new <- setdiff(reports(workflow), before)
    if (length(new)) {
      f <- file.path(PROJ, ".runs", "reports", new[1])
      cat(readLines(f), sep = "\n")
      return(invisible(f))
    }
    if (difftime(Sys.time(), t0, units = "mins") > minutes) {
      message("not finished yet")
      return(invisible(NULL))
    }
    Sys.sleep(15)
  }
}

# Reload the tab (after the app was restarted on the same port).
reload <- function() {
  tab()$Page$reload()
  Sys.sleep(5)
  idle()
}

# A contig's original rows, to prove that edits and the recheck never touch them.
orig_snap <- function(u) {
  k <- function(sql) q(sql, u$ID, u$path, u$scaffold)
  w <- " WHERE ID = ? AND path = ? AND scaffold = ?"
  list(
    asm = k(paste0("SELECT sequence, length, blast_accession, blast_species FROM assemblies", w)),
    ann = k(paste0("SELECT element, start, stop, strand, partial_start, partial_stop, status, flags,
      model_coverage FROM annotations", w, " ORDER BY element")),
    an = k(paste0("SELECT n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings,
      itsx_origin, itsx_scores, annotate_switch FROM annotate", w)),
    scr = k(paste0("SELECT screen_verdict, screen_reason FROM screen", w)),
    chi = k(paste0("SELECT label, score FROM chimera", w)),
    hits = k(paste0("SELECT COUNT(*) AS n FROM blast_hits", w)))
}
```

Start the app (its own process, from the repository source; pipeline runs launch from the project folder), wait until it answers, then start Chrome with remote debugging on the app's URL. Keep both PIDs; never stop them with `pkill -f` (the pattern also matches your own shell):

```bash
P=$HOME/Documents/DJMCompBio/RiboPilot_test/recheck_t24
(nohup Rscript /tmp/rp-t24/app.R /home/dmacguig/Documents/GitHub/RiboPilot "$P" 8791 > /tmp/rp-t24/app.log 2>&1 & echo $! > /tmp/rp-t24/app.pid)
for i in $(seq 60); do curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8791 | grep -q 200 && break; sleep 2; done
CHROME=/home/dmacguig/.local/share/flatpak/app/com.google.Chrome/current/active/files/extra/chrome
mkdir -p /dev/shm/rp_t24_chrome
(nohup "$CHROME" --headless=new --disable-gpu --no-first-run --no-default-browser-check --remote-debugging-port=9222 --user-data-dir=/dev/shm/rp_t24_chrome --window-size=1600,1000 http://127.0.0.1:8791 > /tmp/rp-t24/chrome.log 2>&1 & echo $! > /tmp/rp-t24/chrome.pid)
sleep 8
Rscript -e 'source("/tmp/rp-t24/cdp.R"); idle(); shot("00_app_open")'
```

Expected: `saved /tmp/rp-t24/shots/00_app_open.png`; read it: the app on the Taxonomy step with three samples. The app stops by itself when its last browser session closes, so keep Chrome running until Step 15. If chromote fails on a socket path, rerun with `TMPDIR=/dev/shm` before `Rscript`.

How to drive: each step below is an R file run with `Rscript /tmp/rp-t24/<file>.R`; each `stopifnot()` is a check, and each `shot()` saves a screenshot you read with the Read tool when a check fails or something looks off. If an id has moved (a click stops with `No visible, enabled element #<id>`), look the id up in `R/app_ui.R`, `R/app_modals.R`, or `R/app_server.R`, use the new one, and note it in the report. If an unexpected confirm appears, read `modal_text()` and take the choice that continues the flow (for example Lock anyway, Export anyway). If a step cannot be scripted at all, ask the user to do that one step in their own browser at http://127.0.0.1:8791 while the headless tab stays open, then continue.

- [ ] **Step 10: Taxonomy, Reference, and Assemble in the app**

Create these four files with the Write tool.

`/tmp/rp-t24/d1_taxonomy.R`:

```r
source("/tmp/rp-t24/cdp.R")
set_mode("Taxonomy")
click("tax_run_btn")
stopifnot(wait_text("Taxonomy assigned|Some samples need a taxon|Taxonomy assignment failed", 540))
cat(modal_text(), "\n")
shot("01_taxonomy_assigned")
close_alert()
t <- q("SELECT ID, taxid, sci_name, match_type FROM samples ORDER BY ID")
print(t)
stopifnot(!anyNA(t$taxid))
select_rows(seq_len(nrow(rows())))
click("tax_lock_btn")
if (grepl("Lock anyway", modal_text())) click_text("Lock anyway")
stopifnot(all(q("SELECT taxonomy_lock FROM samples")$taxonomy_lock == 1))
```

`/tmp/rp-t24/d2_reference.R`:

```r
source("/tmp/rp-t24/cdp.R")
set_mode("Reference")
click("ref_suggest_btn")
stopifnot(wait_text("Proposed reference groups|Suggestion failed", 300))
cat(modal_text(), "\n")
click("ref_diff_apply")
pick("ref_view", "Groups")
g <- rows()
print(g)
select_rows(seq_len(nrow(g)))
set_input("ref_row_action", 1L)
click("ref_act_build")
m <- modal_text()
cat(m, "\n")
if (grepl("Run from App", m)) start_run("ref_nf_start", "reference")
shot("02_reference")
print(q("SELECT ref_id, reference_switch, source, seeds_db FROM reference"))
```

`/tmp/rp-t24/d3_assemble_start.R`:

```r
source("/tmp/rp-t24/cdp.R")
set_mode("Assemble")
click("asmb_run_btn")
cat(modal_text(), "\n")
close_alert()
stopifnot(grepl("Run from App", modal_text()))
shot("03_update_window")
start_run("nf_start", "pipeline")
```

`/tmp/rp-t24/d4_assemble_wait.R`:

```r
source("/tmp/rp-t24/cdp.R")
if (is.null(wait_run("pipeline"))) quit(status = 3)
shot("04_pipeline_report")
close_alert()
click("nf_close")
print(q("SELECT ID, assemble_switch FROM assemble ORDER BY ID"))
a <- q("SELECT ID, path, scaffold, annotate_switch, n_features, itsx_warnings, itsx_origin
  FROM annotate ORDER BY ID, path, scaffold")
print(a)
stopifnot(all(q("SELECT assemble_switch FROM assemble")$assemble_switch == 2))
# S6.3: after the main pipeline NULL only means never checked
stopifnot(!anyNA(a$itsx_warnings[!is.na(a$n_features)]))
stopifnot(q("SELECT COUNT(*) AS n FROM review_units")$n == 0)
set_mode("Assemble")
r <- rows()
select_rows(seq_len(nrow(r)))
click("asmb_lock_btn")
if (grepl("Lock anyway", modal_text())) click_text("Lock anyway")
stopifnot(all(q("SELECT assemble_lock FROM assemble")$assemble_lock == 1))
```

Run, in order: `Rscript /tmp/rp-t24/d1_taxonomy.R` (a few minutes: the first Assign Taxonomy downloads and indexes the NCBI taxonomy), `Rscript /tmp/rp-t24/d2_reference.R`, `Rscript /tmp/rp-t24/d3_assemble_start.R`, then `Rscript /tmp/rp-t24/d4_assemble_wait.R` with a 10-minute Bash timeout, repeated while it exits with status 3 (`not finished yet`; the run goes on in the app).
Expected: d1, all three samples get a taxid and are locked. d2, one group (Syngnathiformes) that ends with `reference_switch` 2, `source` `existing`, and a `seeds_db` under `~/ribopilot_ref_dbs` (a toast says the group was wired). If d2 printed `Run from App` instead, the database had to be built: run `Rscript -e 'source("/tmp/rp-t24/cdp.R"); if (!is.null(wait_run("reference"))) { close_alert(); click("ref_nf_close") }'` until it prints the report. d3, the Update window, the pipeline run started. d4, the pipeline report (about 5-10 minutes), every sample at `assemble_switch` 2, every annotated contig with a stored ITSx value (Task 5), no review copies (the main pipeline never makes one), and all three samples locked (CHIM01 asks first; the script answers Lock anyway). If d4 stops after its report was found, finish its remaining lines by hand instead of rerunning it.

- [ ] **Step 11: Trim on Annotation Review**

Create `/tmp/rp-t24/d5_trim.R` with the Write tool. It cuts 50 bp into the first gene of SRR22396758's five-feature contig and keeps the right end exactly at the last gene end:

```r
source("/tmp/rp-t24/cdp.R")
set_mode("Annotate")
r <- rows()
print(r[, intersect(c("ID", "path", "scaffold", "annotate_switch", "check_state", "n_features", "length"), names(r))])
i <- which(r$ID == "SRR22396758" & r$n_features == "5")[1]
stopifnot(!is.na(i))
u <- list(ID = r$ID[i], path = as.integer(r$path[i]), scaffold = as.integer(r$scaffold[i]))
f <- q("SELECT element, start, stop, strand FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?
  ORDER BY start", u$ID, u$path, u$scaffold)
u$first <- f[1, ]
u$from <- f$start[1] + 50L
u$to <- max(f$stop)
saveRDS(u, "/tmp/rp-t24/unit.rds")
saveRDS(orig_snap(u), "/tmp/rp-t24/orig.rds")
stopifnot(q("SELECT COUNT(*) AS n FROM review_units")$n == 0)

# cut 50 bp into the first gene, keep the right end exactly at the last gene end
set_input("ann_details", i)
stopifnot(grepl(u$ID, modal_text(), fixed = TRUE))
click("trim_btn")
set_input("trim_mode", "custom")
set_input("trim_from", u$from)
set_input("trim_to", u$to)
shot("05_trim_window")
click("trim_apply")
shot("06_details_after_trim")
if (grepl("Feature model", modal_text(), fixed = TRUE)) click("annotate_close")

k <- function(sql) q(sql, u$ID, u$path, u$scaffold)
w <- " WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?"
c1 <- k(paste0("SELECT r.length, r.trim_start, r.trim_end, r.edited, r.check_switch, a.length AS orig_length
  FROM review_units r JOIN assemblies a ON a.ID = r.ID AND a.path = r.path AND a.scaffold = r.scaffold", w))
print(c1)
stopifnot(nrow(c1) == 1, c1$edited == 1, c1$check_switch == 1, c1$trim_start == u$from,
  c1$trim_end == u$to, c1$length == u$to - u$from + 1L)
p <- k(paste0("SELECT r.element, r.start, r.strand, r.partial_start_user, r.partial_stop_user
  FROM review_annotations r", w, " ORDER BY r.start"))
print(p)
cut <- p[p$element == u$first$element, ]
stopifnot(cut$start == 1, sum(!is.na(c(p$partial_start_user, p$partial_stop_user))) == 1,
  if (cut$strand == "-") identical(cut$partial_stop_user, 1L) else identical(cut$partial_start_user, 1L))
stopifnot(identical(orig_snap(u), readRDS("/tmp/rp-t24/orig.rds")))

# the edit shows as not run, and Lock is refused until Update
set_mode("Annotate")
r <- rows()
j <- which(r$ID == u$ID & r$path == u$path & r$scaffold == u$scaffold)
print(r[j, intersect(c("annotate_switch", "check_state", "length", "n_flags", "itsx_warnings"), names(r))])
stopifnot(r$annotate_switch[j] == "1")
shot("07_not_run")
select_rows(j)
click("annot_lock_btn")
txt <- page_text()
shot("08_lock_refused")
stopifnot(grepl("Press Update first.", txt, fixed = TRUE))
close_alert()
stopifnot(k("SELECT annotate_lock FROM annotate r WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?")$annotate_lock == 0)
```

Run: `Rscript /tmp/rp-t24/d5_trim.R`
Expected: no error. The trim made the contig's review copy (S5.1): `edited` 1, `check_switch` 1, the window in original coordinates, the cut gene starting at 1 with the one user partial in its own orientation and nothing marked at the right end, which was trimmed exactly at a gene end (Review Focus 2); the original rows are unchanged (S4.1); the State column reads not run (the hidden `check_state` column, when printed, says `not run`); and Lock refuses it with "Press Update first." (S6.2). Read `07_not_run.png` and `08_lock_refused.png`: the check cells of the trimmed row show the state instead of values (S8).

- [ ] **Step 12: Update on Annotation Review**

Create `/tmp/rp-t24/d6_recheck_start.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
saveRDS(file.info(file.path(PROJ, ".nextflow", "history"))[, c("size", "mtime")], "/tmp/rp-t24/hist0.rds")
set_mode("Annotate")
click("annot_run_btn")
m <- modal_text()
cat(m, "\n")
stopifnot(grepl("RiboPilot recheck - 1 contig(s) to update", m, fixed = TRUE),
  !grepl("Resume previous run?", m, fixed = TRUE), !grepl("-resume", m, fixed = TRUE))
shot("09_recheck_window")
start_run("rc_nf_start", "recheck")
```

Create `/tmp/rp-t24/d7_recheck_wait.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
if (is.null(wait_run("recheck"))) quit(status = 3)
shot("10_recheck_report")
close_alert()
click("rc_nf_close")
u <- readRDS("/tmp/rp-t24/unit.rds")
k <- function(sql) q(sql, u$ID, u$path, u$scaffold)
w <- " WHERE ID = ? AND path = ? AND scaffold = ?"
r <- k(paste0("SELECT check_switch, check_note, checked_at, itsx_warnings, itsx_origin, n_features,
  n_flags, edited, length, blast_source FROM review_units", w))
print(r)
stopifnot(nrow(r) == 1, r$check_switch == 2, is.na(r$check_note), !is.na(r$checked_at),
  !is.na(r$itsx_warnings), r$edited == 1)
stopifnot(r$n_features == k(paste0("SELECT COUNT(*) AS n FROM review_annotations", w))$n,
  k(paste0("SELECT COUNT(*) AS n FROM review_screen", w))$n == 1)
stopifnot(q("SELECT COUNT(*) AS n FROM review_units WHERE check_switch IN (3, 4)")$n == 0)
stopifnot(identical(orig_snap(u), readRDS("/tmp/rp-t24/orig.rds")))
stopifnot(dir.exists(file.path(PROJ, ".recheck")), dir.exists(file.path(PROJ, "work_recheck")),
  identical(file.info(file.path(PROJ, ".nextflow", "history"))[, c("size", "mtime")],
    readRDS("/tmp/rp-t24/hist0.rds")))
set_mode("Annotate")
rr <- rows()
j <- which(rr$ID == u$ID & rr$path == u$path & rr$scaffold == u$scaffold)
print(rr[j, intersect(c("annotate_switch", "check_state", "length", "n_flags", "itsx_warnings",
  "screen_verdict", "chimera_label"), names(rr))])
stopifnot(rr$annotate_switch[j] == "2")
shot("11_after_update")
```

Run: `Rscript /tmp/rp-t24/d6_recheck_start.R`, then `Rscript /tmp/rp-t24/d7_recheck_wait.R` with a 10-minute Bash timeout, repeated while it exits with status 3.
Expected: the recheck window says 1 contig to update and has no Resume box (S7.4); the run (a few minutes: barrnap, both ITSx runs, BLAST, screen, alignment, curation, chimera in `ribopilot:recheck`) ends with a report and no failed contig; the copy is done with a stored ITSx value, its own features counted, and a screen row (S7.3); the originals are still unchanged (the recheck writes only review tables); the run used `<project>/.recheck` and `work_recheck` and left the main pipeline's `.nextflow/history` untouched (S7.4); and the State column reads done. Print `r$blast_source` in the report: it shows whether the recheck searched the panel (`local`) or mapped the original hits. If a contig failed, `check_note` holds the reason: report it with `/tmp/rp-t24/shots/10_recheck_report.png` and the newest `.runs/nextflow/recheck_*.nextflow.log`.

- [ ] **Step 13: Lock and Export**

Create `/tmp/rp-t24/d8_lock.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
u <- readRDS("/tmp/rp-t24/unit.rds")
set_mode("Annotate")
select_rows(seq_len(nrow(rows())))
click("annot_lock_btn")
m <- modal_text()
cat(m, "\n")
shot("12_lock")
if (grepl("Lock anyway", m)) click_text("Lock anyway")
close_alert()
l <- q("SELECT an.ID, an.path, an.scaffold, an.annotate_lock, an.problematic, v.check_switch
  FROM annotate an JOIN review_units_v v ON v.ID = an.ID AND v.path = an.path AND v.scaffold = an.scaffold
  ORDER BY an.ID, an.path, an.scaffold")
print(l)
stopifnot(l$annotate_lock[l$ID == u$ID & l$path == u$path & l$scaffold == u$scaffold] == 1)
# D2: nothing that is not done got locked, unless excluded
stopifnot(!any(l$annotate_lock == 1 & l$check_switch != 2 & !(l$problematic %in% "yes")))
```

Create `/tmp/rp-t24/d9_export.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
set_mode("Export")
r <- rows()
print(r[, intersect(c("ID", "path", "scaffold", "export_group"), names(r))])
select_rows(seq_len(nrow(r)))
click("assign_group_btn")
set_input("group_name", "t24")
click("assign_group_go")
if (grepl("Mixed assemblers", modal_text())) click("group_keep_one")
stopifnot(all(q("SELECT export_group FROM annotate WHERE annotate_lock = 1")$export_group == "t24"))
saveRDS(list(
  an = q("SELECT ID, path, scaffold, itsx_warnings, itsx_origin, itsx_scores, n_flags, curate_flags
    FROM annotate ORDER BY ID, path, scaffold"),
  ru = q("SELECT ID, path, scaffold, itsx_warnings, check_switch, CAST(rev AS TEXT) AS rev
    FROM review_units ORDER BY ID, path, scaffold")), "/tmp/rp-t24/before_export.rds")
click("export_btn")
set_input("export_group_sel", "t24")
m <- modal_text()
stopifnot(grepl("Export data", m, fixed = TRUE), !grepl("Final ITSx check", m, fixed = TRUE))
shot("13_export_window")
click("export_go", timeout = 300)
m <- modal_text()
cat(m, "\n")
if (grepl("Leave them out", m)) {
  shot("14_export_flagged")
  click("export_screen_skip", timeout = 300)
}
if (grepl("Replace it", modal_text())) click("export_overwrite_go", timeout = 300)
stopifnot(wait_text("Export complete", 300))
cat(modal_text(), "\n")
shot("15_export_done")
```

Create `/tmp/rp-t24/d10_export_files.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
u <- readRDS("/tmp/rp-t24/unit.rds")
g <- file.path(PROJ, "out", "export", "t24")
print(list.files(g))
stopifnot(!dir.exists(file.path(g, "itsx_check")), !file.exists(file.path(g, "t24_itsx_check.csv")))

# the summary's itsx_check column is the stored value, NA written as not_run(not checked)
s <- utils::read.csv(file.path(g, "t24_summary.csv"), colClasses = "character", na.strings = character(0))
v <- q("SELECT ID, path, scaffold, itsx_warnings, length FROM review_units_v")
s$key <- paste(s$ID, s$path, s$scaffold)
v$key <- paste(v$ID, v$path, v$scaffold)
m <- merge(s[, c("key", "seqid", "itsx_check", "length")], v[, c("key", "itsx_warnings", "length")],
  by = "key", suffixes = c(".csv", ".db"))
print(m)
stopifnot(nrow(m) == nrow(s),
  identical(m$itsx_check, ifelse(is.na(m$itsx_warnings), "not_run(not checked)", m$itsx_warnings)),
  all(as.integer(m$length.csv) == m$length.db))

# the trimmed contig is exported as its copy, with the cut end marked partial
sid <- m$seqid[m$key == paste(u$ID, u$path, u$scaffold)]
stopifnot(length(sid) == 1)
fa <- readLines(file.path(g, "t24.fasta"))
h <- grep("^>", fa)
at <- h[startsWith(fa[h], paste0(">", sid, " "))]
end <- c(h, length(fa) + 1L)[match(at, h) + 1L] - 1L
seq_fa <- paste(fa[(at + 1L):end], collapse = "")
seq_db <- q("SELECT sequence FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?",
  u$ID, u$path, u$scaffold)$sequence
stopifnot(identical(toupper(seq_fa), toupper(seq_db)))
tb <- readLines(file.path(g, "t24.tbl"))
fh <- grep("^>Feature ", tb)
b0 <- fh[tb[fh] == paste(">Feature", sid)]
b1 <- c(fh, length(tb) + 1L)[match(b0, fh) + 1L] - 1L
block <- tb[b0:b1]
cat(block, sep = "\n")
stopifnot(any(grepl("(^|\t)[<>]1(\t|$)", block)))
stopifnot(grepl("not_run", paste(readLines(file.path(g, "README.txt")), collapse = "\n"), fixed = TRUE))

# export writes no check result
b <- readRDS("/tmp/rp-t24/before_export.rds")
stopifnot(identical(b$an, q("SELECT ID, path, scaffold, itsx_warnings, itsx_origin, itsx_scores, n_flags,
    curate_flags FROM annotate ORDER BY ID, path, scaffold")),
  identical(b$ru, q("SELECT ID, path, scaffold, itsx_warnings, check_switch, CAST(rev AS TEXT) AS rev
    FROM review_units ORDER BY ID, path, scaffold")))
cat("export checks passed\n")
```

Run, in order: `Rscript /tmp/rp-t24/d8_lock.R`, `Rscript /tmp/rp-t24/d9_export.R`, `Rscript /tmp/rp-t24/d10_export_files.R`.
Expected: d8, the trimmed contig and every other done contig lock (CHIM01 asks first and is locked anyway), and no contig that is not done is locked (D2). d9, the Export window has no "Final ITSx check" box (D4); CHIM01 is offered as flagged and left out; the done dialog shows the flagged, clean, and not checked counts (S9); read `13_export_window.png` and `15_export_done.png`. d10 prints `export checks passed`: no `itsx_check/` folder or `t24_itsx_check.csv`; the summary's `itsx_check` column equals the stored values; the trimmed contig's record is the copy's sequence with its cut gene end marked partial (`<1` or `>1`) in the `.tbl` (Review Focus 2, S9); the README explains `not_run`; and export changed no check result (S9).

- [ ] **Step 14: An upgraded project in the app**

Point the same browser tab at the app opened on the upgraded copy of 2026_10_01 (version 3 plus every step, nine locked contigs; Review Focus 5). Create `/tmp/rp-t24/d11_upgraded.R` with the Write tool:

```r
source("/tmp/rp-t24/cdp.R")
reload()
set_mode("Annotate")
r <- rows()
print(r[, intersect(c("ID", "path", "scaffold", "annotate_lock", "annotate_switch", "check_state"), names(r))])
stopifnot(nrow(r) > 0, all(r$annotate_switch[r$annotate_lock == "1"] == "2"),
  !grepl("Error:", page_text(), fixed = TRUE))
shot("16_upgraded_annotation_review")
set_mode("Export")
stopifnot(!grepl("Error:", page_text(), fixed = TRUE))
shot("17_upgraded_export")
```

Run:

```bash
kill "$(cat /tmp/rp-t24/app.pid)"; sleep 3
(nohup Rscript /tmp/rp-t24/app.R /home/dmacguig/Documents/GitHub/RiboPilot /tmp/rp-t24/upgrade/2026_10_01 8791 > /tmp/rp-t24/app2.log 2>&1 & echo $! > /tmp/rp-t24/app.pid)
for i in $(seq 60); do curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8791 | grep -q 200 && break; sleep 2; done
RP_PROJ=/tmp/rp-t24/upgrade/2026_10_01 Rscript /tmp/rp-t24/d11_upgraded.R
```

Expected: Annotation Review lists the copy's contigs, every locked one reads done (S6.1 on originals: `annotate_switch` 2 with a stored ITSx value), and neither step shows an `Error:` output. Read both screenshots.

- [ ] **Step 15: Stop the app and the browser**

Run:

```bash
kill "$(cat /tmp/rp-t24/app.pid)" "$(cat /tmp/rp-t24/chrome.pid)"; sleep 3
rm -rf /dev/shm/rp_t24_chrome
ss -ltn | grep -E ':(8791|9222)\b'
git -C /home/dmacguig/Documents/GitHub/RiboPilot status --short
```

Expected: no listener on either port; git status as in Step 1. Keep `/tmp/rp-t24/shots/`, the test project, and the `ribopilot:recheck` image for the user.

- [ ] **Step 16: Report to the user (no commit)**

This task changes no file in the repository, so there is nothing to commit, and nothing is merged or pushed. Tell the user, in plain words:

1. Static checks (Step 2) and the branch state (Step 1), including whether Task 23 touched Troubleshooting.Rmd.
2. The upgrade of the three project copies (Step 3): versions, view row counts, copies made, and the step's messages.
3. The suite (Step 4): tests run, failures (none expected), and every skip with its reason.
4. The mutation table (Steps 5-6): each name with CAUGHT, MISSED, NOT FOUND (and what you broke by hand instead), or untested.
5. The app drive (Steps 8-14): each checkpoint passed or failed, `blast_source` of the rechecked contig, any id you had to look up, and the screenshot folder `/tmp/rp-t24/shots/`.
6. The image note: `ribopilot:recheck` was built for this check (remove it with `docker rmi ribopilot:recheck` when done). Every other image is now stale for this branch: the bin scripts call `RiboPilot:::review_follow_original`, `RiboPilot:::recheck_unit`, and `RiboPilot:::itsx_store_not_run`, so with an older image the main pipeline's per-contig R steps fail and every recheck contig fails. Rebuild `ribopilot:devel` (`docker/deploy-local.sh devel`) before running the local projects that name it (2026_10_01, 2026_10_08) on this branch, and publish `macguigand/ribopilot:<version>` with `docker/deploy-dockerhub.sh` at release (the user's call).
7. Screenshots to retake for the docs (S13): `vignettes/figures/get-started/annotation-table.png` (the Annotation Review toolbar now has Update, and the State column can read waiting for Assemble Update), `vignettes/figures/get-started/export-window.png` (no Final ITSx check box), `vignettes/figures/get-started/export-done.png` (flagged, clean, and not checked counts), plus any `missing` figure from Step 2. The drive's screenshots are headless captures for reference, not doc-ready.
8. What is left behind: the test project `~/Documents/DJMCompBio/RiboPilot_test/recheck_t24`, the scratch folder `/tmp/rp-t24`, and the upgraded copies in `/tmp/rp-t24/upgrade`.
9. Next step, only when the user agrees (project workflow): `git switch devel && git merge --no-ff feat/contig-recheck`, then `git branch -d feat/contig-recheck`. Never push.

