# Metadata import (GEOME, GBIF, NCBI) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port MitoPilot's specimen metadata feature into RiboPilot, same method, in three stages.

**Architecture:** Copy MitoPilot's feature files into `R/` with internal names unchanged; wire them into RiboPilot's single-table app, its `samples`-based export, and its schema-step migration. Network goes through three mockable functions; tests run on recorded fixtures.

**Tech Stack:** R package, SQLite through DBI/RSQLite, httr2, xml2, jsonlite, shiny, reactable, testthat 3e.

**Spec:** `dev/specs/2026-09-30-metadata-import-design.md` (read it first; it is the authority on behaviour and deviations).

**MitoPilot source:** `/home/dmacguig/Documents/GitHub/MitoPilot` (called `$MP` below). Read-only.

## Global constraints

- ASCII only in code and docs. No em dashes. Minimal comments.
- Internal function names, table names, column names, and key formats stay exactly as in MitoPilot.
- `NAMESPACE` is hand-maintained: add `export()` lines by hand, never run `devtools::document()` for NAMESPACE. `man/` pages come from `roxygen2::roxygenise()` (it skips NAMESPACE).
- Schema: every table or column change is in `new_db()` AND the one new step in `db_migrations()`. `tests/testthat/fixtures/schema_v1.sql` never changes.
- Wording for people: "contig" not "scaffold"; "test project" not "practice project".
- No git commands that change state inside a task unless the task says so. No push.
- Every task ends with the named test files passing: `Rscript -e 'devtools::test(filter = "<pattern>", reporter = "summary")'`.
- After each stage, the full suite: `Rscript -e 'devtools::test(reporter = "summary")'` with 0 failures and no new skips beyond the baseline 23.

## Review focus

- A mapping file with `GEOME_BCID`, `GBIF_ID`, or `BioSample` columns must be accepted by `read_mapping()` (not treated as reserved) and must not break `test-mapping-checks.R`'s reserved-names test.
- A sample whose identifier is invalid must never trigger a network call and must show a failed status with the source's message.
- A fetch that fails after a previous success must keep the old records and only change status.
- A ticked export column whose name collides with a computed token (`length`, `partial`, `path`) is shadowed by the computed one; the header must still render.
- A project on schema version 1 must upgrade to version 2 through `backwards_compatibility()` and then match a fresh `new_db()` table for table and column for column.

---

## Stage 1

### Task 1: dependencies, schema, DDL, migration step

**Files:**
- Modify: `DESCRIPTION` (Imports: add `httr2`, `xml2`; Suggests: add `curl`)
- Modify: `R/db_schema.R` (seven DDL constants)
- Modify: `R/init_db.R` (`new_db()` creates tables and columns)
- Modify: `R/db_migrate.R` (`db_migrations()` step 1)
- Modify: `tests/testthat/helper-app-fixture.R` (`fixture_fill_schema()` list, `samples` columns)
- Test: `tests/testthat/test-backwards-compat.R`, `tests/testthat/test-mapping-checks.R`

**Interfaces:**
- Produces: `META_RECORDS_DDL`, `META_STATUS_DDL`, `META_EXPORT_FIELDS_DDL`, `META_CSV_MAP_DDL`, `META_OPTIONS_DDL`, `META_LINKS_DDL`, `META_VIEW_FIELDS_DDL` (character(1) each, `CREATE TABLE IF NOT EXISTS`); `META_SAMPLE_COLS <- c("GEOME_BCID", "GBIF_ID", "BioSample")`; `meta_apply_schema(con)` which runs all DDL and adds any missing `samples` column with `ALTER TABLE samples ADD COLUMN <col> TEXT`.

- [ ] Add the seven DDL constants (verbatim SQL from the spec) and `META_SAMPLE_COLS` and `meta_apply_schema()` to `R/db_schema.R`.
- [ ] In `new_db()`, call `meta_apply_schema(con)` right after `ensure_export_opts(con)` and before `record_schema_event(con, "created", ...)`.
- [ ] In `db_migrations()`, return `list(list(note = "sample metadata (GEOME/GBIF/NCBI) tables and columns", fn = meta_apply_schema))`.
- [ ] Add the constants to `fixture_fill_schema()`'s DDL list and the three columns via its `add_col("samples", ...)` path so app fixtures match.
- [ ] Add `httr2`, `xml2` to Imports and `curl` to Suggests in DESCRIPTION; confirm both are installed (`Rscript -e 'library(httr2); library(xml2)'`).
- [ ] Run `devtools::test(filter = "backwards-compat|mapping-checks|app-table|app-export")`. Expected: parity test passes with version 2; nothing else changes.

### Task 2: API layer (normalizers, fetch chains, HTTP)

**Files:**
- Create: `R/gbif_api.R`, `R/geome_api.R`, `R/ncbi_api.R` (copied from `$MP/R/`)
- Create: `tests/testthat/helper-gbif.R`, `helper-geome.R`, `helper-ncbi.R`, `tests/testthat/fixtures/gbif/*`, `fixtures/geome/*`, `fixtures/ncbi/*` (copied from `$MP/tests/testthat/`)
- Create: `tests/testthat/test-gbif-api.R`, `test-geome-api.R`, `test-ncbi-api.R` (copied, package name swapped)
- Modify: `NAMESPACE` (export `geome_normalize_bcid`, `gbif_normalize_id`, `ncbi_normalize_id`)

**Interfaces:**
- Produces: `.gbif_get(path)`, `.geome_get(path, query = list())`, `.ncbi_get(endpoint, query)`; `.gbif_fetch_chain(ref, cache)`, `.geome_fetch_chain(ref, cache)`, `.ncbi_fetch_chain(ref, cache)` returning `data.frame(level, depth, ref, field, value)`; the three normalizers; `.check_ncbi_reachable()`.

- [ ] Copy the three api files. Change the user agent string to `RiboPilot (https://github.com/dmacguigan/RiboPilot)` and the eutils `tool` value to `RiboPilot`.
- [ ] In `.ncbi_get`, replace `key <- Sys.getenv("ENTREZ_KEY")` with `key <- ncbi_meta_key()`, and add `ncbi_meta_key <- function() { k <- tryCatch(project_ncbi_key("."), error = function(e) ""); if (!nzchar(k)) k <- Sys.getenv("ENTREZ_KEY", ""); k }` in `R/ncbi_api.R`. Keep the rate gap logic.
- [ ] Confirm every helper the files call exists in RiboPilot (`%||%`, `%|NA|%`); if MitoPilot's `%nin%` or `.lst()` is used, copy it into `R/utils.R` (check `test-no-duplicate-definitions.R` afterwards).
- [ ] Copy fixtures and helpers; in the helper files replace `MitoPilot` with `RiboPilot`. Copy the three api test files and swap `MitoPilot:::`/`.package = "MitoPilot"` for RiboPilot. Adjust the `.ncbi_get` key test to set `ENTREZ_KEY` with no project on disk (it must still see the env key).
- [ ] Add the three exports to NAMESPACE. Run `roxygen2::roxygenise()` and keep only the new `man/*.Rd` files plus wording-only diffs.
- [ ] Run `devtools::test(filter = "gbif-api|geome-api|ncbi-api|no-duplicate")`. Expected: all pass; live tests skip offline.

### Task 3: storage and console fetch functions

**Files:**
- Create: `R/meta_db.R`, `R/gbif_db.R`, `R/geome_db.R`, `R/ncbi_db.R`, `R/meta_link.R` (copied; linking code lands now but stays off until stage 2 wires the option; the `remove_metadata` function lives in `meta_link.R`)
- Modify: `R/meta_db.R` (`.meta_ensure_tables` becomes a call to `meta_apply_schema(con)`; drop the legacy `geome_*` migration block)
- Modify: `NAMESPACE` (export `fetch_geome`, `fetch_gbif`, `fetch_biosample`, `remove_metadata`)
- Create: `tests/testthat/test-meta-db.R`, `test-geome-db.R`, `test-gbif-db.R`, `test-ncbi-db.R`, `test-meta-remove.R` (copied, adapted: drop the legacy-migration test and any `new_db_userAsmb` cases)

**Interfaces:**
- Consumes: Task 1 schema, Task 2 chains.
- Produces: `META_SOURCES`, `.meta_set_ref(con, id, source, value)`, `.meta_fetch_into(con, source, ids, cache, link)`, `.meta_fetch_project(path, source, ids, refs, link)`, `.meta_fetch_new(con, mapping, cols, link)`, `.meta_take_cols(mapping, ...)`, `.meta_drop`, `.meta_remove`, `.meta_link_enabled`, `.meta_set_link_enabled`, `fetch_geome()`, `fetch_gbif()`, `fetch_biosample()`, `remove_metadata()`.

- [ ] Copy the files; make the edits above. Each console function opens `file.path(path, ".sqlite")` as in MitoPilot but through `open_project_db()` so busy timeouts match.
- [ ] Copy and adapt the tests. Where a test builds a project with MitoPilot's `new_db()`, use RiboPilot's `new_db()` with `check_reads = FALSE` and the dummy seeds/labels fasta pattern from `test-backwards-compat.R`'s `bc_project()`; those cases that depend on Task 4 wiring (`new_db(mapping_geome = ...)`) move to Task 4's test file.
- [ ] Run `devtools::test(filter = "meta-db|geome-db|gbif-db|ncbi-db|meta-remove")`. Expected: pass.

### Task 4: setup wiring (mapping file, new_project, add_samples, update_sample_metadata)

**Files:**
- Modify: `R/init_db.R` (`read_mapping()`, `check_mapping()`, `new_db()` args and `.meta_fetch_new` call)
- Modify: `R/init_project.R` (`new_project()` args pass-through, reachability pre-flight)
- Modify: `R/init_test_project.R` (pass `...` through, already does)
- Modify: `R/add_samples.R`, `R/update_sample_metadata.R` (args, `.meta_take_cols`, `.meta_fetch_new` / `.meta_sync_changed`)
- Test: `tests/testthat/test-geome-init.R`, `test-geome-update.R` (copied, adapted), `test-mapping-checks.R` (extend)

**Interfaces:**
- Consumes: Task 3.
- Produces: the seven new arguments on `new_db()`, `new_project()`, `add_samples()`, `update_sample_metadata()` with MitoPilot defaults; `.meta_sync_changed(con, old, new, ...)`.

- [ ] Port the metadata parts of `$MP/R/init_checks.R` `check_mapping` (invalid-ID warnings) into RiboPilot's `check_mapping()`; port the `check_project_init` reachability block into `new_project()` after `check_container_image()`, using `.check_resource()` copied into `R/utils.R` if missing.
- [ ] `read_mapping()`: the three columns are matched ignoring case like `ID`/`Taxon` and renamed to the canonical names.
- [ ] `new_db()`: add the arguments; after the samples upsert, `if (isTRUE(link_sources)) .meta_set_link_enabled(con, TRUE)`; then `.meta_fetch_new(con, mapping, list(GEOME = fetch_geome, GBIF = fetch_gbif, NCBI = fetch_biosample), link = link_sources)`; roxygen `@param` text copied from MitoPilot.
- [ ] `add_samples()` and `update_sample_metadata()`: same arguments; `update_sample_metadata()` keeps the Taxon re-queue behaviour and adds `.meta_sync_changed`.
- [ ] Adapt `test-geome-init.R` and `test-geome-update.R`; add a `test-mapping-checks.R` case that a mapping file with `GEOME_BCID` and `BioSample` columns is accepted and the values land normalized in `samples`.
- [ ] Run `devtools::test(filter = "geome-init|geome-update|mapping-checks|gbif-db|ncbi-db")`. Expected: pass.

### Task 5: export fields, combos, header data, summary

**Files:**
- Create: `R/meta_export.R`, `R/geome_export.R`, `R/gbif_export.R`, `R/ncbi_export.R` (copied)
- Modify: `R/export_opts.R` (`export_header_data()` joins ticked columns via `.meta_join`)
- Modify: `R/export.R` (summary rows get the ticked columns; `<group>_summary.csv` header extended)
- Test: `test-geome-export.R`, `test-gbif-export.R`, `test-ncbi-export.R` (copied), `test-export-headers.R` (extend: a ticked `geome_lat_lon` works as `{geome_lat_lon}`; an unticked field errors as "not in the database")

**Interfaces:**
- Produces: `GEOME_COMBOS`, `GBIF_COMBOS`, `NCBI_COMBOS`, `meta_export_cols(con, ids, keys)`, `.meta_join(dat, con)`, `meta_field_summary(con, source)`, `.meta_save_fields(con, keys)`, `.meta_key_col(key)`.

- [ ] Copy the four files unchanged.
- [ ] In `export_header_data()`, after the per-unit rows are built and before the NA blanking, `rows <- .meta_join(rows, con)` (join by `ID`); keep the "common columns only" reduce so every row has the same names.
- [ ] In `export.R`, join the same columns onto the summary data frame so they appear in `<group>_summary.csv` after the existing columns.
- [ ] Copy and adapt the three export test files; extend `test-export-headers.R`.
- [ ] Run `devtools::test(filter = "geome-export|gbif-export|ncbi-export|export-headers|export-opts|app-export")`. Expected: pass.

### Task 6: viewer modal, Metadata column, Set Export Metadata, token chips

**Files:**
- Create: `R/app_specimen.R`, `R/app_export_tokens.R` (copied)
- Create: `inst/app/www/specimen/geome_g.png`, `gbif_leaf.png`, `ncbi_helix.png`, `inst/extdata/countries.csv` (copied)
- Modify: `inst/app/www/custom.css` (append the `mp-spec-*`, `mp-meta-*`, `mp-token-*`, and logo rules from `$MP/inst/app/www/custom.css`)
- Modify: `R/app_server.R` (Metadata column in Taxonomy/Assemble/Annotate/Export data and colDefs; `specimen_viewer_server` wiring; `specimen_fields` handlers; Export modal token chips)
- Modify: `R/app_ui.R` (Export toolbar: "Set Export Metadata" button)
- Modify: `R/app_modals.R` (`export_modal()`: replace the "Available columns" details block with `export_token_ui(...)`)
- Test: `tests/testthat/test-specimen-app.R`, `test-export-tokens.R`, `test-export-metadata-cols.R` (copied, adapted), `test-app-export.R` (extend for chips), `test-app-table.R` (Metadata column present)

**Interfaces:**
- Consumes: Tasks 3 and 5.
- Produces: `specimen_col_def(input_id, sticky)`, `.specimen_status_join(db, con)`, `specimen_viewer_server(id, open, on_change, con)`, `specimen_fields_modal(...)`, `export_token_groups(data, sample_cols, ticked_keys)`, `export_token_ui(groups, ns)`.

- [ ] Copy `app_specimen.R` and `app_export_tokens.R`. Replace `session$userData$con` with a `con` argument (or `app_con()` with `on.exit(DBI::dbDisconnect(con))`). Replace every `gargoyle` call with a `reactiveVal` bump passed in (`bump = function() ...`). Replace `mp_confirm` calls with RiboPilot's `mp_confirm` (same name, `R/app_alerts.R`).
- [ ] `app_server.R`: in `all_rows()` for modes Taxonomy, Assemble, Annotate, Export, run `.specimen_status_join()` and add a `specimen` colDef after Taxon (`specimen_col_def("specimen_open")`); wire `specimen_viewer_server("specimen", open = reactive(input$specimen_open), on_change = <bump the table>)`. Add the `specimen_fields` observers ported from `$MP/R/app_export.R:223-282` (snapshot the Export modal, open the fields modal, reopen on close).
- [ ] `app_ui.R`: add `btn("specimen_fields", "Set Export Metadata", ...)` to `export_ctrls` with the MitoPilot title.
- [ ] `export_modal()`: drop the "Available columns" details block and render `export_token_ui(export_token_groups(...))` in its place; the chips insert into `fasta_header` / `fasta_header_marker` with the same JS as MitoPilot (copy the inline script).
- [ ] Copy the CSS rules and the three logos and `countries.csv`.
- [ ] Copy and adapt the tests; extend `test-app-export.R` and `test-app-table.R`.
- [ ] Run `devtools::test(filter = "specimen-app|export-tokens|export-metadata-cols|app-export|app-table|app-headless")`. Expected: pass.

### Task 7: stage 1 close

- [ ] Full suite: `devtools::test(reporter = "summary")` with 0 failures.
- [ ] `roxygen2::roxygenise()`; check `git diff man/` is only new pages and wording.
- [ ] Commit on `feat/metadata-import`: `feat: sample metadata from GEOME, GBIF, and NCBI (stage 1)`.

## Stage 2

### Task 8: link-following

**Files:**
- Modify: `R/meta_link.R` (already copied in Task 3; now wired), `R/app_specimen.R` (the "Follow links between databases" checkbox and the second Refresh-all pass)
- Test: `tests/testthat/test-meta-links.R` (copied, adapted)

- [ ] Confirm `.meta_link_sample`, `.meta_link_candidates`, `.gbif_find_voucher`, `.geome_fastq_biosample` are present and unchanged; wire the checkbox to `.meta_set_link_enabled` shared across tabs.
- [ ] Copy and adapt `test-meta-links.R`.
- [ ] Run `devtools::test(filter = "meta-links|specimen-app")`. Expected: pass.

### Task 9: conflicts, Compare tab, set_metadata_columns, export warning

**Files:**
- Create: `R/specimen_reconcile.R` (copied)
- Modify: `R/app_specimen.R` (Compare tab and "Mapfile columns..." panel, already in the copied file; confirm enabled), `R/app_server.R` (before an export, `specimen_template_concepts()` + `specimen_export_warnings()` + `mp_confirm("Sample metadata disagrees", "Export anyway")` ported from `$MP/R/app_export.R:1487-1560`)
- Modify: `NAMESPACE` (export `set_metadata_columns`)
- Test: `tests/testthat/test-specimen-conflicts.R` (copied, adapted), `test-app-export.R` (extend: a conflicting `{geome_lat_lon}` asks first)

- [ ] Copy the file; `app_sys("extdata", "countries.csv")` must resolve (RiboPilot has `app_sys()`; confirm its signature).
- [ ] Wire the export confirm in the export flow in `app_server.R`.
- [ ] Copy and adapt tests; extend `test-app-export.R`.
- [ ] Run `devtools::test(filter = "specimen-conflicts|specimen-app|app-export")`. Expected: pass.

### Task 10: stage 2 close

- [ ] Full suite, 0 failures. `roxygenise()`. Commit: `feat: metadata links and conflict checks (stage 2)`.

## Stage 3

### Task 11: Metadata column picker

**Files:**
- Create: `R/meta_view.R`, `R/app_meta_view.R` (copied)
- Modify: `R/app_ui.R` (the `Metadata` button in both `mp-filter-row` blocks, and a Taxonomy filter row if none exists), `R/app_server.R` (`meta_view_setup`, `meta_view_join` on the four views' data, `meta_view_col_defs` merged into `cols`)
- Test: `tests/testthat/test-meta-view.R` (copied, adapted), `test-app-table.R` (extend: a shown field appears as a column)

- [ ] Copy the files; replace gargoyle with the `reactiveVal` bump.
- [ ] Wire button, modal, save, and join.
- [ ] Copy and adapt the test; extend `test-app-table.R`.
- [ ] Run `devtools::test(filter = "meta-view|app-table|specimen-app")`. Expected: pass.

### Task 12: docs

**Files:**
- Create: `vignettes/Sample-Metadata.Rmd` (from `$MP/vignettes/Specimen-Metadata.Rmd`; RiboPilot names; drop figure lines; add a final section "Screenshots to add" listing the eight MitoPilot figure subjects)
- Modify: `_pkgdown.yml` (article under Usage after `Your-Own-Project`; reference section "Sample metadata" listing the new exports), `vignettes/Your-Own-Project.Rmd` (mapping-file section: the three ID columns), `vignettes/GenBank-Submission.Rmd` ("Required source modifiers": fetched tokens), `NEWS.md` (new section "Sample metadata")
- Test: `tests/testthat/test-docs-site.R`

- [ ] Write the pages; `pkgdown::check_pkgdown()` clean; `pkgdown::build_site(install = FALSE, lazy = FALSE)` clean; link check over `docs/articles`.
- [ ] Run `devtools::test(filter = "docs-site")`. Expected: pass.

### Task 13: stage 3 close

- [ ] Full suite, 0 failures. Commit: `feat: metadata column picker and guide (stage 3)`.
