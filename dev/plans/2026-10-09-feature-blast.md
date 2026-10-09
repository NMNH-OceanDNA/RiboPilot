# Feature BLAST Implementation Plan

**Goal:** BLAST each core rDNA feature of every kept contig against the bundled panel, in the main and recheck workflows, and show problems in Annotation Review.

**Spec:** `dev/specs/2026-10-09-feature-blast-design.md`

**Branch:** `feat/feature-blast` (off `devel`).

## Done (main session)

- Task 1, schema: `FEATURE_BLAST_DDL`, `FEATURE_BLAST_HITS_DDL`, review mirrors, `FEATURE_BLAST_VIEWS_DDL`
  (`review_feature_blast_v`, `review_feature_blast_hits_v`), `feature_blast_apply_schema()`,
  migration step 5, `UNIT_EVIDENCE_TABLES`, `review_invalidate()` clears the mirrors.
- Task 2, core: `R/feature_blast.R`
  - `feature_blast_unit(con, id, path, scaffold, dbarg = NULL, taxdump_dir = NULL, threads = 1L, workdir, review = FALSE)` (exported)
  - `feature_blast_warning_text(fb)` -> data.frame(ID, path, scaffold, feature_blast_warnings)
  - statuses: ok, no_hit, non_rdna, off_target, too_small, not_run
- Task 4, recheck: `recheck_unit()` calls `feature_blast_unit(review = TRUE)` after the chimera check; failures swallowed.

Deviations from spec:
- non_rdna fires only on hit class "other"; "uninformative" titles (assembly records) are judged by taxonomy.
- No coordinate stale marker: `review_invalidate()` already deletes a copy's feature BLAST rows on any edit,
  and the table gates check cells until the recheck is done.

## Task 3: Main workflow wiring (subagent, Sonnet)

Files: `inst/nextflow/bin/rp_feature_blast.R` (new, executable), `inst/nextflow/modules/pipeline_r.nf`,
`inst/nextflow/main.nf`, tests in `tests/testthat/test-nextflow.R`, `test-nf-bin-contract.R`, `test-review-follow.R`.

- Script modeled on `rp_chimera.R`: args `<db> <id> <path> <scaffold> [dbarg] [taxdump_dir]`; calls
  `feature_blast_unit(con, id, path, scaffold, dbarg, taxdump_dir, threads = RP_THREADS, workdir = getwd()/feature_blast_work)`,
  then `RiboPilot:::review_follow_original()`; prints `feature_blast: <id> <path> <scaffold> - <n> features`.
- Process `feature_blast` in `pipeline_r.nf` next to `chimera`: `cache false`, `errorStrategy 'ignore'`,
  same DBARG block as `chimera`.
- `main.nf`: include it and call `feature_blast(reconcile.out)`; not joined into triage.
- Tests: follow existing contract tests for every rp_*.R script and process wiring.

## Task 5: Annotation Review UI (subagent, Sonnet)

Files: `R/app_server.R`, `R/app_modals.R`, `R/constants.R`, tests `tests/testthat/test-app-*.R`.

- Table: in the Annotate branch of `all_rows`, after the query, read
  `SELECT ID, path, scaffold, element, status, lineage FROM review_feature_blast_v`, run
  `feature_blast_warning_text()`, left-join by (ID, path, scaffold) into `d$feature_blast_warnings`
  ("" for units with none). Column def `feature_blast_warnings`: header "Feature BLAST warnings" with
  `rt_help_header` and a new `SCAFFOLD_COL_HELP$feature_blast` entry, `rt_gated(rt_longtext())`,
  minWidth ~180; add to the `Check details` column group after `chimera_label`.
- Modal: new `rp_details("Feature BLAST", ...)` panel in `annotate_modal()` after the Chimera check panel,
  help text, `reactable::reactableOutput("feature_blast_tbl")`. Server output reads
  `review_feature_blast_v` and `review_feature_blast_hits_v` for `sel_unit()`, gated like the chimera panel
  (only when the check state is "done"). Columns: Feature, Length, Status (readable label; low confidence
  note for low_conf = 1), Call (species), Match rank, %ID, Cov, E. Expandable row `details` showing the
  stored top hits (accession with `rt_ncbi_link()`, title, %ID, Cov, Bits, class). Empty state text when no rows.
- Status labels: ok "OK", no_hit "No hit", non_rdna "Non-rDNA hit", off_target "Off-target",
  too_small "No BLAST, feature too small", not_run "Not run".

## Task 6: Docs (subagent, Haiku or Sonnet)

- `vignettes/Test-Project-Annotate.Rmd`: short subsection on the Feature BLAST warnings column and panel,
  length rules (under 50 bp skipped, 50 to 99 bp low confidence), phylum off-target rule, tie collapse. Say "contig".
- `NEWS.md` under "Quality checks": one bullet.

## Task 7: Verify

- Full relevant test files; headless app drive (chromote) to see column and panel render.
