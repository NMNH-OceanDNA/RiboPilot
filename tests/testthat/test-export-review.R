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

test_that("the user's 5'/3' choices on a copy override curation's in the feature table", {
  p <- mk_review_export()
  con <- p$con
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 51L, 487L)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, itsx_warnings = ''")
  tbl <- function() readLines(file.path(review_bundle(p$dir, group = "g1"), "samples", "S1.tbl"))
  base <- tbl()
  expect_equal(base[c(2, 4)], c("150\t>1\trRNA", "<437\t351\trRNA"))
  # curation's marks stand until the user says otherwise
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 0 WHERE element = '18S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_stop_user = 0 WHERE element = '28S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 1 WHERE element = '28S'")
  expect_equal(tbl()[c(2, 4)], c("<150\t1\trRNA", "437\t351\trRNA"))
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
