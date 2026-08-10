# Export review report: feature-sequence gathering, per-feature alignment, and a
# render smoke test of inst/export_report.Rmd.

# Three units. All annotate 18S and 28S; only S1 and S2 annotate ITS1, so the
# gatherer has something to skip. Sequences differ so identity is not 100%.
mk_report_project <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT,
                       problematic TEXT, partial TEXT, reviewed TEXT, chimeric INT,
                       n_features INT, missing TEXT, export_group TEXT)")
  set.seed(42)
  base <- paste(sample(c("A", "C", "G", "T"), 400, replace = TRUE), collapse = "")
  ids <- c("S1", "S2", "S3")
  for (k in seq_along(ids)) {
    id <- ids[k]
    seq <- base
    # a handful of substitutions per sample so the alignment is informative
    for (pos in seq(10 * k, 400, by = 40)) {
      substr(seq, pos, pos) <- c("A", "C", "G", "T")[(k + pos) %% 4 + 1]
    }
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES(?,1,1,?,?,0)",
      params = list(id, seq, nchar(seq)))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES(?,1,1,1,'no','no','yes',0,2,'',NULL)",
      params = list(id))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,1,1,'18S','rRNA',1,200)", params = list(id))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,1,1,'28S','rRNA',251,400)", params = list(id))
    if (id != "S3") {
      DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,1,1,'ITS1','misc_RNA',201,250)",
        params = list(id))
    }
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,?,NULL,NULL)",
      params = list(id, paste("Genus species", k)))
  }
  dir
}

report_units <- function() {
  data.frame(ID = c("S1", "S2", "S3"), path = 1L, scaffold = 1L,
    seqid = c("S1", "S2", "S3"), stringsAsFactors = FALSE)
}

open_project <- function(dir, env = parent.frame()) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con), envir = env)
  con
}

test_that("gather_feature_seqs slices the annotated span from each unit", {
  dir <- mk_report_project()
  con <- open_project(dir)
  full <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID='S1'")$sequence[1]

  s18 <- gather_feature_seqs(con, report_units(), "18S")
  expect_s4_class(s18, "DNAStringSet")
  expect_equal(names(s18), c("S1", "S2", "S3"))
  expect_equal(as.character(s18[["S1"]]), substr(full, 1, 200))
  expect_true(all(Biostrings::width(s18) == 200))

  its <- gather_feature_seqs(con, report_units(), "ITS1")
  expect_equal(names(its), c("S1", "S2"))
  expect_equal(as.character(its[["S1"]]), substr(full, 201, 250))

  expect_equal(length(gather_feature_seqs(con, report_units(), "5.8S")), 0L)
})

test_that("feature_alignment returns NULL for fewer than two sequences", {
  dir <- mk_report_project()
  con <- open_project(dir)
  one <- gather_feature_seqs(con, report_units()[1, , drop = FALSE], "18S")
  expect_equal(length(one), 1L)
  expect_null(feature_alignment(one))
  expect_null(feature_alignment(Biostrings::DNAStringSet(character(0))))
})

test_that("feature_alignment aligns, orders, and summarises three sequences", {
  dir <- mk_report_project()
  con <- open_project(dir)
  seqs <- gather_feature_seqs(con, report_units(), "18S")
  res <- feature_alignment(seqs)
  expect_type(res, "list")
  expect_equal(res$n, 3L)
  expect_equal(res$min_len, 200L)
  expect_equal(res$max_len, 200L)
  expect_true(res$mean_identity > 80 && res$mean_identity <= 100)
  # all rows padded to one alignment width, and every input kept exactly once
  expect_equal(length(unique(Biostrings::width(res$alignment))), 1L)
  expect_setequal(names(res$alignment), c("S1", "S2", "S3"))
})

test_that("report_elements lists every annotated element in canonical order", {
  dir <- mk_report_project()
  con <- open_project(dir)
  expect_equal(report_elements(con, report_units()), c("18S", "ITS1", "28S"))
})

test_that("export_report renders the Rmd with an alignment section", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("msaR")
  skip_on_cran()
  dir <- mk_report_project()
  con <- open_project(dir)
  out <- file.path(dir, "report.html")
  export_report(con, report_units(), out, "testgroup", alignments = TRUE)
  expect_true(file.exists(out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_true(grepl("Feature alignments", html, fixed = TRUE))
  expect_true(grepl("testgroup", html, fixed = TRUE))
  # rendering must not leave intermediates next to the report
  expect_equal(list.files(dirname(out), pattern = "\\.(Rmd|md|knit|utf8)"), character(0))
})
