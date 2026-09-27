# The real main.nf, launched with the app's own command and environment on the
# practice project's two bundled samples, with the tools faked
# (helper-fake-tools.R). Every WF1 step must leave its rows behind.

test_that("main.nf runs every WF1 step on the bundled samples", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
  tools <- local_fake_pipeline_tools()
  p <- file.path(withr::local_tempdir(), "proj")
  suppressMessages(suppressWarnings(RiboPilot::new_test_project(p, n = 0, executor = "conda",
    custom_seeds_db = tools$scaffold, custom_labels_db = tools$scaffold, Rproj = FALSE)))
  db <- file.path(p, ".sqlite")
  con <- RiboPilot:::open_project_db(db)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  DBI::dbDisconnect(con)
  no_panels <- withr::local_tempdir()

  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 600)
  expect_equal(res$status, 0L, info = res$stderr)

  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)
  expect_equal(q("SELECT trimmed_reads FROM preprocess ORDER BY ID")$trimmed_reads, c(2000L, 2000L))
  expect_equal(q("SELECT assemble_switch FROM assemble ORDER BY ID")$assemble_switch, c(2L, 2L))
  u <- q("SELECT an.n_features, an.n_flags, asm.triage FROM annotate an
    JOIN assemblies asm USING (ID, path, scaffold)")
  expect_equal(u$n_features, c(5L, 5L))
  expect_false(anyNA(u$n_flags))
  expect_false(anyNA(u$triage))
  for (tbl in c("blast_searched", "screen", "chimera", "coverage")) {
    expect_equal(q(paste("SELECT COUNT(*) AS n FROM", tbl))$n, 2L, info = tbl)
  }
  expect_gt(q("SELECT COUNT(*) AS n FROM blast_ref_alignment")$n, 0L)
})
