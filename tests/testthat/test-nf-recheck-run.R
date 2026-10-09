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
