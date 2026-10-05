# The real main.nf, launched with the app's own command and environment, with
# the tools faked (helper-fake-tools.R). The first test below runs every WF1
# step on the bundled samples (reads -> GetOrganelle) plus the re-queue arm
# (annotate-only re-run); further tests below exercise the depth-gate, the
# fallback-assembler branch (routing only: no spades.py is faked, so no real
# contig comes out of it), and a user-supplied assembly with find_rrna off.
# Not covered: a successful fallback assembly (needs real SPAdes + a working
# local BLAST db).

test_that("main.nf runs every WF1 step on the bundled samples", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  # Out-of-range Nextflow gets pinned (nf_env()) to NF_MAX_SUPPORTED, which the
  # launcher downloads if not already cached -- skip rather than hang offline.
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

  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 600, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)

  con <- RiboPilot:::open_project_db(db)
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

  # Order, not just presence: curate() reads blast_ref_alignment, so it cannot
  # have finished before that table's rows landed. Catches e.g.
  # curate(reconcile.out) swapped in for curate(blast_ref_align.out) in main.nf.
  expect_gte(q("SELECT MAX(time_stamp) AS t FROM annotate")$t,
    q("SELECT MAX(time_stamp) AS t FROM blast_ref_alignment")$t)
  # chimera reads coverage (depth/clips) and blast_hits, so it cannot have
  # started before both barriers closed.
  chim_start <- q("SELECT MIN(time_stamp) AS t FROM chimera")$t
  expect_gte(chim_start, q("SELECT MAX(time_stamp) AS t FROM coverage")$t)
  expect_gte(chim_start, q("SELECT MAX(time_stamp) AS t FROM blast_searched")$t)

  # Re-queue arm: annotate-only re-run, no reassembly/preprocess. assemblies
  # table has no time_stamp of its own for annotate_switch, so use annotate's.
  asm_before <- q("SELECT time_stamp FROM assemble ORDER BY ID")$time_stamp
  pre_before <- q("SELECT trimmed_reads FROM preprocess ORDER BY ID")$trimmed_reads
  unit <- q("SELECT ID, path, scaffold, time_stamp FROM annotate LIMIT 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(unit$ID, unit$path, unit$scaffold))
  DBI::dbDisconnect(con)

  res2 <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 600, cleanup_tree = TRUE)
  expect_equal(res2$status, 0L, info = res2$stderr)

  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  requeued <- q(sprintf("SELECT annotate_switch, time_stamp FROM annotate WHERE ID = '%s' AND path = %d AND scaffold = %d",
    unit$ID, unit$path, unit$scaffold))
  expect_equal(requeued$annotate_switch, 2L)  # back from the forced 1: proves it was re-processed
  expect_gte(requeued$time_stamp, unit$time_stamp)
  expect_equal(q("SELECT time_stamp FROM assemble ORDER BY ID")$time_stamp, asm_before)
  expect_equal(q("SELECT trimmed_reads FROM preprocess ORDER BY ID")$trimmed_reads, pre_before)
})

test_that("main.nf gates a sample below min_depth without assembling it", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
  tools <- local_fake_pipeline_tools()
  p <- file.path(withr::local_tempdir(), "proj")
  # fake fastp always reports trimmed_reads = 2000 (both mates); set min_depth
  # above that so below_min_depth fires for both bundled samples.
  suppressMessages(suppressWarnings(RiboPilot::new_test_project(p, n = 0, executor = "conda",
    custom_seeds_db = tools$scaffold, custom_labels_db = tools$scaffold, Rproj = FALSE,
    min_depth = 2500, fetch_metadata = FALSE)))
  db <- file.path(p, ".sqlite")
  con <- RiboPilot:::open_project_db(db)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  DBI::dbDisconnect(con)
  no_panels <- withr::local_tempdir()

  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)

  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)
  a <- q("SELECT assemble_switch, assemble_notes, paths, scaffolds FROM assemble ORDER BY ID")
  expect_equal(a$assemble_switch, c(3L, 3L))
  expect_true(all(grepl("^below min_depth", a$assemble_notes)))
  expect_true(all(is.na(a$paths)))
  expect_equal(q("SELECT COUNT(*) AS n FROM assemblies")$n, 0L)
})

test_that("main.nf routes assembler = 'bait' through the fallback arm, not GetOrganelle", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
  # No faked spades.py, so a real fallback assembly can't be produced here;
  # this only proves the fallback branch ran instead of GetOrganelle.
  tools <- local_fake_pipeline_tools()
  p <- file.path(withr::local_tempdir(), "proj")
  suppressMessages(suppressWarnings(RiboPilot::new_test_project(p, n = 0, executor = "conda",
    custom_seeds_db = tools$scaffold, custom_labels_db = tools$scaffold, Rproj = FALSE, fetch_metadata = FALSE)))
  db <- file.path(p, ".sqlite")
  con <- RiboPilot:::open_project_db(db)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  DBI::dbExecute(con, "UPDATE assemble_opts SET assembler = 'bait'")
  DBI::dbDisconnect(con)
  no_panels <- withr::local_tempdir()

  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)

  for (id in c("CHIM00", "CHIM01")) {
    expect_false(file.exists(file.path(p, "out", id, "assemble", "get_org.log.txt")), info = id)
  }
  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  a <- DBI::dbGetQuery(con, "SELECT assemble_switch, assemble_notes FROM assemble ORDER BY ID")
  expect_equal(a$assemble_switch, c(3L, 3L))
  expect_true(all(grepl("^failed assembly", a$assemble_notes)))
})

test_that("main.nf finds rDNA in a user-supplied assembly (find_rrna off)", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
  tools <- local_fake_pipeline_tools()
  proj_dir <- withr::local_tempdir()
  p <- file.path(proj_dir, "proj")
  asmb_dir <- file.path(proj_dir, "asmb")
  dir.create(asmb_dir)
  file.copy(tools$scaffold, file.path(asmb_dir, "SUPPLIED_assembly.fa"))
  map_fn <- file.path(proj_dir, "mapping.csv")
  writeLines(c("ID,Taxon,Assembly", "SUPPLIED,Upeneus parvus,SUPPLIED_assembly.fa"), map_fn)

  suppressMessages(suppressWarnings(RiboPilot::new_project(path = p, mapping_fn = map_fn,
    assembly_path = asmb_dir, executor = "conda", custom_seeds_db = tools$scaffold,
    custom_labels_db = tools$scaffold, Rproj = FALSE)))
  db <- file.path(p, ".sqlite")
  con <- RiboPilot:::open_project_db(db)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  # find_rrna = 0: take the supplied file as given (stage_assembly), skipping
  # the find_rrna_scan/pick module entirely.
  DBI::dbExecute(con, "UPDATE assemble_opts SET find_rrna = 0")
  DBI::dbDisconnect(con)
  no_panels <- withr::local_tempdir()

  res <- processx::run("nextflow", RiboPilot:::nf_command(p, resume = FALSE), wd = p,
    env = c(RiboPilot:::nf_env(), RIBOPILOT_DB_DIR = no_panels,
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)

  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)
  a <- q("SELECT find_rrna_notes FROM assemble WHERE ID = 'SUPPLIED'")
  expect_match(a$find_rrna_notes, "using the supplied assembly as given", fixed = TRUE)
  u <- q("SELECT n_features FROM annotate WHERE ID = 'SUPPLIED'")
  expect_equal(nrow(u), 1L)
  expect_equal(u$n_features, 5L)
})
