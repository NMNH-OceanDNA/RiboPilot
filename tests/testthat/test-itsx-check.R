# Stock ITSx check: parsers over stock-ITSx output, the pure flag function,
# storing the result, and what export reports from it.

pos_clean <- "S1\t15120 bp.\tSSU: 1-6832\tITS1: 6833-7243\t5.8S: 7244-7401\tITS2: 7402-7928\tLSU: 7929-15120\t"
pos_chim <- "S2\t3000 bp.\tSSU: 1-1800\tITS1: 1801-2000\t5.8S: 2001-2160\tITS2: 2161-2400\tLSU: 2401-2900\tChimeric! "
pos_no58 <- "S3\t3000 bp.\tSSU: 1-1800\tITS1: 1801-2000\t5.8S: Not found\tITS2: 2161-2400\tLSU: 2401-2900\tBroken or partial sequence, no 5.8S! ITS region too long! "
res_clean <- "S1\t15120\tM\t0\t4\t2.9e-16\t49.9\t49.9\t6788\t7973\t1_SSU\t4_LSU\t\t1_SSU,2_5.8,3_End,4_LSU\tF: 4 1e-08 33.15, M: 4 2.9e-16 49.9"
res_reorder <- "S4\t3000\tF\t0\t3\t1e-10\t40\t40\t100\t2900\t2_5.8\t1_SSU\tChimeric\t2_5.8,1_SSU,4_LSU\tF: 3 1e-10 40"

mk_its <- function(lines) RiboPilot:::parse_itsx(lines)
no_res <- RiboPilot:::parse_itsx_results(character(0))
no_prob <- RiboPilot:::parse_itsx_problematic(character(0))

test_that("parse_itsx_notes pulls every trailing warning off a positions line", {
  n <- RiboPilot:::parse_itsx_notes(c(pos_clean, pos_chim, pos_no58))
  expect_equal(n$seqid, c("S1", "S2", "S3"))
  expect_equal(n$notes[1], "")
  expect_equal(n$notes[2], "Chimeric!")
  expect_equal(n$notes[3], "Broken or partial sequence, no 5.8S!; ITS region too long!")
  expect_equal(nrow(RiboPilot:::parse_itsx_notes(character(0))), 0)
})

test_that("parse_itsx_results reads origin and the domain-order chimera column", {
  r <- RiboPilot:::parse_itsx_results(c(res_clean, res_reorder))
  expect_equal(r$origin, c("M", "F"))
  expect_equal(r$chimeric, c(FALSE, TRUE))
  expect_equal(r$domain_order[1], "1_SSU,2_5.8,3_End,4_LSU")
  expect_equal(nrow(no_res), 0)
})

test_that("parse_itsx_problematic keeps the free-text reason", {
  p <- RiboPilot:::parse_itsx_problematic("S5\tNo ITS1 sequence; LSU seems to overlap 5.8S, sequence may be chimeric")
  expect_equal(p$seqid, "S5")
  expect_match(p$problem, "may be chimeric")
})

test_that("itsx_expected_origins maps a lineage to consistent ITSx groups", {
  expect_equal(RiboPilot:::itsx_expected_origins("Metazoa;Chordata;Actinopteri"), "M")
  expect_setequal(RiboPilot:::itsx_expected_origins("Viridiplantae;Streptophyta;Magnoliopsida"), c("T", "B", "L"))
  expect_null(RiboPilot:::itsx_expected_origins("Bacteria;Proteobacteria"))
  expect_null(RiboPilot:::itsx_expected_origins(NA_character_))
})

test_that("a clean record that agrees with RiboPilot gets no flags", {
  ann <- data.frame(element = c("18S", "ITS1", "5.8S", "ITS2", "28S"),
    start = c(4998L, 6833L, 7244L, 7402L, 7933L), stop = c(6831L, 7243L, 7401L, 7928L, 11991L))
  fl <- RiboPilot:::itsx_check_flags("S1", strrep("A", 15120), ann, mk_its(pos_clean),
    RiboPilot:::parse_itsx_notes(pos_clean), RiboPilot:::parse_itsx_results(res_clean), no_prob,
    lineage = "Metazoa;Chordata")
  expect_equal(fl$flags, character(0))
  expect_equal(fl$note, "")
})

test_that("chimera is flagged from the positions note, the chimeric FASTA, or the results column", {
  its <- mk_its(c(pos_chim, pos_clean))
  notes <- RiboPilot:::parse_itsx_notes(c(pos_chim, pos_clean))
  expect_true("itsx_chimeric" %in% RiboPilot:::itsx_check_flags("S2", NA, NULL, its, notes, no_res, no_prob)$flags)
  # domain-order chimera: only in the chimeric FASTA / results, no note
  res <- RiboPilot:::parse_itsx_results(res_reorder)
  expect_true("itsx_chimeric" %in% RiboPilot:::itsx_check_flags("S4", NA, NULL, its, notes, res, no_prob)$flags)
  expect_true("itsx_chimeric" %in% RiboPilot:::itsx_check_flags("S1", NA, NULL, its, notes, no_res, no_prob,
    chimeric_ids = "S1")$flags)
  expect_false("itsx_chimeric" %in% RiboPilot:::itsx_check_flags("S1", NA, NULL, its, notes, no_res, no_prob)$flags)
})

test_that("missing 5.8S and over-long ITS map to misassembly flags", {
  its <- mk_its(pos_no58)
  fl <- RiboPilot:::itsx_check_flags("S3", NA, NULL, its, RiboPilot:::parse_itsx_notes(pos_no58), no_res, no_prob)
  expect_true(all(c("itsx_no_5.8S", "itsx_ITS_too_long") %in% fl$flags))
  expect_match(fl$note, "itsx_no_5.8S;itsx_ITS_too_long", fixed = TRUE)
})

test_that("an N-run inside the ITSx ITS span is a gap flag; outside it is not", {
  its <- mk_its(pos_clean)
  notes <- RiboPilot:::parse_itsx_notes(pos_clean)
  inside <- paste0(strrep("A", 7000), strrep("N", 30), strrep("A", 15120 - 7030))
  outside <- paste0(strrep("A", 100), strrep("N", 30), strrep("A", 15120 - 130))
  expect_true("itsx_gap_in_ITS" %in% RiboPilot:::itsx_check_flags("S1", inside, NULL, its, notes, no_res, no_prob)$flags)
  expect_false("itsx_gap_in_ITS" %in% RiboPilot:::itsx_check_flags("S1", outside, NULL, its, notes, no_res, no_prob)$flags)
})

test_that("ITS boundaries far from RiboPilot's annotation are reported per element", {
  its <- mk_its(pos_clean)
  notes <- RiboPilot:::parse_itsx_notes(pos_clean)
  ann <- data.frame(element = c("ITS1", "5.8S", "ITS2"),
    start = c(6833L, 7244L, 7500L), stop = c(7243L, 7401L, 7928L))  # ITS2 start off by 98
  fl <- RiboPilot:::itsx_check_flags("S1", NA, ann, its, notes, no_res, no_prob)
  expect_equal(fl$flags, "itsx_boundary_shift(ITS2:98bp)")
  ann$start[3] <- 7402L + 10L  # within tolerance
  expect_equal(RiboPilot:::itsx_check_flags("S1", NA, ann, its, notes, no_res, no_prob)$flags, character(0))
})

test_that("an ITSx origin outside the sample's lineage is flagged, unknown lineages are not", {
  its <- mk_its(pos_clean)
  notes <- RiboPilot:::parse_itsx_notes(pos_clean)
  res <- RiboPilot:::parse_itsx_results(sub("\tM\t", "\tF\t", res_clean, fixed = TRUE))
  fl <- RiboPilot:::itsx_check_flags("S1", NA, NULL, its, notes, res, no_prob, lineage = "Metazoa;Chordata")
  expect_equal(fl$flags, "itsx_origin(score: Fungi 33.15 vs Metazoa 49.9)")
  expect_equal(RiboPilot:::itsx_check_flags("S1", NA, NULL, its, notes, res, no_prob, lineage = "Bacteria")$flags,
    character(0))
})

test_that("an origin flag shows the scores and calls ties", {
  sc <- RiboPilot:::itsx_origin_scores
  expect_equal(sc("F: 1 5.4e-19 57.6, M: 1 5.4e-19 57.6", "F", "M"), "score: Fungi 57.6 vs Metazoa 57.6, tie")
  expect_equal(sc("F: 1 5e-19 57.6, M: 1 5e-18 56.2", "F", "M"), "score: Fungi 57.6 vs Metazoa 56.2, near tie")
  expect_equal(sc("F: 4 1e-20 60, M: 2 1e-08 30", "F", "M"), "score: Fungi 60 vs Metazoa 30")
  expect_equal(sc("F: 1 5e-19 57.6", "F", "M"), "")
  expect_equal(sc("", "F", "M"), "")
  res_ssu <- "S7\t1855\tF\t0\t1\t5.4e-19\t57.6\t14.4\t1811\t1855\t1_SSU\t***\t\t1_SSU\tF: 1 5.4e-19 57.6, M: 1 5.4e-19 57.6"
  pos_ssu <- "S7\t1855 bp.\tSSU: 1-1855\tITS1: Not found\t5.8S: Not found\tITS2: Not found\tLSU: Not found\tBroken or partial sequence, no 5.8S! "
  fl <- RiboPilot:::itsx_check_flags("S7", NA, NULL, mk_its(pos_ssu), RiboPilot:::parse_itsx_notes(pos_ssu),
    RiboPilot:::parse_itsx_results(res_ssu), no_prob, lineage = "Metazoa;Chordata")
  expect_true("itsx_origin(score: Fungi 57.6 vs Metazoa 57.6, tie)" %in% fl$flags)
})

test_that("a record ITSx did not detect is reported as such and nothing else", {
  fl <- RiboPilot:::itsx_check_flags("S9", "ACGT", NULL, mk_its(pos_clean),
    RiboPilot:::parse_itsx_notes(pos_clean), no_res, no_prob, no_detections = "S9")
  expect_equal(fl$flags, "itsx_not_detected")
  # absent from every ITSx file counts the same
  expect_equal(RiboPilot:::itsx_check_flags("S8", "ACGT", NULL, mk_its(pos_clean),
    RiboPilot:::parse_itsx_notes(pos_clean), no_res, no_prob)$flags, "itsx_not_detected")
})

# Minimal project for itsx_check_one
mk_itsx_project <- function() {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT,
                       problematic TEXT, partial TEXT, reviewed TEXT, chimeric INT,
                       n_features INT, missing TEXT, export_group TEXT)")
  seq <- paste0(strrep("A", 200), strrep("C", 200))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?,?,'linear',0)", params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,1,'no','no','yes',0,2,'',NULL)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',1,200)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',201,400)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Genus species',NULL,NULL,'Metazoa;Chordata')")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1)")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_origin TEXT")
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_scores TEXT")
  dir
}

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

test_that("singularity_image_ref turns a Docker Hub name into something singularity can open", {
  dir <- withr::local_tempdir()
  img <- "macguigand/ribopilot:0.0.0.9016"
  # nothing cached: pull from Docker Hub
  expect_equal(RiboPilot:::singularity_image_ref(img, dir, character(0)), paste0("docker://", img))
  # a URI or an existing file is left alone
  expect_equal(RiboPilot:::singularity_image_ref("docker://x/y:1", dir, character(0)), "docker://x/y:1")
  sif <- file.path(dir, "rp.sif"); file.create(sif)
  expect_equal(RiboPilot:::singularity_image_ref(sif, dir, character(0)), sif)
  # Nextflow's cached SIF under work/singularity wins
  wk <- file.path(dir, "work", "singularity"); dir.create(wk, recursive = TRUE)
  cached <- file.path(wk, "macguigand-ribopilot-0.0.0.9016.img"); file.create(cached)
  expect_equal(RiboPilot:::singularity_image_ref(img, dir, character(0)), cached)
  # explicit cacheDir in .config is checked instead
  cd <- file.path(dir, "cache"); dir.create(cd)
  cfg <- c("singularity {", "  enabled = true", paste0("  cacheDir = '", cd, "'"), "}")
  expect_equal(RiboPilot:::singularity_image_ref(img, dir, cfg), paste0("docker://", img))
  file.create(file.path(cd, "macguigand-ribopilot-0.0.0.9016.img"))
  expect_equal(RiboPilot:::singularity_image_ref(img, dir, cfg), file.path(cd, "macguigand-ribopilot-0.0.0.9016.img"))
})

test_that("itsx_check_one stores stock-ITSx flags on the unit's annotate row", {
  dir <- mk_itsx_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_warnings TEXT")
  pre <- file.path(dir, "unit.itsx_check")
  writeLines("S1.1.1\t400 bp.\tSSU: 1-200\tLSU: 201-400\tBroken or partial sequence, no 5.8S! ",
    paste0(pre, ".positions.txt"))
  writeLines("S1.1.1\t400\tM\t1\t2\tSSU\t\t\t\t\t\t\t\tSSU,LSU", paste0(pre, ".extraction.results"))
  writeLines("ITSx run finished at Thu Oct  8 12:00:00 2026.", paste0(pre, ".summary.txt"))
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, pre)
  expect_equal(out$itsx_check, "itsx_no_5.8S")
  expect_equal(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate")$itsx_warnings, "itsx_no_5.8S")
  # no ITSx output at all: recorded as not run, never an error
  out <- RiboPilot:::itsx_check_one(con, "S1", 1L, 1L, file.path(dir, "missing"))
  expect_equal(out$itsx_check, "not_run(ITSx did not finish)")
  expect_match(DBI::dbGetQuery(con, "SELECT itsx_warnings FROM annotate")$itsx_warnings, "^not_run\\(")
})

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
