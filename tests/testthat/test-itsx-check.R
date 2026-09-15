# Final pre-export ITSx check: parsers over stock-ITSx output, the pure flag
# function, and the export_bundle wiring.

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
  expect_equal(fl$flags, "itsx_origin(Fungi)")
  expect_equal(RiboPilot:::itsx_check_flags("S1", NA, NULL, its, notes, res, no_prob, lineage = "Bacteria")$flags,
    character(0))
})

test_that("a record ITSx did not detect is reported as such and nothing else", {
  fl <- RiboPilot:::itsx_check_flags("S9", "ACGT", NULL, mk_its(pos_clean),
    RiboPilot:::parse_itsx_notes(pos_clean), no_res, no_prob, no_detections = "S9")
  expect_equal(fl$flags, "itsx_not_detected")
  # absent from every ITSx file counts the same
  expect_equal(RiboPilot:::itsx_check_flags("S8", "ACGT", NULL, mk_its(pos_clean),
    RiboPilot:::parse_itsx_notes(pos_clean), no_res, no_prob)$flags, "itsx_not_detected")
})

test_that("run_itsx_check reports not-run when the binary is missing", {
  out <- RiboPilot:::run_itsx_check("x.fa", tempfile(), bin = "")
  expect_false(out$ok)
  expect_match(out$message, "not found")
})

# Minimal project for the export wiring (same shape as test-export-units.R)
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
  dir
}

test_that("export_bundle writes a not_run itsx_check column when ITSx is unavailable", {
  skip_if_not_installed("withr")
  dir <- mk_itsx_project()
  withr::local_envvar(PATH = tempdir())   # no ITSx here
  withr::local_envvar(HOME = tempdir())   # and no conda env to fall back to
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

test_that("a real ITSx run over a genuine rDNA record produces a clean check", {
  skip_if_not_installed("withr")
  bin <- RiboPilot:::itsx_check_bin()
  skip_if(!nzchar(bin), "ITSx not available")
  proj <- Sys.getenv("RIBOPILOT_TEST_PROJECT", "")
  skip_if(!nzchar(proj), "RIBOPILOT_TEST_PROJECT not set")
  src <- file.path(proj, ".sqlite")
  skip_if(!file.exists(src), "test project database not present")
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
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT, problematic TEXT,
    partial TEXT, reviewed TEXT, chimeric INT, n_features INT, missing TEXT, export_group TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('REAL','Fish sp.','Metazoa;Chordata')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('REAL',1,1,?,?,'linear',0)", params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('REAL',1,1,1,'no','no','yes',0,5,'',NULL)")
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
  expect_equal(chk$itsx_origin, "Metazoa")
  expect_match(chk$itsx_ITS1, "^[0-9]+-[0-9]+$")
  expect_true(file.exists(file.path(g, "itsx_check", "itsx.positions.txt")))
})
