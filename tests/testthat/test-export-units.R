# Export at (ID, path, scaffold) grain: unit selection reaches the bundle, and
# SeqIDs do not change between a partial and a whole-project export.

# Minimal project: S1 has two exportable units, S2 has one.
mk_project <- function() {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT,
                       problematic TEXT, partial TEXT, reviewed TEXT, chimeric INT,
                       n_features INT, missing TEXT, export_group TEXT)")
  seq <- paste0(strrep("A", 200), strrep("C", 200))
  for (u in list(c("S1", 1, 1), c("S1", 1, 2), c("S2", 1, 1))) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES(?,?,?,?,?,'linear',0)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3]), seq, nchar(seq)))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES(?,?,?,1,'no','no','yes',0,2,'',NULL)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'18S','rRNA',1,200)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'28S','rRNA',201,400)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
  }
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Genus speciesA',NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL)")
  dir
}

bundle <- function(dir, ...) {
  suppressWarnings(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = FALSE, report = FALSE, zip = FALSE, itsx_check = FALSE, ...))
}

test_that("export_bundle honours a per-unit selection", {
  dir <- mk_project()
  # one scaffold of a two-scaffold sample: the sibling must not come along
  g <- bundle(dir, units = data.frame(ID = "S1", path = 1L, scaffold = 2L))
  fsa <- basename(Sys.glob(file.path(g, "samples", "*.fasta")))
  expect_equal(fsa, "S1_p1_s2.fasta")

  # sample-grain ids still export every unit of the sample
  g2 <- bundle(dir, ids = "S1", label = "byid")
  expect_setequal(basename(Sys.glob(file.path(g2, "samples", "*.fasta"))),
    c("S1_p1_s1.fasta", "S1_p1_s2.fasta"))
})

test_that("SeqIDs are stable between a partial and a whole-project export", {
  dir <- mk_project()
  one <- bundle(dir, units = data.frame(ID = "S1", path = 1L, scaffold = 1L), label = "one")
  all <- bundle(dir, label = "all")
  # S1 has two exportable units in the DB, so it is suffixed in BOTH bundles.
  # Counting within the filtered subset instead would emit a bare "S1" here.
  expect_true(file.exists(file.path(one, "samples", "S1_p1_s1.fasta")))
  expect_true(file.exists(file.path(all, "samples", "S1_p1_s1.fasta")))
  # S2 has a single unit, so it keeps its bare ID
  expect_true(file.exists(file.path(all, "samples", "S2.fasta")))
})

test_that("ignored and problematic units are excluded from the bundle", {
  dir <- mk_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE ID='S1' AND scaffold=2")
  DBI::dbDisconnect(con)
  g <- bundle(dir, label = "kept")
  # S1 now has one exportable unit, so it loses the suffix too
  expect_setequal(basename(Sys.glob(file.path(g, "samples", "*.fasta"))), c("S1.fasta", "S2.fasta"))
})

# Add a screen table carrying a paralog pair for S1, as screen_sample() would.
add_screen_paralogs <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version,
       variant_class,variant_group,variant_identity)
     VALUES ('S1',1,1,'on_target',3,'paralog','g1',NULL)")
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version,
       variant_class,variant_group,variant_identity)
     VALUES ('S1',1,2,'on_target',3,'paralog','g2',92.51)")
  invisible(dir)
}

test_that("export records the variant fact for a paralog and does NOT gate", {
  # a paralogous unit still exports: the fact is recorded, the choice stays the
  # curator's (doc 12 section 7). Warning on export is the behaviour this
  # sub-project exists to remove, so any warning here is a failure.
  dir <- mk_project()
  add_screen_paralogs(dir)
  warns <- character(0)
  g <- withCallingHandlers(
    suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
      markers = FALSE, report = FALSE, zip = FALSE, itsx_check = FALSE, ids = "S1")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_equal(warns, character(0))

  smry <- utils::read.csv(Sys.glob(file.path(g, "*_summary.csv")), stringsAsFactors = FALSE)
  expect_true(all(c("variant_class", "variant_group", "variant_identity", "n_variants")
    %in% names(smry)))
  expect_equal(sort(smry$variant_class), c("paralog", "paralog"))
  expect_setequal(smry$variant_group, c("g1", "g2"))
  expect_equal(unique(smry$n_variants), 2L)
  # both units are still in the bundle
  expect_setequal(basename(Sys.glob(file.path(g, "samples", "*.fasta"))),
    c("S1_p1_s1.fasta", "S1_p1_s2.fasta"))
  # and the variant fact never leaks into the submission-risk gate
  # (read.csv turns the empty risk string into NA)
  risk <- as.character(smry$submission_risk)
  expect_true(all(is.na(risk) | !nzchar(risk)))
})

test_that("n_variants counts comparable units, not every screen row", {
  # S1 unit 1 is an on-target singleton; unit 2 is a contaminant with no variant
  # class. n_variants must be 1, which counting screen rows would report as 2.
  dir <- mk_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version,variant_class,variant_group)
     VALUES ('S1',1,1,'on_target',3,'singleton','g1')")
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version)
     VALUES ('S1',1,2,'off_target',3)")
  DBI::dbDisconnect(con)

  g <- bundle(dir, ids = "S1", label = "counts")
  smry <- utils::read.csv(Sys.glob(file.path(g, "*_summary.csv")), stringsAsFactors = FALSE)
  expect_equal(unique(smry$n_variants), 1L)
  expect_equal(smry$variant_class[smry$scaffold == 1], "singleton")
  expect_true(is.na(smry$variant_class[smry$scaffold == 2]))
})

test_that("bundle file contents are exactly as specified", {
  dir <- mk_project()
  g <- bundle(dir, ids = "S2", label = "exact")
  seq <- paste0(strrep("A", 200), strrep("C", 200))
  expect_equal(readLines(file.path(g, "samples", "S2.fasta")),
    c(paste(">S2 [organism=Genus speciesB] [moltype=genomic DNA] Genus speciesB",
        "18S ribosomal RNA gene and 28S ribosomal RNA gene, partial sequence."), seq))
  expect_equal(readLines(file.path(g, "samples", "S2.tbl")),
    c(">Feature S2", "<1\t200\trRNA", "\t\t\tproduct\t18S ribosomal RNA",
      "201\t>400\trRNA", "\t\t\tproduct\t28S ribosomal RNA"))
  expect_equal(readLines(file.path(g, "exact.fasta")),
    readLines(file.path(g, "samples", "S2.fasta")))
  expect_equal(readLines(file.path(g, "exact_summary.csv")),
    c(paste0("\"seqid\",\"ID\",\"path\",\"scaffold\",\"organism\",\"length\",\"n_features\",",
        "\"missing\",\"reviewed\",\"partial\",\"problematic\",\"variant_class\",",
        "\"variant_group\",\"variant_identity\",\"n_variants\",\"submission_risk\""),
      "\"S2\",\"S2\",1,1,\"Genus speciesB\",400,2,\"\",\"yes\",\"no\",\"no\",NA,NA,NA,NA,\"\""))
})

test_that("an empty Taxon is flagged, and never rendered as the word NA", {
  dir <- mk_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE samples SET Taxon = NULL WHERE ID = 'S2'")
  DBI::dbExecute(con, "UPDATE samples SET Taxon = '  ' WHERE ID = 'S1'")
  DBI::dbDisconnect(con)

  expect_warning(g <- suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
    markers = TRUE, report = FALSE, zip = FALSE, itsx_check = FALSE, label = "notax")), "missing_organism")
  # glue renders an NA field as the characters "NA", which reads like a real
  # organism name. The header template cannot drop a field, so the value is
  # blanked and the unit is flagged instead.
  fsa <- readLines(file.path(g, "samples", "S2.fasta"))
  expect_match(fsa[1], "[organism=]", fixed = TRUE)
  expect_false(any(grepl("organism=NA", readLines(file.path(g, "notax.fasta")), fixed = TRUE)))
  expect_false(any(grepl("organism=NA",
    readLines(Sys.glob(file.path(g, "markers", "18S", "S2_18S.fasta"))), fixed = TRUE)))

  smry <- utils::read.csv(file.path(g, "notax_summary.csv"), stringsAsFactors = FALSE)
  expect_true(all(grepl("missing_organism", smry$submission_risk)))
})

test_that("set_export_group is atomic", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, export_group TEXT,
                       CHECK (ID <> 'BAD' OR export_group IS NULL))")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('OK',1,1,NULL)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('BAD',1,1,NULL)")
  expect_error(RiboPilot:::set_export_group(con, c("OK", "BAD"), "g1"))
  # the first unit's assignment must not survive the second one's failure
  got <- DBI::dbGetQuery(con, "SELECT ID, export_group FROM annotate ORDER BY ID")
  expect_true(all(is.na(got$export_group)))

  RiboPilot:::set_export_group(con, c("OK"), "g1")
  expect_equal(DBI::dbGetQuery(con, "SELECT export_group FROM annotate WHERE ID='OK'")$export_group, "g1")
})

test_that("a project with no screen table still exports, with NA variant columns", {
  dir <- mk_project()
  g <- bundle(dir, ids = "S2", label = "noscreen")
  smry <- utils::read.csv(Sys.glob(file.path(g, "*_summary.csv")), stringsAsFactors = FALSE)
  expect_true("variant_class" %in% names(smry))
  expect_true(all(is.na(smry$variant_class)))
  expect_true(all(is.na(smry$n_variants)))
})

test_that("set_export_group takes one group per unit", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, export_group TEXT)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('A',1,1,NULL)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('B',1,1,NULL)")
  u <- data.frame(ID = c("A", "B"), path = 1L, scaffold = 1L)
  RiboPilot:::set_export_group(con, u, c("g-getorganelle", "g-bait"))
  got <- DBI::dbGetQuery(con, "SELECT ID, export_group FROM annotate ORDER BY ID")
  expect_equal(got$export_group, c("g-getorganelle", "g-bait"))
  expect_error(RiboPilot:::set_export_group(con, u, c("x", "y", "z")))
})

test_that("sample_assembler follows assemble_opts and defaults to getorganelle", {
  dir <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, assembler TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('A','default'),('B','default_bait'),('C','gone')")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default','getorganelle'),('default_bait','bait')")
  DBI::dbDisconnect(con)
  expect_equal(RiboPilot:::sample_assembler(dir, c("B", "A", "C", "D")),
    c("bait", "getorganelle", "getorganelle", "getorganelle"))
})

test_that("a group name that is not a safe folder name is refused before anything is deleted", {
  dir <- mk_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_error(RiboPilot:::set_export_group(con, "S1", "Batch 1"), "Try \"Batch_1\"", fixed = TRUE)
  for (g in c("..", ".", "../S1", "18S/ITS", "-x")) {
    expect_error(RiboPilot:::set_export_group(con, "S1", g), "cannot be used", info = g)
  }
  expect_true(all(is.na(DBI::dbGetQuery(con, "SELECT export_group FROM annotate")$export_group)))
  RiboPilot:::set_export_group(con, "S1", "batch_2026.01-a")
  # a bad label must not reach the recursive delete of out/export/<label>
  dir.create(file.path(dir, "out", "keep"), recursive = TRUE)
  expect_error(bundle(dir, label = ".."), "cannot be used")
  expect_true(dir.exists(file.path(dir, "out", "keep")))
})

test_that("the zip lands in out/export from a relative project path, and is replaced", {
  skip_if(!nzchar(Sys.which("zip")), "zip not on PATH")
  dir <- mk_project()
  withr::local_dir(dir)
  z <- function(...) suppressWarnings(suppressMessages(export_bundle(".", label = "z",
    markers = FALSE, report = FALSE, itsx_check = FALSE, zip = TRUE, ...)))
  z()
  zf <- file.path("out", "export", "z.zip")
  expect_true(file.exists(zf))
  n_all <- nrow(utils::unzip(zf, list = TRUE))
  z(ids = "S2")
  expect_lt(nrow(utils::unzip(zf, list = TRUE)), n_all)
})
