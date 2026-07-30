# The assembly fallback is a per-sample assembler choice, not a workflow. The
# load-bearing part is that repointing one sample must not change the assembler
# for every other sample sharing its options preset.

mk_fallback_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (
    assemble_opts TEXT NOT NULL, cpus INTEGER, memory INTEGER, getOrganelle TEXT,
    seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INTEGER,
    max_scaffolds INTEGER, min_assembly_length INTEGER,
    assembler TEXT DEFAULT 'getorganelle', PRIMARY KEY (assemble_opts))")
  DBI::dbExecute(con, "INSERT INTO assemble_opts
    (assemble_opts, cpus, memory, getOrganelle, seeds_db, assembler)
    VALUES ('default', 4, 16, '-F anonym', '/seeds.fa', 'getorganelle')")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, missing TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('FAILED',3,0,'default')")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('GOOD',2,0,'default')")
  con
}

test_that("set_sample_assembler derives a preset instead of editing the shared one", {
  con <- mk_fallback_db(); on.exit(DBI::dbDisconnect(con))

  nm <- RiboPilot:::set_sample_assembler(con, "FAILED", "bait")
  expect_equal(nm, "default_bait")

  # the derived preset carries the new assembler and copies everything else
  d <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'default_bait'")
  expect_equal(nrow(d), 1L)
  expect_equal(d$assembler, "bait")
  expect_equal(d$cpus, 4L)
  expect_equal(d$getOrganelle, "-F anonym")
  expect_equal(d$seeds_db, "/seeds.fa")

  # the SOURCE preset is untouched: this is the whole point
  src <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'default'")
  expect_equal(src$assembler, "getorganelle")

  # only the targeted sample moves
  a <- DBI::dbGetQuery(con, "SELECT ID, assemble_opts FROM assemble ORDER BY ID")
  expect_equal(a$assemble_opts[a$ID == "FAILED"], "default_bait")
  expect_equal(a$assemble_opts[a$ID == "GOOD"], "default")
})

test_that("set_sample_assembler is idempotent and does not duplicate the preset", {
  con <- mk_fallback_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_sample_assembler(con, "FAILED", "bait")
  RiboPilot:::set_sample_assembler(con, "FAILED", "bait")
  n <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemble_opts WHERE assemble_opts = 'default_bait'")$n
  expect_equal(n, 1L)
})

test_that("bait and denovo get separate presets", {
  con <- mk_fallback_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_sample_assembler(con, "FAILED", "bait")
  RiboPilot:::set_sample_assembler(con, "GOOD", "denovo")
  got <- DBI::dbGetQuery(con,
    "SELECT assemble_opts, assembler FROM assemble_opts ORDER BY assemble_opts")
  expect_equal(got$assemble_opts, c("default", "default_bait", "default_denovo"))
  expect_equal(got$assembler, c("getorganelle", "bait", "denovo"))
})

test_that("fallback_candidates finds QC-gate failures and missing anchors", {
  con <- mk_fallback_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::fallback_candidates(con), "FAILED")   # assemble_switch = 3

  DBI::dbExecute(con, "INSERT INTO annotate VALUES('GOOD','18S')")
  expect_setequal(RiboPilot:::fallback_candidates(con), c("FAILED", "GOOD"))

  # a locked sample is never a candidate
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  expect_length(RiboPilot:::fallback_candidates(con), 0L)
})

test_that("nf_command no longer offers a fallback workflow", {
  expect_error(RiboPilot:::nf_command(tempdir(), "fallback"), "arg")
  expect_false(any(grepl("--method", RiboPilot:::nf_command(tempdir(), "assemble"), fixed = TRUE)))
})
