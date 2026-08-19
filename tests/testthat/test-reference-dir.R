# Reference-folder resolution order and the cross-project manifest.

mk_cfg_db <- function(ref_dir = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con,
    "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
     VALUES (1, ?, 1000, 25000)", params = list(ref_dir))
  con
}

test_that("the project setting wins over the environment variable", {
  d <- withr::local_tempdir()
  con <- mk_cfg_db(d); on.exit(DBI::dbDisconnect(con))
  withr::local_envvar(RIBOPILOT_REF_DIR = "/somewhere/else")
  expect_equal(RiboPilot:::ref_config(con)$ref_dir, path.expand(d))
})

test_that("an unset project setting falls back to the environment variable", {
  d <- withr::local_tempdir()
  con <- mk_cfg_db(NA_character_); on.exit(DBI::dbDisconnect(con))
  withr::local_envvar(RIBOPILOT_REF_DIR = d)
  expect_equal(RiboPilot:::ref_config(con)$ref_dir, path.expand(d))
})

test_that("set_ref_config writes only the arguments it is given", {
  con <- mk_cfg_db(NA_character_); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_ref_config(con, min_records = 200)
  cfg <- DBI::dbGetQuery(con, "SELECT * FROM reference_config")
  expect_equal(cfg$min_records[1], 200)
  expect_equal(cfg$max_records[1], 25000)
  expect_true(is.na(cfg$ref_dir[1]))
})

test_that("ref_config reports the stored panel prefix, empty when unset", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(ref_config(con)$blastdb, "")
  set_ref_config(con, blastdb = "/data/blastdb/2026.08.1/ribopilot_rrna")
  expect_equal(ref_config(con)$blastdb, "/data/blastdb/2026.08.1/ribopilot_rrna")
})

test_that("set_ref_config clears the panel prefix when given an empty string", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  set_ref_config(con, blastdb = "/data/p/ribopilot_rrna")
  set_ref_config(con, blastdb = "")
  expect_equal(ref_config(con)$blastdb, "")
  expect_true(is.na(DBI::dbGetQuery(con,
    "SELECT blastdb FROM reference_config WHERE rowid_key = 1")$blastdb[1]))
})

test_that("set_ref_config leaves the panel prefix alone when not supplied", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  set_ref_config(con, blastdb = "/data/p/ribopilot_rrna")
  set_ref_config(con, min_records = 500)
  expect_equal(ref_config(con)$blastdb, "/data/p/ribopilot_rrna")
})

test_that("validate_ref_dir rejects a path that is not a writable directory", {
  f <- withr::local_tempfile()
  writeLines("x", f)
  expect_false(RiboPilot:::validate_ref_dir(f)$ok)
  expect_true(RiboPilot:::validate_ref_dir(withr::local_tempdir())$ok)
})

test_that("the manifest round-trips a built database", {
  d <- withr::local_tempdir()
  expect_equal(nrow(RiboPilot:::read_ref_manifest(d)), 0L)
  RiboPilot:::append_ref_manifest(d, clade = "asteroidea", taxid = 7588L, rank = "class",
    query = "\"Asteroidea\"[Organism]", n_records = 4102L,
    seeds = file.path(d, "seeds", "asteroidea_rrna_seeds.fasta"),
    labels = file.path(d, "labels", "asteroidea_rrna_labels.fasta"))
  m <- RiboPilot:::read_ref_manifest(d)
  expect_equal(nrow(m), 1L)
  expect_equal(m$clade[1], "asteroidea")
  expect_equal(m$taxid[1], 7588L)
  expect_equal(m$n_records[1], 4102L)
})

test_that("appending onto a malformed manifest preserves it as a .bak and warns", {
  d <- withr::local_tempdir()
  p <- RiboPilot:::ref_manifest_path(d)
  garbage <- "this is not\ta valid manifest\njust noise"
  writeLines(garbage, p)
  expect_warning(
    RiboPilot:::append_ref_manifest(d, "asteroidea", 7588L, "class", "q", 10L, "s1", "l1"))
  bak <- paste0(p, ".bak")
  expect_true(file.exists(bak))
  expect_equal(readLines(bak), strsplit(garbage, "\n")[[1]])
  m <- RiboPilot:::read_ref_manifest(d)
  expect_equal(nrow(m), 1L)
  expect_equal(m$clade[1], "asteroidea")
})

test_that("appending twice keeps one row per clade, the newest", {
  d <- withr::local_tempdir()
  RiboPilot:::append_ref_manifest(d, "asteroidea", 7588L, "class", "q", 10L, "s1", "l1")
  RiboPilot:::append_ref_manifest(d, "asteroidea", 7588L, "class", "q", 20L, "s2", "l2")
  m <- RiboPilot:::read_ref_manifest(d)
  expect_equal(nrow(m), 1L)
  expect_equal(m$n_records[1], 20L)
  expect_equal(m$seeds[1], "s2")
})

test_that("an empty panel prefix is valid and means remote", {
  r <- validate_blastdb_prefix("")
  expect_true(r$ok); expect_equal(r$msg, "")
  expect_true(validate_blastdb_prefix(NULL)$ok)
})

test_that("a non-scalar prefix is rejected rather than erroring", {
  expect_false(validate_blastdb_prefix(c("a", "b"))$ok)
})

test_that("a prefix with no volume files is rejected", {
  r <- validate_blastdb_prefix(file.path(tempdir(), "nope", "ribopilot_rrna"))
  expect_false(r$ok)
  expect_match(r$msg, "No BLAST database")
})

test_that("a prefix with volume files is accepted without needing write access", {
  d <- tempfile(); dir.create(d)
  on.exit({Sys.chmod(d, "0755"); unlink(d, recursive = TRUE)})
  file.create(file.path(d, "ribopilot_rrna.nin"))
  Sys.chmod(d, "0555")  # read-only, as on a shared HPC mount
  skip_if(unname(file.access(d, 2)) == 0, "cannot drop write permission here")
  r <- validate_blastdb_prefix(file.path(d, "ribopilot_rrna"))
  expect_true(r$ok)
  expect_equal(r$msg, file.path(d, "ribopilot_rrna"))
})

test_that("a whitespace prefix is rejected, since blastn cannot express it", {
  d <- tempfile("has space"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  file.create(file.path(d, "ribopilot_rrna.nin"))
  r <- validate_blastdb_prefix(file.path(d, "ribopilot_rrna"))
  expect_false(r$ok)
  expect_match(r$msg, "space")
})
