# Persisting clade groups, wiring a built database into assemble_opts, and
# finishing background builds.

mk_ref_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT, ref_id TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT PRIMARY KEY, cpus INT,
    memory INT, getOrganelle TEXT, seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INT,
    max_scaffolds INT, min_assembly_length INT, assembler TEXT, bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, NULL, 1000, 25000)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default',4,16,'--flag','old_s.fa','old_l.fa',NULL,10,10,1000,'getorganelle','--local','--isolate')")
  for (id in c("S1", "S2")) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,NULL,NULL,NULL,1,2,NULL)", params = list(id))
    DBI::dbExecute(con, "INSERT INTO assemble VALUES(?,2,0,'default','default',NULL)", params = list(id))
  }
  con
}

grp <- data.frame(ref_id = "asteroidea", clade_taxid = 7588L, clade_name = "Asteroidea",
  clade_rank = "class", n_records = 4102L, n_samples = 2L, stringsAsFactors = FALSE)
mem <- data.frame(ID = c("S1", "S2"), ref_id = "asteroidea", stringsAsFactors = FALSE)

test_that("write_reference_groups stores groups and stamps membership", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(nrow(r), 1L)
  expect_equal(r$clade_name[1], "Asteroidea")
  expect_equal(r$reference_switch[1], 1L)
  expect_equal(r$query[1], "\"Asteroidea\"[Organism]")
  expect_setequal(DBI::dbGetQuery(con, "SELECT ref_id FROM samples")$ref_id, "asteroidea")
})

test_that("re-suggesting keeps a built group's database and state", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  RiboPilot:::write_reference_groups(con, grp, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 2L)
  expect_equal(r$seeds_db[1], "s.fa")
})

test_that("re-suggesting resets a group whose clade changed", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  moved <- grp; moved$clade_taxid <- 7586L
  RiboPilot:::write_reference_groups(con, moved, mem)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 1L)
  expect_true(is.na(r$seeds_db[1]))
})

test_that("wire_reference_db creates a preset and moves the group's samples", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  n <- RiboPilot:::wire_reference_db(con, "asteroidea", "/refs/a_seeds.fa", "/refs/a_labels.fa")
  expect_equal(n, 2L)
  o <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'asteroidea'")
  expect_equal(nrow(o), 1L)
  expect_equal(o$seeds_db[1], "/refs/a_seeds.fa")
  expect_equal(o$labels_db[1], "/refs/a_labels.fa")
  expect_equal(o$cpus[1], 4)              # copied from the samples' previous preset
  expect_equal(o$getOrganelle[1], "--flag")
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts,
    "asteroidea")
  expect_equal(DBI::dbGetQuery(con, "SELECT reference_switch FROM reference")$reference_switch[1], 2L)
})

test_that("wiring re-queues unlocked samples and leaves locked ones alone", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'S2'")
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa")
  a <- DBI::dbGetQuery(con, "SELECT ID, assemble_switch, assemble_opts FROM assemble ORDER BY ID")
  expect_equal(a$assemble_switch[a$ID == "S1"], 1L)
  expect_equal(a$assemble_switch[a$ID == "S2"], 2L)
  expect_equal(a$assemble_opts[a$ID == "S2"], "default")
})

test_that("skip leaves the samples' preset untouched", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::skip_reference_group(con, "asteroidea")
  r <- DBI::dbGetQuery(con, "SELECT reference_switch, source FROM reference")
  expect_equal(r$reference_switch[1], 5L)
  expect_equal(r$source[1], "skipped")
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("set_reference_clade splits a subset into its own group", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::set_reference_clade(con, ref_id = NULL, clade_taxid = 41692L,
    clade_name = "Forcipulatida", clade_rank = "order", n_records = 900L,
    member_ids = "S1")
  r <- DBI::dbGetQuery(con, "SELECT ref_id, n_samples FROM reference ORDER BY ref_id")
  expect_setequal(r$ref_id, c("asteroidea", "forcipulatida"))
  expect_equal(r$n_samples[r$ref_id == "forcipulatida"], 1L)
  expect_equal(r$n_samples[r$ref_id == "asteroidea"], 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT ref_id FROM samples WHERE ID = 'S1'")$ref_id, "forcipulatida")
})

test_that("start_reference_build records the job and poll finishes it", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = refs)
  RiboPilot:::write_reference_groups(con, grp, mem)

  seeds <- file.path(refs, "seeds", "asteroidea_rrna_seeds.fasta")
  labels <- file.path(refs, "labels", "asteroidea_rrna_labels.fasta")
  launcher <- function(cmd, args, log) {
    dir.create(dirname(seeds), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(labels), recursive = TRUE, showWarnings = FALSE)
    writeLines(">x", seeds); writeLines(">18S", labels)
    writeLines("done", log)
    4242L
  }
  RiboPilot:::start_reference_build(con, proj, "asteroidea", launcher = launcher)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 4L)
  expect_equal(r$pid[1], 4242L)

  n <- RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) FALSE)
  expect_equal(n, 1L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 2L)
  expect_equal(r$seeds_db[1], seeds)
  expect_equal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble WHERE ID='S1'")$assemble_opts,
    "asteroidea")
  expect_equal(nrow(RiboPilot:::read_ref_manifest(refs)), 1L)
})

test_that("a build that produced no FASTAs is recorded as failed", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = withr::local_tempdir())
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::start_reference_build(con, proj, "asteroidea",
    launcher = function(cmd, args, log) { writeLines("esearch: command not found", log); 7L })
  RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) FALSE)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 3L)
  expect_true(nzchar(r$fail_reason[1]))
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("a still-running build is left alone", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::set_ref_config(con, ref_dir = withr::local_tempdir())
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::start_reference_build(con, withr::local_tempdir(), "asteroidea",
    launcher = function(cmd, args, log) { writeLines("working", log); 11L })
  expect_equal(RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) TRUE), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT reference_switch FROM reference")$reference_switch[1], 4L)
})

test_that("reference_pending lists everything not done or skipped", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  expect_equal(RiboPilot:::reference_pending(con), "asteroidea")
  RiboPilot:::skip_reference_group(con, "asteroidea")
  expect_equal(length(RiboPilot:::reference_pending(con)), 0L)
})
