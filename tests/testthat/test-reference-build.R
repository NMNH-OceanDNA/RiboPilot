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

test_that("one build run covers every group and shares its pid", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = withr::local_tempdir())
  two <- rbind(grp, data.frame(ref_id = "gastropoda", clade_taxid = 6448L,
    clade_name = "Gastropoda", clade_rank = "class", n_records = 5000L, n_samples = 1L,
    stringsAsFactors = FALSE))
  RiboPilot:::write_reference_groups(con, two,
    data.frame(ID = c("S1", "S2"), ref_id = c("asteroidea", "gastropoda"),
      stringsAsFactors = FALSE))

  seen <- new.env()
  launcher <- function(cmd, args, log) { seen$cmd <- cmd; seen$args <- args; 4242L }
  RiboPilot:::start_reference_build(con, proj, c("asteroidea", "gastropoda"),
    launcher = launcher)
  r <- DBI::dbGetQuery(con, "SELECT ref_id, pid, reference_switch, log_file FROM reference ORDER BY ref_id")
  expect_equal(nrow(r), 2L)
  expect_true(all(r$pid == 4242L))
  expect_true(all(r$reference_switch == 4L))
  expect_equal(length(unique(r$log_file)), 1L)
  expect_equal(seen$cmd, "nextflow")
  expect_true("-entry" %in% seen$args)
  expect_equal(seen$args[which(seen$args == "-entry") + 1], "BUILD")
  man <- utils::read.delim(seen$args[which(seen$args == "--clades") + 1], stringsAsFactors = FALSE)
  expect_setequal(man$ref_id, c("asteroidea", "gastropoda"))
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

test_that("wire_reference_db errors when no base preset exists, without moving samples", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  DBI::dbExecute(con, "DELETE FROM assemble_opts WHERE assemble_opts = 'default'")
  expect_error(RiboPilot:::wire_reference_db(con, "asteroidea", "s.fa", "l.fa"))
  a <- DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")
  expect_setequal(a$assemble_opts, "default")
  r <- DBI::dbGetQuery(con, "SELECT reference_switch FROM reference")
  expect_equal(r$reference_switch[1], 1L)
})

test_that("stale FASTAs predating the launch are not mistaken for a fresh build", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  proj <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = refs)
  RiboPilot:::write_reference_groups(con, grp, mem)

  seeds <- file.path(refs, "seeds", "asteroidea_rrna_seeds.fasta")
  labels <- file.path(refs, "labels", "asteroidea_rrna_labels.fasta")
  dir.create(dirname(seeds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(labels), recursive = TRUE, showWarnings = FALSE)
  writeLines(">old", seeds); writeLines(">old18S", labels)
  old_time <- Sys.time() - 3600
  Sys.setFileTime(seeds, old_time); Sys.setFileTime(labels, old_time)

  # overwrite: without it the leftovers would (correctly) be reused instead of built
  RiboPilot:::start_reference_build(con, proj, "asteroidea", overwrite = TRUE,
    launcher = function(cmd, args, log) { writeLines("build failed before writing output", log); 99L })
  n <- RiboPilot:::poll_reference_builds(con, alive_fn = function(pid) FALSE)
  expect_equal(n, 1L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 3L)
  expect_true(nzchar(r$fail_reason[1]))
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("set_reference_clade with no members creates a group that survives the count refresh", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  key <- RiboPilot:::set_reference_clade(con, ref_id = NULL, clade_taxid = 41692L,
    clade_name = "Forcipulatida", clade_rank = "order", n_records = 900L)
  expect_equal(key, "forcipulatida")
  r <- DBI::dbGetQuery(con, "SELECT ref_id, n_samples FROM reference WHERE ref_id = 'forcipulatida'")
  expect_equal(nrow(r), 1L)
  expect_equal(r$n_samples[1], 0L)
})

test_that("an already-built database is reused instead of launched", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  refs <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = refs)
  RiboPilot:::write_reference_groups(con, grp, mem)
  p <- RiboPilot:::reference_db_paths(refs, "asteroidea")
  dir.create(dirname(p$seeds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(p$labels), recursive = TRUE, showWarnings = FALSE)
  writeLines(">x", p$seeds); writeLines(">18S", p$labels)

  launched <- FALSE
  RiboPilot:::start_reference_build(con, withr::local_tempdir(), "asteroidea",
    launcher = function(cmd, args, log) { launched <<- TRUE; 1L })
  expect_false(launched)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 2L)
  expect_equal(r$source[1], "existing")
  expect_equal(r$seeds_db[1], p$seeds)
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts,
    "asteroidea")
})

test_that("overwrite launches a build even when the database exists", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  refs <- withr::local_tempdir()
  RiboPilot:::set_ref_config(con, ref_dir = refs)
  RiboPilot:::write_reference_groups(con, grp, mem)
  p <- RiboPilot:::reference_db_paths(refs, "asteroidea")
  dir.create(dirname(p$seeds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(p$labels), recursive = TRUE, showWarnings = FALSE)
  writeLines(">x", p$seeds); writeLines(">18S", p$labels)

  launched <- FALSE
  RiboPilot:::start_reference_build(con, withr::local_tempdir(), "asteroidea", overwrite = TRUE,
    launcher = function(cmd, args, log) { launched <<- TRUE; writeLines("go", log); 55L })
  expect_true(launched)
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(r$reference_switch[1], 4L)
  expect_equal(r$pid[1], 55L)
})

test_that("changing a group's clade re-keys the group and its samples", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  key <- RiboPilot:::set_reference_clade(con, ref_id = "asteroidea", clade_taxid = 7586L,
    clade_name = "Echinodermata", clade_rank = "phylum", n_records = 9000L,
    member_ids = c("S1", "S2"))
  expect_equal(key, "echinodermata")
  r <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  expect_equal(nrow(r), 1L)
  expect_equal(r$ref_id[1], "echinodermata")
  expect_equal(r$clade_name[1], "Echinodermata")
  expect_equal(r$clade_taxid[1], 7586L)
  expect_setequal(DBI::dbGetQuery(con, "SELECT ref_id FROM samples")$ref_id, "echinodermata")

  # re-applying the same clade keeps the key it already has
  expect_equal(RiboPilot:::set_reference_clade(con, ref_id = "echinodermata",
    clade_taxid = 7586L, clade_name = "Echinodermata", clade_rank = "phylum"), "echinodermata")
})

test_that("wiring never takes over a preset the user already made", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('asteroidea',8,32,'--mine','mine_s.fa','mine_l.fa',NULL,10,10,1000,'getorganelle','--local','--isolate')")
  RiboPilot:::write_reference_groups(con, grp, mem)
  RiboPilot:::wire_reference_db(con, "asteroidea", "/refs/s.fa", "/refs/l.fa")

  mine <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'asteroidea'")
  expect_equal(mine$seeds_db[1], "mine_s.fa")
  made <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'asteroidea_ref'")
  expect_equal(nrow(made), 1L)
  expect_equal(made$seeds_db[1], "/refs/s.fa")
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts,
    "asteroidea_ref")
  expect_equal(DBI::dbGetQuery(con, "SELECT preset_name FROM reference")$preset_name[1],
    "asteroidea_ref")

  # re-wiring reuses the stored name rather than inventing another one
  RiboPilot:::wire_reference_db(con, "asteroidea", "/refs/s2.fa", "/refs/l2.fa")
  names_now <- DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble_opts")$assemble_opts
  expect_setequal(names_now, c("default", "asteroidea", "asteroidea_ref"))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT seeds_db FROM assemble_opts WHERE assemble_opts = 'asteroidea_ref'")$seeds_db[1],
    "/refs/s2.fa")
  expect_equal(DBI::dbGetQuery(con, "SELECT seeds_db FROM assemble_opts WHERE assemble_opts = 'asteroidea'")$seeds_db[1],
    "mine_s.fa")
})

test_that("use_existing_reference errors when a file is missing", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_reference_groups(con, grp, mem)
  tmp <- withr::local_tempdir()
  seeds <- file.path(tmp, "seeds.fa")
  writeLines(">x", seeds)
  labels <- file.path(tmp, "labels.fa")  # not created
  expect_error(RiboPilot:::use_existing_reference(con, "asteroidea", seeds, labels))
  expect_setequal(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble")$assemble_opts, "default")
})

test_that("use_existing_reference wires an existing database and leaves a locked sample alone", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'S2'")
  RiboPilot:::write_reference_groups(con, grp, mem)
  tmp <- withr::local_tempdir()
  seeds <- file.path(tmp, "seeds.fa"); writeLines(">x", seeds)
  labels <- file.path(tmp, "labels.fa"); writeLines(">18S", labels)
  n <- RiboPilot:::use_existing_reference(con, "asteroidea", seeds, labels)
  expect_equal(n, 1L)
  r <- DBI::dbGetQuery(con, "SELECT source, seeds_db FROM reference")
  expect_equal(r$source[1], "existing")
  expect_equal(r$seeds_db[1], seeds)
  a <- DBI::dbGetQuery(con, "SELECT ID, assemble_opts FROM assemble ORDER BY ID")
  expect_equal(a$assemble_opts[a$ID == "S1"], "asteroidea")
  expect_equal(a$assemble_opts[a$ID == "S2"], "default")
})
