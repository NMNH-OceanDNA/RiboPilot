# The reference tab's table shape and modal content.

# The wired paths are real files: a "Ready" row is only Ready while its FASTAs
# are on disk, so a fixture pointing at made-up paths would read as broken.
tab_db_paths <- function(dir = tempdir()) {
  p <- list(seeds = file.path(dir, "a_seeds.fa"), labels = file.path(dir, "a_labels.fa"))
  for (f in p) writeLines(c(">x", "ACGT"), f)
  p
}

mk_tab_db <- function(paths = tab_db_paths()) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, NULL, 1000, 25000, NULL, NULL, NULL)")
  DBI::dbExecute(con, "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank,
    query, n_samples, n_records, reference_switch, source, seeds_db, labels_db, time_stamp)
    VALUES ('asteroidea', 7588, 'Asteroidea', 'class', '\"Asteroidea\"[Organism]', 4, 4102,
      2, 'built', ?, ?, 1)", params = list(paths$seeds, paths$labels))
  con
}

test_that("reference_table exposes ref_id as ID and the seed path as database", {
  p <- tab_db_paths()
  con <- mk_tab_db(p); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::reference_table(con)
  expect_equal(d$ID[1], "asteroidea")
  expect_equal(d$clade_name[1], "Asteroidea")
  expect_equal(d$database[1], p$seeds)
  expect_true("reference_switch" %in% names(d))
  # the state and the row's verb are rendered as words, not as icons the user
  # has to decode; `source` was never explained anywhere, so it is gone
  expect_equal(d$status[1], "Ready")
  expect_equal(d$action[1], "Change")
  expect_false("source" %in% names(d))
})

test_that("each reference state gets its own wording", {
  con <- mk_tab_db(); on.exit(DBI::dbDisconnect(con))
  for (s in c(1L, 3L, 4L, 5L)) {
    DBI::dbExecute(con, "UPDATE reference SET reference_switch = ?", params = list(s))
    d <- RiboPilot:::reference_table(con)
    expect_equal(d$status[1], c(`1` = "Needs a database", `3` = "Build failed",
      `4` = "Building", `5` = "Left on current database")[[as.character(s)]])
  }
  # a skipped group is still changeable, and a failed one says so
  DBI::dbExecute(con, "UPDATE reference SET reference_switch = 3")
  expect_equal(RiboPilot:::reference_table(con)$action[1], "Try again")
  DBI::dbExecute(con, "UPDATE reference SET reference_switch = 1")
  expect_equal(RiboPilot:::reference_table(con)$action[1], "Set up database")
})

# The reference folder is shared and can be deleted long after a group was wired
# to it. The row used to keep reading "Ready" with a dead path, offering only
# "Change", and the pipeline then assembled nothing while reporting success.
test_that("a Ready group whose FASTAs are gone reads as missing, not Ready", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  p <- tab_db_paths(dir)
  con <- mk_tab_db(p); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::reference_table(con)$status[1], "Ready")
  unlink(p$seeds)
  d <- RiboPilot:::reference_table(con)
  expect_equal(d$status[1], "Database missing")
  expect_equal(d$action[1], "Set up database")
  expect_match(d$fail_reason[1], "no longer")
  # the labels FASTA counts too, and an empty file is as dead as a missing one
  writeLines(c(">x", "ACGT"), p$seeds); unlink(p$labels); file.create(p$labels)
  expect_equal(RiboPilot:::reference_table(con)$status[1], "Database missing")
  # the labels path never reaches the table: only the seed path is shown
  expect_false("labels_db" %in% names(RiboPilot:::reference_table(con)))
})

test_that("the ladder modal shows every rung with coverage and counts", {
  ladder <- data.frame(level = 2:5, rank = c("phylum", "class", "order", "family"),
    taxid = c(7586L, 7588L, 41692L, 7607L),
    name = c("Echinodermata", "Asteroidea", "Forcipulatida", "Asteriidae"),
    n_samples = c(4L, 4L, 3L, 3L), n_records = c(31500L, 4102L, 900L, 340L),
    stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_ladder_modal("asteroidea", "Asteroidea", ladder, 7588L, 4L)
  html <- as.character(m)
  expect_true(grepl("Asteriidae", html, fixed = TRUE))
  expect_true(grepl("Echinodermata", html, fixed = TRUE))
  expect_true(grepl("4,102", html, fixed = TRUE))
  expect_true(grepl("ref_rung", html, fixed = TRUE))
  expect_true(grepl("ref_search", html, fixed = TRUE))
})

test_that("the ladder modal marks the current clade as chosen", {
  ladder <- data.frame(level = 3L, rank = "class", taxid = 7588L, name = "Asteroidea",
    n_samples = 4L, n_records = 4102L, stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_ladder_modal("asteroidea", "Asteroidea", ladder, 7588L, 4L)
  expect_true(grepl("current", as.character(m), fixed = TRUE))
})

test_that("the existing-database modal lists manifest entries", {
  man <- data.frame(clade = "gastropoda", taxid = 6448L, rank = "class",
    query = "q", n_records = 5000L, built = "2026-07-01 09:00:00",
    seeds = "/refs/seeds/gastropoda_rrna_seeds.fasta",
    labels = "/refs/labels/gastropoda_rrna_labels.fasta", stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("gastropoda", html, fixed = TRUE))
  expect_true(grepl("5,000", html, fixed = TRUE))
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
})

test_that("the existing-database modal still offers typed paths with an empty manifest", {
  man <- RiboPilot:::read_ref_manifest(withr::local_tempdir())
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
  expect_true(grepl("ref_existing_labels", html, fixed = TRUE))
})

test_that("the existing-database modal renders every clade, even with a two-row manifest and a missing rank", {
  # %|NA|% (R/utils.R) is not vectorised; a manifest with more than one row used
  # to raise "the condition has length > 1" as soon as $rank or $n_records held
  # more than one value, or any NA.
  man <- data.frame(clade = c("gastropoda", "bivalvia"), taxid = c(6448L, 6544L),
    rank = c("class", NA_character_), query = c("q1", "q2"), n_records = c(5000L, NA_integer_),
    built = c("2026-07-01 09:00:00", "2026-07-02 10:00:00"),
    seeds = c("/refs/seeds/gastropoda_rrna_seeds.fasta", "/refs/seeds/bivalvia_rrna_seeds.fasta"),
    labels = c("/refs/labels/gastropoda_rrna_labels.fasta", "/refs/labels/bivalvia_rrna_labels.fasta"),
    stringsAsFactors = FALSE)
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("gastropoda", html, fixed = TRUE))
  expect_true(grepl("bivalvia", html, fixed = TRUE))
  expect_true(grepl("5,000", html, fixed = TRUE))
  expect_true(grepl("records unknown", html, fixed = TRUE))
  expect_false(grepl("NA records", html, fixed = TRUE))
})

# ---- reference_handlers(): the ladder-apply split behavior ----
#
# Two sea stars start in one "asteroidea" group. STAR (S1) is Forcipulatida /
# Asteriidae; STAR2 (S2) is Valvatida / Asterinidae -- they share only
# kingdom/phylum/class. Picking Forcipulatida (order, taxid 41692) covers only
# S1, so S2 must be split into a group of its own, not silently kept inside (or
# dumped out of) a clade it does not belong to.
#
# reference_handlers() is not yet mounted by app_server() (Task 7), so it is
# driven directly with shiny::testServer() and a stand-in `deps`. The chosen
# taxid's record count is pre-seeded into the reference folder's cache, so the
# cache-only lookup always answers.
mk_split_db <- function(dir) {
  ref_dir <- file.path(dir, "refs")
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  # pre-seeded so no count ever has to be fetched: 33208/7586/7588 are the ranks
  # S1 and S2 still share, 41692 is the Forcipulatida rung the tests pick, and
  # 7607/41695/41693 complete every rung at or above family that Suggest asks about
  utils::write.table(
    data.frame(taxid = c(33208L, 7586L, 7588L, 41692L, 7607L, 41695L, 41693L),
      n_records = c(900000L, 60000L, 4102L, 900L, 340L, 700L, 200L)),
    file.path(ref_dir, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  dbfile <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, ?, 1000, 25000, NULL, NULL, NULL)", params = list(ref_dir))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, Taxon TEXT, sci_name TEXT,
    lineage TEXT, lineage_taxids TEXT, ref_id TEXT, user_group TEXT, ref_source TEXT,
    time_stamp INTEGER)")
  DBI::dbExecute(con, "INSERT INTO samples (ID, lineage, lineage_taxids, ref_id) VALUES('S1',
    'Metazoa;Echinodermata;Asteroidea;Forcipulatida;Asteriidae;Asterias;Asterias rubens',
    '33208;7586;7588;41692;7607;7608;7610', 'asteroidea')")
  DBI::dbExecute(con, "INSERT INTO samples (ID, lineage, lineage_taxids, ref_id) VALUES('S2',
    'Metazoa;Echinodermata;Asteroidea;Valvatida;Asterinidae;Patiria;Patiria miniata',
    '33208;7586;7588;41695;41693;104852;46514', 'asteroidea')")
  DBI::dbExecute(con, "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank, query,
    n_samples, n_records, reference_switch, time_stamp)
    VALUES ('asteroidea', 7588, 'Asteroidea', 'class', '\"Asteroidea\"[Organism]', 2, 4102, 1, 1)")
  DBI::dbDisconnect(con)
  dbfile
}

mk_split_deps <- function(dbfile, dir) {
  list(con = function() DBI::dbConnect(RSQLite::SQLite(), dbfile),
    project = function() dir,
    refresh = shiny::reactiveVal(0),
    selected_ids = function() character(0),
    table_data = shiny::reactive(data.frame()))
}

test_that("a partial-coverage pick splits the uncovered sample into a new group", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)

  # reference_handlers() is not a module server (its formals are input, output,
  # session, deps -- not id, ...), so shiny::testServer() cannot take `args`;
  # `deps` is instead captured by the wrapper via closure. Binding open_ladder
  # as a LOCAL of the wrapper (not just returning it) is what makes it visible
  # to the test expression below, which testServer() evaluates against the
  # wrapper's own call frame.
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    open_ladder <- h$open_ladder
  }, {
    open_ladder("asteroidea")
    session$setInputs(ref_rung = "41692")   # Forcipulatida: covers only S1
    session$flushReact()
    session$setInputs(ref_clade_apply = 1)
    session$flushReact()
    # the apply runs on a later() tick so the "Applying..." alert can paint
    while (!later::loop_empty()) later::run_now()
  })

  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  s <- DBI::dbGetQuery(con, "SELECT ID, ref_id FROM samples ORDER BY ID")
  # the group is re-keyed on the clade it now holds, so nothing is left keyed
  # on the clade the user moved away from
  expect_equal(s$ref_id[s$ID == "S1"], "forcipulatida")
  # S2 is no longer inside the Forcipulatida group it does not belong to
  expect_false(s$ref_id[s$ID == "S2"] == "forcipulatida")

  ref <- DBI::dbGetQuery(con, "SELECT ref_id, clade_name, clade_rank, n_samples FROM reference ORDER BY ref_id")
  expect_equal(nrow(ref), 2L)
  expect_false("asteroidea" %in% ref$ref_id)
  orig <- ref[ref$ref_id == "forcipulatida", ]
  expect_equal(orig$clade_name, "Forcipulatida")
  expect_equal(orig$n_samples, 1L)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT query FROM reference WHERE ref_id = 'forcipulatida'")$query, "txid41692[Organism:exp]")
  split_key <- s$ref_id[s$ID == "S2"]
  spl <- ref[ref$ref_id == split_key, ]
  expect_equal(nrow(spl), 1L)
  expect_equal(spl$n_samples, 1L)
  # S2's own group is keyed on ITS shared lineage (deepest = itself), not on
  # the ref_id the edited group just took
  expect_false(spl$ref_id == "forcipulatida")
})

test_that("a clade covering none of the group's samples is refused, not applied to everyone", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    open_ladder <- h$open_ladder
  }, {
    open_ladder("asteroidea")
    session$setInputs(ref_rung = "999999999")   # not in either sample's lineage
    session$flushReact()
    session$setInputs(ref_clade_apply = 1)
    session$flushReact()
    # the apply runs on a later() tick so the "Applying..." alert can paint
    while (!later::loop_empty()) later::run_now()
  })

  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  s <- DBI::dbGetQuery(con, "SELECT ID, ref_id FROM samples ORDER BY ID")
  expect_equal(s$ref_id, c("asteroidea", "asteroidea"))   # untouched
  ref <- DBI::dbGetQuery(con, "SELECT clade_name, clade_taxid, n_samples FROM reference WHERE ref_id = 'asteroidea'")
  expect_equal(ref$clade_name, "Asteroidea")
  expect_equal(ref$clade_taxid, 7588L)
  expect_equal(ref$n_samples, 2L)
})

test_that("Suggest counts only uncached clades and then walks offline", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  calls <- 0L
  withr::local_options(list(ribopilot.count_runner = function(project, args, log = NULL) {
    calls <<- calls + 1L; 0L
  }))
  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_suggest_btn = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
    # Suggest only previews now; Apply is what writes
    session$setInputs(ref_diff_apply = 1)
    session$flushReact()
  })
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  ref <- DBI::dbGetQuery(con, "SELECT ref_id, time_stamp FROM reference")
  expect_gt(nrow(ref), 0L)
  # the fixture pre-inserts a row stamped 1; only a real write_reference_groups()
  # re-stamps it, so this fails if Suggest no-opped or threw into the error alert
  expect_true(all(ref$time_stamp > 1L))
  expect_equal(calls, 0L)   # every rung was already cached
})

test_that("Suggest sends the uncached rungs to the runner", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  # drop two rungs from the cache so a count run is required
  cache <- utils::read.delim(file.path(dir, "refs", "counts.tsv"))
  utils::write.table(cache[!(cache$taxid %in% c(7607L, 41693L)), , drop = FALSE],
    file.path(dir, "refs", "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  sent <- NULL
  withr::local_options(list(ribopilot.count_runner = function(project, args, log = NULL) {
    man <- args[which(args == "--clades") + 1L]
    sent <<- utils::read.delim(man)
    0L
  }))
  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_suggest_btn = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  expect_setequal(sent$taxid, c(7607L, 41693L))
  expect_true(all(grepl("^txid[0-9]+\\[Organism:exp\\]$", sent$query)))
  # the stub runner produced no count files, so the two dropped rungs are still
  # missing from the cache, which is what proves the walk stayed offline
  after <- RiboPilot:::read_counts_cache(file.path(dir, "refs"))
  expect_false(any(c(7607L, 41693L) %in% after$taxid))
})

test_that("opening the ladder modal never fetches a count for an uncached rung", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  # 7588 (Asteroidea) is a rung the modal's group_ladder() asks about; drop it
  cache <- utils::read.delim(file.path(dir, "refs", "counts.tsv"))
  utils::write.table(cache[cache$taxid != 7588L, , drop = FALSE],
    file.path(dir, "refs", "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    open_ladder <- h$open_ladder
  }, {
    open_ladder("asteroidea")
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  after <- RiboPilot:::read_counts_cache(file.path(dir, "refs"))
  expect_false(7588L %in% after$taxid)
})

test_that("apply-clade counts a cached clade without calling the runner", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  calls <- 0L
  withr::local_options(list(ribopilot.count_runner = function(project, args, log = NULL) {
    calls <<- calls + 1L; 0L
  }))
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    open_ladder <- h$open_ladder
  }, {
    open_ladder("asteroidea")
    session$setInputs(ref_rung = "41692")
    session$flushReact()
    session$setInputs(ref_clade_apply = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  expect_equal(calls, 0L)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  n <- DBI::dbGetQuery(con, "SELECT n_records FROM reference WHERE ref_id = 'forcipulatida'")
  expect_equal(n$n_records, 900L)
})

# ---- the row action: which groups a click acts on ----
#
# The action cell replaced a bar of verbs, so the rule that decides its scope is
# the whole point: a click on a row inside the current selection acts on the
# selection, a click anywhere else acts on that row alone. Neither case may be
# refused after the click, which is what the old "select exactly one" bar
# buttons did.
mk_two_group_deps <- function(dbfile, dir, selected) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile)
  DBI::dbExecute(con, "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank,
    query, n_samples, n_records, reference_switch, time_stamp)
    VALUES ('echinoidea', 7625, 'Echinoidea', 'class', '\"Echinoidea\"[Organism]', 0, 3110, 1, 1)")
  tbl <- RiboPilot:::reference_table(con)
  DBI::dbDisconnect(con)
  deps <- mk_split_deps(dbfile, dir)
  deps$table_data <- shiny::reactive(tbl)
  deps$selected_ids <- function() selected
  deps
}

skipped_ids <- function(dbfile) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "SELECT ref_id FROM reference WHERE reference_switch = 5 ORDER BY ref_id")$ref_id
}

test_that("clicking a row inside the selection acts on the whole selection", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_two_group_deps(dbfile, dir, selected = c("asteroidea", "echinoidea"))

  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_row_action = "1")   # row 1 = asteroidea, and it is ticked
    session$flushReact()
    session$setInputs(ref_act_skip = 1)
    session$flushReact()
  })

  expect_equal(skipped_ids(dbfile), c("asteroidea", "echinoidea"))
})

test_that("clicking a row outside the selection acts on that row alone", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_two_group_deps(dbfile, dir, selected = "echinoidea")

  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_row_action = "1")   # row 1 = asteroidea, NOT ticked
    session$flushReact()
    session$setInputs(ref_act_skip = 1)
    session$flushReact()
  })

  # the ticked-but-unclicked group is left alone rather than swept in
  expect_equal(skipped_ids(dbfile), "asteroidea")
})

# ---- reference_run_modal(): the Nextflow-run modal, mirroring run_modal() ----

test_that("the reference run modal carries the command, progress, and controls", {
  html <- as.character(RiboPilot:::reference_run_modal(2, "/lab/refs"))
  for (id in c("ref_api_key", "ref_refseq", "ref_nf_resume", "ref_nf_cmd_text",
               "ref_nf_prog_all", "ref_nf_start", "ref_nf_stop", "ref_nf_close",
               "ref_gears")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  expect_true(grepl("/lab/refs", html, fixed = TRUE))
  expect_true(grepl("2 clade(s) to build", html, fixed = TRUE))
})

test_that("the headless reference modal offers the job script instead of a local run", {
  html <- as.character(RiboPilot:::reference_run_modal(1, "/lab/refs", headless = TRUE,
    sched = list(name = "Slurm", cmd = "sbatch"), script = "#!/bin/bash\nnextflow run"))
  expect_true(grepl("ref_submit_job", html, fixed = TRUE))
  expect_true(grepl("ref_save_script", html, fixed = TRUE))
  expect_false(grepl("ref_nf_start", html, fixed = TRUE))
})

test_that("headless Build offers a BUILD job script and saves it, running nothing", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  withr::local_options(list(RiboPilot.headless = TRUE))

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    session$setInputs(ref_save_script = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })

  f <- Sys.glob(file.path(dir, "jobs", "ribopilot_ref_*.sh"))
  expect_true(file.exists(f))
  script <- paste(readLines(f), collapse = "\n")
  # -log is a top-level option, so it sits between the binary and `run`
  expect_match(script, "'nextflow' '-log'", fixed = TRUE)
  expect_match(script, "'run'", fixed = TRUE)
  expect_match(script, "'-entry' 'BUILD'", fixed = TRUE)
  expect_true(grepl("--clades", script, fixed = TRUE))
  # the job reads the manifest long after the modal closed, so it must be on disk
  man <- utils::read.delim(RiboPilot:::reference_manifest_path(dir), stringsAsFactors = FALSE)
  expect_equal(man$ref_id, "asteroidea")
  # nothing was launched: the group is still pending, with no pid
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  r <- DBI::dbGetQuery(con, "SELECT reference_switch, pid FROM reference WHERE ref_id = 'asteroidea'")
  expect_equal(r$reference_switch[1], 1L)
  expect_true(is.na(r$pid[1]))
})

test_that("the saved job script carries the API key and RefSeq settings typed in the modal", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  withr::local_options(list(RiboPilot.headless = TRUE))
  # the stale copy below is built out of band, so it only matches what the
  # handler generated if both get the same submission timestamp
  testthat::local_mocked_bindings(
    job_basename = function(...) "ribopilot_ref_2026-08-12_09-45-03", .package = "RiboPilot")

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    # a real browser holds the script the modal was rendered with, generated
    # before the settings inputs existed; that stale copy must not be saved
    session$setInputs(ref_submit_script_text = RiboPilot:::build_reference_submit_script(
      dir, RiboPilot:::reference_manifest_path(dir), file.path(dir, "refs"),
      sched = RiboPilot:::detect_scheduler()))
    session$setInputs(ref_api_key = "KEY123", ref_refseq = TRUE)
    session$flushReact()
    session$setInputs(ref_save_script = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })

  script <- paste(readLines(Sys.glob(file.path(dir, "jobs", "ribopilot_ref_*.sh"))), collapse = "\n")
  expect_match(script, "'--ncbi_api_key' 'KEY123'", fixed = TRUE)
  expect_match(script, "'--refseq' 'true'", fixed = TRUE)
})

test_that("an emptied script textarea falls back to a generated script", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  withr::local_options(list(RiboPilot.headless = TRUE))

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    session$setInputs(ref_submit_script_text = "")
    session$setInputs(ref_save_script = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })

  script <- paste(readLines(Sys.glob(file.path(dir, "jobs", "ribopilot_ref_*.sh"))), collapse = "\n")
  expect_match(script, "'-entry' 'BUILD'", fixed = TRUE)
})

test_that("the in-app Build run pins NXF_VER and passes what the project config needs", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  writeLines(c(RiboPilot:::container_engine_block("docker"),
    "process { executor = 'local' }",
    "params {", "  rawDir = '/data/raw'", "}"), file.path(dir, ".config"))

  seen <- new.env()
  # a stand-in for the processx handle: nothing is spawned, and the run reads as
  # already finished so the streaming observer runs its exit branch too
  fake <- list(get_pid = function() 424242L, is_alive = function() FALSE,
    read_output_lines = function() c("executor >  local (1)", "run finished ok"),
    read_error_lines = function() character(0), kill = function() TRUE)
  testthat::local_mocked_bindings(
    nf_pin_version = function(...) "25.10.6",
    nf_spawn = function(cmd, args, wd, env) {
      seen$cmd <- cmd; seen$args <- args; seen$env <- env; fake
    }, .package = "RiboPilot")

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    session$setInputs(ref_nf_start = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
    session$flushReact()
  })

  expect_equal(seen$cmd, "nextflow")
  expect_equal(unname(seen$env[["NXF_VER"]]), "25.10.6")
  expect_equal(seen$args[which(seen$args == "-entry") + 1], "BUILD")
  expect_equal(seen$args[which(seen$args == "--projectDir") + 1], path.expand(dir))
  expect_false("--rawDir" %in% seen$args)

  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  log <- DBI::dbGetQuery(con, "SELECT log_file FROM reference WHERE ref_id = 'asteroidea'")$log_file[1]
  expect_true(!is.na(log) && file.exists(log))
  expect_true(any(grepl("run finished ok", readLines(log), fixed = TRUE)))
})

test_that("a hand-edited script is written as edited", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  withr::local_options(list(RiboPilot.headless = TRUE))

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    session$setInputs(ref_submit_script_text = "#!/bin/bash\n#SBATCH -p bigmem\nnextflow run mine")
    session$setInputs(ref_save_script = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })

  script <- paste(readLines(Sys.glob(file.path(dir, "jobs", "ribopilot_ref_*.sh"))), collapse = "\n")
  expect_true(grepl("#SBATCH -p bigmem", script, fixed = TRUE))
})

test_that("the deferred apply runs with a session but no reactive context", {
  # The app defers apply_clade through later(); withReactiveDomain() restores
  # the session but NOT a reactive context, so any reactiveVal read inside it
  # throws "operation not allowed without an active reactive context".
  # testServer() evaluates its expression inside isolate(), which hides that,
  # so the call below deliberately happens outside one.
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  handlers <- NULL; sess <- NULL
  shiny::testServer(function(input, output, session) {
    handlers <<- RiboPilot:::reference_handlers(input, output, session, deps)
    sess <<- session
  }, { session$flushReact() })

  expect_no_error(
    shiny::withReactiveDomain(sess, handlers$apply_clade(7588L, "asteroidea")))

  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT clade_name FROM reference")$clade_name,
    "Asteroidea")
})

test_that("the run modal prefills the key it is given, not the environment", {
  withr::local_envvar(c(NCBI_API_KEY = "from_env"))
  html <- as.character(RiboPilot:::reference_run_modal(1, "/lab/refs", api_key = "from_project"))
  expect_true(grepl("from_project", html, fixed = TRUE))
  expect_false(grepl("from_env", html, fixed = TRUE))
})

test_that("headless Build defers the submit so the wait popup paints first", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  deps$selected_ids <- function() "asteroidea"
  withr::local_options(list(RiboPilot.headless = TRUE))

  # a fake scheduler that records every call, so one click cannot submit twice
  bin <- withr::local_tempdir()
  calls <- file.path(bin, "calls.txt")
  writeLines(c("#!/bin/sh", sprintf('echo "$@" >> %s', shQuote(calls)),
    'echo "Your job 12345 has been submitted"'), file.path(bin, "qsub"))
  Sys.chmod(file.path(bin, "qsub"), "0755")
  withr::local_path(bin, action = "prefix")
  skip_if(is.null(RiboPilot:::detect_scheduler()), "no scheduler on PATH")

  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::reference_handlers(input, output, session, deps)
    do_build <- h$do_build
  }, {
    do_build("asteroidea")
    session$flushReact()
    session$setInputs(ref_submit_job = 1)
    session$flushReact()
    # the submit is queued behind the "Submitting..." popup, not run inline: a
    # blocking qsub here is what left the button looking dead and got it clicked
    # a second time
    expect_false(file.exists(calls))
    while (!later::loop_empty()) later::run_now()
  })

  expect_true(file.exists(calls))
  expect_length(readLines(calls), 1L)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  # the side effects that used to run inline still run after the deferral
  expect_equal(DBI::dbGetQuery(con,
    "SELECT reference_switch FROM reference WHERE ref_id = 'asteroidea'")$reference_switch, 4)
})

test_that("Suggest writes nothing until Apply", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  withr::local_options(list(ribopilot.count_runner = function(project, args, log = NULL) 0L))
  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_suggest_btn = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  # the fixture's row is stamped 1; a write would have re-stamped it
  expect_equal(DBI::dbGetQuery(con, "SELECT time_stamp FROM reference")$time_stamp, 1L)
})

test_that("the sample view names each sample's group and who set it", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), mk_split_db(dir)); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET ref_source = 'manual', user_group = 'A' WHERE ID = 'S1'")
  d <- RiboPilot:::reference_sample_table(con)
  expect_equal(d$ID, c("S1", "S2"))
  expect_equal(d$ref_group, c("Asteroidea", "Asteroidea"))
  expect_equal(d$set_by, c("you", "auto"))
  expect_equal(d$user_group, c("A", ""))
  expect_equal(nrow(RiboPilot:::reference_members(con)), 2L)
})
