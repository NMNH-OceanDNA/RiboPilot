# The reference tab's table shape and modal content.

mk_tab_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, NULL, 1000, 25000)")
  DBI::dbExecute(con, "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank,
    query, n_samples, n_records, reference_switch, source, seeds_db, labels_db, time_stamp)
    VALUES ('asteroidea', 7588, 'Asteroidea', 'class', '\"Asteroidea\"[Organism]', 4, 4102,
      2, 'built', '/refs/a_seeds.fa', '/refs/a_labels.fa', 1)")
  con
}

test_that("reference_table exposes ref_id as ID and the seed path as database", {
  con <- mk_tab_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::reference_table(con)
  expect_equal(d$ID[1], "asteroidea")
  expect_equal(d$clade_name[1], "Asteroidea")
  expect_equal(d$database[1], "/refs/a_seeds.fa")
  expect_true("reference_switch" %in% names(d))
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

test_that("the build modal names the folder it will write to", {
  m <- RiboPilot:::reference_build_modal(c("asteroidea", "actinopteri"), "/lab/refs")
  html <- as.character(m)
  expect_true(grepl("/lab/refs", html, fixed = TRUE))
  expect_true(grepl("2", html, fixed = TRUE))
  expect_true(grepl("ref_api_key", html, fixed = TRUE))
  expect_true(grepl("ref_conda", html, fixed = TRUE))
  expect_true(grepl("ref_refseq", html, fixed = TRUE))
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
# taxid's record count is pre-seeded into the reference folder's cache so
# cached_count_fn() never falls through to genbank_rrna_count(), which would
# spawn a process.
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
  DBI::dbExecute(con, "INSERT INTO reference_config VALUES(1, ?, 1000, 25000)", params = list(ref_dir))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, lineage TEXT, lineage_taxids TEXT, ref_id TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',
    'Metazoa;Echinodermata;Asteroidea;Forcipulatida;Asteriidae;Asterias;Asterias rubens',
    '33208;7586;7588;41692;7607;7608;7610', 'asteroidea')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2',
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
  withr::local_options(list(ribopilot.count_runner = function(project, args) {
    calls <<- calls + 1L; 0L
  }))
  shiny::testServer(function(input, output, session) {
    RiboPilot:::reference_handlers(input, output, session, deps)
  }, {
    session$setInputs(ref_suggest_btn = 1)
    session$flushReact()
    while (!later::loop_empty()) later::run_now()
  })
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile); on.exit(DBI::dbDisconnect(con))
  expect_gt(nrow(DBI::dbGetQuery(con, "SELECT * FROM reference")), 0L)
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
  withr::local_options(list(ribopilot.count_runner = function(project, args) {
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
  expect_true(all(grepl("[Organism]", sent$query, fixed = TRUE)))
})

test_that("apply-clade counts a cached clade without calling the runner", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dbfile <- mk_split_db(dir)
  deps <- mk_split_deps(dbfile, dir)
  calls <- 0L
  withr::local_options(list(ribopilot.count_runner = function(project, args) {
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
