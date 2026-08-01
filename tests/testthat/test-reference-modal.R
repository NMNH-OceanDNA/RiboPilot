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
  expect_true(grepl("4102", html, fixed = TRUE))
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
  expect_true(grepl("5000", html, fixed = TRUE))
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
})

test_that("the existing-database modal still offers typed paths with an empty manifest", {
  man <- RiboPilot:::read_ref_manifest(withr::local_tempdir())
  m <- RiboPilot:::reference_existing_modal(man, "asteroidea")
  html <- as.character(m)
  expect_true(grepl("ref_existing_seeds", html, fixed = TRUE))
  expect_true(grepl("ref_existing_labels", html, fixed = TRUE))
})
