# Writing a sample's taxon: lineage path, alternatives, and screen follow the taxon.

mk_tax_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, taxid INT, sci_name TEXT, rank TEXT,
    lineage TEXT, lineage_taxids TEXT, match_type TEXT, tax_alternatives TEXT, taxonomy_switch INT,
    taxonomy_lock INT, reference_support INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT,
    start INT, stop INT, strand TEXT)")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  # the re-screen also reads annotate and the review copies
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, n_features INT,
    itsx_warnings TEXT, itsx_origin TEXT, itsx_scores TEXT)")
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon, taxonomy_switch, taxonomy_lock)
    VALUES ('S1', 'Morus sp.', 1, 0)")
  con
}
# three plant accessions on S1's only unit
add_plant_hits <- function(con) {
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES ('S1', 1, 1, 'ACGT', 4)")
  for (k in 1:3) DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp,
      accession, lineage, lineage_taxids, bitscore, hit_class) VALUES ('S1', 1, 1, ?, 1, ?, ?, ?, ?, 'rdna')",
    params = list(k, paste0("ACC", k),
      "Viridiplantae;Streptophyta;Magnoliopsida;Rosales;Moraceae;Morus;Morus alba",
      "33090;35493;3398;3744;3487;3497;3498", 1000 - k))
}

test_that("resolving leaves a homonym unassigned and stores its alternatives", {
  con <- mk_tax_db()
  mt <- RiboPilot:::resolve_samples_taxonomy(con, "S1", local_taxdump())
  expect_equal(unname(mt), "ambiguous")
  row <- DBI::dbGetQuery(con, "SELECT * FROM samples WHERE ID = 'S1'")
  expect_true(is.na(row$taxid))
  expect_true(is.na(row$lineage_taxids))
  expect_equal(row$match_type, "ambiguous")
  expect_true(grepl("(Viridiplantae)", row$tax_alternatives, fixed = TRUE))
  expect_true(grepl("(Metazoa)", row$tax_alternatives, fixed = TRUE))
})

test_that("changing a taxon rewrites lineage_taxids in the same update", {
  con <- mk_tax_db()
  dir <- local_taxdump()
  RiboPilot:::assign_samples_taxid(con, "S1", 3497L, dir)
  RiboPilot:::assign_samples_taxid(con, "S1", 37577L, dir)
  row <- DBI::dbGetQuery(con, "SELECT * FROM samples WHERE ID = 'S1'")
  expect_equal(row$lineage_taxids, "33208;7711;8782;30410;30446;37577;")
  expect_true(is.na(row$tax_alternatives))
})

test_that("a taxon change re-screens the sample against its new lineage", {
  con <- mk_tax_db()
  dir <- local_taxdump()
  add_plant_hits(con)
  RiboPilot:::assign_samples_taxid(con, "S1", 3497L, dir)
  v <- DBI::dbGetQuery(con, "SELECT screen_verdict FROM screen WHERE ID = 'S1'")$screen_verdict
  expect_equal(v, "on_target")
  RiboPilot:::assign_samples_taxid(con, "S1", 37577L, dir)
  v <- DBI::dbGetQuery(con, "SELECT screen_verdict FROM screen WHERE ID = 'S1'")$screen_verdict
  expect_equal(v, "off_target")
})

test_that("re-resolving to a name that matches nothing clears the old lineage", {
  con <- mk_tax_db()
  dir <- local_taxdump()
  RiboPilot:::assign_samples_taxid(con, "S1", 37577L, dir)
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'Zzzzyx qqqqq' WHERE ID = 'S1'")
  mt <- RiboPilot:::resolve_samples_taxonomy(con, "S1", dir)
  expect_equal(unname(mt), "unresolved")
  row <- DBI::dbGetQuery(con, "SELECT * FROM samples WHERE ID = 'S1'")
  expect_true(is.na(row$lineage_taxids))
})

test_that("the resolve summary counts samples by confidence", {
  expect_equal(RiboPilot:::tax_count_text(c("exact", "exact", "fuzzy_species", "ambiguous", NA)),
    "Checked 5 samples: 2 exact, 1 fuzzy species, 1 ambiguous, 1 unresolved. Samples marked ambiguous or unresolved have no taxon yet: click a sample's ID to pick one, then lock.")
  expect_equal(RiboPilot:::tax_count_text(c("exact", "exact")),
    "Checked 2 samples: 2 exact. Review and lock.")
})

test_that("Assign in the taxon search leaves a taxonomy-locked sample unchanged", {
  d <- withr::local_tempdir()
  db <- file.path(d, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbWriteTable(con, "samples", data.frame(ID = c("S1", "S2"), taxonomy_lock = c(1L, 0L)))
  DBI::dbDisconnect(con)
  seen <- new.env()
  testthat::local_mocked_bindings(
    mp_toast = function(text, ...) seen$toast <- text,
    assign_samples_taxid = function(con, ids, tid, td) { seen$ids <- ids; length(ids) },
    .package = "RiboPilot")
  deps <- list(con = function() DBI::dbConnect(RSQLite::SQLite(), db), refresh = shiny::reactiveVal(0),
    taxdump = function() d, samples_for = function(ids) data.frame(ID = ids, Taxon = "", sci_name = NA, tax_alternatives = NA))
  shiny::testServer(function(input, output, session) {
    session$userData$open <- RiboPilot:::taxonomy_handlers(input, output, session, deps)
  }, {
    session$userData$open("S1")
    session$setInputs(tax_pick = "7955", tax_apply = 1)
    expect_match(seen$toast, "Left 1 locked sample unchanged")
    expect_null(seen$ids)
    session$userData$open(c("S1", "S2"))
    session$setInputs(tax_apply = 2)
    for (i in 1:5) later::run_now(0.2)
    expect_equal(seen$ids, "S2")
  })
})

PLANT_LIN <- "Viridiplantae;Streptophyta;Magnoliopsida;Rosales;Moraceae;Morus;Morus alba"
PLANT_T <- "33090;35493;3398;3744;3487;3497;3498"
BIRD_LIN <- "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;unclassified Morus species"
BIRD_T <- "33208;7711;8782;30410;30446;37577;"
ITSX_SC <- "M: 4 2.9e-16 49.9, T: 4 1e-08 33.15"

# S1 has just moved from the plant Morus to the bird Morus. Its four contigs
# and S2's one hold plant hits and were checked under the plant lineage;
# S1.1.4 was never annotated. Copies: S1.1.1 done and locked, holding bird
# hits; S1.1.2 not run; S1.1.3 running; S1.1.4 on hold; S2.1.1 not run.
rescreen_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE samples SET lineage = ?, lineage_taxids = ? WHERE ID = 'S1'",
    params = list(BIRD_LIN, BIRD_T))
  units <- data.frame(ID = c("S1", "S1", "S1", "S1", "S2"), scaffold = c(1:4, 1L))
  for (i in seq_len(nrow(units))) {
    p <- list(units$ID[i], units$scaffold[i])
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence)
      VALUES (?, 1, ?, 'linear', 100, ?)", params = c(p, strrep("ACGT", 25)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        annotate_opts, curate_opts, n_features, n_flags, curate_flags, itsx_warnings, itsx_origin,
        itsx_scores)
      VALUES (?, 1, ?, 2, 0, 'default', 'default', 1, 99, 'stale', '', 'T', ?)", params = c(p, ITSX_SC))
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length, strand)
      VALUES (?, 1, ?, '18S', 1, 90, 90, '+')", params = p)
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, lineage, lineage_taxids,
      is_nuclear_rrna) VALUES (?, 1, ?, ?, ?, 1)", params = c(p, PLANT_LIN, PLANT_T))
    for (k in 1:3) DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp,
        lineage, lineage_taxids, bitscore, hit_class) VALUES (?, 1, ?, ?, 1, ?, ?, ?, 'rdna')",
      params = c(p, k, PLANT_LIN, PLANT_T, 1000 - k))
  }
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M',
    itsx_warnings = 'itsx_chimeric;itsx_origin(score: Metazoa 49.9 vs Tracheophyta 33.15)'
    WHERE ID = 'S1' AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'not_run(ITSx did not finish)'
    WHERE ID = 'S1' AND scaffold = 3")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL, itsx_origin = NULL,
    itsx_warnings = 'itsx_origin(Tracheophyta)' WHERE ID = 'S1' AND scaffold = 4")

  for (i in seq_len(nrow(units))) RiboPilot:::review_copy(con, units$ID[i], 1L, units$scaffold[i])
  DBI::dbExecute(con, "UPDATE review_units SET rev = 1000, time_stamp = 1, n_flags = 99,
    curate_flags = 'stale', itsx_warnings = '', itsx_origin = 'T', itsx_scores = ?, check_switch = 1",
    params = list(ITSX_SC))
  for (s in 1:4) {
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = ? WHERE ID = 'S1' AND scaffold = ?",
      params = list(c(2L, 1L, 4L, 0L)[s], s))
  }
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE ID = 'S1' AND scaffold = 1")
  DBI::dbExecute(con, "DELETE FROM review_blast_hits WHERE ID = 'S1' AND scaffold = 1")
  for (k in 1:3) DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank,
      hsp, lineage, lineage_taxids, bitscore, hit_class) VALUES ('S1', 1, 1, ?, 1, ?, ?, ?, 'rdna')",
    params = list(k, "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;Morus bassanus",
      paste0(BIRD_T, "37578"), 1000 - k))
}

test_that("a taxon change re-curates the contigs, refreshes done copies, and requeues running ones", {
  con <- chimera_test_con(c("S1", "S2"))
  on.exit(DBI::dbDisconnect(con))
  rescreen_seed(con)
  RiboPilot:::rescreen_samples(con, "S1")
  q <- function(sql) DBI::dbGetQuery(con, sql)

  # originals: screened against the new lineage, curated when annotated
  expect_equal(q("SELECT screen_verdict FROM screen ORDER BY ID, scaffold")$screen_verdict,
    rep("off_target", 4))
  an <- q("SELECT ID, scaffold, curate_flags, itsx_warnings FROM annotate ORDER BY ID, scaffold")
  expect_equal(an$curate_flags[4:5], c("stale", "stale"))
  for (i in 1:3) {
    expect_true("contaminant" %in% strsplit(an$curate_flags[i], "; ")[[1]], info = an$curate_flags[i])
  }
  expect_equal(an$itsx_warnings, c("itsx_chimeric",
    "itsx_origin(score: Tracheophyta 33.15 vs Metazoa 49.9)", "not_run(ITSx did not finish)",
    "itsx_origin(Tracheophyta)", ""))

  # copies: the done one re-scored from its own hits; queued and running ones requeued
  ru <- q("SELECT ID, scaffold, check_switch, CAST(rev AS TEXT) AS rev, time_stamp, curate_flags,
    itsx_warnings FROM review_units ORDER BY ID, scaffold")
  expect_equal(ru$check_switch, c(2L, 1L, 4L, 0L, 1L))
  expect_equal(ru$rev, c("1000", "1001", "1001", "1000", "1000"))
  expect_equal(ru$time_stamp > 1, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_false(identical(ru$curate_flags[1], "stale"))
  expect_equal(ru$curate_flags[2:5], rep("stale", 4))
  expect_equal(ru$itsx_warnings, c("itsx_origin(score: Tracheophyta 33.15 vs Metazoa 49.9)",
    "", "", "", ""))
  expect_equal(q("SELECT ID, scaffold, screen_verdict FROM review_screen"),
    data.frame(ID = "S1", scaffold = 1L, screen_verdict = "on_target"))
})

test_that("one failing copy refresh does not abort a taxon change", {
  con <- chimera_test_con(c("S1", "S2"))
  on.exit(DBI::dbDisconnect(con))
  rescreen_seed(con)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2 WHERE ID = 'S1' AND scaffold = 2")
  real <- RiboPilot:::screen_one
  testthat::local_mocked_bindings(screen_one = function(con, id, path = 1L, scaffold = 1L,
                                                        force = FALSE, review = FALSE) {
    if (review && scaffold == 1L) stop("bad copy")
    real(con, id, path, scaffold, force = force, review = review)
  }, .package = "RiboPilot")
  expect_no_error(RiboPilot:::rescreen_samples(con, "S1"))
  ru <- DBI::dbGetQuery(con, "SELECT scaffold, CAST(rev AS TEXT) AS rev, curate_flags
    FROM review_units WHERE ID = 'S1' ORDER BY scaffold")
  expect_false(identical(ru$curate_flags[2], "stale"))
  expect_equal(ru$rev[3], "1001")
})
