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
