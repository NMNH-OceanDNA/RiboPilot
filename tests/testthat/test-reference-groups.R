# Clade grouping: walk each sample's lineage up until the clade has enough
# GenBank records, then merge samples that land on the same clade.

# Rank-aligned kingdom..species lineages, as samples.lineage / lineage_taxids
# hold them (7 fields, ";"-joined, empty taxid at a filled rank).
FISH <- list(
  names = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
  taxids = "33208;7711;186623;7952;1002673;7954;7955")
FISH2 <- list(
  names = "Metazoa;Chordata;Actinopteri;Perciformes;Percidae;Perca;Perca flavescens",
  taxids = "33208;7711;186623;8111;8168;8166;8167")
STAR <- list(
  names = "Metazoa;Echinodermata;Asteroidea;Forcipulatida;Asteriidae;Asterias;Asterias rubens",
  taxids = "33208;7586;7588;41692;7607;7608;7610")
STAR2 <- list(
  names = "Metazoa;Echinodermata;Asteroidea;Valvatida;Asterinidae;Patiria;Patiria miniata",
  taxids = "33208;7586;7588;41695;41693;104852;46514")

mk_ladder_db <- function(rows) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT, ref_id TEXT)")
  for (nm in names(rows)) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,?,?,?,1,2,NULL)",
      params = list(nm, NA_integer_, rows[[nm]]$names, rows[[nm]]$taxids))
  }
  con
}

# counts by taxid; anything absent counts as 0 records
counter <- function(map) function(taxid, name) {
  v <- map[[as.character(taxid)]]
  if (is.null(v)) 0L else as.integer(v)
}

test_that("sample_ladders yields one row per resolved rank, deepest last", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_equal(nrow(l), 7L)
  expect_equal(l$rank[1], "kingdom")
  expect_equal(l$rank[7], "species")
  expect_equal(l$name[3], "Actinopteri")
  expect_equal(l$taxid[3], 186623L)
})

test_that("a rank with no taxid is dropped from the ladder", {
  gappy <- list(names = "Metazoa;Chordata;unclassified Chordata;Cypriniformes;Danionidae;Danio;Danio rerio",
    taxids = "33208;7711;;7952;1002673;7954;7955")
  con <- mk_ladder_db(list(F1 = gappy)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_false("class" %in% l$rank)
  expect_equal(nrow(l), 6L)
})

test_that("fish and sea stars land in two groups", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2, S1 = STAR, S2 = STAR2))
  on.exit(DBI::dbDisconnect(con))
  # families and orders are record-poor here; the two classes clear the floor
  cf <- counter(list("186623" = 30000, "7588" = 4102, "7711" = 900000, "7586" = 60000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 50000)
  expect_equal(nrow(res$groups), 2L)
  expect_setequal(res$groups$clade_name, c("Actinopteri", "Asteroidea"))
  expect_equal(res$members$ref_id[res$members$ID == "F1"], "actinopteri")
  expect_equal(res$members$ref_id[res$members$ID == "S2"], "asteroidea")
  expect_equal(res$groups$n_samples[res$groups$clade_name == "Actinopteri"], 2L)
})

test_that("a record-rich family stops the walk at family", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("1002673" = 5000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(res$groups$clade_name, "Danionidae")
  expect_equal(res$groups$clade_rank, "family")
})

test_that("the ceiling stops the walk before an over-broad clade", {
  con <- mk_ladder_db(list(F1 = FISH)); on.exit(DBI::dbDisconnect(con))
  # nothing reaches the floor until kingdom, which is far over the ceiling
  cf <- counter(list("33208" = 900000, "7711" = 400000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  # the walk keeps the best clade it saw rather than one above the ceiling
  expect_equal(res$groups$clade_rank, "phylum")
  expect_equal(res$groups$clade_name, "Chordata")
})

test_that("a nested group merges into its parent when the parent is under the ceiling", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2)); on.exit(DBI::dbDisconnect(con))
  # Danionidae alone clears the floor, Perciformes only does at class level
  cf <- counter(list("1002673" = 2000, "186623" = 20000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(nrow(res$groups), 1L)
  expect_equal(res$groups$clade_name, "Actinopteri")
  expect_equal(res$groups$n_samples, 2L)
})

test_that("a nested group stays separate when the parent is over the ceiling", {
  con <- mk_ladder_db(list(F1 = FISH, F2 = FISH2)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("1002673" = 2000, "186623" = 90000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(nrow(res$groups), 2L)
})

test_that("group_ladder reports coverage for the shared lineage", {
  con <- mk_ladder_db(list(S1 = STAR, S2 = STAR2)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  g <- RiboPilot:::group_ladder(l, c("S1", "S2"))
  expect_equal(g$name[nrow(g)], "Asteroidea")
  expect_equal(g$n_samples[nrow(g)], 2L)
  expect_true(all(is.na(g$n_records)))
})

test_that("ladder_coverage names the samples under a taxid", {
  con <- mk_ladder_db(list(F1 = FISH, S1 = STAR, S2 = STAR2)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  expect_setequal(RiboPilot:::ladder_coverage(l, 7588L), c("S1", "S2"))
  expect_setequal(RiboPilot:::ladder_coverage(l, 41692L), "S1")
})

test_that("samples with no usable lineage are reported, not grouped", {
  con <- mk_ladder_db(list(F1 = FISH))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO samples VALUES('X1',NULL,NULL,NULL,0,1,NULL)")
  l <- RiboPilot:::sample_ladders(con)
  expect_false("X1" %in% l$ID)
  res <- RiboPilot:::suggest_reference_groups(l, counter(list("1002673" = 5000)))
  expect_false("X1" %in% res$members$ID)
})

test_that("ref_id_from_name makes a filename-safe slug", {
  expect_equal(RiboPilot:::ref_id_from_name("Asteroidea"), "asteroidea")
  expect_equal(RiboPilot:::ref_id_from_name("Danio rerio"), "danio_rerio")
})
