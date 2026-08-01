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
  # nothing ever reaches the floor (100000); phylum and kingdom both clear it
  # by count but are over the ceiling (25000), so neither can be the fallback
  cf <- counter(list("186623" = 0, "7952" = 0, "1002673" = 0,
    "7711" = 50000, "33208" = 60000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 100000, max_records = 25000)
  # the walk keeps the deepest clade under the ceiling rather than one above it
  expect_equal(res$groups$clade_rank, "class")
  expect_equal(res$groups$clade_name, "Actinopteri")
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

test_that("group_ladder returns the empty frame when members share no rung", {
  # kingdom is unresolved (dropped) for both, and they diverge at phylum, so
  # nothing survives the unanimity filter
  gappy_fish <- list(
    names = "unclassified Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
    taxids = ";7711;186623;7952;1002673;7954;7955")
  gappy_star <- list(
    names = "unclassified Metazoa;Echinodermata;Asteroidea;Forcipulatida;Asteriidae;Asterias;Asterias rubens",
    taxids = ";7586;7588;41692;7607;7608;7610")
  con <- mk_ladder_db(list(X = gappy_fish, Y = gappy_star)); on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  g <- RiboPilot:::group_ladder(l, c("X", "Y"))
  expect_equal(nrow(g), 0L)
  expect_named(g, c("level", "rank", "taxid", "name", "n_samples", "n_records"))
})

test_that("a chain of merges folds to a fixed point rather than stranding a sample", {
  # A's own family clears the floor; B's family fails but its order (shared
  # with A) clears it; C's family and order fail but its class (shared with
  # both) clears it. A naive single pass leaves A merged into B's order while
  # B merges past it into C's class, stranding A in a group nobody chose.
  fish_b <- list(
    names = "Metazoa;Chordata;Actinopteri;Cypriniformes;FakeCyprinidFamily;FakeGenus;Fake species",
    taxids = "33208;7711;186623;7952;9001;9002;9003")
  con <- mk_ladder_db(list(A = FISH, B = fish_b, C = FISH2)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("1002673" = 2000, "9001" = 0, "7952" = 5000,
    "8168" = 0, "8111" = 0, "186623" = 40000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 50000)
  expect_equal(nrow(res$groups), 1L)
  expect_equal(res$groups$clade_name, "Actinopteri")
  expect_equal(res$groups$n_samples, 3L)
  expect_setequal(res$members$ref_id, "actinopteri")
})

test_that("groups are deduped by clade taxid, and colliding slugs are disambiguated", {
  collide1 <- list(
    names = "Metazoa;Chordata;Actinopteri;Cypriniformes;Foo bar;FooGenus;Foo species",
    taxids = "33208;7711;186623;7952;111111;111112;111113")
  collide2 <- list(
    names = "Metazoa;Chordata;Actinopteri;Cypriniformes;Foo-bar;BarGenus;Bar species",
    taxids = "33208;7711;186623;7952;222222;222223;222224")
  con <- mk_ladder_db(list(A = collide1, B = collide2)); on.exit(DBI::dbDisconnect(con))
  cf <- counter(list("111111" = 2000, "222222" = 2000))
  res <- RiboPilot:::suggest_reference_groups(RiboPilot:::sample_ladders(con), cf,
    min_records = 1000, max_records = 25000)
  expect_equal(nrow(res$groups), 2L)
  expect_equal(length(unique(res$groups$ref_id)), 2L)
  expect_setequal(res$groups$clade_taxid, c(111111L, 222222L))
  a_ref <- res$members$ref_id[res$members$ID == "A"]
  b_ref <- res$members$ref_id[res$members$ID == "B"]
  expect_false(identical(a_ref, b_ref))
  expect_equal(res$groups$clade_taxid[res$groups$ref_id == a_ref], 111111L)
  expect_equal(res$groups$clade_taxid[res$groups$ref_id == b_ref], 222222L)
})

test_that("candidate_rungs covers every rung at or above the start rank, once each", {
  con <- mk_ladder_db(list(F1 = FISH, S1 = STAR, S2 = STAR2))
  on.exit(DBI::dbDisconnect(con))
  l <- RiboPilot:::sample_ladders(con)
  r <- RiboPilot:::candidate_rungs(l, start_rank = "family")
  expect_true(all(c("ref_id", "taxid", "clade_name", "query") %in% names(r)))
  expect_false(any(duplicated(r$taxid)))
  expect_true(all(r$taxid %in% l$taxid[l$rank %in%
    c("kingdom", "phylum", "class", "order", "family")]))
  expect_false(any(r$taxid %in% l$taxid[l$rank %in% c("genus", "species")]))
  expect_true(grepl("[Organism]", r$query[1], fixed = TRUE))
  expect_equal(r$ref_id[r$taxid == 7588L], "asteroidea")
})

test_that("candidate_rungs on an empty ladder returns the empty shape", {
  con <- mk_ladder_db(list()); on.exit(DBI::dbDisconnect(con))
  r <- RiboPilot:::candidate_rungs(RiboPilot:::sample_ladders(con))
  expect_equal(nrow(r), 0L)
  expect_true(all(c("ref_id", "taxid", "clade_name", "query") %in% names(r)))
})
