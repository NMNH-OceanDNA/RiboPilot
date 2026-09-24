# Rank-aware taxonomy primitives. All pure: no taxonkit, no taxdump, no network.
# Lineage strings below are real `taxonkit reformat -F -t` output shapes.

# Danio rerio
FISH_N <- "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio"
FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
# Gallus gallus: same kingdom + phylum, diverges at class
BIRD_N <- "Metazoa;Chordata;Aves;Galliformes;Phasianidae;Gallus;Gallus gallus"
BIRD_T <- "33208;7711;8782;8976;9005;9030;9031"
# Escherichia coli: shares nothing
BACT_N <- "Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"
BACT_T <- "3379134;1224;1236;91347;543;561;562"
# A congener of the fish
CONGENER_T <- "33208;7711;186623;7952;2743709;7954;7998"
# Keratoisis: genus-level, species rank FILLED by -F so its taxid is empty
CORAL_N <- "Metazoa;Cnidaria;Anthozoa;Scleralcyonacea;Keratoisididae;Keratoisis;unclassified Keratoisis species"
CORAL_T <- "33208;6073;6101;3028851;2931962;243285;"

test_that("placeholder and uninformative taxa are recognised", {
  f <- RiboPilot:::is_placeholder_taxon
  expect_true(f("unclassified Keratoisis species"))
  expect_true(f("uncultured eukaryote"))
  expect_true(f("environmental samples"))
  expect_true(f("unidentified marine sponge"))
  expect_true(f(NA_character_))
  expect_true(f(""))
  expect_false(f("Danio rerio"))
  expect_false(f("Keratoisis"))
})

test_that("rank-aligned parsing keeps seven slots and blanks filled ranks", {
  tx <- RiboPilot:::parse_rank_taxids(CORAL_T)
  expect_length(tx, 7)
  expect_equal(unname(tx[["genus"]]), 243285L)
  expect_true(is.na(tx[["species"]]))          # -F filled it, taxid is empty

  nm <- RiboPilot:::parse_rank_names(CORAL_N)
  expect_equal(unname(nm[["genus"]]), "Keratoisis")
  expect_true(is.na(nm[["species"]]))          # "unclassified ..." is not a name
  expect_equal(unname(nm[["kingdom"]]), "Metazoa")
})

test_that("rank_distance counts rank steps, not string prefixes", {
  expect_equal(rank_distance(FISH_T, FISH_T), 0)        # same species
  expect_equal(rank_distance(CONGENER_T, FISH_T), 1)    # same genus
  expect_equal(rank_distance(BIRD_T, FISH_T), 5)        # shares phylum Chordata
  expect_equal(rank_distance(BACT_T, FISH_T), 7)        # shares nothing
  expect_equal(rank_distance(CORAL_T, FISH_T), 6)       # shares kingdom Metazoa only
})

test_that("a metazoan contaminant is distinguishable from the host", {
  # The case lineages_conflict() structurally could not fire on: both are
  # animals, so a shared-token test reports "clean".
  expect_equal(rank_distance(CORAL_T, FISH_T), 6)
  expect_false(rank_at_least(rank_distance(CORAL_T, FISH_T), "phylum"))
  expect_true(rank_at_least(rank_distance(BIRD_T, FISH_T), "phylum"))
})

test_that("an unassigned intermediate rank does not break deeper matching", {
  # Same lineage as the coral but with class and order unresolved. The old
  # positional-prefix comparison stopped at the first gap and reported kingdom;
  # the taxid-aligned version still sees the shared family.
  gappy <- "33208;6073;;;2931962;243285;"
  expect_equal(rank_distance(gappy, CORAL_T), 1)   # agrees down to genus
  expect_equal(RiboPilot:::deepest_shared_rank(
    RiboPilot:::parse_rank_taxids(gappy),
    RiboPilot:::parse_rank_taxids(CORAL_T)), 6)
})

test_that("rank labels and thresholds read the way a biologist expects", {
  expect_equal(rank_distance_label(0), "species")
  expect_equal(rank_distance_label(1), "genus")
  expect_equal(rank_distance_label(6), "kingdom")
  expect_equal(rank_distance_label(7), "none")
  expect_true(is.na(rank_distance_label(NA)))

  expect_true(rank_at_least(2, "family"))
  expect_false(rank_at_least(3, "family"))
  expect_true(rank_at_least(6, "kingdom"))
  expect_false(rank_at_least(7, "kingdom"))
  expect_error(rank_at_least(1, "tribe"))
})

test_that("bitscore-weighted LCA survives one mislabelled accession", {
  # Three real fish hits plus one bogus bacterial record with a lower bitscore.
  paths <- list(FISH_T, CONGENER_T, FISH_T, BACT_T)
  bs <- c(2000, 1950, 1900, 400)
  lca <- bitscore_weighted_lca(paths, bs, names_paths = list(FISH_N, FISH_N, FISH_N, BACT_N))
  expect_equal(lca$rank, "genus")          # bacterium excluded by the bitscore window
  expect_equal(lca$taxid, 7954L)
  expect_equal(lca$n_hits, 3L)
})

test_that("LCA backs off to a shallower rank when top hits genuinely disagree", {
  lca <- bitscore_weighted_lca(list(FISH_T, BIRD_T), c(1000, 990))
  expect_equal(lca$rank, "phylum")
  expect_equal(lca$n_hits, 2L)

  none <- bitscore_weighted_lca(list(FISH_T, BACT_T), c(1000, 1000))
  expect_equal(none$rank, "none")
  expect_true(is.na(none$taxid))
})

test_that("LCA ignores ranks that are not resolved in every selected hit", {
  # coral has no species taxid; the LCA must not claim a species-level call
  lca <- bitscore_weighted_lca(list(CORAL_T, CORAL_T), c(500, 500))
  expect_equal(lca$rank, "genus")
})

test_that("degenerate LCA inputs do not error", {
  e <- bitscore_weighted_lca(list(FISH_T), NA_real_)
  expect_equal(e$n_hits, 0L)
  expect_equal(e$rank, "none")
})
