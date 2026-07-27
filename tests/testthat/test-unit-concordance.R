# Sibling concordance: is a scaffold the odd one out among its own siblings?
# Fires where lineages_conflict() structurally cannot (metazoan-in-metazoan),
# and must go SILENT on sparse clades where every scaffold's best hit differs.

# kingdom;phylum;class;order;family;genus;species
tx <- function(k, p, c, o, f, g, s) paste(k, p, c, o, f, g, s, sep = ";")
nm <- function(k, p, c, o, f, g, s) paste(k, p, c, o, f, g, s, sep = ";")

# Notropis species, same genus, different species
FISH <- function(sp) tx(33208, 7711, 186623, 7952, 2743709, 7954, sp)
FISHN <- function(sp) nm("Metazoa", "Chordata", "Actinopteri", "Cypriniformes",
  "Leuciscidae", "Notropis", paste("Notropis sp", sp))
# a copepod: same kingdom, different phylum
COPE  <- tx(33208, 6656, 6854, 6843, 1, 2, 3)
COPEN <- nm("Metazoa", "Arthropoda", "Copepoda", "Cyclopoida", "Cyclopidae", "Cyclops", "Cyclops sp")

mk <- function(taxids, lineages, ignore = 0) {
  n <- length(taxids)
  data.frame(path = rep(1L, n), scaffold = seq_len(n),
    screen_taxids = taxids, screen_lineage = lineages,
    ignore = rep(ignore, length.out = n), stringsAsFactors = FALSE)
}

test_that("a metazoan contaminant among metazoan siblings is flagged", {
  d <- mk(c(FISH(101), FISH(102), FISH(103), FISH(104), COPE),
          c(FISHN(1), FISHN(2), FISHN(3), FISHN(4), COPEN))
  r <- RiboPilot:::unit_concordance(d)
  expect_equal(r$odd_one_out, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  # species differ among the fish, so the majority resolves at genus
  expect_equal(unique(r$sib_rank), "genus")
  expect_equal(unique(r$sib_taxon), "Notropis")
})

test_that("a sparse clade with no two hits alike raises nothing", {
  # Five scaffolds of one octocoral, each best-hitting a different genus in the
  # same class. A genus-level concordance test would flag all five.
  d <- mk(
    c(tx(33208, 6073, 6101, 11, 21, 31, 41), tx(33208, 6073, 6101, 12, 22, 32, 42),
      tx(33208, 6073, 6101, 13, 23, 33, 43), tx(33208, 6073, 6101, 14, 24, 34, 44),
      tx(33208, 6073, 6101, 15, 25, 35, 45)),
    rep(nm("Metazoa", "Cnidaria", "Anthozoa", "o", "f", "g", "s"), 5))
  r <- RiboPilot:::unit_concordance(d)
  expect_false(any(r$odd_one_out))
  expect_equal(unique(r$sib_rank), "class")   # they agree, just not deeply
})

test_that("fewer than three resolved siblings never votes", {
  d <- mk(c(FISH(101), COPE), c(FISHN(1), COPEN))
  r <- RiboPilot:::unit_concordance(d)
  expect_false(any(r$odd_one_out))
  expect_true(all(is.na(r$sib_rank)))

  # two resolved + one unresolved is still only two votes
  d2 <- mk(c(FISH(101), COPE, NA), c(FISHN(1), COPEN, NA))
  expect_false(any(RiboPilot:::unit_concordance(d2)$odd_one_out))
})

test_that("a tie is not a majority", {
  d <- mk(c(FISH(101), FISH(101), COPE, COPE), c(FISHN(1), FISHN(1), COPEN, COPEN))
  r <- RiboPilot:::unit_concordance(d)
  # 2-2 at every rank below kingdom: no strict majority, so no accusation
  expect_false(any(r$odd_one_out))
  expect_equal(unique(r$sib_rank), "kingdom")  # all four are Metazoa
})

test_that("units with no resolved hit are never accused", {
  d <- mk(c(FISH(101), FISH(102), FISH(103), NA, ""),
          c(FISHN(1), FISHN(2), FISHN(3), NA, ""))
  r <- RiboPilot:::unit_concordance(d)
  expect_false(any(r$odd_one_out))          # the NA rows abstain, they do not fail
  expect_equal(nrow(r), 5)
})

test_that("concordance returns one row per input unit", {
  expect_equal(nrow(RiboPilot:::unit_concordance(mk(character(0), character(0)))), 0)
  d <- mk(c(FISH(1), FISH(2), FISH(3)), c(FISHN(1), FISHN(2), FISHN(3)))
  expect_equal(nrow(RiboPilot:::unit_concordance(d)), 3)
})
