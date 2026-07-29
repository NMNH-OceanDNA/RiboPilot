# The screen's taxonomy channel. It answers "do this unit's hits agree on a
# taxon, and is that taxon the host's". It needs NO identity threshold, which is
# why it survives the envelope being deferred (see doc 08 section 3.2.1).

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
FUNGI_T <- "4751;4890;147545;5042;1131492;5052;5062"
COPE_T <- "33208;6656;6854;6843;1;2;3"

# one hit row as blast_hits stores it
hit <- function(tax, bs, cls = "rdna", lin = "K;p;c;o;f;g;s") {
  data.frame(hit_class = cls, bitscore = bs, lineage_taxids = tax,
    lineage = lin, stringsAsFactors = FALSE)
}
hits <- function(...) do.call(rbind, list(...))

test_that("the constants are exactly the documented strings", {
  expect_equal(RiboPilot:::SCREEN_VERDICTS,
    c("on_target", "off_target", "no_reference", "uncertain"))
  expect_equal(RiboPilot:::SCREEN_CLASSES,
    c("bacterial", "archaeal", "fungal", "algal_plant", "other_euk",
      "metazoan_other", "human"))
  expect_equal(RiboPilot:::SCREEN_NO_REF_CODES,
    c("no_hits", "no_taxdump", "unresolved_sample_lineage", "unresolved_hit_lineages"))
  expect_equal(sum(unlist(RiboPilot:::SCREEN_WEIGHTS)), 1)
})

test_that("rank depth runs from species down to nothing", {
  expect_equal(RiboPilot:::screen_rank_depth("species"), 1)
  expect_equal(RiboPilot:::screen_rank_depth("kingdom"), 1 / 7, tolerance = 1e-6)
  expect_equal(RiboPilot:::screen_rank_depth("none"), 0)
  expect_equal(RiboPilot:::screen_rank_depth(NA_character_), 0)
  # monotone: a deeper rank is never worth less
  d <- vapply(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    RiboPilot:::screen_rank_depth, numeric(1))
  expect_false(is.unsorted(d))
})

test_that("each confidence term is computed independently and weighted", {
  tm <- RiboPilot:::screen_confidence_terms(lca_n_hits = 10, lca_rank = "species",
    margin = 500, top_bitscore = 500)
  expect_equal(tm$support, 1)
  expect_equal(tm$resolution, 1)
  expect_equal(tm$separation, 1)
  expect_equal(tm$confidence, 1)

  # support saturates at 10 hits and scales below it
  expect_equal(RiboPilot:::screen_confidence_terms(5, "species", 500, 500)$support, 0.5)
  expect_equal(RiboPilot:::screen_confidence_terms(50, "species", 500, 500)$support, 1)
  # a NA margin contributes nothing rather than erroring
  tm2 <- RiboPilot:::screen_confidence_terms(10, "species", NA_real_, 500)
  expect_equal(tm2$separation, 0)
  expect_equal(tm2$confidence, 0.6)
})

test_that("a concordant unit is on_target", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_equal(r$host_concordant, 1L)
  expect_true(is.na(r$class))
})

test_that("a discordant unit is off_target and names its subclass", {
  h <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "off_target")
  expect_equal(r$host_concordant, 0L)
  expect_equal(r$class, "fungal")
})

test_that("placeholder taxa are excluded from the vote entirely", {
  # uncultured records dominate core_nt 18S and resolve to superkingdom only,
  # which every rank rule would read as maximum taxonomic distance
  h <- hits(hit(FISH_T, 900, lin = "Eukaryota;;;;;;uncultured eukaryote"),
            hit(FISH_T, 890, lin = "Eukaryota;;;;;;environmental sample"),
            hit(FISH_T, 880, lin = "Eukaryota;;;;;;uncultured organism"))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_hit_lineages")
})

test_that("no hits at all is no_reference, not off_target", {
  r <- RiboPilot:::screen_channel_taxonomy(hits()[0, ], FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "no_hits")
})

test_that("an unresolvable sample lineage is no_reference with its own code", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, NA_character_)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_sample_lineage")
})

test_that("unanimity is maximal separation, not missing separation", {
  # a unit whose hits ALL agree with the host has no opposing hit to measure a
  # margin against. Scoring that 0 would put the cleanest possible scaffold
  # below the uncertainty gate, which is backwards.
  tm <- RiboPilot:::screen_confidence_terms(3, "species", NA_real_, 900,
    unanimous = TRUE)
  expect_equal(tm$separation, 1)
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_gte(r$confidence, 0.6)
})

test_that("a resolvable but weakly supported unit is uncertain, with evidence kept", {
  # one hit resolving only to kingdom: unanimous, but almost nothing behind it
  h <- hits(hit("4751;;;;;;", 100, lin = "Fungi;;;;;;"))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "uncertain")
  expect_lt(r$confidence, 0.6)
  # the evidence is still there, which is the point of recording uncertain
  expect_false(is.na(r$lca_taxid))
  expect_equal(r$host_concordant, 0L)
})

test_that("the reason names every term with its value and its weight", {
  h <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  for (k in names(RiboPilot:::SCREEN_WEIGHTS)) expect_match(r$reason, k, fixed = TRUE)
  expect_match(r$reason, "0.30", fixed = TRUE)
  expect_match(r$reason, "0.40", fixed = TRUE)
})

test_that("no identity threshold appears anywhere in the scoring core", {
  # the envelope was deferred after measurement (doc 08 section 3.2.1); a pident
  # comparison here would be a defect, so pin its absence. Deparse the functions
  # rather than reading R/screen.R, because under R CMD check the package is
  # installed and the source tree is not present.
  src <- paste(c(deparse(RiboPilot:::screen_channel_taxonomy),
                 deparse(RiboPilot:::screen_confidence_terms)), collapse = " ")
  expect_false(grepl("pident", src, fixed = TRUE))
})
