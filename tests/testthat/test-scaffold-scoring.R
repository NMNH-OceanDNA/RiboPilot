# Candidate-path scoring and auto-select. This is what decides which assembly the
# curator is shown first, so the sub-scores are pinned individually: a composite
# that happens to rank correctly can still be built from a broken term.

test_that("parse_gfa_depth reads depth as RC/LN per segment", {
  gfa <- withr::local_tempfile(fileext = ".gfa")
  writeLines(c(
    "H\tVN:Z:1.0",
    "S\t101\tACGTACGTAC\tLN:i:10\tRC:i:1000",
    "S\t102\tACGTA\tLN:i:5\tRC:i:100",
    "L\t101\t+\t102\t+\t0M"), gfa)
  d <- RiboPilot:::parse_gfa_depth(gfa)
  expect_equal(d$segment, c("101", "102"))
  expect_equal(d$depth, c(100, 20))
})

test_that("parse_gfa_depth falls back to sequence length when LN is absent", {
  gfa <- withr::local_tempfile(fileext = ".gfa")
  writeLines(c("S\t1\tACGTACGTAC\tRC:i:50"), gfa)
  expect_equal(RiboPilot:::parse_gfa_depth(gfa)$depth, 5)
})

test_that("parse_gfa_depth returns an empty frame for a GFA with no segments", {
  gfa <- withr::local_tempfile(fileext = ".gfa")
  writeLines("H\tVN:Z:1.0", gfa)
  expect_equal(nrow(RiboPilot:::parse_gfa_depth(gfa)), 0)
})

test_that("scaffold_depth is the length-weighted mean over traversed segments", {
  gd <- data.frame(segment = c("11", "22"), length = c(100, 900),
    rc = c(1, 1), depth = c(10, 110), stringsAsFactors = FALSE)
  # weighted: (100*10 + 900*110) / 1000
  expect_equal(RiboPilot:::scaffold_depth("scaffold_1--11-,22+", gd), 100)
  # a header naming no known segment has no depth, rather than depth 0
  expect_true(is.na(RiboPilot:::scaffold_depth("scaffold_2--99+", gd)))
})

test_that("taxonomic_match scores the deepest agreeing rank", {
  f <- RiboPilot:::taxonomic_match
  lin <- "Metazoa;Chordata;Actinopteri;Perciformes;Sebastidae;Sebastes;Sebastes alutus"
  expect_equal(f(lin, lin), 1)                                  # species
  expect_equal(f("Metazoa;Chordata;Actinopteri;;;;", lin), 3 / 7)  # class
  expect_equal(f("Viridiplantae;Streptophyta;;;;;", lin), 0)
  expect_equal(f(NA_character_, lin), 0)
  expect_equal(f(lin, NA_character_), 0)
  expect_equal(f(toupper(lin), lin), 1)                         # case-insensitive
})

test_that("a gap at one rank no longer discards agreement at every deeper rank", {
  # The prefix walk this replaced broke at the missing class and scored 2/7,
  # calling a congeneric hit almost unrelated. Under-described invertebrates hit
  # this constantly, which is exactly backwards.
  lin <- "Metazoa;Chordata;Actinopteri;Perciformes;Sebastidae;Sebastes;Sebastes alutus"
  gap <- "Metazoa;Chordata;;Perciformes;Sebastidae;Sebastes;Sebastes alutus"
  expect_equal(RiboPilot:::taxonomic_match(gap, lin), 1)
})

test_that("taxonomic_match ignores taxonkit's unclassified placeholders", {
  # `reformat -F` writes "unclassified <parent>" into a rank it had to fill; two
  # unrelated taxa must not be scored as agreeing because both say unclassified.
  a <- "Metazoa;Cnidaria;unclassified Cnidaria;unclassified Cnidaria;X;Y;Y z"
  b <- "Metazoa;Chordata;unclassified Chordata;unclassified Chordata;P;Q;Q r"
  expect_equal(RiboPilot:::taxonomic_match(a, b), 1 / 7)   # agrees at kingdom only
})

test_that("taxonomic_match compares taxids when asked", {
  f <- RiboPilot:::taxonomic_match
  a <- "33208;7711;186623;7952;2743709;7954;7955"
  b <- "33208;7711;186623;7952;2743709;7954;9999"
  expect_equal(f(a, b, use_taxids = TRUE), 6 / 7)          # same genus
  expect_equal(f(a, a, use_taxids = TRUE), 1)
  expect_equal(f("33208;;;;;;", "9999;;;;;;", use_taxids = TRUE), 0)
})

test_that("score_candidates prefers taxids over names when both sides have them", {
  # names disagree (a synonym), taxids agree: the taxid answer must win
  cand <- data.frame(path = 1, depth = 1, blast_pident = 90,
    blast_lineage = "Metazoa;Chordata;;;;;Synonym name",
    blast_lineage_taxids = "33208;7711;186623;7952;2743709;7954;7955",
    anchors = "18S", length = 8000, stringsAsFactors = FALSE)
  by_name <- score_candidates(cand,
    sample_lineage = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio")
  by_taxid <- score_candidates(cand,
    sample_lineage = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
    sample_taxids = "33208;7711;186623;7952;2743709;7954;7955")
  expect_equal(by_name$s_tax, round(2 / 7, 3))
  expect_equal(by_taxid$s_tax, 1)
})

test_that("score_candidates ranks, flags one auto_select, and returns sorted rows", {
  cand <- data.frame(
    path = 1:3,
    depth = c(10, 100, 50),
    blast_pident = c(80, 99, 90),
    blast_lineage = c("A; B", "A; B; C", "A"),
    anchors = c("18S", "18S,5.8S,28S", "18S,28S"),
    length = c(2000, 8000, 8000),
    stringsAsFactors = FALSE)
  out <- score_candidates(cand, sample_lineage = "A; B; C")
  expect_equal(out$rank, 1:3)                 # returned in rank order
  expect_equal(sum(out$auto_select), 1L)
  expect_equal(out$path[1], 2L)               # best on every term
  expect_true(all(c("s_tax", "s_identity", "s_anchors", "s_depth", "s_length",
    "score", "rank", "auto_select") %in% names(out)))
})

test_that("score_candidates sub-scores are each on [0, 1]", {
  cand <- data.frame(path = 1:2, depth = c(10, 100), blast_pident = c(50, 100),
    blast_lineage = c("A", "A; B; C"), anchors = c("", "18S,5.8S,28S"),
    length = c(500, 8000), stringsAsFactors = FALSE)
  out <- score_candidates(cand, sample_lineage = "A; B; C")
  for (col in c("s_tax", "s_identity", "s_anchors", "s_depth", "s_length")) {
    expect_true(all(out[[col]] >= 0 & out[[col]] <= 1), info = col)
  }
  expect_equal(out$s_anchors[out$path == 2], 1)
  expect_equal(out$s_anchors[out$path == 1], 0)
})

test_that("s_depth is relative to the deepest candidate, not absolute", {
  mk <- function(dep) data.frame(path = 1:2, depth = dep, blast_pident = c(90, 90),
    blast_lineage = c("A", "A"), anchors = c("18S", "18S"), length = c(8000, 8000),
    stringsAsFactors = FALSE)
  shallow <- score_candidates(mk(c(5, 10)), sample_lineage = "A")
  deep <- score_candidates(mk(c(500, 1000)), sample_lineage = "A")
  # the deepest always scores 1, so absolute coverage carries no information here
  expect_equal(shallow$s_depth, deep$s_depth)
  expect_equal(max(shallow$s_depth), 1)
})

test_that("s_length is 1 inside the expected window and decays linearly outside", {
  sc <- function(L, expected = c(4000, 12000)) {
    out <- score_candidates(
      data.frame(path = 1, depth = 1, blast_pident = 90, blast_lineage = "A",
        anchors = "18S", length = L, stringsAsFactors = FALSE),
      sample_lineage = "A", expected_length = expected)
    out$s_length
  }
  expect_equal(sc(4000), 1)
  expect_equal(sc(8000), 1)
  expect_equal(sc(12000), 1)
  expect_equal(sc(2000), 0)          # a full 50% below the floor
  expect_equal(sc(3000), 0.5)
  expect_equal(sc(18000), 0)         # a full 50% above the ceiling
  expect_equal(sc(15000), 0.5)
  expect_equal(sc(NA_integer_), 0)
})

test_that("score_candidates handles an empty candidate set", {
  empty <- data.frame(path = integer(0), depth = numeric(0), blast_pident = numeric(0),
    blast_lineage = character(0), anchors = character(0), length = integer(0),
    stringsAsFactors = FALSE)
  expect_equal(nrow(score_candidates(empty)), 0)
})

test_that("missing evidence scores zero rather than propagating NA", {
  cand <- data.frame(path = 1:2, depth = c(NA_real_, NA_real_),
    blast_pident = c(NA_real_, 90), blast_lineage = c(NA_character_, "A"),
    anchors = c(NA_character_, "18S"), length = c(NA_integer_, 8000L),
    stringsAsFactors = FALSE)
  out <- score_candidates(cand, sample_lineage = "A")
  expect_false(any(is.na(out$score)))
  expect_equal(out$score[out$path == 1], 0)
})

test_that("weights change the ranking they are supposed to change", {
  cand <- data.frame(path = 1:2, depth = c(100, 1), blast_pident = c(50, 100),
    blast_lineage = c("A", "A"), anchors = c("18S", "18S"), length = c(8000, 8000),
    stringsAsFactors = FALSE)
  depth_led <- score_candidates(cand, sample_lineage = "A",
    weights = c(tax = 0, identity = 0, anchors = 0, depth = 1, length = 0))
  ident_led <- score_candidates(cand, sample_lineage = "A",
    weights = c(tax = 0, identity = 1, anchors = 0, depth = 0, length = 0))
  expect_equal(depth_led$path[1], 1L)
  expect_equal(ident_led$path[1], 2L)
})
