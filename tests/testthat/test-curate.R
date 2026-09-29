# Curation flag logic (pure): structural, length, identity, truncation,
# contamination, chimera.

good_feats <- function() {
  f <- data.frame(
    element = c("18S", "ITS1", "5.8S", "ITS2", "28S"),
    type = c("rRNA", "misc_RNA", "rRNA", "misc_RNA", "rRNA"),
    start = c(1, 1801, 2201, 2361, 2701), stop = c(1800, 2200, 2360, 2700, 6000),
    strand = "+", model_coverage = c(100, NA, 100, NA, 100), stringsAsFactors = FALSE)
  f$length <- f$stop - f$start + 1
  f
}
FISH <- "Metazoa; Chordata; Actinopteri; Gadus"

test_that("a clean feature model curates with no flags", {
  r <- RiboPilot:::curate_rrna_core(good_feats(), RiboPilot:::rrna_params("actinopterygii"),
    sample_lineage = FISH, hit_lineage = FISH, hit_is_nuclear_rrna = 1L)
  expect_equal(r$n_flags, 0)
  expect_true(all(r$features$status == "ok"))
})

test_that("missing anchors and out-of-order are flagged", {
  f <- good_feats()[good_feats()$element != "28S", ]
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())
  expect_true("missing_28S" %in% r$flags)

  f2 <- good_feats()
  f2$start <- rev(f2$start); f2$stop <- rev(f2$stop)  # scramble order
  r2 <- RiboPilot:::curate_rrna_core(f2, RiboPilot:::rrna_params())
  expect_true("out_of_order" %in% r2$flags)
})

test_that("wrong_strand fires on a minus-strand feature (not dead code)", {
  f <- good_feats(); f$strand <- "-"
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    sample_lineage = FISH, hit_lineage = FISH)
  flat <- unlist(strsplit(r$features$flags[nzchar(r$features$flags)], ","))
  expect_true("wrong_strand" %in% flat)
})

test_that("truncated anchor is flagged from model coverage", {
  f <- good_feats(); f$model_coverage[f$element == "18S"] <- 70
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())
  fl <- r$features$flags[r$features$element == "18S"]
  expect_true(grepl("truncated", fl))
})

test_that("a truncated feature is partial only at the end that meets no neighbour", {
  f <- good_feats(); f$model_coverage[f$element == "18S"] <- 70
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())
  x <- r$features[r$features$element == "18S", ]
  expect_equal(c(x$partial_start, x$partial_stop), c(1L, NA))  # 3' end meets ITS1
  expect_true(all(is.na(r$features$partial_start[r$features$element != "18S"])))

  # meets neighbours at both ends: cannot tell which is short, so both
  f$model_coverage[f$element == "18S"] <- 100; f$model_coverage[f$element == "5.8S"] <- 50
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())
  x <- r$features[r$features$element == "5.8S", ]
  expect_equal(c(x$partial_start, x$partial_stop), c(1L, 1L))

  # a gap wider than the snap distance does not count as meeting; strand flips ends
  f <- good_feats()[good_feats()$element %in% c("18S", "ITS1"), ]
  f$model_coverage[1] <- 70; f$start[2] <- 1807; f$strand <- "-"
  x <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())$features[1, ]
  expect_equal(c(x$partial_start, x$partial_stop), c(1L, 1L))
  f$start[2] <- 1806
  x <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())$features[1, ]
  expect_equal(c(x$partial_start, x$partial_stop), c(NA, 1L))

  # a feature nested inside the truncated one is not a neighbour
  f <- good_feats()[good_feats()$element %in% c("18S", "ITS1"), ]
  f$model_coverage[1] <- 70; f$start[2] <- 1500; f$stop[2] <- 1600
  x <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())$features[1, ]
  expect_equal(c(x$partial_start, x$partial_stop), c(1L, 1L))

  # a neighbour overlapping by a few bp still meets it
  f$start[2] <- 1798; f$stop[2] <- 2200
  x <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())$features[1, ]
  expect_equal(c(x$partial_start, x$partial_stop), c(1L, NA))
})

test_that("curate_one fills partial ends once and keeps the curator's choice", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations")
  f <- good_feats()
  f$model_coverage[f$element == "28S"] <- 60
  DBI::dbAppendTable(con, "annotations", cbind(ID = "S1", path = 1L, scaffold = 1L, f))
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  p <- function() DBI::dbGetQuery(con,
    "SELECT partial_start s, partial_stop e FROM annotations WHERE element = '28S'")
  expect_equal(unlist(p()), c(s = NA, e = 1L))
  DBI::dbExecute(con, "UPDATE annotations SET partial_stop = 0 WHERE element = '28S'")
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  expect_equal(unlist(p()), c(s = NA, e = 0L))
})

test_that("low_identity respects coverage and per-element floors", {
  f <- good_feats()
  # partial reference (low coverage): low pid must NOT flag
  r1 <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    ref_pid = c("18S" = 3), ref_cov = c("18S" = 0.1))
  expect_false(grepl("low_identity", r1$features$flags[r1$features$element == "18S"]))
  # well-covered, genuinely divergent 28S: flag
  r2 <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    ref_pid = c("28S" = 60), ref_cov = c("28S" = 0.9))
  expect_true(grepl("low_identity", r2$features$flags[r2$features$element == "28S"]))
})

test_that("contamination checks use the screen, not silence", {
  f <- good_feats()
  OCTO <- "Metazoa; Cnidaria; Anthozoa; Octocorallia; Xenia"
  SYMB <- "Eukaryota; Sar; Alveolata; Dinophyceae; Symbiodinium"
  # symbiont hit -> contaminant
  r1 <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    sample_lineage = OCTO, hit_lineage = SYMB, hit_is_nuclear_rrna = 1L)
  expect_true("contaminant" %in% r1$flags)
  # unknown lineage -> contaminant_unchecked, not clean
  r2 <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    sample_lineage = OCTO, hit_lineage = NA, hit_is_nuclear_rrna = 1L)
  expect_true("contaminant_unchecked" %in% r2$flags)
  # top hit is not nuclear rDNA -> not_nuclear_rdna
  r3 <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params(),
    sample_lineage = OCTO, hit_lineage = OCTO, hit_is_nuclear_rrna = 0L)
  expect_true("not_nuclear_rdna" %in% r3$flags)
})

test_that("lineages_conflict and lineage_known behave", {
  expect_true(RiboPilot:::lineages_conflict("Metazoa; Cnidaria", "Bacteria; Escherichia"))
  # two animals in different phyla: the case the shared-token test could not fire on
  expect_true(RiboPilot:::lineages_conflict("Metazoa; Cnidaria", "Metazoa; Chordata"))
  expect_false(RiboPilot:::lineages_conflict("Metazoa", NA))       # unknown -> no conflict
  expect_true(RiboPilot:::lineage_known("Metazoa; Cnidaria"))
  expect_false(RiboPilot:::lineage_known("Metazoa"))               # nothing below kingdom
  expect_false(RiboPilot:::lineage_known(NA))
  expect_false(RiboPilot:::lineage_known(""))
})

test_that("the extracted ruleset resolver keeps curate_one's precedence", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, curate_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE curate_opts(curate_opts TEXT, ruleset TEXT)")
  DBI::dbExecute(con, "INSERT INTO curate_opts VALUES('default','metazoa_generic')")
  DBI::dbExecute(con, "INSERT INTO curate_opts VALUES('fish','actinopterygii')")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,'default')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Eukaryota; Metazoa; Cnidaria; Scyphozoa')")

  # explicit argument wins over everything
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1, "asteroidea")$ruleset,
    "asteroidea")
  # a generic curate_opts falls through to the lineage
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1)$ruleset, "scyphozoa")
  # a non-generic curate_opts is honoured as-is
  DBI::dbExecute(con, "UPDATE annotate SET curate_opts='fish'")
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1)$ruleset, "actinopterygii")
  # a lineage matching no clade records the fallback rather than hiding it
  DBI::dbExecute(con, "UPDATE annotate SET curate_opts='default'")
  DBI::dbExecute(con, "UPDATE samples SET lineage='Eukaryota; Metazoa; Porifera'")
  r <- RiboPilot:::resolve_ruleset(con, "S1", 1, 1)
  expect_equal(r$ruleset, "metazoa_generic")
  expect_equal(r$note, "no_clade_ruleset")
})

test_that("a gap between neighbouring features is flagged; one between non-neighbours is not", {
  f <- good_feats(); f$start[f$element == "ITS1"] <- 1805
  r <- RiboPilot:::curate_rrna_core(f, RiboPilot:::rrna_params())
  expect_true("gap" %in% r$flags)
  expect_match(r$features$flags[r$features$element == "ITS1"], "gap")
  f2 <- good_feats()[good_feats()$element %in% c("18S", "5.8S", "28S"), ]
  expect_false("gap" %in% RiboPilot:::curate_rrna_core(f2, RiboPilot:::rrna_params())$flags)
})

test_that("curation writes its flags to curate_flags and keeps the curator's note", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = 'check 3prime end with PCR'")
  RiboPilot::curate_one(con, "S1", 1L, 1L)
  a <- DBI::dbGetQuery(con, "SELECT annotate_notes, curate_flags, n_flags FROM annotate WHERE ID = 'S1'")
  expect_equal(a$annotate_notes, "check 3prime end with PCR")
  expect_match(a$curate_flags, "missing_18S", fixed = TRUE)
  expect_gt(a$n_flags, 0)
})
