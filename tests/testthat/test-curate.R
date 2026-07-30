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
    sample_lineage = FISH, hit_lineage = FISH, hit_is_nuclear_rrna = 1L, chimeric = 0L)
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

test_that("chimeric flag is surfaced", {
  r <- RiboPilot:::curate_rrna_core(good_feats(), RiboPilot:::rrna_params(),
    sample_lineage = FISH, hit_lineage = FISH, hit_is_nuclear_rrna = 1L, chimeric = 1L)
  expect_true("chimeric" %in% r$flags)
})

test_that("lineages_conflict and lineage_known behave", {
  expect_true(RiboPilot:::lineages_conflict("Metazoa; Cnidaria", "Bacteria; Escherichia"))
  expect_false(RiboPilot:::lineages_conflict("Metazoa; Cnidaria", "Metazoa; Chordata"))
  expect_false(RiboPilot:::lineages_conflict("Metazoa", NA))       # unknown -> no conflict
  expect_true(RiboPilot:::lineage_known("Metazoa; Cnidaria"))
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
