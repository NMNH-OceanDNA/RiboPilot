# The two contamination detectors that used to be structurally unable to fire:
# the screen's taxonomy channel, which discarded every non-rDNA hit before
# forming a verdict, and lineages_conflict(), which asked for fully disjoint
# token sets and so could never separate two animals.
#
# Both only ever FLAG. Nothing here arms a drop.

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
FISH_N <- "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio"
CORAL_T <- "33208;6073;6101;3028851;37509;6134;6135"
CORAL_N <- "Metazoa;Cnidaria;Octocorallia;Malacalcyonacea;Alcyoniidae;Sinularia;Sinularia flexibilis"
BIRD_T <- "33208;7711;8782;9126;9205;9207;9208"
BIRD_N <- "Metazoa;Chordata;Aves;Passeriformes;Corvidae;Corvus;Corvus corax"
CONGENER_T <- "33208;7711;186623;7952;2743709;7954;7956"
CONGENER_N <- "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio kyathit"
BACT_T <- "3379134;1224;1236;91347;543;561;562"
BACT_N <- "Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"

ct_hit <- function(tax, bs, cls = "rdna", lin = "K;p;c;o;f;g;s") {
  data.frame(hit_class = cls, bitscore = bs, lineage_taxids = tax,
    lineage = lin, hit_rank = NA_integer_, stringsAsFactors = FALSE)
}

# --- defect 1: the screen's taxonomy channel ---------------------------------

test_that("a unit whose best hits are bacterial 16S is off_target, not no_reference", {
  # classify_hit_title() calls a 16S defline `other` by design (doc 10 s1), and
  # that used to remove every row before the LCA, so a symbiont's rRNA operon
  # contig was reported to the curator as UNJUDGED.
  h <- rbind(ct_hit(BACT_T, 900, "other", BACT_N),
             ct_hit(BACT_T, 890, "other", BACT_N),
             ct_hit(BACT_T, 880, "other", BACT_N))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "off_target")
  expect_equal(r$host_concordant, 0L)
  expect_true(is.na(r$reason_code))
  # and the curator can see the evidence was not nuclear rDNA
  expect_equal(r$hit_class, "other")
  expect_match(r$reason, "non-rDNA records", fixed = TRUE)
})

test_that("a genuinely unscreened unit still reports no_reference", {
  # the distinction the change must preserve: "could not judge" vs "judged, and
  # it is wrong". Each of these is a missing channel, not contrary evidence.
  no_hits <- RiboPilot:::screen_channel_taxonomy(ct_hit(FISH_T, 900)[0, ], FISH_T)
  expect_equal(no_hits$verdict, "no_reference")
  expect_equal(no_hits$reason_code, "no_hits")

  # whole-genome / chromosome records: the title describes the container, so
  # `uninformative` stays excluded and the unit remains unjudged
  asm <- rbind(ct_hit(FISH_T, 900, "uninformative", FISH_N),
               ct_hit(FISH_T, 890, "uninformative", FISH_N),
               ct_hit(FISH_T, 880, "uninformative", FISH_N))
  r <- RiboPilot:::screen_channel_taxonomy(asm, FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_hit_lineages")
  expect_true(is.na(r$hit_class))

  # an unresolved sample lineage keeps its own code
  expect_equal(RiboPilot:::screen_channel_taxonomy(
    ct_hit(FISH_T, 900), NA_character_)$reason_code, "unresolved_sample_lineage")
})

test_that("a host unit is unaffected by non-rDNA hits joining the vote", {
  h <- rbind(ct_hit(FISH_T, 900, "rdna", FISH_N),
             ct_hit(FISH_T, 890, "rdna", FISH_N),
             ct_hit(FISH_T, 880, "rdna", FISH_N),
             ct_hit(FISH_T, 300, "other", FISH_N))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_equal(r$hit_class, "rdna")
})

test_that("placeholder hit lineages are still dropped before the vote", {
  h <- rbind(ct_hit(BACT_T, 900, "other", "Bacteria;;;;;;uncultured bacterium"),
             ct_hit(BACT_T, 890, "other", "Bacteria;;;;;;environmental sample"),
             ct_hit(BACT_T, 880, "other", "Bacteria;;;;;;uncultured organism"))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_hit_lineages")
})

# --- defect 2: lineages_conflict() -------------------------------------------

test_that("two metazoans in different phyla now conflict", {
  # every metazoan lineage contains Metazoa, so the shared-token test read a
  # coral on a fish as clean no matter how distant
  expect_true(RiboPilot:::lineages_conflict(FISH_N, CORAL_N, FISH_T, CORAL_T))
  expect_true(RiboPilot:::lineages_conflict(FISH_N, BACT_N, FISH_T, BACT_T))
  # and the names-only path reaches the same answer
  expect_true(RiboPilot:::lineages_conflict(FISH_N, CORAL_N))
})

test_that("ordinary taxonomic distance within a phylum does NOT conflict", {
  # a bird against a fish is class distance and a congener is genus distance:
  # both are the normal state of an under-represented reference database
  expect_false(RiboPilot:::lineages_conflict(FISH_N, BIRD_N, FISH_T, BIRD_T))
  expect_false(RiboPilot:::lineages_conflict(FISH_N, CONGENER_N, FISH_T, CONGENER_T))
  expect_false(RiboPilot:::lineages_conflict(FISH_N, FISH_N, FISH_T, FISH_T))
})

test_that("an unassigned intermediate rank does not manufacture a conflict", {
  # the sparse-clade case: class and order unknown, phylum shared
  sparse_t <- "33208;7711;;;2743709;7954;7955"
  sparse_n <- "Metazoa;Chordata;unclassified Chordata;unclassified Chordata;Danionidae;Danio;Danio rerio"
  expect_false(RiboPilot:::lineages_conflict(FISH_N, sparse_n, FISH_T, sparse_t))
})

test_that("placeholder names create neither a spurious agreement nor a conflict", {
  ph_n <- "unclassified;uncultured eukaryote;;;;;"
  ph_t <- ";;;;;;"
  # unusable on one side -> unknown, which the caller reports as unchecked
  expect_false(RiboPilot:::lineages_conflict(FISH_N, ph_n, FISH_T, ph_t))
  expect_false(RiboPilot:::lineage_known(ph_n))
  # two placeholders do not agree with each other either
  expect_false(RiboPilot:::lineages_conflict(ph_n, ph_n, ph_t, ph_t))
  # a shared placeholder at one rank cannot rescue two different phyla
  a <- "Metazoa;Chordata;unclassified;unclassified;unclassified;unclassified;unclassified"
  b <- "Metazoa;Cnidaria;unclassified;unclassified;unclassified;unclassified;unclassified"
  expect_true(RiboPilot:::lineages_conflict(a, b))
})

test_that("a taxid slot whose name is a placeholder is not compared", {
  # taxonkit leaves a filled rank's taxid empty, but a real taxid can still sit
  # under a junk name ("environmental samples"); it must not create agreement
  a <- "Metazoa;Chordata;environmental samples;;;;"
  b <- "Metazoa;Cnidaria;environmental samples;;;;"
  expect_true(RiboPilot:::lineages_conflict(a, b, "33208;7711;48479;;;;", "33208;6073;48479;;;;"))
})

test_that("the contaminant flag fires on a different-phylum hit and only flags", {
  feats <- data.frame(
    element = c("18S", "ITS1", "5.8S", "ITS2", "28S"),
    type = c("rRNA", "misc_RNA", "rRNA", "misc_RNA", "rRNA"),
    start = c(1L, 1801L, 2001L, 2201L, 2401L),
    stop = c(1800L, 2000L, 2160L, 2400L, 5400L),
    length = c(1800L, 200L, 160L, 200L, 3000L),
    strand = "+", stringsAsFactors = FALSE)
  r <- RiboPilot:::curate_rrna_core(feats, RiboPilot:::rrna_params(),
    sample_lineage = FISH_N, hit_lineage = CORAL_N,
    sample_taxids = FISH_T, hit_taxids = CORAL_T, hit_is_nuclear_rrna = 1L)
  expect_true("contaminant" %in% r$flags)
  expect_false("contaminant_unchecked" %in% r$flags)

  clean <- RiboPilot:::curate_rrna_core(feats, RiboPilot:::rrna_params(),
    sample_lineage = FISH_N, hit_lineage = BIRD_N,
    sample_taxids = FISH_T, hit_taxids = BIRD_T, hit_is_nuclear_rrna = 1L)
  expect_false("contaminant" %in% clean$flags)
})
