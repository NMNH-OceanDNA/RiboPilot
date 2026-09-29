# The screen's taxonomy channel. It answers "do this unit's hits agree on a
# taxon, and is that taxon the host's". It needs NO identity threshold, which is
# why it survives the envelope being deferred (see doc 08 section 3.2.1).

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
FUNGI_T <- "4751;4890;147545;5042;1131492;5052;5062"
COPE_T <- "33208;6656;6854;6843;1;2;3"

# one hit row as blast_hits stores it. hit_rank is the per-subject rank, so two
# rows sharing one are two HSPs of ONE accession; the default gives every row
# its own accession.
hit <- function(tax, bs, cls = "rdna", lin = "K;p;c;o;f;g;s", hr = NA_integer_) {
  data.frame(hit_class = cls, bitscore = bs, lineage_taxids = tax,
    lineage = lin, hit_rank = hr, stringsAsFactors = FALSE)
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

test_that("finding A: a congener hit is concordant, not discordant, on the margin", {
  # one exact species match plus two same-genus relatives: all three share the
  # sample's lineage above species, so none of them opposes the LCA. Testing
  # the sample's DEEPEST taxid instead of ANY shared taxid used to count the
  # congeners as discordant and score this uncertain.
  cong1 <- "33208;7711;186623;7952;2743709;7954;90001"
  cong2 <- "33208;7711;186623;7952;2743709;7954;90002"
  h <- hits(hit(FISH_T, 900), hit(cong1, 500), hit(cong2, 480))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_gte(r$confidence, 0.6)
})

test_that("finding B: a lone hit cannot decide a scaffold on its own", {
  # ONE fungal hit is discordant and, being alone, unopposed by definition; a
  # count floor keeps that from reaching off_target unopposed
  h1 <- hits(hit(FUNGI_T, 900))
  r1 <- RiboPilot:::screen_channel_taxonomy(h1, FISH_T)
  expect_equal(r1$verdict, "uncertain")

  # three discordant hits clear the floor and still reach off_target
  h3 <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  r3 <- RiboPilot:::screen_channel_taxonomy(h3, FISH_T)
  expect_equal(r3$verdict, "off_target")
})

test_that("finding C: class_for matches the whole taxid path, not rank slot 1", {
  BACT_T <- "2;1224;1236;91347;543;561;562"
  expect_equal(RiboPilot:::screen_class_for(BACT_T, NA_character_), "bacterial")
  expect_equal(RiboPilot:::screen_class_for(FUNGI_T, NA_character_), "fungal")
  expect_equal(RiboPilot:::screen_class_for(FISH_T, NA_character_), "metazoan_other")
  # taxon label wins even though the path itself is metazoan
  expect_equal(RiboPilot:::screen_class_for(FISH_T, "Homo sapiens"), "human")

  # the call site, not just the helper: drive bacterial and fungal through
  # screen_channel_taxonomy() too
  hb <- hits(hit(BACT_T, 900), hit(BACT_T, 890), hit(BACT_T, 880))
  rb <- RiboPilot:::screen_channel_taxonomy(hb, FISH_T)
  expect_equal(rb$verdict, "off_target")
  expect_equal(rb$class, "bacterial")

  hf <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  rf <- RiboPilot:::screen_channel_taxonomy(hf, FISH_T)
  expect_equal(rf$verdict, "off_target")
  expect_equal(rf$class, "fungal")
})

test_that("a different-phylum contaminant is not waved through as concordant", {
  # sharing only the kingdom slot means "both are animals", which every
  # metazoan contaminant satisfies. Three exact host hits plus two hits that
  # share nothing below kingdom (a real different phylum) must split the
  # vote, not read as unanimous.
  DIFFPHYLUM_T <- "33208;6656;6854;6843;1;2;3"
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880),
            hit(DIFFPHYLUM_T, 500), hit(DIFFPHYLUM_T, 480))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "uncertain")
  expect_false(is.na(r$margin))
  expect_lt(r$terms$separation, 1)
})

test_that("three HSPs of ONE accession are one vote, not three", {
  # blast_hits stores one row per HSP and blast.nf runs blastn with
  # -max_hsps 5, so an accession casts up to five votes unless the HSPs are
  # collapsed first. COUNT(DISTINCT hit_rank) here is 1: whatever this unit
  # knows, it knows from a single GenBank subject, which is not enough to
  # satisfy either the LCA's min_hits or the unanimity floor.
  h <- hits(hit(FUNGI_T, 880, hr = 1), hit(FUNGI_T, 900, hr = 1),
            hit(FUNGI_T, 890, hr = 1))
  expect_equal(length(unique(h$hit_rank)), 1L)

  kept <- RiboPilot:::screen_one_vote_per_subject(h)
  expect_equal(nrow(kept), 1L)
  expect_equal(kept$bitscore, 900)   # the best-bitscore HSP is the one kept

  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$lca_n_hits, 1L)     # one distinct accession, one vote
  expect_equal(r$terms$separation, 0)
  expect_lt(r$confidence, 0.6)
  expect_equal(r$verdict, "uncertain")
})

test_that("three DISTINCT accessions still produce a confident verdict", {
  # the HSP collapse must not disable the feature: the same three bitscores
  # spread across three subjects are three votes and do decide the unit
  h <- hits(hit(FUNGI_T, 900, hr = 1), hit(FUNGI_T, 890, hr = 2),
            hit(FUNGI_T, 880, hr = 3))
  expect_equal(length(unique(h$hit_rank)), 3L)
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$lca_n_hits, 3L)
  expect_equal(r$verdict, "off_target")
  expect_gte(r$confidence, 0.6)
})

test_that("two distinct accessions do not clear the gate on unanimity alone", {
  # unanimity sets separation to 1, worth 0.40 of the 0.60 gate, so two rows
  # agreeing at species rank used to reach 0.76 with nothing else behind them
  h <- hits(hit(FUNGI_T, 900, hr = 1), hit(FUNGI_T, 890, hr = 2))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$lca_n_hits, 2L)
  expect_equal(r$terms$separation, 0)
  expect_lt(r$confidence, 0.6)
  expect_equal(r$verdict, "uncertain")
})

test_that("three hits of one congener are on_target, not off_target", {
  # the congener's own species is absent from the sample's path, so testing
  # `lca_taxid %in% samp` calls a clean unit off_target at 0.790. That is the
  # too-strict rule the margin already discarded, and it fires hardest on the
  # under-described taxa this tool exists to serve.
  cong <- "33208;7711;186623;7952;2743709;7954;90001"
  h <- hits(hit(cong, 900, hr = 1), hit(cong, 890, hr = 2), hit(cong, 880, hr = 3))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$lca_taxid, 90001L)          # the LCA is NOT on the sample's path
  expect_false(r$lca_taxid %in% RiboPilot:::parse_rank_taxids(FISH_T))
  expect_equal(r$host_concordant, 1L)        # but its ancestry is
  expect_equal(r$verdict, "on_target")
})

test_that("a same-phylum different-class hit is still concordant", {
  # phylum is the shallowest thing the concordance test accepts, and for a
  # sparse clade it may be the nearest reference GenBank has
  otherclass <- "33208;7711;80000;80001;80002;80003;80004"
  h <- hits(hit(otherclass, 900, hr = 1), hit(otherclass, 890, hr = 2),
            hit(otherclass, 880, hr = 3))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$host_concordant, 1L)
  expect_equal(r$verdict, "on_target")
})

test_that("a kingdom-only sample lineage is no_reference, not off_target", {
  # dropping the shallowest slot leaves nothing to compare at phylum or
  # deeper. Read as discordant, every hit becomes unanimous and clean host
  # data scores off_target at 0.790.
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, "33208;;;;;;")
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_sample_lineage")
  expect_true(is.na(r$host_concordant))
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

mk_screen_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage_taxids TEXT, reference_support INT)")
  # sequence/length are in the base assemblies DDL and the variant pass reads
  # them, so the fixture must carry them or it tests a schema that cannot exist
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT,
    sequence TEXT, length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
    element TEXT, start INT, stop INT, strand TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',?,NULL)", params = list(FISH_T))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,2,NULL,NULL)")
  # a second sample, present so per-sample driver calls can be pinned as scoped
  # to their own sample rather than accidentally scanning the whole project
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2',?,NULL)", params = list(FISH_T))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S2',1,1,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S2',1,2,NULL,NULL)")
  con
}
add_hits <- function(con, scaf, tax, n = 3, bs = 900,
                      lin = "Eukaryota;Chordata;Actinopteri;o;f;g;s", id = "S1",
                      rd = rank_distance(tax, FISH_T)) {
  for (i in seq_len(n)) {
    DBI::dbExecute(con,
      "INSERT INTO blast_hits
         (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class,rank_distance)
       VALUES (?,1,?,?,1,?,?,?,'rdna',?)",
      params = list(id, scaf, i, bs - i, tax, lin, rd))
  }
}

test_that("screen_one writes one row carrying the verdict and its evidence", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE scaffold=1")
  expect_equal(nrow(r), 1L)
  expect_equal(r$screen_verdict, "on_target")
  expect_equal(r$host_concordant, 1L)
  expect_equal(r$screen_version, RiboPilot:::SCREEN_VERSION)
  expect_true(!is.na(r$screen_reason) && nzchar(r$screen_reason))
})

test_that("a unit searched again is screened again", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "DELETE FROM blast_hits")
  add_hits(con, 1, FUNGI_T, lin = "Eukaryota;Ascomycota;Sordariomycetes;o;f;g;s")
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "new", 3L)
  screen_sample(con, "S1", variants = FALSE)
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_verdict v FROM screen WHERE ID='S1' AND scaffold=1")$v,
    "off_target")
})

test_that("a unit with no hits gets a no_reference row rather than none", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  screen_one(con, "S1", 1, 2)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict FROM screen WHERE scaffold=2")
  expect_equal(r$screen_verdict, "no_reference")
})

test_that("a stored row at the current version is left alone unless forced", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "UPDATE screen SET screen_reason='SENTINEL' WHERE scaffold=1")
  screen_one(con, "S1", 1, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_reason r FROM screen WHERE scaffold=1")$r,
    "SENTINEL")
  screen_one(con, "S1", 1, 1, force = TRUE)
  expect_false(identical(
    DBI::dbGetQuery(con, "SELECT screen_reason r FROM screen WHERE scaffold=1")$r, "SENTINEL"))
})

test_that("a row stored at a lower version is recomputed without force", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "UPDATE screen SET screen_version=0, screen_reason='OLD' WHERE scaffold=1")
  screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_version v, screen_reason s FROM screen WHERE scaffold=1")
  expect_equal(r$v, RiboPilot:::SCREEN_VERSION)
  expect_false(identical(r$s, "OLD"))
})

test_that("a stored no_reference row is recomputed once evidence arrives, without force", {
  # no_reference records MISSING evidence, and the evidence routinely arrives
  # later: BLAST disabled on the first pass, or the sample's lineage resolved
  # afterwards. Skipping on version alone made the verdict permanent.
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = NULL WHERE ID = 'S1'")
  screen_one(con, "S1", 1, 1)
  r0 <- DBI::dbGetQuery(con,
    "SELECT screen_verdict v, screen_version ver FROM screen WHERE scaffold=1")
  expect_equal(r0$v, "no_reference")
  expect_equal(r0$ver, RiboPilot:::SCREEN_VERSION)   # already current, so only
                                                    # the verdict can trigger it

  add_hits(con, 1, FISH_T)
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = ? WHERE ID = 'S1'",
    params = list(FISH_T))
  screen_one(con, "S1", 1, 1)
  r1 <- DBI::dbGetQuery(con,
    "SELECT screen_verdict v, lca_n_hits n FROM screen WHERE scaffold=1")
  expect_equal(r1$v, "on_target")
  expect_equal(r1$n, 3L)
})

test_that("screen_sample covers every unit", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T); add_hits(con, 2, FISH_T)
  n <- screen_sample(con, "S1")
  expect_equal(n, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen")$n, 2L)
})

test_that("screen_sample touches only its own sample, not the whole project", {
  # a prior sub-project shipped an ungated whole-project rescan here, caught
  # only by measurement (O(units^2) on a live run); this pins scope by
  # construction rather than relying on that happening again. S1 and S2 are
  # deliberately ASYMMETRIC -- different unit counts, and S2 owns a (path,
  # scaffold) tuple S1 has none of -- so an unfiltered query cannot land on
  # the right rows by re-iterating a shared shape.
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S2',2,5,NULL,NULL)")
  add_hits(con, 1, FISH_T, id = "S1")
  add_hits(con, 2, FISH_T, id = "S1")
  add_hits(con, 1, FISH_T, id = "S2")
  add_hits(con, 2, FISH_T, id = "S2")

  screen_sample(con, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen WHERE ID='S1'")$n, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen WHERE ID='S2'")$n, 0L)

  s1_before <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  screen_sample(con, "S2")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen WHERE ID='S2'")$n, 3L)
  s1_after <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(s1_before, s1_after)
})

test_that("the scoring core reads its weights from the params sublist", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  base <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  # a ruleset that weighted support alone would score this differently; if the
  # core ignored the argument, both calls would agree
  tuned <- RiboPilot:::screen_channel_taxonomy(h, FISH_T,
    params = list(support = 1, resolution = 0, separation = 0, min_confidence = 0.6))
  expect_false(isTRUE(all.equal(base$confidence, tuned$confidence)))
  expect_equal(tuned$confidence, 0.3)  # lca_n_hits 3 -> support 0.3, weight 1
})

test_that("SCREEN_VERSION is 5 so earlier rows recompute", {
  # bumped by 3c twice: once for the variant columns, once when the coverage and
  # shared-element guards changed what an isomer means; then again when `other`
  # hits joined the vote
  expect_equal(RiboPilot:::SCREEN_VERSION, 5L)
})

test_that("an unanswered local search suggests a remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 0L, source = "local")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM screen")
  expect_equal(r$screen_verdict[1], "no_reference")
  expect_match(r$screen_reason[1], "bundled")
  expect_match(r$screen_reason[1], "NCBI")
})

test_that("an unanswered remote search does not suggest another remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 0L, source = "remote")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_reason FROM screen")
  expect_false(grepl("remote search", r$screen_reason[1]))
  # pins the wording itself, not just the absence of "remote search": a clause
  # later widened to fire on remote too would still trip this
  expect_equal(r$screen_reason[1], "no BLAST hits are stored for this scaffold")
})

test_that("a local search whose hits all fail to resolve suggests a remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  # every hit's deepest lineage field is a placeholder, so all are filtered
  # out: reason_code "unresolved_hit_lineages", not the zero-hit short circuit
  add_hits(con, 1, FISH_T, lin = "Eukaryota;Chordata;Actinopteri;o;f;g;uncultured organism")
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 3L, source = "local")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM screen")
  expect_equal(r$screen_verdict[1], "no_reference")
  expect_match(r$screen_reason[1], "bundled")
  expect_match(r$screen_reason[1], "NCBI")
})

test_that("a remote search whose hits all fail to resolve does not suggest another remote one", {
  con <- mk_screen_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  add_hits(con, 1, FISH_T, lin = "Eukaryota;Chordata;Actinopteri;o;f;g;uncultured organism")
  RiboPilot:::mark_unit_searched(con, "S1", 1, 1, "abc", 3L, source = "remote")
  RiboPilot:::screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict, screen_reason FROM screen")
  expect_equal(r$screen_verdict[1], "no_reference")
  expect_false(grepl("remote search", r$screen_reason[1]))
  expect_equal(r$screen_reason[1], "no hit resolved to a usable lineage")
})
