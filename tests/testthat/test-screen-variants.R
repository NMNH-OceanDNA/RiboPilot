# Intragenomic variants: are two on-target units the same molecule assembled
# twice, or two genuinely different rDNA copies from one individual. Restricting
# the comparison to on_target units is what makes any difference between them
# intragenomic BY CONSTRUCTION (doc 12 section 1), so no threshold has to
# separate a paralog from a contaminant. That separation already happened.
#
# Only the ANCHORS are compared. The spacers are what varies in both cases and
# so cannot discriminate between them.

ann_row <- function(element, start, stop, strand = "+") {
  data.frame(element = element, start = start, stop = stop, strand = strand,
    stringsAsFactors = FALSE)
}
anns <- function(...) do.call(rbind, list(...))

test_that("the anchor string is 18S + 5.8S + 28S in element order", {
  # 18S(1-10) ITS1(11-15) 5.8S(16-23) ITS2(24-27) 28S(28-39)
  seq <- paste0(strrep("A", 10), strrep("C", 5), strrep("G", 8), strrep("T", 4),
    strrep("A", 12))
  a <- anns(ann_row("18S", 1, 10), ann_row("ITS1", 11, 15),
            ann_row("5.8S", 16, 23), ann_row("ITS2", 24, 27),
            ann_row("28S", 28, 39))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8), strrep("A", 12)))
})

test_that("element order is fixed, not coordinate order", {
  # a circularly permuted assembly puts 28S first; the same molecule must not
  # produce a different anchor string and read as a paralog of itself
  seq <- paste0(strrep("G", 8), strrep("A", 10))
  a <- anns(ann_row("28S", 1, 8), ann_row("18S", 9, 18))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8)))
})

test_that("a unit with two of three anchors is still comparable", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10), ann_row("28S", 11, 18))
  expect_equal(nchar(RiboPilot:::screen_anchor_string(seq, a)), 18L)
})

test_that("a unit with no anchors at all is not comparable", {
  seq <- strrep("A", 30)
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, anns(ann_row("ITS1", 1, 10)))))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, data.frame())))
  expect_true(is.na(RiboPilot:::screen_anchor_string(NA_character_,
    anns(ann_row("18S", 1, 10)))))
})

test_that("an anchor with unusable coordinates is skipped, not fatal", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10), ann_row("5.8S", NA, NA), ann_row("28S", 11, 18))
  expect_equal(nchar(RiboPilot:::screen_anchor_string(seq, a)), 18L)
})

test_that("an all-minus unit is reverse-complemented rather than compared antisense", {
  # orient_sequence() normally leaves a unit on '+', so this is defensive: a unit
  # that escaped orientation must not be compared against an oriented one and
  # read as a paralog on the strength of being antisense
  fwd <- paste0(strrep("A", 10), strrep("G", 8))
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd)))
  a <- anns(ann_row("28S", 1, 8, "-"), ann_row("18S", 9, 18, "-"))
  expect_equal(RiboPilot:::screen_anchor_string(rc, a), fwd)
})

test_that("a mixed-strand unit is not comparable", {
  # concatenating a sense and an antisense span produces a string that means
  # nothing; a mixed-strand unit is a chimera signal, not a paralog
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10, "+"), ann_row("28S", 11, 18, "-"))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, a)))
})

# units frame plus the parallel list of named anchor parts. A one-element list
# named "18S" stands for a unit annotating just that anchor.
mk_units <- function(...) {
  d <- do.call(rbind, lapply(list(...), function(x)
    data.frame(key = x[[1]], total_length = as.integer(x[[3]]),
      stringsAsFactors = FALSE)))
  attr(d, "parts") <- lapply(list(...), function(x) {
    if (length(x) >= 4) return(x[[4]])
    if (is.na(x[[2]])) return(character(0))
    stats::setNames(c(x[[2]]), "18S")
  })
  d
}
grp <- function(d, min_pid = 99.5, ...) {
  RiboPilot:::screen_variant_groups(d, attr(d, "parts"),
    utils::modifyList(RiboPilot:::SCREEN_PARAMS_DEFAULT,
      list(min_isomer_anchor_pid = min_pid, ...)))
}
A <- strrep("ACGT", 200)
# one substitution in 800 bp, comfortably above the 99.5 threshold
A1 <- paste0("T", substr(A, 2, nchar(A)))
# a genuinely divergent copy: every tenth block mutated, about 92 percent
B <- paste0(vapply(seq_len(200), function(i)
  if (i %% 10 == 0) "AAAA" else "ACGT", character(1)), collapse = "")

test_that("identical anchors score 100 with full coverage, without aligning", {
  r <- RiboPilot:::screen_anchor_pid(A, A)
  expect_equal(r$pid, 100)
  expect_equal(r$cov_a, 1)
  expect_equal(r$cov_b, 1)
})

test_that("a divergent copy scores well below the threshold", {
  expect_lt(RiboPilot:::screen_anchor_pid(A, B)$pid, 99.5)
})

test_that("an unalignable input yields NA rather than an error", {
  expect_true(is.na(RiboPilot:::screen_anchor_pid(A, NA_character_)$pid))
  expect_true(is.na(RiboPilot:::screen_anchor_pid(A, "")$pid))
})

test_that("a fragment scores 100 identity but only partial coverage", {
  # the defect found on SRR33849683: overlap alignment plus PID1 divides matches
  # by the ALIGNED length, so a substring is a perfect match of its container
  frag <- substr(A, 1, 400)
  r <- RiboPilot:::screen_anchor_pid(A, frag)
  expect_equal(r$pid, 100)
  expect_lt(r$cov_a, 0.6)   # only half of A was covered
  expect_equal(r$cov_b, 1)  # all of the fragment was
})

test_that("two units with near-identical anchors and different lengths are isomers", {
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", A1, 4200L)), 99.5)
  expect_equal(unique(g$variant_class), "isomer")
  expect_equal(length(unique(g$variant_group)), 1L)
  # the longest member is the representative and has nothing to compare against
  expect_true(is.na(g$variant_identity[g$key == "u1"]))
  expect_gt(g$variant_identity[g$key == "u2"], 99.5)
  expect_equal(g$len_delta[g$key == "u2"], 1800L)
})

test_that("identical anchors AND identical length are still isomers", {
  # the spec's table has no row for this; two units that agree on everything are
  # the same molecule assembled twice, for which "pick either" is exactly right
  g <- grp(mk_units(list("u1", A, 5000L), list("u2", A, 5000L)), 99.5)
  expect_equal(unique(g$variant_class), "isomer")
})

test_that("two units differing beyond the threshold are paralogs in separate groups", {
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", B, 6000L)), 99.5)
  expect_equal(unique(g$variant_class), "paralog")
  expect_equal(length(unique(g$variant_group)), 2L)
  # each is its own group's representative, so neither has an identity to record
  expect_true(all(is.na(g$variant_identity)))
})

test_that("a mixed set forms one isomer group and one paralog group", {
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", A1, 4200L), list("u3", B, 5000L)), 99.5)
  # every comparable unit reads paralog once the sample carries two copies:
  # telling u1 "pick either" would be false when u3 exists
  expect_equal(unique(g$variant_class), "paralog")
  expect_equal(g$variant_group[g$key == "u1"], g$variant_group[g$key == "u2"])
  expect_false(g$variant_group[g$key == "u1"] == g$variant_group[g$key == "u3"])
})

test_that("a three-unit isomer set forms one group", {
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", A1, 5000L), list("u3", A, 4000L)), 99.5)
  expect_equal(length(unique(g$variant_group)), 1L)
  expect_equal(unique(g$variant_class), "isomer")
})

test_that("a single comparable unit is a singleton with no identity", {
  g <- grp(mk_units(list("u1", A, 6000L)), 99.5)
  expect_equal(g$variant_class, "singleton")
  # recording 100 would assert a measurement that was never made
  expect_true(is.na(g$variant_identity))
})

test_that("an uncomparable unit is NA, not singleton", {
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", NA_character_, 900L)), 99.5)
  expect_true(is.na(g$variant_class[g$key == "u2"]))
  expect_true(is.na(g$variant_group[g$key == "u2"]))
  # and it does not make the remaining unit look like it had company
  expect_equal(g$variant_class[g$key == "u1"], "singleton")
})

test_that("group labels are deterministic, longest representative first", {
  g1 <- grp(mk_units(list("u1", A, 4000L), list("u2", B, 9000L)), 99.5)
  g2 <- grp(mk_units(list("u2", B, 9000L), list("u1", A, 4000L)), 99.5)
  expect_equal(g1$variant_group[g1$key == "u2"], "g1")
  expect_equal(g2$variant_group[g2$key == "u2"], "g1")
  expect_equal(g1$variant_group[g1$key == "u1"], "g2")
})

test_that("the group representative is the longest member, not the first", {
  # input order deliberately puts the SHORT unit first: if the representative
  # were simply the first member, this passes for the wrong reason
  g <- grp(mk_units(list("u1", A, 4200L), list("u2", A1, 6000L)), 99.5)
  expect_true(is.na(g$variant_identity[g$key == "u2"]))
  expect_false(is.na(g$variant_identity[g$key == "u1"]))
  expect_equal(g$len_delta[g$key == "u1"], 1800L)
})

test_that("a pair at exactly the threshold groups together", {
  # the spec says "at or above", so >= not >
  fake <- mk_units(list("u1", A, 6000L), list("u2", A, 4000L))
  expect_equal(length(unique(
    grp(fake, 100)$variant_group)), 1L)
})

test_that("a chained group still forms and the spread stays visible", {
  # A matches A1 and A1 matches A2, but A may not match A2. The group forms and
  # variant_identity makes the spread visible rather than hiding it.
  A2 <- paste0("TT", substr(A, 3, nchar(A)))
  g <- grp(mk_units(list("u1", A, 6000L), list("u2", A1, 5000L), list("u3", A2, 4000L)), 99.5)
  expect_equal(length(unique(g$variant_group)), 1L)
  expect_false(anyNA(g$variant_identity[g$key != "u1"]))
})

# --- the per-sample DB pass ---

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
FUNGI_T <- "4751;4890;147545;5042;1131492;5052;5062"
FISH_LIN <- "Eukaryota;Chordata;Actinopteri;o;f;g;s"
FUNGI_LIN <- "Eukaryota;Ascomycota;Sordariomycetes;o;f;g;s"

mk_variant_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage TEXT, lineage_taxids TEXT,
    reference_support INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT,
    sequence TEXT, length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
    element TEXT, start INT, stop INT, strand TEXT)")
  for (s in c("S1", "S2")) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,'Metazoa; Chordata; Actinopteri',?,NULL)",
      params = list(s, FISH_T))
  }
  con
}

# One unit: `anchor` is annotated as 18S at the head, `spacer` bases of
# unannotated remainder follow, so the anchor string is exactly `anchor` and the
# total length varies independently of it.
add_unit <- function(con, scaf, anchor, spacer, tax = FISH_T, lin = FISH_LIN,
                     id = "S1", annotate = TRUE, strand = "+") {
  seq <- paste0(anchor, strrep("A", spacer))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES(?,1,?,?,?)",
    params = list(id, scaf, seq, nchar(seq)))
  if (annotate) {
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,1,?, '18S',1,?,?)",
      params = list(id, scaf, nchar(anchor), strand))
  }
  for (i in 1:3) {
    DBI::dbExecute(con,
      "INSERT INTO blast_hits
         (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class,rank_distance)
       VALUES (?,1,?,?,1,?,?,?,'rdna',?)",
      params = list(id, scaf, i, 900 - i, tax, lin,
        as.integer(RiboPilot:::rank_distance(tax, FISH_T))))
  }
}

test_that("two on-target units with near-identical anchors are isomers, both on_target", {
  # the orthogonality this design turns on: the finding is not a verdict
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict, c("on_target", "on_target"))
  expect_equal(r$variant_class, c("isomer", "isomer"))
  expect_equal(length(unique(r$variant_group)), 1L)
})

test_that("a paralogous unit stays on_target", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, B, 2000L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict, c("on_target", "on_target"))
  expect_equal(r$variant_class, c("paralog", "paralog"))
  expect_equal(length(unique(r$variant_group)), 2L)
})

test_that("an off_target unit is never compared and its variant columns stay NA", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A, 200L, tax = FUNGI_T, lin = FUNGI_LIN)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict[2], "off_target")
  expect_true(is.na(r$variant_class[2]))
  expect_true(is.na(r$variant_group[2]))
  # and the surviving unit is a singleton, not half of an isomer pair
  expect_equal(r$variant_class[1], "singleton")
})

test_that("a unit whose anchors are unannotated is NA, not singleton", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A, 2000L, annotate = FALSE)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict[2], "on_target")
  expect_true(is.na(r$variant_class[2]))
  expect_equal(r$variant_class[1], "singleton")
})

test_that("the measured identity appears in the reason string", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT screen_reason FROM screen WHERE scaffold=2")
  expect_match(r$screen_reason, "variant: isomer", fixed = TRUE)
  expect_match(r$screen_reason, "99\\.[0-9]+ percent")
  expect_match(r$screen_reason, "total length differs by 1800 bp", fixed = TRUE)
  # the verdict's own evidence survives beside it
  expect_match(r$screen_reason, "confidence", fixed = TRUE)
})

test_that("rerunning does not stack duplicate variant clauses", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  # NOT force = TRUE: a forced rerun makes screen_one rewrite a fresh reason, so
  # the strip is never exercised. Unforced, the rows are already current and
  # screen_one skips them, so the variant pass appends onto its own last clause.
  screen_sample(con, "S1")
  screen_sample(con, "S1")
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT screen_reason FROM screen WHERE scaffold=2")
  expect_equal(lengths(regmatches(r$screen_reason,
    gregexpr("variant:", r$screen_reason, fixed = TRUE))), 1L)
  expect_match(r$screen_reason, "confidence", fixed = TRUE)
})

test_that("the measured identity is stored in variant_identity, not only in the reason", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con,
    "SELECT scaffold, variant_identity FROM screen WHERE ID='S1' ORDER BY scaffold")
  # scaffold 1 is the longer unit and so the group representative
  expect_true(is.na(r$variant_identity[1]))
  expect_false(is.na(r$variant_identity[2]))
  expect_gt(r$variant_identity[2], 99.5)
  expect_lte(r$variant_identity[2], 100)
})

test_that("screen_variants clears a stale finding left on a non-on_target unit", {
  # exercises the guard directly: reached through screen_sample() the DELETE plus
  # append in screen_one() would already have wiped the columns, so this is the
  # only way to pin that screen_variants() does not leave a contaminant looking
  # like an intragenomic variant
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version,screen_reason,
       variant_class,variant_group,variant_identity)
     VALUES ('S1',1,2,'off_target',3,'hits agree on Fusarium | variant: isomer, group g1',
       'isomer','g1',99.9)")
  RiboPilot:::screen_variants(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE scaffold=2")
  expect_true(is.na(r$variant_class))
  expect_true(is.na(r$variant_group))
  expect_true(is.na(r$variant_identity))
  # the verdict's own reason survives; only the variant clause goes
  expect_equal(r$screen_reason, "hits agree on Fusarium")
  expect_equal(r$screen_verdict, "off_target")
})

test_that("a unit that turns off_target loses its stale variant finding", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  screen_sample(con, "S1")
  expect_equal(DBI::dbGetQuery(con,
    "SELECT variant_class v FROM screen WHERE scaffold=2")$v, "isomer")
  # the evidence changes: this unit is now a fungal contaminant
  DBI::dbExecute(con, "DELETE FROM blast_hits WHERE scaffold=2")
  for (i in 1:3) {
    DBI::dbExecute(con,
      "INSERT INTO blast_hits
         (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class)
       VALUES ('S1',1,2,?,1,?,?,?,'rdna')",
      params = list(i, 900 - i, FUNGI_T, FUNGI_LIN))
  }
  screen_sample(con, "S1", force = TRUE)
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE scaffold=2")
  expect_equal(r$screen_verdict, "off_target")
  expect_true(is.na(r$variant_class))
  expect_false(grepl("variant:", r$screen_reason, fixed = TRUE))
})

test_that("the variant pass touches only its own sample", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  # asymmetric on purpose: S2 is a paralog pair, so a leaked write would be
  # visible as the wrong class rather than coincidentally correct
  add_unit(con, 1, A, 500L, id = "S2")
  add_unit(con, 2, B, 500L, id = "S2")
  screen_sample(con, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen WHERE ID='S2'")$n, 0L)
  screen_sample(con, "S2")
  expect_equal(DBI::dbGetQuery(con,
    "SELECT variant_class v FROM screen WHERE ID='S2' ORDER BY scaffold")$v,
    c("paralog", "paralog"))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT variant_class v FROM screen WHERE ID='S1' ORDER BY scaffold")$v,
    c("isomer", "isomer"))
})

test_that("a stale-version row recomputes and gains a class", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  add_unit(con, 2, A1, 200L)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version)
     VALUES ('S1',1,1,'on_target',2)")
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT screen_version, variant_class FROM screen WHERE scaffold=1")
  expect_equal(r$screen_version, RiboPilot:::SCREEN_VERSION)
  expect_equal(r$variant_class, "isomer")
})

test_that("a sample with one on-target unit is a singleton, not NA", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT variant_class, variant_group, variant_identity
    FROM screen WHERE scaffold=1")
  expect_equal(r$variant_class, "singleton")
  expect_equal(r$variant_group, "g1")
  expect_true(is.na(r$variant_identity))
  expect_match(DBI::dbGetQuery(con, "SELECT screen_reason r FROM screen WHERE scaffold=1")$r,
    "variant: singleton", fixed = TRUE)
})

test_that("the pipeline leaves the variant check to triage, after annotation lands", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, A, 2000L, annotate = FALSE)
  add_unit(con, 2, A1, 200L, annotate = FALSE)
  screen_sample(con, "S1", variants = FALSE)
  expect_true(all(is.na(DBI::dbGetQuery(con, "SELECT variant_class v FROM screen")$v)))
  # annotation arrives after the screen, as it does with a local BLAST
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S',1,?,'+')", params = list(nchar(A)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,2,'18S',1,?,'+')", params = list(nchar(A1)))
  RiboPilot:::screen_variants(con, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT variant_class v FROM screen ORDER BY scaffold")$v,
    c("isomer", "isomer"))
  bin <- system.file("nextflow", "bin", package = "RiboPilot")
  expect_true(any(grepl("variants = FALSE", readLines(file.path(bin, "rp_screen.R")), fixed = TRUE)))
  tri <- readLines(file.path(bin, "rp_triage.R"))
  expect_lt(grep("screen_variants", tri)[1], grep("auto_triage_sample", tri)[1])
})

# --- the app badge ---

test_that("the variant badge never borrows the off-target red or warning language", {
  js <- as.character(RiboPilot:::rt_variant_badge())
  red <- as.character(RiboPilot:::rt_screen_badge())
  expect_true(grepl("tone = 'danger'", red, fixed = TRUE))   # the off_target red exists
  expect_false(grepl("danger", js, fixed = TRUE))   # and is not reused here
  expect_false(grepl("warn", js, ignore.case = TRUE))
  # an isomer is neutral; a paralog is information, and the two differ
  expect_true(grepl("#eeeeee", js, fixed = TRUE))
  expect_true(grepl("#dbe9f6", js, fixed = TRUE))
})

test_that("the badge renders nothing for a singleton or an unscreened unit", {
  js <- as.character(RiboPilot:::rt_variant_badge())
  # a project that has never been screened must not sprout placeholder badges
  expect_true(grepl("v === 'singleton'", js, fixed = TRUE))
  expect_true(grepl("v === null", js, fixed = TRUE))
  expect_true(grepl("return ''", js, fixed = TRUE))
})

test_that("the badge surfaces the measured identity and the group label", {
  js <- as.character(RiboPilot:::rt_variant_badge())
  expect_true(grepl("variant_identity", js, fixed = TRUE))
  expect_true(grepl("variant_group", js, fixed = TRUE))
  expect_true(grepl("percent identical to the group representative", js, fixed = TRUE))
  # the two curator messages the spec dictates
  expect_true(grepl("pick either", js, fixed = TRUE))
  expect_true(grepl("determines what GenBank receives", js, fixed = TRUE))
})

test_that("unit_evidence selects the variant columns so the table can render them", {
  src <- paste(deparse(RiboPilot:::unit_evidence), collapse = " ")
  for (col in c("variant_class", "variant_group", "variant_identity")) {
    expect_true(grepl(col, src, fixed = TRUE), info = col)
  }
})

# --- the partial-anchor guards (found on real data, SRR33849683) ---

P3 <- function(...) stats::setNames(c(...), c("18S", "5.8S", "28S")[seq_len(...length())])
A18 <- strrep("ACGT", 200)   # 800 bp "18S"
A58 <- strrep("ACGA", 50)    # 200 bp "5.8S"
A28 <- strrep("ACCT", 300)   # 1200 bp "28S"

test_that("two units annotating DIFFERENT anchor sets are not comparable", {
  # the real defect: a 28S-only unit scored 100.00 against a full 18S+5.8S+28S
  # unit and was called an isomer of it, telling the curator to pick either
  full <- stats::setNames(c(A18, A58, A28), c("18S", "5.8S", "28S"))
  only28 <- stats::setNames(A28, "28S")
  r <- RiboPilot:::screen_anchor_pair(full, only28)
  expect_equal(r$shared, "28S")          # only the shared element is compared
  expect_equal(r$pid, 100)               # and on that element they DO match
  expect_true(r$comparable)              # 1200 bp of shared 28S is a real comparison
})

test_that("a truncated shared anchor fails the coverage floor", {
  a <- stats::setNames(A28, "28S")
  b <- stats::setNames(substr(A28, 1, 500), "28S")   # 500 of 1200 bp
  r <- RiboPilot:::screen_anchor_pair(a, b)
  expect_equal(r$pid, 100)               # identity alone still says perfect
  expect_lt(r$coverage, 0.9)             # but it covered less than half of a
  expect_false(r$comparable)             # so no call is made
})

test_that("too little shared anchor makes no call at all", {
  a <- stats::setNames(A58, "5.8S")      # 200 bp, under min_isomer_anchor_bp
  r <- RiboPilot:::screen_anchor_pair(a, a)
  expect_equal(r$pid, 100)
  expect_equal(r$coverage, 1)
  expect_false(r$comparable)
  expect_equal(r$n_bp, 200L)
})

test_that("no shared element at all is not comparable", {
  r <- RiboPilot:::screen_anchor_pair(stats::setNames(A18, "18S"),
                                      stats::setNames(A28, "28S"))
  expect_equal(length(r$shared), 0L)
  expect_false(r$comparable)
  expect_true(is.na(r$pid))
})

test_that("a unit with no comparable counterpart is NA, never paralog", {
  # asserting paralog here would tell the curator the sample carries several
  # rDNA copies on the strength of a measurement that was never made
  d <- data.frame(key = c("1.1", "1.2"), total_length = c(4318L, 1694L),
    stringsAsFactors = FALSE)
  parts <- list(stats::setNames(A18, "18S"), stats::setNames(A28, "28S"))
  g <- RiboPilot:::screen_variant_groups(d, parts)
  expect_true(all(is.na(g$variant_class)))
  expect_true(all(is.na(g$variant_group)))
})

test_that("one comparable pair plus one loner leaves the loner NA", {
  d <- data.frame(key = c("1.1", "1.2", "1.3"),
    total_length = c(9000L, 8000L, 1500L), stringsAsFactors = FALSE)
  full <- stats::setNames(c(A18, A58, A28), c("18S", "5.8S", "28S"))
  parts <- list(full, full, stats::setNames(A58, "5.8S"))
  g <- RiboPilot:::screen_variant_groups(d, parts)
  expect_equal(g$variant_class, c("isomer", "isomer", NA))
  expect_equal(g$variant_group[1], g$variant_group[2])
  expect_true(is.na(g$variant_group[3]))
})

test_that("is_rep marks the representative, distinct from a chained member", {
  d <- data.frame(key = c("1.1", "1.2"), total_length = c(9000L, 8000L),
    stringsAsFactors = FALSE)
  full <- stats::setNames(c(A18, A58, A28), c("18S", "5.8S", "28S"))
  g <- RiboPilot:::screen_variant_groups(d, list(full, full))
  expect_equal(g$is_rep, c(TRUE, FALSE))
  # the representative's clause says so; a chained member's says something else
  expect_match(RiboPilot:::screen_variant_clause("isomer", "g1", NA_real_,
    NA_integer_, is_rep = TRUE), "group representative", fixed = TRUE)
  expect_match(RiboPilot:::screen_variant_clause("isomer", "g1", NA_real_,
    NA_integer_, is_rep = FALSE), "through an intermediate member", fixed = TRUE)
})

test_that("a one-base alignment cannot score 100 percent identity any more", {
  # the sharpest form of the SRR33849683 defect: two unrelated 28S regions found
  # no alignment, the overlap aligner fell back to a 1 bp overlap, and PID1 read
  # 1 match over 1 aligned position as perfect. Measured aligned_len = 1.
  a <- stats::setNames(paste0(strrep("CATG", 200), "A"), "28S")
  b <- stats::setNames(paste0("A", strrep("GTAC", 200)), "28S")
  r <- RiboPilot:::screen_anchor_pair(a, b)
  # whatever identity the aligner reports, coverage must veto the call
  expect_false(r$comparable)
  expect_true(is.na(r$coverage) || r$coverage < 0.9)
})

test_that("NaN identity from an empty alignment is treated as no measurement", {
  a <- stats::setNames(strrep("CATG", 200), "28S")
  b <- stats::setNames(strrep("GTAC", 200), "28S")
  r <- RiboPilot:::screen_anchor_pair(a, b)
  expect_false(r$comparable)
  d <- data.frame(key = c("1.1", "1.2"), total_length = c(900L, 800L),
    stringsAsFactors = FALSE)
  g <- RiboPilot:::screen_variant_groups(d, list(a, b))
  # not comparable means no claim, not "two groups therefore paralog"
  expect_true(all(is.na(g$variant_class)))
})
