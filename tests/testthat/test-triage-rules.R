# The rule registry is the whole triage decision. Each rule is tested in
# isolation, and the 2x2 over the two evidence channels is tested exhaustively,
# because a gap there is a scaffold that silently keeps or silently drops.

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
COPE_T <- "33208;6656;6854;6843;1;2;3"

# One unit row shaped like unit_triage_state() output.
u <- function(...) {
  r <- list(...)
  data.frame(
    path = 1L, scaffold = r$scaffold %||% 1L, length = 5000L, ignore = 0L,
    triage = NA_character_, triage_reason = NA_character_,
    triage_source = NA_character_, triage_rule = NA_character_,
    depth = 100, screen_species = r$sp %||% NA_character_,
    screen_taxids = r$tx %||% NA_character_, rank_distance = NA_integer_,
    is_nuclear_rrna = r$nuc %||% NA_integer_, sample_taxids = FISH_T,
    sci_name = "Danio rerio", n_features = r$feats %||% 5L,
    annotated_n = r$ann %||% 5L, annotate_lock = 0L, reviewed = "no",
    stringsAsFactors = FALSE)
}
us <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(seq_along(rows), function(i) {
    do.call(u, c(rows[[i]], list(scaffold = i)))
  }))
}

test_that("the rule keys and their order are fixed", {
  expect_equal(RiboPilot:::triage_rule_keys(),
    c("not_annotated", "no_screen", "no_feat_non_rdna",
      "no_feat_rdna_hit", "feat_non_rdna_hit", "odd_one_out"))
})

test_that("annotation that never ran is flagged, never dropped", {
  # ann = NA means annotate.n_features is NULL: reconcile never reached this
  # unit. annotate_tools carries errorStrategy 'ignore', so this is a tool
  # failure, and the old rule dropped it as if it were evidence.
  d <- RiboPilot:::triage_rules(us(list(ann = NA, feats = 0, nuc = 0, sp = "Homo sapiens")))
  expect_equal(d$new_rule, "not_annotated")
  expect_true(is.na(d$new_triage))
  expect_match(d$new_reason, "annotation has not run")
})

test_that("a scaffold with no usable hit is flagged, never dropped", {
  d <- RiboPilot:::triage_rules(us(list(feats = 5, nuc = NA)))
  expect_equal(d$new_rule, "no_screen")
  expect_true(is.na(d$new_triage))
  expect_match(d$new_reason, "GenBank screen has not run")
})

test_that("the 2x2 over the two channels is exhaustive", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T)))
  expect_equal(d$new_rule[1], "no_feat_non_rdna")
  expect_equal(d$new_triage[1], "drop")
  expect_equal(d$new_rule[2], "no_feat_rdna_hit")
  expect_true(is.na(d$new_triage[2]))
  expect_equal(d$new_rule[3], "feat_non_rdna_hit")
  expect_true(is.na(d$new_triage[3]))
  expect_true(is.na(d$new_rule[4]))
  expect_equal(d$new_triage[4], "keep")
})

test_that("taxonomic distance alone never drops a scaffold", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Cyclops sp.", tx = COPE_T)))
  expect_equal(d$new_rule[4], "odd_one_out")
  expect_false(identical(d$new_triage[4], "drop"))
})

test_that("depth never appears in any reason", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T)))
  # %|NA|% is scalar only (R/utils.R:10), so drop the NAs instead of coalescing
  expect_false(any(grepl("depth", stats::na.omit(d$new_reason), ignore.case = TRUE)))
})

test_that("turning a rule off moves its units to the next MATCHING rule", {
  # the odd unit matches both no_feat_non_rdna and odd_one_out, so turning the
  # first off must fall through to the second rather than to keep
  frame <- function() us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 0, nuc = 0, sp = "Cyclops sp.", tx = COPE_T))
  d0 <- RiboPilot:::triage_rules(frame())
  expect_equal(d0$new_rule[4], "no_feat_non_rdna")
  expect_equal(d0$new_triage[4], "drop")

  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "off"))
  d1 <- RiboPilot:::triage_rules(frame(), act)
  expect_equal(d1$new_rule[4], "odd_one_out")
  expect_true(is.na(d1$new_triage[4]))

  # when NO other rule matches, turning the only matching rule off IS a keep
  d2 <- RiboPilot:::triage_rules(
    us(list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T)), act)
  expect_true(is.na(d2$new_rule))
  expect_equal(d2$new_triage, "keep")
})

test_that("an action a rule cannot support is softened to flag", {
  act <- RiboPilot:::normalize_triage_actions(
    list(odd_one_out = "drop", no_screen = "drop", not_annotated = "drop"))
  expect_equal(act$odd_one_out, "flag")
  expect_equal(act$no_screen, "flag")
  expect_equal(act$not_annotated, "flag")

  # off is always legal, and drop is honoured on the combination rules
  act2 <- RiboPilot:::normalize_triage_actions(
    list(odd_one_out = "off", no_feat_rdna_hit = "drop"))
  expect_equal(act2$odd_one_out, "off")
  expect_equal(act2$no_feat_rdna_hit, "drop")
})

test_that("setting every combination rule off or flag yields no drops", {
  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "flag"))
  d <- RiboPilot:::triage_rules(us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T)), act)
  expect_false(any(d$new_triage %in% "drop"))
})

test_that("counts partition on Decided and overlap on Matching", {
  d <- us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T))
  cnt <- RiboPilot:::triage_rule_counts(d)
  expect_equal(nrow(cnt), 6L)
  expect_equal(sum(cnt$n_decided), 3L)          # the clean unit fires nothing
  expect_true(all(cnt$n_matching >= cnt$n_decided))

  # a rule set to off decides nothing but still matches
  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "off"))
  cnt2 <- RiboPilot:::triage_rule_counts(d, act)
  off <- cnt2[cnt2$key == "no_feat_non_rdna", ]
  expect_equal(off$n_decided, 0L)
  expect_equal(off$n_matching, 1L)
})

test_that("an empty frame returns empty columns rather than erroring", {
  d <- RiboPilot:::triage_rules(us()[0, ])
  expect_equal(nrow(d), 0L)
  expect_true(all(c("new_triage", "new_reason", "new_rule") %in% names(d)))
})
