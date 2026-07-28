# The modal is the only place a curator sees what triage will do, so its
# contents are pinned: one control per rule, drop offered only where the
# evidence supports it, and a legend that says what the two counts mean.

test_that("the modal renders one control per rule, in registry order", {
  m <- RiboPilot:::triage_modal(3, RiboPilot:::triage_actions_default())
  html <- as.character(m)
  for (k in RiboPilot:::triage_rule_keys()) {
    expect_match(html, paste0("triage_rule_", k), fixed = TRUE)
  }
  pos <- vapply(RiboPilot:::triage_rule_keys(),
    function(k) regexpr(paste0("triage_rule_", k), html, fixed = TRUE)[1], numeric(1))
  expect_false(is.unsorted(pos))
})

test_that("drop is offered only on the combination rules", {
  m <- RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default())
  html <- as.character(m)
  for (def in RiboPilot:::TRIAGE_RULE_DEFS) {
    blk <- RiboPilot:::triage_rule_control(def, def$default)
    has_drop <- grepl("Drop", as.character(blk), fixed = TRUE)
    expect_equal(has_drop, "drop" %in% def$actions)
  }
  expect_true(nzchar(html))
})

test_that("the legend explains both counts without jargon", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_match(html, "first one that matches decides whether to drop or flag", fixed = TRUE)
  expect_match(html, "even if a higher-priority rule acts on it first", fixed = TRUE)
})

test_that("the modal no longer offers a depth control", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_false(grepl("triage_depth_frac", html, fixed = TRUE))
  expect_false(grepl("deepest scaffold", html, fixed = TRUE))
})
