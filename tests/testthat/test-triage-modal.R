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
  # asserts against the RENDERED modal, not against triage_rule_control()'s own
  # inputs: comparing a helper's output to the argument it was given passes
  # whatever the modal actually shows
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  keys <- RiboPilot:::triage_rule_keys()
  pos <- vapply(keys, function(k)
    regexpr(paste0('id="triage_rule_', k, '"'), html, fixed = TRUE)[1], numeric(1))
  expect_false(any(pos < 0))
  ends <- c(pos[-1] - 1, nchar(html))
  for (i in seq_along(keys)) {
    seg <- substr(html, pos[i], ends[i])
    expect_true(grepl('value="off"', seg, fixed = TRUE), info = keys[i])
    expect_true(grepl('value="flag"', seg, fixed = TRUE), info = keys[i])
    expect_equal(grepl('value="drop"', seg, fixed = TRUE),
      "drop" %in% RiboPilot:::TRIAGE_RULE_DEFS[[i]]$actions, info = keys[i])
  }
})

test_that("the legend explains both counts without jargon", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_match(html, "first one that matches decides whether to drop or flag", fixed = TRUE)
  expect_match(html, "even if a higher-priority rule acts on it first", fixed = TRUE)
})

test_that("the modal no longer offers a depth control", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_false(grepl("triage_depth_frac", html, fixed = TRUE))
  expect_false(grepl("deepest contig", html, fixed = TRUE))
})
