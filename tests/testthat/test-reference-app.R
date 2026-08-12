# The Reference mode is reachable from the UI and its rows carry the state icon.

test_that("the mode picker offers Reference between Taxonomy and Assemble", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  expect_true(grepl("Reference", html, fixed = TRUE))
  expect_lt(regexpr("Reference", html, fixed = TRUE),
    regexpr("Annotation Review", html, fixed = TRUE))
})

test_that("the Reference control bar keeps only the two tab-wide actions", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  for (id in c("ref_suggest_btn", "ref_dir_btn")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  # the per-group verbs moved into the row; a bar button for them would bring
  # back the guess-the-selection-rule problem they were moved to fix
  for (id in c("ref_edit_btn", "ref_existing_btn", "ref_skip_btn", "ref_build_btn")) {
    expect_false(grepl(id, html, fixed = TRUE), info = id)
  }
})

test_that("the empty-table message can wrap without stacking on itself", {
  # reactable's own stylesheet sets `.rt-no-data { line-height: 0 }`, which only
  # works for its one-word default. Every mode's message here is a sentence, and
  # the Reference one runs to three lines, so without this override each wrapped
  # line is drawn on top of the last.
  css <- paste(readLines(RiboPilot:::app_sys("app/www/custom.css"), warn = FALSE), collapse = "\n")
  # drop comments first: the one above this rule quotes the very declaration
  # being overridden, and would otherwise match instead of the rule
  css <- gsub("(?s)/\\*.*?\\*/", "", css, perl = TRUE)
  rule <- regmatches(css, regexpr("\\.rt-no-data\\s*\\{[^}]*\\}", css))
  expect_length(rule, 1L)
  expect_match(rule, "line-height:\\s*[1-9]")
})

test_that("the row action modal names what skipping actually costs", {
  html <- as.character(RiboPilot:::reference_action_modal("asteroidea", "Asteroidea"))
  for (id in c("ref_act_build", "ref_act_existing", "ref_act_clade", "ref_act_skip")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  expect_true(grepl("Asteroidea", html, fixed = TRUE))
  expect_true(grepl("fish", html, fixed = TRUE))
})

test_that("the clade option is dropped for a multi-group selection", {
  html <- as.character(RiboPilot:::reference_action_modal(c("a", "b"), "2 selected groups"))
  expect_false(grepl("ref_act_clade", html, fixed = TRUE))
  expect_true(grepl("ref_act_build", html, fixed = TRUE))
})
