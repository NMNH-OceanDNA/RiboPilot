# The Reference mode is reachable from the UI and its rows carry the state icon.

test_that("the mode picker offers Reference between Taxonomy and Assemble", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  expect_true(grepl("Reference", html, fixed = TRUE))
  expect_lt(regexpr("Reference", html, fixed = TRUE),
    regexpr("Annotation Review", html, fixed = TRUE))
})

test_that("the Reference control bar exposes every action", {
  html <- as.character(RiboPilot:::app_ui(NULL))
  for (id in c("ref_suggest_btn", "ref_edit_btn", "ref_existing_btn", "ref_skip_btn",
               "ref_build_btn", "ref_dir_btn")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
})

test_that("skipped rows get their own state icon", {
  expect_true("5" %in% names(RiboPilot:::REFERENCE_STATE_ICONS))
})
