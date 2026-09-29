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

test_that("the empty-table message stays readable and clear of the header", {
  # Two reactable defaults assume its own one-word default. `line-height: 0`
  # draws each wrapped line on top of the last, and absolute centring on
  # `.rt-tbody-no-data` puts the block on the header row, because an empty body
  # has no height. Both overrides have to survive for the message to be legible.
  css <- paste(readLines(RiboPilot:::app_sys("app/www/custom.css"), warn = FALSE), collapse = "\n")
  # drop comments first: the one above this rule quotes the very declarations
  # being overridden, and would otherwise match instead of the rule
  css <- gsub("(?s)/\\*.*?\\*/", "", css, perl = TRUE)
  rule <- regmatches(css, regexpr("\\.rt-no-data\\s*\\{[^}]*\\}", css))
  expect_length(rule, 1L)
  expect_match(rule, "line-height:\\s*[1-9]")
  expect_match(rule, "position:\\s*static")
  # tables with rows keep an empty message node, which must take no space
  expect_match(css, "\\.rt-tbody-no-data \\.rt-no-data\\s*\\{")
})

test_that("the settings modal groups the folder and the record window", {
  html <- as.character(RiboPilot:::reference_settings_modal(
    list(ref_dir = "/lab/refs", min_records = 1000L, max_records = 25000L)))
  expect_true(grepl("Reference settings", html, fixed = TRUE))
  expect_true(grepl("Reference database folder", html, fixed = TRUE))
  expect_true(grepl("/lab/refs", html, fixed = TRUE))
  for (id in c("ref_dir_path", "ref_min", "ref_max", "ref_dir_apply")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  # the measure is named once, in the heading, so the two fields can be one word
  # each instead of repeating "GenBank records for a suggested clade" twice
  expect_true(grepl("GenBank records per suggested clade", html, fixed = TRUE))
  expect_false(grepl("Minimum GenBank records", html, fixed = TRUE))
  expect_false(grepl("Maximum GenBank records", html, fixed = TRUE))
})

test_that("the settings modal exposes a panel path input carrying the current value", {
  # distinct from the field's placeholder text (R/app_reference.R), so this can
  # only pass via the value= wiring, not the placeholder rendering unconditionally
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000,
    blastdb = "/mnt/shared/panels/myclade_rrna")
  html <- as.character(RiboPilot:::reference_settings_modal(cfg))
  expect_match(html, "ref_blastdb_path")
  expect_match(html, "/mnt/shared/panels/myclade_rrna", fixed = TRUE)
})

test_that("the settings modal renders when no panel is configured", {
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000, blastdb = "")
  expect_match(as.character(RiboPilot:::reference_settings_modal(cfg)), "ref_blastdb_path")
})

test_that("the settings modal names the panel a run last searched against", {
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000, blastdb = "")
  html <- as.character(RiboPilot:::reference_settings_modal(cfg, "ribopilot_rrna 2026.08.1"))
  expect_match(html, "ribopilot_rrna 2026.08.1", fixed = TRUE)
})

test_that("the settings modal still renders, with no provenance line, when nothing was recorded", {
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000, blastdb = "")
  html <- as.character(RiboPilot:::reference_settings_modal(cfg))
  expect_match(html, "ref_blastdb_path")
  expect_false(grepl("Last run searched against", html, fixed = TRUE))
})

test_that("the row action modal names what skipping actually costs", {
  html <- as.character(RiboPilot:::reference_action_modal("asteroidea", "Asteroidea"))
  for (id in c("ref_act_build", "ref_act_existing", "ref_act_clade", "ref_act_skip")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
  expect_true(grepl("Asteroidea", html, fixed = TRUE))
  expect_true(grepl("their current database", html, fixed = TRUE))
  # a stack of short buttons opts out of the app's 600px modal floor
  expect_true(grepl("rp-modal-narrow", html, fixed = TRUE))
  # the hint sits under a bare button, which has no bottom margin for
  # opts_help()'s default negative top margin to absorb
  expect_false(grepl("margin-top: -", html, fixed = TRUE))
})

test_that("the clade option is dropped for a multi-group selection", {
  html <- as.character(RiboPilot:::reference_action_modal(c("a", "b"), "2 selected groups"))
  expect_false(grepl("ref_act_clade", html, fixed = TRUE))
  expect_true(grepl("ref_act_build", html, fixed = TRUE))
})
