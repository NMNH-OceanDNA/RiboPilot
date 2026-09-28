test_that("mp_js_attr escapes quotes and angle brackets", {
  expect_equal(mp_js_attr("it's <b>"), "it&#39;s &lt;b&gt;")
})

test_that("rt_dynamicIcon only names the glyph when labels are supplied", {
  bare <- as.character(rt_dynamicIcon(c(`0` = "fa fa-lock-open")))
  expect_true(grepl("var labels = {};", bare, fixed = TRUE))
  expect_false(grepl("reactable-btth", bare, fixed = TRUE))

  named <- as.character(rt_dynamicIcon(
    c(`0` = "fa fa-lock-open"),
    labels = c(`0` = "Unlocked")
  ))
  expect_true(grepl("role='img'", named, fixed = TRUE))
  expect_true(grepl("Unlocked", named, fixed = TRUE))
})

test_that("rt_dynamicIcon with an inputId sends its row to that input", {
  expect_true(grepl("var clickid = '';", as.character(rt_lock_icon()), fixed = TRUE))
  js <- as.character(rt_lock_icon("lock_row"))
  expect_true(grepl("var clickid = 'lock_row';", js, fixed = TRUE))
  expect_true(grepl("cellInfo.index + 1", js, fixed = TRUE))
})

test_that("a skippable confirm stops asking once its box is ticked", {
  shown <- 0L
  local_mocked_bindings(confirmSweetAlert = function(...) shown <<- shown + 1L,
    .package = "shinyWidgets")
  s <- shiny::MockShinySession$new()
  expect_false(mp_confirm("x", "T", "body", "Go", skippable = TRUE, session = s))
  s$setInputs(x_skip = TRUE)
  expect_true(mp_confirm("x", "T", "body", "Go", skippable = TRUE, session = s))
  expect_false(mp_confirm("x", "T", "body", "Go", session = s))
  expect_equal(shown, 2L)
})

test_that("rt_icon_bttn_text keeps the cell value as the render gate", {
  js <- as.character(rt_icon_bttn_text("id", "fa fa-eye", "view",
                                       label = "Coverage", title = "Show it"))
  # The empty-value early return must survive relabelling, or buttons appear
  # on rows with nothing behind them.
  expect_true(grepl("value==='') {", js, fixed = TRUE))
  expect_lt(regexpr("return;", js, fixed = TRUE),
            regexpr("value = 'Coverage';", js, fixed = TRUE))
  expect_true(grepl("title='Show it'", js, fixed = TRUE))
  expect_true(grepl("type='button'", js, fixed = TRUE))
})

test_that("rt_bool_badge uses pill classes and never a green no", {
  js <- as.character(rt_bool_badge())
  expect_true(grepl("mp-pill", js, fixed = TRUE))
  expect_true(grepl("'yes': 'success'", js, fixed = TRUE))
  expect_true(grepl("'no': 'neutral'", js, fixed = TRUE))
  expect_true(grepl("not set", js, fixed = TRUE))

  hidden <- as.character(rt_bool_badge(hide_no = TRUE))
  expect_true(grepl("var hide = ['no']", hidden, fixed = TRUE))
})

test_that("rt_header returns plain text without a tip", {
  expect_identical(rt_header("Gene Order"), "Gene Order")
  expect_true(grepl(
    "mp-th-tip",
    as.character(rt_header("Gene Order", "How the genes are ordered")),
    fixed = TRUE
  ))
})

test_that("rt_longtext escapes the tooltip and skips empty cells", {
  js <- as.character(rt_longtext())
  expect_true(grepl("mp-abbr", js, fixed = TRUE))
  expect_true(grepl("if (raw === '') return ''", js, fixed = TRUE))
})

test_that("rt_link renders locked cells as text", {
  js <- as.character(rt_link("assemble-pre_opts", lock_col = "assemble_lock"))
  expect_match(js, "mp-locked-cell")
  expect_match(js, "\\['assemble_lock'\\]")
  expect_false(grepl("mp-locked-cell", as.character(rt_link("x"))))
})

test_that("mp_col_header names columns and carries tips from the registry", {
  expect_match(as.character(RiboPilot:::mp_col_header("assemble_switch")), "mp-th-tip")
  expect_identical(RiboPilot:::mp_col_header("ID"), "ID")
  expect_equal(names(RiboPilot:::mp_state_icon_classes("assemble")), c("0", "1", "4", "2", "3"))
})
