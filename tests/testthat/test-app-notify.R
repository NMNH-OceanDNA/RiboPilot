# Center alert vs bottom-right toast: the mapping from the app's message types
# to sweetalert types and titles.

capture_alert <- function(expr) {
  got <- NULL
  testthat::local_mocked_bindings(
    sendSweetAlert = function(...) {
      got <<- list(...)
      invisible(NULL)
    },
    .package = "shinyWidgets")
  force(expr)
  got
}

test_that("message maps to a success alert", {
  a <- capture_alert(rp_alert("Updated 3 sample(s).", session = NULL))
  expect_equal(a$type, "success")
  expect_equal(a$title, "Done")
  expect_equal(a$text, "Updated 3 sample(s).")
  expect_false(a$html)
})

test_that("warning and error keep their own type and title", {
  w <- capture_alert(rp_alert("No samples selected.", "warning", session = NULL))
  expect_equal(w$type, "warning")
  expect_equal(w$title, "Heads up")
  e <- capture_alert(rp_alert("Export failed.", "error", session = NULL))
  expect_equal(e$type, "error")
  expect_equal(e$title, "Error")
})

test_that("an unknown type falls back rather than erroring", {
  a <- capture_alert(rp_alert("hi", "not-a-type", session = NULL))
  expect_equal(a$type, "success")
})

test_that("an explicit title wins", {
  a <- capture_alert(rp_alert("hi", "error", title = "Header cannot be used", session = NULL))
  expect_equal(a$title, "Header cannot be used")
})

test_that("a tag body is sent as html, so one alert can carry several lines", {
  a <- capture_alert(rp_alert(shiny::tagList(shiny::tags$p("one"), shiny::tags$p("two")),
    "warning", session = NULL))
  expect_true(a$html)
})
