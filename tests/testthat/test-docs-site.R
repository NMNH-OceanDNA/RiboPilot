test_that("site pages are ASCII and follow the app route", {
  root <- testthat::test_path("../..")
  vig <- file.path(root, "vignettes")
  if (!dir.exists(vig)) skip("package sources not available")
  pages <- c(list.files(vig, pattern = "[.]Rmd$", full.names = TRUE),
    file.path(root, "README.Rmd"))
  for (f in pages) {
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    expect_false(any(grepl("[^\\x01-\\x7F]", txt, perl = TRUE)), info = basename(f))
    expect_false(any(grepl("^\\s*RiboPilot\\(", txt)), info = basename(f))
    expect_false(any(grepl("PASTE_YOUR_KEY_HERE", txt, fixed = TRUE)), info = basename(f))
  }
  gs <- readLines(file.path(vig, "RiboPilot.Rmd"), warn = FALSE)
  expect_false(any(grepl("run_assembly\\(|run_fallback\\(|export_project\\(|conda_env", gs)))
})

test_that("the reference index names every Rd topic and only those", {
  root <- testthat::test_path("../..")
  skip_if_not(file.exists(file.path(root, "_pkgdown.yml")), "package sources not available")
  skip_if_not_installed("yaml")
  listed <- unlist(lapply(yaml::read_yaml(file.path(root, "_pkgdown.yml"))$reference, `[[`, "contents"))
  rd <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  aliases <- unlist(lapply(rd, function(f) sub("^\\\\alias\\{(.*)\\}$", "\\1",
    grep("^\\\\alias\\{", readLines(f, warn = FALSE), value = TRUE))))
  expect_equal(setdiff(listed, aliases), character(0))
  expect_equal(setdiff(sub("[.]Rd$", "", basename(rd)), listed), character(0))
})
