# Guard: the packaging script's asset names must agree with install_ribopilot_db().

test_that("the packaging script emits the asset names the installer requests", {
  skip_on_cran()
  p <- "../../dev/refdb/package_release.sh"
  skip_if_not(file.exists(p), "dev/ is Rbuildignored; this runs from the source tree only")
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_match(txt, 'ribopilot_rrna-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'decoys-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'db-${VER}', fixed = TRUE)
  # the other half of the contract: the names the installer actually asks for
  fn <- paste(deparse(body(install_ribopilot_db)), collapse = "\n")
  expect_match(fn, 'ribopilot_rrna-%s.tar.gz', fixed = TRUE)
  expect_match(fn, 'decoys-%s.tar.gz', fixed = TRUE)
  expect_match(fn, '"db-"', fixed = TRUE)
})
