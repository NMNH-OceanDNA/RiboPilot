# Guard: the packaging script's asset names must agree with install_ribopilot_db().

test_that("the packaging script emits the asset names the installer requests", {
  skip_on_cran()
  s <- readLines("../../dev/refdb/package_release.sh", warn = FALSE)
  txt <- paste(s, collapse = "\n")
  expect_match(txt, 'ribopilot_rrna-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'decoys-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'db-${VER}', fixed = TRUE)
  # the installer builds exactly these
  expect_equal(sprintf("ribopilot_rrna-%s.tar.gz", "1.2.3"), "ribopilot_rrna-1.2.3.tar.gz")
  expect_equal(paste0("db-", "1.2.3"), "db-1.2.3")
})
