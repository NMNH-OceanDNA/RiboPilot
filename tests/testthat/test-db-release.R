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

test_that("the Dockerfile and deploy script agree with the installer on asset names and RIBOPILOT_DB_DIR", {
  skip_on_cran()
  dockerfile <- "../../docker/Dockerfile"
  deploy <- "../../docker/deploy-dockerhub.sh"
  skip_if_not(file.exists(dockerfile), "docker/ is Rbuildignored; this runs from the source tree only")
  skip_if_not(file.exists(deploy), "docker/ is Rbuildignored; this runs from the source tree only")

  df <- paste(readLines(dockerfile, warn = FALSE), collapse = "\n")
  expect_match(df, "ADD docker/ribopilot_rrna-${RP_DB_VERSION}.tar.gz", fixed = TRUE)
  expect_match(df, "ADD docker/decoys-${RP_DB_VERSION}.tar.gz", fixed = TRUE)
  expect_match(df, "ENV RIBOPILOT_DB_DIR=/opt/ribopilot/blastdb", fixed = TRUE)

  ds <- paste(readLines(deploy, warn = FALSE), collapse = "\n")
  expect_match(ds, "ribopilot_rrna-${RP_DB_VERSION}.tar.gz", fixed = TRUE)
  expect_match(ds, "decoys-${RP_DB_VERSION}.tar.gz", fixed = TRUE)
  expect_match(ds, "RIBOPILOT_DB_DIR", fixed = TRUE)
})
