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

test_that("the packager, the installer and discovery demand the same panel files", {
  # three places used to disagree, so a panel could package and install cleanly
  # and still be invisible to discovery: remote fallback, no message
  skip_on_cran()
  p <- "../../dev/refdb/package_release.sh"
  skip_if_not(file.exists(p), "dev/ is Rbuildignored; this runs from the source tree only")
  txt <- readLines(p, warn = FALSE)
  i <- grep("the panel is incomplete", txt)
  expect_length(i, 1)
  j <- max(grep("^for e in ", txt[seq_len(i)]))
  exts <- strsplit(trimws(sub("^for e in (.*); do$", "\\1", txt[j])), " +")[[1]]
  expect_setequal(exts, RiboPilot:::PANEL_REQUIRED_EXT)

  empty <- withr::local_tempdir()
  miss <- setdiff(RiboPilot:::panel_missing_parts(empty), "manifest.json")
  expect_setequal(sub("^ribopilot_rrna\\.", "", miss), RiboPilot:::PANEL_REQUIRED_EXT)
})

test_that("the deploy scripts can build an image with no panel at all", {
  # the release the scripts stage from may not exist yet; that must not block
  # every image build, and a panel-less image must not advertise a panel
  skip_on_cran()
  dockerfile <- "../../docker/Dockerfile"
  skip_if_not(file.exists(dockerfile), "docker/ is Rbuildignored; this runs from the source tree only")
  df <- readLines(dockerfile, warn = FALSE)
  expect_true(any(grepl("^FROM db-\\$\\{RP_DB\\} AS final", df)))
  # RIBOPILOT_DB_DIR is set only in the stage that actually carries the panel
  db_env <- grep("^ENV RIBOPILOT_DB_DIR=", df)
  expect_length(db_env, 1)
  expect_gt(db_env, max(grep("^FROM base AS db-1", df)))

  for (s in c("../../docker/deploy-local.sh", "../../docker/deploy-dockerhub.sh")) {
    ds <- paste(readLines(s, warn = FALSE), collapse = "\n")
    expect_match(ds, "--no-db", fixed = TRUE, info = s)
    expect_match(ds, "RP_DB=${RP_DB}", fixed = TRUE, info = s)
    # gh, not curl: the repo is private, so an unauthenticated fetch 404s
    expect_match(ds, "gh release download", fixed = TRUE, info = s)
    expect_false(grepl("curl -fsSL", ds, fixed = TRUE), info = s)
  }
})

test_that("the deploy script checksum-verifies the db panel tarballs before building", {
  skip_on_cran()
  deploy <- "../../docker/deploy-dockerhub.sh"
  skip_if_not(file.exists(deploy), "docker/ is Rbuildignored; this runs from the source tree only")
  ds <- paste(readLines(deploy, warn = FALSE), collapse = "\n")
  # fetches the published checksum alongside the tarball
  expect_match(ds, ".sha256", fixed = TRUE)
  # and actually runs the verification, not just downloads it
  expect_match(ds, '-c "$sumfile"', fixed = TRUE)
})
