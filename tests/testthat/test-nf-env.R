# How the pipeline's R processes find RiboPilot.
#
# rp_env.R exists so a stale installed copy cannot silently shadow the code
# being developed: both are version 0.0.0.9000, so nothing else detects the
# mismatch. It never fired. nf_command() passed `<root>/inst`, which has no
# DESCRIPTION, so the source branch was skipped and every pipeline R process
# quietly used the installed package instead.
#
# Under R CMD check the tests run against the INSTALLED package, so there is no
# source checkout to find and the source-loading path cannot be exercised at
# all. Those tests skip rather than assert something only true in a dev session.

is_source_checkout <- function(d) {
  file.exists(file.path(d, "DESCRIPTION")) &&
    length(list.files(file.path(d, "R"), pattern = "\\.[Rr]$")) > 0
}

test_that("the pipeline is pointed at a package source root, not inst/", {
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!is_source_checkout(src), "installed package, not a source checkout")

  dir <- withr::local_tempdir()
  file.create(file.path(dir, ".sqlite"))
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(dir, ".config"))
  args <- RiboPilot:::nf_command(dir, "assemble")
  passed <- args[which(args == "--ribopilot_src") + 1]
  expect_length(passed, 1)
  # what rp_env.R checks before it will load from source at all
  expect_true(file.exists(file.path(passed, "DESCRIPTION")))
  expect_gt(length(list.files(file.path(passed, "R"), pattern = "\\.[Rr]$")), 0)
  expect_equal(args[which(args == "--ribopilot_version") + 1],
    as.character(utils::packageVersion("RiboPilot")))
})

test_that("an installed package passes no source folder", {
  skip_if_not_installed("withr")
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "lib", "RiboPilot"), recursive = TRUE)
  file.create(file.path(tmp, "lib", "RiboPilot", "DESCRIPTION"))
  testthat::local_mocked_bindings(app_sys = function(...) file.path(tmp, "lib", "RiboPilot"),
    .package = "RiboPilot")
  expect_equal(RiboPilot:::ribopilot_src_default(), "")

  dir <- withr::local_tempdir()
  file.create(file.path(dir, ".sqlite"))
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(dir, ".config"))
  args <- RiboPilot:::nf_command(dir)
  expect_false("--ribopilot_src" %in% args)
})

test_that("an installed library directory is not mistaken for a source checkout", {
  skip_on_cran()
  env_r <- system.file("nextflow", "bin", "rp_env.R", package = "RiboPilot")
  skip_if(!nzchar(env_r), "rp_env.R not found")
  # an installed package has a DESCRIPTION and an R/ directory too, but R/ holds
  # the lazy-load database, not sources. Loading it as source leaves an empty
  # namespace, so every pipeline process dies on "could not find function".
  lib <- withr::local_tempdir()
  file.create(file.path(lib, "DESCRIPTION"))
  dir.create(file.path(lib, "R"))
  file.create(file.path(lib, "R", "RiboPilot.rdb"))

  # evaluate only the definitions at the top of rp_env.R, not the load it does
  env <- new.env()
  for (e in as.list(parse(env_r))) {
    eval(e, env)
    if (exists("rp_is_source_dir", envir = env, inherits = FALSE)) break
  }
  expect_false(env$rp_is_source_dir(lib))

  # and it still recognises a real checkout, when this session has one
  src <- RiboPilot:::ribopilot_src_default()
  if (is_source_checkout(src)) expect_true(env$rp_is_source_dir(src))
})

test_that("rp_env.R loads from source when handed either the root or inst/", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  env_r <- system.file("nextflow", "bin", "rp_env.R", package = "RiboPilot")
  skip_if(!nzchar(env_r), "rp_env.R not found")
  root <- RiboPilot:::ribopilot_src_default()
  skip_if(!is_source_checkout(root), "installed package, not a source checkout")

  for (src in c(root, file.path(root, "inst"))) {
    out <- suppressWarnings(system2("Rscript", c("-e",
      shQuote(sprintf("source('%s'); cat(pkgload::is_dev_package('RiboPilot'))", env_r))),
      env = paste0("RIBOPILOT_SRC=", src), stdout = TRUE, stderr = FALSE))
    expect_true(any(grepl("TRUE", out)), info = src)
  }
})

test_that("rp_env.R says so when it falls back to the installed package", {
  skip_on_cran()
  env_r <- system.file("nextflow", "bin", "rp_env.R", package = "RiboPilot")
  skip_if(!nzchar(env_r), "rp_env.R not found")
  skip_if_not_installed("RiboPilot")

  # the container case: RIBOPILOT_SRC points somewhere this process cannot see.
  # The message has to survive rp_env.R's own suppressMessages(), or "my edit
  # had no effect" looks exactly like a bug in the edit.
  out <- suppressWarnings(system2("Rscript", c("-e", shQuote(sprintf("source('%s')", env_r))),
    env = "RIBOPILOT_SRC=/nonexistent/ribopilot", stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("not visible from this process", out)))
  expect_true(any(grepl("installed RiboPilot", out)))
})
