# How the pipeline's R processes find RiboPilot.
#
# rp_env.R exists so a stale installed copy cannot silently shadow the code
# being developed: both are version 0.0.0.9000, so nothing else detects the
# mismatch. It never fired. nf_command() passed `<root>/inst`, which has no
# DESCRIPTION, so the source branch was skipped and every pipeline R process
# quietly used the installed package instead.

test_that("the pipeline is pointed at a package source root, not inst/", {
  src <- RiboPilot:::ribopilot_src_default()
  expect_true(file.exists(file.path(src, "DESCRIPTION")))

  dir <- withr::local_tempdir()
  file.create(file.path(dir, ".sqlite"))
  args <- RiboPilot:::nf_command(dir, "assemble")
  passed <- args[which(args == "--ribopilot_src") + 1]
  expect_length(passed, 1)
  # what rp_env.R checks before it will load from source at all
  expect_true(file.exists(file.path(passed, "DESCRIPTION")))
  expect_true(dir.exists(file.path(passed, "R")))
})

test_that("rp_env.R loads from source when handed either the root or inst/", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  env_r <- system.file("nextflow", "bin", "rp_env.R", package = "RiboPilot")
  skip_if(!nzchar(env_r), "rp_env.R not found")
  root <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(root, "DESCRIPTION")), "not a source checkout")

  for (src in c(root, file.path(root, "inst"))) {
    out <- suppressWarnings(system2("Rscript", c("-e",
      shQuote(sprintf("source('%s'); cat(pkgload::is_dev_package('RiboPilot'))", env_r))),
      env = paste0("RIBOPILOT_SRC=", src), stdout = TRUE, stderr = FALSE))
    expect_true(any(grepl("TRUE", out)), info = src)
  }
})
