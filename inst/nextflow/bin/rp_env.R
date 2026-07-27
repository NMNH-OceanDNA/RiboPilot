# Shared bootstrap sourced by the WF1/WF2 R-script processes.
#
# RIBOPILOT_SRC WINS over an installed package. nf_command() always passes
# --ribopilot_src for the running package, and preferring the installed copy
# meant a stale install silently shadowed it: the development version is
# 0.0.0.9000 and so is every older install, so nothing detects the mismatch and
# the pipeline runs old code with new data. That surfaces as an argument-
# signature error at best and wrong results at worst.
suppressMessages({
  src <- Sys.getenv("RIBOPILOT_SRC")
  if (nzchar(src) && file.exists(file.path(src, "DESCRIPTION"))) {
    pkgload::load_all(src, quiet = TRUE, attach_testthat = FALSE)
  } else if (requireNamespace("RiboPilot", quietly = TRUE)) {
    library(RiboPilot)
  } else {
    stop("RiboPilot not installed and RIBOPILOT_SRC not set (or not a package source dir)")
  }
})
