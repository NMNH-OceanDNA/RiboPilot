# Shared bootstrap sourced by the WF1/WF2 R-script processes.
#
# RIBOPILOT_SRC WINS over an installed package. nf_command() always passes
# --ribopilot_src for the running package, and preferring the installed copy
# meant a stale install silently shadowed it: the development version is
# 0.0.0.9000 and so is every older install, so nothing detects the mismatch and
# the pipeline runs old code with new data. That surfaces as an argument-
# signature error at best and wrong results at worst.
#
# The source directory is often invisible from inside a container, which is not
# an error: the image ships an installed RiboPilot and that is what should run.
# It IS worth saying out loud, because "my edit had no effect" is otherwise
# indistinguishable from a bug in the edit.

# A package source checkout, not an installed library directory. Installed
# packages also carry a DESCRIPTION, and pkgload cannot load one.
rp_is_source_dir <- function(d) {
  nzchar(d) && dir.exists(d) &&
    file.exists(file.path(d, "DESCRIPTION")) && dir.exists(file.path(d, "R"))
}

suppressMessages({
  src <- Sys.getenv("RIBOPILOT_SRC")
  # accept either the package root or its inst/ (what system.file() returns)
  cand <- Filter(rp_is_source_dir, c(src, dirname(src)))
  if (length(cand) > 0 && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(cand[1], quiet = TRUE, attach_testthat = FALSE)
  } else if (requireNamespace("RiboPilot", quietly = TRUE)) {
    if (nzchar(src)) {
      message("RIBOPILOT_SRC='", src, "' is not usable here",
        if (!dir.exists(src)) " (not visible from this process)" else
          if (!requireNamespace("pkgload", quietly = TRUE)) " (no pkgload)" else
            " (not a package source directory)",
        "; using the installed RiboPilot instead. Uncommitted package edits will",
        " NOT be in this run -- rebuild the image, or run without a container.")
    }
    library(RiboPilot)
  } else {
    stop("RiboPilot not installed and RIBOPILOT_SRC not set (or not a package source dir)")
  }
})
