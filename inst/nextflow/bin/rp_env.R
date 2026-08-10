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

# Positional Rscript arguments: a missing one is NA, NA flows into a query that
# matches nothing, and the script still exits 0, so Nextflow reads a no-op as
# success. Require the mandatory count before any work starts.
rp_args <- function(a, n, usage) {
  ok <- length(a) >= n && !anyNA(a[seq_len(n)]) && all(nzchar(a[seq_len(n)]))
  if (!ok) {
    stop("usage: ", usage, "\n  need ", n, " argument(s), got ", length(a),
      call. = FALSE)
  }
  invisible(a)
}

# A package source checkout, not an installed library directory. Installed
# packages also carry a DESCRIPTION, and pkgload cannot load one.
rp_is_source_dir <- function(d) {
  nzchar(d) && dir.exists(d) &&
    file.exists(file.path(d, "DESCRIPTION")) && dir.exists(file.path(d, "R"))
}

rp_src <- Sys.getenv("RIBOPILOT_SRC")
# accept either the package root or its inst/ (what system.file() returns)
rp_cand <- Filter(rp_is_source_dir, c(rp_src, dirname(rp_src)))
rp_have_pkgload <- requireNamespace("pkgload", quietly = TRUE)

suppressMessages({
  if (length(rp_cand) > 0 && rp_have_pkgload) {
    pkgload::load_all(rp_cand[1], quiet = TRUE, attach_testthat = FALSE)
  } else if (requireNamespace("RiboPilot", quietly = TRUE)) {
    library(RiboPilot)
  } else {
    stop("RiboPilot not installed and RIBOPILOT_SRC not set (or not a package source dir)")
  }
})

# Outside the suppressMessages() above, or the whole point of saying it is lost.
if (length(rp_cand) == 0 && nzchar(rp_src)) {
  message("RIBOPILOT_SRC='", rp_src, "' is not usable here",
    if (!dir.exists(rp_src)) " (not visible from this process)"
    else if (!rp_have_pkgload) " (no pkgload)"
    else " (not a package source directory)",
    "; using the installed RiboPilot instead. Uncommitted package edits are NOT",
    " in this run -- rebuild the image, or run without a container.")
}
