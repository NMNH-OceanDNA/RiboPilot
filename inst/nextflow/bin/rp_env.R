# Shared bootstrap sourced by the WF1/WF2 R-script processes.
#
# RIBOPILOT_SRC wins over an installed package. nf_command() passes
# --ribopilot_src only from a source checkout, so a stale install cannot shadow
# the code being developed.
#
# Otherwise the installed RiboPilot runs (in a container, the image's). The app
# warns when a project's image is for another RiboPilot version.

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

# A package source checkout, not an installed library directory. An installed
# package carries a DESCRIPTION and an R/ directory too, but that R/ holds the
# lazy-load database rather than sources, so the .R files are what tell them
# apart. Loading an installed directory as source leaves an empty namespace and
# every process dies on "could not find function".
rp_is_source_dir <- function(d) {
  nzchar(d) && dir.exists(d) && file.exists(file.path(d, "DESCRIPTION")) &&
    length(list.files(file.path(d, "R"), pattern = "\\.[Rr]$")) > 0
}

rp_src <- Sys.getenv("RIBOPILOT_SRC")
# accept either the package root or its inst/ (what system.file() returns)
rp_cand <- Filter(rp_is_source_dir, c(rp_src, dirname(rp_src)))
rp_have_pkgload <- requireNamespace("pkgload", quietly = TRUE)
rp_loaded_src <- length(rp_cand) > 0 && rp_have_pkgload

suppressMessages({
  if (rp_loaded_src) {
    pkgload::load_all(rp_cand[1], quiet = TRUE, attach_testthat = FALSE)
  } else if (requireNamespace("RiboPilot", quietly = TRUE)) {
    library(RiboPilot)
  } else {
    stop("RiboPilot not installed and RIBOPILOT_SRC not set (or not a package source dir)")
  }
})

# Outside the suppressMessages() above, or the whole point of saying it is lost.
if (!rp_loaded_src && nzchar(rp_src)) {
  message("RIBOPILOT_SRC='", rp_src, "' is not usable here",
    if (!dir.exists(rp_src)) " (not visible from this process)"
    else if (length(rp_cand) == 0) " (not a package source directory)"
    else " (no pkgload)",
    "; using the installed RiboPilot instead. Uncommitted package edits are NOT",
    " in this run -- rebuild the image, or run without a container.")
}
