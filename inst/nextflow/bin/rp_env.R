# Shared bootstrap sourced by the WF1/WF2 R-script processes.
#
# RIBOPILOT_SRC wins over an installed package. nf_command() passes
# --ribopilot_src only from a source checkout, so a stale install cannot shadow
# the code being developed.
#
# Otherwise the installed RiboPilot runs (in a container, the image's). The
# scripts in this folder come from the RiboPilot that started the run, and the
# two must be the same version and build.

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

# NULL when the scripts in `bin` and the RiboPilot in `img` are the same build of
# version `want`, else a message saying what differs
rp_mismatch <- function(bin, img, want, have) {
  f <- list.files(bin, pattern = "\\.(R|sh)$")
  same <- identical(f, list.files(img, pattern = "\\.(R|sh)$")) &&
    identical(unname(tools::md5sum(file.path(bin, f))), unname(tools::md5sum(file.path(img, f))))
  if (!nzchar(want)) want <- have
  if (want == have && same) return(NULL)
  ctr <- file.exists("/.dockerenv") || nzchar(Sys.getenv("SINGULARITY_NAME")) ||
    nzchar(Sys.getenv("APPTAINER_NAME"))
  fix <- if (!ctr) {
    paste0(" Install the same RiboPilot in the R library this task uses (",
      dirname(dirname(dirname(img))), ").")
  } else if (want != have) {
    paste0(" Set container in the project's .config to the image made from RiboPilot ", want,
      " (macguigand/ribopilot:", want, ", or your own build of the same code).")
  } else {
    paste0(" Build an image from the RiboPilot that starts the run (docker/deploy-local.sh in its",
      " source folder) and set container in the project's .config to it.")
  }
  paste0("Version mismatch: this run was started by RiboPilot ", want, ", but this task loaded RiboPilot ",
    have, if (want == have) " built from different code", ".", fix, " Then run again.")
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

if (!rp_loaded_src && exists("rp_file")) {
  rp_why <- rp_mismatch(dirname(rp_file), system.file("nextflow", "bin", package = "RiboPilot"),
    Sys.getenv("RIBOPILOT_VERSION"), as.character(utils::packageVersion("RiboPilot")))
  if (!is.null(rp_why)) {
    message(rp_why)
    quit(save = "no", status = 1)
  }
}
