#' Session cache for the Nextflow version probe
#'
#' The crayon import is declared here because the startup banner calls
#' `crayon::` from inside glue template strings, which R CMD check's static
#' analysis cannot see; without it the package looks unused.
#' dbplyr is likewise never called by name: dplyr::tbl() on a DBI connection
#' loads it to translate, so it is a real runtime dependency that static
#' analysis cannot see either.
#' @importFrom crayon bold
#' @importFrom dbplyr sql
#' @noRd
.nf_probe_cache <- new.env(parent = emptyenv())

#' Installed Nextflow version, probed at most once per R session
#'
#' A cold Nextflow launcher downloads its dependencies before printing a
#' banner, so probing on every attach can block for tens of seconds.
#' @noRd
nf_version_session <- function() {
  if (!exists("version", envir = .nf_probe_cache, inherits = FALSE)) {
    assign("version", nf_installed_version(), envir = .nf_probe_cache)
  }
  get("version", envir = .nf_probe_cache, inherits = FALSE)
}

.onAttach <- function(libname, pkgname) {
  version <- nf_version_session()
  status <- nf_version_status(version)
  # Refuse to load on an unsupported-old Nextflow, but never during package
  # build/check/dev-load (that final load-test would fail on CI/install).
  if (status == "too_old" && !in_r_build_context()) {
    stop(glue::glue(
      "Nextflow {version} is below the minimum supported version ",
      "({NF_MIN_SUPPORTED}). {nf_update_hint()}"), call. = FALSE)
  }
  glue::glue("Welcome to {crayon::bold(crayon::white('{RiboPilot}'))}!") |> packageStartupMessage()
  switch(status,
    missing = {
      packageStartupMessage("A Nextflow installation is needed to run the RiboPilot pipeline.")
      glue::glue("Please install Nextflow (>= {NF_MIN_SUPPORTED}) from ",
        "{crayon::underline('https://www.nextflow.io/')}") |> packageStartupMessage()
    },
    too_old = packageStartupMessage(glue::glue(
      "Nextflow {version} is below the supported range ({nf_supported_label()}). ",
      "{nf_update_hint()}")),
    too_new = packageStartupMessage(glue::glue(
      "Using Nextflow {version} (newer than the verified range {nf_supported_label()}); ",
      "runs will be pinned to {NF_MAX_SUPPORTED} via NXF_VER.")),
    supported = glue::glue("Using Nextflow {version} (supported {nf_supported_label()})") |>
      packageStartupMessage()
  )
  invisible()
}

# Column names used non-standardly inside dplyr / ggplot2 aes(); R CMD check
# cannot see that they are data columns rather than missing globals.
utils::globalVariables(c(
  ".data", "ID", "R1", "R2", "depth", "element", "mid", "pos", "start"
))
