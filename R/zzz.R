#' Session cache for the Nextflow version probe
#'
#' dbplyr is never called by name: dplyr::tbl() on a DBI connection loads it
#' to translate, so it is a real runtime dependency that R CMD check's static
#' analysis cannot see.
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
  # Hydra sessions lack Java, qsub, and ~/bin on PATH; fix before probing Nextflow
  tryCatch(ensure_hydra_setup(), error = function(e) {
    packageStartupMessage("Hydra setup failed: ", conditionMessage(e),
      ". Run hydra_setup() by hand.")
  })
  version <- nf_version_session()
  status <- nf_version_status(version)
  packageStartupMessage("Welcome to RiboPilot!")
  switch(status,
    missing = if (nf_on_path()) {
      packageStartupMessage("A 'nextflow' executable was found but 'nextflow -version' did not report a version.")
      packageStartupMessage("Ensure Java 17+ is on PATH in this R session (or set NXF_JAVA_HOME).")
    } else {
      packageStartupMessage("A Nextflow installation is needed to run the RiboPilot pipeline.")
      packageStartupMessage("Please install Nextflow (>= ", NF_MIN_SUPPORTED, ") from https://www.nextflow.io/")
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
  ".data", "ID", "R1", "R2", "depth", "element", "mid", "pos", "start",
  "label", "taxon", "pident", "x", "src"
))
