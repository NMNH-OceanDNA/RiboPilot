.onAttach <- function(libname, pkgname) {
  version <- nf_installed_version()
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
