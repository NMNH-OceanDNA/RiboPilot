# Nextflow version compatibility for the RiboPilot pipeline.
#
# The supported range is aligned with MitoPilot so a single Nextflow install
# serves both packages. RiboPilot's bounds are NOT the same hard breakages
# MitoPilot has, though:
#   * Lower bound 24.10.0: a tested-good floor. RiboPilot does not rely on
#     `workflow.failOnIgnore` (it sets no failOnIgnore; its errorStrategy
#     'ignore' steps are genuinely best-effort), so this is a validated-range
#     floor rather than a known failure below it.
#   * Upper bound (broken from 26.0.0): PRECAUTIONARY. Nextflow 26 removed
#     `NF.isDsl2()`, which breaks the nf-sqldb plugin -- but RiboPilot uses no
#     nf-sqldb (processes open the project SQLite by path). RiboPilot is simply
#     unverified on the 26.x line, so runs are pinned down rather than trusted.

#' Minimum supported Nextflow version.
#' @noRd
NF_MIN_SUPPORTED <- "24.10.0"

#' Latest verified-good Nextflow release (pin target when the installed version
#' is outside the supported range).
#' @noRd
NF_MAX_SUPPORTED <- "25.10.6"

#' First Nextflow release line RiboPilot has not verified (>= this is "too_new").
#' @noRd
NF_BROKEN_FROM <- "26.0.0"

#' Human-readable supported range, for messages.
#' @noRd
nf_supported_label <- function() {
  paste0(NF_MIN_SUPPORTED, " - 25.10.x")
}

#' Is a `nextflow` executable on PATH?
#'
#' Lets callers tell "not installed" apart from "installed but `nextflow
#' -version` failed" (e.g. Java missing from the session), so the guidance can
#' point at the real problem.
#' @param nextflow The nextflow binary name/path.
#' @return `TRUE` if found on PATH.
#' @noRd
nf_on_path <- function(nextflow = "nextflow") {
  nzchar(Sys.which(nextflow))
}

#' Raw output of `nextflow -version`
#'
#' stderr is merged and the exit status ignored: a fresh launcher prints
#' "Downloading nextflow dependencies..." (and may exit non-zero) before the
#' banner, and some builds print the banner on stderr.
#' @param nextflow The nextflow binary name/path.
#' @noRd
nf_version_output <- function(nextflow = "nextflow") {
  suppressWarnings(tryCatch(
    system2(nextflow, args = "-version", stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  ))
}

#' Extract the `X.Y.Z` version from `nextflow -version` output, dropping any
#' `-edge` / build suffix.
#' @return A version string like "25.10.6", or `NA_character_` if not parseable.
#' @noRd
nf_parse_version <- function(out) {
  if (length(out) == 0) {
    return(NA_character_)
  }
  m <- stringr::str_match(out, "version\\s+([0-9]+\\.[0-9]+\\.[0-9]+)")
  v <- stats::na.omit(m[, 2])
  if (length(v) == 0) NA_character_ else v[1]
}

#' Parse the installed Nextflow version
#' @param nextflow The nextflow binary name/path.
#' @noRd
nf_installed_version <- function(nextflow = "nextflow") {
  nf_parse_version(nf_version_output(nextflow))
}

#' Concise error output from a failed `nextflow -version`
#'
#' Surfaced in the "found but no version" message so the real cause (most often
#' "Cannot find Java or it's a wrong version") is visible. Keeps the informative
#' lines (ERROR / NOTE / Java), trimmed short.
#' @param out Output lines from [nf_version_output()].
#' @return A short string, or "" if nothing useful was captured.
#' @noRd
nf_probe_error <- function(out) {
  if (length(out) == 0) {
    return("")
  }
  hit <- out[grepl("ERROR|NOTE|Java|JAVA", out, ignore.case = TRUE)]
  if (length(hit) == 0) hit <- out
  paste(utils::head(trimws(hit), 4), collapse = "\n")
}

#' Classify a Nextflow version against the supported range
#'
#' Pure: pass a version string to test without a live Nextflow.
#' @param version Version string (defaults to the installed Nextflow).
#' @return One of "missing", "too_old", "supported", "too_new".
#' @noRd
nf_version_status <- function(version = nf_installed_version()) {
  if (length(version) == 0 || is.na(version)) {
    return("missing")
  }
  v <- numeric_version(version)
  if (v < numeric_version(NF_MIN_SUPPORTED)) {
    return("too_old")
  }
  if (v >= numeric_version(NF_BROKEN_FROM)) {
    return("too_new")
  }
  "supported"
}

#' Choose the NXF_VER value for a given installed version (pure)
#'
#' In-range installs pin to themselves so cluster nodes without internet do not
#' download a different build; the original string is preserved (a
#' `numeric_version` round-trip would drop a leading zero, e.g. "25.04.8" ->
#' "25.4.8", which NXF_VER rejects). Out-of-range installs pin to the latest
#' verified-good release.
#' @param version Installed version string.
#' @return A version string for `NXF_VER`, or `NA_character_` if missing.
#' @noRd
nf_pin_for <- function(version) {
  if (length(version) == 0 || is.na(version)) {
    return(NA_character_)
  }
  if (nf_version_status(version) == "supported") {
    return(version)
  }
  NF_MAX_SUPPORTED
}

#' NXF_VER value to pin for a run (detects, then [nf_pin_for()]).
#' @param nextflow The nextflow binary name/path.
#' @return A version string for `NXF_VER`, or `NA_character_` if Nextflow missing.
#' @noRd
nf_pin_version <- function(nextflow = "nextflow") {
  nf_pin_for(nf_installed_version(nextflow))
}

#' Shell prefix that pins `NXF_VER` for a `nextflow` invocation
#'
#' Returns `"NXF_VER='<ver>' "` to prepend to a `nextflow ...` command string
#' (used by the CLI drivers, which run via `bash -c`), or `""` when Nextflow is
#' missing so nothing is pinned.
#' @param nextflow The nextflow binary name/path.
#' @noRd
nxf_ver_prefix <- function(nextflow = "nextflow") {
  pin <- nf_pin_version(nextflow)
  if (is.na(pin)) "" else paste0("NXF_VER=", shQuote(pin), " ")
}

#' Update hint shown when the installed Nextflow is out of range.
#' @noRd
nf_update_hint <- function() {
  "Update Nextflow, e.g. mamba install -n ribopilot 'bioconda::nextflow>=24.10,<26'"
}

#' Are we inside R's install / build / check / dev-load machinery?
#'
#' `R CMD INSTALL` / `check` attach the namespace in a load-test, so an
#' unconditional `stop()` in `.onAttach` would break installation on any machine
#' without a compatible Nextflow (including CI). `devtools::load_all()` likewise
#' runs `.onAttach`. These env vars are set only by that machinery.
#' @noRd
in_r_build_context <- function() {
  nzchar(Sys.getenv("R_INSTALL_PKG")) ||
    nzchar(Sys.getenv("R_CMD")) ||
    nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
    nzchar(Sys.getenv("DEVTOOLS_LOAD"))
}

#' Check the installed Nextflow version at a pipeline entry point
#'
#' Stops (or warns) when the installed Nextflow is below the supported range,
#' warns when above it (runs are pinned via `NXF_VER`). Called from the app
#' launcher and the `run_*()` pipeline drivers.
#'
#' @param context Short label for where the check runs (e.g. "run_assembly"),
#'   included in messages.
#' @param nextflow The nextflow binary name/path to probe.
#' @param on_too_old "stop" (default) errors out; "warn" only warns.
#' @return Invisibly, the version status string.
#' @noRd
check_nextflow_version <- function(context = NULL, nextflow = "nextflow",
                                   on_too_old = c("stop", "warn")) {
  on_too_old <- match.arg(on_too_old)
  probe <- nf_version_output(nextflow)
  version <- nf_parse_version(probe)
  status <- nf_version_status(version)
  ctx <- if (is.null(context)) "" else paste0(" (", context, ")")

  switch(status,
    missing = warning(
      if (nf_on_path(nextflow)) {
        err <- nf_probe_error(probe)
        glue::glue(
          "A 'nextflow' executable was found but 'nextflow -version' did not ",
          "report a version{ctx}.",
          if (nzchar(err)) "\nNextflow reported:\n{err}" else "",
          "\nEnsure Java 17+ is on PATH in this R session (or set NXF_JAVA_HOME); ",
          "if you just installed Nextflow, run `nextflow -version` once in a ",
          "terminal to finish its download. RiboPilot needs Nextflow ",
          "{nf_supported_label()}."
        )
      } else {
        glue::glue(
          "Nextflow not found{ctx}. RiboPilot needs Nextflow ",
          "{nf_supported_label()}. Install from https://www.nextflow.io/"
        )
      },
      call. = FALSE
    ),
    too_old = {
      msg <- glue::glue(
        "Nextflow {version} is below the minimum supported version ",
        "({NF_MIN_SUPPORTED}){ctx}. {nf_update_hint()}"
      )
      if (on_too_old == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    },
    too_new = warning(
      glue::glue(
        "Nextflow {version} is newer than RiboPilot has verified{ctx} ",
        "(supported {nf_supported_label()}). Runs will be pinned to Nextflow ",
        "{NF_MAX_SUPPORTED} via NXF_VER. To avoid re-downloading an engine, ",
        "install a supported version: {nf_update_hint()}"
      ),
      call. = FALSE
    ),
    supported = invisible(NULL)
  )
  invisible(status)
}
