# Assembly-fallback driver. Queues samples that failed the GetOrganelle (B1) QC
# gate (or an explicit set) for a reference-guided baiting / de-novo fallback,
# then launches the self-contained main_fallback.nf workflow.

#' Samples that failed the assembly QC gate (candidates for a fallback)
#'
#' Failed assembly (`assemble_switch = 3`) or annotated but missing a conserved
#' anchor (18S / 28S). @param con Open DBI connection. @return character IDs.
#' @noRd
fallback_candidates <- function(con) {
  DBI::dbGetQuery(con,
    "SELECT a.ID FROM assemble a LEFT JOIN annotate an ON a.ID = an.ID
     WHERE a.assemble_lock = 0
       AND (a.assemble_switch = 3
            OR an.missing LIKE '%18S%' OR an.missing LIKE '%28S%')")$ID
}

#' Run an assembly fallback (Option B2 baiting / Option A de-novo)
#'
#' Queues samples for the fallback workflow and launches it. By default it
#' targets samples that failed the primary GetOrganelle assembly QC gate (no
#' assembly, or a missing 18S/28S anchor); pass `ids` to force a specific set
#' (e.g. an Option A whole-assembly cross-check of a good sample).
#'
#' @param path Project directory (contains `.sqlite` and `.config`).
#' @param ids Sample IDs to run (default: the QC-gate failures).
#' @param method "bait" (B2 reference-guided baiting, default) or "denovo"
#'   (Option A, assemble all reads).
#' @param nextflow Path to the `nextflow` binary.
#' @param work Nextflow work directory (default `<path>/work`).
#' @param conda_env Optional conda env providing the tools; its `bin/` is
#'   prepended to PATH for the run.
#' @param ribopilot_src Package source for the R-script processes.
#' @return (invisibly) the Nextflow exit status (or 0 if nothing queued).
#' @export
run_fallback <- function(path = ".", ids = NULL, method = c("bait", "denovo"),
                         nextflow = "nextflow", work = file.path(path, "work"),
                         conda_env = NULL, ribopilot_src = NULL) {
  method <- match.arg(method)
  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  path <- normalizePath(path)
  stopifnot(file.exists(file.path(path, ".sqlite")), file.exists(file.path(path, ".config")))
  check_nextflow_version("run_fallback", nextflow = nextflow)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  if (is.null(ids)) ids <- fallback_candidates(con)
  if (length(ids) == 0) {
    DBI::dbDisconnect(con)
    message("No samples need an assembly fallback (none failed the QC gate).")
    return(invisible(0L))
  }
  for (id in ids) {
    DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 5 WHERE ID = ?", params = list(id))
  }
  DBI::dbDisconnect(con)

  message("Running Nextflow assembly fallback (", method, ") for ", length(ids), " sample(s)...")
  args <- nf_command(path, "fallback", resume = TRUE, work = work,
    ribopilot_src = ribopilot_src, method = method)
  nf_cmd <- paste0(nxf_ver_prefix(nextflow), paste(c(shQuote(nextflow), shQuote(args)), collapse = " "))
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) {
    warning("Nextflow (fallback) exited with status ", status)
  } else {
    message("Fallback complete (assemble -> score -> stitch -> BLAST + coverage) written to the database. ",
      "Review and lock assemblies in the app, then run annotation (WF2).")
  }
  invisible(status)
}
