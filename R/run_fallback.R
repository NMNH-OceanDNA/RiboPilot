# Assembly-fallback driver. Points samples that failed the GetOrganelle (B1) QC
# gate at the baiting / de-novo assembler, then runs the ordinary assembly
# workflow. The fallback used to be a separate Nextflow workflow; it is now a
# per-sample assembler inside WF1 (dev/design/14_ASSEMBLER_SWAP.md).

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

#' Point one sample at an assemble_opts preset using the requested assembler
#'
#' Derives a new preset rather than editing the sample's current one in place:
#' presets are shared, so editing would change the assembler for every other
#' sample using it. The derived preset copies every other column from its
#' source, and is reused on a second call rather than duplicated.
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param method "bait" or "denovo".
#' @return (invisibly) the preset name the sample now points at.
#' @noRd
set_sample_assembler <- function(con, id, method) {
  src <- DBI::dbGetQuery(con,
    "SELECT assemble_opts FROM assemble WHERE ID = ?", params = list(id))$assemble_opts
  if (length(src) == 0 || is.na(src[1])) src <- "default" else src <- src[1]
  derived <- paste0(src, "_", method)

  exists_already <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemble_opts WHERE assemble_opts = ?",
    params = list(derived))$n > 0
  if (!exists_already) {
    row <- DBI::dbGetQuery(con,
      "SELECT * FROM assemble_opts WHERE assemble_opts = ?", params = list(src))
    if (nrow(row) == 0) {
      row <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts LIMIT 1")
    }
    if (nrow(row) > 0) {
      row$assemble_opts <- derived
      row$assembler <- method
      DBI::dbAppendTable(con, "assemble_opts", row)
    }
  } else {
    DBI::dbExecute(con, "UPDATE assemble_opts SET assembler = ? WHERE assemble_opts = ?",
      params = list(method, derived))
  }

  DBI::dbExecute(con, "UPDATE assemble SET assemble_opts = ? WHERE ID = ?",
    params = list(derived, id))
  invisible(derived)
}

#' Run an assembly fallback (Option B2 baiting / Option A de-novo)
#'
#' Points samples at the baiting / de-novo assembler and runs the ordinary
#' assembly workflow. By default it targets samples that failed the primary
#' GetOrganelle assembly QC gate (no assembly, or a missing 18S/28S anchor);
#' pass `ids` to force a specific set (e.g. an Option A whole-assembly
#' cross-check of a good sample).
#'
#' Each targeted sample is repointed at a derived `assemble_opts` preset named
#' `<current>_<method>`, so other samples sharing its original preset are
#' unaffected.
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
  path <- normalizePath(path)
  stopifnot(file.exists(file.path(path, ".sqlite")), file.exists(file.path(path, ".config")))

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  if (is.null(ids)) ids <- fallback_candidates(con)
  if (length(ids) == 0) {
    DBI::dbDisconnect(con)
    message("No samples need an assembly fallback (none failed the QC gate).")
    return(invisible(0L))
  }
  for (id in ids) {
    set_sample_assembler(con, id, method)
    DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 1 WHERE ID = ?", params = list(id))
  }
  DBI::dbDisconnect(con)

  message("Queued ", length(ids), " sample(s) for the ", method, " assembler.")
  run_assembly(path = path, nextflow = nextflow, work = work,
    conda_env = conda_env, ribopilot_src = ribopilot_src)
}
