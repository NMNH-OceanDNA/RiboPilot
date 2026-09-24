# Shiny app helpers.

#' Attach external resources (CSS/JS/favicon), mirroring MitoPilot styling
#' @noRd
add_external_resources <- function() {
  shiny::addResourcePath("www", app_sys("app", "www"))
  shiny::tags$head(
    shiny::tags$title("RiboPilot"),
    shiny::tags$link(rel = "shortcut icon", href = "www/favicon.ico"),
    shiny::tags$link(rel = "stylesheet", href = "www/custom.css"),
    shiny::tags$script(src = "www/custom.js"),
    shinyjs::useShinyjs()
  )
}

#' Enumerate Nextflow task work directories for one sample
#'
#' Parses the project's .nextflow.log* files (runs started from the app) and
#' nextflow_logs/*.nextflow.log (submitted jobs) for completed-task records,
#' which Nextflow writes for every attempt (success or failure), each with its
#' own work dir. Mirrors MitoPilot's work-dir browser.
#' @noRd
find_workdirs <- function(project_dir, sample_id) {
  empty <- data.frame(process = character(0), status = character(0),
    modified = character(0), workdir = character(0), exists = logical(0),
    stringsAsFactors = FALSE)
  if (is.null(project_dir) || is.null(sample_id) || !nzchar(sample_id)) return(empty)
  logs <- c(
    list.files(project_dir, pattern = "^\\.nextflow\\.log", full.names = TRUE, all.files = TRUE),
    list.files(file.path(project_dir, JOB_DIRS$nextflow), pattern = "\\.nextflow\\.log",
      full.names = TRUE))
  if (length(logs) == 0) return(empty)
  lines <- unlist(lapply(logs, function(f) tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))))
  th <- grep("Task completed > TaskHandler\\[", lines, value = TRUE)
  if (length(th) == 0) return(empty)

  # perl = TRUE: TRE reads \s inside a bracket expression as a literal backslash
  # and 's', so [^]\s] silently truncated every path at its first "s".
  grab1 <- function(x, pat) vapply(regmatches(x, regexec(pat, x, perl = TRUE)),
    function(v) if (length(v) >= 2) v[2] else NA_character_, character(1))
  name <- grab1(th, "name:\\s*(.*?);")
  exit <- grab1(th, "exit:\\s*(.*?);")
  workdir <- grab1(th, "workDir:\\s*([^\\]\\s]+)")
  pt <- regmatches(name, regexec("([^:[:space:]]+)\\s*\\(([^)]+)\\)\\s*$", name))
  process <- vapply(pt, function(v) if (length(v) >= 3) v[2] else NA_character_, character(1))
  sample <- sub("\\.[0-9]+$", "", vapply(pt, function(v) if (length(v) >= 3) v[3] else NA_character_, character(1)))

  keep <- !is.na(workdir) & !is.na(process) & sample == sample_id
  if (!any(keep)) return(empty)
  out <- data.frame(process = process[keep], exit = exit[keep],
    workdir = trimws(workdir[keep]), stringsAsFactors = FALSE)
  out <- out[!duplicated(out$workdir), , drop = FALSE]
  ex <- trimws(out$exit)
  out$status <- ifelse(is.na(ex) | !nzchar(ex), "unknown",
    ifelse(ex == "0", "success", "failed"))
  out$exists <- dir.exists(out$workdir)
  mt <- file.info(out$workdir)$mtime
  out$modified <- ifelse(is.na(mt), "-", format(mt, "%Y-%m-%d %H:%M"))
  out[, c("process", "status", "modified", "workdir", "exists")]
}

#' Show a directory to the user, by whatever means this session actually has
#'
#' Two different situations, and the desktop assumption is wrong in one of them.
#' On a desktop session the R process and the user are on the same machine, so
#' the platform's file manager is the right answer. In an RStudio Server session,
#' in cluster mode, or when the browser is on another computer, the R process is
#' on a host the user has no desktop on: `xdg-open` there reaches nothing, and
#' RStudio Server intercepts it and serves a URL that 404s. The Files pane is
#' the only file browser that session has in front of the user, so navigate
#' that instead.
#'
#' @return invisibly `TRUE` if a browser was opened, `FALSE` if this session has
#'   no way to show a directory (the caller then has to say where it is).
#' @noRd
open_path <- function(path) {
  if (length(path) != 1 || is.na(path) || !nzchar(path)) return(invisible(FALSE))
  if (app_headless() || !isTRUE(getOption("RiboPilot.local_browser", TRUE))) {
    return(invisible(rstudio_show_dir(path)))
  }
  cmd <- if (Sys.info()[["sysname"]] == "Darwin") "open" else
    if (.Platform$OS.type == "windows") "explorer" else "xdg-open"
  try(system2(cmd, shQuote(path), stdout = FALSE, stderr = FALSE, wait = FALSE), silent = TRUE)
  invisible(TRUE)
}

#' Point the RStudio Files pane at a directory
#'
#' Guarded on every step: rstudioapi is a Suggests, `isAvailable()` is FALSE when
#' the app was launched from a plain R session, and `filesPaneNavigate()` is not
#' in older RStudio. A remote session without RStudio has no file browser at all,
#' which is a real answer rather than an error.
#' @return `TRUE` if the pane was navigated.
#' @noRd
rstudio_show_dir <- function(path) {
  if (!requireNamespace("rstudioapi", quietly = TRUE)) return(FALSE)
  ok <- tryCatch(
    rstudioapi::isAvailable() && rstudioapi::hasFun("filesPaneNavigate"),
    error = function(e) FALSE)
  if (!isTRUE(ok)) return(FALSE)
  tryCatch({
    rstudioapi::filesPaneNavigate(path)
    TRUE
  }, error = function(e) FALSE)
}

#' Active project directory (from the app options, else cwd)
#' @noRd
app_project <- function() {
  shiny::getShinyOption("project") %||% getwd()
}

#' Open a connection to the project database
#'
#' Goes through [open_project_db()] so app connections get the same busy timeout
#' and journal mode as the pipeline's: the app reads and writes the same file
#' while Nextflow tasks are running, and a bare dbConnect has a busy timeout of 0.
#' @noRd
app_con <- function() {
  open_project_db(file.path(app_project(), ".sqlite"))
}

#' Query the project database and disconnect
#' @noRd
app_query <- function(sql, params = NULL) {
  con <- app_con()
  on.exit(DBI::dbDisconnect(con))
  if (is.null(params)) DBI::dbGetQuery(con, sql) else DBI::dbGetQuery(con, sql, params = params)
}

#' Run an app database write, reporting a failure instead of ending the session
#'
#' An error inside an observer closes the Shiny session, and a write can fail
#' while a pipeline task holds the project database.
#' @param expr The write.
#' @return TRUE when `expr` ran, FALSE after an error alert.
#' @noRd
app_write <- function(expr) {
  tryCatch({
    expr
    TRUE
  }, error = function(e) {
    busy <- grepl("locked|busy", conditionMessage(e), ignore.case = TRUE)
    rp_alert(if (busy) {
      paste("The project database is busy, most likely because a pipeline run is",
        "writing to it, so this change was not saved. Wait a moment and try again.")
    } else {
      paste("This change was not saved:", conditionMessage(e))
    }, type = "error", title = if (busy) "Database busy" else "Not saved")
    FALSE
  })
}

#' Is the remote GenBank BLAST switched on for this sample?
#'
#' Lets the modal tell "BLAST was never asked for" apart from "BLAST was asked
#' for and has not come back". Defaults to TRUE when the option row is missing,
#' matching the pipeline's own default (blast.nf).
#' @param id Sample ID.
#' @noRd
blast_enabled <- function(id) {
  if (is.null(id) || !nzchar(id)) {
    return(TRUE)
  }
  d <- app_query(
    "SELECT COALESCE(bo.run_blast, 1) AS run_blast FROM assemble a
       JOIN blast_opts bo ON a.blast_opts = bo.blast_opts WHERE a.ID = ?", params = list(id))
  if (nrow(d) == 0) TRUE else !identical(as.integer(d$run_blast[1]), 0L)
}

#' Per-unit evidence for one sample: assembly, depth, screen hit, taxonomy
#'
#' Everything the reviewer needs at the moment of the keep/drop decision, in one
#' query. All of it is already written by WF1 and none of it was rendered
#' anywhere: `blast_screen` (the unfiltered top core_nt hit, the only column that
#' can say "this is not the host") appeared in no app file, `coverage.mean_depth`
#' was read only by the modal plot, and `blast_ref_candidates` kept up to 50
#' ranked hits of which one was ever shown.
#'
#' @param id Sample ID.
#' @return data.frame with one row per (path, scaffold), ordered path/scaffold,
#'   including a `uid` key, relative depth, and the sibling-concordance columns
#'   from [unit_concordance()]. Empty data.frame when the sample has no units.
#' @noRd
unit_evidence <- function(id) {
  if (is.null(id) || !nzchar(id)) {
    return(data.frame())
  }
  d <- app_query(
    # topology is deliberately not selected: GetOrganelle calls an rDNA result
    # "circular" when it disentangles the assembly graph as a cycle, which a
    # TANDEM ARRAY does by construction. Nuclear rDNA is not a circular
    # molecule, and showing the word next to a scaffold said it was.
    "SELECT a.path, a.scaffold, a.length, COALESCE(a.ignore,0) AS ignore,
       a.triage, a.triage_reason, a.triage_source, a.source_contig,
       a.blast_species, a.blast_accession, a.blast_pident, a.blast_qcovs, a.blast_source,
       ROUND(COALESCE(c.median_depth, c.mean_depth), 1) AS depth,
       CASE WHEN c.median_depth IS NULL THEN 'mean' ELSE 'median' END AS depth_stat,
       ROUND(c.mean_depth, 1) AS depth_mean,
       c.depths AS depths_raw,
       b.species AS screen_species, b.accession AS screen_accession, b.title AS screen_title,
       b.pident AS screen_pident, b.qcovs AS screen_qcovs, b.bitscore AS screen_bitscore,
       b.lineage AS screen_lineage, b.lineage_taxids AS screen_taxids,
       b.rank_distance, b.is_nuclear_rrna,
       sc.screen_verdict, sc.screen_class, sc.screen_confidence, sc.screen_reason, sc.lca_taxon,
       sc.variant_class, sc.variant_group, sc.variant_identity,
       ROUND(ps.score, 3) AS path_score,
       s.lineage_taxids AS sample_taxids, s.sci_name,
       ch.label AS chimera_label, ch.reason AS chimera_reason,
       (SELECT an.itsx_warnings FROM annotate an WHERE an.ID = a.ID AND an.path = a.path
          AND an.scaffold = a.scaffold LIMIT 1) AS itsx_warnings,
       (SELECT GROUP_CONCAT(n.element) FROM
          (SELECT element FROM annotations WHERE ID = a.ID AND path = a.path
             AND scaffold = a.scaffold ORDER BY start) n) AS elements
     FROM assemblies a
     LEFT JOIN coverage c    ON c.ID = a.ID AND c.path = a.path AND c.scaffold = a.scaffold
     LEFT JOIN blast_screen b ON b.ID = a.ID AND b.path = a.path AND b.scaffold = a.scaffold
     LEFT JOIN screen sc      ON sc.ID = a.ID AND sc.path = a.path AND sc.scaffold = a.scaffold
     -- scored per path, shown per scaffold: the whole path is the thing scored
     LEFT JOIN path_scores ps ON ps.ID = a.ID AND ps.path = a.path
     LEFT JOIN samples s      ON s.ID = a.ID
     LEFT JOIN chimera ch     ON ch.ID = a.ID AND ch.path = a.path AND ch.scaffold = a.scaffold
     WHERE a.ID = ? ORDER BY a.path, a.scaffold", params = list(id))
  if (nrow(d) == 0) {
    return(d)
  }
  d$uid <- paste(id, d$path, d$scaffold, sep = "|")
  d$unit <- sprintf("p%s.s%s", d$path, d$scaffold)
  d$depth_track <- vapply(d$depths_raw %||% rep(NA_character_, nrow(d)),
    depth_track_string, character(1), USE.NAMES = FALSE)
  d$depths_raw <- NULL
  # Median per position, not mean: an rDNA scaffold's conserved 18S/28S cores
  # recruit multi-mapping reads that the spacers do not, so a mean is inflated in
  # proportion to how much core the scaffold contains -- which is composition,
  # not abundance. Pre-v10 rows have no median and fall back to the mean until
  # backfill_coverage_median() runs; depth_stat records which is in play.
  #
  # Either way this is a WITHIN-SAMPLE ratio only. coverage.nf maps one path at
  # a time, so depths compare within a path, not across paths.
  mx <- suppressWarnings(max(d$depth, na.rm = TRUE))
  d$depth_max <- if (is.finite(mx) && mx > 0) mx else NA_real_
  d$shares <- rank_distance_label(d$rank_distance)
  cbind(d, unit_concordance(d))
}

#' Sibling concordance: is this scaffold the odd one out among its own siblings?
#'
#' The check that fires exactly where the shipped `lineages_conflict()`
#' structurally cannot: a copepod, epibiont or gut-content contaminant on a
#' metazoan sample shares "Metazoa" with the host, so a shared-token test reports
#' it clean. Comparing a sample's scaffolds to EACH OTHER needs no absolute
#' threshold at all.
#'
#' Concordance is taken at the DEEPEST RANK AT WHICH A STRICT MAJORITY EXISTS,
#' not at a fixed rank. That keeps the copepod-on-fish catch (a majority resolves
#' at phylum: many Chordata vs one Arthropoda) while going silent on the sparse
#' clades the tool exists for, where a genus-level test would flag every unit
#' because each scaffold's best hit is a different genus.
#'
#' @param d Unit rows carrying `screen_taxids` (rank-aligned hit taxid paths).
#' @return data.frame(sib_rank, sib_taxon, odd_one_out) with `nrow(d)` rows.
#' @noRd
unit_concordance <- function(d) {
  n <- nrow(d)
  out <- data.frame(sib_rank = rep(NA_character_, n), sib_taxon = rep(NA_character_, n),
    odd_one_out = rep(FALSE, n), stringsAsFactors = FALSE)
  tx <- d$screen_taxids %||% rep(NA_character_, n)
  # only units with a resolved hit lineage can vote
  votes <- which(!is.na(tx) & nzchar(tx))
  if (length(votes) < 3) {
    return(out)
  }
  mat <- vapply(tx[votes], parse_rank_taxids, integer(length(RRNA_LINEAGE_RANKS)))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(RRNA_LINEAGE_RANKS))
  names_mat <- vapply(d$screen_lineage[votes] %||% rep(NA_character_, length(votes)),
    parse_rank_names, character(length(RRNA_LINEAGE_RANKS)))
  if (is.null(dim(names_mat))) names_mat <- matrix(names_mat, nrow = length(RRNA_LINEAGE_RANKS))

  # deepest rank (species is deepest) with a STRICT majority among resolved votes
  for (p in rev(seq_along(RRNA_LINEAGE_RANKS))) {
    v <- mat[p, ]
    v <- v[!is.na(v)]
    if (length(v) < 3) next
    tab <- table(v)
    top <- names(tab)[which.max(tab)]
    if (max(tab) <= length(v) / 2) next          # strict majority only
    keep <- !is.na(mat[p, ]) & mat[p, ] == as.integer(top)
    out$sib_rank[votes] <- RRNA_LINEAGE_RANKS[p]
    out$sib_taxon[votes] <- names_mat[p, which(keep)[1]]
    out$odd_one_out[votes] <- !keep
    break
  }
  out
}

#' Put text on the user's clipboard from the server side
#'
#' The Clipboard API needs a secure context, which a Shiny app served over plain
#' http on anything but localhost is not, so the textarea + execCommand path is
#' a real fallback rather than legacy-browser politeness. jsonlite does the
#' escaping: a FASTA blob is full of newlines and would otherwise break the
#' generated statement.
#' @param text Character scalar to copy.
#' @noRd
copy_to_clipboard <- function(text) {
  js_str <- as.character(jsonlite::toJSON(text, auto_unbox = TRUE))
  shinyjs::runjs(sprintf("
    (function(t){
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(t);
      } else {
        var a = document.createElement('textarea');
        a.value = t; a.style.position = 'fixed'; a.style.opacity = '0';
        document.body.appendChild(a); a.focus(); a.select();
        try { document.execCommand('copy'); } finally { document.body.removeChild(a); }
      }
    })(%s);", js_str))
  invisible(text)
}

#' Did this sample supply its own assembly?
#'
#' Drives the two pieces of the app that only make sense for one: the source
#' contig column and the list of contigs the rDNA search left behind. FALSE for
#' a project whose mapping file never had an `Assembly` column.
#'
#' @param id Sample ID.
#' @return TRUE when the sample names an assembly file.
#' @noRd
sample_has_user_assembly <- function(id) {
  a <- tryCatch(
    app_query("SELECT Assembly FROM samples WHERE ID = ?", params = list(id))$Assembly,
    error = function(e) NULL)
  length(a) > 0 && !is.na(a[1]) && nzchar(a[1])
}

#' Contigs the rDNA search found but did not carry forward
#'
#' The kept contigs become assembly units and are shown in the scaffold table;
#' these are the ones that never got there, which is the only place a user can
#' see that they existed at all.
#'
#' @param id Sample ID.
#' @return data.frame(contig, length, rrna_features, reason), empty when the
#'   sample dropped nothing.
#' @noRd
rejected_rrna_contigs <- function(id) {
  d <- tryCatch(app_query(
    "SELECT contig, length, rrna_features, reason FROM rrna_candidates
     WHERE ID = ? AND COALESCE(selected,0) = 0 ORDER BY rank", params = list(id)),
    error = function(e) NULL)
  d %||% data.frame(contig = character(0), length = integer(0),
    rrna_features = character(0), reason = character(0), stringsAsFactors = FALSE)
}
