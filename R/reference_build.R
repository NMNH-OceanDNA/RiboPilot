# Clade groups in the project database: how they are stored, built, and wired
# into the assembly presets the pipeline actually reads.
#
# A build is a detached process rather than an in-app wait: downloading tens of
# thousands of GenBank records takes minutes, and the app has to stay usable.
# The database row is the only state that survives a closed app, so the row
# carries the pid and log path and the finishing logic is a poll, not a callback.

#' Expected seed/label paths for a clade in a reference folder
#' @return list(seeds, labels).
#' @noRd
reference_db_paths <- function(dir, ref_id) {
  prefix <- paste0(ref_id, "_rrna")
  list(
    seeds = file.path(path.expand(dir), "seeds", paste0(prefix, "_seeds.fasta")),
    labels = file.path(path.expand(dir), "labels", paste0(prefix, "_labels.fasta")))
}

#' Store a suggestion, keeping what an unchanged group already achieved
#'
#' Re-suggesting must not throw away a database somebody waited ten minutes for,
#' so a group whose clade taxid is unchanged keeps its state and paths; a group
#' pointing somewhere new starts over.
#'
#' @param con Open project DB connection.
#' @param groups data.frame from [suggest_reference_groups()]$groups.
#' @param members data.frame(ID, ref_id).
#' @return (invisibly) the number of groups written.
#' @noRd
write_reference_groups <- function(con, groups, members) {
  old <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM reference"), error = function(e) NULL)
  DBI::dbExecute(con, "DELETE FROM reference")
  now <- as.integer(Sys.time())
  for (i in seq_len(nrow(groups))) {
    keep <- NULL
    if (!is.null(old) && nrow(old) > 0) {
      j <- which(old$ref_id == groups$ref_id[i] &
          !is.na(old$clade_taxid) & old$clade_taxid == groups$clade_taxid[i])
      if (length(j) > 0) keep <- old[j[1], ]
    }
    DBI::dbExecute(con,
      "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank, query, n_samples,
         n_records, reference_switch, source, seeds_db, labels_db, log_file, pid, fail_reason,
         time_stamp)
       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
      params = list(groups$ref_id[i], as.integer(groups$clade_taxid[i]), groups$clade_name[i],
        groups$clade_rank[i], clade_query(groups$clade_name[i]),
        as.integer(groups$n_samples[i]), as.integer(groups$n_records[i]),
        if (is.null(keep)) 1L else as.integer(keep$reference_switch),
        if (is.null(keep)) NA_character_ else keep$source,
        if (is.null(keep)) NA_character_ else keep$seeds_db,
        if (is.null(keep)) NA_character_ else keep$labels_db,
        if (is.null(keep)) NA_character_ else keep$log_file,
        if (is.null(keep)) NA_integer_ else keep$pid,
        if (is.null(keep)) NA_character_ else keep$fail_reason,
        now))
  }
  DBI::dbExecute(con, "UPDATE samples SET ref_id = NULL")
  if (nrow(members) > 0) {
    upd <- DBI::dbSendStatement(con, "UPDATE samples SET ref_id = ? WHERE ID = ?")
    on.exit(DBI::dbClearResult(upd), add = TRUE)
    DBI::dbBind(upd, list(members$ref_id, members$ID))
  }
  invisible(nrow(groups))
}

#' Recount `reference.n_samples` from current membership, dropping empty groups
#'
#' @param protect ref_id(s) to keep even at n_samples = 0 (a group just created
#'   or targeted by this call, which may not have members yet).
#' @noRd
refresh_reference_counts <- function(con, protect = character(0)) {
  DBI::dbExecute(con,
    "UPDATE reference SET n_samples =
       (SELECT COUNT(*) FROM samples s WHERE s.ref_id = reference.ref_id)")
  if (length(protect) > 0) {
    ph <- paste(rep("?", length(protect)), collapse = ",")
    DBI::dbExecute(con,
      sprintf("DELETE FROM reference WHERE COALESCE(n_samples, 0) = 0 AND ref_id NOT IN (%s)", ph),
      params = as.list(protect))
  } else {
    DBI::dbExecute(con, "DELETE FROM reference WHERE COALESCE(n_samples, 0) = 0")
  }
  invisible(NULL)
}

#' Point a set of samples at a chosen clade, creating the group if needed
#'
#' @param con Open project DB connection.
#' @param ref_id Group key, or NULL to derive one from `clade_name`.
#' @param clade_taxid,clade_name,clade_rank,n_records The chosen NCBI node.
#' @param member_ids Sample IDs to move into the group.
#' @return (invisibly) the group key.
#' @noRd
set_reference_clade <- function(con, ref_id = NULL, clade_taxid, clade_name, clade_rank,
                               n_records = NA_integer_, member_ids = character(0)) {
  key <- ref_id %||% ref_id_from_name(clade_name)
  exists_row <- nrow(DBI::dbGetQuery(con, "SELECT 1 FROM reference WHERE ref_id = ? LIMIT 1",
    params = list(key))) > 0
  now <- as.integer(Sys.time())
  if (exists_row) {
    # a clade change invalidates whatever was built for the old one
    DBI::dbExecute(con,
      "UPDATE reference SET clade_taxid = ?, clade_name = ?, clade_rank = ?, query = ?,
         n_records = ?, reference_switch = 1, source = NULL, seeds_db = NULL, labels_db = NULL,
         fail_reason = NULL, time_stamp = ?
       WHERE ref_id = ? AND (clade_taxid IS NULL OR clade_taxid != ?)",
      params = list(as.integer(clade_taxid), clade_name, clade_rank, clade_query(clade_name),
        as.integer(n_records), now, key, as.integer(clade_taxid)))
  } else {
    DBI::dbExecute(con,
      "INSERT INTO reference (ref_id, clade_taxid, clade_name, clade_rank, query, n_samples,
         n_records, reference_switch, time_stamp)
       VALUES (?,?,?,?,?,0,?,1,?)",
      params = list(key, as.integer(clade_taxid), clade_name, clade_rank,
        clade_query(clade_name), as.integer(n_records), now))
  }
  if (length(member_ids) > 0) {
    upd <- DBI::dbSendStatement(con, "UPDATE samples SET ref_id = ? WHERE ID = ?")
    DBI::dbBind(upd, list(rep(key, length(member_ids)), member_ids))
    DBI::dbClearResult(upd)
  }
  refresh_reference_counts(con, protect = key)
  invisible(key)
}

#' Mark a group skipped (its samples keep the database they already have)
#' @noRd
skip_reference_group <- function(con, ref_id) {
  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 5, source = 'skipped', fail_reason = NULL,
       time_stamp = ? WHERE ref_id = ?",
    params = list(as.integer(Sys.time()), ref_id))
  invisible(NULL)
}

#' Give a group's samples an assemble_opts preset pointing at its database
#'
#' The preset is named after the group and copies every other setting from the
#' preset the group's samples are on, so switching database changes the database
#' and nothing else. Locked samples keep their preset: their stored results were
#' produced with the old one.
#'
#' @param con Open project DB connection.
#' @param ref_id Group key.
#' @param seeds,labels Absolute paths to the seed and label FASTAs.
#' @param source "built" or "existing".
#' @return (invisibly) the number of samples moved.
#' @noRd
wire_reference_db <- function(con, ref_id, seeds, labels, source = "built") {
  # resolve the base preset FIRST: without one, there is nothing to point at the
  # new database, and moving samples or marking the group done would leave them
  # on a preset with no seeds/labels at all.
  base <- DBI::dbGetQuery(con,
    "SELECT o.* FROM assemble_opts o JOIN assemble a ON a.assemble_opts = o.assemble_opts
     JOIN samples s ON s.ID = a.ID WHERE s.ref_id = ? LIMIT 1", params = list(ref_id))
  if (nrow(base) == 0) {
    base <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts WHERE assemble_opts = 'default'")
  }
  if (nrow(base) == 0) {
    stop("No base assemble_opts preset found for reference group '", ref_id, "'.", call. = FALSE)
  }
  row <- base[1, , drop = FALSE]
  row$assemble_opts <- ref_id
  row$seeds_db <- seeds
  row$labels_db <- labels
  DBI::dbExecute(con, "DELETE FROM assemble_opts WHERE assemble_opts = ?", params = list(ref_id))
  DBI::dbAppendTable(con, "assemble_opts", row)

  ids <- DBI::dbGetQuery(con,
    "SELECT s.ID FROM samples s JOIN assemble a ON a.ID = s.ID
     WHERE s.ref_id = ? AND COALESCE(a.assemble_lock, 0) = 0",
    params = list(ref_id))$ID
  if (length(ids) > 0) {
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbExecute(con,
      sprintf("UPDATE assemble SET assemble_opts = ?, assemble_switch = 1, time_stamp = ?
               WHERE ID IN (%s)", ph),
      params = c(list(ref_id, as.integer(Sys.time())), as.list(ids)))
  }
  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 2, source = ?, seeds_db = ?, labels_db = ?,
       fail_reason = NULL, pid = NULL, time_stamp = ? WHERE ref_id = ?",
    params = list(source, seeds, labels, as.integer(Sys.time()), ref_id))
  invisible(length(ids))
}

#' Point a group at seed/label FASTAs that already exist
#' @return (invisibly) the number of samples moved.
#' @noRd
use_existing_reference <- function(con, ref_id, seeds, labels) {
  seeds <- path.expand(seeds); labels <- path.expand(labels)
  for (f in c(seeds, labels)) {
    if (!file.exists(f)) stop("Reference file not found: ", f, call. = FALSE)
  }
  wire_reference_db(con, ref_id, seeds, labels, source = "existing")
}

#' Launch a detached command, appending its output to `log`
#'
#' Split out so tests can substitute a launcher that never spawns anything.
#' @return the pid.
#' @noRd
processx_launcher <- function(cmd, args, log) {
  p <- processx::process$new(cmd, args, stdout = log, stderr = "2>&1", cleanup = FALSE)
  p$get_pid()
}

#' Is a pid still running?
#' @noRd
process_alive <- function(pid) {
  if (is.na(pid)) return(FALSE)
  isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(as.integer(pid))), error = function(e) FALSE))
}

#' Start a background build for one clade group
#'
#' @param con Open project DB connection.
#' @param project Project directory (the log lives in `<project>/.reference`).
#' @param ref_id Group key.
#' @param api_key,refseq,conda_env,overwrite Passed to [build_reference_db()].
#' @param launcher function(cmd, args, log) -> pid.
#' @return (invisibly) the log path.
#' @noRd
start_reference_build <- function(con, project, ref_id, api_key = NULL, refseq = FALSE,
                                  conda_env = NULL, overwrite = FALSE,
                                  launcher = processx_launcher) {
  row <- DBI::dbGetQuery(con, "SELECT * FROM reference WHERE ref_id = ?", params = list(ref_id))
  if (nrow(row) == 0) stop("No reference group '", ref_id, "'.", call. = FALSE)
  cfg <- ref_config(con)
  logdir <- file.path(project, ".reference")
  dir.create(logdir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(logdir, paste0(ref_id, ".log"))
  unlink(log)

  call <- sprintf(
    "RiboPilot::build_reference_db(query = %s, clade = %s, out_dir = %s, api_key = %s, refseq = %s, conda_env = %s, overwrite = %s)",
    deparse(row$query[1]), deparse(ref_id), deparse(cfg$ref_dir),
    if (is.null(api_key) || !nzchar(api_key)) "NULL" else deparse(api_key),
    if (isTRUE(refseq)) "TRUE" else "FALSE",
    if (is.null(conda_env) || !nzchar(conda_env)) "NULL" else deparse(conda_env),
    if (isTRUE(overwrite)) "TRUE" else "FALSE")

  # captured before the launcher runs, so any FASTA this build writes has an
  # mtime >= launch_time; anything already on disk with an older mtime is a
  # stale leftover from a previous build, not evidence this one succeeded.
  launch_time <- as.integer(Sys.time())
  pid <- launcher(file.path(R.home("bin"), "Rscript"), c("-e", call), log)

  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 4, pid = ?, log_file = ?, fail_reason = NULL,
       time_stamp = ? WHERE ref_id = ?",
    params = list(as.integer(pid), log, launch_time, ref_id))
  invisible(log)
}

#' Last non-empty lines of a build log, for the failure message
#' @noRd
log_tail <- function(log, n = 3) {
  if (is.null(log) || is.na(log) || !file.exists(log)) return("")
  x <- tryCatch(readLines(log, warn = FALSE), error = function(e) character(0))
  x <- trimws(x[nzchar(trimws(x))])
  if (length(x) == 0) return("")
  paste(utils::tail(x, n), collapse = " | ")
}

#' Finish any background build whose process has exited
#'
#' @param con Open project DB connection.
#' @param alive_fn function(pid) -> logical.
#' @return (invisibly) the number of rows changed.
#' @noRd
poll_reference_builds <- function(con, alive_fn = process_alive) {
  running <- tryCatch(DBI::dbGetQuery(con,
    "SELECT * FROM reference WHERE reference_switch = 4"), error = function(e) data.frame())
  if (nrow(running) == 0) return(invisible(0L))
  cfg <- ref_config(con)
  n <- 0L
  for (i in seq_len(nrow(running))) {
    if (isTRUE(alive_fn(running$pid[i]))) next
    ref_id <- running$ref_id[i]
    p <- reference_db_paths(cfg$ref_dir, ref_id)
    launch_time <- suppressWarnings(as.integer(running$time_stamp[i]))
    fresh <- function(f) {
      if (!file.exists(f)) return(FALSE)
      if (is.na(launch_time)) return(TRUE)
      mt <- suppressWarnings(as.integer(file.info(f)$mtime))
      !is.na(mt) && mt >= launch_time
    }
    if (fresh(p$seeds) && fresh(p$labels)) {
      wire_reference_db(con, ref_id, p$seeds, p$labels, source = "built")
      append_ref_manifest(cfg$ref_dir, clade = ref_id, taxid = running$clade_taxid[i],
        rank = running$clade_rank[i], query = running$query[i],
        n_records = running$n_records[i], seeds = p$seeds, labels = p$labels)
    } else {
      DBI::dbExecute(con,
        "UPDATE reference SET reference_switch = 3, pid = NULL, fail_reason = ?, time_stamp = ?
         WHERE ref_id = ?",
        params = list(paste0("Build produced no fresh seed/label FASTAs. ",
          log_tail(running$log_file[i])), as.integer(Sys.time()), ref_id))
    }
    n <- n + 1L
  }
  invisible(n)
}

#' Groups that are neither built nor skipped
#' @noRd
reference_pending <- function(con) {
  tryCatch(DBI::dbGetQuery(con,
    "SELECT ref_id FROM reference WHERE COALESCE(reference_switch, 0) NOT IN (2, 5)")$ref_id,
    error = function(e) character(0))
}
