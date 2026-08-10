# Clade groups in the project database: how they are stored, built, and wired
# into the assembly presets the pipeline actually reads.
#
# A build is a detached process rather than an in-app wait: downloading tens of
# thousands of GenBank records takes minutes, and the app has to stay usable.
# The database row is the only state that survives a closed app, so the row
# carries the pid and log path and the finishing logic is a poll, not a callback.

#' Expected seed/label paths for a clade in a reference folder
#'
#' The single definition of the naming scheme: the builder writes these paths
#' and the build poll looks for them, so both must read it from here.
#' @return list(prefix, seeds, labels).
#' @noRd
reference_db_paths <- function(dir, ref_id) {
  prefix <- paste0(ref_id, "_rrna")
  list(
    prefix = prefix,
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
         n_records, reference_switch, source, seeds_db, labels_db, log_file, pid, preset_name,
         fail_reason, time_stamp)
       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
      params = list(groups$ref_id[i], as.integer(groups$clade_taxid[i]), groups$clade_name[i],
        groups$clade_rank[i], clade_query(groups$clade_name[i]),
        as.integer(groups$n_samples[i]), as.integer(groups$n_records[i]),
        if (is.null(keep)) 1L else as.integer(keep$reference_switch),
        if (is.null(keep)) NA_character_ else keep$source,
        if (is.null(keep)) NA_character_ else keep$seeds_db,
        if (is.null(keep)) NA_character_ else keep$labels_db,
        if (is.null(keep)) NA_character_ else keep$log_file,
        if (is.null(keep)) NA_integer_ else keep$pid,
        if (is.null(keep)) NA_character_ else keep$preset_name,
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

#' Slug a clade name into a ref_id no existing group is already using
#'
#' The derived slug for a clade can coincide with a group that already exists
#' (or with the very row being rewritten), which would silently merge or
#' overwrite it; the taxid suffix keeps the key its own.
#'
#' @param con Open project DB connection.
#' @param name,taxid The clade.
#' @param avoid Extra keys the result must not equal.
#' @noRd
unique_split_ref_id <- function(con, name, taxid, avoid = character(0)) {
  slug <- ref_id_from_name(name)
  existing <- tryCatch(DBI::dbGetQuery(con, "SELECT ref_id FROM reference")$ref_id,
    error = function(e) character(0))
  if (!(slug %in% existing) && !(slug %in% avoid)) return(slug)
  paste(slug, taxid, sep = "_")
}

#' Point a set of samples at a chosen clade, creating the group if needed
#'
#' The ref_id is the clade's identity everywhere downstream (output file names,
#' the folder manifest, the preset name), so a group whose clade actually
#' changes is re-keyed on the new clade and its samples follow.
#'
#' @param con Open project DB connection.
#' @param ref_id Group key, or NULL to derive one from `clade_name`.
#' @param clade_taxid,clade_name,clade_rank,n_records The chosen NCBI node.
#' @param member_ids Sample IDs to move into the group.
#' @return (invisibly) the group key, which is NOT `ref_id` when the clade changed.
#' @noRd
set_reference_clade <- function(con, ref_id = NULL, clade_taxid, clade_name, clade_rank,
                               n_records = NA_integer_, member_ids = character(0)) {
  key <- ref_id %||% ref_id_from_name(clade_name)
  cur <- DBI::dbGetQuery(con, "SELECT clade_taxid FROM reference WHERE ref_id = ? LIMIT 1",
    params = list(key))
  exists_row <- nrow(cur) > 0
  now <- as.integer(Sys.time())
  if (exists_row) {
    old_taxid <- suppressWarnings(as.integer(cur$clade_taxid[1]))
    if (is.na(old_taxid) || old_taxid != as.integer(clade_taxid)) {
      # a clade change invalidates whatever was built for the old one, and the
      # key has to move too or the new clade's database would be written under
      # the old clade's name
      new_key <- unique_split_ref_id(con, clade_name, as.integer(clade_taxid))
      DBI::dbExecute(con,
        "UPDATE reference SET ref_id = ?, clade_taxid = ?, clade_name = ?, clade_rank = ?,
           query = ?, n_records = ?, reference_switch = 1, source = NULL, seeds_db = NULL,
           labels_db = NULL, preset_name = NULL, fail_reason = NULL, time_stamp = ?
         WHERE ref_id = ?",
        params = list(new_key, as.integer(clade_taxid), clade_name, clade_rank,
          clade_query(clade_name), as.integer(n_records), now, key))
      DBI::dbExecute(con, "UPDATE samples SET ref_id = ? WHERE ref_id = ?",
        params = list(new_key, key))
      key <- new_key
    }
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

#' Name for a group's preset, never taking over one the user already made
#'
#' A preset already named after the clade may belong to the user and be used by
#' samples outside the group, so it is left alone and a distinct name is used.
#' The choice is stored on the reference row so re-wiring reuses it.
#' @noRd
reference_preset_name <- function(con, ref_id) {
  stored <- tryCatch(DBI::dbGetQuery(con,
    "SELECT preset_name FROM reference WHERE ref_id = ?", params = list(ref_id))$preset_name,
    error = function(e) character(0))
  if (length(stored) > 0 && !is.na(stored[1]) && nzchar(stored[1])) return(stored[1])
  taken <- tryCatch(DBI::dbGetQuery(con, "SELECT assemble_opts FROM assemble_opts")$assemble_opts,
    error = function(e) character(0))
  if (!(ref_id %in% taken)) return(ref_id)
  cand <- paste0(ref_id, "_ref")
  i <- 2L
  while (cand %in% taken) {
    cand <- paste0(ref_id, "_ref", i)
    i <- i + 1L
  }
  cand
}

#' Give a group's samples an assemble_opts preset pointing at its database
#'
#' The preset is named after the group (see [reference_preset_name()]) and
#' copies every other setting from the preset the group's samples are on, so
#' switching database changes the database and nothing else. Locked samples keep
#' their preset: their stored results were produced with the old one.
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
  preset <- reference_preset_name(con, ref_id)
  row <- base[1, , drop = FALSE]
  row$assemble_opts <- preset
  row$seeds_db <- seeds
  row$labels_db <- labels
  DBI::dbExecute(con, "DELETE FROM assemble_opts WHERE assemble_opts = ?", params = list(preset))
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
      params = c(list(preset, as.integer(Sys.time())), as.list(ids)))
  }
  DBI::dbExecute(con,
    "UPDATE reference SET reference_switch = 2, source = ?, seeds_db = ?, labels_db = ?,
       preset_name = ?, fail_reason = NULL, pid = NULL, time_stamp = ? WHERE ref_id = ?",
    params = list(source, seeds, labels, preset, as.integer(Sys.time()), ref_id))
  invisible(length(ids))
}

#' Point a group at seed/label FASTAs that already exist
#' @return (invisibly) the number of samples moved.
#' @noRd
use_existing_reference <- function(con, ref_id, seeds, labels) {
  seeds <- path.expand(seeds); labels <- path.expand(labels)
  for (f in c(seeds, labels)) {
    if (!file.exists(f)) stop("Reference file not found: ", f, call. = FALSE)
    if (!usable_fasta(f)) {
      stop("Reference file is empty or not a FASTA: ", f, call. = FALSE)
    }
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

#' Spawn a Nextflow process with piped output
#'
#' Indirected so a test can capture the command and environment without any
#' process being started.
#' @noRd
nf_spawn <- function(cmd, args, wd, env) {
  if (!nzchar(Sys.which(cmd))) stop("Nextflow not found on PATH.", call. = FALSE)
  processx::process$new(cmd, args, stdout = "|", stderr = "|", wd = wd, env = env)
}

#' Is a pid still running?
#' @noRd
process_alive <- function(pid) {
  if (is.na(pid)) return(FALSE)
  isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(as.integer(pid))), error = function(e) FALSE))
}

#' Where a project's BUILD clade manifest lives
#'
#' Fixed, so the command shown in the modal is the command that runs.
#' @noRd
reference_manifest_path <- function(project) {
  file.path(project, ".reference", "clades_build.tsv")
}

#' Reference rows for a set of group keys, erroring on an unknown key
#' @noRd
reference_rows <- function(con, ref_ids) {
  ref_ids <- unique(as.character(ref_ids))
  ph <- paste(rep("?", length(ref_ids)), collapse = ",")
  rows <- DBI::dbGetQuery(con,
    sprintf("SELECT * FROM reference WHERE ref_id IN (%s)", ph), params = as.list(ref_ids))
  gone <- setdiff(ref_ids, rows$ref_id)
  if (length(gone) > 0) stop("No reference group '", gone[1], "'.", call. = FALSE)
  rows
}

#' Is a built FASTA actually usable?
#'
#' Existence is not enough: a build that dies partway leaves the empty files it
#' created behind, and treating those as a database wires a preset that recruits
#' no reads and reports success while doing it.
#' @noRd
usable_fasta <- function(f) {
  if (!file.exists(f)) return(FALSE)
  sz <- suppressWarnings(file.info(f)$size)
  if (is.na(sz) || sz <= 0) return(FALSE)
  first <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) character(0))
  length(first) > 0 && startsWith(trimws(first[1]), ">")
}

#' Group keys whose seed/label FASTAs are not both usable yet
#' @noRd
reference_missing_dbs <- function(ref_dir, ref_ids) {
  ref_ids <- unique(as.character(ref_ids))
  ref_ids[!vapply(ref_ids, function(r) {
    p <- reference_db_paths(ref_dir, r)
    usable_fasta(p$seeds) && usable_fasta(p$labels)
  }, logical(1))]
}

#' Write the BUILD manifest for a set of reference rows
#' @noRd
write_build_manifest <- function(project, rows) {
  write_clade_manifest(reference_manifest_path(project),
    data.frame(ref_id = rows$ref_id, taxid = as.integer(rows$clade_taxid),
      clade_name = rows$clade_name, query = rows$query, stringsAsFactors = FALSE))
}

#' Append run output to the log a build's rows point at
#'
#' The in-app run pipes Nextflow's output into R for the progress pane, so the
#' log the rows record only exists if the stream is teed into it.
#' @noRd
append_run_log <- function(log, lines) {
  if (is.null(log) || length(log) != 1 || is.na(log) || !nzchar(log)) return(invisible(NULL))
  if (length(lines) == 0) return(invisible(NULL))
  cat(paste0(cli::ansi_strip(lines), "\n"), file = log, sep = "", append = TRUE)
  invisible(NULL)
}

#' Start one background build run covering a set of clade groups
#'
#' One Nextflow run handles every group, so all launched rows share the pid, the
#' log, and the launch time; [poll_reference_builds()] still judges each row on
#' its own FASTAs.
#'
#' @param con Open project DB connection.
#' @param project Project directory (the log lives in `<project>/.reference`).
#' @param ref_ids Group keys.
#' @param api_key,refseq,resume Passed to [nf_reference_command()].
#' @param overwrite Rebuild groups whose database already exists.
#' @param launcher function(cmd, args, log) -> pid.
#' @return (invisibly) the log path, or NULL when every group was reused.
#' @noRd
start_reference_build <- function(con, project, ref_ids, api_key = NULL, refseq = FALSE,
                                  resume = TRUE, overwrite = FALSE,
                                  launcher = processx_launcher) {
  ref_ids <- unique(as.character(ref_ids))
  if (length(ref_ids) == 0) return(invisible(NULL))
  rows <- reference_rows(con, ref_ids)
  cfg <- ref_config(con)

  # Without overwrite a group whose database is already on disk would be rebuilt
  # for nothing; wire what is there and keep it out of the run.
  keep <- rep(TRUE, nrow(rows))
  for (i in seq_len(nrow(rows))) {
    have <- reference_db_paths(cfg$ref_dir, rows$ref_id[i])
    if (!isTRUE(overwrite) && usable_fasta(have$seeds) && usable_fasta(have$labels)) {
      wire_reference_db(con, rows$ref_id[i], have$seeds, have$labels, source = "existing")
      keep[i] <- FALSE
    }
  }
  rows <- rows[keep, , drop = FALSE]
  if (nrow(rows) == 0) return(invisible(NULL))

  rundir <- file.path(project, ".reference")
  dir.create(rundir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(rundir, paste0("build_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  unlink(log)
  file.create(log)

  man <- write_build_manifest(project, rows)
  args <- nf_reference_command(project, "BUILD", man, cfg$ref_dir, cfg$ref_dir,
    api_key = api_key %||% "", refseq = isTRUE(refseq), resume = isTRUE(resume))

  # captured before the launcher runs, so any FASTA this build writes has an
  # mtime >= launch_time; anything already on disk with an older mtime is a
  # stale leftover from a previous build, not evidence this one succeeded.
  launch_time <- as.integer(Sys.time())
  pid <- launcher("nextflow", args, log)

  ph <- paste(rep("?", nrow(rows)), collapse = ",")
  DBI::dbExecute(con,
    sprintf("UPDATE reference SET reference_switch = 4, pid = ?, log_file = ?, fail_reason = NULL,
       time_stamp = ? WHERE ref_id IN (%s)", ph),
    params = c(list(as.integer(pid), log, launch_time), as.list(rows$ref_id)))
  invisible(log)
}

#' Last non-empty lines of a build log, for the failure message
#' @noRd
log_tail <- function(log, n = 3) {
  if (is.null(log) || is.na(log) || !file.exists(log)) return("")
  x <- tryCatch(readLines(log, warn = FALSE), error = function(e) character(0))
  x <- trimws(cli::ansi_strip(x))
  x <- x[nzchar(x)]
  # an in-app run tees the redrawing progress board into the log; its rows would
  # otherwise crowd the error that actually ended the run out of the tail
  board <- grepl(nf_proc_rx, x, perl = TRUE) | grepl("^executor\\s*>", x)
  if (any(!board)) x <- x[!board]
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
      if (!usable_fasta(f)) return(FALSE)
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
