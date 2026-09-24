# Per-sample pipeline steps, called by the WF1 Nextflow R-script processes
# (bin/rp_score.R, bin/rp_reconcile.R). Each takes an open DBI
# connection + a sample ID, does one stage, and writes to the database. Kept
# separate from the project-wide drivers so a Nextflow process can run one stage
# for one sample.

#' Score one sample's candidate paths and write path_scores (+ auto_select)
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param out_dir Pipeline output dir (assembly fastas + gfa live under
#'   `<out_dir>/<id>/assemble`).
#' @param taxdump_dir Optional pinned taxdump for hit-lineage resolution.
#' @return (invisibly) the scored data.frame.
#' @export
score_one <- function(con, id, out_dir, taxdump_dir = NULL) {
  td <- if (have_taxdump(taxdump_dir)) taxdump_dir else NULL
  seeds_db <- DBI::dbGetQuery(con,
    "SELECT ao.seeds_db FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(id))$seeds_db[1]
  blast_db <- make_blast_db(local_ref(seeds_db))

  adir <- file.path(out_dir, id, "assemble")
  fastas <- Sys.glob(file.path(adir, paste0(id, "_assembly_*.fasta")))
  fastas <- fastas[!grepl("_assembly_0\\.fasta$", fastas)]
  if (length(fastas) == 0) {
    return(invisible(NULL))
  }
  gfa <- Sys.glob(file.path(adir, "*.gfa"))
  tax <- DBI::dbGetQuery(con,
    "SELECT lineage, lineage_taxids FROM samples WHERE ID = ?", params = list(id))

  scored <- score_assembly(
    path_fastas = fastas,
    gfa_path = if (length(gfa) > 0) gfa[1] else NA_character_,
    sample_lineage = tax$lineage[1] %||% NA_character_,
    sample_taxids = tax$lineage_taxids[1] %||% NA_character_,
    blast_db = blast_db, taxdump_dir = td, with_staxids = FALSE)
  write_path_scores(con, id, scored)
  tryCatch(blast_scaffolds_local(con, id, blast_db), error = function(e) NULL)
  invisible(scored)
}

#' Per-scaffold BLAST (local seed DB) -> assemblies.blast_* best hit
#'
#' Gives each assembly unit its own top hit for review (accession / species /
#' percent identity), so per-scaffold hits are visible without a remote BLAST.
#' Best-effort: a no-op when blastn / the seed DB is unavailable.
#' @param con Open DBI connection. @param id Sample ID. @param blast_db BLAST db prefix.
#' @noRd
blast_scaffolds_local <- function(con, id, blast_db) {
  if (is.null(blast_db) || is.na(blast_db) || !nzchar(blast_db)) {
    return(invisible())
  }
  scaf <- DBI::dbGetQuery(con,
    "SELECT path, scaffold, sequence FROM assemblies WHERE ID = ? ORDER BY path, scaffold", params = list(id))
  if (nrow(scaf) == 0) {
    return(invisible())
  }
  fa <- tempfile(fileext = ".fasta"); on.exit(unlink(fa), add = TRUE)
  writeLines(paste0(">", scaf$path, ".", scaf$scaffold, "\n", scaf$sequence), fa)
  bh <- tryCatch(blast_best_hits(fa, blast_db), error = function(e) data.frame())
  if (nrow(bh) == 0) {
    return(invisible())
  }
  for (i in seq_len(nrow(bh))) {
    ps <- strsplit(bh$qseqid[i], "\\.")[[1]]
    if (length(ps) != 2) next
    DBI::dbExecute(con,
      # blast_source records which search named the species: the remote path
      # writes 'core_nt', so a hit left unlabelled here reads as a GenBank call.
      "UPDATE assemblies SET blast_accession = ?, blast_species = ?, blast_pident = ?,
         blast_source = 'seeds'
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(sub("_.*", "", bh$sseqid[i]), sseqid_species(bh$sseqid[i]),
        round(bh$pident[i], 1), id, as.integer(ps[1]), as.integer(ps[2])))
  }
  invisible()
}

#' Reconcile barrnap + ITSx outputs for one sample; write annotations
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param gff barrnap GFF file (or NULL).
#' @param positions ITSx positions file (or NULL).
#' @return (invisibly) the feature data.frame.
#' @export
reconcile_one <- function(con, id, path = 1L, scaffold = 1L, gff = NULL, positions = NULL) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  bar <- if (!is.null(gff) && file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
  pos_lines <- if (!is.null(positions) && file.exists(positions)) readLines(positions, warn = FALSE) else character(0)
  its <- parse_itsx(pos_lines)
  feat <- reconcile_features(bar, its)
  chimeric <- as.integer(itsx_is_chimeric(pos_lines))
  ts <- as.integer(Sys.time())

  # one transaction: a failed append after a committed delete would leave the
  # unit with no annotations at all
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    if (nrow(feat) > 0) {
      df <- data.frame(ID = id, path = path, scaffold = scaffold, element = feat$element, type = feat$type,
        start = feat$start, stop = feat$stop, length = feat$length, strand = feat$strand,
        model_coverage = feat$model_cov %||% NA_real_,
        tool = feat$tool, edited = 0L, notes = NA_character_, time_stamp = ts, stringsAsFactors = FALSE)
      DBI::dbAppendTable(con, "annotations", df)
    }
  })
  missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
  DBI::dbExecute(con,
    "UPDATE annotate SET n_features = ?, missing = ?, chimeric = ?, annotate_switch = 2, time_stamp = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(nrow(feat), paste(missing, collapse = ","), chimeric, ts, id, path, scaffold))
  invisible(feat)
}

#' Filesystems whose locking cannot support SQLite WAL
#' @noRd
NETWORK_FS <- c("nfs", "nfs4", "lustre", "gpfs", "beegfs", "panfs",
                "cifs", "smb2", "afs", "fuseblk")

#' Classify the filesystem holding `path` as "local", "network", or "unknown"
#'
#' Linux only (`stat -f`); everything else is "unknown" and treated as local.
#' @noRd
filesystem_type <- function(path) {
  if (!identical(unname(Sys.info()[["sysname"]]), "Linux")) return("unknown")
  fs <- tryCatch(
    suppressWarnings(system2("stat", c("-f", "-c", "%T", shQuote(path)),
      stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(fs) == 0 || is.na(fs[1]) || !nzchar(fs[1])) return("unknown")
  if (tolower(fs[1]) %in% NETWORK_FS) "network" else "local"
}

#' Cache of filesystem type by directory, so the app does not fork `stat` on
#' every query (it opens a fresh connection per query).
#' @noRd
.fs_cache <- new.env(parent = emptyenv())

#' Journal mode for a project DB
#'
#' WAL coordinates readers/writers through an mmap'd `-shm` file and requires
#' every process to be on one host, so it is unsupported on the network
#' filesystems HPC projects live on. Use the rollback journal there.
#' Override with `RIBOPILOT_SQLITE_JOURNAL` ("wal" or "delete").
#' @noRd
sqlite_journal_mode <- function(db) {
  ov <- tolower(trimws(Sys.getenv("RIBOPILOT_SQLITE_JOURNAL", "")))
  if (nzchar(ov)) return(ov)
  dir <- dirname(db)
  fs <- .fs_cache[[dir]]
  if (is.null(fs)) {
    fs <- filesystem_type(dir)
    assign(dir, fs, envir = .fs_cache)
  }
  if (identical(fs, "network")) "delete" else "wal"
}

#' Open a project DB connection with a busy timeout and a filesystem-appropriate
#' journal mode
#'
#' The busy timeout is set first: `journal_mode` takes an exclusive lock, so
#' setting it before the timeout would run that one statement with no timeout.
#' @param db Path to the project SQLite file.
#' @export
open_project_db <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
  DBI::dbExecute(con, sprintf("PRAGMA journal_mode=%s;", sqlite_journal_mode(db)))
  con
}

#' Open a SQLite database read-only
#'
#' The shared reference store may live on a path the running user cannot write,
#' which is the recommended cluster deployment. `journal_mode` takes an
#' exclusive lock and fails there, so this deliberately sets only busy_timeout.
#'
#' @param path Path to an existing SQLite file.
#' @return An open read-only DBI connection.
#' @noRd
open_sqlite_ro <- function(path) {
  if (!file.exists(path)) {
    stop("no such database: ", path, call. = FALSE)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
  con
}

#' Was this sample's last assembly stopped by the max_paths / max_scaffolds gate
#'
#' Recognised from the state and the note the gate writes, so it stays true
#' across a restart and needs no extra column.
#' @noRd
gate_failed <- function(con, id) {
  r <- tryCatch(DBI::dbGetQuery(con,
    "SELECT assemble_switch, COALESCE(assemble_notes,'') AS notes FROM assemble WHERE ID = ?",
    params = list(id)), error = function(e) data.frame())
  if (nrow(r) == 0) return(FALSE)
  isTRUE(r$assemble_switch[1] == 3L) && grepl("max_paths|max_scaffolds", r$notes[1])
}

#' Settle a sample's state against the work still queued for it
#'
#' The sample state is written by the pipeline as the sample moves through WF1:
#' loading an assembly opens it at 4 (in progress) and triage, the last
#' per-sample step, closes it at 2 (done). This recomputes it from what is
#' actually pending, so a sample also closes when the remaining work is taken
#' away rather than finished -- units held at 0 or on ignored scaffolds are a
#' deliberate skip, not unfinished work.
#'
#' Only 2 and 4 are rewritten. A sample the user has held (0), queued to re-run
#' (1), or whose assembly failed (3) is saying something this cannot improve on.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return (invisibly) the state now stored.
#' @noRd
refresh_sample_state <- function(con, id) {
  pending <- tryCatch(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM annotate an
       LEFT JOIN assemblies asm
         ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
      WHERE an.ID = ? AND an.annotate_switch IN (1,4) AND COALESCE(asm.ignore,0) = 0",
    params = list(id))$n[1], error = function(e) NA_integer_)
  if (is.na(pending)) return(invisible(NA_integer_))
  st <- if (pending > 0) 4L else 2L
  DBI::dbExecute(con,
    "UPDATE assemble SET assemble_switch = ? WHERE ID = ? AND assemble_switch IN (2,4)",
    params = list(st, id))
  invisible(st)
}

#' Per-unit tables whose rows describe one assembled sequence
#' @noRd
UNIT_EVIDENCE_TABLES <- c("annotations", "coverage", "screen", "chimera", "chimera_windows",
  "blast_hits", "blast_screen", "blast_searched", "blast_ref_candidates",
  "blast_ref_alignment", "assembly_backup")

#' Forget everything known about a unit whose sequence was replaced
#'
#' A re-assembly can give a surviving (path, scaffold) key a different sequence.
#' Its gene coordinates, BLAST hits, screen verdict, chimera call, and depth then
#' describe another molecule, so they are deleted, and the unit's lock and
#' curator marks are cleared so it is annotated and reviewed again.
#' @noRd
reset_unit_evidence <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  for (tbl in intersect(UNIT_EVIDENCE_TABLES, DBI::dbListTables(con))) {
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
      params = p)
  }
  DBI::dbExecute(con,
    "UPDATE annotate SET annotate_switch = 1, annotate_lock = 0, reviewed = 'no',
       partial = 'no', problematic = 'no', export_group = NULL
     WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
}

#' Load one sample's published assembly fastas into the database
#'
#' Parses `out/<ID>/assemble/<ID>_assembly_*.fasta` (headers
#' `>ID.path.scaffold topology`), writes per-scaffold rows into `assemblies`,
#' and updates the per-sample `assemble` summary.
#'
#' @param con Open DBI connection.
#' @param out_dir The pipeline `out/` directory.
#' @param ID Sample ID.
#' @return (invisibly) the number of scaffolds loaded.
#' @export
load_assemblies <- function(con, out_dir, ID) {
  fas <- Sys.glob(file.path(out_dir, ID, "assemble", paste0(ID, "_assembly_*.fasta")))
  fas <- fas[!grepl("_assembly_0\\.fasta$", fas)] # skip the "no assembly" sentinel
  ts <- as.integer(Sys.time())

  if (length(fas) == 0) {
    DBI::dbExecute(con,
      "UPDATE assemble SET paths=NULL, scaffolds=NULL, length=NULL, topology=NULL,
       assemble_switch=3, assemble_notes='failed assembly', time_stamp=? WHERE ID=?",
      params = list(ts, ID))
    return(invisible(0L))
  }

  rows <- list()
  n_bad <- 0L
  for (fa in fas) {
    seqs <- Biostrings::readDNAStringSet(fa)
    for (i in seq_along(seqs)) {
      hdr <- names(seqs)[i]
      parts <- strsplit(hdr, "\\s+")[[1]]
      # path and scaffold are the LAST two dot-separated fields: the ID may
      # itself contain dots, and a positional read would take them for the keys
      idpath <- strsplit(parts[1], "\\.")[[1]]
      np <- length(idpath)
      pth <- if (np >= 3) suppressWarnings(as.integer(idpath[np - 1])) else NA_integer_
      scf <- if (np >= 3) suppressWarnings(as.integer(idpath[np])) else NA_integer_
      if (is.na(pth) || is.na(scf)) {
        n_bad <- n_bad + 1L
        next
      }
      topo <- if (length(parts) > 1) parts[2] else NA_character_
      # whatever the source called this contig, kept so a supplied assembly can
      # be cross-referenced back to the user's own file
      src <- if (length(parts) > 2) paste(parts[-(1:2)], collapse = " ") else NA_character_
      s <- as.character(seqs[[i]])
      rows[[length(rows) + 1]] <- data.frame(
        ID = ID, path = pth, scaffold = scf,
        topology = topo, length = nchar(s), length_raw = nchar(s),
        sequence = s, source_contig = src, ignore = 0L, edited = 0L, time_stamp = ts,
        stringsAsFactors = FALSE
      )
    }
  }
  if (n_bad > 0) {
    warning(sprintf("%s: skipped %d assembly record(s) with an unparsable header", ID, n_bad),
      call. = FALSE)
  }
  if (length(rows) == 0) {
    DBI::dbExecute(con,
      "UPDATE assemble SET paths=NULL, scaffolds=NULL, length=NULL, topology=NULL,
       assemble_switch=3, assemble_notes='failed assembly', time_stamp=? WHERE ID=?",
      params = list(ts, ID))
    return(invisible(0L))
  }
  df <- do.call(rbind, rows)

  # max_paths / max_scaffolds are GATES, not caps: a sample over either one is
  # failed here and nothing downstream runs for it. max_paths used to silently
  # keep its top N paths by length, and max_scaffolds was read by nothing at all,
  # so a 47-scaffold assembly went on to cost one serial remote BLAST per
  # scaffold. Checked before the per-scaffold orientation pass, which shells out
  # to an anchor caller for every scaffold. NA or <= 0 means no limit.
  lim <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ao.max_paths, ao.max_scaffolds FROM assemble a
       JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(ID)), error = function(e) data.frame())
  # queried separately so a database without the column loses only the filter,
  # not the two gates as well
  ml <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ao.min_assembly_length FROM assemble a
       JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(ID)), error = function(e) data.frame())
  over <- function(n, limit) {
    !is.null(limit) && length(limit) && !is.na(limit) && limit > 0 && n > limit
  }

  # min_assembly_length is a FILTER, not a gate: a short contig is noise, not a
  # reason to fail the sample. Parked (as a triage drop, below) rather than
  # discarded, so the sequence is still there if the curator lowers the limit,
  # and counted out of the gate here so junk fragments cannot trip max_scaffolds.
  minlen <- if (nrow(ml) > 0) ml$min_assembly_length[1] else NA_integer_
  short <- if (!is.na(minlen) && minlen > 0) df$length < minlen else rep(FALSE, nrow(df))
  short_note <- if (any(short)) sprintf("%d scaffold(s) under min_assembly_length (%d bp)",
    sum(short), as.integer(minlen)) else NULL

  kept <- df[!short, , drop = FALSE]
  n_paths_raw <- length(unique(kept$path))
  # per path, not summed: 3 paths of 2 scaffolds is not a 6-scaffold assembly
  max_per_path_raw <- if (nrow(kept) > 0) max(table(kept$path)) else 0L
  gate <- c(
    # nothing left to work with is a failure, not a silently empty sample
    if (nrow(kept) == 0) short_note else NULL,
    if (over(n_paths_raw, lim$max_paths[1]))
      sprintf("%d paths exceeds max_paths (%d)", n_paths_raw, lim$max_paths[1]),
    if (over(max_per_path_raw, lim$max_scaffolds[1]))
      sprintf("%d scaffolds in a path exceeds max_scaffolds (%d)",
        max_per_path_raw, lim$max_scaffolds[1]))

  if (length(gate) > 0) {
    # Store the sequences anyway: raising the limit and re-running should not
    # mean re-assembling, and the curator needs to see what tripped the gate.
    # ignore = 1 is what actually stops the work -- every downstream selection
    # (remote BLAST, coverage, the annotate unit query, export) filters on it.
    df$ignore <- 1L
    # one transaction: a failed append after a committed delete would destroy the
    # sample's assemblies with nothing to put back
    DBI::dbWithTransaction(con, {
      DBI::dbExecute(con, "DELETE FROM assemblies WHERE ID = ?", params = list(ID))
      DBI::dbAppendTable(con, "assemblies", df)
      DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 0 WHERE ID = ?", params = list(ID))
      DBI::dbExecute(con,
        "UPDATE assemble SET paths=?, scaffolds=?, length=?, topology=?,
         assemble_switch=3, assemble_notes=?, time_stamp=? WHERE ID=?",
        params = list(length(unique(df$path)), nrow(df),
          paste(sort(unique(df$length), decreasing = TRUE), collapse = ";"),
          paste(sort(unique(stats::na.omit(df$topology))), collapse = ";"),
          paste(gate, collapse = "; "), ts, ID))
    })
    return(invisible(nrow(df)))
  }

  # Orient each kept scaffold to the sense (5'->3') strand before storing, so
  # per-scaffold annotation/coverage/export all see plus-strand features. Stitching
  # used to normalize orientation once per sample; with independent per-scaffold
  # annotation it must happen per scaffold. No-op (unchanged) when no rRNA anchor
  # caller is on PATH or the scaffold has no anchors; length is preserved by RC.
  for (i in seq_len(nrow(df))) {
    or <- tryCatch(orient_rrna_sense(df$sequence[i]), error = function(e) NULL)
    if (!is.null(or) && isTRUE(or$flipped)) df$sequence[i] <- or$seq
  }

  # Carry the curator's keep/drop decisions across a re-assembly. The rows are
  # replaced wholesale below, and opts_apply re-queues a sample whenever its
  # assemble_opts change, so without this every ignore flag is silently reset.
  prev <- tryCatch(DBI::dbGetQuery(con,
    "SELECT path, scaffold, COALESCE(ignore,0) AS ignore FROM assemblies WHERE ID = ?",
    params = list(ID)), error = function(e) data.frame())
  # ...unless the last run was a gate failure, which ignored every scaffold
  # wholesale. Those are not curator decisions, and carrying them would leave the
  # sample empty after raising the limit and re-running -- the exact recovery the
  # gate is supposed to allow.
  if (isTRUE(gate_failed(con, ID))) prev <- data.frame()
  if (nrow(prev) > 0) {
    i <- match(paste(df$path, df$scaffold), paste(prev$path, prev$scaffold))
    df$ignore <- ifelse(is.na(i), 0L, as.integer(prev$ignore[i]))
  }
  old <- DBI::dbGetQuery(con, "SELECT path, scaffold, sequence FROM assemblies WHERE ID = ?",
    params = list(ID))
  j <- match(paste(df$path, df$scaffold), paste(old$path, old$scaffold))
  changed <- which(!is.na(j) & (is.na(old$sequence[j]) | old$sequence[j] != df$sequence))

  # replace this sample's assemblies, in one transaction: a failed append after a
  # committed delete would destroy them with nothing to put back
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM assemblies WHERE ID = ?", params = list(ID))
    DBI::dbAppendTable(con, "assemblies", df)
    for (i in changed) reset_unit_evidence(con, ID, df$path[i], df$scaffold[i])
  })

  # Per-sample summary, described over the scaffolds that passed the length
  # filter. One basis for every column, or the row contradicts itself: a path of
  # nothing but junk fragments would otherwise inflate the counts while dropping
  # out of the path total, and the length list would advertise contigs the sample
  # is not using.
  n_paths <- n_paths_raw
  max_per_path <- max_per_path_raw      # drives the "disconnected contigs" note
  n_scaffolds <- nrow(kept)             # total units, not the biggest path's count
  length_str <- paste(sort(unique(kept$length), decreasing = TRUE), collapse = ";")
  topo_str <- paste(sort(unique(stats::na.omit(kept$topology))), collapse = ";")
  notes <- paste(c(
    if (n_paths > 1) "Unable to resolve single assembly from reads"
    else if (max_per_path > 1) "Output contains disconnected contigs" else NULL,
    short_note),
    collapse = "; ")
  DBI::dbExecute(con,
    "UPDATE assemble SET paths=?, scaffolds=?, length=?, topology=?,
     assemble_switch=4, assemble_notes=?, time_stamp=? WHERE ID=?",
    params = list(n_paths, n_scaffolds, length_str, topo_str, notes, ts, ID))

  # Seed one annotate unit per non-ignored scaffold, inheriting the sample's opts
  # from any existing annotate row (the init placeholder or a prior run), then
  # drop annotate units whose scaffold no longer exists after this (re)assembly.
  DBI::dbExecute(con,
    "INSERT OR IGNORE INTO annotate
       (ID, path, scaffold, annotate_switch, annotate_lock, annotate_opts, curate_opts, reviewed)
     SELECT a.ID, a.path, a.scaffold, 1, 0,
            COALESCE(o.annotate_opts,'default'), COALESCE(o.curate_opts,'default'), 'no'
     FROM assemblies a
     LEFT JOIN (SELECT annotate_opts, curate_opts FROM annotate WHERE ID = ? LIMIT 1) o
     WHERE a.ID = ? AND COALESCE(a.ignore,0) = 0",
    params = list(ID, ID))
  # Drop annotate units (and their orphaned annotations/coverage) whose scaffold
  # no longer exists after this (re)assembly. Deliberately NOT gated on
  # assemblies.ignore: a curator-dropped scaffold still exists, and treating it
  # as an orphan here would delete on the next WF1 run exactly the annotation
  # that reconcile_ignore_unit now takes care to preserve.
  for (tbl in c("annotate", "annotations", "coverage")) {
    DBI::dbExecute(con, sprintf(
      "DELETE FROM %s WHERE ID = ? AND NOT EXISTS
         (SELECT 1 FROM assemblies a WHERE a.ID = %s.ID AND a.path = %s.path
            AND a.scaffold = %s.scaffold)",
      tbl, tbl, tbl, tbl), params = list(ID))
  }
  # The length filter is recorded as a triage decision, not by writing
  # assemblies.ignore here. Triage owns that column: it derives it from the
  # verdict, and "undecided" means ignore = 0. Setting the flag directly meant
  # triage overturned this at the end of the same run -- it saw a scaffold with
  # no annotation (which is exactly what parking causes), called that undecided,
  # and un-parked it. Going through set_unit_triage() also puts the reason in
  # front of the curator and holds the unit's annotate row.
  if (any(short)) {
    reason <- sprintf("shorter than min_assembly_length (%d bp)", as.integer(minlen))
    for (i in which(short)) {
      set_unit_triage(con, ID, df$path[i], df$scaffold[i], "drop",
        reason = reason, source = "auto", rule = "too_short")
    }
  }

  # A (re)assembly invalidates any pre-trim undo snapshots for this sample.
  tryCatch(DBI::dbExecute(con, "DELETE FROM assembly_backup WHERE ID = ?", params = list(ID)),
    error = function(e) NULL)
  invisible(nrow(df))
}

#' Trim one assembly unit (ID,path,scaffold) to [from, to], shifting its
#' annotations and coverage to match. Features that fall entirely outside the
#' kept range are dropped; overlapping features are clamped. Used by the Annotate
#' tab's "Trim assembly" control.
#'
#' @param con Open DBI connection.
#' @param id,path,scaffold Assembly unit.
#' @param from,to 1-based inclusive kept range in the unit's coordinates.
#' @return (invisibly) the new unit length.
#' @noRd
trim_assembly <- function(con, id, path = 1L, scaffold = 1L, from, to) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  s <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(s) == 0 || is.na(s$sequence[1]) || !nzchar(s$sequence[1])) {
    stop("No sequence for unit ", id, ".", path, ".", scaffold, call. = FALSE)
  }
  n <- nchar(s$sequence[1])
  from <- max(1L, as.integer(from)); to <- min(n, as.integer(to))
  if (from > to) stop("Invalid trim range.", call. = FALSE)

  # Snapshot the original (pre-trim) unit state once, so [restore_assembly()] can
  # undo this. Only the first trim writes it, so restore always returns the full unit.
  ensure_assembly_backup(con)
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assembly_backup WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))$n == 0) {
    ann0 <- DBI::dbGetQuery(con, "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    cov0 <- DBI::dbGetQuery(con, "SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    DBI::dbExecute(con,
      "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length, annotations, coverage, time_stamp)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(id, path, scaffold, s$sequence[1], n,
        as.character(jsonlite::toJSON(ann0, dataframe = "rows", na = "null")),
        if (nrow(cov0)) cov0$depths[1] else NA_character_, as.integer(Sys.time())))
  }
  new_seq <- substr(s$sequence[1], from, to); new_len <- nchar(new_seq)

  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(new_seq, new_len, id, path, scaffold))

  ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  for (i in seq_len(nrow(ann))) {
    ns <- ann$start[i] - from + 1; nt <- ann$stop[i] - from + 1
    if (nt < 1 || ns > new_len) {
      DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = list(id, path, scaffold, ann$element[i]))
      next
    }
    ns <- max(1L, ns); nt <- min(new_len, nt)
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = list(ns, nt, nt - ns + 1, id, path, scaffold, ann$element[i]))
  }

  cov <- DBI::dbGetQuery(con, "SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  dep <- if (nrow(cov) > 0) parse_depths(cov$depths) else numeric(0)
  # `from` is clamped to the sequence, which can outrun a shorter depth track
  # (truncated `samtools depth -a`, or a re-assembly coverage never re-ran on)
  if (length(dep) > 0 && from <= length(dep)) {
    write_coverage_summary(con, id, path, scaffold, dep[seq(from, min(length(dep), to))])
  }
  DBI::dbExecute(con, "UPDATE coverage SET clips = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  chimera_shift(con, id, path, scaffold, from - 1L, new_len)
  invisible(new_len)
}

#' Rewrite one unit's coverage row from a depth track
#'
#' median_depth is what the scaffold table and triage read, so it has to move
#' with the track, not just mean/min/max.
#' @param con Open DBI connection. @param id,path,scaffold Assembly unit.
#' @param dep Numeric depths in position order.
#' @noRd
write_coverage_summary <- function(con, id, path, scaffold, dep) {
  if (length(dep) == 0) return(invisible(FALSE))
  DBI::dbExecute(con,
    "UPDATE coverage SET depths = ?, mean_depth = ?, median_depth = ?, min_depth = ?, max_depth = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(paste(dep, collapse = ","), round(mean(dep), 2),
      round(stats::median(dep), 2), as.integer(min(dep)), as.integer(max(dep)),
      id, path, scaffold))
  invisible(TRUE)
}

#' Ensure the assembly_backup table exists (idempotent).
#' @noRd
ensure_assembly_backup <- function(con) {
  DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)
}

#' Reconcile a unit's annotate row after its assemblies.ignore flag changed
#'
#' Ignoring a scaffold drops it from annotation/export, so its annotate unit (and
#' any annotations/coverage) are removed; un-ignoring re-seeds the annotate unit
#' and queues it for (re-)annotation. Mirrors MitoPilot's ignore reconciliation.
#' @param con Open DBI connection. @param id,path,scaffold Assembly unit.
#' @param ignore New ignore flag (1 = ignored, 0 = kept).
#' @noRd
reconcile_ignore_unit <- function(con, id, path, scaffold, ignore) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  if (isTRUE(as.integer(ignore) == 1L)) {
    # HOLD, do not delete. WF1 annotates every unit as it is assembled, so the
    # rows discarded here are completed work, and "drop it, look again, put it
    # back" is the normal motion when a divergent scaffold might be a real
    # intragenomic array copy rather than a contaminant. The unit is already
    # invisible to the app, to annotation and to export -- all three filter
    # COALESCE(asm.ignore,0) = 0 -- so deleting bought nothing.
    DBI::dbExecute(con,
      "UPDATE annotate SET annotate_switch = 0 WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
  } else {
    o <- DBI::dbGetQuery(con,
      "SELECT annotate_opts, curate_opts FROM annotate WHERE ID = ? LIMIT 1", params = list(id))
    DBI::dbExecute(con,
      "INSERT OR IGNORE INTO annotate
         (ID, path, scaffold, annotate_switch, annotate_lock, annotate_opts, curate_opts, reviewed)
       VALUES (?, ?, ?, 1, 0, ?, ?, 'no')",
      params = list(id, path, scaffold,
        if (nrow(o)) o$annotate_opts[1] else "default",
        if (nrow(o)) o$curate_opts[1] else "default"))
    # Re-queue only when the unit has never been annotated. Feature rows are not
    # the test on their own: a scaffold annotation legitimately found nothing on
    # has none, and re-queueing it every time triage touches it holds the sample
    # at "in progress" forever. annotate.n_features is written (as 0) whenever
    # reconcile_one has run, so it, not the row count, marks annotation as done.
    n <- tryCatch(DBI::dbGetQuery(con,
      "SELECT COUNT(*) n FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))$n[1], error = function(e) 0L)
    done <- tryCatch(DBI::dbGetQuery(con,
      "SELECT COUNT(*) n FROM annotate
        WHERE ID = ? AND path = ? AND scaffold = ? AND n_features IS NOT NULL",
      params = list(id, path, scaffold))$n[1], error = function(e) 0L)
    DBI::dbExecute(con,
      sprintf("UPDATE annotate SET annotate_switch = %d WHERE ID = ? AND path = ? AND scaffold = ?",
        if (isTRUE(n > 0) || isTRUE(done > 0)) 2L else 1L),
      params = list(id, path, scaffold))
  }
  # holding the last queued unit finishes the sample; releasing one reopens it
  refresh_sample_state(con, id)
  invisible()
}

#' Restore a sample's assembly to its original (pre-trim) state
#'
#' Reverses one or more [trim_assembly()] calls by writing back the snapshot
#' taken at the first trim: the full-length sequence, the original annotation
#' rows, and the untrimmed coverage track. The snapshot is left in place so a
#' fresh trim + restore returns to the same original.
#' @param con Open project DB connection. @param id,path,scaffold Assembly unit.
#' @return (invisibly) the restored length.
#' @noRd
restore_assembly <- function(con, id, path = 1L, scaffold = 1L) {
  ensure_assembly_backup(con)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  bak <- DBI::dbGetQuery(con,
    "SELECT sequence, length, annotations, coverage FROM assembly_backup WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(bak) == 0 || is.na(bak$sequence[1])) {
    stop("No original assembly on record for this unit (nothing to restore).", call. = FALSE)
  }
  cur <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$sequence[1]
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(bak$sequence[1], bak$length[1], id, path, scaffold))

  # one transaction: a failed restore after a committed delete would leave the
  # unit with neither the trimmed annotations nor the snapshot ones
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    if (!is.na(bak$annotations[1]) && nzchar(bak$annotations[1])) {
      ann <- jsonlite::fromJSON(bak$annotations[1])
      if (is.data.frame(ann) && nrow(ann) > 0) {
        ann$ID <- id; ann$path <- path; ann$scaffold <- scaffold
        DBI::dbAppendTable(con, "annotations", ann)
      }
    }
  })

  write_coverage_summary(con, id, path, scaffold, parse_depths(bak$coverage))
  off <- if (is.na(cur)) -1L else as.integer(regexpr(cur, bak$sequence[1], fixed = TRUE))
  if (off > 1L) chimera_shift(con, id, path, scaffold, 1L - off, bak$length[1])
  invisible(bak$length[1])
}
