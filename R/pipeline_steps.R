# Per-sample pipeline steps, called by the WF1 Nextflow R-script processes
# (bin/rp_score.R, bin/rp_stitch.R, bin/rp_reconcile.R). Each takes an open DBI
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
  lineage <- DBI::dbGetQuery(con, "SELECT lineage FROM samples WHERE ID = ?", params = list(id))$lineage[1]

  scored <- score_assembly(
    path_fastas = fastas,
    gfa_path = if (length(gfa) > 0) gfa[1] else NA_character_,
    sample_lineage = lineage %||% NA_character_,
    blast_db = blast_db, taxdump_dir = td, with_staxids = FALSE)
  write_path_scores(con, id, scored)
  invisible(scored)
}

#' Stitch one sample's selected path; write the `stitched` table + return fasta
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param fasta_out Optional path to also write the stitched sequence as fasta.
#' @return (invisibly) the stitched sequence.
#' @export
stitch_one <- function(con, id, fasta_out = NULL) {
  seeds_env <- new.env()
  seeds_cache <- list(get = function(db) {
    if (is.na(db) || !nzchar(db)) {
      return(NULL)
    }
    ref <- local_ref(db)
    if (is.na(ref)) {
      return(NULL)
    }
    if (is.null(seeds_env[[ref]])) seeds_env[[ref]] <- Biostrings::readDNAStringSet(ref)
    seeds_env[[ref]]
  })
  st <- stitch_selected(con, id, seeds_cache)
  if (is.null(st)) {
    return(invisible(NULL))
  }
  # Normalize orientation to sense here, at the single point where the stitched
  # sequence is written, so WF1 coverage (computed after stitch) and WF2
  # annotation both see the same sense sequence. Deterministic, so re-stitching
  # in WF2 reproduces it -- coverage never needs reversing.
  or <- orient_rrna_sense(st$seq)
  st$seq <- or$seq
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM stitched WHERE ID = ?", params = list(id))
  DBI::dbExecute(con,
    "INSERT INTO stitched (ID, path, sequence, length, reference_accession, n_gaps, n_redundant, time_stamp)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(id, st$path, st$seq, nchar(st$seq), st$ref_acc, st$n_gaps, st$n_redundant, ts))
  if (!is.null(fasta_out)) {
    writeLines(c(paste0(">", id), st$seq), fasta_out)
  }
  invisible(st$seq)
}

#' Reconcile barrnap + ITSx outputs for one sample; write annotations
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param gff barrnap GFF file (or NULL).
#' @param positions ITSx positions file (or NULL).
#' @return (invisibly) the feature data.frame.
#' @export
reconcile_one <- function(con, id, gff = NULL, positions = NULL) {
  # Ensure any columns added after this project was created exist before writing
  # (e.g. annotations.model_coverage); reconcile runs before curate, which also
  # migrates, so an older project would otherwise fail here on the first WF2 run.
  ensure_curate_schema(con)
  bar <- if (!is.null(gff) && file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
  pos_lines <- if (!is.null(positions) && file.exists(positions)) readLines(positions, warn = FALSE) else character(0)
  its <- parse_itsx(pos_lines)
  feat <- reconcile_features(bar, its)
  chimeric <- as.integer(itsx_is_chimeric(pos_lines))
  ts <- as.integer(Sys.time())

  DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ?", params = list(id))
  if (nrow(feat) > 0) {
    df <- data.frame(ID = id, element = feat$element, type = feat$type,
      start = feat$start, stop = feat$stop, length = feat$length, strand = feat$strand,
      model_coverage = feat$model_cov %||% NA_real_,
      tool = feat$tool, edited = 0L, notes = NA_character_, time_stamp = ts, stringsAsFactors = FALSE)
    DBI::dbAppendTable(con, "annotations", df)
  }
  missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
  DBI::dbExecute(con,
    "UPDATE annotate SET n_features = ?, missing = ?, chimeric = ?, annotate_switch = 2, time_stamp = ? WHERE ID = ?",
    params = list(nrow(feat), paste(missing, collapse = ","), chimeric, ts, id))
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
  if (identical(filesystem_type(dirname(db)), "network")) "delete" else "wal"
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
  for (fa in fas) {
    seqs <- Biostrings::readDNAStringSet(fa)
    for (i in seq_along(seqs)) {
      hdr <- names(seqs)[i]
      parts <- strsplit(hdr, "\\s+")[[1]]
      idpath <- strsplit(parts[1], "\\.")[[1]] # ID, path, scaffold
      topo <- if (length(parts) > 1) parts[2] else NA_character_
      s <- as.character(seqs[[i]])
      rows[[length(rows) + 1]] <- data.frame(
        ID = ID, path = as.integer(idpath[2]), scaffold = as.integer(idpath[3]),
        topology = topo, length = nchar(s), length_raw = nchar(s),
        sequence = s, ignore = 0L, edited = 0L, time_stamp = ts,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)

  # Enforce max_paths BEFORE scoring. A tandem repeat like the rDNA cistron makes
  # GetOrganelle disentangle the graph into many isomer paths (up to its
  # --max-paths-num cap); scoring every one runs barrnap per path and dominates
  # the runtime. Keep only the max_paths paths with the greatest total assembled
  # length (the most complete repeat units). max_paths NA/<=0 means no cap.
  n_paths_raw <- length(unique(df$path))
  n_dropped <- 0L
  mp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ao.max_paths FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(ID))$max_paths[1], error = function(e) NA_integer_)
  if (!is.na(mp) && mp > 0 && n_paths_raw > mp) {
    tot <- stats::aggregate(length ~ path, data = df, FUN = sum)
    keep <- tot$path[order(-tot$length, tot$path)][seq_len(mp)]
    df <- df[df$path %in% keep, , drop = FALSE]
    n_dropped <- n_paths_raw - mp
  }

  # replace this sample's assemblies
  DBI::dbExecute(con, "DELETE FROM assemblies WHERE ID = ?", params = list(ID))
  DBI::dbAppendTable(con, "assemblies", df)

  # per-sample summary
  n_paths <- length(unique(df$path))
  n_scaffolds <- max(table(df$path))
  length_str <- paste(sort(unique(df$length), decreasing = TRUE), collapse = ";")
  topo_str <- paste(sort(unique(stats::na.omit(df$topology))), collapse = ";")
  cap_note <- if (n_dropped > 0)
    sprintf("kept top %d of %d paths (max_paths)", mp, n_paths_raw) else ""
  notes <- paste(c(
    if (n_paths > 1) "Unable to resolve single assembly from reads"
    else if (n_scaffolds > 1) "Output contains disconnected contigs" else NULL,
    if (nzchar(cap_note)) cap_note else NULL), collapse = "; ")
  DBI::dbExecute(con,
    "UPDATE assemble SET paths=?, scaffolds=?, length=?, topology=?,
     assemble_switch=2, assemble_notes=?, time_stamp=? WHERE ID=?",
    params = list(n_paths, n_scaffolds, length_str, topo_str, notes, ts, ID))
  invisible(nrow(df))
}

#' Stitch a sample's selected path against its top-hit reference
#' @noRd
stitch_selected <- function(con, id, seeds_cache) {
  sel <- DBI::dbGetQuery(con,
    "SELECT path, hit_id FROM path_scores WHERE ID = ? AND auto_select = 1", params = list(id))
  if (nrow(sel) == 0) {
    p <- DBI::dbGetQuery(con, "SELECT MIN(path) p FROM assemblies WHERE ID = ?", params = list(id))$p
    hit_id <- NA_character_
  } else {
    p <- sel$path[1]
    hit_id <- sel$hit_id[1]
  }
  if (is.na(p)) {
    return(NULL)
  }
  scaf <- DBI::dbGetQuery(con,
    "SELECT scaffold, sequence FROM assemblies WHERE ID = ? AND path = ? ORDER BY scaffold",
    params = list(id, p))
  if (nrow(scaf) == 0) {
    return(NULL)
  }
  seeds_db <- DBI::dbGetQuery(con,
    "SELECT ao.seeds_db FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(id))$seeds_db[1]

  seeds <- seeds_cache$get(seeds_db)
  # Only use the top BLAST hit as a reference; do NOT fall back to the single
  # longest seed of an arbitrary species, which orders scaffolds by an unrelated
  # sequence.
  refseq <- if (!is.null(seeds) && !is.na(hit_id) && hit_id %in% names(seeds)) {
    seeds[hit_id]
  } else {
    NULL
  }
  ref_acc <- if (!is.null(refseq)) sub("_.*", "", names(refseq)[1]) else NA_character_

  # Multi-scaffold: order by each scaffold's own anchors (reference-free). This
  # is where reference placement scrambles order for sparse clades whose only
  # GenBank hit is a single-element (e.g. 28S-only) record: an 18S scaffold
  # cannot align to a 28S reference, so it is appended out of order.
  if (nrow(scaf) > 1) {
    seqs <- stats::setNames(scaf$sequence, as.character(scaf$scaffold))
    oc <- tryCatch(order_scaffolds_by_content(seqs), error = function(e) NULL)
    if (!is.null(oc) && oc$n_anchored > 0) {
      return(list(seq = oc$sequence, ref_acc = ref_acc, path = p,
        n_gaps = 0L, n_redundant = NA_integer_))
    }
  }

  scaf_fa <- tempfile(fileext = ".fasta")
  writeLines(paste0(">", scaf$scaffold, "\n", scaf$sequence), scaf_fa)
  if (is.null(refseq)) {
    return(list(seq = paste(scaf$sequence, collapse = ""), ref_acc = NA_character_,
      path = p, n_gaps = NA_integer_, n_redundant = NA_integer_))
  }
  ref_fa <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(refseq, ref_fa)
  st <- tryCatch(stitch_scaffolds(scaf_fa, ref_fa), error = function(e) NULL)
  if (is.null(st)) {
    return(list(seq = paste(scaf$sequence, collapse = ""), ref_acc = ref_acc,
      path = p, n_gaps = NA_integer_, n_redundant = NA_integer_))
  }
  list(seq = st$sequence, ref_acc = ref_acc,
    path = p, n_gaps = st$n_gaps, n_redundant = st$n_redundant)
}

#' Trim a sample's stitched assembly to [from, to] (stitched coordinates),
#' shifting its annotations and coverage to match. Features that fall entirely
#' outside the kept range are dropped; overlapping features are clamped. Used by
#' the Annotate tab's "Trim assembly" control.
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param from,to 1-based inclusive kept range in stitched coordinates.
#' @return (invisibly) the new stitched length.
#' @noRd
trim_assembly <- function(con, id, from, to) {
  s <- DBI::dbGetQuery(con, "SELECT sequence FROM stitched WHERE ID = ?", params = list(id))
  if (nrow(s) == 0 || is.na(s$sequence[1]) || !nzchar(s$sequence[1])) {
    stop("No stitched sequence for ", id, call. = FALSE)
  }
  n <- nchar(s$sequence[1])
  from <- max(1L, as.integer(from)); to <- min(n, as.integer(to))
  if (from > to) stop("Invalid trim range.", call. = FALSE)
  new_seq <- substr(s$sequence[1], from, to); new_len <- nchar(new_seq)

  DBI::dbExecute(con, "UPDATE stitched SET sequence = ?, length = ? WHERE ID = ?",
    params = list(new_seq, new_len, id))

  ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ?", params = list(id))
  for (i in seq_len(nrow(ann))) {
    ns <- ann$start[i] - from + 1; nt <- ann$stop[i] - from + 1
    if (nt < 1 || ns > new_len) {
      DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND element = ?",
        params = list(id, ann$element[i]))
      next
    }
    ns <- max(1L, ns); nt <- min(new_len, nt)
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND element = ?",
      params = list(ns, nt, nt - ns + 1, id, ann$element[i]))
  }

  cov <- DBI::dbGetQuery(con, "SELECT depths FROM coverage WHERE ID = ?", params = list(id))
  if (nrow(cov) > 0 && !is.na(cov$depths[1]) && nzchar(cov$depths[1])) {
    dep <- as.numeric(strsplit(cov$depths[1], ",")[[1]])
    dep <- dep[from:min(length(dep), to)]
    DBI::dbExecute(con,
      "UPDATE coverage SET depths = ?, mean_depth = ?, min_depth = ?, max_depth = ? WHERE ID = ?",
      params = list(paste(dep, collapse = ","), mean(dep), min(dep), max(dep), id))
  }
  invisible(new_len)
}
