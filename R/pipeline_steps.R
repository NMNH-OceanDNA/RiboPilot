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
  bar <- if (!is.null(gff) && file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
  its <- if (!is.null(positions) && file.exists(positions)) parse_itsx(readLines(positions, warn = FALSE)) else data.frame()
  feat <- reconcile_features(bar, its)
  ts <- as.integer(Sys.time())

  DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ?", params = list(id))
  if (nrow(feat) > 0) {
    df <- data.frame(ID = id, element = feat$element, type = feat$type,
      start = feat$start, stop = feat$stop, length = feat$length, strand = feat$strand,
      tool = feat$tool, edited = 0L, notes = NA_character_, time_stamp = ts, stringsAsFactors = FALSE)
    DBI::dbAppendTable(con, "annotations", df)
  }
  missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
  DBI::dbExecute(con,
    "UPDATE annotate SET n_features = ?, missing = ?, annotate_switch = 2, time_stamp = ? WHERE ID = ?",
    params = list(nrow(feat), paste(missing, collapse = ","), ts, id))
  invisible(feat)
}

#' Open a project DB connection with WAL + busy timeout (safe for concurrent
#' writes from parallel Nextflow processes)
#' @export
open_project_db <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "PRAGMA journal_mode=WAL;")
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
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

  # replace this sample's assemblies
  DBI::dbExecute(con, "DELETE FROM assemblies WHERE ID = ?", params = list(ID))
  DBI::dbAppendTable(con, "assemblies", df)

  # per-sample summary
  n_paths <- length(unique(df$path))
  n_scaffolds <- max(table(df$path))
  length_str <- paste(sort(unique(df$length), decreasing = TRUE), collapse = ";")
  topo_str <- paste(sort(unique(stats::na.omit(df$topology))), collapse = ";")
  notes <- if (n_paths > 1) "Unable to resolve single assembly from reads"
  else if (n_scaffolds > 1) "Output contains disconnected contigs" else ""
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
  refseq <- NULL
  if (!is.null(seeds)) {
    refseq <- if (!is.na(hit_id) && hit_id %in% names(seeds)) seeds[hit_id] else seeds[which.max(Biostrings::width(seeds))]
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
    return(list(seq = paste(scaf$sequence, collapse = ""), ref_acc = sub("_.*", "", names(refseq)[1]),
      path = p, n_gaps = NA_integer_, n_redundant = NA_integer_))
  }
  list(seq = st$sequence, ref_acc = sub("_.*", "", names(refseq)[1]),
    path = p, n_gaps = st$n_gaps, n_redundant = st$n_redundant)
}
