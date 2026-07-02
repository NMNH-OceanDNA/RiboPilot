# Annotation driver: for each sample, stitch the auto-selected path into one
# ordered sequence, annotate its nrDNA feature model, and store the stitched
# sequence + features in the database.

#' Annotate assembled samples for a project
#'
#' For each sample with a scored/auto-selected path (falls back to path 1),
#' stitches that path against its top-BLAST-hit reference, runs [annotate_seq()],
#' and writes the `stitched`, `annotations`, and `annotate` tables.
#'
#' @param path Project directory.
#' @param out_dir Pipeline output dir (default `<path>/out`) - unused directly
#'   but kept for symmetry; sequences come from the `assemblies` table.
#' @param con Optional open DBI connection.
#' @param cpu CPUs for ITSx.
#' @return (invisibly) NULL.
#' @export
annotate_project <- function(path = ".", out_dir = file.path(path, "out"), con = NULL, cpu = 2) {
  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
    on.exit(DBI::dbDisconnect(con), add = TRUE)
  }
  ts <- as.integer(Sys.time())

  ids <- DBI::dbGetQuery(con, "SELECT DISTINCT ID FROM assemblies")$ID
  if (length(ids) == 0) {
    message("No assembled samples to annotate.")
    return(invisible())
  }

  # seed DB per sample -> reference sequences (cache by path)
  seed_cache <- list()
  seeds_for <- function(seeds_db) {
    if (is.na(seeds_db) || !nzchar(seeds_db)) {
      return(NULL)
    }
    ref <- local_ref(seeds_db)
    if (is.na(ref)) {
      return(NULL)
    }
    if (is.null(seed_cache[[ref]])) {
      seed_cache[[ref]] <<- Biostrings::readDNAStringSet(ref)
    }
    seed_cache[[ref]]
  }

  for (id in ids) {
    # auto-selected path (fall back to lowest path number)
    sel <- DBI::dbGetQuery(con,
      "SELECT path, hit_id FROM path_scores WHERE ID = ? AND auto_select = 1", params = list(id))
    if (nrow(sel) == 0) {
      p <- DBI::dbGetQuery(con, "SELECT MIN(path) p FROM assemblies WHERE ID = ?", params = list(id))$p
      hit_id <- NA_character_
    } else {
      p <- sel$path[1]
      hit_id <- sel$hit_id[1]
    }
    if (is.na(p)) next

    scaf <- DBI::dbGetQuery(con,
      "SELECT scaffold, sequence FROM assemblies WHERE ID = ? AND path = ? ORDER BY scaffold",
      params = list(id, p))
    if (nrow(scaf) == 0) next

    seeds_db <- DBI::dbGetQuery(con,
      "SELECT ao.seeds_db FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
      params = list(id))$seeds_db[1]
    seeds <- seeds_for(seeds_db)

    # reference = top-hit seed, else longest seed
    refseq <- NULL
    if (!is.null(seeds)) {
      if (!is.na(hit_id) && hit_id %in% names(seeds)) {
        refseq <- seeds[hit_id]
      } else {
        refseq <- seeds[which.max(Biostrings::width(seeds))]
      }
    }

    # write scaffolds + reference, stitch
    scaf_fa <- tempfile(fileext = ".fasta")
    writeLines(paste0(">", scaf$scaffold, "\n", scaf$sequence), scaf_fa)
    stitched <- NULL
    ref_acc <- NA_character_
    if (!is.null(refseq)) {
      ref_fa <- tempfile(fileext = ".fasta")
      Biostrings::writeXStringSet(refseq, ref_fa)
      ref_acc <- sub("_.*", "", names(refseq)[1])
      st <- tryCatch(stitch_scaffolds(scaf_fa, ref_fa), error = function(e) NULL)
      if (!is.null(st)) stitched <- st
    }
    # fall back: no reference -> concatenate scaffolds as-is
    seq_to_annotate <- if (!is.null(stitched)) stitched$sequence else paste(scaf$sequence, collapse = "")
    n_gaps <- if (!is.null(stitched)) stitched$n_gaps else NA_integer_
    n_red <- if (!is.null(stitched)) stitched$n_redundant else NA_integer_

    DBI::dbExecute(con, "DELETE FROM stitched WHERE ID = ?", params = list(id))
    DBI::dbExecute(con,
      "INSERT INTO stitched (ID, path, sequence, length, reference_accession, n_gaps, n_redundant, time_stamp)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(id, p, seq_to_annotate, nchar(seq_to_annotate), ref_acc, n_gaps, n_red, ts))

    # annotate
    ann_fa <- tempfile(fileext = ".fasta")
    writeLines(c(paste0(">", id), seq_to_annotate), ann_fa)
    feat <- tryCatch(annotate_seq(ann_fa, cpu = cpu), error = function(e) data.frame())

    DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ?", params = list(id))
    if (nrow(feat) > 0) {
      df <- data.frame(ID = id, element = feat$element, type = feat$type,
        start = feat$start, stop = feat$stop, length = feat$length, strand = feat$strand,
        tool = feat$tool, edited = 0L, notes = NA_character_, time_stamp = ts,
        stringsAsFactors = FALSE)
      DBI::dbAppendTable(con, "annotations", df)
    }
    missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
    DBI::dbExecute(con,
      "UPDATE annotate SET n_features = ?, missing = ?, annotate_switch = 2, time_stamp = ? WHERE ID = ?",
      params = list(nrow(feat), paste(missing, collapse = ","), ts, id))
    message("Annotated ", id, ": ", nrow(feat), " features",
      if (length(missing)) paste0(" (missing: ", paste(missing, collapse = ","), ")") else "")
  }
  invisible()
}
