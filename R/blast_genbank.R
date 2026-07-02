# WF1 remote-GenBank BLAST: parse blastn hits into ranked reference candidates
# and set each sample's representative reference. Mirrors MitoPilot's
# blast_genbank result handling, adapted to RiboPilot's sqlite3/R-script model
# and targeting nuclear rRNA + ITS references.

#' Parse remote blastn hits into ranked reference candidates
#'
#' Reads a blastn `-outfmt "6 saccver stitle pident qcovs evalue"` file for one
#' sample, dedups by accession (best pident*qcovs), ranks, writes
#' `blast_ref_candidates`, and sets the rank-1 hit as the sample's reference
#' (`stitched.reference_accession` + `assemble.blast_*`).
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param hits_file blastn outfmt-6 output (or a "NO_SIGNIFICANT_HITS" sentinel).
#' @return (invisibly) the ranked candidate data.frame.
#' @export
parse_blast_candidates <- function(con, id, hits_file) {
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM blast_ref_candidates WHERE ID = ?", params = list(id))

  lines <- if (file.exists(hits_file)) readLines(hits_file, warn = FALSE) else character(0)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0 || any(grepl("NO_SIGNIFICANT_HITS", lines))) {
    DBI::dbExecute(con,
      "UPDATE assemble SET assemble_notes = TRIM(COALESCE(assemble_notes,'') || ' [blast] no significant GenBank hits') WHERE ID = ?",
      params = list(id))
    return(invisible(data.frame()))
  }

  p <- do.call(rbind, lapply(strsplit(lines, "\t"), function(x) {
    if (length(x) < 5) return(NULL)
    data.frame(accession = x[1], title = x[2],
      pident = suppressWarnings(as.numeric(x[3])),
      qcovs = suppressWarnings(as.numeric(x[4])),
      evalue = suppressWarnings(as.numeric(x[5])), stringsAsFactors = FALSE)
  }))
  if (is.null(p) || nrow(p) == 0) return(invisible(data.frame()))

  p$species <- extract_species(p$title)
  p$score <- (p$pident %||% 0) * (p$qcovs %||% 0)
  # dedup by accession, keep best-scoring hit
  p <- p[order(-p$score), ]
  p <- p[!duplicated(p$accession), ]
  p <- p[order(-p$score), ]
  p$rank <- seq_len(nrow(p))

  df <- data.frame(ID = id, rank = p$rank, accession = p$accession, species = p$species,
    title = p$title, pident = p$pident, qcovs = p$qcovs, evalue = p$evalue,
    time_stamp = ts, stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_ref_candidates", df)

  top <- p[1, ]
  DBI::dbExecute(con,
    "UPDATE assemble SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?, blast_evalue = ? WHERE ID = ?",
    params = list(top$accession, top$species, top$pident, top$qcovs, top$evalue, id))
  DBI::dbExecute(con,
    "UPDATE stitched SET reference_accession = ? WHERE ID = ?",
    params = list(top$accession, id))
  invisible(df)
}

#' Best-effort species name from a GenBank defline title
#' @noRd
extract_species <- function(title) {
  # NCBI deflines usually start "Genus species ..."; keep the first two words.
  vapply(title, function(t) {
    if (is.na(t) || !nzchar(t)) return(NA_character_)
    w <- strsplit(trimws(t), "\\s+")[[1]]
    paste(utils::head(w, 2), collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}
