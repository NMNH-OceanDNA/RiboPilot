# WF1 reference fetch: store an NCBI reference record (sequence + rRNA/ITS
# features) fetched via edirect efetch, so the app can draw a reference synteny
# map and per-feature sample-vs-reference alignments.

#' Map an NCBI feature product string to a RiboPilot element name
#' @noRd
ref_product_to_element <- function(product, key) {
  p <- tolower(product %||% "")
  if (grepl("18s|small subunit", p)) return("18S")
  if (grepl("5\\.8s", p)) return("5.8S")
  if (grepl("28s|26s|large subunit", p)) return("28S")
  if (grepl("transcribed spacer 1|its1|its 1", p)) return("ITS1")
  if (grepl("transcribed spacer 2|its2|its 2", p)) return("ITS2")
  NA_character_
}

#' Parse an NCBI feature table (efetch -format ft) into rRNA/ITS features
#'
#' @param ft_lines Character vector (lines of the feature-table file).
#' @return data.frame(element, type, start, stop, strand); rRNA/ITS rows only.
#' @noRd
parse_feature_table <- function(ft_lines) {
  feats <- list()
  for (ln in ft_lines) {
    if (grepl("^>Feature", ln)) next
    parts <- strsplit(ln, "\t")[[1]]
    # coordinate line: "<start> <stop> <key>" (start/stop may carry </> partials)
    if (grepl("^[<>0-9]", ln) && length(parts) >= 3 &&
        grepl("^[<>]?[0-9]+$", parts[1]) && grepl("^[<>]?[0-9]+$", parts[2])) {
      s <- as.integer(gsub("[<>]", "", parts[1]))
      e <- as.integer(gsub("[<>]", "", parts[2]))
      feats[[length(feats) + 1]] <- list(key = parts[3], product = NA_character_,
        start = min(s, e), stop = max(s, e), strand = if (s > e) "-" else "+")
    } else if (length(parts) >= 2 && parts[length(parts) - 1] == "product") {
      if (length(feats) > 0) feats[[length(feats)]]$product <- parts[length(parts)]
    }
  }
  if (length(feats) == 0) return(data.frame())
  df <- do.call(rbind, lapply(feats, function(f) {
    el <- ref_product_to_element(f$product, f$key)
    if (is.na(el)) return(NULL)
    data.frame(element = el, type = f$key, start = f$start, stop = f$stop,
      strand = f$strand, stringsAsFactors = FALSE)
  }))
  if (is.null(df)) return(data.frame())
  df[!duplicated(df$element), , drop = FALSE]
}

#' TRUE if the file's first non-blank line is a FASTA header
#' @noRd
is_fasta_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(FALSE)
  ls <- readLines(path, warn = FALSE)
  ls <- ls[nzchar(trimws(ls))]
  length(ls) > 0 && grepl("^>", ls[1])
}

#' Read a single-record FASTA into a bare sequence string
#'
#' Strips only whitespace and digits; IUPAC ambiguity codes are kept so the
#' sequence stays the same length as the GenBank record the feature-table
#' coordinates refer to.
#' @noRd
read_fasta_seq <- function(path) {
  if (!is_fasta_file(path)) return("")
  ls <- readLines(path, warn = FALSE)
  paste(gsub("[[:space:][:digit:]]", "", ls[!grepl("^>", ls)]), collapse = "")
}

#' Store one fetched NCBI reference (sequence + features) in the project DB
#'
#' @param con Open project DB connection.
#' @param accession NCBI accession.
#' @param fasta_file efetch `-format fasta` output.
#' @param ft_file efetch `-format ft` output.
#' @param title Optional defline title (for species).
#' @return (invisibly) the reference sequence length.
#' @export
write_blast_ref <- function(con, accession, fasta_file, ft_file, title = NA_character_) {
  ts <- as.integer(Sys.time())
  # efetch can answer HTTP 200 with an XML error body; that is not a record.
  if (!is_fasta_file(fasta_file)) return(invisible(0L))
  seq <- read_fasta_seq(fasta_file)
  if (!nzchar(seq)) return(invisible(0L))
  if (is.na(title) && file.exists(fasta_file)) {
    h <- grep("^>", readLines(fasta_file, warn = FALSE), value = TRUE)
    if (length(h)) title <- sub("^>\\S+\\s*", "", h[1])
  }
  species <- extract_species(title)

  ft <- if (file.exists(ft_file)) parse_feature_table(readLines(ft_file, warn = FALSE)) else data.frame()
  # one transaction: a failed append after a committed delete would drop the
  # stored record and its features with nothing to put back
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM blast_ref_sequences WHERE accession = ?", params = list(accession))
    DBI::dbAppendTable(con, "blast_ref_sequences", data.frame(accession = accession,
      species = species, title = title, taxid = NA_character_, lineage = NA_character_,
      sequence = seq, length = nchar(seq), time_stamp = ts, stringsAsFactors = FALSE))
    DBI::dbExecute(con, "DELETE FROM blast_ref_annotations WHERE accession = ?", params = list(accession))
    if (nrow(ft) > 0) {
      DBI::dbAppendTable(con, "blast_ref_annotations", data.frame(accession = accession,
        element = ft$element, type = ft$type, start = ft$start, stop = ft$stop,
        strand = ft$strand, time_stamp = ts, stringsAsFactors = FALSE))
    }
  })
  invisible(nchar(seq))
}
