# WF1 remote-GenBank BLAST: parse blastn hits into ranked reference candidates
# and set each sample's representative reference. Mirrors MitoPilot's
# blast_genbank result handling, adapted to RiboPilot's sqlite3/R-script model
# and targeting nuclear rRNA + ITS references.

#' Is a BLAST subject title a nuclear rDNA record?
#'
#' Client-side equivalent of the restrictive `entrez_query` that used to be
#' pushed to NCBI. Applied to the hit list rather than the search space, so the
#' search itself can stay unrestricted.
#' @param title Character vector of subject titles.
#' @noRd
is_nuclear_rrna_title <- function(title) {
  t <- ifelse(is.na(title), "", title)
  grepl(NUCLEAR_RRNA_TITLE_RE, t, ignore.case = TRUE) &
    !grepl(ORGANELLAR_TITLE_RE, t, ignore.case = TRUE)
}

#' Parse remote blastn hits into ranked reference candidates + a contamination
#' screen
#'
#' Reads a blastn `-outfmt "6 saccver stitle pident qcovs evalue staxids"` file
#' for one sample and splits it two ways:
#'
#' * **Reference panel** - hits whose title is nuclear rDNA. Deduped by
#'   accession (best pident*qcovs), ranked into `blast_ref_candidates`, rank-1
#'   set as the sample's reference (`stitched.reference_accession` +
#'   `assemble.blast_*`). Same records the restricted search returned before.
#' * **Contamination screen** - the top hit over *all* of core_nt, unfiltered,
#'   into `blast_screen`. This is the only signal that can distinguish "the host"
#'   from "a symbiont" or "not nuclear rDNA at all", so it is deliberately taken
#'   before the nuclear-rDNA filter.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param hits_file blastn outfmt-6 output (or a "NO_SIGNIFICANT_HITS" sentinel).
#' @param taxdump_dir Optional pinned taxdump for resolving the screen hit's
#'   lineage. Without it the lineage is left NA, which curation flags rather
#'   than treats as clean.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
parse_blast_candidates <- function(con, id, hits_file, taxdump_dir = NULL) {
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM blast_ref_candidates WHERE ID = ?", params = list(id))
  DBI::dbExecute(con, "DELETE FROM blast_screen WHERE ID = ?", params = list(id))

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
      evalue = suppressWarnings(as.numeric(x[5])),
      # staxids may be absent (older outfmt) or ";"-joined for merged records.
      staxid = if (length(x) >= 6) {
        suppressWarnings(as.integer(strsplit(x[6], ";")[[1]][1]))
      } else {
        NA_integer_
      },
      stringsAsFactors = FALSE)
  }))
  if (is.null(p) || nrow(p) == 0) return(invisible(data.frame()))

  p$species <- extract_species(p$title)
  p$score <- (p$pident %||% 0) * (p$qcovs %||% 0)
  # dedup by accession, keep best-scoring hit
  p <- p[order(-p$score), ]
  p <- p[!duplicated(p$accession), ]
  p <- p[order(-p$score), ]

  # --- contamination screen: top hit over everything, before any filtering ---
  write_blast_screen(con, id, p[1, ], taxdump_dir = taxdump_dir, ts = ts)

  # --- reference panel: nuclear rDNA only ---
  ref <- p[is_nuclear_rrna_title(p$title), , drop = FALSE]
  if (nrow(ref) == 0) {
    DBI::dbExecute(con,
      "UPDATE assemble SET assemble_notes = TRIM(COALESCE(assemble_notes,'') || ' [blast] no nuclear rDNA hits among GenBank matches') WHERE ID = ?",
      params = list(id))
    return(invisible(data.frame()))
  }
  ref$rank <- seq_len(nrow(ref))

  df <- data.frame(ID = id, rank = ref$rank, accession = ref$accession,
    species = ref$species, title = ref$title, pident = ref$pident,
    qcovs = ref$qcovs, evalue = ref$evalue, time_stamp = ts,
    stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_ref_candidates", df)

  top <- ref[1, ]
  DBI::dbExecute(con,
    "UPDATE assemble SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?, blast_evalue = ? WHERE ID = ?",
    params = list(top$accession, top$species, top$pident, top$qcovs, top$evalue, id))
  DBI::dbExecute(con,
    "UPDATE stitched SET reference_accession = ? WHERE ID = ?",
    params = list(top$accession, id))
  invisible(df)
}

#' Write the contamination-screen row for a sample's top unfiltered hit
#'
#' Resolves the hit's lineage from its staxid when a taxdump is available,
#' falling back to the parsed species name. Leaves `lineage` NA when neither
#' works; curation treats an unknown lineage as a flag, not as clean.
#' @noRd
write_blast_screen <- function(con, id, top, taxdump_dir = NULL, ts = as.integer(Sys.time())) {
  lin <- NA_character_
  if (!is.na(top$staxid) && have_taxdump(taxdump_dir)) {
    lin <- tryCatch({
      dd <- if (file.exists(file.path(taxdump_dir, "names.dmp"))) {
        taxdump_dir
      } else {
        extract_taxdump(taxdump_dir)
      }
      tk_taxid_info(top$staxid, dd)$lineage[1]
    }, error = function(e) NA_character_)
  }
  if ((is.na(lin) || !nzchar(lin)) && !is.na(top$species) && nzchar(top$species)) {
    lin <- tryCatch(
      resolve_taxonomy(top$species, taxdump_dir = taxdump_dir)$lineage[1],
      error = function(e) NA_character_)
  }
  DBI::dbAppendTable(con, "blast_screen", data.frame(
    ID = id, accession = top$accession, title = top$title, species = top$species,
    staxid = top$staxid, lineage = lin %||% NA_character_,
    pident = top$pident, qcovs = top$qcovs, evalue = top$evalue,
    is_nuclear_rrna = as.integer(is_nuclear_rrna_title(top$title)),
    time_stamp = ts, stringsAsFactors = FALSE))
  invisible(lin)
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
