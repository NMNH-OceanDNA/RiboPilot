# WF1 remote-GenBank BLAST: parse blastn hits into ranked reference candidates
# and set each sample's representative reference. Mirrors MitoPilot's
# blast_genbank result handling, adapted to RiboPilot's sqlite3/R-script model
# and targeting nuclear rRNA + ITS references.

#' Does this title name a nuclear rDNA record?
#'
#' Thin wrapper over [classify_hit_title()] so the binary predicate and the
#' three-state classifier can never disagree. `uninformative` reads as FALSE
#' here, which is what every existing caller already assumed.
#' @param title Character vector of subject titles.
#' @noRd
is_nuclear_rrna_title <- function(title) {
  classify_hit_title(title) == "rdna"
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
#' @param path,scaffold Assembly unit within the sample.
#' @param hits_file blastn outfmt-6 output (or a "NO_SIGNIFICANT_HITS" sentinel).
#' @param taxdump_dir Optional pinned taxdump for resolving the screen hit's
#'   lineage. Without it the lineage is left NA, which curation flags rather
#'   than treats as clean.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
parse_blast_candidates <- function(con, id, path = 1L, scaffold = 1L, hits_file, taxdump_dir = NULL) {
  ts <- as.integer(Sys.time())
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  # rp_blast.R runs against the project DB directly, so guarantee the v9 schema
  # here rather than assuming the user migrated; idempotent, as in curate_one().
  tryCatch(migrate_db(con), error = function(e) NULL)
  for (tbl in c("blast_hits", "blast_ref_candidates", "blast_screen")) {
    DBI::dbExecute(con,
      sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
      params = list(id, path, scaffold))
  }

  h <- parse_blast_hits_file(hits_file)
  if (nrow(h) == 0) {
    return(invisible(data.frame()))
  }
  h$species <- extract_species(h$title)
  h$is_nuclear_rrna <- as.integer(is_nuclear_rrna_title(h$title))

  # Rank SUBJECTS by their best bitscore. The old pident*qcovs score
  # systematically preferred short high-identity hits: a 400 bp conserved-18S-core
  # HSP against an unrelated organism scores 99*95, while a full-length repeat
  # from a genuine relative scores 95*95. Bitscore scales with alignment length.
  best <- stats::aggregate(rank_score ~ accession, data = h, FUN = max)
  best <- best[order(-best$rank_score), , drop = FALSE]
  h$hit_rank <- match(h$accession, best$accession)
  h <- h[order(h$hit_rank, -h$bitscore, h$qstart), , drop = FALSE]
  h$hsp <- stats::ave(seq_len(nrow(h)), h$hit_rank, FUN = seq_along)

  # --- resolve every distinct staxid once, rank-aligned ---
  samp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  if (is.na(samp %||% NA) || !nzchar(samp %||% "")) {
    tryCatch({
      ensure_sample_lineage_taxids(con, ids = id, taxdump_dir = taxdump_dir)
      samp <- DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?",
        params = list(id))$lineage_taxids[1]
    }, error = function(e) NULL)
  }
  h$lineage <- NA_character_
  h$lineage_taxids <- NA_character_
  lin <- tryCatch(tk_ranked_lineage(h$staxid, taxdump_dir), error = function(e) data.frame())
  if (nrow(lin) > 0) {
    i <- match(h$staxid, lin$taxid)
    h$lineage <- lin$lineage[i]
    h$lineage_taxids <- lin$lineage_taxids[i]
  }
  h$rank_distance <- if (!is.na(samp %||% NA) && nzchar(samp %||% "")) {
    vapply(h$lineage_taxids, function(x) {
      if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, samp))
    }, integer(1), USE.NAMES = FALSE)
  } else {
    NA_integer_
  }

  # --- full retained distribution ---
  hits <- data.frame(ID = id, path = path, scaffold = scaffold,
    hit_rank = h$hit_rank, hsp = h$hsp, accession = h$accession, title = h$title,
    species = h$species, staxid = h$staxid, lineage = h$lineage,
    lineage_taxids = h$lineage_taxids, rank_distance = h$rank_distance,
    pident = h$pident, qcovs = h$qcovs, evalue = h$evalue, bitscore = h$bitscore,
    aln_len = h$aln_len, qstart = h$qstart, qend = h$qend, sstart = h$sstart,
    send = h$send, is_nuclear_rrna = h$is_nuclear_rrna, time_stamp = ts,
    stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_hits", hits)

  # one row per subject, best HSP first, for the two derived tables
  subj <- h[!duplicated(h$hit_rank), , drop = FALSE]

  # --- contamination screen: top subject over everything, before any filtering ---
  write_blast_screen(con, id, path, scaffold, subj[1, ], ts = ts)

  # --- reference panel: nuclear rDNA only ---
  ref <- subj[subj$is_nuclear_rrna == 1L, , drop = FALSE]
  if (nrow(ref) == 0) {
    return(invisible(data.frame()))
  }
  ref$rank <- seq_len(nrow(ref))

  df <- data.frame(ID = id, path = path, scaffold = scaffold, rank = ref$rank, accession = ref$accession,
    species = ref$species, title = ref$title, pident = ref$pident,
    qcovs = ref$qcovs, evalue = ref$evalue, bitscore = ref$bitscore, time_stamp = ts,
    stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_ref_candidates", df)

  top <- ref[1, ]
  DBI::dbExecute(con,
    "UPDATE assemblies SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?,
       blast_evalue = ?, blast_source = 'core_nt'
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(top$accession, top$species, top$pident, top$qcovs, top$evalue, id, path, scaffold))
  invisible(df)
}

#' Parse a blastn outfmt-6 hits file into one row per HSP
#'
#' Tolerates both the pre-v9 six-column format and the current
#' `saccver stitle pident qcovs evalue staxids bitscore length qstart qend
#' sstart send`, so a re-parse of an older run still works (with the extra
#' columns NA and ranking falling back to pident*qcovs).
#' @param hits_file Path to the blastn output.
#' @return data.frame, one row per HSP; empty when there are no hits.
#' @noRd
parse_blast_hits_file <- function(hits_file) {
  empty <- data.frame()
  lines <- if (file.exists(hits_file)) readLines(hits_file, warn = FALSE) else character(0)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0 || any(grepl("NO_SIGNIFICANT_HITS", lines))) {
    return(empty)
  }
  num <- function(x) suppressWarnings(as.numeric(x))
  int <- function(x) suppressWarnings(as.integer(x))
  rows <- lapply(strsplit(lines, "\t", fixed = TRUE), function(x) {
    if (length(x) < 5) {
      return(NULL)
    }
    g <- function(i) if (length(x) >= i) x[i] else NA_character_
    data.frame(
      accession = x[1], title = x[2],
      pident = num(x[3]), qcovs = num(x[4]), evalue = num(x[5]),
      # staxids may be absent (older outfmt) or ";"-joined for merged records
      staxid = int(strsplit(g(6) %|NA|% "", ";")[[1]][1]),
      bitscore = num(g(7)), aln_len = int(g(8)),
      qstart = int(g(9)), qend = int(g(10)),
      sstart = int(g(11)), send = int(g(12)),
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(empty)
  }
  d <- do.call(rbind, rows)
  d <- d[!is.na(d$accession) & nzchar(d$accession), , drop = FALSE]
  if (nrow(d) == 0) {
    return(empty)
  }
  # Rank on bitscore when the search provided it; otherwise fall back to the
  # pre-v9 pident*qcovs so an old hits file still parses.
  # %|NA|% is scalar-only; these are vectors.
  fallback <- ifelse(is.na(d$pident), 0, d$pident) * ifelse(is.na(d$qcovs), 0, d$qcovs)
  d$rank_score <- ifelse(is.na(d$bitscore), fallback, d$bitscore)
  d
}

#' Write the contamination-screen row for a unit's top unfiltered hit
#'
#' `top` is one row of the parsed hit table, with its lineage already resolved
#' in the single batched call made by [parse_blast_candidates()]. Leaves
#' `lineage` NA when the hit carried no usable taxid; curation treats an unknown
#' lineage as a flag, not as clean.
#' @noRd
write_blast_screen <- function(con, id, path, scaffold, top, ts = as.integer(Sys.time())) {
  cols <- DBI::dbListFields(con, "blast_screen")
  row <- data.frame(
    ID = id, path = as.integer(path), scaffold = as.integer(scaffold),
    accession = top$accession, title = top$title, species = top$species,
    staxid = top$staxid, lineage = top$lineage %||% NA_character_,
    lineage_taxids = top$lineage_taxids %||% NA_character_,
    rank_distance = top$rank_distance %||% NA_integer_,
    pident = top$pident, qcovs = top$qcovs, evalue = top$evalue,
    bitscore = top$bitscore %||% NA_real_,
    is_nuclear_rrna = as.integer(top$is_nuclear_rrna %||% is_nuclear_rrna_title(top$title)),
    time_stamp = ts, stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_screen", row[, intersect(names(row), cols), drop = FALSE])
  invisible(row$lineage)
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
