# WF1 remote-GenBank BLAST: parse blastn hits into ranked reference candidates
# and set each sample's representative reference. Mirrors MitoPilot's
# blast_genbank result handling, adapted to RiboPilot's sqlite3/R-script model
# and targeting nuclear rRNA + ITS references.

#' Floors a subject must clear to become a reference the sample is aligned to.
#'
#' Applied to the reference panel only, never to the search or to the screen's
#' hit distribution. rDNA anchors are conserved well above 85 percent within a
#' phylum. Below-half query coverage usually means a fragment rather than a
#' comparable record, but REF_CANDIDATE_MIN_ALN_BP below admits the exception:
#' a long absolute alignment against a full rDNA repeat, where low qcovs comes
#' from recovered spacer rather than a short match.
#' @noRd
REF_CANDIDATE_MIN_PIDENT <- 85
REF_CANDIDATE_MIN_QCOVS <- 50

#' Absolute aligned length that admits a subject regardless of query coverage.
#'
#' qcovs against a full rDNA repeat measures recovered intergenic spacer, not
#' match quality, so it falls as an assembly improves. 1500 bp is most of an 18S
#' and well under a 28S. NOT yet calibrated outside Actinopterygii: confirm on a
#' clade with shorter elements before treating it as settled.
#' @noRd
REF_CANDIDATE_MIN_ALN_BP <- 1500

#' Does this subject qualify as a reference the sample is aligned to?
#'
#' @param x data.frame with hit_class, pident, qcovs, aln_len.
#' @return logical vector, one element per row.
#' @noRd
ref_candidate_ok <- function(x) {
  aln <- if ("aln_len" %in% names(x)) x$aln_len else rep(NA_integer_, nrow(x))
  !is.na(x$hit_class) & x$hit_class == "rdna" &
    !is.na(x$pident) & x$pident >= REF_CANDIDATE_MIN_PIDENT &
    ((!is.na(x$qcovs) & x$qcovs >= REF_CANDIDATE_MIN_QCOVS) |
      (!is.na(aln) & aln >= REF_CANDIDATE_MIN_ALN_BP))
}

#' Parse remote blastn hits into ranked reference candidates + a contamination
#' screen
#'
#' Reads a blastn `-outfmt "6 saccver stitle pident qcovs evalue staxids"` file
#' for one sample and splits it two ways:
#'
#' * **Reference panel** - hits whose title is nuclear rDNA. Deduped by
#'   accession (best pident*qcovs), ranked into `blast_ref_candidates`, rank-1
#'   written to the unit's `assemblies.blast_*` columns. Same records the restricted search returned before.
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
  clear_unit <- function() {
    DBI::dbExecute(con,
      "DELETE FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
  }

  h <- parse_blast_hits_file(hits_file)
  if (nrow(h) == 0) {
    DBI::dbWithTransaction(con, clear_unit())
    # blast_hits is now empty for this unit; let derive_blast_products clear
    # the derived tables (and any now-stale assemblies.blast_*) too.
    return(derive_blast_products(con, id, path, scaffold, ts = ts))
  }
  h$species <- extract_species(h$title)
  h$hit_class <- classify_hit_title(h$title)
  h$is_nuclear_rrna <- as.integer(h$hit_class == "rdna")

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
    send = h$send, slen = h$slen, is_nuclear_rrna = h$is_nuclear_rrna, hit_class = h$hit_class,
    time_stamp = ts, stringsAsFactors = FALSE)
  DBI::dbWithTransaction(con, {
    clear_unit()
    DBI::dbAppendTable(con, "blast_hits", hits)
  })
  derive_blast_products(con, id, path, scaffold, ts = ts)
}

#' Derive a unit's reference panel and contamination screen from stored hits
#'
#' Reads `blast_hits` rather than a hits file, so thresholds can change without
#' re-running a remote search. Replaces the unit's rows in
#' `blast_ref_candidates` and `blast_screen` and updates `assemblies.blast_*`.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param ts Timestamp to stamp the written rows with.
#' @return (invisibly) the ranked candidate data.frame.
#' @export
derive_blast_products <- function(con, id, path = 1L, scaffold = 1L,
                                  ts = as.integer(Sys.time())) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  clear <- function() {
    for (tbl in c("blast_ref_candidates", "blast_screen")) {
      DBI::dbExecute(con,
        sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", tbl),
        params = list(id, path, scaffold))
    }
  }
  clear_assembly_ref <- function() {
    DBI::dbExecute(con,
      "UPDATE assemblies SET blast_accession = NULL, blast_species = NULL, blast_pident = NULL,
         blast_qcovs = NULL, blast_evalue = NULL, blast_source = NULL
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
  }
  h <- DBI::dbGetQuery(con,
    "SELECT * FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp",
    params = list(id, path, scaffold))
  if (nrow(h) == 0) {
    DBI::dbWithTransaction(con, { clear(); clear_assembly_ref() })
    return(invisible(data.frame()))
  }

  subj <- h[!duplicated(h$hit_rank), , drop = FALSE]
  ref <- subj[ref_candidate_ok(subj), , drop = FALSE]
  df <- NULL
  if (nrow(ref) > 0) {
    ref$rank <- seq_len(nrow(ref))
    df <- data.frame(ID = id, path = path, scaffold = scaffold, rank = ref$rank,
      accession = ref$accession, species = ref$species, title = ref$title,
      pident = ref$pident, qcovs = ref$qcovs, evalue = ref$evalue,
      bitscore = ref$bitscore, time_stamp = ts, stringsAsFactors = FALSE)
  }

  DBI::dbWithTransaction(con, {
    clear()
    write_blast_screen(con, id, path, scaffold, subj[1, ], hits = h, ts = ts)
    if (!is.null(df)) {
      DBI::dbAppendTable(con, "blast_ref_candidates", df)
      top <- ref[1, ]
      searched_source <- DBI::dbGetQuery(con,
        "SELECT source FROM blast_searched WHERE ID = ? AND path = ? AND scaffold = ?",
        params = list(id, path, scaffold))$source
      # a missing row means remote.
      blast_source <- if (length(searched_source) > 0 && identical(searched_source[1], "local")) {
        "local_rrna"
      } else {
        "core_nt"
      }
      DBI::dbExecute(con,
        "UPDATE assemblies SET blast_accession = ?, blast_species = ?, blast_pident = ?, blast_qcovs = ?,
           blast_evalue = ?, blast_source = ?
         WHERE ID = ? AND path = ? AND scaffold = ?",
        params = list(top$accession, top$species, top$pident, top$qcovs, top$evalue,
          blast_source, id, path, scaffold))
    } else {
      clear_assembly_ref()
    }
  })
  invisible(df %||% data.frame())
}

#' Parse a blastn outfmt-6 hits file into one row per HSP
#'
#' Tolerates both the pre-v9 six-column format and the current
#' `saccver stitle pident qcovs evalue staxids bitscore length qstart qend
#' sstart send slen qseqid` (slen parsed at column 13; qseqid at column 14),
#' so a re-parse of an older run still works (with the extra columns NA and
#' ranking falling back to pident*qcovs).
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
      sstart = int(g(11)), send = int(g(12)), slen = int(g(13)),
      qseqid = g(14),
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
#'
#' `hits` is the unit's full hit frame, carrying at least `hit_class` and
#' `bitscore`. It is required: the flag is a property of the whole distribution,
#' and any single-row fallback would restore the hit-number-one reading this
#' replaced.
#' @noRd
write_blast_screen <- function(con, id, path, scaffold, top, hits,
                               ts = as.integer(Sys.time())) {
  row <- data.frame(
    ID = id, path = as.integer(path), scaffold = as.integer(scaffold),
    accession = top$accession, title = top$title, species = top$species,
    staxid = top$staxid, lineage = top$lineage %||% NA_character_,
    lineage_taxids = top$lineage_taxids %||% NA_character_,
    rank_distance = top$rank_distance %||% NA_integer_,
    pident = top$pident, qcovs = top$qcovs, evalue = top$evalue,
    bitscore = top$bitscore %||% NA_real_,
    is_nuclear_rrna = consensus_hit_class(hits$hit_class, hits$bitscore),
    time_stamp = ts, stringsAsFactors = FALSE)
  DBI::dbAppendTable(con, "blast_screen", row)
  invisible(row$lineage)
}

#' Re-derive stored BLAST products for every unit that has hits
#'
#' Replays [derive_blast_products()] over `blast_hits`, so a changed threshold or
#' a fixed filter reaches existing projects without re-running a remote search.
#'
#' Changing which candidate ranks 1 invalidates `blast_ref_alignment`
#' (`compute_ref_alignment()` always aligns against the rank-1 candidate) and any
#' curation derived from it; this function does not mark those stale.
#' `inst/nextflow/bin/rp_rederive.R` is a manual tool wired into no `.nf`
#' process, so after running it the operator must re-run the reference-fetch and
#' alignment steps by hand.
#'
#' @param con Open project DB connection.
#' @param ids Optional character vector of sample IDs; default is every sample.
#' @return (invisibly) the number of units re-derived.
#' @export
rederive_blast_products <- function(con, ids = NULL) {
  sql <- "SELECT DISTINCT ID, path, scaffold FROM blast_hits"
  units <- if (is.null(ids)) {
    DBI::dbGetQuery(con, paste(sql, "ORDER BY ID, path, scaffold"))
  } else {
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbGetQuery(con, sprintf("%s WHERE ID IN (%s) ORDER BY ID, path, scaffold", sql, ph),
      params = as.list(ids))
  }
  if (nrow(units) == 0) return(invisible(0L))
  for (i in seq_len(nrow(units))) {
    derive_blast_products(con, units$ID[i], units$path[i], units$scaffold[i])
  }
  invisible(nrow(units))
}

#' Record that a unit was searched, and by which search
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param query_md5 md5 of the query sequence and the search settings, as blast.nf computes it.
#' @param n_hits Number of hits the search returned.
#' @param source Which search produced the hits, "local" or "remote".
#' @param db The database that answered: "core_nt", or the local panel's name and release.
#' @return (invisibly) TRUE.
#' @export
mark_unit_searched <- function(con, id, path = 1L, scaffold = 1L, query_md5,
                               n_hits = 0L, source = "remote", db = "") {
  source <- match.arg(source, c("remote", "local"))
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO blast_searched
       (ID, path, scaffold, query_md5, n_hits, source, db, time_stamp)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
    params = list(id, as.integer(path), as.integer(scaffold), query_md5,
      as.integer(n_hits), source, db, as.integer(Sys.time())))
  # new hits: the stored verdict was judged on the old ones
  if (DBI::dbExistsTable(con, "screen")) {
    DBI::dbExecute(con, "DELETE FROM screen WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, as.integer(path), as.integer(scaffold)))
  }
  invisible(TRUE)
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
