# Per-feature BLAST: each core feature of a contig searched on its own against
# the bundled local panel, so a feature with no hit, a non-rDNA hit, or a hit
# in another phylum shows up even when the whole-contig BLAST looks fine.

FEATURE_BLAST_ELEMENTS <- c("18S", "ITS1", "5.8S", "ITS2", "28S")
FEATURE_BLAST_MIN_LEN <- 50L
FEATURE_BLAST_LOWCONF_LEN <- 100L
FEATURE_BLAST_EVALUE <- 1e-5
FEATURE_BLAST_KEEP_HITS <- 10L

#' blastn task for a feature length: NA (too small), "blastn", or "megablast"
#' @noRd
feature_blast_task <- function(len) {
  len <- suppressWarnings(as.integer(len))
  ifelse(is.na(len) | len < FEATURE_BLAST_MIN_LEN, NA_character_,
    ifelse(len < FEATURE_BLAST_LOWCONF_LEN, "blastn", "megablast"))
}

#' Core features of a unit with their sequences cut from the contig
#' @return data.frame(element, start, stop, length, sequence), in
#'   FEATURE_BLAST_ELEMENTS order.
#' @noRd
feature_blast_targets <- function(con, id, path, scaffold, review = FALSE) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  seq <- DBI::dbGetQuery(con, sprintf("SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("assemblies", review)), params = p)$sequence[1]
  f <- DBI::dbGetQuery(con, sprintf("SELECT element, start, stop FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("annotations", review)), params = p)
  f <- f[f$element %in% FEATURE_BLAST_ELEMENTS & !is.na(f$start) & !is.na(f$stop), , drop = FALSE]
  f <- f[order(match(f$element, FEATURE_BLAST_ELEMENTS)), , drop = FALSE]
  if (is.na(seq %||% NA) || !nrow(f)) {
    return(data.frame(element = character(0), start = integer(0), stop = integer(0),
      length = integer(0), sequence = character(0), stringsAsFactors = FALSE))
  }
  f$sequence <- substring(seq, f$start, f$stop)
  f$length <- nchar(f$sequence)
  rownames(f) <- NULL
  f[, c("element", "start", "stop", "length", "sequence")]
}

#' One call per feature from a parsed hit frame
#'
#' `hits` comes from parse_blast_hits_file() with qseqid = element, plus
#' lineage and lineage_taxids columns. Near-tied hits (within CHIMERA_TIE_FRAC
#' of the top bitscore) collapse to the clade they share, as the chimera
#' windows do; the feature is off-target only when no tied hit (or the best
#' hit, without a tie) is in the sample's phylum, and every phylum is known.
#' @return list(calls = one row per target, hits = top hits per element).
#' @noRd
feature_blast_calls <- function(hits, targets, sample_taxids = NA_character_) {
  has_sample <- !is.na(sample_taxids) && nzchar(sample_taxids)
  rd_of <- function(lt) vapply(lt, function(x)
    if (is.na(x) || !has_sample) NA_integer_ else as.integer(rank_distance(x, sample_taxids)),
    integer(1), USE.NAMES = FALSE)
  if (nrow(hits)) {
    hits <- hits[!is.na(hits$evalue) & hits$evalue <= FEATURE_BLAST_EVALUE, , drop = FALSE]
    hits <- hits[order(hits$qseqid, -hits$bitscore), , drop = FALSE]
    hits <- hits[!duplicated(hits[c("qseqid", "accession")]), , drop = FALSE]
    hits$hit_class <- classify_hit_title(hits$title)
  }
  task <- feature_blast_task(targets$length)
  calls <- lapply(seq_len(nrow(targets)), function(k) {
    t <- targets[k, ]
    row <- data.frame(element = t$element, start = t$start, stop = t$stop, length = t$length,
      task = task[k], low_conf = as.integer(!is.na(task[k]) & task[k] == "blastn"),
      status = NA_character_, note = NA_character_, accession = NA_character_,
      species = NA_character_, staxid = NA_integer_, lineage = NA_character_,
      lineage_taxids = NA_character_, n_tied = NA_integer_, match_rank = NA_character_,
      rank_distance = NA_integer_, pident = NA_real_, qcovs = NA_real_, evalue = NA_real_,
      bitscore = NA_real_, stringsAsFactors = FALSE)
    if (is.na(task[k])) {
      row$status <- "too_small"
      row$note <- "no BLAST, feature too small"
      return(row)
    }
    h <- if (nrow(hits)) hits[hits$qseqid == t$element, , drop = FALSE] else hits
    if (!nrow(h)) {
      row$status <- "no_hit"
      return(row)
    }
    b <- h[1, ]
    row[c("accession", "staxid", "lineage", "lineage_taxids", "pident", "qcovs", "evalue", "bitscore")] <-
      b[c("accession", "staxid", "lineage", "lineage_taxids", "pident", "qcovs", "evalue", "bitscore")]
    row$species <- extract_species(b$title)
    if (b$hit_class == "other") {
      row$status <- "non_rdna"
      return(row)
    }
    hl <- h[!is.na(h$lineage_taxids), , drop = FALSE]
    tied <- hl[0, ]
    if (nrow(hl) >= 2) {
      lca <- bitscore_weighted_lca(hl$lineage_taxids, hl$bitscore, frac = CHIMERA_TIE_FRAC, min_hits = 1L)
      if (lca$n_hits >= 2 && !identical(lca$rank, "species")) {
        tied <- hl[seq_len(lca$n_hits), , drop = FALSE]
        row$lineage <- chimera_cut_lineage(hl$lineage[1], lca$rank)
        row$lineage_taxids <- chimera_cut_lineage(hl$lineage_taxids[1], lca$rank)
        row$n_tied <- lca$n_hits
        p <- strsplit(row$lineage %|NA|% "", ";", fixed = TRUE)[[1]]
        clade <- utils::tail(c(NA_character_, p[nzchar(p)]), 1)
        if (!is.na(clade)) row$species <- sprintf("%s (%d tied hits)", clade, lca$n_hits)
      }
    }
    row$rank_distance <- rd_of(row$lineage_taxids)
    row$match_rank <- rank_distance_label(row$rank_distance)
    # off target only when every candidate's phylum is known and none is the sample's
    phy <- function(x) vapply(x, function(v) parse_rank_taxids(v)[2], integer(1), USE.NAMES = FALSE)
    cand <- if (nrow(tied)) tied$lineage_taxids else b$lineage_taxids
    sp <- if (has_sample) phy(sample_taxids) else NA_integer_
    cp <- if (any(is.na(cand))) NA_integer_ else phy(cand)
    row$status <- if (!has_sample) {
      row$note <- "sample lineage unknown"
      "ok"
    } else if (is.na(b$lineage_taxids)) {
      row$note <- "hit has no taxonomy"
      "ok"
    } else if (is.na(sp) || anyNA(cp)) {
      row$note <- "phylum unknown"
      "ok"
    } else if (!sp %in% cp) "off_target" else "ok"
    row
  })
  calls <- if (length(calls)) do.call(rbind, calls) else NULL
  keep <- if (nrow(hits)) {
    hits$hit_rank <- stats::ave(seq_len(nrow(hits)), hits$qseqid, FUN = seq_along)
    k <- hits[hits$hit_rank <= FEATURE_BLAST_KEEP_HITS, , drop = FALSE]
    data.frame(element = k$qseqid, hit_rank = k$hit_rank, accession = k$accession, title = k$title,
      staxid = k$staxid, lineage = k$lineage, pident = k$pident, qcovs = k$qcovs, evalue = k$evalue,
      bitscore = k$bitscore, hit_class = k$hit_class, stringsAsFactors = FALSE)
  } else NULL
  list(calls = calls, hits = keep)
}

#' @noRd
feature_blast_run <- function(query, dbarg, task, out, threads = 1L, blastn = "blastn") {
  status <- system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query), "-task", task,
    "-outfmt", shQuote("6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"),
    "-evalue", FEATURE_BLAST_EVALUE, "-max_target_seqs", "50", "-max_hsps", "1", "-num_threads", threads),
    stdout = out, stderr = FALSE)
  if (!identical(as.integer(status), 0L)) stop(sprintf("blastn -task %s exited with status %s", task, status))
  parse_blast_hits_file(out)
}

#' BLAST each core feature of a unit and store the calls
#'
#' Features under FEATURE_BLAST_MIN_LEN are not searched. With no `dbarg`, or
#' when blastn fails, searchable features are stored as not_run with a note.
#' @param review TRUE reads and writes the contig's review copy.
#' @return (invisibly) the stored calls.
#' @export
feature_blast_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                               threads = 1L, workdir = tempfile("feature_blast_"), review = FALSE,
                               blastn = "blastn") {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  targets <- feature_blast_targets(con, id, path, scaffold, review)
  task <- feature_blast_task(targets$length)
  hits <- data.frame()
  dc <- character(0)
  fail <- if (!isTRUE(nzchar(dbarg, keepNA = TRUE))) "no local BLAST panel" else NULL
  if (is.null(fail) && any(!is.na(task))) {
    dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
    run <- function(sel, tk) {
      if (!any(sel)) return(data.frame())
      q <- file.path(workdir, paste0(tk, ".fa"))
      writeLines(as.vector(rbind(paste0(">", targets$element[sel]), targets$sequence[sel])), q)
      feature_blast_run(q, dbarg, tk, file.path(workdir, paste0(tk, ".txt")), threads, blastn)
    }
    fail <- tryCatch({
      mega <- !is.na(task) & task == "megablast"
      h1 <- run(mega, "megablast")
      seen <- if (nrow(h1)) unique(h1$qseqid) else character(0)
      h2 <- run(mega & !targets$element %in% seen, "dc-megablast")
      h3 <- run(!is.na(task) & task == "blastn", "blastn")
      hits <- do.call(rbind, Filter(nrow, list(h1, h2, h3))) %||% data.frame()
      dc <- if (nrow(h2)) unique(h2$qseqid) else character(0)
      NULL
    }, error = function(e) paste("BLAST failed:", conditionMessage(e)))
  }
  if (nrow(hits)) {
    lin <- tryCatch(tk_ranked_lineage(hits$staxid, taxdump_dir), error = function(e) data.frame())
    i <- if (nrow(lin)) match(hits$staxid, lin$taxid) else rep(NA_integer_, nrow(hits))
    hits$lineage <- if (nrow(lin)) lin$lineage[i] else NA_character_
    hits$lineage_taxids <- if (nrow(lin)) lin$lineage_taxids[i] else NA_character_
  }
  tryCatch(ensure_sample_lineage_taxids(con, ids = id, taxdump_dir = taxdump_dir), error = function(e) NULL)
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?",
    params = list(id))$lineage_taxids[1], error = function(e) NA_character_)
  res <- feature_blast_calls(hits, targets, samp %||% NA_character_)
  calls <- res$calls
  if (!is.null(calls)) calls$task[calls$element %in% dc] <- "dc-megablast"
  if (!is.null(calls) && !is.null(fail)) {
    calls$status[calls$status != "too_small"] <- "not_run"
    calls$note[calls$status == "not_run"] <- fail
    res$hits <- NULL
  }
  ts <- as.integer(Sys.time())
  key <- function(d) cbind(data.frame(ID = rep(id, nrow(d)), path = path, scaffold = scaffold,
    stringsAsFactors = FALSE), d)
  fb <- review_tbl("feature_blast", review)
  fbh <- review_tbl("feature_blast_hits", review)
  with_savepoint(con, {
    for (t in c(fb, fbh)) DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", t),
      params = list(id, path, scaffold))
    if (!is.null(calls)) {
      calls$db_source <- dbarg %||% NA_character_
      calls$time_stamp <- ts
      DBI::dbAppendTable(con, fb, key(calls))
    }
    if (!is.null(res$hits)) DBI::dbAppendTable(con, fbh, key(res$hits))
  }, name = "rp_feature_blast")
  invisible(calls)
}

#' The Annotation Review "Feature BLAST warnings" text per unit
#'
#' Lists no hit, non-rDNA, off-target, too-small, and weak (found only by the
#' dc-megablast retry) features in element order; "" when a unit has none.
#' @param fb Rows of feature_blast (or review_feature_blast_v) with ID, path,
#'   scaffold, element, status, lineage, and optionally task.
#' @return data.frame(ID, path, scaffold, feature_blast_warnings).
#' @noRd
feature_blast_warning_text <- function(fb) {
  if (is.null(fb$task)) fb$task <- rep(NA_character_, nrow(fb))
  fb$status[fb$status %in% "ok" & fb$task %in% "dc-megablast"] <- "weak"
  fb <- fb[fb$status %in% c("no_hit", "non_rdna", "off_target", "too_small", "weak"), , drop = FALSE]
  if (!nrow(fb)) {
    return(data.frame(ID = character(0), path = integer(0), scaffold = integer(0),
      feature_blast_warnings = character(0), stringsAsFactors = FALSE))
  }
  fb <- fb[order(fb$ID, fb$path, fb$scaffold, match(fb$element, FEATURE_BLAST_ELEMENTS)), , drop = FALSE]
  taxon <- chimera_rank_name(fb$lineage, "phylum")
  deepest <- vapply(strsplit(ifelse(is.na(fb$lineage), "", fb$lineage), ";", fixed = TRUE),
    function(p) utils::tail(c(NA_character_, p[nzchar(p)]), 1), character(1))
  taxon <- ifelse(is.na(taxon), deepest, taxon)
  msg <- paste0(fb$element, ": ", c(no_hit = "no hit", non_rdna = "non-rDNA hit",
    too_small = "no BLAST, feature too small", off_target = "off-target", weak = "weak hit")[fb$status],
    ifelse(fb$status == "off_target" & !is.na(taxon), paste0(" (", taxon, ")"), ""))
  k <- paste(fb$ID, fb$path, fb$scaffold, sep = "\r")
  out <- fb[!duplicated(k), c("ID", "path", "scaffold"), drop = FALSE]
  out$feature_blast_warnings <- as.vector(tapply(msg, factor(k, levels = unique(k)), paste, collapse = "; "))
  rownames(out) <- NULL
  out
}
