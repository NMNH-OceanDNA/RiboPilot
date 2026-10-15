# Paths alignment panel: homologous contigs across a sample's assembly paths,
# aligned with DECIPHER and shown in the sequence viewer's alignment mode.

#' Homologous sets of contigs: connected components of "shares an anchor"
#' @param feats data.frame(path, scaffold, element), one row per annotation.
#' @param link NULL (any shared anchor links) or function(ua, ub, element)
#'   giving TRUE when that pair's copies of the element are homologous.
#' @return list of list(label, units = data.frame(path, scaffold)); only sets
#'   of two or more contigs, largest first.
#' @noRd
paths_homology_sets <- function(feats, link = NULL) {
  f <- feats[feats$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (!nrow(f)) return(list())
  f$unit <- paste(f$path, f$scaffold)
  units <- unique(f$unit)
  comp <- stats::setNames(seq_along(units), units)
  for (el in unique(f$element)) {
    members <- unique(f$unit[f$element == el])
    if (is.null(link)) {
      comp[comp %in% comp[members]] <- min(comp[members])
      next
    }
    if (length(members) < 2) next
    for (pr in utils::combn(members, 2, simplify = FALSE)) {
      if (comp[pr[1]] != comp[pr[2]] && isTRUE(link(pr[1], pr[2], el))) {
        comp[comp %in% comp[pr]] <- min(comp[pr])
      }
    }
  }
  sets <- split(names(comp), comp)
  sets <- Filter(function(u) length(u) >= 2, sets)
  if (!length(sets)) return(list())
  out <- lapply(sets, function(u) {
    pu <- do.call(rbind, strsplit(u, " ", fixed = TRUE))
    units <- data.frame(path = as.integer(pu[, 1]), scaffold = as.integer(pu[, 2]))
    units <- units[order(units$path, units$scaffold), , drop = FALSE]
    rownames(units) <- NULL
    anchors <- SCREEN_ANCHOR_ELEMENTS[SCREEN_ANCHOR_ELEMENTS %in% f$element[f$unit %in% u]]
    list(label = sprintf("%s (%d contigs)", paste(anchors, collapse = "-"), nrow(units)), units = units)
  })
  n <- vapply(out, function(s) nrow(s$units), 1L)
  first <- vapply(out, function(s) s$units$path[1] * 1e6 + s$units$scaffold[1], 1)
  unname(out[order(-n, first)])
}

#' Link function for paths_homology_sets(): two copies of an element link when
#' their ends-free alignment covers REF_ALN_MIN_COV of the shorter copy over at
#' least REF_ALN_MIN_BP bases
#' @noRd
paths_overlap_link <- function(con, id) {
  asm <- DBI::dbGetQuery(con, "SELECT path, scaffold, sequence FROM assemblies WHERE ID = ?", params = list(id))
  n <- DBI::dbGetQuery(con, "SELECT path, scaffold, element, start, stop FROM annotations WHERE ID = ?",
    params = list(id))
  seqs <- stats::setNames(asm$sequence, paste(asm$path, asm$scaffold))
  n$unit <- paste(n$path, n$scaffold)
  memo <- list()
  piece <- function(u, el) {
    r <- n[n$unit == u & n$element == el, , drop = FALSE]
    if (!nrow(r) || is.na(seqs[u])) return(NA_character_)
    substr(seqs[[u]], r$start[1], r$stop[1])
  }
  function(ua, ub, element) {
    key <- paste(c(sort(c(ua, ub)), element), collapse = " ")
    if (!is.null(memo[[key]])) return(memo[[key]])
    a <- piece(ua, element); b <- piece(ub, element)
    ok <- !is.na(a) && !is.na(b) && min(nchar(a), nchar(b)) >= REF_ALN_MIN_BP &&
      tryCatch({
        aln <- pairwise_nt(a, b, "overlap")
        s <- if (nchar(a) <= nchar(b)) as.character(pwalign::pattern(aln)) else as.character(pwalign::subject(aln))
        isTRUE(alignment_coverage(s, min(nchar(a), nchar(b))) >= REF_ALN_MIN_COV &&
          nchar(gsub("-", "", s)) >= REF_ALN_MIN_BP)
      }, error = function(e) FALSE)
    memo[[key]] <<- ok
    ok
  }
}

#' Align the contigs of one set
#' @param seqs Named character vector (unit keys as names).
#' @return Aligned DNAStringSet or NULL.
#' @noRd
paths_align <- function(seqs) {
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  if (length(seqs) < 2) return(NULL)
  # ponytail: no cap on contig count; DECIPHER on a few 10 kb contigs takes seconds
  dna <- Biostrings::DNAStringSet(toupper(seqs))
  tryCatch(DECIPHER::AlignSeqs(dna, verbose = FALSE), error = function(e) NULL)
}

#' Viewer payload for an alignment: one row per contig, features in columns
#' @noRd
paths_payload <- function(aln, feats, labels) {
  seqs <- as.character(aln)
  rows <- lapply(names(seqs), function(u) {
    s <- seqs[[u]]
    cols <- which(strsplit(s, "", fixed = TRUE)[[1]] != "-")
    f <- feats[feats$unit == u & !is.na(feats$start) & !is.na(feats$stop), , drop = FALSE]
    list(unit = u, label = unname(labels[u]) %|NA|% u, seq = s,
      features = lapply(seq_len(nrow(f)), function(i) list(
        row = f$element[i], gene = f$element[i],
        pos1 = as.integer(cols[f$start[i]]), pos2 = as.integer(cols[f$stop[i]]),
        dir = if (identical(f$strand[i], "-")) "-" else "+",
        partial5 = isTRUE(f$partial_start[i] == 1), partial3 = isTRUE(f$partial_stop[i] == 1),
        color = unname(RRNA_ELEMENT_COLORS[f$element[i]] %|NA|% "#888888"))))
  })
  list(len = as.integer(nchar(seqs[1])), rows = rows)
}

#' Columns where the rows disagree (gaps count)
#' @noRd
paths_diff_count <- function(rows) {
  if (length(rows) < 2) return(0L)
  s <- vapply(rows, `[[`, "", "seq")
  f <- regexpr("[^-]", s); l <- nchar(sub("-+$", "", s))
  m <- do.call(rbind, strsplit(s, "", fixed = TRUE))
  sum(vapply(seq_len(ncol(m)), function(j) {
    x <- m[f <= j & j <= l, j]
    length(x) > 1 && length(unique(x)) > 1
  }, TRUE))
}

#' Row labels for a set's contigs: "p<path>.s<scaffold>" plus triage when set
#' @return Named character vector keyed by "p<path>.s<scaffold>".
#' @noRd
paths_labels <- function(con, id, units) {
  keys <- sprintf("p%d.s%d", units$path, units$scaffold)
  a <- DBI::dbGetQuery(con, "SELECT path, scaffold, triage FROM assemblies WHERE ID = ?", params = list(id))
  tr <- a$triage[match(paste(units$path, units$scaffold), paste(a$path, a$scaffold))]
  stats::setNames(trimws(paste(keys, ifelse(is.na(tr), "", tr))), keys)
}

#' Align one homologous set and send it to the paths viewer
#' @return The payload, or NULL when the set could not be aligned.
#' @noRd
paths_send <- function(session, id, set) {
  con <- app_con(); on.exit(DBI::dbDisconnect(con))
  u <- set$units
  keys <- sprintf("p%d.s%d", u$path, u$scaffold)
  seqs <- stats::setNames(rep(NA_character_, nrow(u)), keys)
  labels <- paths_labels(con, id, u)
  feats <- list()
  for (i in seq_len(nrow(u))) {
    a <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, u$path[i], u$scaffold[i]))
    if (!nrow(a)) next
    seqs[i] <- a$sequence[1]
    f <- DBI::dbGetQuery(con, "SELECT element, start, stop, strand, partial_start, partial_stop
      FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(id, u$path[i], u$scaffold[i]))
    if (nrow(f)) {
      e <- partial_ends(f, nchar(a$sequence[1]), c(FALSE, FALSE))
      f$partial_start <- as.integer(e$p5); f$partial_stop <- as.integer(e$p3)
      feats[[keys[i]]] <- cbind(unit = keys[i], f, stringsAsFactors = FALSE)
    }
  }
  aln <- paths_align(seqs)
  if (is.null(aln)) return(NULL)
  feats <- if (length(feats)) do.call(rbind, feats) else
    data.frame(unit = character(), element = character(), start = numeric(), stop = numeric(),
      strand = character(), partial_start = numeric(), partial_stop = numeric())
  pl <- paths_payload(aln, feats, labels)
  pl$id <- "seqv_paths-canvas"; pl$unit <- id
  session$sendCustomMessage("rpseq", pl)
  pl
}

#' Summary line above the paths viewer
#' @noRd
paths_summary_text <- function(pl) {
  if (is.null(pl)) return("Could not align these contigs.")
  sprintf("%d contigs, %s columns, %s differ", length(pl$rows),
    format(pl$len, big.mark = ","), format(paths_diff_count(pl$rows), big.mark = ","))
}
