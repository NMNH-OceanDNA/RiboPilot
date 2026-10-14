# Finding nuclear rDNA contigs inside a user-supplied assembly.
# See dev/design/20_USER_ASSEMBLIES.md.

#' Select nuclear rDNA contigs from a barrnap scan
#'
#' Given the merged barrnap GFF for one sample's supplied assembly, keep the
#' contigs that carry a nuclear rDNA anchor. A contig qualifies on an 18S or a
#' 28S hit; a lone 5.8S hit does not, because 5.8S is short enough that a single
#' hit is as likely to be a degraded pseudogene copy as a real repeat unit.
#'
#' Only contigs barrnap hit at all are scored, so the evidence returned is the
#' handful of rRNA-bearing contigs rather than every contig in a draft genome.
#'
#' @param gff_lines Character vector, the merged barrnap GFF.
#' @param lengths Optional named numeric of contig lengths, names = contig ids.
#' @param max_contigs Most contigs to keep, best first (default 10). Fed from
#'   the sample's `assemble_opts.max_scaffolds`, which is also the gate
#'   [load_assemblies()] applies afterwards.
#' @return list with `candidates` (contig ids, best first) and `evidence` (one
#'   row per contig scored: contig, length, rrna_features, rrna_bases, rank,
#'   selected, reason)
#' @export
select_rrna_contigs <- function(gff_lines, lengths = NULL, max_contigs = 10) {
  empty <- data.frame(
    contig = character(0), length = integer(0), rrna_features = character(0),
    rrna_bases = integer(0), rank = integer(0), selected = integer(0),
    reason = character(0), stringsAsFactors = FALSE
  )
  feats <- if (length(gff_lines)) parse_barrnap(gff_lines) else data.frame()
  if (nrow(feats) == 0) {
    return(list(candidates = character(0), evidence = empty))
  }

  feats$bases <- pmax(feats$stop - feats$start + 1L, 0L)
  by_contig <- split(feats, feats$seqid)
  ev <- do.call(rbind, lapply(names(by_contig), function(ctg) {
    f <- by_contig[[ctg]]
    el <- unique(f$element)
    data.frame(
      contig = ctg,
      length = if (!is.null(lengths) && ctg %in% names(lengths)) {
        as.integer(lengths[[ctg]])
      } else {
        NA_integer_
      },
      rrna_features = paste(el[order(match(el, c("18S", "5.8S", "28S")))], collapse = ","),
      rrna_bases = as.integer(sum(f$bases)),
      has_anchor = any(el %in% c("18S", "28S")),
      stringsAsFactors = FALSE
    )
  }))

  # best first: most rRNA sequence, then longest contig
  len_key <- ev$length
  len_key[is.na(len_key)] <- 0L
  ev <- ev[order(-ev$rrna_bases, -len_key), , drop = FALSE]
  ev$rank <- seq_len(nrow(ev))
  ev$selected <- 0L
  ev$reason <- NA_character_
  ev$reason[!ev$has_anchor] <- "no 18S or 28S anchor"

  ok <- which(ev$has_anchor)
  cap <- if (is.na(max_contigs) || max_contigs <= 0) length(ok) else max_contigs
  keep <- utils::head(ok, cap)
  ev$selected[keep] <- 1L
  if (length(ok) > length(keep)) {
    ev$reason[setdiff(ok, keep)] <- sprintf("beyond the %d best rDNA contigs", as.integer(cap))
  }
  ev$has_anchor <- NULL
  rownames(ev) <- NULL

  list(candidates = ev$contig[ev$selected == 1L], evidence = ev)
}

#' Write the rDNA contig search evidence for one sample
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param evidence data.frame from [select_rrna_contigs()].
#' @return (invisibly) rows written.
#' @export
write_rrna_candidates <- function(con, id, evidence) {
  DBI::dbExecute(con, "DELETE FROM rrna_candidates WHERE ID = ?", params = list(id))
  if (is.null(evidence) || nrow(evidence) == 0) {
    return(invisible(0L))
  }
  evidence$ID <- id
  evidence$time_stamp <- as.integer(Sys.time())
  DBI::dbAppendTable(con, "rrna_candidates", evidence[, c(
    "ID", "contig", "length", "rrna_features", "rrna_bases", "rank",
    "selected", "reason", "time_stamp"
  )])
  invisible(nrow(evidence))
}
