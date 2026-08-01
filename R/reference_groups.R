# Turning frozen sample taxonomy into clade reference-database groups.
#
# A seed database is only useful in a narrow band: too narrow and GenBank has
# nothing to seed from, too broad and the seeds stop being specific to the
# organism being assembled. The walk therefore climbs each sample's lineage
# and stops at the first clade that clears a record floor, even when that
# clade is itself over the ceiling; the ceiling instead governs the fallback
# taken when nothing ever clears the floor, and which merges into an ancestor
# clade are allowed. All of this is pure given a count function, which is
# what makes it testable without touching NCBI.

#' Rank-aligned ladders for every sample with a usable lineage
#'
#' One row per resolved rank per sample, kingdom first and species last. Ranks
#' the taxonomy could not resolve carry no taxid and are dropped, so a sample
#' with an unassigned class simply has no class rung.
#'
#' @param con Open project DB connection.
#' @return data.frame(ID, level, rank, taxid, name).
#' @noRd
sample_ladders <- function(con) {
  empty <- data.frame(ID = character(0), level = integer(0), rank = character(0),
    taxid = integer(0), name = character(0), stringsAsFactors = FALSE)
  s <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ID, lineage, lineage_taxids FROM samples
     WHERE lineage_taxids IS NOT NULL AND TRIM(lineage_taxids) != ''"),
    error = function(e) NULL)
  if (is.null(s) || nrow(s) == 0) return(empty)
  rows <- lapply(seq_len(nrow(s)), function(i) {
    tx <- parse_rank_taxids(s$lineage_taxids[i])
    nm <- parse_rank_names(s$lineage[i])
    keep <- which(!is.na(tx) & !is.na(nm))
    if (length(keep) == 0) return(NULL)
    data.frame(ID = s$ID[i], level = as.integer(keep), rank = RRNA_LINEAGE_RANKS[keep],
      taxid = unname(tx[keep]), name = unname(nm[keep]), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Filename-safe group key from a clade name
#' @noRd
ref_id_from_name <- function(name) {
  x <- tolower(trimws(name %||% ""))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

#' Sample IDs whose lineage passes through a taxid
#' @param ladders From [sample_ladders()]. @param taxid NCBI taxid.
#' @noRd
ladder_coverage <- function(ladders, taxid) {
  unique(ladders$ID[ladders$taxid == as.integer(taxid)])
}

#' The lineage shared by a set of samples, with coverage and optional counts
#'
#' @param ladders From [sample_ladders()].
#' @param ids Sample IDs in the group.
#' @param count_fn Optional function(taxid, name) -> integer record count.
#' @return data.frame(level, rank, taxid, name, n_samples, n_records), deepest last.
#' @noRd
group_ladder <- function(ladders, ids, count_fn = NULL) {
  l <- ladders[ladders$ID %in% ids, , drop = FALSE]
  empty <- data.frame(level = integer(0), rank = character(0), taxid = integer(0),
    name = character(0), n_samples = integer(0), n_records = integer(0),
    stringsAsFactors = FALSE)
  if (nrow(l) == 0) return(empty)
  n_ids <- length(unique(l$ID))
  # a rung belongs to the shared ladder only when every member passes through
  # it; a rung only some members reach is not part of the group's shared path
  agg <- stats::aggregate(list(n_samples = l$ID), by = list(taxid = l$taxid),
    FUN = function(x) length(unique(x)))
  key <- l[!duplicated(l$taxid), c("level", "rank", "taxid", "name")]
  out <- merge(key, agg, by = "taxid")
  out <- out[out$n_samples == n_ids, , drop = FALSE]
  out <- out[order(out$level), , drop = FALSE]
  if (nrow(out) == 0) {
    out$n_records <- integer(0)
    rownames(out) <- NULL
    return(out[, c("level", "rank", "taxid", "name", "n_samples", "n_records"), drop = FALSE])
  }
  out$n_records <- NA_integer_
  if (!is.null(count_fn)) {
    out$n_records <- vapply(seq_len(nrow(out)),
      function(i) as.integer(count_fn(out$taxid[i], out$name[i]) %|NA|% NA_integer_), integer(1))
  }
  rownames(out) <- NULL
  out[, c("level", "rank", "taxid", "name", "n_samples", "n_records"), drop = FALSE]
}

#' Walk one sample's ladder up to the first clade clearing the record floor
#'
#' Returns the deepest rung at or above `start_rank` whose count reaches
#' `min_records`, and stops there even when that count is over `max_records`:
#' once a clade has enough evidence, climbing further only makes the seed
#' database less specific. If no rung ever reaches `min_records`, the deepest
#' rung whose count stays at or under `max_records` is kept instead, so the
#' fallback never lands on a clade that is far too broad; if every rung is
#' over the ceiling, the last rung considered is returned anyway.
#' @return one-row data.frame(level, rank, taxid, name, n_records), or NULL.
#' @noRd
walk_to_clade <- function(rungs, count_fn, min_records, max_records, start_rank) {
  start <- match(start_rank, RRNA_LINEAGE_RANKS)
  cand <- rungs[rungs$level <= start, , drop = FALSE]
  cand <- cand[order(-cand$level), , drop = FALSE]   # deepest first, walking up
  if (nrow(cand) == 0) return(NULL)
  last <- NULL
  last_any <- NULL
  for (i in seq_len(nrow(cand))) {
    n <- suppressWarnings(as.integer(count_fn(cand$taxid[i], cand$name[i])))
    if (is.na(n)) next
    row <- data.frame(level = cand$level[i], rank = cand$rank[i], taxid = cand$taxid[i],
      name = cand$name[i], n_records = n, stringsAsFactors = FALSE)
    last_any <- row
    if (n <= max_records) last <- row
    if (n >= min_records) return(row)
  }
  last %||% last_any
}

#' Suggest clade groups for a project's samples
#'
#' @param ladders From [sample_ladders()].
#' @param count_fn function(taxid, name) -> integer GenBank record count.
#' @param min_records,max_records Record window bounding the walk.
#' @param start_rank Rank the walk starts from (default "family").
#' @return list(groups = data.frame(ref_id, clade_taxid, clade_name, clade_rank,
#'   n_records, n_samples), members = data.frame(ID, ref_id)).
#' @noRd
suggest_reference_groups <- function(ladders, count_fn, min_records = 1000,
                                     max_records = 25000, start_rank = "family") {
  ids <- unique(ladders$ID)
  empty <- list(
    groups = data.frame(ref_id = character(0), clade_taxid = integer(0),
      clade_name = character(0), clade_rank = character(0), n_records = integer(0),
      n_samples = integer(0), stringsAsFactors = FALSE),
    members = data.frame(ID = character(0), ref_id = character(0), stringsAsFactors = FALSE))
  if (length(ids) == 0) return(empty)

  picks <- lapply(ids, function(id) {
    p <- walk_to_clade(ladders[ladders$ID == id, , drop = FALSE], count_fn,
      min_records, max_records, start_rank)
    if (is.null(p)) return(NULL)
    p$ID <- id
    p
  })
  picks <- Filter(Negate(is.null), picks)
  if (length(picks) == 0) return(empty)
  picks <- do.call(rbind, picks)

  # Fold a group into an ancestor group when one exists and stays under the
  # ceiling: two fish families are better served by one class database than by
  # two family databases that each miss the other's reads. Only clades that
  # some sample's own walk already landed on are candidates, so this never
  # invents a new clade or re-queries count_fn. A fold can itself land on a
  # clade that has since folded further (A -> order because B picked order,
  # then B -> class because C picked class), so the pass repeats to a fixed
  # point rather than running once; level strictly decreases on every change,
  # so this always terminates.
  repeat {
    chosen <- unique(picks[, c("level", "rank", "taxid", "name", "n_records")])
    changed <- FALSE
    for (i in seq_len(nrow(picks))) {
      anc <- ladders$taxid[ladders$ID == picks$ID[i]]
      up <- chosen[chosen$taxid %in% anc & chosen$level < picks$level[i] &
          chosen$n_records <= max_records, , drop = FALSE]
      if (nrow(up) > 0) {
        up <- up[order(-up$level), , drop = FALSE][1, ]
        if (!identical(picks$taxid[i], up$taxid)) changed <- TRUE
        picks[i, c("level", "rank", "taxid", "name", "n_records")] <- up
      }
    }
    if (!changed) break
  }

  # Groups are keyed on clade taxid, not on the derived slug: two distinct
  # clades can slugify to the same ref_id ("Foo bar" / "Foo-bar"), and
  # collapsing on the slug would silently merge unrelated groups.
  slug <- ref_id_from_name(picks$name)
  uniq <- unique(data.frame(taxid = picks$taxid, slug = slug, stringsAsFactors = FALSE))
  dup <- uniq$slug[duplicated(uniq$slug) | duplicated(uniq$slug, fromLast = TRUE)]
  uniq$ref_id <- ifelse(uniq$slug %in% dup, paste(uniq$slug, uniq$taxid, sep = "_"), uniq$slug)
  picks$ref_id <- uniq$ref_id[match(picks$taxid, uniq$taxid)]

  members <- picks[, c("ID", "ref_id")]
  rownames(members) <- NULL
  g <- picks[!duplicated(picks$taxid), , drop = FALSE]
  groups <- data.frame(
    ref_id = g$ref_id, clade_taxid = as.integer(g$taxid), clade_name = g$name,
    clade_rank = g$rank, n_records = as.integer(g$n_records),
    n_samples = as.integer(table(members$ref_id)[g$ref_id]),
    stringsAsFactors = FALSE)
  rownames(groups) <- NULL
  list(groups = groups, members = members)
}
