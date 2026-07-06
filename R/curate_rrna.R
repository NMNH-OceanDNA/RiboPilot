# Reference-guided curation of the nuclear rDNA feature model. Follows
# MitoPilot's design: a shared core (curate_rrna_core) driven by per-clade
# params (params_rrna.R), plus thin clade wrappers. The core is pure (no DB);
# curate_one() is the DB driver called by the WF1/WF2 Nextflow curate process,
# right after reconcile writes the annotations table.
#
# The core never base-aligns spacers. It checks the structural model (anchors
# present, ordered, one strand, sane lengths) and, when reference data is
# available, flags divergent anchors and contaminant lineages. Everything else
# is left to the human reviewer.

#' Curate one reconciled nrDNA feature model against a clade ruleset (pure)
#'
#' @param features data.frame(element, type, start, stop, length, strand, ...)
#'   as produced by [reconcile_features()] / stored in the `annotations` table.
#' @param params A ruleset from [rrna_params()].
#' @param ref_pid Optional named numeric vector (element -> percent identity to
#'   the top reference), from the `blast_ref_alignment` table. Missing elements
#'   are skipped (identity check is opportunistic).
#' @param sample_lineage,hit_lineage Optional lineage strings for the contaminant
#'   check (sample's assigned lineage vs. top BLAST-hit lineage).
#' @return list(features, flags, n_flags, notes) where `features` gains `status`
#'   ("ok"/"flag") and `flags` (comma-separated per-feature codes), and the
#'   sample-level `flags`/`notes` summarize missing anchors, ordering, etc.
#' @export
curate_rrna_core <- function(features, params = rrna_params(),
                             ref_pid = NULL, sample_lineage = NA, hit_lineage = NA) {
  sample_flags <- character(0)
  feat_flags <- list()
  add_feat_flag <- function(el, code) feat_flags[[el]] <<- c(feat_flags[[el]], code)

  has_feats <- !is.null(features) && nrow(features) > 0
  present <- if (has_feats) features$element else character(0)

  # 1. Anchors + spacers present ----
  for (g in RRNA_GENES) if (!(g %in% present)) sample_flags <- c(sample_flags, paste0("missing_", g))
  for (s in RRNA_SPACERS) if (!(s %in% present)) sample_flags <- c(sample_flags, paste0("missing_", s))

  if (has_feats) {
    ord <- order(features$start)
    fo <- features[ord, , drop = FALSE]

    # 2. Order: observed elements must follow the canonical 5'->3' sequence ----
    canon <- RRNA_ELEMENT_ORDER[RRNA_ELEMENT_ORDER %in% fo$element]
    if (!identical(fo$element, canon)) sample_flags <- c(sample_flags, "out_of_order")

    # 3. Overlaps between adjacent features ----
    if (nrow(fo) > 1) {
      for (i in 2:nrow(fo)) {
        if (!is.na(fo$start[i]) && !is.na(fo$stop[i - 1]) && fo$start[i] <= fo$stop[i - 1]) {
          add_feat_flag(fo$element[i], "overlap")
          sample_flags <- c(sample_flags, "overlap")
        }
      }
    }

    # 4. Strand: everything should be on '+' after orientation ----
    if ("strand" %in% names(features)) {
      for (i in seq_len(nrow(features))) {
        st <- features$strand[i]
        if (!is.na(st) && !identical(st, "+")) {
          add_feat_flag(features$element[i], "wrong_strand")
          sample_flags <- c(sample_flags, "wrong_strand")
        }
      }
    }

    # 5. Per-feature length tolerance ----
    for (i in seq_len(nrow(features))) {
      el <- features$element[i]
      rng <- params$lengths[[el]]
      len <- features$length[i]
      if (is.null(rng) || is.na(len)) next
      if (len < rng[1] || len > rng[2]) add_feat_flag(el, "length")
    }

    # 6. Reference identity of conserved anchors (opportunistic) ----
    if (!is.null(ref_pid) && length(ref_pid) > 0) {
      for (g in RRNA_GENES) {
        pid <- if (g %in% names(ref_pid)) ref_pid[[g]] else NA_real_
        if (!is.na(pid) && pid < params$min_anchor_pid) {
          add_feat_flag(g, "low_identity")
        }
      }
    }
  }

  # 7. Contaminant lineage check ----
  if (isTRUE(params$check_contaminant) && lineages_conflict(sample_lineage, hit_lineage)) {
    sample_flags <- c(sample_flags, "contaminant")
  }

  # Assemble per-feature status ----
  if (has_feats) {
    features$flags <- vapply(features$element, function(el) {
      fl <- feat_flags[[el]]
      if (is.null(fl)) "" else paste(unique(fl), collapse = ",")
    }, character(1))
    features$status <- ifelse(nzchar(features$flags), "flag", "ok")
  } else {
    features <- data.frame()
  }

  sample_flags <- unique(sample_flags)
  notes <- if (length(sample_flags)) paste(sample_flags, collapse = "; ") else NA_character_
  list(features = features, flags = sample_flags, n_flags = length(sample_flags), notes = notes)
}

#' Thin actinopterygii curation wrapper.
#' @param features,ref_pid,sample_lineage,hit_lineage see [curate_rrna_core()].
#' @return see [curate_rrna_core()].
#' @export
curate_actinopterygii_rrna <- function(features, ref_pid = NULL,
                                       sample_lineage = NA, hit_lineage = NA) {
  curate_rrna_core(features, rrna_params("actinopterygii"),
    ref_pid = ref_pid, sample_lineage = sample_lineage, hit_lineage = hit_lineage)
}

#' Do two lineage strings share no taxon in common (gross contaminant)?
#' Conservative: only TRUE when both are present and, after dropping ubiquitous
#' high-rank tokens, they share zero taxa. Catches wrong-kingdom hits (e.g. a
#' fungal hit on a fish) without flagging same-clade samples.
#' @noRd
LINEAGE_STOPWORDS <- c("cellular organisms", "eukaryota", "opisthokonta")
lineages_conflict <- function(a, b) {
  toks <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
    t <- tolower(trimws(unlist(strsplit(x, "[;,]"))))
    unique(setdiff(t[nchar(t) > 3], LINEAGE_STOPWORDS))
  }
  ta <- toks(a); tb <- toks(b)
  if (length(ta) == 0 || length(tb) == 0) return(FALSE)
  length(intersect(ta, tb)) == 0
}

#' Ensure the curation columns/tables exist (idempotent migration)
#'
#' Adds `curate_opts` (table + row), the `annotate.curate_opts`/`annotate.n_flags`
#' columns, and the `annotations.status`/`annotations.flags` columns to a project
#' DB created before curation existed. New projects get these in [new_db()].
#' @param con Open DBI connection.
#' @export
ensure_curate_schema <- function(con) {
  has_col <- function(tbl, col) col %in% DBI::dbListFields(con, tbl)
  has_tbl <- function(tbl) tbl %in% DBI::dbListTables(con)

  if (!has_tbl("curate_opts")) {
    DBI::dbExecute(con,
      "CREATE TABLE curate_opts (
        curate_opts TEXT NOT NULL,
        ruleset TEXT,
        min_anchor_pid REAL,
        check_contaminant INTEGER,
        PRIMARY KEY (curate_opts)
      );")
  }
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM curate_opts")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO curate_opts (curate_opts, ruleset, min_anchor_pid, check_contaminant)
       VALUES ('default', 'metazoa_generic', 85, 1)")
  }
  if (has_tbl("annotate")) {
    if (!has_col("annotate", "curate_opts")) {
      DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN curate_opts TEXT DEFAULT 'default'")
    }
    if (!has_col("annotate", "n_flags")) {
      DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN n_flags INTEGER")
    }
  }
  if (has_tbl("annotations")) {
    if (!has_col("annotations", "status")) {
      DBI::dbExecute(con, "ALTER TABLE annotations ADD COLUMN status TEXT")
    }
    if (!has_col("annotations", "flags")) {
      DBI::dbExecute(con, "ALTER TABLE annotations ADD COLUMN flags TEXT")
    }
  }
  invisible()
}

#' Curate one sample's annotations and write results back to the DB
#'
#' Reads the reconciled `annotations`, the sample's `curate_opts` ruleset, any
#' per-feature reference identities (`blast_ref_alignment`), and the sample vs.
#' hit lineages, runs [curate_rrna_core()], then updates each feature's
#' `status`/`flags` and the per-sample `annotate` summary (`n_flags`,
#' `annotate_notes`) in one transaction.
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param ruleset Optional ruleset override (else read from `curate_opts`).
#' @return (invisibly) the curate_rrna_core() result.
#' @export
curate_one <- function(con, id, ruleset = NULL) {
  ensure_curate_schema(con)

  feats <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop, length, strand FROM annotations WHERE ID = ?",
    params = list(id))

  # Ruleset: explicit arg > sample's curate_opts row > default.
  if (is.null(ruleset)) {
    ruleset <- tryCatch(DBI::dbGetQuery(con,
      "SELECT co.ruleset FROM annotate an JOIN curate_opts co ON an.curate_opts = co.curate_opts
       WHERE an.ID = ?", params = list(id))$ruleset[1], error = function(e) NA_character_)
  }
  params <- rrna_params(if (length(ruleset) && !is.na(ruleset)) ruleset else "metazoa_generic")

  # Per-feature reference identity (opportunistic; empty when BLAST/ref absent).
  ref_pid <- NULL
  rp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT element, pid FROM blast_ref_alignment WHERE ID = ?", params = list(id)),
    error = function(e) data.frame())
  if (nrow(rp) > 0) ref_pid <- stats::setNames(rp$pid, rp$element)

  sample_lineage <- DBI::dbGetQuery(con, "SELECT lineage FROM samples WHERE ID = ?",
    params = list(id))$lineage[1]
  hit_lineage <- tryCatch(DBI::dbGetQuery(con,
    "SELECT blast_lineage FROM path_scores WHERE ID = ? AND auto_select = 1",
    params = list(id))$blast_lineage[1], error = function(e) NA_character_)

  res <- curate_rrna_core(feats, params, ref_pid = ref_pid,
    sample_lineage = sample_lineage, hit_lineage = hit_lineage)

  ts <- as.integer(Sys.time())
  DBI::dbWithTransaction(con, {
    if (nrow(res$features) > 0) {
      for (i in seq_len(nrow(res$features))) {
        DBI::dbExecute(con,
          "UPDATE annotations SET status = ?, flags = ? WHERE ID = ? AND element = ?",
          params = list(res$features$status[i],
            if (nzchar(res$features$flags[i])) res$features$flags[i] else NA_character_,
            id, res$features$element[i]))
      }
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET n_flags = ?, annotate_notes = ?, time_stamp = ? WHERE ID = ?",
      params = list(res$n_flags, res$notes, ts, id))
  })
  invisible(res)
}
