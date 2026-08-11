# Scoring driver: after assembly, score each sample's candidate paths and write
# the results (with the auto-selected path) to the path_scores table.
#
# Degrades gracefully: if the project has no pinned taxdump (so BLAST-hit
# lineages can't be resolved) or the sample's lineage is unset, the taxonomy
# term is simply 0 and the other terms (identity, anchors, depth, length) still
# rank the candidates.

#' Resolve a seeds/BLAST reference to a local file (download a URL to temp)
#' @noRd
local_ref <- function(ref) {
  if (is.na(ref) || !nzchar(ref)) {
    return(NA_character_)
  }
  if (grepl("^https?://", ref)) {
    dest <- tempfile(fileext = ".fasta")
    if (utils::download.file(ref, dest, quiet = TRUE, mode = "wb") != 0) {
      return(NA_character_)
    }
    return(dest)
  }
  if (file.exists(ref)) ref else NA_character_
}

#' Build a nucleotide BLAST db from a fasta (returns db prefix or NA)
#' @noRd
make_blast_db <- function(fasta) {
  if (is.na(fasta) || !nzchar(unname(Sys.which("makeblastdb")))) {
    return(NA_character_)
  }
  db <- tempfile("rp_scoredb_")
  ok <- system2("makeblastdb", c("-in", shQuote(fasta), "-dbtype", "nucl", "-out", shQuote(db)),
    stdout = FALSE, stderr = FALSE)
  if (ok == 0) db else NA_character_
}

#' Score assembled candidate paths for a project and write `path_scores`
#'
#' @param path Project directory.
#' @param out_dir Pipeline output dir (default `<path>/out`).
#' @param taxdump_dir Pinned taxdump (default `<path>/.taxonomy`); optional.
#' @param blast_db Optional prebuilt BLAST db prefix for taxonomy of hits. If
#'   NULL, one is built from each sample's `seeds_db` (with species-labeled
#'   headers).
#' @param con Optional open DBI connection (else opened from `<path>/.sqlite`).
#' @return (invisibly) NULL.
#' @noRd
score_project <- function(path = ".", out_dir = file.path(path, "out"),
                          taxdump_dir = file.path(path, ".taxonomy"),
                          blast_db = NULL, con = NULL) {
  if (is.null(con)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
    on.exit(DBI::dbDisconnect(con), add = TRUE)
  }
  td <- if (have_taxdump(taxdump_dir)) taxdump_dir else NULL

  info <- DBI::dbGetQuery(con, "
    SELECT a.ID, s.lineage, s.lineage_taxids, ao.seeds_db
    FROM assemble a
    JOIN samples s        ON a.ID = s.ID
    JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
    WHERE EXISTS (SELECT 1 FROM assemblies x WHERE x.ID = a.ID)")
  if (nrow(info) == 0) {
    message("No assembled samples to score.")
    return(invisible())
  }

  # cache one BLAST db per distinct seeds_db
  db_cache <- list()
  get_db <- function(seeds) {
    if (!is.null(blast_db)) {
      return(blast_db)
    }
    key <- seeds %||% ""
    if (is.null(db_cache[[key]])) {
      db_cache[[key]] <<- make_blast_db(local_ref(seeds))
    }
    db_cache[[key]]
  }

  for (i in seq_len(nrow(info))) {
    id <- info$ID[i]
    adir <- file.path(out_dir, id, "assemble")
    fastas <- Sys.glob(file.path(adir, paste0(id, "_assembly_*.fasta")))
    fastas <- fastas[!grepl("_assembly_0\\.fasta$", fastas)]
    if (length(fastas) == 0) next
    gfa <- Sys.glob(file.path(adir, "*.gfa"))
    scored <- score_assembly(
      path_fastas = fastas,
      gfa_path = if (length(gfa) > 0) gfa[1] else NA_character_,
      sample_lineage = info$lineage[i] %||% NA_character_,
      sample_taxids = info$lineage_taxids[i] %||% NA_character_,
      blast_db = get_db(info$seeds_db[i]),
      taxdump_dir = td,
      with_staxids = FALSE
    )
    write_path_scores(con, id, scored)
    message("Scored ", id, ": ", nrow(scored), " path(s), auto-select = path ",
      scored$path[scored$auto_select][1])
  }
  invisible()
}
