# Wire real assembly output into the scaffold/path scorer.
#
# Turns a set of GetOrganelle candidate path fastas + the assembly graph into a
# scored candidate table by extracting, per path: length, graph depth, anchor
# completeness (barrnap/Rfam), and the top BLAST hit + its resolved lineage.
# Then calls score_candidates().

#' Distinct rRNA anchors (18S/5.8S/28S) present in a fasta, via barrnap
#'
#' Uses the bioconda barrnap fork's eukaryote Rfam models (`--kingdom fun`;
#' SSU/LSU/5_8S covariance models). Returns the distinct anchor names found.
#'
#' @param fasta Path to a nucleotide fasta.
#' @param kingdom barrnap kingdom (default "fun" = eukaryote Rfam models here).
#' @return character vector subset of RRNA_GENES.
#' @noRd
assembly_anchors <- function(fasta, kingdom = "fun") {
  bin <- unname(Sys.which("barrnap"))
  if (!nzchar(bin)) {
    return(character(0))
  }
  gff <- suppressWarnings(system2(bin, c("--kingdom", kingdom, "--quiet", shQuote(fasta)),
    stdout = TRUE, stderr = FALSE))
  gff <- gff[!startsWith(gff, "#")]
  names <- unique(sub(".*Name=([^;]+).*", "\\1", gff[grepl("Name=", gff)]))
  # barrnap names like "18S_rRNA", "28S_rRNA", "5_8S_rRNA"
  map <- c("18S_rRNA" = "18S", "5_8S_rRNA" = "5.8S", "5.8S_rRNA" = "5.8S", "28S_rRNA" = "28S")
  out <- unname(map[names])
  intersect(RRNA_GENES, stats::na.omit(out))
}

#' Best BLAST hit per query against a nucleotide db
#'
#' @param query_fasta Path to query fasta.
#' @param db BLAST nucleotide db prefix.
#' @param with_staxids If TRUE, request `staxids` (nt-style dbs carry taxids);
#'   otherwise the subject id is used to derive taxonomy (custom seed DBs whose
#'   ids are `ACCESSION_Genus_species`).
#' @return data.frame(qseqid, sseqid, pident, bitscore, staxid) best-per-query,
#'   or empty.
#' @noRd
blast_best_hits <- function(query_fasta, db, with_staxids = FALSE) {
  bin <- unname(Sys.which("blastn"))
  if (!nzchar(bin)) {
    stop("blastn not found (see environment.yml).", call. = FALSE)
  }
  fmt <- if (with_staxids) "6 qseqid sseqid pident bitscore staxids" else "6 qseqid sseqid pident bitscore"
  out <- suppressWarnings(system2(bin, c(
    "-query", shQuote(query_fasta), "-db", shQuote(db),
    "-max_target_seqs", "5", "-evalue", "1e-10", "-outfmt", shQuote(fmt)
  ), stdout = TRUE, stderr = FALSE))
  parse_blast_best(out, with_staxids)
}

#' Parse blastn tabular output down to the best hit per query
#' @noRd
parse_blast_best <- function(out, with_staxids = FALSE) {
  if (length(out) == 0) {
    return(data.frame())
  }
  f <- strsplit(out, "\t")
  df <- data.frame(
    qseqid = vapply(f, `[`, character(1), 1),
    sseqid = vapply(f, `[`, character(1), 2),
    pident = as.numeric(vapply(f, `[`, character(1), 3)),
    bitscore = as.numeric(vapply(f, `[`, character(1), 4)),
    staxid = if (with_staxids) sub(";.*", "", vapply(f, function(x) x[5] %||% NA, character(1))) else NA_character_,
    stringsAsFactors = FALSE
  )
  # best hit per query by bitscore
  df <- df[order(df$qseqid, -df$bitscore), ]
  df[!duplicated(df$qseqid), ]
}

#' Derive "Genus species" from a custom seed-DB subject id (ACC_Genus_species)
#' @noRd
sseqid_species <- function(sseqid) {
  sp <- sub("^[^_]+_", "", sseqid) # drop accession
  gsub("_", " ", sp)
}

#' Score GetOrganelle candidate paths from real assembly output
#'
#' @param path_fastas Character vector of per-path fasta files.
#' @param gfa_path The assembly graph `.gfa` (for depth). May be NA.
#' @param sample_lineage The sample's resolved rank-aligned lineage.
#' @param sample_taxids The sample's rank-aligned lineage taxids. Preferred over
#'   the names when present: taxids carry no homonyms.
#' @param blast_db BLAST nucleotide db prefix (e.g. the Actinopterygii seeds db,
#'   or an nt-style db).
#' @param taxdump_dir Pinned taxdump dir (for resolving hit lineages).
#' @param with_staxids TRUE if `blast_db` carries taxids; else taxonomy is taken
#'   from the seed-DB subject id species.
#' @param expected_length,weights Passed to [score_candidates()].
#' @return Scored candidate data.frame (see [score_candidates()]).
#' @export
score_assembly <- function(path_fastas,
                           gfa_path = NA_character_,
                           sample_lineage = NA_character_,
                           sample_taxids = NA_character_,
                           blast_db = NULL,
                           taxdump_dir = NULL,
                           with_staxids = FALSE,
                           expected_length = c(4000, 12000),
                           weights = DEFAULT_SCORE_WEIGHTS) {
  gd <- if (!is.na(gfa_path) && file.exists(gfa_path)) parse_gfa_depth(gfa_path) else NULL

  per_path <- lapply(seq_along(path_fastas), function(i) {
    fa <- path_fastas[i]
    seqs <- readLines(fa, warn = FALSE)
    headers <- sub("^>", "", seqs[startsWith(seqs, ">")])
    len <- sum(nchar(seqs[!startsWith(seqs, ">")]))
    depth <- if (!is.null(gd) && length(headers) > 0) {
      mean(vapply(headers, scaffold_depth, numeric(1), gfa_depth = gd), na.rm = TRUE)
    } else {
      NA_real_
    }
    anchors <- paste(assembly_anchors(fa), collapse = ",")
    bh <- if (!is.null(blast_db)) blast_best_hits(fa, blast_db, with_staxids) else data.frame()
    top <- if (nrow(bh) > 0) bh[which.max(bh$bitscore), ] else NULL
    data.frame(
      path = i,
      length = len,
      depth = depth,
      anchors = anchors,
      blast_pident = if (!is.null(top)) top$pident else NA_real_,
      hit_id = if (!is.null(top)) top$sseqid else NA_character_,
      hit_staxid = if (!is.null(top)) top$staxid else NA_character_,
      hit_species = if (!is.null(top) && !with_staxids) sseqid_species(top$sseqid) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  cand <- do.call(rbind, per_path)

  # resolve top-hit lineages
  cand$blast_lineage <- NA_character_
  cand$blast_lineage_taxids <- NA_character_
  if (with_staxids && have_taxdump(taxdump_dir) && any(!is.na(cand$hit_staxid))) {
    info <- tk_taxid_info(as.integer(stats::na.omit(cand$hit_staxid)),
      if (file.exists(file.path(taxdump_dir, "names.dmp"))) taxdump_dir else extract_taxdump(taxdump_dir))
    i <- match(cand$hit_staxid, info$taxid)
    cand$blast_lineage <- info$lineage[i]
    cand$blast_lineage_taxids <- info$lineage_taxids[i]
  } else if (!with_staxids && have_taxdump(taxdump_dir) && any(!is.na(cand$hit_species))) {
    lin <- resolve_taxonomy(stats::na.omit(unique(cand$hit_species)), taxdump_dir = taxdump_dir)
    i <- match(cand$hit_species, lin$Taxon)
    cand$blast_lineage <- lin$lineage[i]
    cand$blast_lineage_taxids <- lin$lineage_taxids[i]
  }

  score_candidates(cand, sample_lineage = sample_lineage, sample_taxids = sample_taxids,
    expected_length = expected_length, weights = weights)
}
