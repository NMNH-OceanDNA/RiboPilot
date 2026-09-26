#!/usr/bin/env Rscript
# Pick the nuclear rDNA contigs out of a user-supplied assembly.
# Args: <db> <id> <assembly.fasta> <out_prefix> <chunk_size> [chunk.gff ...]
# Writes <out_prefix>_assembly_1.fasta (assembler header contract), the merged
# scan rrna_scan.gff, the rrna_candidates evidence, an assemble note, and
# status.txt (ok|empty).
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 5, "rp_find_rrna.R <db> <id> <assembly.fasta> <out_prefix> <chunk_size> [chunk.gff ...]")

local({
  db <- a[1]; id <- a[2]; fa <- a[3]; prefix <- a[4]; gffs <- a[-(1:5)]
  con <- open_project_db(db)
  on.exit(DBI::dbDisconnect(con))

  seqs <- Biostrings::readDNAStringSet(fa)
  # barrnap reports the first whitespace token of the header as seqid
  contigs <- sub("\\s.*$", "", names(seqs))
  lengths <- stats::setNames(Biostrings::width(seqs), contigs)

  cap <- tryCatch(DBI::dbGetQuery(con,
    "SELECT ao.max_scaffolds FROM assemble a
       JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(id))$max_scaffolds[1], error = function(e) NA_integer_)
  if (length(cap) != 1 || is.na(cap)) cap <- 10L

  gff_lines <- as.character(unlist(lapply(gffs, readLines, warn = FALSE)))
  writeLines(gff_lines, "rrna_scan.gff")
  lost <- ceiling(length(seqs) / as.integer(a[5])) - length(gffs)
  res <- select_rrna_contigs(gff_lines, lengths = lengths, max_contigs = cap)
  write_rrna_candidates(con, id, res$evidence)

  ts <- as.integer(Sys.time())
  if (lost > 0 || length(res$candidates) == 0) {
    n_hit <- nrow(res$evidence)
    note <- if (lost > 0) {
      sprintf("rRNA search incomplete: %d of %d part(s) of the supplied assembly could not be scanned",
        lost, lost + length(gffs))
    } else if (n_hit == 0) {
      "no rRNA genes found in the supplied assembly"
    } else {
      sprintf("no rDNA contig confirmed: %d contig(s) had rRNA hits but no 18S or 28S anchor", n_hit)
    }
    # only find_rrna_notes: the app shows the two note columns as one cell, so
    # writing the same sentence to both printed it twice. Nothing assembled
    # here, so any assemble_notes left by an earlier run is stale.
    DBI::dbExecute(con,
      "UPDATE assemble SET paths=NULL, scaffolds=NULL, length=NULL, topology=NULL,
       assemble_switch=3, assemble_notes=NULL, find_rrna_notes=?, time_stamp=? WHERE ID=?",
      params = list(note, ts, id))
    # a re-run that found nothing must not leave last run's units queued
    DBI::dbExecute(con,
      "UPDATE annotate SET annotate_switch = 0
       WHERE ID = ? AND COALESCE(annotate_lock,0) = 0 AND COALESCE(reviewed,'no') = 'no'",
      params = list(id))
    writeLines("empty", "status.txt")
    cat("no rDNA contigs:", id, "\n")
    return(invisible())
  }

  # Same header contract format_assembly.sh writes, so load_assemblies() and
  # everything downstream treat this exactly like an assembled path.
  keep <- match(res$candidates, contigs)
  out <- seqs[keep]
  names(out) <- sprintf("%s.1.%d linear %s", id, seq_along(keep), res$candidates)
  Biostrings::writeXStringSet(out, sprintf("%s_assembly_1.fasta", prefix))

  # its own column: load_assemblies() overwrites assemble_notes right after this
  DBI::dbExecute(con,
    "UPDATE assemble SET find_rrna_notes = ?, time_stamp = ? WHERE ID = ?",
    params = list(sprintf("%d rDNA contig(s) found in the supplied assembly",
      length(res$candidates)), ts, id))
  writeLines("ok", "status.txt")
  cat("rDNA contigs found:", id, length(res$candidates), "\n")
})
