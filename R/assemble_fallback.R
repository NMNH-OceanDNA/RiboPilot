# Assembly fallbacks for samples where GetOrganelle (B1) fails the QC gate.
#
# Option B2 (reference-guided baiting): map reads to the clade reference rRNA
# with bowtie2, keep the concordantly-aligned pairs, assemble that baited pool
# with SPAdes, then keep the contigs that BLAST to the label DB (true rDNA).
# Option A (whole-assembly cross-check): skip baiting and assemble all reads.
#
# Both emit the same `<id>_assembly_<n>.fasta` header format the GetOrganelle
# path produces (">id.path.scaffold topology ..."), which is what lets this be a
# per-sample assembler swap inside WF1 rather than a workflow of its own: the
# whole tail (load -> score -> annotate -> reconcile -> curate -> screen ->
# triage) runs unchanged.

#' Bait reads by mapping to a reference with bowtie2 (concordant pairs only)
#'
#' @param r1,r2 Paired FASTQ(.gz) paths.
#' @param reference Reference rRNA FASTA (clade seeds / top hit).
#' @param out_prefix Prefix for the baited-read output pair.
#' @param threads Threads for bowtie2.
#' @return Character vector `c(baited_R1, baited_R2)`, or NULL if baiting failed
#'   or recruited no pairs.
#' @noRd
bait_reads <- function(r1, r2, reference, out_prefix, threads = 4) {
  if (!nzchar(unname(Sys.which("bowtie2"))) || !nzchar(unname(Sys.which("bowtie2-build")))) {
    return(NULL)
  }
  if (is.null(reference) || is.na(reference) || !file.exists(reference)) {
    return(NULL)
  }
  idx <- tempfile("rp_baitidx_")
  if (system2("bowtie2-build", c("--quiet", "--threads", threads, shQuote(reference), shQuote(idx)),
    stdout = FALSE, stderr = FALSE) != 0) {
    return(NULL)
  }
  ok <- system2("bowtie2", c("-x", shQuote(idx), "-1", shQuote(r1), "-2", shQuote(r2),
    "--no-unal", "--al-conc-gz", shQuote(paste0(out_prefix, ".%.fastq.gz")),
    "--local", "--threads", threads, "-S", shQuote(tempfile(fileext = ".sam"))),
    stdout = FALSE, stderr = FALSE)
  b1 <- paste0(out_prefix, ".1.fastq.gz"); b2 <- paste0(out_prefix, ".2.fastq.gz")
  if (ok != 0 || !file.exists(b1) || !file.exists(b2) || file.info(b1)$size < 50) {
    return(NULL)
  }
  c(b1, b2)
}

#' Assemble a read pair with SPAdes; return the contigs FASTA path (or NA)
#' @noRd
run_spades <- function(r1, r2, out_dir, threads = 4) {
  if (!nzchar(unname(Sys.which("spades.py")))) {
    return(NA_character_)
  }
  ok <- system2("spades.py", c("--isolate", "-1", shQuote(r1), "-2", shQuote(r2),
    "-o", shQuote(out_dir), "-t", threads), stdout = FALSE, stderr = FALSE)
  contigs <- file.path(out_dir, "contigs.fasta")
  # --isolate can reject very small/odd libraries; retry with the default pipeline
  if (ok != 0 || !file.exists(contigs)) {
    system2("spades.py", c("-1", shQuote(r1), "-2", shQuote(r2),
      "-o", shQuote(out_dir), "-t", threads), stdout = FALSE, stderr = FALSE)
  }
  if (file.exists(contigs)) contigs else NA_character_
}

#' Keep the rDNA contigs (BLAST hit vs the label DB) and write assembly fastas
#'
#' Filters `contigs_fa` to sequences long enough and with a BLAST hit against
#' `labels_db`, then writes them as one candidate path in the standard
#' `<id>_assembly_1.fasta` header format (scaffolds ordered by length).
#'
#' @param contigs_fa Assembler contigs FASTA.
#' @param labels_db Label DB FASTA (single-locus rRNA records).
#' @param id Sample ID.
#' @param out_fa Output `<id>_assembly_1.fasta` path.
#' @param min_contig Minimum contig length to consider.
#' @param min_pident,max_evalue BLAST thresholds to call a contig rDNA.
#' @return Number of contigs kept (0 if none).
#' @export
extract_rdna_contigs <- function(contigs_fa, labels_db, id, out_fa,
                                 min_contig = 300, min_pident = 70, max_evalue = 1e-10) {
  if (is.na(contigs_fa) || !file.exists(contigs_fa)) {
    return(0L)
  }
  seqs <- Biostrings::readDNAStringSet(contigs_fa)
  seqs <- seqs[Biostrings::width(seqs) >= min_contig]
  if (length(seqs) == 0) {
    return(0L)
  }
  db <- make_blast_db(local_ref(labels_db))
  keep_names <- names(seqs)
  if (!is.na(db) && nzchar(unname(Sys.which("blastn")))) {
    qfa <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(seqs, qfa)
    out <- suppressWarnings(system2("blastn", c("-query", shQuote(qfa), "-db", shQuote(db),
      "-outfmt", shQuote("6 qseqid pident evalue"), "-max_target_seqs", "1", "-evalue", "1e-5"),
      stdout = TRUE, stderr = FALSE))
    hits <- if (length(out)) utils::read.table(text = out, sep = "\t",
      col.names = c("qseqid", "pident", "evalue"), stringsAsFactors = FALSE) else data.frame()
    if (nrow(hits) > 0) {
      hits <- hits[hits$pident >= min_pident & hits$evalue <= max_evalue, ]
      keep_names <- intersect(names(seqs), unique(hits$qseqid))
    } else {
      keep_names <- character(0)
    }
  }
  seqs <- seqs[names(seqs) %in% keep_names]
  if (length(seqs) == 0) {
    return(0L)
  }
  seqs <- seqs[order(Biostrings::width(seqs), decreasing = TRUE)]
  lines <- character(0)
  for (k in seq_along(seqs)) {
    lines <- c(lines,
      sprintf(">%s.1.%d linear %s", id, k, names(seqs)[k]),
      as.character(seqs[[k]]))
  }
  writeLines(lines, out_fa)
  length(seqs)
}

#' Run an assembly fallback for one sample (Option B2 baiting, or Option A)
#'
#' Writes `<out_dir>/<id>_assembly_<n>.fasta` (a `_0` sentinel if nothing is
#' recovered), matching the GetOrganelle output contract so the rest of the
#' pipeline is unchanged.
#'
#' @param id Sample ID.
#' @param r1,r2 Trimmed paired reads.
#' @param reference Clade reference rRNA FASTA (seeds) for baiting.
#' @param labels_db Label DB FASTA for the rDNA contig filter.
#' @param out_dir Directory to write the assembly fastas into (the process cwd).
#' @param method "bait" (B2, default) or "denovo" (A, assemble all reads).
#' @param threads Threads for mapping/assembly.
#' @return (invisibly) number of contigs recovered.
#' @export
fallback_assemble_one <- function(id, r1, r2, reference, labels_db, out_dir = ".",
                                  method = c("bait", "denovo"), threads = 4) {
  method <- match.arg(method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_fa <- file.path(out_dir, paste0(id, "_assembly_1.fasta"))

  reads <- c(r1, r2)
  if (method == "bait") {
    reference <- local_ref(reference)  # resolve a URL/local seeds DB to a file
    baited <- bait_reads(r1, r2, reference, tempfile("rp_baited_"), threads = threads)
    if (is.null(baited)) {
      message("Baiting recovered no reads for ", id, "; falling back to all reads.")
    } else {
      reads <- baited
    }
  }

  sp_dir <- file.path(tempdir(), paste0("rp_spades_", id, "_", as.integer(Sys.time())))
  contigs <- run_spades(reads[1], reads[2], sp_dir, threads = threads)
  n <- extract_rdna_contigs(contigs, labels_db, id, out_fa)
  unlink(sp_dir, recursive = TRUE)

  if (n == 0L) {
    writeLines(">No assembly found", file.path(out_dir, paste0(id, "_assembly_0.fasta")))
    message("Fallback (", method, ") recovered no rDNA contigs for ", id, ".")
  } else {
    message("Fallback (", method, ") recovered ", n, " rDNA contig(s) for ", id, ".")
  }
  invisible(n)
}
