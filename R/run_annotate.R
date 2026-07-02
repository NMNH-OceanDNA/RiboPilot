# Annotation driver (WF2): separate from assembly (WF1) so the assembly can be
# validated first. Prepares each sample's stitched selected path (R), runs the
# annotation Nextflow workflow (barrnap + ITSx per sample), then reconciles the
# tool outputs into the feature model and loads them.

#' Stitch a sample's selected path against its top-hit reference
#' @return list(seq, ref_acc, path, n_gaps, n_redundant) or NULL
#' @noRd
stitch_selected <- function(con, id, seeds_cache) {
  sel <- DBI::dbGetQuery(con,
    "SELECT path, hit_id FROM path_scores WHERE ID = ? AND auto_select = 1", params = list(id))
  if (nrow(sel) == 0) {
    p <- DBI::dbGetQuery(con, "SELECT MIN(path) p FROM assemblies WHERE ID = ?", params = list(id))$p
    hit_id <- NA_character_
  } else {
    p <- sel$path[1]
    hit_id <- sel$hit_id[1]
  }
  if (is.na(p)) {
    return(NULL)
  }
  scaf <- DBI::dbGetQuery(con,
    "SELECT scaffold, sequence FROM assemblies WHERE ID = ? AND path = ? ORDER BY scaffold",
    params = list(id, p))
  if (nrow(scaf) == 0) {
    return(NULL)
  }
  seeds_db <- DBI::dbGetQuery(con,
    "SELECT ao.seeds_db FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(id))$seeds_db[1]

  seeds <- seeds_cache$get(seeds_db)
  refseq <- NULL
  if (!is.null(seeds)) {
    refseq <- if (!is.na(hit_id) && hit_id %in% names(seeds)) seeds[hit_id] else seeds[which.max(Biostrings::width(seeds))]
  }

  scaf_fa <- tempfile(fileext = ".fasta")
  writeLines(paste0(">", scaf$scaffold, "\n", scaf$sequence), scaf_fa)
  if (is.null(refseq)) {
    return(list(seq = paste(scaf$sequence, collapse = ""), ref_acc = NA_character_,
      path = p, n_gaps = NA_integer_, n_redundant = NA_integer_))
  }
  ref_fa <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(refseq, ref_fa)
  st <- tryCatch(stitch_scaffolds(scaf_fa, ref_fa), error = function(e) NULL)
  if (is.null(st)) {
    return(list(seq = paste(scaf$sequence, collapse = ""), ref_acc = sub("_.*", "", names(refseq)[1]),
      path = p, n_gaps = NA_integer_, n_redundant = NA_integer_))
  }
  list(seq = st$sequence, ref_acc = sub("_.*", "", names(refseq)[1]),
    path = p, n_gaps = st$n_gaps, n_redundant = st$n_redundant)
}

#' Run the annotation workflow (WF2) for a project
#'
#' Stitches each sample's auto-selected path, runs barrnap + ITSx via Nextflow,
#' and reconciles + loads the feature model. Run after validating assemblies
#' (WF1 / [run_assembly()]).
#'
#' @param path Project directory.
#' @param nextflow Path to the nextflow binary.
#' @param work Nextflow work dir (default `<path>/work`).
#' @param conda_env Optional conda env providing the tools (see [run_assembly()]).
#' @param cpu CPUs per annotate task.
#' @return (invisibly) the annotation samplesheet.
#' @export
run_annotate <- function(path = ".", nextflow = "nextflow", work = file.path(path, "work"),
                         conda_env = NULL, cpu = 4) {
  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  path <- normalizePath(path)
  main_nf <- app_sys("nextflow", "main_annotate.nf")
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  ids <- DBI::dbGetQuery(con, "SELECT DISTINCT ID FROM assemblies")$ID
  if (length(ids) == 0) {
    message("No assembled samples to annotate. Run assembly (WF1) first.")
    return(invisible())
  }

  # seed cache
  seeds_env <- new.env()
  seeds_cache <- list(get = function(db) {
    if (is.na(db) || !nzchar(db)) {
      return(NULL)
    }
    ref <- local_ref(db)
    if (is.na(ref)) {
      return(NULL)
    }
    if (is.null(seeds_env[[ref]])) seeds_env[[ref]] <- Biostrings::readDNAStringSet(ref)
    seeds_env[[ref]]
  })

  # ---- prep: stitch selected paths, write fastas + samplesheet ----
  in_dir <- file.path(path, "annotate_in")
  dir.create(in_dir, showWarnings = FALSE)
  ts <- as.integer(Sys.time())
  rows <- list()
  for (id in ids) {
    st <- stitch_selected(con, id, seeds_cache)
    if (is.null(st)) next
    fa <- file.path(in_dir, paste0(id, ".fasta"))
    writeLines(c(paste0(">", id), st$seq), fa)
    DBI::dbExecute(con, "DELETE FROM stitched WHERE ID = ?", params = list(id))
    DBI::dbExecute(con,
      "INSERT INTO stitched (ID, path, sequence, length, reference_accession, n_gaps, n_redundant, time_stamp)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
      params = list(id, st$path, st$seq, nchar(st$seq), st$ref_acc, st$n_gaps, st$n_redundant, ts))
    rows[[length(rows) + 1]] <- data.frame(ID = id, stitched = fa, cpus = cpu, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) {
    message("Nothing to annotate.")
    return(invisible())
  }
  ss <- do.call(rbind, rows)
  ss_path <- file.path(path, "annotate_samplesheet.tsv")
  utils::write.table(ss, ss_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # ---- run Nextflow WF2 ----
  message("Running Nextflow annotation for ", nrow(ss), " sample(s)...")
  nf_cmd <- paste(c(
    shQuote(nextflow), "run", shQuote(main_nf),
    "-work-dir", shQuote(work), "-resume",
    "--samplesheet", shQuote(ss_path), "--publishDir", "out"
  ), collapse = " ")
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) warning("Nextflow (annotate) exited with status ", status)

  # ---- load: parse tool outputs, reconcile, write annotations ----
  for (id in ss$ID) {
    adir <- file.path(path, "out", id, "annotate")
    gff <- file.path(adir, paste0(id, ".barrnap.gff"))
    pos <- file.path(adir, paste0(id, ".itsx.positions.txt"))
    bar <- if (file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
    its <- if (file.exists(pos)) parse_itsx(readLines(pos, warn = FALSE)) else data.frame()
    feat <- reconcile_features(bar, its)

    DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ?", params = list(id))
    if (nrow(feat) > 0) {
      df <- data.frame(ID = id, element = feat$element, type = feat$type,
        start = feat$start, stop = feat$stop, length = feat$length, strand = feat$strand,
        tool = feat$tool, edited = 0L, notes = NA_character_, time_stamp = ts, stringsAsFactors = FALSE)
      DBI::dbAppendTable(con, "annotations", df)
    }
    missing <- setdiff(c(RRNA_GENES, RRNA_SPACERS), feat$element)
    DBI::dbExecute(con,
      "UPDATE annotate SET n_features = ?, missing = ?, annotate_switch = 2, time_stamp = ? WHERE ID = ?",
      params = list(nrow(feat), paste(missing, collapse = ","), ts, id))
    message("Annotated ", id, ": ", nrow(feat), " features",
      if (length(missing)) paste0(" (missing: ", paste(missing, collapse = ","), ")") else "")
  }
  invisible(ss)
}
