# Assembly driver: build a samplesheet from the project database, run the
# Nextflow assembly pipeline (preprocess + GetOrganelle), and load the published
# outputs back into the database (assemblies + assemble summary).
#
# DB I/O is driver-side (R), not in Nextflow, so no nf-sqldb plugin is needed.

#' Resolve a conda env (name or path) to its `bin/` directory
#' @noRd
resolve_conda_bin <- function(env) {
  # already a bin dir
  if (dir.exists(env) && basename(env) == "bin") {
    return(normalizePath(env))
  }
  # a path to an env prefix
  if (dir.exists(file.path(env, "bin"))) {
    return(normalizePath(file.path(env, "bin")))
  }
  # a bare env name: search common conda roots
  roots <- c(
    file.path(dirname(dirname(Sys.getenv("CONDA_PREFIX", ""))), "envs"),
    file.path(path.expand("~"), c("miniconda3", "miniforge3", "mambaforge", "anaconda3"), "envs")
  )
  for (r in roots) {
    cand <- file.path(r, env, "bin")
    if (dir.exists(cand)) {
      return(normalizePath(cand))
    }
  }
  stop("Could not locate conda env '", env, "'. Pass a full path to the env or its bin/.", call. = FALSE)
}

#' Read a filled `params.X = '...'` value from the project .config
#' @noRd
config_param <- function(config_lines, key) {
  hit <- grep(paste0("^\\s*", key, "\\s*="), config_lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  val <- sub(paste0("^\\s*", key, "\\s*=\\s*"), "", hit[1])
  val <- sub("//.*$", "", val)
  trimws(gsub("['\"]", "", val))
}

#' Run the assembly pipeline for a project
#'
#' Queries the database for samples ready to assemble, writes a Nextflow
#' samplesheet, runs the pipeline, and loads the results.
#'
#' @param path Project directory (contains `.sqlite` and `.config`).
#' @param nextflow Path to the `nextflow` binary (default: on PATH).
#' @param work Nextflow work directory (default `<path>/work`).
#' @param load (logical) Load published outputs into the database afterwards.
#' @param score (logical) After loading, score candidate paths and write
#'   `path_scores` (see [score_project()]).
#' @param conda_env Optional conda environment providing the pipeline tools
#'   (nextflow, fastp, GetOrganelle, ...). Give the env name (looked up under
#'   `~/miniconda3|miniforge3|mambaforge/envs`, or `$CONDA_PREFIX`'s sibling) or a
#'   full path to the env or its `bin/`. Its `bin/` is prepended to `PATH` for
#'   the run, so you don't have to `conda activate` before starting RStudio.
#' @return (invisibly) the samplesheet data.frame.
#' @export
run_assembly <- function(path = ".", nextflow = "nextflow", work = file.path(path, "work"),
                         load = TRUE, score = TRUE, conda_env = NULL) {
  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  path <- normalizePath(path)
  db <- file.path(path, ".sqlite")
  cfg <- file.path(path, ".config")
  main_nf <- app_sys("nextflow", "main.nf")
  stopifnot(file.exists(db), file.exists(cfg), file.exists(main_nf))

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  ss <- DBI::dbGetQuery(con, "
    SELECT a.ID, p.R1, p.R2, ao.cpus, ao.memory, po.fastp AS fastp,
           ao.getOrganelle AS getOrganelle, ao.seeds_db, ao.labels_db
    FROM assemble a
    JOIN preprocess p     ON a.ID = p.ID
    JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
    JOIN pre_opts po      ON p.pre_opts = po.pre_opts
    WHERE a.assemble_switch IN (1, 4) AND a.assemble_lock = 0")
  if (nrow(ss) == 0) {
    message("No samples ready to assemble.")
    return(invisible(ss))
  }

  # expand any '~' in local DB paths (Nextflow/GetOrganelle treat '~' literally)
  expand_local <- function(x) ifelse(grepl("^https?://", x), x, path.expand(x))
  ss$seeds_db <- expand_local(ss$seeds_db)
  ss$labels_db <- expand_local(ss$labels_db)

  ss_path <- file.path(path, "samplesheet.tsv")
  utils::write.table(ss, ss_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  cfg_lines <- readLines(cfg)
  raw_dir <- config_param(cfg_lines, "rawDir")
  min_depth <- config_param(cfg_lines, "minDepth")

  message("Running Nextflow assembly for ", nrow(ss), " sample(s)...")
  # Launch Nextflow FROM the project dir so launchDir (hence publishDir out/) is
  # the project, not the caller's cwd.
  nf_cmd <- paste(c(
    shQuote(nextflow), "run", shQuote(main_nf), "-c", shQuote(cfg),
    "-work-dir", shQuote(work), "-resume",
    "--samplesheet", shQuote(ss_path),
    "--rawDir", shQuote(raw_dir),
    "--minDepth", min_depth
  ), collapse = " ")
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) {
    warning("Nextflow exited with status ", status)
  }

  if (isTRUE(load)) {
    out_dir <- file.path(path, "out")
    for (id in ss$ID) load_assemblies(con, out_dir, id)
    message("Loaded assemblies into the database.")
    if (isTRUE(score)) {
      score_project(path = path, out_dir = out_dir, con = con)
    }
  }
  message("Assembly (WF1) complete. Review/validate in the app, then run ",
    "annotation with run_annotate().")
  invisible(ss)
}

#' Load one sample's published assembly fastas into the database
#'
#' Parses `out/<ID>/assemble/<ID>_assembly_*.fasta` (headers
#' `>ID.path.scaffold topology`), writes per-scaffold rows into `assemblies`,
#' and updates the per-sample `assemble` summary.
#'
#' @param con Open DBI connection.
#' @param out_dir The pipeline `out/` directory.
#' @param ID Sample ID.
#' @return (invisibly) the number of scaffolds loaded.
#' @export
load_assemblies <- function(con, out_dir, ID) {
  fas <- Sys.glob(file.path(out_dir, ID, "assemble", paste0(ID, "_assembly_*.fasta")))
  fas <- fas[!grepl("_assembly_0\\.fasta$", fas)] # skip the "no assembly" sentinel
  ts <- as.integer(Sys.time())

  if (length(fas) == 0) {
    DBI::dbExecute(con,
      "UPDATE assemble SET paths=NULL, scaffolds=NULL, length=NULL, topology=NULL,
       assemble_switch=3, assemble_notes='failed assembly', time_stamp=? WHERE ID=?",
      params = list(ts, ID))
    return(invisible(0L))
  }

  rows <- list()
  for (fa in fas) {
    seqs <- Biostrings::readDNAStringSet(fa)
    for (i in seq_along(seqs)) {
      hdr <- names(seqs)[i]
      parts <- strsplit(hdr, "\\s+")[[1]]
      idpath <- strsplit(parts[1], "\\.")[[1]] # ID, path, scaffold
      topo <- if (length(parts) > 1) parts[2] else NA_character_
      s <- as.character(seqs[[i]])
      rows[[length(rows) + 1]] <- data.frame(
        ID = ID, path = as.integer(idpath[2]), scaffold = as.integer(idpath[3]),
        topology = topo, length = nchar(s), length_raw = nchar(s),
        sequence = s, ignore = 0L, edited = 0L, time_stamp = ts,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, rows)

  # replace this sample's assemblies
  DBI::dbExecute(con, "DELETE FROM assemblies WHERE ID = ?", params = list(ID))
  DBI::dbAppendTable(con, "assemblies", df)

  # per-sample summary
  n_paths <- length(unique(df$path))
  n_scaffolds <- max(table(df$path))
  length_str <- paste(sort(unique(df$length), decreasing = TRUE), collapse = ";")
  topo_str <- paste(sort(unique(stats::na.omit(df$topology))), collapse = ";")
  notes <- if (n_paths > 1) "Unable to resolve single assembly from reads"
  else if (n_scaffolds > 1) "Output contains disconnected contigs" else ""
  DBI::dbExecute(con,
    "UPDATE assemble SET paths=?, scaffolds=?, length=?, topology=?,
     assemble_switch=2, assemble_notes=?, time_stamp=? WHERE ID=?",
    params = list(n_paths, n_scaffolds, length_str, topo_str, notes, ts, ID))
  invisible(nrow(df))
}
