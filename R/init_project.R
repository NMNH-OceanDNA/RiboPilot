#' Initialize a new RiboPilot project
#'
#' Mirrors MitoPilot's `new_project()`. Creates the project directory, copies
#' the mapping file, initializes the SQLite tracking database, and writes a
#' `.config` Nextflow config from the chosen executor template.
#'
#' @param path Path to the project directory (default = current working dir).
#' @param mapping_fn Path to a comma-separated mapping file. Minimally an `ID`
#'   column (unique sample id), a `Taxon` column, and `R1`/`R2` columns naming
#'   the paired read files relative to `data_path` (or an `Assembly` column).
#'   Column names are matched ignoring case. Extra metadata columns are
#'   preserved.
#' @param mapping_id Column in the mapping file with unique ids (default "ID").
#' @param data_path Directory where the raw read data is located. Required
#'   when any sample lists reads; not needed when every sample supplies its own
#'   assembly.
#' @param assembly_path Directory holding the assembly files named in the
#'   mapping file's optional `Assembly` column. A sample with an assembly is
#'   searched for its nuclear rDNA contigs instead of being assembled from
#'   reads, and its `R1`/`R2` may be left blank. See
#'   dev/design/20_USER_ASSEMBLIES.md.
#' @param min_depth Minimum number of reads left after trimming, counting both
#'   reads of a pair (2000000 is 1000000 pairs). A sample below it is not
#'   assembled and is marked Failed with a note.
#' @param executor Executor for the Nextflow pipeline. Ready to use as-is:
#'   "docker" (default) and "local" run the RiboPilot container on this
#'   machine, and "NMNH_Hydra" and "NOAA_SEDNA" are ready-made Singularity
#'   profiles for those clusters. "conda" uses a self-installed conda env (not
#'   supported; warns). For any other cluster, build a profile with
#'   generate_config() and pass its name here. A bare scheduler name ("slurm",
#'   "sge", "pbs", "lsf") is an error: it has no queue or cluster options yet.
#' @param container Container to use for pipeline execution.
#' @param custom_seeds_db Path/URL to a custom GetOrganelle nr seeds database.
#' @param custom_labels_db Path/URL to a custom GetOrganelle nr labels database.
#' @param ref_dir Folder for built reference databases (default: the
#'   `RIBOPILOT_REF_DIR` environment variable, else [ribopilot_ref_dir()]).
#' @param config Optional path to an existing custom Nextflow config file.
#' @param profile_dir Directory searched for saved cluster profiles.
#' @param ncbi_api_key Optional NCBI API key (raises request rate limits).
#' @param check_reads (logical) Check that the read files named in the mapping
#'   file exist under `data_path` (default TRUE). Set FALSE when the data will
#'   be staged later.
#' @param check_assemblies (logical) Check that the assembly files named in the
#'   mapping file exist under `assembly_path` (default TRUE).
#' @param mapping_geome Name of the mapping-file column holding GEOME BCIDs
#'   (optional). Stored as `GEOME_BCID`. See `vignette("Sample-Metadata")`.
#'   Passed to [new_db()].
#' @param fetch_geome Fetch GEOME metadata for samples with a BCID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_geome()] later.
#'   Passed to [new_db()].
#' @param mapping_gbif Name of the mapping-file column holding GBIF occurrence
#'   IDs (optional). Stored as `GBIF_ID`. See `vignette("Sample-Metadata")`.
#'   Passed to [new_db()].
#' @param fetch_gbif Fetch GBIF metadata for samples with a GBIF ID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_gbif()] later.
#'   Passed to [new_db()].
#' @param mapping_ncbi Name of the mapping-file column holding NCBI
#'   BioSample or SRA accessions (optional). Must not be the sample ID column.
#'   Stored as `BioSample`. See `vignette("Sample-Metadata")`.
#'   Passed to [new_db()].
#' @param fetch_ncbi Fetch NCBI metadata for samples with a BioSample
#'   value during setup (default TRUE). Set FALSE when offline and run
#'   [fetch_ncbi()] later.
#'   Passed to [new_db()].
#' @param link_sources Follow links between GEOME, GBIF, and NCBI records to
#'   fill in IDs a sample does not have yet (default FALSE). Saved as the project
#'   setting. See `vignette("Sample-Metadata")`.
#'   Passed to [new_db()].
#' @param Rproj (logical) Initialize/open an RStudio project (default TRUE; no
#'   effect outside RStudio).
#' @param force (logical) Recreate an existing project database (default
#'   FALSE). Without it, `new_project()` stops when the folder already holds a
#'   project.
#' @param ... Additional default processing parameters passed to [new_db()].
#'
#' @export
new_project <- function(
    path = ".",
    mapping_fn = NULL,
    mapping_id = "ID",
    data_path = NULL,
    assembly_path = NULL,
    min_depth = 2000000,
    executor = c("docker", "conda", "local", "slurm", "pbs", "sge", "lsf", "NMNH_Hydra", "NOAA_SEDNA"),
    container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
    custom_seeds_db = NULL,
    custom_labels_db = NULL,
    ref_dir = NULL,
    config = NULL,
    profile_dir = ribopilot_config_dir(),
    ncbi_api_key = NULL,
    check_reads = TRUE,
    check_assemblies = TRUE,
    mapping_geome = "GEOME_BCID",
    fetch_geome = TRUE,
    mapping_gbif = "GBIF_ID",
    fetch_gbif = TRUE,
    mapping_ncbi = "BioSample",
    fetch_ncbi = TRUE,
    link_sources = FALSE,
    Rproj = TRUE,
    force = FALSE,
    ...) {
  # Every check runs before anything is written ----
  path <- check_project_path(path, force)
  if (length(data_path) == 1) {
    data_path <- normalizePath(data_path, mustWork = FALSE)
  }
  if (length(assembly_path) == 1) {
    assembly_path <- normalizePath(assembly_path, mustWork = FALSE)
  }
  stop_unsafe_paths(c(data_path = data_path %||% NA, assembly_path = assembly_path %||% NA,
    custom_seeds_db = custom_seeds_db %||% NA, custom_labels_db = custom_labels_db %||% NA))

  # Mapping file ----
  mapping <- read_mapping(mapping_fn, mapping_id, list(...)[["mapping_taxon"]] %||% "Taxon") |>
    check_mapping(mapping_id, data_path, assembly_path, check_reads, check_assemblies,
      mapping_geome, mapping_gbif, mapping_ncbi)
  if (any(!is.na(mapping$R1)) && length(data_path) != 1) {
    stop("Samples list read files (R1, R2), but no data_path was given. ",
      "Pass data_path = the folder that holds them.", call. = FALSE)
  }

  # Executor ----
  executor <- executor[1]
  if (is.null(config) && (is.null(executor) || !nzchar(executor))) {
    stop("Invalid executor.")
  }
  if (!is.null(executor) && executor %in% CONTAINER_EXECUTORS) {
    check_container_image(container)
  }
  preflight <- list(
    list(mapping_geome, fetch_geome, geome_normalize_bcid, "https://api.geome-db.org/docs/geomeAPI.json", "GEOME"),
    list(mapping_gbif, fetch_gbif, gbif_normalize_id, "https://api.gbif.org/v1/enumeration/country", "GBIF"),
    list(mapping_ncbi, fetch_ncbi, ncbi_normalize_id, NULL, "NCBI"))
  for (pf in preflight) {
    if (pf[[1]] %in% names(mapping) && isTRUE(pf[[2]]) && any(!is.na(pf[[3]](mapping[[pf[[1]]]])))) {
      if (is.null(pf[[4]])) {
        ok <- tryCatch({ .ncbi_get("einfo", list()); TRUE }, error = function(e) FALSE)
        if (!ok) warning("NCBI: not reachable right now: ", NCBI_EUTILS, call. = FALSE)
      } else if (identical(url_status(pf[[4]]), "unknown")) {
        warning(pf[[5]], ": not reachable right now: ", pf[[4]], call. = FALSE)
      }
    }
  }
  if (identical(executor, "conda")) {
    warning("executor = \"conda\" is not supported: it needs every pipeline tool ",
      "installed by hand. Use executor = \"docker\" (or a Singularity profile on a cluster).",
      call. = FALSE)
  }

  # Config file ----
  config <- config %||% resolve_config(executor, profile_dir = profile_dir)
  if (!file.exists(config)) {
    stop("Config file not found.")
  }
  cfg_lines <- readLines(config)
  # For container executors, bind-mount every directory the pipeline reads that
  # is not the project folder: the reads, the supplied assemblies, and the
  # folders holding any LOCAL seed/label reference DBs (Nextflow stages those as
  # symlinks, which resolve inside the container only if their target is
  # mounted). The templates mount the project folder and nothing else, so a
  # project with no reads has no unmountable path to name: "-v NA:NA" is not a
  # mount docker or singularity accepts.
  extra_mounts <- ""
  engine <- config_container_engine(cfg_lines)
  if (!is.na(engine)) {
    local_ref_dir <- function(p) {
      if (!is.null(p) && length(p) == 1 && nzchar(p) &&
          !grepl("^https?://", p) && file.exists(p)) {
        dirname(normalizePath(p))
      } else {
        NULL
      }
    }
    dirs <- unique(c(
      if (length(data_path) == 1 && nzchar(data_path)) data_path,
      if (length(assembly_path) == 1 && nzchar(assembly_path)) assembly_path,
      local_ref_dir(custom_seeds_db), local_ref_dir(custom_labels_db)))
    dirs <- setdiff(dirs, normalizePath(path, mustWork = FALSE))
    if (length(dirs)) {
      mounts <- if (engine == "singularity") {
        paste0("-B ", dirs)
      } else {
        paste0("-v ", dirs, ":", dirs)
      }
      extra_mounts <- paste0(" ", paste(mounts, collapse = " "))
    }
  }
  out_lines <- cfg_lines |>
    fill_config(list(
      CONTAINER_ID = container,
      # 'NA' sentinel rather than an empty string for a project with no reads:
      # rawDir is only ever pasted onto a read file name, which never happens
      # for a sample that supplied its own assembly, and a visible sentinel says
      # so where an empty string would read as a missing setting.
      RAW_DIR = if (length(data_path) == 1 && nzchar(data_path)) data_path else "NA",
      ASMB_DIR = if (length(assembly_path) == 1) assembly_path else "NA",
      MIN_DEPTH = format(min_depth, scientific = FALSE),
      NCBI_API_KEY = ncbi_api_key %||% "",
      EXTRA_MOUNTS = extra_mounts,
      HOST_CPUS = ps::ps_cpu_count(),
      HOST_MEMORY = floor(ps::ps_system_memory()$total / 1024^3)
    ))
  # Scheduler templates carry tokens only generate_config() fills (queue,
  # cluster options, container engine). Writing them through produces a config
  # Nextflow cannot parse, so fail here instead of at run time.
  left <- config_tokens(out_lines)
  if (length(left)) {
    stop(glue::glue(
      "Config template '{basename(config)}' still has unfilled placeholders: ",
      "{paste(left, collapse = ', ')}. Build a profile with generate_config() ",
      "first, then pass its name as `executor`."), call. = FALSE)
  }

  # Write ----
  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  mapping_out <- file.path(path, "mapping.csv")
  if (!identical(normalizePath(mapping_fn), mapping_out)) {
    file.copy(mapping_fn, mapping_out, overwrite = TRUE)
  }
  db <- file.path(path, ".sqlite")
  if (file.exists(db)) {
    message("Overwriting existing database")
    file.remove(db)
  }
  done <- FALSE
  on.exit(if (!done) unlink(db), add = TRUE)
  # column renames were already reported by the check above
  withCallingHandlers(suppressMessages(new_db(
    db_path = db,
    mapping_fn = mapping_out,
    mapping_id = mapping_id,
    data_path = data_path,
    assembly_path = assembly_path,
    check_reads = check_reads,
    check_assemblies = check_assemblies,
    seeds_db = custom_seeds_db,
    labels_db = custom_labels_db,
    ref_dir = ref_dir,
    mapping_geome = mapping_geome,
    fetch_geome = fetch_geome,
    mapping_gbif = mapping_gbif,
    fetch_gbif = fetch_gbif,
    mapping_ncbi = mapping_ncbi,
    fetch_ncbi = fetch_ncbi,
    link_sources = link_sources,
    ...
  )), ribopilot_meta_id_warning = function(w) invokeRestart("muffleWarning"))
  config_out <- file.path(path, ".config")
  writeLines(out_lines, config_out)
  # the config carries the NCBI API key, so keep it owner-only
  Sys.chmod(config_out, "0600")
  done <- TRUE

  message("Project initialized successfully.")
  message("Please open and review the .config file to ensure all required options are specified.")
  if (Rproj && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    rstudioapi::initializeProject(path)
    rstudioapi::openProject(path, newSession = TRUE)
  }
}
