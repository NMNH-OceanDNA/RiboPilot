#' Initialize a new RiboPilot project
#'
#' Mirrors MitoPilot's `new_project()`. Creates the project directory, copies
#' the mapping file, initializes the SQLite tracking database, and writes a
#' `.config` Nextflow config from the chosen executor template.
#'
#' @param path Path to the project directory (default = current working dir).
#' @param mapping_fn Path to a mapping CSV. Minimally an `ID` column (unique
#'   sample id), a `Taxon` column, and `R1`/`R2` columns naming the paired read
#'   inputs. Extra metadata columns are preserved.
#' @param mapping_id Column in the mapping file with unique ids (default "ID").
#' @param data_path Directory where the raw read data is located.
#' @param min_depth Minimum trimmed paired reads to proceed with assembly.
#' @param executor Executor for the Nextflow pipeline. Ready to use as-is:
#'   "conda" (default) for a local conda env, and "docker"/"local" for the
#'   RiboPilot container. Every scheduler ("slurm", "sge", "pbs", "lsf",
#'   "awsbatch", "NMNH_Hydra", "NOAA_SEDNA") is a template, not a usable
#'   executor: it needs a queue and cluster options that only [generate_config()]
#'   supplies, so build a profile there first and pass its name here. Passing a
#'   bare scheduler name is an error rather than a config Nextflow cannot parse.
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
#' @param Rproj (logical) Initialize/open an RStudio project (default TRUE).
#' @param force (logical) Recreate an existing project database (default FALSE).
#' @param ... Additional default processing parameters passed to [new_db()].
#'
#' @export
new_project <- function(
    path = ".",
    mapping_fn = NULL,
    mapping_id = "ID",
    data_path = NULL,
    min_depth = 2000000,
    executor = c("conda", "docker", "local", "awsbatch", "slurm", "pbs", "sge", "lsf", "NMNH_Hydra", "NOAA_SEDNA"),
    container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
    custom_seeds_db = NULL,
    custom_labels_db = NULL,
    ref_dir = NULL,
    config = NULL,
    profile_dir = ribopilot_config_dir(),
    ncbi_api_key = NULL,
    check_reads = TRUE,
    Rproj = TRUE,
    force = FALSE,
    ...) {
  # Create directory ----
  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  path <- normalizePath(path)

  if (length(data_path) == 1) {
    data_path <- normalizePath(data_path)
  }

  # Mapping file ----
  if (is.null(mapping_fn) || !file.exists(mapping_fn)) {
    stop("A mapping file is required to initialize a new project")
  }
  mapping_out <- file.path(path, "mapping.csv")
  if (!identical(normalizePath(mapping_fn), mapping_out)) {
    file.copy(mapping_fn, mapping_out, overwrite = TRUE)
  }

  # Executor ----
  executor <- executor[1]
  if (is.null(config) && (is.null(executor) || !nzchar(executor))) {
    stop("Invalid executor.")
  }
  if (!is.null(executor) && executor %in% CONTAINER_EXECUTORS) {
    check_container_image(container)
  }

  # RStudio project (optional) ----
  if (Rproj && !isFALSE(Sys.getenv("RSTUDIO", FALSE))) {
    if (isFALSE(requireNamespace("rstudioapi", quietly = TRUE))) {
      message("package 'rstudioapi' not available. Skipping RStudio project initialization.")
    } else {
      rstudioapi::initializeProject(path)
      on.exit(rstudioapi::openProject(path, newSession = TRUE))
    }
  }

  # Database ----
  db <- file.path(path, ".sqlite")
  if (file.exists(db) && !force) {
    message("Database already exists. Use force = TRUE to overwrite (old data will be lost).")
    return(invisible())
  }
  if (file.exists(db) && force) {
    message("Overwriting existing database")
    file.remove(db)
  }

  new_db(
    db_path = db,
    mapping_fn = mapping_out,
    mapping_id = mapping_id,
    data_path = data_path,
    check_reads = check_reads,
    seeds_db = custom_seeds_db,
    labels_db = custom_labels_db,
    ref_dir = ref_dir,
    ...
  )

  # Config file ----
  config <- config %||% resolve_config(executor, profile_dir = profile_dir)
  if (!file.exists(config)) {
    stop("Config file not found.")
  }
  cfg_lines <- readLines(config)
  # For container executors, bind-mount the directories holding any LOCAL
  # seed/label reference DBs so Nextflow's staged symlinks resolve inside the
  # container (the project + raw dirs are already mounted).
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
    ref_dirs <- unique(c(local_ref_dir(custom_seeds_db), local_ref_dir(custom_labels_db)))
    mounted <- normalizePath(c(path, if (length(data_path)) data_path), mustWork = FALSE)
    ref_dirs <- setdiff(ref_dirs, mounted)
    if (length(ref_dirs)) {
      mounts <- if (engine == "singularity") {
        paste0("-B ", ref_dirs)
      } else {
        paste0("-v ", ref_dirs, ":", ref_dirs)
      }
      extra_mounts <- paste0(" ", paste(mounts, collapse = " "))
    }
  }
  out_lines <- cfg_lines |>
    fill_config(list(
      CONTAINER_ID = container,
      RAW_DIR = data_path,
      ASMB_DIR = "NA",
      MIN_DEPTH = format(min_depth, scientific = FALSE),
      GENETIC_CODE = "1",
      NCBI_API_KEY = ncbi_api_key %||% "",
      EXTRA_MOUNTS = extra_mounts
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
  config_out <- file.path(path, ".config")
  writeLines(out_lines, config_out)
  # the config carries the NCBI API key, so keep it owner-only
  Sys.chmod(config_out, "0600")

  message("Project initialized successfully.")
  message("Please open and review the .config file to ensure all required options are specified.")
}
