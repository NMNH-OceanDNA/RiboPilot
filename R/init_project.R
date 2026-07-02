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
#' @param executor Executor for the Nextflow pipeline: a built-in template
#'   ("local" (default), "awsbatch", "slurm", "sge", "pbs", "lsf", "NMNH_Hydra",
#'   "NOAA_SEDNA") or a saved profile from [generate_config()].
#' @param container Container to use for pipeline execution.
#' @param custom_seeds_db Path/URL to a custom GetOrganelle nr seeds database.
#' @param custom_labels_db Path/URL to a custom GetOrganelle nr labels database.
#' @param config Optional path to an existing custom Nextflow config file.
#' @param profile_dir Directory searched for saved cluster profiles.
#' @param ncbi_api_key Optional NCBI API key (raises request rate limits).
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
    executor = c("local", "awsbatch", "slurm", "sge", "pbs", "lsf", "NMNH_Hydra", "NOAA_SEDNA"),
    container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
    custom_seeds_db = NULL,
    custom_labels_db = NULL,
    config = NULL,
    profile_dir = ribopilot_config_dir(),
    ncbi_api_key = NULL,
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
    seeds_db = custom_seeds_db,
    labels_db = custom_labels_db,
    ...
  )

  # Config file ----
  config <- config %||% resolve_config(executor, profile_dir = profile_dir)
  if (!file.exists(config)) {
    stop("Config file not found.")
  }
  readLines(config) |>
    fill_config(list(
      CONTAINER_ID = container,
      RAW_DIR = data_path,
      ASMB_DIR = "NA",
      MIN_DEPTH = format(min_depth, scientific = FALSE),
      GENETIC_CODE = "1",
      NCBI_API_KEY = ncbi_api_key %||% ""
    )) |>
    writeLines(file.path(path, ".config"))

  message("Project initialized successfully.")
  message("Please open and review the .config file to ensure all required options are specified.")
}
