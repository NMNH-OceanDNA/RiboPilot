#' Initialize a test project
#'
#' Sets up a test project using the same fish samples as MitoPilot and fetches
#' the associated data. `n` limits the number of samples for faster setup.
#'
#' NOTE: unlike MitoPilot, RiboPilot's target locus (nuclear rDNA) is NOT present
#' in MitoPilot's pre-filtered (mito-mapped) read files. Until RiboPilot ships
#' its own rDNA-filtered subsets, use `full_size = TRUE` to download the full
#' reads from ENA.
#'
#' @param path Path to the test project directory (default = current working
#'   dir). Created if missing.
#' @param n How many samples to include (default = Inf, all).
#' @param full_size (logical) Download full raw reads from ENA (default = FALSE).
#'   Required for meaningful rDNA assembly until filtered subsets are bundled.
#' @param executor Executor for the Nextflow pipeline (default "local").
#' @param container Container for the pipeline.
#' @param Rproj (logical) Initialize/open an RStudio project (default TRUE).
#' @param force (logical) Recreate an existing project database (default FALSE).
#' @param ... Additional arguments passed to [new_project()] / [new_db()].
#'
#' @export
new_test_project <- function(
    path = ".",
    n = Inf,
    full_size = FALSE,
    executor = "local",
    container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
    Rproj = TRUE,
    force = FALSE,
    ...) {
  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  path <- normalizePath(path)

  mapping <- app_sys(file.path("test_data", "mapping_test.csv")) |>
    utils::read.csv() |>
    dplyr::slice_head(n = n)
  readr::write_csv(mapping, file.path(path, "mapping.csv"), quote = "none", na = "")

  # Data ----
  data_dir <- file.path(path, "data")
  if (!dir.exists(data_dir)) dir.create(data_dir)

  if (full_size) {
    message("Fetching full test data from ENA. This will take a while...")
    purrr::pwalk(mapping, function(...) {
      cur <- list(...)
      message(glue::glue("{cur$ID} - {cur$Taxon}"))
      acc <- cur$ID
      pre <- stringr::str_sub(acc, 1, 6)
      suf <- stringr::str_extract(acc, "..$") |> stringr::str_pad(3, "left", "0")
      for (r in c("1", "2")) {
        fn <- file.path(data_dir, glue::glue("{acc}_R{r}.fastq.gz"))
        ok <- function() {
          glue::glue("-t {fn} >/dev/null 2>&1 && echo complete || echo incomplete") |>
            system2("gzip", args = _, stdout = TRUE)
        }
        status <- if (file.exists(fn)) ok() else "incomplete"
        while (status == "incomplete") {
          glue::glue(
            "curl http://ftp.sra.ebi.ac.uk/vol1/fastq/{pre}/{suf}/{acc}/{acc}_{r}.fastq.gz",
            "--silent -o {fn}",
            .sep = " "
          ) |> system()
          status <- ok()
        }
      }
    })
  } else {
    have_files <- purrr::every(
      c(mapping$R1, mapping$R2),
      ~ nzchar(app_sys(file.path("test_data", .x)))
    )
    if (!have_files) {
      stop(
        "Pre-filtered rDNA test reads are not bundled yet.\n",
        "Re-run with full_size = TRUE to download full reads from ENA.",
        call. = FALSE
      )
    }
    purrr::pwalk(mapping, function(...) {
      cur <- list(...)
      message(glue::glue("{cur$ID} - {cur$Taxon}"))
      file.copy(app_sys(file.path("test_data", cur$R1)), data_dir)
      file.copy(app_sys(file.path("test_data", cur$R2)), data_dir)
    })
  }

  # Initialize project ----
  new_project(
    path = path,
    mapping_fn = file.path(path, "mapping.csv"),
    mapping_id = "ID",
    data_path = file.path(path, "data", ""),
    container = container,
    min_depth = 500,
    executor = executor,
    Rproj = Rproj,
    force = force,
    ...
  )
}
