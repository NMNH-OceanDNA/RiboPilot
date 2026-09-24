#' Add samples to an existing project
#'
#' Appends new samples to a project's database, seeding every per-sample
#' table (samples, preprocess, assemble, annotate) the way [new_db()] does
#' for a new sample. Runs the same read_mapping()/check_mapping() checks as
#' [new_project()].
#'
#' @param project Path to the project directory (default = current working
#'   directory).
#' @param mapping_fn Path to a comma-separated mapping file for the new
#'   samples only. Same shape as [new_project()]'s mapping file: an `ID`
#'   column (none already in the project), a `Taxon` column, and `R1`/`R2`
#'   columns naming paired read files relative to `data_path` (or an
#'   `Assembly` column). Column names are matched ignoring case. Extra
#'   metadata columns are added to the `samples` table if not already
#'   present.
#' @param mapping_id Column in the mapping file with unique ids (default
#'   "ID").
#' @param mapping_taxon Column in the mapping file with a taxonomic
#'   identifier (default "Taxon").
#' @param data_path Directory where the new samples' raw read data is
#'   located. Required when any new sample lists reads.
#' @param assembly_path Directory holding the new samples' assembly files
#'   named in the mapping file's optional `Assembly` column.
#' @param check_reads (logical) Check that the read files named in the
#'   mapping file exist under `data_path` (default TRUE).
#' @param check_assemblies (logical) Check that the assembly files named in
#'   the mapping file exist under `assembly_path` (default TRUE).
#'
#' @export
add_samples <- function(
    project = ".",
    mapping_fn = NULL,
    mapping_id = "ID",
    mapping_taxon = "Taxon",
    data_path = NULL,
    assembly_path = NULL,
    check_reads = TRUE,
    check_assemblies = TRUE) {
  project <- normalizePath(project, mustWork = FALSE)
  db_path <- file.path(project, ".sqlite")
  if (!file.exists(db_path)) {
    stop("No RiboPilot project database found in ", project, call. = FALSE)
  }
  if (length(data_path) == 1) data_path <- normalizePath(data_path, mustWork = FALSE)
  if (length(assembly_path) == 1) assembly_path <- normalizePath(assembly_path, mustWork = FALSE)
  stop_unsafe_paths(c(data_path = data_path %||% NA, assembly_path = assembly_path %||% NA))

  mapping <- read_mapping(mapping_fn, mapping_id, mapping_taxon) |>
    check_mapping(mapping_id, data_path, assembly_path, check_reads, check_assemblies)
  if (any(!is.na(mapping$R1)) && length(data_path) != 1) {
    stop("Samples list read files (R1, R2), but no data_path was given. ",
      "Pass data_path = the folder that holds them.", call. = FALSE)
  }

  con <- open_project_db(db_path)
  on.exit(DBI::dbDisconnect(con))

  existing <- DBI::dbReadTable(con, "samples")$ID
  dup <- mapping[[mapping_id]][tolower(mapping[[mapping_id]]) %in% tolower(existing)]
  if (length(dup)) {
    stop("Already in the project (use update_sample_metadata() to change them): ",
      paste(truncate_ids(dup), collapse = ", "), call. = FALSE)
  }

  mapping <- mapping |>
    dplyr::mutate(
      ID = .data[[mapping_id]],
      Taxon = .data[[mapping_taxon]],
      taxid = NA_integer_,
      sci_name = NA_character_,
      rank = NA_character_,
      match_type = NA_character_,
      lineage = NA_character_,
      lineage_taxids = NA_character_,
      ref_id = NA_character_,
      user_group = NA_character_,
      ref_source = NA_character_,
      taxonomy_switch = 1L,
      taxonomy_lock = 0L,
      time_stamp = NA_integer_
    )

  new_cols <- setdiff(colnames(mapping), DBI::dbListFields(con, "samples"))
  for (col in new_cols) {
    glue::glue_sql("ALTER TABLE samples ADD COLUMN {col}", col = col, .con = con) |>
      DBI::dbExecute(con, statement = _)
  }

  dplyr::tbl(con, "samples") |>
    dplyr::rows_insert(mapping, in_place = TRUE, copy = TRUE, by = "ID", conflict = "ignore")

  dplyr::tbl(con, "preprocess") |>
    dplyr::rows_insert(
      mapping |>
        dplyr::select(ID, R1, R2) |>
        dplyr::mutate(
          pre_opts = "default",
          reads = NA_real_,
          trimmed_reads = NA_real_,
          mean_length = NA_real_,
          time_stamp = NA_integer_
        ),
      in_place = TRUE, copy = TRUE, by = "ID", conflict = "ignore"
    )

  dplyr::tbl(con, "assemble") |>
    dplyr::rows_insert(
      mapping |>
        dplyr::select(ID) |>
        dplyr::mutate(
          length = NA_character_,
          topology = NA_character_,
          paths = NA_integer_,
          scaffolds = NA_integer_,
          assemble_notes = NA_character_,
          find_rrna_notes = NA_character_,
          assemble_switch = 1,
          assemble_lock = 0,
          assemble_opts = "default",
          blast_opts = "default",
          time_stamp = NA_integer_
        ),
      in_place = TRUE, copy = TRUE, by = "ID", conflict = "ignore"
    )

  dplyr::tbl(con, "annotate") |>
    dplyr::rows_insert(
      data.frame(
        ID = mapping$ID, path = 1L, scaffold = 1L,
        annotate_switch = 1L, annotate_lock = 0L,
        annotate_opts = "default", curate_opts = "default", reviewed = "no"
      ),
      in_place = TRUE, copy = TRUE, by = c("ID", "path", "scaffold"), conflict = "ignore"
    )

  message("Added ", nrow(mapping), " sample(s) to ", project, ": ",
    paste(mapping$ID, collapse = ", "))
  invisible(mapping$ID)
}
