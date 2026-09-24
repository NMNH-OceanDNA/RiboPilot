#' Update metadata for existing samples
#'
#' Updates metadata columns for samples already in a project's database
#' (adding new columns if needed). Never touches R1, R2, Assembly, or any
#' pipeline state or lock column, except that a sample whose Taxon changes
#' has its resolved taxonomy cleared, unlocked, and queued for the Taxonomy
#' step again. Use [add_samples()] to change what a sample reads or
#' assembles from. Runs the same read_mapping()/check_mapping()
#' checks as [new_project()].
#'
#' @param project Path to the project directory (default = current working
#'   directory).
#' @param mapping_fn Path to a comma-separated mapping file. Same shape as
#'   [new_project()]'s mapping file. Every `ID` must already be in the
#'   project.
#' @param mapping_id Column in the mapping file with unique ids (default
#'   "ID").
#' @param mapping_taxon Column in the mapping file with a taxonomic
#'   identifier (default "Taxon").
#' @param data_path Directory holding the read files named in the mapping
#'   file (checked, not stored).
#' @param assembly_path Directory holding the assembly files named in the
#'   mapping file (checked, not stored).
#' @param check_reads (logical) Check that the read files named in the
#'   mapping file exist under `data_path` (default TRUE).
#' @param check_assemblies (logical) Check that the assembly files named in
#'   the mapping file exist under `assembly_path` (default TRUE).
#'
#' @export
update_sample_metadata <- function(
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

  con <- open_project_db(db_path)
  on.exit(DBI::dbDisconnect(con))

  existing <- DBI::dbReadTable(con, "samples")$ID
  missing <- mapping[[mapping_id]][!tolower(mapping[[mapping_id]]) %in% tolower(existing)]
  if (length(missing)) {
    stop("Not in the project (use add_samples() to add them): ",
      paste(truncate_ids(missing), collapse = ", "), call. = FALSE)
  }

  mapping <- mapping |>
    dplyr::mutate(ID = .data[[mapping_id]], Taxon = .data[[mapping_taxon]]) |>
    dplyr::select(-dplyr::any_of(c("R1", "R2", "Assembly")))
  mapping$ID <- existing[match(tolower(mapping$ID), tolower(existing))]
  old <- DBI::dbGetQuery(con, "SELECT ID, Taxon FROM samples")
  retax <- mapping$ID[!paste(mapping$ID, mapping$Taxon) %in% paste(old$ID, old$Taxon)]

  new_cols <- setdiff(colnames(mapping), DBI::dbListFields(con, "samples"))
  for (col in new_cols) {
    glue::glue_sql("ALTER TABLE samples ADD COLUMN {col}", col = col, .con = con) |>
      DBI::dbExecute(con, statement = _)
  }

  dplyr::tbl(con, "samples") |>
    dplyr::rows_upsert(mapping, in_place = TRUE, copy = TRUE, by = "ID")
  if (length(retax)) {
    DBI::dbExecute(con, paste("UPDATE samples SET taxid = NULL, sci_name = NULL, rank = NULL,
      match_type = NULL, tax_alternatives = NULL, lineage = NULL, lineage_taxids = NULL,
      taxonomy_switch = 1, taxonomy_lock = 0 WHERE ID IN (",
      paste(rep("?", length(retax)), collapse = ","), ")"), params = as.list(retax))
    message("Taxon changed for ", paste(truncate_ids(retax), collapse = ", "),
      ": their taxonomy was cleared and unlocked. Run the Taxonomy step again for them.")
  }

  message("Updated metadata for ", nrow(mapping), " sample(s) in ", project, ": ",
    paste(mapping$ID, collapse = ", "))
  invisible(mapping$ID)
}
