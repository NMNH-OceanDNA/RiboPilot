geome_calls <- new.env()

geome_fixture_get <- function(path, query = list()) {
  geome_calls[[path]] <- (geome_calls[[path]] %||% 0L) + 1L
  f <- testthat::test_path("fixtures", "geome",
                           paste0(gsub("[^A-Za-z0-9.]", "_", path), ".json"))
  if (!file.exists(f)) stop("BCID not found in GEOME", call. = FALSE)
  jsonlite::fromJSON(f, simplifyVector = FALSE)
}

geome_reset_calls <- function() rm(list = ls(geome_calls), envir = geome_calls)

# A project built by new_project() from a mapping data frame; reads are empty files.
meta_project <- function(mapping, envir = parent.frame(), quiet = TRUE, ...) {
  dir <- withr::local_tempdir(.local_envir = envir)
  data_dir <- file.path(dir, "data")
  dir.create(data_dir)
  file.create(file.path(data_dir, c(mapping$R1, mapping$R2)))
  map_fn <- file.path(dir, "mapping.csv")
  utils::write.csv(mapping, map_fn, row.names = FALSE)
  proj <- file.path(dir, "proj")
  np <- function() suppressMessages(new_project(path = proj, mapping_fn = map_fn,
    data_path = data_dir, executor = "conda", Rproj = FALSE, custom_seeds_db = map_fn,
    custom_labels_db = map_fn, ...))
  if (quiet) suppressWarnings(np()) else np()
  list(proj = proj, data_dir = data_dir, dir = dir)
}

meta_q <- function(proj, sql) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, sql)
}

meta_reads <- function(ids) {
  data.frame(ID = ids, Taxon = "Fish a", R1 = paste0(ids, "_1.fq"), R2 = paste0(ids, "_2.fq"))
}
