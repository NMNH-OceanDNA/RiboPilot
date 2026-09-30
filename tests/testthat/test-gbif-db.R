gbif_csv <- function(dir, col = "GBIF_ID", ids = c("6186461308", ""), extra = list()) {
  m <- data.frame(ID = c("s1", "s2"), Taxon = "Fish a", R1 = c("a_1.fq", "b_1.fq"),
                  R2 = c("a_2.fq", "b_2.fq"))
  m[[col]] <- ids
  for (nm in names(extra)) m[[nm]] <- extra[[nm]]
  f <- file.path(dir, "mapping.csv")
  utils::write.csv(m, f, row.names = FALSE)
  f
}

gbif_q <- function(d, sql) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, sql)
}

gbif_mem_db <- function() {
  f <- withr::local_tempfile(fileext = ".sqlite", .local_envir = parent.frame())
  con <- DBI::dbConnect(RSQLite::SQLite(), f)
  withr::defer(DBI::dbDisconnect(con), envir = parent.frame())
  DBI::dbWriteTable(con, "samples", data.frame(ID = c("s1", "s2"), Taxon = "x"))
  con
}

test_that(".meta_ensure_tables adds GBIF_ID and GBIF values normalize", {
  con <- gbif_mem_db()
  .meta_ensure_tables(con)
  expect_true(all(c("GEOME_BCID", "GBIF_ID") %in% DBI::dbListFields(con, "samples")))
  expect_equal(.meta_store_value("GBIF", c(" https://www.gbif.org/occurrence/6186461308 ", "junk", "", NA)),
               c("6186461308", "junk", NA, NA))
})

test_that(".meta_fetch_into stores GBIF levels next to GEOME and a refetch failure keeps data", {
  local_mocked_bindings(.gbif_get = gbif_fixture_get, .geome_get = geome_fixture_get)
  con <- gbif_mem_db()
  .meta_fetch_into(con, "GEOME", "s1", "ark:/21547/CYB2REEDY")
  .meta_fetch_into(con, "GBIF", "s1", "6186461308")
  lv <- DBI::dbGetQuery(con, "SELECT DISTINCT level FROM meta_records WHERE source = 'GBIF'")$level
  expect_setequal(lv, c("Occurrence", "Dataset", "Organization"))
  expect_setequal(DBI::dbGetQuery(con, "SELECT source FROM meta_status")$source, c("GEOME", "GBIF"))
  n <- DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM meta_records WHERE source = 'GBIF'")$n
  local_mocked_bindings(.gbif_get = function(path) stop("could not reach GBIF (timeout)", call. = FALSE))
  expect_warning(.meta_fetch_into(con, "GBIF", "s1", "6186461308"), "GBIF fetch failed for s1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM meta_records WHERE source = 'GBIF'")$n, n)
  st <- DBI::dbGetQuery(con, "SELECT status, message FROM meta_status WHERE source = 'GBIF'")
  expect_equal(st$status, "failed")
  expect_match(st$message, "could not reach GBIF")
  expect_gt(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM meta_records WHERE source = 'GEOME'")$n, 0)
})

test_that("a non-numeric GBIF ID is marked failed without calling GBIF", {
  local_mocked_bindings(.gbif_get = function(...) stop("should not be called"))
  con <- gbif_mem_db()
  res <- suppressWarnings(.meta_fetch_into(con, "GBIF", "s1", "abc"))
  expect_equal(res$status, "failed")
  expect_match(res$message, "not a GBIF occurrence ID")
})

test_that("fetch_gbif sets, fetches, and clears IDs without touching GEOME", {
  local_mocked_bindings(.gbif_get = gbif_fixture_get)
  d <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  DBI::dbWriteTable(con, "samples", data.frame(ID = c("s1", "s2"), Taxon = "x", GEOME_BCID = "ark:/1/A"))
  .meta_ensure_tables(con)
  DBI::dbExecute(con, "INSERT INTO meta_status VALUES ('s1', 'GEOME', 'ark:/1/A', 'ok', NULL, 1)")
  DBI::dbDisconnect(con)
  fetch_gbif(d, ids = c("s1", "s2"),
             gbifs = c("https://www.gbif.org/occurrence/6186461308", "2336663130"))
  expect_equal(gbif_q(d, "SELECT GBIF_ID FROM samples ORDER BY ID")$GBIF_ID, c("6186461308", "2336663130"))
  fetch_gbif(d, ids = "s1", gbifs = "")
  expect_true(is.na(gbif_q(d, "SELECT GBIF_ID FROM samples WHERE ID = 's1'")$GBIF_ID))
  expect_equal(gbif_q(d, "SELECT ID FROM meta_status WHERE source = 'GBIF'")$ID, "s2")
  expect_equal(gbif_q(d, "SELECT ID FROM meta_status WHERE source = 'GEOME'")$ID, "s1")
  expect_error(fetch_gbif(d, ids = c("s1", "s2"), gbifs = "1"), "ids and gbifs must be the same length")
  expect_message(res <- fetch_gbif(d, ids = "s1"), "No samples with a GBIF ID to fetch")
})
