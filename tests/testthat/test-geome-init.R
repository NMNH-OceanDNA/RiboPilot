geome_mapping <- function(dir, bcid_col = "GEOME_BCID", bcids = c("ark:/21547/CYB2REEDY", "")) {
  m <- data.frame(ID = c("s1", "s2"), Taxon = "Fish a", R1 = c("a_1.fq", "b_1.fq"),
                  R2 = c("a_2.fq", "b_2.fq"))
  m[[bcid_col]] <- bcids
  f <- file.path(dir, "mapping.csv")
  utils::write.csv(m, f, row.names = FALSE)
  f
}

geome_new_db <- function(d, f, ...) {
  new_db(db_path = file.path(d, ".sqlite"), mapping_fn = f, seeds_db = f, labels_db = f,
    check_reads = FALSE, ...)
}

test_that("check_mapping warns on a non-ARK BCID and errors on a missing custom column", {
  m <- data.frame(ID = "s1", Taxon = "x", R1 = "a", R2 = "b", Bcid = "junk")
  expect_warning(RiboPilot:::check_mapping(m, check_reads = FALSE, check_assemblies = FALSE,
    mapping_geome = "Bcid"), "GEOME")
  expect_error(RiboPilot:::check_mapping(m, check_reads = FALSE, check_assemblies = FALSE,
    mapping_geome = "Nope"), "Nope")
  expect_error(RiboPilot:::check_mapping(m, check_reads = FALSE, check_assemblies = FALSE,
    mapping_geome = "ID"), "can't also be the GEOME")
})

test_that("new_db stores a renamed BCID column and fetches it", {
  local_mocked_bindings(.geome_get = geome_fixture_get)
  d <- withr::local_tempdir()
  geome_new_db(d, geome_mapping(d, "Bcid"), mapping_geome = "Bcid")
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  s <- DBI::dbGetQuery(con, "SELECT * FROM samples ORDER BY ID")
  expect_false("Bcid" %in% names(s))
  expect_equal(s$GEOME_BCID, c("ark:/21547/CYB2REEDY", NA))
  expect_equal(DBI::dbGetQuery(con, "SELECT ID, status FROM meta_status")$ID, "s1")
})

test_that("new_db fetch_geome = FALSE stores BCIDs without calling GEOME", {
  local_mocked_bindings(.geome_get = function(...) stop("should not be called"))
  d <- withr::local_tempdir()
  geome_new_db(d, geome_mapping(d), fetch_geome = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM meta_status")$n, 0L)
  expect_true("meta_export_fields" %in% DBI::dbListTables(con))
})

test_that("new_db without any BCID column still creates GEOME tables", {
  d <- withr::local_tempdir()
  f <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "s1", Taxon = "x", R1 = "a_1.fq", R2 = "a_2.fq"), f, row.names = FALSE)
  geome_new_db(d, f)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_true("GEOME_BCID" %in% DBI::dbListFields(con, "samples"))
})

test_that("new_project passes the metadata arguments through and checks reachability once", {
  local_mocked_bindings(.geome_get = geome_fixture_get)
  seen <- character()
  local_mocked_bindings(url_status = function(u) { seen <<- c(seen, u); "ok" })
  m <- meta_reads(c("s1", "s2"))
  m$Bcid <- c("ark:/21547/CYB2REEDY", "")
  p <- meta_project(m, mapping_geome = "Bcid")
  expect_length(seen, 1L)
  expect_match(seen, "geome-db")
  expect_equal(meta_q(p$proj, "SELECT GEOME_BCID FROM samples ORDER BY ID")$GEOME_BCID,
    c("ark:/21547/CYB2REEDY", NA))
  expect_equal(meta_q(p$proj, "SELECT ID FROM meta_status")$ID, "s1")
})

test_that("the reachability check warns, never stops, and is skipped when fetch is off", {
  local_mocked_bindings(url_status = function(u) "unknown")
  m <- meta_reads(c("s1", "s2"))
  m$GEOME_BCID <- c("ark:/21547/CYB2REEDY", "")
  d <- withr::local_tempdir()
  dd <- file.path(d, "data"); dir.create(dd)
  file.create(file.path(dd, c(m$R1, m$R2)))
  fn <- file.path(d, "mapping.csv")
  utils::write.csv(m, fn, row.names = FALSE)
  run <- function(...) new_project(path = file.path(d, "proj"), mapping_fn = fn, data_path = dd,
    executor = "conda", Rproj = FALSE, custom_seeds_db = fn, custom_labels_db = fn, force = TRUE, ...)
  w <- testthat::capture_warnings(suppressMessages(run(fetch_geome = FALSE)))
  expect_false(any(grepl("not reachable", w)))
  local_mocked_bindings(.geome_get = function(...) stop("offline"))
  w <- testthat::capture_warnings(suppressMessages(run()))
  expect_true(any(grepl("GEOME: not reachable right now", w)))
})
