ncbi_proj <- function(ids = c("SRR21844202", "s2"), envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbWriteTable(con, "samples", data.frame(ID = ids, Taxon = "Fundulus majalis"))
  DBI::dbDisconnect(con)
  dir
}
q <- function(dir, sql) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, sql)
}

test_that(".meta_ensure_tables adds a BioSample column", {
  dir <- ncbi_proj()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  .meta_ensure_tables(con)
  expect_true("BioSample" %in% DBI::dbListFields(con, "samples"))
})

test_that("fetch_ncbi sets IDs and stores records", {
  local_mocked_bindings(.ncbi_get = ncbi_fixture_get)
  dir <- ncbi_proj()
  res <- fetch_ncbi(dir, ids = "s2", ncbi_ids = "samn63236902")
  expect_equal(res$status, "ok")
  expect_equal(q(dir, "SELECT BioSample FROM samples WHERE ID = 's2'")$BioSample, "SAMN63236902")
  r <- q(dir, "SELECT DISTINCT level FROM meta_records WHERE ID = 's2' AND source = 'NCBI'")
  expect_setequal(r$level, c("BioSample", "BioProject"))
})

test_that("a failed refresh keeps earlier NCBI data", {
  local_mocked_bindings(.ncbi_get = ncbi_fixture_get)
  dir <- ncbi_proj()
  fetch_ncbi(dir, ids = "s2", ncbi_ids = "SAMN63236902")
  n <- q(dir, "SELECT COUNT(*) n FROM meta_records WHERE ID = 's2'")$n
  local_mocked_bindings(.ncbi_get = function(...) stop("could not reach NCBI (offline)", call. = FALSE))
  expect_warning(fetch_ncbi(dir, ids = "s2"), "offline")
  expect_equal(q(dir, "SELECT COUNT(*) n FROM meta_records WHERE ID = 's2'")$n, n)
})

test_that(".meta_take_cols copies source columns and refuses the ID column", {
  m <- data.frame(ID = c("SRR1", "SRR2"), Taxon = "x")
  expect_error(.meta_take_cols(m, c(NCBI = "ID")), "can't also be the NCBI ID column")
  m1 <- data.frame(SampleID = c("SRR1", "SRR2"), Taxon = "x")
  m1$ID <- m1$SampleID
  expect_error(.meta_take_cols(m1, c(NCBI = "SampleID"), keep = "SampleID"), "SampleID")
  m2 <- data.frame(run = c("SRR1", "SRR2"), ID = c("SRR1", "SRR2"), Taxon = "x")
  out2 <- .meta_take_cols(m2, c(NCBI = "run"))
  expect_equal(out2$BioSample, c("SRR1", "SRR2"))
  m4 <- data.frame(ID = c("s1", "s2"), BioSample = c("SRR1", "junk"), Taxon = "x")
  expect_equal(.meta_take_cols(m4, c(NCBI = "BioSample"))$BioSample, c("SRR1", "junk"))
})

ncbi_mapping <- function(dir) {
  m <- data.frame(ID = c("SRR21844202", "s2"), Taxon = "Fundulus majalis",
                  R1 = c("a_1.fq", "b_1.fq"), R2 = c("a_2.fq", "b_2.fq"),
                  BioSample = c("SRR21844202", NA))
  f <- file.path(dir, "mapping.csv")
  utils::write.csv(m, f, row.names = FALSE)
  f
}

test_that(".meta_take_cols lets two sources share one column", {
  m2 <- data.frame(ID = "s1", Taxon = "x", Acc = "123")
  out2 <- .meta_take_cols(m2, c(GBIF = "Acc", NCBI = "Acc"))
  expect_equal(out2$GBIF_ID, "123")
  expect_equal(out2$BioSample, "123")
  expect_false("Acc" %in% names(out2))
})

ncbi_mapping <- function(dir) {
  m <- data.frame(ID = c("SRR21844202", "s2"), Taxon = "Fundulus majalis",
                  R1 = c("a_1.fq", "b_1.fq"), R2 = c("a_2.fq", "b_2.fq"),
                  BioSample = c("SRR21844202", NA))
  f <- file.path(dir, "mapping.csv")
  utils::write.csv(m, f, row.names = FALSE)
  f
}

ncbi_new_db <- function(d, f, ...) {
  new_db(db_path = file.path(d, ".sqlite"), mapping_fn = f, seeds_db = f, labels_db = f,
    check_reads = FALSE, ...)
}

test_that("check_mapping validates BioSample values and refuses the ID column", {
  m <- data.frame(ID = "s1", Taxon = "x", R1 = "a", R2 = "b", Acc = "junk")
  chk <- function(...) RiboPilot:::check_mapping(m, check_reads = FALSE, check_assemblies = FALSE, ...)
  expect_warning(chk(mapping_ncbi = "Acc"), "BioSample")
  expect_error(chk(mapping_ncbi = "Nope"), "Nope")
  expect_error(chk(mapping_ncbi = "ID"), "can't also be the NCBI ID column")
  m$Acc <- "SAMN00000001"
  expect_no_warning(chk(mapping_ncbi = "Acc"))
  m$BioSample <- "SAMN00000002"
  expect_error(chk(mapping_ncbi = "Acc"), "BioSample")
})

test_that("new_db stores BioSample from its column and fetches it", {
  local_mocked_bindings(.ncbi_get = ncbi_fixture_get)
  d <- withr::local_tempdir()
  ncbi_new_db(d, ncbi_mapping(d))
  s <- q(d, "SELECT ID, BioSample FROM samples ORDER BY ID")
  expect_equal(s$ID, c("SRR21844202", "s2"))
  expect_equal(s$BioSample, c("SRR21844202", NA))
  expect_equal(q(d, "SELECT ID, status FROM meta_status WHERE source = 'NCBI'")$status, "ok")
})

test_that("new_db fetch_ncbi = FALSE stores IDs without calling NCBI", {
  local_mocked_bindings(.ncbi_get = function(...) stop("should not be called"))
  d <- withr::local_tempdir()
  ncbi_new_db(d, ncbi_mapping(d), fetch_ncbi = FALSE)
  expect_equal(q(d, "SELECT COUNT(*) n FROM meta_status")$n, 0L)
})

test_that("new_project warns when NCBI is unreachable and the check uses a real request", {
  seen <- NULL
  local_mocked_bindings(.ncbi_get = function(endpoint, query) stop("could not reach NCBI (timeout)", call. = FALSE))
  m <- meta_reads(c("s1", "s2"))
  m$BioSample <- c("SAMN00000001", NA)
  w <- testthat::capture_warnings(meta_project(m, quiet = FALSE))
  expect_true(any(grepl("NCBI: not reachable right now", w, fixed = TRUE)))
})
