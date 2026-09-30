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

test_that("fetch_biosample sets IDs and stores records", {
  local_mocked_bindings(.ncbi_get = ncbi_fixture_get)
  dir <- ncbi_proj()
  res <- fetch_biosample(dir, ids = "s2", biosamples = "samn63236902")
  expect_equal(res$status, "ok")
  expect_equal(q(dir, "SELECT BioSample FROM samples WHERE ID = 's2'")$BioSample, "SAMN63236902")
  r <- q(dir, "SELECT DISTINCT level FROM meta_records WHERE ID = 's2' AND source = 'NCBI'")
  expect_setequal(r$level, c("BioSample", "BioProject"))
})

test_that("a failed refresh keeps earlier NCBI data", {
  local_mocked_bindings(.ncbi_get = ncbi_fixture_get)
  dir <- ncbi_proj()
  fetch_biosample(dir, ids = "s2", biosamples = "SAMN63236902")
  n <- q(dir, "SELECT COUNT(*) n FROM meta_records WHERE ID = 's2'")$n
  local_mocked_bindings(.ncbi_get = function(...) stop("could not reach NCBI (offline)", call. = FALSE))
  expect_warning(fetch_biosample(dir, ids = "s2"), "offline")
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
