bs_recs <- function(...) {
  f <- list(...)
  data.frame(level = "BioSample", depth = 1L, ref = "SAMN1", field = names(f),
             value = unlist(f, use.names = FALSE))
}

test_that("NCBI combos pass GenBank values through and drop missing text", {
  r <- bs_recs(collection_date = "2024-06-27", geo_loc_name = "Canada: Nova Scotia",
               specimen_voucher = "ROM:12345", collected_by = "R. Martin",
               identified_by = "not collected", sex = "Female", dev_stage = "missing: control sample",
               lat_lon = "44.5 N 63.1 W", accession = "SAMN1", organism = "Zoarces americanus")
  f <- function(k) NCBI_COMBOS[[k]]$fn(r)
  expect_equal(f("collection_date"), "2024-06-27")
  expect_equal(f("geo_loc_name"), "Canada: Nova Scotia")
  expect_equal(f("specimen_voucher"), "ROM:12345")
  expect_true(is.na(f("identified_by")))
  expect_equal(f("sex"), "female")
  expect_true(is.na(f("dev_stage")))
  expect_equal(f("lat_lon"), "44.5 N 63.1 W")
  expect_equal(f("biosample"), "SAMN1")
  expect_true(is.na(f("bioproject")))
})

test_that("NCBI voucher falls back to GEOME's genbankSpecimenVoucher", {
  f <- NCBI_COMBOS$specimen_voucher$fn
  expect_equal(f(bs_recs(genbankSpecimenVoucher = "USNM:Fish:454742")), "USNM:Fish:454742")
  expect_equal(f(bs_recs(specimen_voucher = "ROM:1", genbankSpecimenVoucher = "USNM:Fish:2")), "ROM:1")
  expect_equal(f(bs_recs(specimen_voucher = "missing", genbankSpecimenVoucher = "USNM:Fish:2")),
               "USNM:Fish:2")
})

test_that("ncbi lat_lon converts decimal pairs and rejects junk", {
  f <- function(x) NCBI_COMBOS$lat_lon$fn(bs_recs(lat_lon = x))
  expect_equal(f("44.5, -63.1"), "44.5 N 63.1 W")
  expect_equal(f("-12.25 130.5"), "12.25 S 130.5 E")
  expect_equal(f("44.5 N, 63.1 W"), "44.5 N 63.1 W")
  expect_true(is.na(f("not collected")))
  expect_true(is.na(f("95, 10")))
})

test_that("bioproject combo takes the first project", {
  r <- rbind(bs_recs(accession = "SAMN1"),
             data.frame(level = "BioProject", depth = 2:3, ref = c("PRJNA2", "PRJNA3"),
                        field = "accession", value = c("PRJNA2", "PRJNA3")))
  expect_equal(NCBI_COMBOS$bioproject$fn(r), "PRJNA2")
  expect_identical(.meta_combos("NCBI"), NCBI_COMBOS)
})

test_that("SRA biosample and bioproject combos read the SRA record", {
  r <- rbind(bs_recs(accession = "SAMN1"),
             data.frame(level = "SRA", depth = 0L, ref = "SRR1", field = c("biosample", "bioproject"),
                        value = c("SAMN9", "PRJNA9")))
  expect_equal(NCBI_COMBOS$biosample_sra$fn(r), "SAMN9")
  expect_equal(NCBI_COMBOS$bioproject_sra$fn(r), "PRJNA9")
  expect_true(is.na(NCBI_COMBOS$biosample_sra$fn(bs_recs(accession = "SAMN1"))))
  expect_true(is.na(NCBI_COMBOS$bioproject_sra$fn(bs_recs(accession = "SAMN1"))))
})

test_that("NCBI concept values split geo_loc_name", {
  r <- bs_recs(geo_loc_name = "Canada: Scotian Shelf, Nova Scotia", organism = "Zoarces americanus")
  expect_equal(.ncbi_concept_value("country", r), "Canada")
  expect_equal(.ncbi_concept_value("locality", r), "Scotian Shelf, Nova Scotia")
  expect_equal(.ncbi_concept_value("taxon", r), "Zoarces americanus")
  expect_true(is.na(.ncbi_concept_value("locality", bs_recs(geo_loc_name = "Canada"))))
})

test_that("missing text in either part of geo_loc_name is blank", {
  r <- bs_recs(geo_loc_name = "USA: missing")
  expect_equal(.ncbi_concept_value("country", r), "USA")
  expect_true(is.na(.ncbi_concept_value("locality", r)))
})

test_that("more BioSample placeholders count as missing", {
  for (v in c("Not Recorded", "not available", "unspecified")) {
    expect_true(is.na(NCBI_COMBOS$sex$fn(bs_recs(sex = v))), info = v)
  }
})
