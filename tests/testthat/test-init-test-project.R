test_that("bundled (Local) test samples are copied, not downloaded", {
  d <- withr::local_tempdir()
  src <- file.path(d, "src"); dir.create(src)
  writeLines("x", file.path(src, "A_R1.fastq.gz")); writeLines("y", file.path(src, "A_R2.fastq.gz"))
  map <- data.frame(ID = "A", Taxon = "t", R1 = "A_R1.fastq.gz", R2 = "A_R2.fastq.gz", Local = TRUE)
  out <- file.path(d, "data"); dir.create(out)
  RiboPilot:::.stage_test_reads(map, out, local_dir = src)
  expect_true(file.exists(file.path(out, "A_R1.fastq.gz")))
  expect_true(file.exists(file.path(out, "A_R2.fastq.gz")))
})

test_that("the bundled chimera samples are present and mapped", {
  map <- read.csv(system.file("test_data", "mapping_test.csv", package = "RiboPilot"))
  chim <- map[map$ID %in% c("CHIM00", "CHIM01"), ]
  expect_equal(nrow(chim), 2)
  expect_true(all(chim$Local))
  for (f in c(chim$R1, chim$R2)) expect_true(file.exists(system.file("test_data", f, package = "RiboPilot")))
})

test_that("n limits downloaded samples but bundled ones are always kept", {
  map <- read.csv(system.file("test_data", "mapping_test.csv", package = "RiboPilot"))
  kept <- RiboPilot:::.select_test_samples(map, 1)
  expect_equal(kept$ID, c(map$ID[!(map$Local %in% TRUE)][1], "CHIM00", "CHIM01"))
  expect_equal(nrow(RiboPilot:::.select_test_samples(map, Inf)), nrow(map))
})

test_that("bundled reads are unambiguously Phred+33", {
  for (id in c("CHIM00", "CHIM01")) {
    f <- system.file("test_data", paste0(id, "_R1.fastq.gz"), package = "RiboPilot")
    q <- readLines(gzfile(f), n = 4)[4]
    expect_lt(min(utf8ToInt(q)), 64L)
  }
})

test_that("the small path fetches filtered reads from the release and skips present ones", {
  d <- withr::local_tempdir()
  out <- file.path(d, "data"); dir.create(out)
  fetched <- character(0)
  fake_fetch <- function(base_url, tag, asset, dir) {
    fetched <<- c(fetched, paste(tag, asset))
    f <- file.path(dir, asset); con <- gzfile(f, "w"); writeLines(c("@r", "A", "+", "?"), con); close(con)
    f
  }
  testthat::local_mocked_bindings(fetch_verified = fake_fetch, .package = "RiboPilot")
  map <- data.frame(ID = "SRR1", Taxon = "t", R1 = "SRR1_R1.fastq.gz", R2 = "SRR1_R2.fastq.gz", Local = FALSE)
  RiboPilot:::.stage_test_reads(map, out, full_size = FALSE)
  expect_equal(fetched, paste(RiboPilot:::RP_TEST_READS_TAG, c("SRR1_R1.fastq.gz", "SRR1_R2.fastq.gz")))
  expect_true(RiboPilot:::.gz_ok(file.path(out, "SRR1_R1.fastq.gz")))
  RiboPilot:::.stage_test_reads(map, out, full_size = FALSE)
  expect_length(fetched, 2)
})

test_that("a failed small-path download names the full_size escape hatch", {
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(fetch_verified = function(...) stop("HTTP 404"), .package = "RiboPilot")
  expect_error(RiboPilot:::.fetch_small_reads("SRR1", "a.gz", "b.gz", d), "full_size = TRUE")
})
