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
