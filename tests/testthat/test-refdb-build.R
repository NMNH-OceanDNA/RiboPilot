skip_if_not(file.exists("../../dev/refdb/build_refdb.sh"))

test_that("build_refdb.sh declares the exact Entrez query the panel is defined by", {
  s <- paste(readLines("../../dev/refdb/build_refdb.sh", warn = FALSE), collapse = "\n")
  expect_match(s, "txid33208[Organism:exp]", fixed = TRUE)
  expect_match(s, "NOT mitochondrion[Filter]", fixed = TRUE)
  expect_match(s, "NOT plastid[Filter]", fixed = TRUE)
  # the filter that silently removed every RefSeq record must never reappear
  expect_false(grepl("efilter", s, fixed = TRUE))
})

test_that("every refdb script sets the strict shell options and C locale", {
  for (f in list.files("../../dev/refdb", pattern = "[.]sh$", full.names = TRUE)) {
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_match(s, "set -Eeuo pipefail", fixed = TRUE, label = basename(f))
    expect_match(s, "LC_ALL=C", fixed = TRUE, label = basename(f))
  }
})

test_that("the docker ARG default agrees with RP_DB_VERSION", {
  skip_if_not(file.exists("../../docker/Dockerfile"))
  v <- sub('.*"(.*)".*', "\\1",
    grep("^RP_DB_VERSION", readLines("../../R/db_dist.R", warn = FALSE), value = TRUE)[1])
  arg <- grep("^ARG RP_DB_VERSION", readLines("../../docker/Dockerfile", warn = FALSE), value = TRUE)[1]
  expect_equal(sub(".*=", "", arg), v)
})
