# GenBank record counts: parsing the builder script's output, and the on-disk
# cache that keeps the suggestion walk down to a handful of network calls.

test_that("parse_record_count reads the builder script's count line", {
  out <- c("GenBank query:", "  (\"Asteroidea\"[Organism]) AND (...)",
           "", "Found 4102 matching nuclear rRNA record(s).")
  expect_equal(RiboPilot:::parse_record_count(out), 4102L)
})

test_that("parse_record_count returns NA when the count line is absent", {
  expect_true(is.na(RiboPilot:::parse_record_count(c("boom", "esearch: command not found"))))
})

test_that("clade_query quotes the organism term", {
  expect_equal(RiboPilot:::clade_query("Asteroidea"), "\"Asteroidea\"[Organism]")
})

test_that("clade_query prefers the taxid, which NCBI cannot rename out from under us", {
  expect_equal(RiboPilot:::clade_query("Asteroidea", 7588L), "txid7588[Organism:exp]")
  expect_equal(RiboPilot:::clade_query("Asteroidea", NA_integer_), "\"Asteroidea\"[Organism]")
  expect_equal(RiboPilot:::clade_query("Asteroidea", NULL), "\"Asteroidea\"[Organism]")
})

test_that("cached counts are queried by taxid when one is known", {
  d <- withr::local_tempdir()
  seen <- NULL
  f <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) { seen <<- query; 12L })
  f(7588L, "Asteroidea")
  expect_equal(seen, "txid7588[Organism:exp]")
})

test_that("cached_count_fn calls the network once per taxid", {
  d <- withr::local_tempdir()
  calls <- 0L
  fake <- function(query, ...) { calls <<- calls + 1L; 4102L }
  f <- RiboPilot:::cached_count_fn(d, count_fn = fake)
  expect_equal(f(7588L, "Asteroidea"), 4102L)
  expect_equal(f(7588L, "Asteroidea"), 4102L)
  expect_equal(calls, 1L)
})

test_that("the count cache survives a new session", {
  d <- withr::local_tempdir()
  f1 <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) 99L)
  f1(7588L, "Asteroidea")
  f2 <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) stop("no network"))
  expect_equal(f2(7588L, "Asteroidea"), 99L)
})

test_that("a failed lookup is not cached", {
  d <- withr::local_tempdir()
  n <- 0L
  f <- RiboPilot:::cached_count_fn(d, count_fn = function(query, ...) {
    n <<- n + 1L
    if (n == 1L) NA_integer_ else 7L
  })
  expect_true(is.na(f(1L, "Nowhere")))
  expect_equal(f(1L, "Nowhere"), 7L)
})

test_that("read_counts_cache returns an empty frame when there is no cache", {
  d <- withr::local_tempdir()
  x <- RiboPilot:::read_counts_cache(d)
  expect_equal(nrow(x), 0L)
  expect_equal(names(x), c("taxid", "n_records"))
})

test_that("cache_count_fn answers from the cache and never calls out", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = 7898L, n_records = 11750L),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  f <- RiboPilot:::cache_count_fn(d)
  expect_equal(f(7898L, "Actinopterygii"), 11750L)
  expect_true(is.na(f(6656L, "Arthropoda")))
})
