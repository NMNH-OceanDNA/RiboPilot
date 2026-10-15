feats_df <- function(...) {
  rows <- list(...)
  data.frame(path = vapply(rows, `[[`, 1, 1), scaffold = vapply(rows, `[[`, 1, 2),
    element = vapply(rows, `[[`, "", 3), stringsAsFactors = FALSE)
}

test_that("contigs sharing an anchor form one set; stray contigs are dropped", {
  f <- feats_df(list(1, 1, "18S"), list(1, 1, "28S"), list(2, 1, "18S"), list(2, 2, "28S"),
    list(3, 1, "ITS1"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(s[[1]]$units, data.frame(path = c(1L, 2L, 2L), scaffold = c(1L, 1L, 2L)), ignore_attr = TRUE)
  expect_equal(s[[1]]$label, "18S-28S (3 contigs)")
})

test_that("three single-contig paths are one set of three", {
  f <- feats_df(list(1, 1, "18S"), list(2, 1, "18S"), list(3, 1, "18S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(nrow(s[[1]]$units), 3L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
})

test_that("disjoint sets come out largest first and singletons vanish", {
  f <- feats_df(list(1, 1, "28S"), list(2, 1, "28S"), list(1, 2, "18S"), list(2, 2, "18S"),
    list(3, 2, "18S"), list(4, 1, "5.8S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 2L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
  expect_equal(s[[2]]$label, "28S (2 contigs)")
  expect_equal(s[[2]]$units$path, c(1L, 2L))
})

test_that("no anchors or empty input gives an empty list", {
  expect_equal(RiboPilot:::paths_homology_sets(feats_df(list(1, 1, "ITS1"), list(2, 1, "ITS2"))), list())
  expect_equal(RiboPilot:::paths_homology_sets(data.frame(path = integer(), scaffold = integer(),
    element = character())), list())
})
