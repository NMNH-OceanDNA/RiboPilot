test_that(".meta_flatten keeps large round numbers in plain digits", {
  out <- .meta_flatten(list(gbifID = 100000, country = "Peru"), "Occurrence", 0, "r")
  expect_equal(out$value, c("100000", "Peru"))
})

test_that(".meta_chr keeps long IDs in plain digits and NA as NA", {
  expect_equal(.meta_chr(c(6186461308, NA, 12, -170.58225)), c("6186461308", NA, "12", "-170.58225"))
  expect_equal(.meta_chr(c("a", NA)), c("a", NA))
})
