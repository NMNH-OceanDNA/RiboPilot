mk_ref_cache <- function() {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(d, "seeds")); dir.create(file.path(d, "blastdb"))
  writeLines(">a\nACGT", file.path(d, "seeds", "x_seeds.fasta"))
  writeLines("taxid\tn_records\n1\t2", file.path(d, "counts.tsv"))
  writeLines("panel", file.path(d, "blastdb", "panel.nsq"))
  d
}

test_that("clear_ref_cache warns about every project and refuses without a person", {
  d <- mk_ref_cache()
  expect_error(expect_message(clear_ref_cache(d), "EVERY RiboPilot project"), "interactive R session")
  expect_true(file.exists(file.path(d, "seeds", "x_seeds.fasta")))
})

test_that("clear_ref_cache deletes only on the typed word, and leaves the panel", {
  d <- mk_ref_cache()
  local_mocked_bindings(interactive = function() TRUE, readline = function(...) "yes")
  expect_length(suppressMessages(clear_ref_cache(d)), 0)
  expect_true(file.exists(file.path(d, "counts.tsv")))

  local_mocked_bindings(readline = function(...) "clear")
  expect_setequal(basename(suppressMessages(clear_ref_cache(d))), c("counts.tsv", "seeds"))
  expect_equal(list.files(d), "blastdb")
  expect_message(clear_ref_cache(d), "Nothing to clear")
})
