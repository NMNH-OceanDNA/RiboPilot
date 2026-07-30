test_that("the bowtie2 baiting command carries the preset's flags", {
  args <- bt2_args("idx", "a_R1.fq", "a_R2.fq", "out", threads = 8,
    opts = "--very-sensitive-local --score-min G,20,8")
  expect_true(all(c("--very-sensitive-local", "--score-min", "G,20,8") %in% args))
  # the flags the pipeline needs are not user-overridable
  expect_true(all(c("--no-unal", "--al-conc-gz") %in% args))
  expect_equal(args[which(args == "--threads") + 1], "8")
})

test_that("baiting flags default to local mode when the preset is empty", {
  expect_true("--local" %in% bt2_args("idx", "1.fq", "2.fq", "out", 4, opts = ""))
  expect_true("--local" %in% bt2_args("idx", "1.fq", "2.fq", "out", 4, opts = NA))
})

test_that("the SPAdes command carries the preset's flags", {
  args <- spades_args("1.fq", "2.fq", "outdir", threads = 6, opts = "--careful -k 21,33")
  expect_true(all(c("--careful", "-k", "21,33") %in% args))
  expect_false("--isolate" %in% args)   # the preset replaces the default, not adds to it
  expect_equal(args[which(args == "-t") + 1], "6")
})

test_that("SPAdes flags default to isolate mode when the preset is empty", {
  expect_true("--isolate" %in% spades_args("1.fq", "2.fq", "o", 4, opts = ""))
})

test_that("the SPAdes retry drops the preset flags", {
  # a preset flag can be what made the first attempt fail, so the retry is bare
  expect_false(any(c("--careful", "--isolate") %in%
    spades_args("1.fq", "2.fq", "o", 4, opts = "--careful", retry = TRUE)))
})
