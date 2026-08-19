# What nf_command() puts in the argument vector, and how the cluster submit
# scripts turn that vector into a shell command line.

mk_cfg_project <- function(dir, lines) {
  writeLines(lines, file.path(dir, ".config"))
  dir
}

# Let bash do the word splitting, so the assertions are about what the scheduler
# will actually hand to Nextflow rather than about the text of the script.
shell_words <- function(script) {
  line <- grep("^'nextflow'", strsplit(script, "\n", fixed = TRUE)[[1]], value = TRUE)
  expect_length(line, 1)
  line <- paste0("printf '%s\\n' ", sub("^'nextflow' ", "", line))
  system2("bash", c("-c", shQuote(line)), stdout = TRUE)
}

test_that("a config missing minDepth does not pass the string NA", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "}"))
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_false("NA" %in% a)
  expect_equal(a[which(a == "--minDepth") + 1], "0")
})

test_that("a minDepth left as an unfilled NA token does not pass the string NA", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "  minDepth = NA", "}"))
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_equal(a[which(a == "--minDepth") + 1], "0")
  expect_false("NA" %in% a)
})

test_that("a config missing rawDir fails at command-build time", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  minDepth = 1000", "}"))
  # --rawDir NA makes docker create a host directory called NA, and every sample
  # then finds no reads hours into the run
  expect_error(RiboPilot:::nf_command(proj, resume = FALSE), "rawDir")
})

test_that("an absent ncbi_api_key is dropped rather than passed as NA", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "}"))
  expect_false("--ncbi_api_key" %in% RiboPilot:::nf_command(proj, resume = FALSE))
})

test_that("the submit script keeps a multi-word argument as one shell word", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  s <- file.path(refs, "x_seeds.fasta"); l <- file.path(refs, "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_cfg_project(proj, c(RiboPilot:::container_engine_block("docker"),
    "params {", "  rawDir = '/data/raw'", "}"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (assemble_opts TEXT, seeds_db TEXT, labels_db TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('a', ?, ?)", params = list(s, l))
  DBI::dbDisconnect(con)

  binds <- RiboPilot:::nf_command(proj, resume = FALSE)
  binds <- binds[which(binds == "-process.containerOptions") + 1]
  expect_match(binds, " ")   # the element this test exists for

  words <- shell_words(RiboPilot:::build_submit_script(proj, resume = FALSE, sched = NULL))
  i <- which(words == "-process.containerOptions")
  expect_length(i, 1)
  expect_equal(words[i + 1], binds)
})

test_that("a project path containing a space survives into the submit script", {
  skip_if_not_installed("withr")
  proj <- file.path(withr::local_tempdir(), "my project")
  dir.create(proj)
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "}"))

  script <- RiboPilot:::build_submit_script(proj, resume = FALSE, sched = NULL)
  words <- shell_words(script)
  expect_true(proj %in% words)
  expect_equal(words[which(words == "--projectDir") + 1], proj)
  expect_true(grepl(paste0("cd ", shQuote(proj)), script, fixed = TRUE))
})

test_that("-ansi-log false stays two arguments", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/data/raw'", "}"))
  words <- shell_words(RiboPilot:::build_submit_script(proj, resume = FALSE, sched = NULL))
  expect_equal(utils::tail(words, 2), c("-ansi-log", "false"))
})

test_that("an unconfigured project passes no --blastdb or --ref_dir flag at all", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/reads'", "}"))
  args <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_false("--blastdb" %in% args)
  expect_false("--ref_dir" %in% args)
  expect_false(any(args == ""))
})

test_that("a configured panel still passes --blastdb with its value", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_cfg_project(proj, c("params {", "  rawDir = '/reads'", "}"))
  con <- RiboPilot:::open_project_db(file.path(proj, ".sqlite"))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
    VALUES (1, NULL, 1000, 25000)")
  RiboPilot:::set_ref_config(con, blastdb = "/panels/2026.08.1/ribopilot_rrna")
  DBI::dbDisconnect(con)
  args <- RiboPilot:::nf_command(proj, resume = FALSE)
  i <- which(args == "--blastdb")
  expect_length(i, 1)
  expect_equal(args[i + 1], "/panels/2026.08.1/ribopilot_rrna")
})

test_that("the reference submit script quotes per argument too", {
  skip_if_not_installed("withr")
  parent <- withr::local_tempdir()
  proj <- file.path(parent, "my project"); dir.create(proj)
  refs <- file.path(parent, "ref dbs"); dir.create(refs)
  mk_cfg_project(proj, c(RiboPilot:::container_engine_block("docker"),
    "params {", "  rawDir = '/data/raw'", "}"))
  manifest <- file.path(proj, "clades.tsv")
  writeLines("clade\tquery", manifest)

  words <- shell_words(RiboPilot:::build_reference_submit_script(proj, manifest, refs,
    sched = NULL))
  expect_equal(words[which(words == "--refDir") + 1], refs)
  expect_equal(words[which(words == "--clades") + 1], manifest)
  expect_equal(words[which(words == "--projectDir") + 1], proj)
  expect_equal(utils::tail(words, 2), c("-ansi-log", "false"))
})
