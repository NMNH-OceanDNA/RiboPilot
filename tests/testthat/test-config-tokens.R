# Config templates carry two rounds of placeholders: scheduler settings filled by
# generate_config(), and project settings filled by new_project(). A template that
# skips the first round produces a config Nextflow cannot parse, so both rounds
# check what is left rather than writing tokens through.

test_that("config_tokens finds unfilled placeholders", {
  expect_equal(RiboPilot:::config_tokens(c("a = '<<QUEUE>>'", "b = 1")), "<<QUEUE>>")
  expect_equal(RiboPilot:::config_tokens("nothing here"), character(0))
  expect_equal(RiboPilot:::config_tokens(c("<<A>>", "<<A>>", "<<B>>")), c("<<A>>", "<<B>>"))
})

test_that("scheduler templates are not fully filled by new_project's token set", {
  # This is the trap: new_project(executor = "slurm") fills only the project
  # tokens, leaving queue / cluster options / container engine behind.
  tmpl <- readLines(RiboPilot:::app_sys("config.slurm"))
  filled <- RiboPilot:::fill_config(tmpl,
    stats::setNames(as.list(rep("x", length(RiboPilot:::PROJECT_CONFIG_TOKENS))),
      RiboPilot:::PROJECT_CONFIG_TOKENS))
  expect_true(length(RiboPilot:::config_tokens(filled)) > 0)
  expect_true("<<QUEUE>>" %in% RiboPilot:::config_tokens(filled))
})

test_that("config.local needs no generate_config round", {
  tmpl <- readLines(RiboPilot:::app_sys("config.local"))
  filled <- RiboPilot:::fill_config(tmpl,
    stats::setNames(as.list(rep("x", length(RiboPilot:::PROJECT_CONFIG_TOKENS))),
      RiboPilot:::PROJECT_CONFIG_TOKENS))
  expect_equal(RiboPilot:::config_tokens(filled), character(0))
})

SCHEDULERS <- c("slurm", "sge", "pbs", "lsf")

test_that("walltime is written as a time directive when asked for", {
  dir <- withr::local_tempdir()
  for (s in SCHEDULERS) {
    suppressMessages(generate_config(name = s, scheduler = s, queue = "q",
      walltime = "24h", profile_dir = dir))
    out <- readLines(file.path(dir, paste0("config.", s)))
    expect_true(any(grepl("^\\s*time = '24h'", out)), info = s)
  }
})

test_that("an unset walltime leaves neither a directive nor a token", {
  dir <- withr::local_tempdir()
  for (s in SCHEDULERS) {
    suppressMessages(generate_config(name = s, scheduler = s, queue = "q",
      profile_dir = dir))
    out <- readLines(file.path(dir, paste0("config.", s)))
    expect_false(any(grepl("time =", out, fixed = TRUE)), info = s)
    expect_false(any(grepl("WALLTIME", out, fixed = TRUE)), info = s)
  }
})

test_that("every scheduler template leaves only project tokens for new_project", {
  dir <- withr::local_tempdir()
  for (s in SCHEDULERS) {
    suppressMessages(generate_config(name = s, scheduler = s, queue = "q",
      walltime = "24h", profile_dir = dir))
    out <- readLines(file.path(dir, paste0("config.", s)))
    expect_true(all(RiboPilot:::config_tokens(out) %in%
      paste0("<<", RiboPilot:::PROJECT_CONFIG_TOKENS, ">>")), info = s)
  }
})

test_that("docker blocks run as the calling user", {
  # otherwise SQLite drops root-owned .sqlite-wal / .sqlite-shm files into the
  # project directory and the Shiny app can no longer open the database
  for (f in c("config.local", "config.docker", "config.awsbatch")) {
    txt <- paste(readLines(RiboPilot:::app_sys(f)), collapse = " ")
    expect_match(txt, "fixOwnership = true", info = f)
    expect_match(txt, "-u \\$(id -u):\\$(id -g)", fixed = TRUE, info = f)
  }
  # saved profiles are generated, not copied, so they need the same treatment
  blk <- RiboPilot:::container_engine_block("docker")
  expect_match(blk, "fixOwnership = true")
  expect_match(blk, "-u \\$(id -u):\\$(id -g)", fixed = TRUE)
})

test_that("generate_config warns about tokens it cannot fill", {
  dir <- withr::local_tempdir()
  # awsbatch carries region / CLI path / batch queue, which neither
  # generate_config() nor new_project() fills.
  expect_warning(
    generate_config(name = "aws_test", scheduler = "awsbatch", profile_dir = dir),
    "by hand")
  expect_no_warning(
    suppressMessages(generate_config(name = "slurm_test", scheduler = "slurm",
      queue = "compute", profile_dir = dir)))
  out <- readLines(file.path(dir, "config.slurm_test"))
  # Only the project-level tokens should be left for new_project().
  expect_true(all(RiboPilot:::config_tokens(out) %in%
    paste0("<<", RiboPilot:::PROJECT_CONFIG_TOKENS, ">>")))
})
