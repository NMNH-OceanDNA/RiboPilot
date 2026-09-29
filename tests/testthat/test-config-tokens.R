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
  for (f in c("config.local", "config.docker")) {
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
  # sge with no parallel environment leaves <<PENV>>, which neither round fills
  expect_warning(
    suppressMessages(generate_config(name = "sge_nopenv", scheduler = "sge", queue = "q",
      penv = NULL, profile_dir = dir)),
    "by hand")
  expect_no_warning(
    suppressMessages(generate_config(name = "slurm_test", scheduler = "slurm",
      queue = "compute", profile_dir = dir)))
  out <- readLines(file.path(dir, "config.slurm_test"))
  # Only the project-level tokens should be left for new_project().
  expect_true(all(RiboPilot:::config_tokens(out) %in%
    paste0("<<", RiboPilot:::PROJECT_CONFIG_TOKENS, ">>")))
})

test_that("awsbatch is no longer offered", {
  expect_false("awsbatch" %in% RiboPilot:::list_configs(profile_dir = withr::local_tempdir())$name)
  expect_error(generate_config("x", scheduler = "awsbatch", profile_dir = withr::local_tempdir()),
    "should be one of")
})

test_that("a local profile refuses a non-Docker engine", {
  expect_error(
    generate_config("x", scheduler = "local", container_engine = "singularity",
      profile_dir = withr::local_tempdir()),
    "runs the Docker image", fixed = TRUE)
})

test_that("generated container blocks mount the project folder only", {
  for (engine in c("docker", "singularity", "apptainer")) {
    blk <- RiboPilot:::container_engine_block(engine)
    expect_false(grepl("rawDir", blk, fixed = TRUE))
    expect_match(blk, "${params.projectDir}", fixed = TRUE)
  }
})

test_that("singularity profiles set a pull timeout and a shared cache", {
  blk <- RiboPilot:::container_engine_block("singularity", cache = "/shared/sif")
  expect_match(blk, "pullTimeout = '60 min'", fixed = TRUE)
  expect_match(blk, "cacheDir = '/shared/sif'", fixed = TRUE)
  expect_false(grepl("pullTimeout", RiboPilot:::container_engine_block("docker"), fixed = TRUE))

  withr::local_envvar(c(NXF_SINGULARITY_CACHEDIR = "/site/sif", NXF_APPTAINER_CACHEDIR = ""))
  dir <- withr::local_tempdir()
  suppressMessages(generate_config("a", scheduler = "slurm", queue = "q", profile_dir = dir))
  out <- readLines(file.path(dir, "config.a"))
  expect_true(any(grepl("cacheDir = '/site/sif'", out, fixed = TRUE)))

  withr::local_envvar(c(NXF_SINGULARITY_CACHEDIR = "", NXF_APPTAINER_CACHEDIR = ""))
  ref <- withr::local_tempdir()
  withr::local_envvar(c(RIBOPILOT_REF_DIR = ref))
  suppressMessages(generate_config("b", scheduler = "slurm", queue = "q", profile_dir = dir))
  out <- readLines(file.path(dir, "config.b"))
  expect_true(any(grepl(paste0("cacheDir = '", file.path(ref, "containers"), "'"),
    out, fixed = TRUE)))

  for (f in c("config.NMNH_Hydra", "config.NOAA_SEDNA")) {
    txt <- paste(readLines(RiboPilot:::app_sys(f)), collapse = " ")
    expect_match(txt, "pullTimeout = '60 min'", fixed = TRUE, info = f)
  }
})

test_that("a generated profile parses for a project with no reads", {
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  pin <- RiboPilot:::nf_pin_version()
  skip_if(is.na(pin), "no supported Nextflow")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  writeLines(c("workflow {", "  println 'parsed'", "}"), "main.nf")
  cases <- list(c("slurm", "singularity"), c("sge", "singularity"), c("pbs", "singularity"),
    c("lsf", "singularity"), c("slurm", "docker"))
  for (cs in cases) {
    nm <- paste(cs, collapse = "_")
    suppressMessages(generate_config(nm, scheduler = cs[1], container_engine = cs[2],
      queue = "q", container_cache = file.path(dir, "cache"), profile_dir = dir))
    cfg <- RiboPilot:::fill_config(readLines(file.path(dir, paste0("config.", nm))),
      list(CONTAINER_ID = "img", RAW_DIR = "NA", ASMB_DIR = dir, MIN_DEPTH = "0",
        GENETIC_CODE = "1", NCBI_API_KEY = "", EXTRA_MOUNTS = "", HOST_CPUS = 1, HOST_MEMORY = 1))
    writeLines(cfg, paste0(nm, ".config"))
    out <- suppressWarnings(system2("nextflow", c("-q", "-log", "nf.log", "run", "main.nf",
      "-c", paste0(nm, ".config"), "-work-dir", "work", "--projectDir", dir),
      stdout = TRUE, stderr = TRUE, env = paste0("NXF_VER=", pin)))
    expect_true(any(grepl("parsed", out, fixed = TRUE)), info = paste(nm, paste(out, collapse = "\n")))
  }
})
