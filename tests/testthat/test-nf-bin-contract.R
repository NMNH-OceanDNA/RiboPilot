# Contract every inst/nextflow/bin entry point has to honour.
#
# These scripts are what Nextflow actually runs, so a script that "cleans up"
# with a top-level on.exit() (a no-op outside a function frame) or that accepts
# a missing positional argument as NA reports success over missing data.

bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
bin_files <- list.files(bin_dir, pattern = "\\.R$", full.names = TRUE)
bin_files <- bin_files[basename(bin_files) != "rp_env.R"]

test_that("the bin scripts are found", {
  expect_true(nzchar(bin_dir))
  expect_gt(length(bin_files), 5)
})

test_that("no script relies on a top-level on.exit()", {
  for (f in bin_files) {
    txt <- readLines(f, warn = FALSE)
    expect_false(any(grepl("^(con <-[^;]*;\\s*)?on\\.exit\\(", txt)), info = basename(f))
  }
})

test_that("every script that opens the project DB closes it from a frame", {
  for (f in bin_files) {
    txt <- readLines(f, warn = FALSE)
    if (!any(grepl("open_project_db", txt))) next
    info <- basename(f)
    expect_true(any(grepl("local\\(\\{", txt)), info = info)
    # the on.exit must be indented, i.e. inside the local() frame
    exits <- grep("on\\.exit\\(DBI::dbDisconnect\\(con\\)\\)", txt, value = TRUE)
    expect_length(exits, 1)
    expect_match(exits, "^\\s+on\\.exit", info = info)
  }
})

test_that("every script validates its argument count before doing work", {
  # rp_dbarg.R accepts zero or more arguments by design (no argument means
  # "discover an installed panel") and touches no project database, so the
  # rp_args() requirement below does not apply to it. Derive the set that
  # actually skips rp_args() from the scripts and compare it to that intent, so
  # a new script that forgets rp_args() is a failure here rather than a silent
  # addition to a hand-written list.
  exempt <- "rp_dbarg.R"
  no_rp_args <- vapply(bin_files, function(f) {
    !any(grepl("^rp_args\\(", readLines(f, warn = FALSE)))
  }, logical(1))
  expect_setequal(basename(bin_files)[no_rp_args], exempt)

  for (f in bin_files) {
    if (basename(f) %in% exempt) next
    txt <- readLines(f, warn = FALSE)
    info <- basename(f)
    n_arg <- grep("^rp_args\\(a, [0-9]+,", txt)
    expect_length(n_arg, 1)
    # before the first call that could touch the database
    first_work <- grep("open_project_db|local\\(\\{|build_reference_db", txt)
    if (length(first_work)) expect_lt(n_arg[1], min(first_work), label = info)
  }
})

test_that("every script fails clearly when rp_env.R cannot be located", {
  for (f in bin_files) {
    txt <- readLines(f, warn = FALSE)
    expect_true(any(grepl("length\\(rp_file\\) != 1", txt)), info = basename(f))
  }
})

test_that("rp_ref_db.R rejects an unrecognised mode instead of assuming build", {
  txt <- readLines(file.path(bin_dir, "rp_ref_db.R"), warn = FALSE)
  expect_true(any(grepl('a\\[1\\] %in% c\\("count", "build"\\)', txt)))
})

test_that("a script run with too few arguments exits non-zero", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  script <- file.path(bin_dir, "rp_screen.R")

  err <- tempfile()
  status <- suppressWarnings(system2("Rscript", shQuote(script),
    env = paste0("RIBOPILOT_SRC=", src), stdout = FALSE, stderr = err))
  expect_gt(status, 0)
  expect_match(paste(readLines(err, warn = FALSE), collapse = " "), "usage:")
})

test_that("environment.yml and the Dockerfile pin the same tools to the same versions", {
  root <- RiboPilot:::ribopilot_src_default()
  env_yml <- file.path(root, "environment.yml")
  dockerfile <- file.path(root, "docker", "Dockerfile")
  skip_if(!file.exists(env_yml) || !file.exists(dockerfile), "not a source checkout")

  # "  - name=version" / "  - name>=version" with an optional trailing comment
  parse_specs <- function(lines, pat) {
    m <- regmatches(lines, regexec(pat, lines))
    m <- Filter(function(x) length(x) == 3, m)
    stats::setNames(vapply(m, `[`, "", 3), vapply(m, `[`, "", 2))
  }
  env_lines <- readLines(env_yml, warn = FALSE)
  env_lines <- env_lines[seq(grep("^dependencies:", env_lines)[1] + 1, length(env_lines))]
  env_specs <- parse_specs(env_lines,
    "^\\s*-\\s*([A-Za-z0-9_.-]+)\\s*([<>=][^ #]*)\\s*(?:#.*)?$")

  # nothing in the tool chain left floating
  unpinned <- grep("^\\s*-\\s*[A-Za-z0-9_.-]+\\s*(#.*)?$", env_lines, value = TRUE)
  expect_length(unpinned, 0)

  docker_lines <- readLines(dockerfile, warn = FALSE)
  docker_specs <- parse_specs(docker_lines, "^\\s*'([A-Za-z0-9_.-]+)([<>=][^']*)'")
  # nothing installed by name alone
  expect_length(grep("^\\s*'[A-Za-z0-9_.-]+'\\s*(\\\\|&&)", docker_lines, value = TRUE), 0)
  # exact versions, not ranges
  expect_true(all(startsWith(docker_specs, "=")))
  tools <- docker_specs[!grepl("^(r|bioconductor)-", names(docker_specs))]
  expect_setequal(names(tools), setdiff(names(env_specs), "nextflow"))
  for (p in names(tools)) expect_equal(unname(env_specs[p]), unname(tools[p]), info = p)
})

test_that("bin scripts reach package internals the way a container will", {
  ns <- asNamespace("RiboPilot")
  ex <- getNamespaceExports("RiboPilot")
  for (f in list.files(bin_dir, pattern = "\\.R$", full.names = TRUE)) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    fc <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
    qualified <- fc$parent %in% pd$parent[pd$token %in% c("NS_GET", "NS_GET_INT")]
    bare <- unique(fc$text[!qualified])
    bad <- bare[vapply(bare, exists, NA, envir = ns, inherits = FALSE) & !bare %in% ex]
    expect_length(bad, 0)
  }
})

test_that("format_assembly.sh numbers paths by GetOrganelle's path number", {
  skip_on_os("windows")
  sh <- file.path(bin_dir, "format_assembly.sh")
  wd <- withr::local_tempdir()
  go <- file.path(wd, "go_out")
  dir.create(go)
  # glob order would put graph1.10 before graph1.2
  for (i in c(1, 2, 10)) {
    writeLines(c(paste0(">scaffold_", i, "--e1"), "ACGT"),
      file.path(go, sprintf("embplant_nr.K105.complete.graph1.%d.path_sequence.fasta", i)))
  }
  status <- withr::with_dir(wd,
    system2("bash", c(shQuote(sh), "SAMP", "go_out"), stdout = FALSE, stderr = FALSE))
  expect_equal(status, 0)
  first <- vapply(1:3, function(n) {
    readLines(file.path(wd, sprintf("SAMP_assembly_%d.fasta", n)), warn = FALSE)[1]
  }, "")
  expect_match(first[1], "scaffold_1--e1")
  expect_match(first[2], "scaffold_2--e1")
  expect_match(first[3], "scaffold_10--e1")
})

test_that("format_assembly.sh writes GetOrganelle's errors into the no-assembly marker", {
  skip_on_os("windows")
  sh <- file.path(bin_dir, "format_assembly.sh")
  wd <- withr::local_tempdir()
  go <- file.path(wd, "go_out")
  dir.create(go)
  writeLines(c(
    "2026-09-23 16:47:12,102 - WARNING: Guessing that you are using too few data for assembling anonym!",
    "2026-09-23 16:47:12,102 - WARNING: GetOrganelle is still trying ...",
    "2026-09-23 16:48:00,020 - ERROR: Error with running SPAdes: == Error ==  system call for: \"['/opt/conda/bin/spades-core', '/w/K115/config.info']\" finished abnormally, OS return value: 21",
    "2026-09-23 16:48:00,030 - ERROR: Assembling failed."), file.path(go, "get_org.log.txt"))
  status <- withr::with_dir(wd,
    system2("bash", c(shQuote(sh), "SAMP", "go_out"), stdout = FALSE, stderr = FALSE))
  expect_equal(status, 0)
  expect_equal(readLines(file.path(wd, "SAMP_assembly_0.fasta"))[1], paste0(
    ">No assembly found: Guessing that you are using too few data for assembling anonym!; ",
    "Error with running SPAdes: == Error ==  finished abnormally, OS return value: 21; ",
    "Assembling failed."))
  # no GetOrganelle log: the plain marker, as before
  unlink(file.path(go, "get_org.log.txt"))
  withr::with_dir(wd, system2("bash", c(shQuote(sh), "SAMP", "go_out"), stdout = FALSE, stderr = FALSE))
  expect_equal(readLines(file.path(wd, "SAMP_assembly_0.fasta"))[1], ">No assembly found")
})
