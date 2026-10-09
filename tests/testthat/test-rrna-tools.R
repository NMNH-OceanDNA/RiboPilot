# inst/nextflow/bin/rp_rrna_tools.sh runs barrnap or pybarrnap, the tuned ITSx
# run, and the stock ITSx check on one contig. The tools are stubs on PATH;
# each stub appends its argv to <name>.calls, one call per line, tab-separated.

rrna_sh <- function() system.file("nextflow", "bin", "rp_rrna_tools.sh", package = "RiboPilot")

mk_tool_stub <- function(dir, name, body = character(0)) {
  calls <- shQuote(file.path(dir, paste0(name, ".calls")))
  p <- file.path(dir, name)
  writeLines(c("#!/usr/bin/env bash",
    sprintf("printf '%%s\\t' \"$@\" >> %s; echo >> %s", calls, calls), body), p)
  Sys.chmod(p, "0755")
  p
}

stub_calls <- function(dir, name) {
  f <- file.path(dir, paste0(name, ".calls"))
  if (!file.exists(f)) return(list())
  lapply(readLines(f), function(l) strsplit(l, "\t", fixed = TRUE)[[1]])
}

run_rrna_sh <- function(dir, args) {
  processx::run("bash", c(rrna_sh(), args), wd = dir, error_on_status = FALSE,
    env = c("current", PATH = paste(dir, Sys.getenv("PATH"), sep = .Platform$path.sep)))
}

mk_unit_fasta <- function(dir) {
  writeLines(c(">S1.1.1", strrep("ACGT", 50)), file.path(dir, "unit.fasta"))
}

ITSX_OUT <- 'while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; esac; shift; done'

test_that("rp_rrna_tools.sh is executable and refuses a wrong argument count", {
  skip_on_os("windows")
  expect_equal(unname(file.access(rrna_sh(), 1)), 0L)
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap")
  mk_tool_stub(dir, "ITSx")
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "--a", "b", "", "2"))
  expect_equal(res$status, 2L)
  expect_match(res$stderr, "usage:", fixed = TRUE)
  expect_length(stub_calls(dir, "pybarrnap"), 0)
  expect_false(file.exists(file.path(dir, "unit.tools_status")))
})

test_that("rp_rrna_tools.sh runs pybarrnap and both ITSx runs with the pipeline's arguments", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap",
    "printf 'S1.1.1\\tbarrnap:0.9\\trRNA\\t1\\t180\\t0\\t+\\t.\\tName=18S_rRNA\\n'")
  mk_tool_stub(dir, "ITSx", c(ITSX_OUT,
    '[ "$out" = unit.itsx ] && echo "S1.1.1 tuned" > "$out.positions.txt"',
    'echo "$out note" >&2', "exit 0"))
  # empty tool, kingdom, and taxon take the defaults; both option strings word-split
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "--allow_single_domain 1e-9,0 -E 1e-3",
    "--evalue 1e-5 --lencutoff 0.5", "2"))
  expect_equal(res$status, 0L)
  expect_equal(stub_calls(dir, "pybarrnap"), list(c("-k", "euk", "-t", "2", "-q",
    "--evalue", "1e-5", "--lencutoff", "0.5", "unit.fasta")))
  expect_equal(stub_calls(dir, "ITSx"), list(
    c("-i", "unit.fasta", "-o", "unit.itsx", "--save_regions", "none", "--graphical", "F",
      "--cpu", "2", "-t", "all", "--preserve", "T", "--allow_single_domain", "1e-9,0", "-E", "1e-3"),
    c("-i", "unit.fasta", "-o", "unit.itsx_check", "--preserve", "T", "--detailed_results", "T",
      "--save_regions", "none", "--graphical", "F", "--cpu", "2")))
  expect_match(readLines(file.path(dir, "unit.barrnap.gff")), "Name=18S_rRNA", fixed = TRUE)
  expect_equal(readLines(file.path(dir, "unit.itsx.positions.txt")), "S1.1.1 tuned")
  expect_equal(readLines(file.path(dir, "unit.itsx_check.positions.txt")), character(0))
  # only the stock run's stderr is kept
  expect_equal(readLines(file.path(dir, "unit.itsx_check.err")), "unit.itsx_check note")
  expect_match(res$stderr, "unit.itsx note", fixed = TRUE)
  expect_equal(readLines(file.path(dir, "unit.tools_status")),
    c("barrnap=0", "itsx=0", "itsx_check=0"))
})

test_that("a barrnap failure sets the exit status after both ITSx runs", {
  skip_on_os("windows")
  for (case in list(list(tool = "barrnap", code = 1L), list(tool = "pybarrnap", code = 137L))) {
    dir <- withr::local_tempdir()
    mk_unit_fasta(dir)
    mk_tool_stub(dir, case$tool, sprintf("exit %d", case$code))
    mk_tool_stub(dir, "ITSx")
    res <- run_rrna_sh(dir, c("unit.fasta", case$tool, "euk", "", "", "--evalue 1e-5", "1"))
    expect_equal(res$status, case$code, info = case$tool)
    expect_equal(stub_calls(dir, case$tool)[[1]], if (case$tool == "barrnap") {
      c("--kingdom", "fun", "--quiet", "--evalue", "1e-5", "unit.fasta")
    } else {
      c("-k", "euk", "-t", "1", "-q", "--evalue", "1e-5", "unit.fasta")
    })
    expect_length(stub_calls(dir, "ITSx"), 2)
    expect_equal(readLines(file.path(dir, "unit.tools_status")),
      c(paste0("barrnap=", case$code), "itsx=0", "itsx_check=0"), info = case$tool)
    expect_true(file.exists(file.path(dir, "unit.itsx.positions.txt")))
    expect_true(file.exists(file.path(dir, "unit.itsx_check.positions.txt")))
  }
})

test_that("an ITSx failure is recorded but never fails the script", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_unit_fasta(dir)
  mk_tool_stub(dir, "pybarrnap")
  mk_tool_stub(dir, "ITSx", c(ITSX_OUT,
    '[ "$out" = unit.itsx ] && exit 2',
    'echo "FATAL ERROR: hmmscan failed" >&2', "exit 3"))
  res <- run_rrna_sh(dir, c("unit.fasta", "", "", "", "", "", "1"))
  expect_equal(res$status, 0L)
  expect_equal(readLines(file.path(dir, "unit.tools_status")),
    c("barrnap=0", "itsx=2", "itsx_check=3"))
  expect_equal(readLines(file.path(dir, "unit.itsx.positions.txt")), character(0))
  expect_equal(readLines(file.path(dir, "unit.itsx_check.positions.txt")), character(0))
  expect_equal(readLines(file.path(dir, "unit.itsx_check.err")), "FATAL ERROR: hmmscan failed")
})

test_that("annotate_tools hands the tools to rp_rrna_tools.sh and keeps its outputs", {
  nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))
  start <- grep("^process annotate_tools \\{", nf)
  expect_length(start, 1)
  body <- nf[start:(start - 1 + which(nf[start:length(nf)] == "}")[1])]
  expect_equal(trimws(grep("^\\s*rp_rrna_tools\\.sh ", body, value = TRUE)),
    'rp_rrna_tools.sh unit.fasta "\\$TOOL" "\\$KINGDOM" "\\$ITSXT" "\\$ITSXO" "\\$RRNAO" ${task.cpus}')
  expect_false(any(grepl("barrnap --kingdom|pybarrnap -k|ITSx -i|tools_status", body)))
  expect_true(any(grepl('path("unit.itsx_check*")', body, fixed = TRUE)))
})

test_that("annotate_tools fails on a barrnap error and retries an out-of-memory exit", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
  if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
  skip_if(!nzchar(Sys.which("sqlite3")), "sqlite3 not installed")
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  con <- chimera_test_con("S1")
  db <- con@dbname
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence)
    VALUES ('S1', 1, 1, 'ACGTACGTAC'), ('S1', 1, 2, 'ACGTACGTAC')")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold, annotate_opts)
    VALUES ('S1', 1, 1, 'default'), ('S1', 1, 2, 'default')")
  DBI::dbDisconnect(con)

  root <- withr::local_tempdir()
  seen <- shQuote(file.path(root, "seen"))
  # contig 1 is killed for memory on its first try; contig 2 always fails
  mk_tool_stub(root, "pybarrnap", c(
    'u=$(head -n 1 "${@: -1}" | cut -c2-)',
    sprintf('echo "$u" >> %s', seen),
    sprintf('n=$(grep -cx "$u" %s)', seen),
    'case "$u.$n" in S1.1.1.1) exit 137 ;; S1.1.2.*) exit 1 ;; esac'))
  mk_tool_stub(root, "ITSx")
  writeLines(c(
    sprintf("include { annotate_tools } from '%s'", file.path(nf_dir, "modules", "pipeline_r.nf")),
    "workflow {",
    "    annotate_tools(Channel.of(tuple('S1', 1, 1, 1, 1), tuple('S1', 1, 2, 1, 1)))",
    "    annotate_tools.out.view { \"EMITTED ${it[0]}.${it[1]}.${it[2]}\" }",
    "}"), file.path(root, "main.nf"))

  res <- processx::run("nextflow", c("run", file.path(root, "main.nf"),
    "-c", file.path(nf_dir, "nextflow.config"), "-work-dir", file.path(root, "work"),
    "-ansi-log", "false", "--db", db), wd = root,
    env = c(RiboPilot:::nf_env(),
      PATH = paste(root, file.path(nf_dir, "bin"), Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE, timeout = 300, cleanup_tree = TRUE)
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(sort(readLines(file.path(root, "seen"))), c("S1.1.1", "S1.1.1", "S1.1.2"))
  expect_match(res$stdout, "EMITTED S1.1.1", fixed = TRUE)
  expect_no_match(res$stdout, "EMITTED S1.1.2", fixed = TRUE)
})
