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

# a stand-in spades.py that fails the way SPAdes does, writing spades.log first
fake_spades <- function(status, env = parent.frame()) {
  bin <- withr::local_tempdir(.local_envir = env)
  writeLines(c("#!/bin/bash", 'out=""',
    'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out=$2; shift; done', 'mkdir -p "$out"',
    'echo "  0:00:01 ... Invalid kmer coverage histogram" > "$out/spades.log"',
    sprintf("echo '== Error ==  system call for: \"[spades-core]\" finished abnormally, OS return value: %d' >> \"$out/spades.log\"", status),
    sprintf("exit %d", status)), file.path(bin, "spades.py"))
  Sys.chmod(file.path(bin, "spades.py"), "0755")
  withr::local_path(bin, action = "prefix", .local_envir = env)
}

test_that("a SPAdes failure is written into the no-assembly marker with its exit status", {
  skip_on_os("windows")
  fake_spades(21L)
  out <- withr::local_tempdir()
  lab <- file.path(out, "labels.fasta"); writeLines(c(">l", "ACGT"), lab)
  n <- suppressMessages(RiboPilot::fallback_assemble_one("S1", "r1.fq", "r2.fq", NA, lab,
    out_dir = out, method = "denovo"))
  expect_equal(as.integer(n), 0L)
  expect_null(attr(n, "retry_status"))
  head <- readLines(file.path(out, "S1_assembly_0.fasta"))[1]
  expect_match(head, "^>No assembly found: SPAdes exit 21: ")
  expect_match(head, "Invalid kmer coverage histogram", fixed = TRUE)
  expect_false(dir.exists(file.path(out, "rp_spades_S1")))
})

test_that("a SPAdes memory kill hands its exit status back for a retry", {
  skip_on_os("windows")
  fake_spades(247L)
  out <- withr::local_tempdir()
  lab <- file.path(out, "labels.fasta"); writeLines(c(">l", "ACGT"), lab)
  n <- suppressMessages(RiboPilot::fallback_assemble_one("S1", "r1.fq", "r2.fq", NA, lab,
    out_dir = out, method = "denovo"))
  expect_equal(attr(n, "retry_status"), 247L)
})

test_that("the fallback script stops asking for a retry after the last attempt", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if_not_installed("pkgload")
  root <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(root, "DESCRIPTION")), "not a source checkout")
  fake_spades(247L)
  d <- withr::local_tempdir()
  lab <- file.path(d, "labels.fasta"); writeLines(c(">l", "ACGT"), lab)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, seeds_db TEXT,
    labels_db TEXT, bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES ('S1', 'default')")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('default', NULL, ?, NULL, NULL)", params = list(lab))
  DBI::dbDisconnect(con)
  withr::local_dir(d)
  run <- function(attempt) suppressWarnings(system2("Rscript",
    c(system.file("nextflow", "bin", "rp_fallback.R", package = "RiboPilot"),
      file.path(d, ".sqlite"), "S1", "r1.fq", "r2.fq", "denovo", "1", attempt),
    env = paste0("RIBOPILOT_SRC=", root), stdout = FALSE, stderr = FALSE))
  expect_equal(run("1"), 247L)
  expect_equal(run("4"), 0L)
  expect_match(readLines("S1_assembly_0.fasta")[1], "SPAdes exit 247")
})

test_that("the fallback refuses to assemble without a usable label database", {
  skip_on_os("windows")
  fake_spades(21L)
  out <- withr::local_tempdir()
  n <- suppressMessages(RiboPilot::fallback_assemble_one("S1", "r1.fq", "r2.fq", NA, NA,
    out_dir = out, method = "denovo"))
  expect_equal(as.integer(n), 0L)
  # SPAdes (stubbed to fail) would have overwritten this with its own exit note
  expect_equal(readLines(file.path(out, "S1_assembly_0.fasta"))[1],
    ">No assembly found: no usable label database to pick out rDNA contigs")
})

test_that("baiting refuses to run without a usable seed database", {
  skip_on_os("windows")
  fake_spades(21L)
  out <- withr::local_tempdir()
  lab <- file.path(out, "labels.fasta"); writeLines(c(">l", "ACGT"), lab)
  n <- suppressMessages(RiboPilot::fallback_assemble_one("S1", "r1.fq", "r2.fq",
    file.path(out, "missing_seeds.fasta"), lab, out_dir = out, method = "bait"))
  expect_equal(as.integer(n), 0L)
  expect_equal(readLines(file.path(out, "S1_assembly_0.fasta"))[1],
    ">No assembly found: no usable seed database to bait reads with")
})

test_that("no contig is kept when the label database cannot be used", {
  d <- withr::local_tempdir()
  contigs <- file.path(d, "contigs.fasta")
  writeLines(c(">NODE_1", strrep("ACGT", 100), ">NODE_2", strrep("TTGA", 100)), contigs)
  out_fa <- file.path(d, "S1_assembly_1.fasta")
  expect_equal(RiboPilot:::extract_rdna_contigs(contigs, NA, "S1", out_fa), 0L)
  expect_false(file.exists(out_fa))
})

test_that("only contigs with a label hit are kept", {
  skip_on_os("windows")
  bin <- withr::local_tempdir()
  writeLines(c("#!/bin/bash", "exit 0"), file.path(bin, "makeblastdb"))
  writeLines(c("#!/bin/bash", "printf 'NODE_1\\t99.0\\t1e-50\\n'"), file.path(bin, "blastn"))
  Sys.chmod(file.path(bin, c("makeblastdb", "blastn")), "0755")
  withr::local_path(bin, action = "prefix")
  d <- withr::local_tempdir()
  contigs <- file.path(d, "contigs.fasta")
  writeLines(c(">NODE_1", strrep("ACGT", 100), ">NODE_2", strrep("TTGA", 100)), contigs)
  lab <- file.path(d, "labels.fasta"); writeLines(c(">l", "ACGT"), lab)
  out_fa <- file.path(d, "S1_assembly_1.fasta")
  expect_equal(RiboPilot:::extract_rdna_contigs(contigs, lab, "S1", out_fa), 1L)
  expect_equal(readLines(out_fa)[1], ">S1.1.1 linear NODE_1")
})
