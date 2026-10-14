# The two rRNA callers take different kingdom vocabularies (pybarrnap:
# bac/arc/euk/all; the bioconda Perl barrnap fork: bac/arc/fun, where "fun" is
# its eukaryote Rfam models). Passing the stored kingdom through untranslated
# makes barrnap exit non-zero and every scaffold gets zero anchors. External
# binaries are faked on PATH so nothing here needs barrnap or pybarrnap.

GFF_18S <- paste("seq1", "barrnap", "rRNA", "10", "1800", "1e-100", "+", ".",
  "Name=18S_rRNA;product=18S ribosomal RNA", sep = "\t")

# A stand-in binary that records the argv it was called with and prints `out`.
mk_fake_tool <- function(dir, name, out = character(0), status = 0L) {
  p <- file.path(dir, name)
  lines <- c(
    "#!/bin/sh",
    sprintf('printf "%%s\\n" "$*" > %s', shQuote(file.path(dir, paste0(name, ".args")))))
  if (length(out) > 0) {
    lines <- c(lines, sprintf('printf "%%s\\n" %s',
      paste(vapply(out, shQuote, character(1)), collapse = " ")))
  }
  if (status != 0L) {
    lines <- c(lines, sprintf('echo "%s: bad kingdom" >&2', name), sprintf("exit %d", status))
  }
  writeLines(lines, p)
  Sys.chmod(p, "0755")
  p
}

fake_args <- function(dir, name) {
  f <- file.path(dir, paste0(name, ".args"))
  if (!file.exists(f)) NA_character_ else readLines(f, warn = FALSE)[1]
}

mk_fasta <- function(dir) {
  f <- file.path(dir, "unit.fasta")
  writeLines(c(">seq1", strrep("ACGT", 100)), f)
  f
}

test_that("rrna_kingdom_for maps the eukaryote intent per tool", {
  f <- RiboPilot:::rrna_kingdom_for
  # the Perl fork has no euk/all; "fun" is its eukaryote Rfam model set
  expect_equal(f("euk", "barrnap"), "fun")
  expect_equal(f("all", "barrnap"), "fun")
  expect_equal(f("euk", "pybarrnap"), "euk")
  expect_equal(f("all", "pybarrnap"), "all")
  # pybarrnap has no "fun"
  expect_equal(f("fun", "pybarrnap"), "euk")
  expect_equal(f("fun", "barrnap"), "fun")
})

test_that("rrna_kingdom_for passes bac/arc through for both tools", {
  f <- RiboPilot:::rrna_kingdom_for
  for (k in c("bac", "arc")) {
    expect_equal(f(k, "barrnap"), k)
    expect_equal(f(k, "pybarrnap"), k)
  }
})

test_that("rrna_kingdom_for falls back to the eukaryote intent when unset", {
  f <- RiboPilot:::rrna_kingdom_for
  expect_equal(f(NA_character_, "barrnap"), "fun")
  expect_equal(f("", "pybarrnap"), "euk")
  expect_equal(f("EUK", "barrnap"), "fun")
})

test_that("a barrnap failure warns instead of returning a silent empty frame", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_fake_tool(dir, "barrnap", status = 2L)
  withr::local_path(dir)
  fa <- mk_fasta(dir)

  expect_warning(out <- RiboPilot:::run_barrnap(fa, kingdom = "euk"),
    "barrnap.*kingdom 'fun'.*failed")
  expect_equal(nrow(out), 0)
})

test_that("an empty result warns too, naming the tool and kingdom", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  bin <- mk_fake_tool(dir, "pybarrnap")
  fa <- mk_fasta(dir)
  expect_warning(RiboPilot:::run_rrna_tool(bin, "pybarrnap", fa, "euk"),
    "pybarrnap.*kingdom 'euk'.*no rRNA records")
})

test_that("a successful call is quiet and parses", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  bin <- mk_fake_tool(dir, "barrnap", out = GFF_18S)
  fa <- mk_fasta(dir)
  expect_silent(gff <- RiboPilot:::run_rrna_tool(bin, "barrnap", fa, "fun"))
  expect_equal(RiboPilot:::parse_barrnap(gff)$element, "18S")
})

test_that("run_barrnap translates euk to fun on the command line", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_fake_tool(dir, "barrnap", out = GFF_18S)
  withr::local_path(dir)
  fa <- mk_fasta(dir)

  expect_equal(RiboPilot:::run_barrnap(fa, kingdom = "euk")$element, "18S")
  expect_match(fake_args(dir, "barrnap"), "--kingdom fun")
})

test_that("call_rrna_anchors leaves euk alone for pybarrnap", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_fake_tool(dir, "pybarrnap", out = GFF_18S)
  withr::local_path(dir)
  fa <- mk_fasta(dir)

  expect_equal(RiboPilot:::call_rrna_anchors(fa, kingdom = "euk")$element, "18S")
  expect_match(fake_args(dir, "pybarrnap"), "-k euk")
})

test_that("rrna_caller prefers pybarrnap and falls back to barrnap", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  mk_fake_tool(dir, "barrnap")
  withr::local_path(dir)
  expect_equal(RiboPilot:::rrna_caller()$tool, "barrnap")

  mk_fake_tool(dir, "pybarrnap")
  expect_equal(RiboPilot:::rrna_caller()$tool, "pybarrnap")
})

test_that("assembly_anchors and the annotation path agree on tool and kingdom", {
  skip_on_os("windows")
  # one PATH override per tool, so the pybarrnap fake cannot shadow the
  # barrnap-only case
  check_agree <- function(tool) {
    dir <- withr::local_tempdir()
    mk_fake_tool(dir, tool, out = GFF_18S)
    withr::local_path(dir)
    fa <- mk_fasta(dir)

    expect_equal(RiboPilot:::assembly_anchors(fa, kingdom = "euk"), "18S")
    scorer_args <- fake_args(dir, tool)
    expect_equal(RiboPilot:::call_rrna_anchors(fa, kingdom = "euk")$element, "18S")
    annot_args <- fake_args(dir, tool)

    expect_equal(scorer_args, annot_args)
    expect_match(scorer_args, if (tool == "barrnap") "fun" else "euk")
  }
  check_agree("pybarrnap")
  check_agree("barrnap")
})

test_that("assembly_anchors is empty and quiet when no caller is installed", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(rrna_caller = function() NULL, .package = "RiboPilot")
  expect_equal(RiboPilot:::assembly_anchors(mk_fasta(dir)), character(0))
})

test_that("rp_rrna_tools.sh translates the kingdom as rrna_kingdom_for does", {
  skip_on_os("windows")
  sh <- system.file("nextflow", "bin", "rp_rrna_tools.sh", package = "RiboPilot")
  kingdom_passed <- function(tool, kingdom) {
    dir <- withr::local_tempdir()
    mk_fake_tool(dir, tool, out = GFF_18S)
    mk_fake_tool(dir, "ITSx")
    withr::local_path(dir)
    mk_fasta(dir)
    status <- withr::with_dir(dir, system2("bash",
      shQuote(c(sh, "unit.fasta", tool, kingdom, "", "", "", "1")), stdout = FALSE, stderr = FALSE))
    expect_equal(status, 0)
    args <- strsplit(fake_args(dir, tool), " ", fixed = TRUE)[[1]]
    args[which(args %in% c("--kingdom", "-k")) + 1]
  }
  for (tool in c("barrnap", "pybarrnap")) {
    for (k in c("euk", "all", "fun", "bac", "arc", "")) {
      expect_equal(kingdom_passed(tool, k), RiboPilot:::rrna_kingdom_for(k, tool),
        info = paste(tool, k))
    }
  }
})

test_that("run_itsx warns when ITSx exits non-zero", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 1L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  expect_warning(RiboPilot:::run_itsx(fa), "ITSx")
})

test_that("run_itsx warns when no positions file is produced", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  expect_warning(RiboPilot:::run_itsx(fa), "positions")
})

test_that("run_itsx passes the requested taxon through to -t", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa, taxon = "M"))
  args <- fake_args(d, "ITSx")
  expect_match(args, "-t M")
})

test_that("ITSX_MAX_BP is the documented HMMER ceiling", {
  expect_equal(RiboPilot:::ITSX_MAX_BP, 100000)
})

test_that("run_itsx omits --complement unless asked", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa))
  expect_false(grepl("--complement", fake_args(d, "ITSx")))
})

test_that("run_itsx passes --complement F when asked", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- mk_fasta(d)
  suppressWarnings(RiboPilot:::run_itsx(fa, complement = FALSE))
  expect_match(fake_args(d, "ITSx"), "--complement F")
})

test_that("annotate_records asks ITSx not to search the complement", {
  d <- withr::local_tempdir()
  mk_fake_tool(d, "pybarrnap", status = 0L)
  mk_fake_tool(d, "ITSx", status = 0L)
  withr::local_path(d, action = "prefix")
  fa <- file.path(d, "refs.fa")
  writeLines(c(">ACC1.1 x", paste(rep("ACGT", 60), collapse = "")), fa)
  suppressWarnings(RiboPilot:::annotate_records(fa))
  expect_match(fake_args(d, "ITSx"), "--complement F")
})
