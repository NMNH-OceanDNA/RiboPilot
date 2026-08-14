# annotate_records() maps tool output back to the right record and the right
# coordinate frame. Most tests below exercise the pieces directly; the last
# block calls annotate_records() itself, with pybarrnap/ITSx faked on PATH
# (same PATH-shim approach as test-rrna-kingdom.R).

test_that("annotate_records demultiplexes barrnap and ITSx output by record", {
  bar <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("18S", "18S"),
    start = c(10L, 20L), stop = c(1810L, 1920L),
    strand = c("+", "+"), model_cov = c(100, 100), stringsAsFactors = FALSE)
  its <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("ITS1", "ITS2"),
    start = c(1811L, 1921L), stop = c(2000L, 2100L), stringsAsFactors = FALSE)
  out <- RiboPilot:::reconcile_records(bar, its)
  expect_equal(sort(unique(out$seqid)), c("ACC1", "ACC2"))
  # the load-bearing property is that a record's coordinates come from its own
  # rows, not from whichever row sorted first
  expect_equal(out$start[out$seqid == "ACC1" & out$element == "18S"], 10L)
  expect_equal(out$start[out$seqid == "ACC2" & out$element == "18S"], 20L)
  expect_true("ITS1" %in% out$element[out$seqid == "ACC1"])
  expect_true("ITS2" %in% out$element[out$seqid == "ACC2"])
})

test_that("minus-strand ITSx coordinates are mapped back to the input frame", {
  # a 1000 bp record; ITSx searched the reverse complement and reported 100-200
  # there, which is 801-901 in the input frame
  flipped <- data.frame(seqid = "ACC1", element = "5.8S",
    start = 100L, stop = 200L, stringsAsFactors = FALSE)
  back <- RiboPilot:::flip_coords(flipped, c(ACC1 = 1000L))
  expect_equal(back$start, 801L)
  expect_equal(back$stop, 901L)
})

test_that("flip_coords leaves records it has no length for untouched", {
  d <- data.frame(seqid = "ACCX", element = "5.8S", start = 5L, stop = 10L,
    stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::flip_coords(d, c(ACC1 = 1000L)), d)
})

test_that("records over the ITSx ceiling are excluded from the ITSx input only", {
  lens <- c(SMALL = 5000L, HUGE = 200000L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "SMALL")
})

test_that("itsx_eligible includes the exact ceiling and excludes one bp over", {
  max_bp <- RiboPilot:::ITSX_MAX_BP
  lens <- c(AT = max_bp, OVER = max_bp + 1L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "AT")
})

test_that("reconciling a minus-strand record in the revcomp frame keeps all five elements", {
  # ACC1 is 5000 bp, stored in non-sense (minus-strand) orientation. barrnap
  # reports anchors in the INPUT frame, descending (28S low ... 18S high, all
  # strand "-"). ITSx was run on the reverse complement and reports its
  # regions in THAT frame (ascending, canonical), unflipped -- this is what
  # annotate_records must do: flip bar into the revcomp frame, reconcile
  # there against ITSx's native output, then flip the finished model back.
  L <- 5000L
  lens <- c(ACC1 = L)
  bar <- data.frame(
    seqid = "ACC1", element = c("28S", "5.8S", "18S"),
    start = c(1L, 3651L, 4001L), stop = c(3450L, 3800L, 5000L),
    strand = "-", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(
    seqid = "ACC1", element = c("ITS1", "5.8S", "ITS2"),
    start = c(1001L, 1201L, 1351L), stop = c(1200L, 1350L, 1550L),
    stringsAsFactors = FALSE)

  bar_rc <- RiboPilot:::flip_coords(bar, lens)
  out_rc <- RiboPilot:::reconcile_records(bar_rc, its)
  out <- RiboPilot:::flip_coords(out_rc, lens)

  expect_setequal(out$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
  expect_true(all(out$start >= 1L & out$stop <= L))
  expect_true(all(out$stop >= out$start))
})

# --- annotate_records() driver, with pybarrnap/ITSx faked on PATH ---

# a canned GFF3 stdout, ignoring argv (mirrors mk_fake_tool in test-rrna-kingdom.R)
mk_fake_pybarrnap <- function(dir, gff_lines) {
  p <- file.path(dir, "pybarrnap")
  writeLines(c("#!/bin/sh",
    sprintf('printf "%%s\\n" %s', paste(vapply(gff_lines, shQuote, character(1)), collapse = " "))), p)
  Sys.chmod(p, "0755")
  p
}

# writes a .positions.txt next to -o, one canned ITS1/5.8S/ITS2 line per
# record found in the -i fasta; ignores actual sequence content, same as the
# fake pybarrnap above -- only the driver's own routing/flip logic is under
# test here, not the tools' real biology
mk_fake_itsx <- function(dir) {
  p <- file.path(dir, "ITSx")
  writeLines(c(
    "#!/bin/sh",
    "INFA=\"\"; PREFIX=\"\"",
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in",
    "    -i) INFA=\"$2\"; shift 2 ;;",
    "    -o) PREFIX=\"$2\"; shift 2 ;;",
    "    *) shift ;;",
    "  esac",
    "done",
    ": > \"${PREFIX}.positions.txt\"",
    "grep '^>' \"$INFA\" | sed 's/^>//; s/ .*//' | while read -r ID; do",
    "  printf '%s\\tITS1: 1801-1850\\t5.8S: 1851-1950\\tITS2: 1951-2000\\n' \"$ID\" >> \"${PREFIX}.positions.txt\"",
    "done",
    "exit 0"), p)
  Sys.chmod(p, "0755")
  p
}

mk_records_fasta <- function(dir) {
  f <- file.path(dir, "records.fasta")
  writeLines(c(
    ">PLUS1", strrep("A", 2000),
    ">MINUS1", strrep("A", 2000),
    ">HUGE1", strrep("A", RiboPilot:::ITSX_MAX_BP + 1L)), f)
  f
}

test_that("annotate_records routes plus/minus groups, strands a minus record '-', and keeps an oversized record anchors-only", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  gff <- c(
    # PLUS1 and HUGE1: barrnap calls 18S "+"
    paste("PLUS1", "barrnap", "rRNA", "10", "1800", "1e-100", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA", sep = "\t"),
    paste("HUGE1", "barrnap", "rRNA", "10", "1800", "1e-100", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA", sep = "\t"),
    # MINUS1: barrnap calls ONLY 5.8S "-" (no 18S/28S anchor). ITSx also
    # calls 5.8S, which is the case reconcile_features() drops barrnap's
    # strand for (IMPORTANT 1); annotate_records must still come back "-".
    paste("MINUS1", "barrnap", "rRNA", "50", "150", "1e-100", "-", ".",
      "Name=5_8S_rRNA;product=5.8S ribosomal RNA", sep = "\t"))
  mk_fake_pybarrnap(dir, gff)
  mk_fake_itsx(dir)
  withr::local_path(dir, action = "prefix")
  fa <- mk_records_fasta(dir)

  out <- RiboPilot:::annotate_records(fa, kingdom = "euk", cpu = 1, taxon = "all")

  # plus/minus split routed each record through the right ITSx call: both
  # PLUS1 (unrevcomp'd) and MINUS1 (revcomp'd, flipped back) picked up ITSx's
  # spacer elements
  expect_true(all(c("ITS1", "5.8S", "ITS2") %in% out$element[out$seqid == "PLUS1"]))
  expect_true(all(c("ITS1", "5.8S", "ITS2") %in% out$element[out$seqid == "MINUS1"]))

  # a minus record comes back strand "-", not the "+" default
  expect_true(all(out$strand[out$seqid == "MINUS1"] == "-"))
  expect_true(all(out$strand[out$seqid == "PLUS1"] == "+"))

  # MINUS1's spacers land at their flipped-back, input-frame coordinates
  # (reversed order vs. PLUS1's raw ITSx coordinates) -- proof the minus
  # group was actually revcomp'd and flipped back, not routed as plus
  m <- out[out$seqid == "MINUS1", ]
  expect_equal(m$start[m$element == "ITS2"], 1L)
  expect_equal(m$stop[m$element == "ITS2"], 50L)
  expect_equal(m$start[m$element == "5.8S"], 51L)
  expect_equal(m$stop[m$element == "5.8S"], 150L)
  expect_equal(m$start[m$element == "ITS1"], 151L)
  expect_equal(m$stop[m$element == "ITS1"], 200L)

  # a record over ITSX_MAX_BP is excluded from ITSx but still appears,
  # anchors only
  h <- out[out$seqid == "HUGE1", ]
  expect_equal(h$element, "18S")
  expect_equal(h$tool, "barrnap")
})
