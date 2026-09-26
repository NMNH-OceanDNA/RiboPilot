# Selecting nuclear rDNA contigs out of a user-supplied assembly.
# See dev/design/20_USER_ASSEMBLIES.md.

gff_line <- function(seqid, name, start, stop) {
  paste(seqid, "barrnap", "rRNA", start, stop, "1e-100", "+", ".",
    paste0("Name=", name, ";product=", name), sep = "\t")
}

test_that("a contig with an 18S anchor is kept and one without is not", {
  gff <- c(
    gff_line("ctg1", "18S_rRNA", 1, 1800),
    gff_line("ctg1", "28S_rRNA", 3000, 6500),
    gff_line("ctg2", "5_8S_rRNA", 100, 250)
  )
  res <- select_rrna_contigs(gff, lengths = c(ctg1 = 12000, ctg2 = 5000))

  expect_equal(res$candidates, "ctg1")
  expect_equal(nrow(res$evidence), 2L)
  expect_equal(res$evidence$selected[res$evidence$contig == "ctg2"], 0L)
  expect_match(res$evidence$reason[res$evidence$contig == "ctg2"], "no 18S or 28S")
  expect_equal(res$evidence$length[res$evidence$contig == "ctg1"], 12000L)
  expect_match(res$evidence$rrna_features[res$evidence$contig == "ctg1"], "18S,28S")
})

test_that("a single supplied rDNA contig passes the scan", {
  res <- select_rrna_contigs(gff_line("rDNA_unit", "28S_rRNA", 500, 4000))
  expect_equal(res$candidates, "rDNA_unit")
  expect_true(is.na(res$evidence$length))
})

test_that("no rRNA hits at all yields no candidates and no evidence", {
  res <- select_rrna_contigs(character(0))
  expect_length(res$candidates, 0L)
  expect_equal(nrow(res$evidence), 0L)

  res <- select_rrna_contigs("# nothing here")
  expect_length(res$candidates, 0L)
  expect_equal(nrow(res$evidence), 0L)
})

test_that("contigs rank by rRNA content and the cap keeps the best", {
  gff <- c(
    gff_line("small", "18S_rRNA", 1, 500),
    gff_line("big", "18S_rRNA", 1, 1800),
    gff_line("big", "28S_rRNA", 2000, 5500),
    gff_line("mid", "18S_rRNA", 1, 1700)
  )
  res <- select_rrna_contigs(gff, max_contigs = 2)

  expect_equal(res$candidates, c("big", "mid"))
  expect_equal(res$evidence$rank, 1:3)
  expect_match(res$evidence$reason[res$evidence$contig == "small"], "beyond the 2 best")
})

test_that("max_contigs of NA or zero means no cap", {
  gff <- c(gff_line("a", "18S_rRNA", 1, 1800), gff_line("b", "28S_rRNA", 1, 3000))
  expect_length(select_rrna_contigs(gff, max_contigs = NA)$candidates, 2L)
  expect_length(select_rrna_contigs(gff, max_contigs = 0)$candidates, 2L)
})

test_that("a scan part that never came back fails the sample instead of shrinking the search", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if_not_installed("pkgload")
  root <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(root, "DESCRIPTION")), "not a source checkout")
  d <- withr::local_tempdir()
  db <- file.path(d, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_opts TEXT, paths INT, scaffolds INT,
    length TEXT, topology TEXT, assemble_switch INT, assemble_notes TEXT, find_rrna_notes TEXT,
    time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, max_scaffolds INT)")
  DBI::dbExecute(con, "CREATE TABLE rrna_candidates(ID TEXT, contig TEXT, length INT,
    rrna_features TEXT, rrna_bases INT, rank INT, selected INT, reason TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, annotate_switch INT, annotate_lock INT,
    reviewed TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default', 10)")
  DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_opts, assemble_switch) VALUES('S1', 'default', 1)")
  withr::local_dir(d)
  s <- strrep("ACGT", 20)
  writeLines(c(">ctg1", s, ">ctg2", s, ">ctg3", s), "asm.fa")
  for (k in 1:3) writeLines(gff_line(paste0("ctg", k), "18S_rRNA", 1, 40), sprintf("chunk_%d.gff", k))
  run <- function(...) suppressWarnings(system2("Rscript",
    c(system.file("nextflow", "bin", "rp_find_rrna.R", package = "RiboPilot"),
      db, "S1", "asm.fa", "S1", "1", ...),
    env = paste0("RIBOPILOT_SRC=", root), stdout = FALSE, stderr = FALSE))
  state <- function() DBI::dbGetQuery(con,
    "SELECT assemble_switch, find_rrna_notes FROM assemble WHERE ID = 'S1'")

  # chunk size 1 over three contigs is three parts; the second never came back
  expect_equal(run("chunk_1.gff", "chunk_3.gff"), 0L)
  expect_equal(state()$assemble_switch, 3L)
  expect_equal(state()$find_rrna_notes,
    "rRNA search incomplete: 1 of 3 part(s) of the supplied assembly could not be scanned")
  expect_equal(readLines("status.txt"), "empty")
  expect_false(file.exists("S1_assembly_1.fasta"))

  expect_equal(run("chunk_1.gff", "chunk_2.gff", "chunk_3.gff"), 0L)
  expect_equal(readLines("status.txt"), "ok")
  expect_equal(state()$find_rrna_notes, "3 rDNA contig(s) found in the supplied assembly")
  expect_length(readLines("rrna_scan.gff"), 3L)
})
