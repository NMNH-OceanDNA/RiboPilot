# Retained BLAST hit distribution: parsing, bitscore ranking, and the two
# tables derived from it. No network, no taxonkit (lineages resolve to NA).

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"

# saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send
hit <- function(acc, title, pident, qcovs, bits, len, qs, qe, staxid = 7955) {
  paste(acc, title, pident, qcovs, "1e-90", staxid, bits, len, qs, qe, 1, len, sep = "\t")
}

mk_blast_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, blast_accession TEXT,
    blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,?)", params = list(FISH_T))
  con
}

write_hits <- function(lines) {
  f <- withr::local_tempfile(.local_envir = parent.frame())
  writeLines(lines, f)
  f
}

test_that("every HSP is retained, not one row per subject", {
  f <- write_hits(c(
    hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99.1, 80, 1500, 1600, 1, 1600),
    hit("AA1", "Danio rerio 18S ribosomal RNA gene", 97.0, 80, 400, 420, 3000, 3420),
    hit("AA2", "Danio rerio 28S ribosomal RNA gene", 98.0, 60, 900, 950, 1, 950)))
  d <- RiboPilot:::parse_blast_hits_file(f)
  expect_equal(nrow(d), 3)
  expect_equal(sum(d$accession == "AA1"), 2)
  expect_equal(d$qstart[2], 3000)
  expect_equal(d$aln_len[1], 1600)
})

test_that("subjects rank by bitscore, not pident * qcovs", {
  # The trap: a short conserved-core hit at 99% to an unrelated organism vs a
  # long hit at 95% to a genuine relative. pident*qcovs prefers the short one.
  f <- write_hits(c(
    hit("SHORT", "Gallus gallus 18S ribosomal RNA gene", 99.5, 95, 500, 400, 1, 400, staxid = 9031),
    hit("LONG",  "Danio rerio 18S ribosomal RNA gene",   95.0, 95, 4000, 5000, 1, 5000)))
  d <- RiboPilot:::parse_blast_hits_file(f)
  expect_equal(d$rank_score, c(500, 4000))          # bitscore, not the product

  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  parse_blast_candidates(con, "S1", 1, 1, f)
  h <- DBI::dbGetQuery(con, "SELECT accession, hit_rank FROM blast_hits ORDER BY hit_rank")
  expect_equal(h$accession[1], "LONG")              # the long genuine hit wins
  scr <- DBI::dbGetQuery(con, "SELECT accession FROM blast_screen")
  expect_equal(scr$accession, "LONG")
})

test_that("the screen takes the top subject unfiltered; the panel is rDNA only", {
  f <- write_hits(c(
    hit("CONTAM", "Saprolegnia parasitica strain X genome assembly", 99.0, 90, 5000, 6000, 1, 6000, staxid = 4769),
    hit("REF",    "Danio rerio 18S ribosomal RNA gene, partial",     92.0, 70, 1200, 1400, 1, 1400)))
  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  parse_blast_candidates(con, "S1", 1, 1, f)

  # screen sees the contaminant: this is the only signal that can say "not the host"
  scr <- DBI::dbGetQuery(con, "SELECT accession, is_nuclear_rrna FROM blast_screen")
  expect_equal(scr$accession, "CONTAM")
  expect_equal(scr$is_nuclear_rrna, 0L)

  # reference panel keeps only the nuclear-rDNA-titled record
  ref <- DBI::dbGetQuery(con, "SELECT accession, rank FROM blast_ref_candidates ORDER BY rank")
  expect_equal(ref$accession, "REF")

  # assemblies.blast_* comes from the panel, and is now attributed
  a <- DBI::dbGetQuery(con, "SELECT blast_accession, blast_source FROM assemblies")
  expect_equal(a$blast_accession, "REF")
  expect_equal(a$blast_source, "core_nt")
})

test_that("an empty or sentinel hits file is handled without error", {
  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  expect_silent(parse_blast_candidates(con, "S1", 1, 1, write_hits("NO_SIGNIFICANT_HITS")))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_hits")$n, 0)
  expect_equal(nrow(RiboPilot:::parse_blast_hits_file(write_hits(character(0)))), 0)
  expect_equal(nrow(RiboPilot:::parse_blast_hits_file("/nonexistent/hits.txt")), 0)
})

test_that("a pre-v9 six-column hits file still parses", {
  old <- write_hits(c(
    paste("AA1", "Danio rerio 18S ribosomal RNA gene", 99.0, 80, "1e-90", 7955, sep = "\t"),
    paste("AA2", "Danio rerio 28S ribosomal RNA gene", 95.0, 90, "1e-80", 7955, sep = "\t")))
  d <- RiboPilot:::parse_blast_hits_file(old)
  expect_equal(nrow(d), 2)
  expect_true(all(is.na(d$bitscore)))
  expect_equal(d$rank_score, c(99 * 80, 95 * 90))   # falls back to the old score

  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  expect_silent(parse_blast_candidates(con, "S1", 1, 1, old))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_hits")$n, 2)
})

test_that("re-parsing a unit replaces its rows rather than duplicating them", {
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99, 80, 1500, 1600, 1, 1600))
  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  parse_blast_candidates(con, "S1", 1, 1, f)
  parse_blast_candidates(con, "S1", 1, 1, f)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_hits")$n, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_candidates")$n, 1)
})

test_that("rank_distance is left NA when no taxdump resolved the hits", {
  # Offline projects must record "not judged", never a spurious distance.
  f <- write_hits(hit("AA1", "Danio rerio 18S ribosomal RNA gene", 99, 80, 1500, 1600, 1, 1600))
  con <- mk_blast_db(); on.exit(DBI::dbDisconnect(con))
  parse_blast_candidates(con, "S1", 1, 1, f, taxdump_dir = NULL)
  h <- DBI::dbGetQuery(con, "SELECT rank_distance, lineage_taxids FROM blast_hits")
  expect_true(is.na(h$rank_distance))
  expect_true(is.na(h$lineage_taxids))
})
