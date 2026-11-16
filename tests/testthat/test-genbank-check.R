# GenBank BLAST check: the pure flag function on the SUB16547976 fixture
# (NCBI's own flags from its 2026-10-15 reply), the grammar, and the unit driver.

gc_fix <- function(f) test_path("fixtures", "genbank_check", f)
gc_hits <- function(f) RiboPilot:::parse_blast_hits_file(gc_fix(f))
gc_split <- function(h) split(h, h$qseqid)
gc_seed <- function(con, seq = "ACGT") {
  DBI::dbExecute(con, "INSERT OR REPLACE INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,?,?)",
    params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT OR REPLACE INTO annotate (ID, path, scaffold) VALUES ('S1',1,1)")
}

test_that("the SUB16547976 fixture reproduces NCBI's two flag lists exactly", {
  its <- gc_split(gc_hits("SUB16547976_ITS_eukaryote_sequences.tsv"))
  rr <- gc_split(gc_hits("SUB16547976_SSU_LSU_eukaryote_rRNA.tsv"))
  want <- read.delim(gc_fix("SUB16547976_ncbi_flags.tsv"), stringsAsFactors = FALSE)
  got <- vapply(want$seqid, function(s) RiboPilot:::genbank_check_calls(its[[s]], rr[[s]]), character(1))
  expect_equal(unname(grepl("genbank_low_coverage(ITS", got, fixed = TRUE)), want$low_coverage == 1)
  expect_equal(unname(grepl("genbank_multi_hsp(ITS", got, fixed = TRUE)), want$hsp_problem == 1)
  expect_false(any(grepl("rRNA", got[want$low_coverage == 0 & want$hsp_problem == 0])))
  expect_equal(unname(got[c("07-01-61_S93", "07-02-111_S335")]), c("", ""))
  expect_match(got[["07-01-65_S97"]],
    "^genbank_low_coverage\\(ITS 70%\\); genbank_multi_hsp\\(ITS 2 HSPs: 1-1855, 2274-2457\\) top hits: ITS OP250806.1 Datnioides microlepis 70%, rRNA MG895428.1 Brachypleura novaezeelandiae 63%$")
})

test_that("genbank_check_calls covers the grammar", {
  mk <- function(acc, spans, qcovs, title = "Danio rerio SSU ITS1") data.frame(accession = acc, title = title,
    pident = 99, qcovs = qcovs, evalue = 0, staxid = 1L, bitscore = rev(seq_along(spans)) * 100,
    aln_len = 10L, qstart = spans, qend = spans + 9L, sstart = 1L, send = 10L, slen = 100L, qseqid = "q",
    stringsAsFactors = FALSE)
  none <- data.frame()
  expect_equal(RiboPilot:::genbank_check_calls(mk("A1", 1L, 100), mk("B1", 1L, 60)), "")
  expect_equal(RiboPilot:::genbank_check_calls(none, none), "genbank_no_hit(ITS); genbank_no_hit(rRNA)")
  expect_match(RiboPilot:::genbank_check_calls(none, mk("B1", 1L, 60)),
    "^genbank_no_hit\\(ITS\\) top hits: ITS none, rRNA B1 Danio rerio 60%$")
  # spans in query order even when HSP rows arrive in bitscore order
  expect_match(RiboPilot:::genbank_check_calls(mk("A1", c(2000L, 1L), 80), mk("B1", 1L, 60)),
    "genbank_multi_hsp\\(ITS 2 HSPs: 1-10, 2000-2009\\)")
  expect_match(RiboPilot:::genbank_check_calls(mk("A1", 1L, 100), mk("B1", c(1L, 50L), 60)),
    "^genbank_multi_hsp\\(rRNA 2 HSPs: 1-10, 50-59\\) top hits")
  expect_equal(RiboPilot:::genbank_check_display(c("", "not_run(x)",
    "genbank_low_coverage(ITS 70%); genbank_multi_hsp(ITS 2 HSPs: 1-2, 3-4); genbank_no_hit(rRNA) top hits: ITS A 1 2 70%, rRNA none")),
    c("", "not_run(x)", "low coverage ITS 70%, ITS 2 HSPs, no rRNA hit"))
  expect_equal(RiboPilot:::genbank_check_db_dir("/x/2026/ribopilot_rrna /x/2026/16S_ribosomal_RNA"), "/x/2026")
  expect_null(RiboPilot:::genbank_check_db_dir(""))
  expect_null(RiboPilot:::genbank_check_db_dir(NULL))
})

test_that("genbank_check_unit stores text, not_run reasons, and the review copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  gc_seed(con)
  d <- withr::local_tempdir()
  for (n in RiboPilot:::GENBANK_CHECK_DBS) file.create(file.path(d, paste0(n, ".nin")))
  gc <- function() DBI::dbGetQuery(con, "SELECT genbank_check FROM annotate WHERE ID = 'S1'")$genbank_check
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = NULL)
  expect_equal(gc(), "not_run(no local BLAST)")
  unlink(file.path(d, "ITS_eukaryote_sequences.nin"))
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(no ITS_eukaryote_sequences database)")
  file.create(file.path(d, "ITS_eukaryote_sequences.nin"))
  testthat::local_mocked_bindings(genbank_blast_run = function(...) stop("boom"), .package = "RiboPilot")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(BLAST failed: boom)")
  its <- gc_hits("SUB16547976_ITS_eukaryote_sequences.tsv"); its <- its[its$qseqid == "07-01-65_S97", ]
  rr <- gc_hits("SUB16547976_SSU_LSU_eukaryote_rRNA.tsv"); rr <- rr[rr$qseqid == "07-01-65_S97", ]
  testthat::local_mocked_bindings(genbank_blast_run = function(query, dbarg, ...)
    if (grepl("ITS_eukaryote", dbarg)) its else rr, .package = "RiboPilot")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_match(gc(), "^genbank_low_coverage\\(ITS 70%\\)")
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = NULL")
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d)
  expect_equal(gc(), "not_run(no sequence)")
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = 'ACGT'")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  RiboPilot:::genbank_check_unit(con, "S1", 1L, 1L, db_dir = d, review = TRUE)
  expect_match(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units")$genbank_check, "^genbank_low_coverage")
  expect_equal(gc(), "not_run(no sequence)")
})

test_that("genbank_check_unit reproduces NCBI on SUB16547976 with real blastn", {
  d <- Sys.getenv("RIBOPILOT_GENBANK_CHECK_DBDIR", "")
  skip_if(!nzchar(d), "RIBOPILOT_GENBANK_CHECK_DBDIR not set")
  skip_if(!nzchar(Sys.which("blastn")), "blastn not on PATH")
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  fa <- readLines(gc_fix("SUB16547976.fasta"))
  ids <- sub("^>(\\S+).*", "\\1", fa[startsWith(fa, ">")]); seqs <- fa[!startsWith(fa, ">")]
  want <- read.delim(gc_fix("SUB16547976_ncbi_flags.tsv"), stringsAsFactors = FALSE)
  for (i in seq_along(ids)) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES (?,1,1,?,?)",
      params = list(ids[i], seqs[i], nchar(seqs[i])))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES (?,1,1)", params = list(ids[i]))
    RiboPilot:::genbank_check_unit(con, ids[i], 1L, 1L, db_dir = d)
  }
  got <- DBI::dbGetQuery(con, "SELECT ID, genbank_check FROM annotate")
  got <- got$genbank_check[match(want$seqid, got$ID)]
  info <- system2("blastdbcmd", c("-db", shQuote(file.path(d, "ITS_eukaryote_sequences")), "-info"), stdout = TRUE)
  message("ITS db build: ", grep("^Date:", info, value = TRUE)[1])
  expect_equal(grepl("genbank_low_coverage(ITS", got, fixed = TRUE), want$low_coverage == 1)
  expect_equal(grepl("genbank_multi_hsp(ITS", got, fixed = TRUE), want$hsp_problem == 1)
})
