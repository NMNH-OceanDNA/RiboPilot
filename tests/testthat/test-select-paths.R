# Wiring real assembly output into the scorer: per-path length, graph depth,
# anchors, and the top BLAST hit. External tools (barrnap, blastn) are mocked so
# these run anywhere; what is pinned is the assembly of the candidate table.

mk_path_fasta <- function(dir, i, headers, seqs) {
  f <- file.path(dir, sprintf("path%d.fasta", i))
  writeLines(as.vector(rbind(paste0(">", headers), seqs)), f)
  f
}

test_that("sseqid_species drops the accession and unescapes underscores", {
  f <- RiboPilot:::sseqid_species
  expect_equal(f("MT12345_Sebastes_alutus"), "Sebastes alutus")
  expect_equal(f("XR999_Homo_sapiens"), "Homo sapiens")
})

test_that("parse_blast_best keeps the top bitscore per query", {
  hits <- RiboPilot:::parse_blast_best(c(
    "q1\tACC_A_a\t90.0\t100",
    "q1\tACC_B_b\t99.0\t500",   # better bitscore, lower down the output
    "q2\tACC_C_c\t80.0\t200"))
  expect_equal(nrow(hits), 2)
  expect_equal(hits$sseqid[hits$qseqid == "q1"], "ACC_B_b")
  expect_equal(hits$bitscore[hits$qseqid == "q1"], 500)
  expect_true(all(is.na(hits$staxid)))
})

test_that("parse_blast_best takes the first taxid when staxids are requested", {
  hits <- RiboPilot:::parse_blast_best("q1\tACC\t99.0\t500\t7955;9606", with_staxids = TRUE)
  expect_equal(hits$staxid, "7955")
})

test_that("parse_blast_best returns an empty frame for no hits", {
  expect_equal(nrow(RiboPilot:::parse_blast_best(character(0))), 0)
})

test_that("score_assembly derives length and depth per path", {
  dir <- withr::local_tempdir()
  # two paths, one clearly longer, both traversing known graph segments
  f1 <- mk_path_fasta(dir, 1, "scaffold_1--11+", strrep("ACGT", 500))    # 2000 bp
  f2 <- mk_path_fasta(dir, 2, "scaffold_1--22+", strrep("ACGT", 2000))   # 8000 bp
  gfa <- file.path(dir, "assembly.gfa")
  writeLines(c("S\t11\tA\tLN:i:100\tRC:i:1000",     # depth 10
               "S\t22\tA\tLN:i:100\tRC:i:10000"), gfa)  # depth 100

  local_mocked_bindings(assembly_anchors = function(...) c("18S", "28S"),
    .package = "RiboPilot")
  out <- score_assembly(c(f1, f2), gfa_path = gfa, sample_lineage = "A; B")

  expect_equal(sort(out$length), c(2000, 8000))
  expect_equal(out$depth[out$path == 1], 10)
  expect_equal(out$depth[out$path == 2], 100)
  expect_equal(out$anchors[out$path == 1], "18S,28S")
  # path 2 wins on both depth and length plausibility
  expect_equal(out$path[1], 2L)
  expect_true(out$auto_select[1])
})

test_that("score_assembly works with no graph and no BLAST db", {
  dir <- withr::local_tempdir()
  f1 <- mk_path_fasta(dir, 1, "scaffold_1", strrep("ACGT", 1500))
  local_mocked_bindings(assembly_anchors = function(...) character(0),
    .package = "RiboPilot")
  out <- score_assembly(f1)
  expect_equal(nrow(out), 1)
  expect_true(is.na(out$depth))
  expect_true(is.na(out$blast_pident))
  expect_true(is.na(out$blast_lineage))
  # missing evidence must not produce an NA score
  expect_false(is.na(out$score))
  expect_true(out$auto_select)
})

test_that("score_assembly carries the top hit's identity into the score", {
  dir <- withr::local_tempdir()
  f1 <- mk_path_fasta(dir, 1, "s1", strrep("ACGT", 2000))
  f2 <- mk_path_fasta(dir, 2, "s2", strrep("ACGT", 2000))
  local_mocked_bindings(
    .package = "RiboPilot",
    assembly_anchors = function(...) "18S",
    blast_best_hits = function(query_fasta, db, with_staxids = FALSE) {
      pid <- if (grepl("path1", query_fasta)) 70 else 99
      data.frame(qseqid = "q", sseqid = "ACC_Genus_species", pident = pid,
        bitscore = 100, staxid = NA_character_, stringsAsFactors = FALSE)
    })
  out <- score_assembly(c(f1, f2), blast_db = "dummy")
  expect_equal(out$blast_pident[out$path == 1], 70)
  expect_equal(out$blast_pident[out$path == 2], 99)
  expect_equal(out$hit_species[out$path == 1], "Genus species")
  # identical on every other term, so the better identity must rank first
  expect_equal(out$path[1], 2L)
})
