# WF1's per-unit R steps, run through their real bin scripts in pipeline order
# on a project made by new_db().

run_bin <- function(script, args, tools) {
  processx::run("Rscript",
    c(file.path(system.file("nextflow", "bin", package = "RiboPilot"), script), args),
    env = c("current", RIBOPILOT_SRC = RiboPilot:::ribopilot_src_default(),
      PATH = paste(tools$bin, Sys.getenv("PATH"), sep = .Platform$path.sep)),
    error_on_status = FALSE)
}

test_that("score, reconcile, curate, and triage run through their bin scripts", {
  skip_on_cran()
  skip_on_os("windows")
  tools <- local_fake_pipeline_tools()
  d <- withr::local_tempdir()
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  db <- file.path(d, ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = tools$scaffold, labels_db = tools$scaffold)
  seq <- readLines(tools$scaffold)[2]
  writeLines(c(">S1.1.1 linear fake", seq), file.path(d, "S1_assembly_1.fasta"))
  unit <- file.path(d, "unit.fasta")
  writeLines(c(">S1.1.1", seq), unit)
  gff <- file.path(d, "unit.gff")
  system2(file.path(tools$bin, "pybarrnap"), c("-q", unit), stdout = gff)
  system2(file.path(tools$bin, "ITSx"), c("-i", unit, "-o", file.path(d, "unit.itsx")))
  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  q <- function(sql) DBI::dbGetQuery(con, sql)

  res <- run_bin("rp_score.R", c(db, "S1", d), tools)
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(q("SELECT assemble_switch FROM assemble WHERE ID = 'S1'")$assemble_switch, 4L)
  expect_equal(nchar(q("SELECT sequence FROM assemblies WHERE ID = 'S1'")$sequence), 6310L)
  expect_equal(q("SELECT auto_select FROM path_scores WHERE ID = 'S1'")$auto_select, 1L)

  res <- run_bin("rp_reconcile.R",
    c(db, "S1", "1", "1", gff, file.path(d, "unit.itsx.positions.txt"), file.path(d, "unit.itsx")), tools)
  expect_equal(res$status, 0L, info = res$stderr)
  a <- q("SELECT annotate_switch, n_features, n_flags FROM annotate WHERE ID = 'S1'")
  expect_equal(a$annotate_switch, 2L)
  expect_equal(a$n_features, 5L)
  expect_equal(sort(q("SELECT element FROM annotations WHERE ID = 'S1'")$element),
    c("18S", "28S", "5.8S", "ITS1", "ITS2"))
  expect_true(is.na(a$n_flags))

  res <- run_bin("rp_curate.R", c(db, "S1", "1", "1"), tools)
  expect_equal(res$status, 0L, info = res$stderr)
  expect_false(is.na(q("SELECT n_flags FROM annotate WHERE ID = 'S1'")$n_flags))

  res <- run_bin("rp_triage.R", c(db, "S1"), tools)
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(q("SELECT assemble_switch FROM assemble WHERE ID = 'S1'")$assemble_switch, 2L)
})

test_that("blast_scaffolds_local names each scaffold's best seed hit", {
  d <- withr::local_tempdir()
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  db <- file.path(d, ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE, seeds_db = fa, labels_db = fa)
  writeLines(c(">S1.1.1 linear a", "ACGTACGTAC", ">S1.1.2 linear b", "TTGGCCAATT"),
    file.path(d, "S1_assembly_1.fasta"))
  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con))
  RiboPilot:::load_assemblies(con, d, "S1")
  local_mocked_bindings(blast_best_hits = function(...) data.frame(qseqid = "1.2",
    sseqid = "AB123_Danio_rerio", pident = 99.44), .package = "RiboPilot")
  RiboPilot:::blast_scaffolds_local(con, "S1", "seeds_db")
  got <- DBI::dbGetQuery(con, "SELECT scaffold, blast_accession, blast_species, blast_pident,
    blast_source FROM assemblies WHERE ID = 'S1' ORDER BY scaffold")
  expect_equal(got$blast_accession, c(NA, "AB123"))
  expect_equal(got$blast_species, c(NA, "Danio rerio"))
  expect_equal(got$blast_pident, c(NA, 99.4))
  expect_equal(got$blast_source, c(NA, "seeds"))
})
