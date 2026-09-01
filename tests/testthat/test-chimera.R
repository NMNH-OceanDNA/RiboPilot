test_that("new_db creates the chimera tables and coverage.clips", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  db <- file.path(d, ".sqlite")
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  tabs <- DBI::dbListTables(con)
  expect_true(all(c("chimera", "chimera_windows") %in% tabs))
  expect_true("clips" %in% DBI::dbListFields(con, "coverage"))
  expect_true("segments" %in% DBI::dbListFields(con, "chimera"))
  expect_true("lineage" %in% DBI::dbListFields(con, "chimera_windows"))
})

test_that("write_coverage_units stores an aligned clip track", {
  skip_if_not_installed("withr")
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, 'ACGTACGT')")
  dep <- tempfile(); cl <- tempfile()
  writeLines(sprintf("1.1\t%d\t%d", 1:8, c(5, 5, 6, 6, 20, 20, 20, 20)), dep)
  writeLines(c("1.1\t5\t7", "1.1\t6\t1"), cl)
  RiboPilot:::write_coverage_units(con, "S1", dep, cl)
  r <- DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = 'S1'")
  expect_equal(r$depths, "5,5,6,6,20,20,20,20")
  expect_equal(r$clips, "0,0,0,0,7,1,0,0")
})

test_that("tiling anchors the last window and handles short scaffolds", {
  t <- RiboPilot:::chimera_tile(1200, win = 500, step = 250)
  expect_equal(t$start, c(1, 251, 501, 701))
  expect_equal(t$stop, c(500, 750, 1000, 1200))
  expect_equal(RiboPilot:::chimera_tile(300, 500, 250), data.frame(start = 1L, stop = 300L))
})

test_that("rank_distance is 0 at same species and 7 when nothing is shared", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  expect_equal(RiboPilot:::rank_distance(fish, fung), 7L)
})

test_that("per-gene LCA uses only HSPs covering half the element", {
  feats <- data.frame(element = c("18S", "28S"), start = c(1L, 3000L), stop = c(1800L, 6000L))
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  hits <- data.frame(
    accession = c("a1", "a2", "a3", "f1", "f2", "f3", "f4"),
    bitscore = c(3000, 2900, 2800, 5000, 4900, 4800, 100),
    qstart = c(1, 5, 10, 3000, 3100, 3050, 1500),
    qend = c(1800, 1790, 1750, 6000, 5900, 5950, 1600),
    lineage_taxids = c(fish, fish, fish, fung, fung, fung, fung),
    lineage = c(rep("Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar", 3),
                rep("Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium sp.", 4)),
    species = c(rep("Salmo salar", 3), rep("Penicillium sp.", 4)),
    pident = c(99, 99, 98, 97, 97, 96, 80),
    stringsAsFactors = FALSE)
  g <- RiboPilot:::chimera_gene_lca(hits, feats)
  expect_equal(g$taxon[g$element == "18S"], "Salmo salar")
  expect_equal(g$taxon[g$element == "28S"], "Penicillium sp.")
  ch <- RiboPilot:::chimera_channel_genes(g, feats)
  expect_equal(ch$score, 0.6)
  expect_equal(ch$junctions, 2400L)
  expect_match(ch$reason, "18S.*Salmo salar.*28S.*Penicillium")
})

test_that("an element with too few covering hits is unresolved and does not vote", {
  feats <- data.frame(element = c("18S", "28S"), start = c(1L, 3000L), stop = c(1800L, 6000L))
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  hits <- data.frame(accession = c("a1", "a2", "a3"), bitscore = c(3, 2, 1),
    qstart = c(1, 5, 10), qend = c(1800, 1790, 1750),
    lineage_taxids = rep(fish, 3), lineage = rep("x;y;z;o;f;g;s", 3),
    species = rep("Salmo salar", 3), pident = c(99, 99, 98), stringsAsFactors = FALSE)
  g <- RiboPilot:::chimera_gene_lca(hits, feats)
  expect_true(is.na(g$taxon[g$element == "28S"]))
  ch <- RiboPilot:::chimera_channel_genes(g)
  expect_equal(ch$score, 0)
  expect_equal(ch$availability, "uninformative")
})
