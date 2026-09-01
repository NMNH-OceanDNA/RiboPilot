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

chim_win <- function(starts, lineage, lineage_taxids, pident = 99, species = "sp") {
  data.frame(start = starts, stop = starts + 499L, accession = paste0("a", starts),
    species = species, staxid = 1L, lineage = lineage, lineage_taxids = lineage_taxids,
    rank_distance = NA_integer_, hit_class = "nuclear_rrna", pident = pident,
    bitscore = 500, db_source = "panel", stringsAsFactors = FALSE)
}

test_that("segments collapse adjacent windows sharing a class", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  w <- rbind(
    chim_win(seq(1, 2001, 250), "A;Chordata;Actinopterygii;o;f;g;Salmo salar", fish, species = "Salmo salar"),
    chim_win(seq(2251, 4001, 250), "F;Ascomycota;Eurotiomycetes;o;f;g;Penicillium sp.", fung, 96, "Penicillium sp."))
  s <- RiboPilot:::chimera_segments(w)
  expect_equal(nrow(s), 2)
  expect_equal(s$taxon, c("Actinopterygii", "Eurotiomycetes"))
  expect_equal(s$start, c(1L, 2251L))
  expect_equal(s$stop, c(2500L, 4500L))
  ch <- RiboPilot:::chimera_channel_windows(s)
  expect_equal(ch$score, 0.5)
  expect_equal(ch$junctions, 2375L)
  expect_match(ch$reason, "Actinopterygii")
  expect_match(ch$reason, "Penicillium")
})

test_that("a sliver segment cannot fire the window channel", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  w <- rbind(
    chim_win(seq(1, 5001, 250), "A;Chordata;Actinopterygii;o;f;g;s", fish),
    chim_win(seq(5251, 5501, 250), "F;Ascomycota;Eurotiomycetes;o;f;g;s", fung))
  ch <- RiboPilot:::chimera_channel_windows(RiboPilot:::chimera_segments(w))
  expect_equal(ch$score, 0)
})

test_that("identity-only mode reports junctions but scores zero", {
  w <- chim_win(seq(1, 3001, 250), NA_character_, NA_character_,
    pident = c(rep(99, 6), rep(80, 7)))
  s <- RiboPilot:::chimera_segments(w)
  ch <- RiboPilot:::chimera_channel_windows(s, identity_only = TRUE)
  expect_equal(ch$score, 0)
  expect_length(ch$junctions, 1)
})

test_that("chimera_window_hits with no hits returns NA rows, one per tile", {
  tiles <- data.frame(start = c(1L, 251L), stop = c(500L, 750L))
  w <- RiboPilot:::chimera_window_hits(data.frame(), tiles)
  expect_equal(nrow(w), nrow(tiles))
  expect_true(all(is.na(w$accession)))
  expect_true(all(is.na(w$pident)))
})

test_that("a taxon segment next to an NA-lineage high-identity segment cannot fire", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  w <- rbind(
    chim_win(seq(1, 2001, 250), "A;Chordata;Actinopterygii;o;f;g;s", fish),
    chim_win(seq(2251, 4001, 250), NA_character_, NA_character_, pident = 95))
  ch <- RiboPilot:::chimera_channel_windows(RiboPilot:::chimera_segments(w))
  expect_equal(ch$score, 0)
})

test_that("ribotyper short output is parsed into codes", {
  rt <- RiboPilot:::chimera_parse_ribotyper(test_path("fixtures/chimera/ribotyper.short.out"))
  expect_equal(rt$target, c("u18S", "u28S", "u5_8S"))
  expect_equal(rt$pass, c(TRUE, FALSE, FALSE))
  expect_equal(rt$codes[[2]], c("InconsistentHits", "LowScoreDifference"))
  ch <- RiboPilot:::chimera_channel_ribotyper(rt, data.frame(target = c("u18S", "u28S"), start = c(1L, 3000L)))
  expect_equal(ch$score, 0.8)
  expect_match(ch$reason, "28S.*InconsistentHits")
})

test_that("a space inside the unexpected_features column does not drop codes", {
  rt <- RiboPilot:::chimera_parse_ribotyper(test_path("fixtures/chimera/ribotyper.short.out"))
  expect_equal(rt$codes[[3]], c("LowScoreDifference", "BothStrands"))
})

test_that("itsx chimeric alone is a weak ribotyper-channel signal", {
  ch <- RiboPilot:::chimera_channel_ribotyper(NULL, NULL, itsx_chimeric = TRUE)
  expect_equal(ch$score, 0.2)
})

test_that("depth steps ignore feature edges", {
  dep <- c(rep(10, 1000), rep(50, 1000))
  feats <- data.frame(element = "18S", start = 1L, stop = 1000L)
  expect_length(RiboPilot:::chimera_depth_steps(dep, feats), 0)
  expect_equal(RiboPilot:::chimera_depth_steps(dep, data.frame(element = character(0), start = integer(0), stop = integer(0))), 1000L)
})

test_that("clip pileups need depth and fraction", {
  dep <- rep(20, 500); cl <- integer(500); cl[300] <- 8; cl[10] <- 8
  expect_equal(RiboPilot:::chimera_clip_pileups(dep, cl), 300L)
  ch <- RiboPilot:::chimera_channel_reads(dep, cl, data.frame(element = character(0), start = integer(0), stop = integer(0)))
  expect_equal(ch$score, 0.7)
  expect_equal(ch$junctions, 300L)
})

chim_ch <- function(score, junctions = integer(0), availability = "available") {
  list(score = score, junctions = junctions, reason = "r", availability = availability)
}

test_that("taxonomy-only evidence caps at suspect", {
  s <- RiboPilot:::chimera_score(list(genes = chim_ch(0.6, 2400L), windows = chim_ch(0.5, 2375L),
    ribotyper = chim_ch(0, availability = "unavailable"), reads = chim_ch(0, availability = "unavailable")))
  expect_equal(s$score, 0.7)
  expect_equal(s$label, "suspect")
})

test_that("corroborated junctions add a bonus and read evidence can reach likely", {
  s <- RiboPilot:::chimera_score(list(genes = chim_ch(0.6, 2400L), windows = chim_ch(0.5, 2375L),
    ribotyper = chim_ch(0), reads = chim_ch(0.7, 2410L)))
  expect_equal(s$score, 0.95)
  expect_equal(s$label, "likely")
  expect_equal(s$junctions, 2375L)
})

test_that("labels follow the thresholds", {
  expect_equal(RiboPilot:::chimera_label(0.29), "clean")
  expect_equal(RiboPilot:::chimera_label(0.3), "suspect")
  expect_equal(RiboPilot:::chimera_label(0.7), "suspect")
  expect_equal(RiboPilot:::chimera_label(0.71), "likely")
})

test_that("chimera_unit writes lineage-bearing rows with mocked tools", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  seq <- paste(rep("ACGT", 1000), collapse = "")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, ?)", params = list(seq))
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = '33208;7711;186623;8030;8016;8017;8018' WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1',1,1,?,?)",
    params = list(paste(rep(20, 4000), collapse = ","), paste(c(rep(0, 2199), 9, rep(0, 1800)), collapse = ",")))
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    lines <- vapply(ids, function(w) {
      st <- as.integer(sub("^w(\\d+)_.*", "\\1", w))
      if (st < 2200) sprintf("NC_1\tSalmo salar 18S\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t1800\t%s", w)
      else sprintf("KX_1\tPenicillium sp. 28S\t97.0\t100\t0\t5073\t880\t500\t1\t500\t1\t500\t4000\t%s", w)
    }, character(1))
    writeLines(lines, out); 0L
  }
  fake_rt <- function(fasta, outdir, ...) {
    writeLines(c("#idx target class strnd p/f unexpected_features", "1 u18S SSU.Eukarya plus PASS -",
      "2 u28S LSU.Eukarya plus PASS -"), file.path(outdir, "rt.ribotyper.short.out")); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = c(8030L, 5073L),
    lineage = c("Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar",
                "Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium sp."),
    lineage_taxids = c("33208;7711;186623;8030;8016;8017;8018", "4751;4890;147545;5042;5043;5044;5045"))
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, chimera_run_ribotyper = fake_rt,
    tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel", taxdump_dir = NULL)
  expect_equal(r$label, "likely")
  w <- DBI::dbGetQuery(con, "SELECT start, lineage FROM chimera_windows WHERE ID='S1' ORDER BY start")
  expect_gt(nrow(w), 10)
  expect_match(w$lineage[1], "Actinopterygii")
  seg <- jsonlite::fromJSON(DBI::dbGetQuery(con, "SELECT segments FROM chimera")$segments)
  expect_equal(seg$taxon, c("Actinopterygii", "Eurotiomycetes"))
  expect_match(DBI::dbGetQuery(con, "SELECT reason FROM chimera")$reason, "Penicillium")
})

test_that("chimera_label_prefix passes through a non-fasta label_db unchanged", {
  expect_equal(RiboPilot:::chimera_label_prefix("labels", tempdir()), "labels")
  expect_null(RiboPilot:::chimera_label_prefix(NULL, tempdir()))
})

test_that("chimera_label_prefix builds a blast db for a labels fasta", {
  skip_if_not_installed("withr")
  skip_if(!nzchar(Sys.which("makeblastdb")))
  d <- withr::local_tempdir()
  fa <- file.path(d, "labels.fasta")
  writeLines(c(">L1", "ACGTACGTACGTACGTACGT"), fa)
  prefix <- RiboPilot:::chimera_label_prefix(fa, d)
  expect_true(file.exists(paste0(prefix, ".nin")))
})

test_that("chimera_unit falls back to label_db when dbarg is NULL, scoring identity-only", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  seq <- paste(rep("ACGT", 1000), collapse = "")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, ?)", params = list(seq))
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = '33208;7711;186623;8030;8016;8017;8018' WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1',1,1,?,?)",
    params = list(paste(rep(20, 4000), collapse = ","), paste(c(rep(0, 2199), 9, rep(0, 1800)), collapse = ",")))
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    lines <- vapply(ids, function(w) {
      st <- as.integer(sub("^w(\\d+)_.*", "\\1", w))
      if (st < 2200) sprintf("NC_1\tSalmo salar 18S\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t1800\t%s", w)
      else sprintf("KX_1\tPenicillium sp. 28S\t97.0\t100\t0\t5073\t880\t500\t1\t500\t1\t500\t4000\t%s", w)
    }, character(1))
    writeLines(lines, out); 0L
  }
  fake_rt <- function(fasta, outdir, ...) {
    writeLines(c("#idx target class strnd p/f unexpected_features", "1 u18S SSU.Eukarya plus PASS -",
      "2 u28S LSU.Eukarya plus PASS -"), file.path(outdir, "rt.ribotyper.short.out")); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = c(8030L, 5073L),
    lineage = c("Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar",
                "Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium sp."),
    lineage_taxids = c("33208;7711;186623;8030;8016;8017;8018", "4751;4890;147545;5042;5043;5044;5045"))
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, chimera_run_ribotyper = fake_rt,
    tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = NULL, taxdump_dir = NULL, label_db = "labels")
  expect_equal(r$db_source, "labels")
  ch <- jsonlite::fromJSON(r$channels)
  expect_equal(ch$windows$availability, "identity_only")
  expect_equal(ch$windows$score, 0)
})
