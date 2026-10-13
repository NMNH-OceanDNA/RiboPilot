test_that("chimera_unit is exported", {
  expect_true("chimera_unit" %in% getNamespaceExports("RiboPilot"))
})

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

test_that("a dropped scaffold keeps its coverage track through a coverage re-run", {
  skip_if_not_installed("withr")
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, ignore) VALUES ('S1', 1, 1, 'ACGT', 0)")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, ignore) VALUES ('S1', 1, 2, 'ACGT', 1)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths) VALUES ('S1', 1, 1, 'old'), ('S1', 1, 2, '9,9')")
  dep <- tempfile()
  writeLines(sprintf("1.1\t%d\t5", 1:4), dep)
  RiboPilot:::write_coverage_units(con, "S1", dep)
  r <- DBI::dbGetQuery(con, "SELECT scaffold, depths FROM coverage WHERE ID = 'S1' ORDER BY scaffold")
  expect_equal(r$depths, c("5,5,5,5", "9,9"))
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
  expect_match(ch$reason, "18S.*Salmo salar.*28S.*Penicillium.*different kingdom")
  expect_true(ch$severe)
})

test_that("the gene channel treats an ancestor as agreement and grades real switches", {
  gt <- function(lt2, pid2 = 99) data.frame(element = c("18S", "5.8S"), taxon = c("A", "B"),
    rank = NA_character_, taxid = NA_integer_, species = "sp", pident = c(99, pid2), n_hits = 3L,
    lineage_taxids = c("33208;7711;186623;;;;", lt2), stringsAsFactors = FALSE)
  anc <- data.frame(element = c("18S", "5.8S"), taxon = c("Chordata", "Pleuronectiformes"),
    rank = NA_character_, taxid = NA_integer_, species = "sp", pident = 99, n_hits = 3L,
    lineage_taxids = c("33208;7711;;;;;", "33208;7711;186623;8256;;;"), stringsAsFactors = FALSE)
  anc_ch <- RiboPilot:::chimera_channel_genes(anc)
  expect_equal(anc_ch$score, 0.3)
  expect_false(anc_ch$severe)
  anc_deep <- data.frame(element = c("18S", "5.8S"), taxon = c("Actinopterygii", "Pleuronectiformes"),
    rank = NA_character_, taxid = NA_integer_, species = "sp", pident = 99, n_hits = 3L,
    lineage_taxids = c("33208;7711;186623;;;;", "33208;7711;186623;8256;;;"), stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::chimera_channel_genes(anc_deep)$score, 0)
  cls <- RiboPilot:::chimera_channel_genes(gt("33208;7711;40674;;;;"))
  expect_equal(cls$score, 0.6)
  expect_false(cls$severe)
  expect_true(RiboPilot:::chimera_channel_genes(gt("4751;4890;147545;;;;"))$severe)
  expect_false(RiboPilot:::chimera_channel_genes(gt("4751;4890;147545;;;;", pid2 = 85))$severe)
})

test_that("a severe switch is found on any pair, not only the first one checked", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  g <- data.frame(element = c("18S", "ITS1", "5.8S", "ITS2", "28S"), taxon = "x",
    rank = NA_character_, taxid = NA_integer_, species = "sp", pident = c(98, 88, 95, 96, 97),
    n_hits = 3L, lineage_taxids = c(fish, fung, fung, fung, fung), stringsAsFactors = FALSE)
  expect_true(RiboPilot:::chimera_channel_genes(g)$severe)
  s <- data.frame(start = c(1L, 1001L, 2001L), stop = c(1000L, 2000L, 3000L),
    taxon = c("Actinopterygii", "Eurotiomycetes", "Chondrichthyes"), rank = "class",
    species = "sp", pident = c(80, 97, 98), n_windows = 4L, db_source = "panel",
    lineage_taxids = c(fish, fung, "33208;7711;7777;;;;"), stringsAsFactors = FALSE)
  expect_true(RiboPilot:::chimera_channel_windows(s)$severe)
})

test_that("genes that agree at the same shallow rank are not a chimera", {
  g <- data.frame(element = c("18S", "28S"), taxon = "Mollusca", rank = NA_character_,
    taxid = NA_integer_, species = "sp", pident = 99, n_hits = 3L,
    lineage_taxids = "33208;6447;;;;;", stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::chimera_channel_genes(g)$score, 0)
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
  expect_true(ch$severe)
})

test_that("window hits tied at the top bitscore resolve to the rank they share", {
  hits <- data.frame(qseqid = c("w1_500", "w1_500", "w251_750"), accession = c("XR_1", "OP_1", "OP_2"),
    title = c("Canis aureus 28S ribosomal RNA", "Paralichthys sp. 28S ribosomal RNA", "Paralichthys sp. 28S ribosomal RNA"),
    staxid = c(1L, 2L, 2L), pident = c(100, 99.8, 99), bitscore = c(924, 918, 900), stringsAsFactors = FALSE)
  fake_lin <- function(taxids, dir) data.frame(taxid = c(1L, 2L),
    lineage = c("Metazoa;Chordata;Mammalia;Carnivora;Canidae;Canis;Canis aureus",
                "Metazoa;Chordata;Actinopteri;Pleuronectiformes;Paralichthyidae;Paralichthys;Paralichthys sp."),
    lineage_taxids = c("33208;7711;40674;33554;9608;9611;44565", "33208;7711;186623;8256;8257;8258;8259"),
    stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  w <- RiboPilot:::chimera_window_hits(hits, data.frame(start = c(1L, 251L), stop = c(500L, 750L)))
  expect_equal(RiboPilot:::chimera_rank_name(w$lineage, "class"), c(NA, "Actinopteri"))
  expect_equal(w$species, c("Chordata (2 tied hits)", "Paralichthys sp."))
  expect_equal(RiboPilot:::chimera_segments(w)$taxon, c("unresolved", "Actinopteri"))
  fish <- "33208;7711;186623;8256;8257;8258;8259"
  s <- RiboPilot:::chimera_segments(rbind(
    chim_win(seq(1, 751, 250), "Metazoa;Chordata;Actinopteri;o;f;g;s", fish),
    chim_win(seq(1001, 1501, 250), "Metazoa;Chordata;;;;;", "33208;7711;;;;;"),
    chim_win(seq(1751, 2501, 250), "Metazoa;Chordata;Actinopteri;o;f;g;s", fish)))
  expect_equal(s$taxon, c("Actinopteri", "unresolved", "Actinopteri"))
  expect_equal(RiboPilot:::chimera_channel_windows(s)$score, 0)
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

test_that("a class-level switch needs four windows a side; a phylum switch needs two", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  frog <- "33208;7711;8292;8342;8343;8344;8345"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  seg <- function(lt, taxon, nw) data.frame(start = c(1L, 401L), stop = c(500L, 900L),
    taxon = c("Actinopteri", taxon), rank = "class", species = "sp", pident = 100,
    n_windows = nw, db_source = "panel", lineage_taxids = c(fish, lt), stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::chimera_channel_windows(seg(frog, "Amphibia", c(4L, 2L)))$score, 0)
  ch <- RiboPilot:::chimera_channel_windows(seg(frog, "Amphibia", c(4L, 4L)))
  expect_equal(ch$score, 0.5)
  expect_false(ch$severe)
  ch <- RiboPilot:::chimera_channel_windows(seg(fung, "Eurotiomycetes", c(4L, 2L)))
  expect_equal(ch$score, 0.5)
  expect_true(ch$severe)
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

test_that("clip pileups need depth and fraction and skip the scaffold ends", {
  dep <- rep(20, 1000); cl <- integer(1000); cl[500] <- 9; cl[100] <- 9; cl[700] <- 7
  expect_equal(RiboPilot:::chimera_clip_pileups(dep, cl), 500L)
  ch <- RiboPilot:::chimera_channel_reads(dep, cl, data.frame(element = character(0), start = integer(0), stop = integer(0)))
  expect_equal(ch$score, 0.7)
  expect_equal(ch$junctions, 500L)
})

test_that("depth steps skip repeats, scaffold ends, and thin coverage", {
  set.seed(1)
  rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
  none <- data.frame(element = character(0), start = integer(0), stop = integer(0))
  dip <- rep(60, 4600); dip[2001:2600] <- 5
  expect_length(RiboPilot:::chimera_depth_steps(dip, none, paste0(rnd(2000), strrep(rnd(100), 6), rnd(2000))), 0)
  expect_true(length(RiboPilot:::chimera_depth_steps(dip, none, rnd(4600))) > 0)
  expect_length(RiboPilot:::chimera_depth_steps(c(rep(60, 4200), rep(5, 400)), none, rnd(4600)), 0)
  expect_length(RiboPilot:::chimera_depth_steps(c(rep(15, 2300), rep(3, 2300)), none, rnd(4600)), 0)
  expect_equal(RiboPilot:::chimera_depth_steps(c(rep(60, 2300), rep(10, 2300)), none, rnd(4600)), 2250L)
})


chim_ch <- function(score, junctions = integer(0), availability = "available") {
  list(score = score, junctions = junctions, reason = "r", availability = availability)
}

test_that("only a severe channel makes a severe label; other evidence is minor", {
  sev <- RiboPilot:::chimera_score(list(genes = chim_ch(0.6, 2400L), windows = c(chim_ch(0.5, 2375L), severe = TRUE),
    ribotyper = chim_ch(0), reads = chim_ch(0.7, 2410L)))
  expect_equal(sev$label, "severe")
  expect_equal(sev$junctions, 2375L)
  reads <- RiboPilot:::chimera_score(list(genes = chim_ch(0), windows = chim_ch(0),
    ribotyper = chim_ch(0), reads = chim_ch(0.4, 3000L)))
  expect_equal(reads$label, "minor")
  itsx <- RiboPilot:::chimera_score(list(genes = chim_ch(0), windows = chim_ch(0),
    ribotyper = chim_ch(0.2), reads = chim_ch(0)))
  expect_equal(itsx$label, "clean")
})

test_that("labels follow the thresholds", {
  expect_equal(RiboPilot:::chimera_label(0.29), "clean")
  expect_equal(RiboPilot:::chimera_label(0.3), "minor")
  expect_equal(RiboPilot:::chimera_label(0.9), "minor")
  expect_equal(RiboPilot:::chimera_label(0.1, severe = TRUE), "severe")
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
  expect_equal(r$label, "severe")
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

test_that("channel 2 drops to no_taxonomy when nothing resolves, even with a real dbarg", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  seq <- paste(rep("ACGT", 1000), collapse = "")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, ?)", params = list(seq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    lines <- vapply(ids, function(w) sprintf("NC_1\thit\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t4000\t%s", w),
      character(1))
    writeLines(lines, out); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = integer(0), lineage = character(0), lineage_taxids = character(0))
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel", taxdump_dir = NULL)
  channels <- jsonlite::fromJSON(r$channels)
  expect_equal(channels$windows$availability, "no_taxonomy")
  expect_equal(channels$windows$score, 0)
})

test_that("scaffold column help and badge exist", {
  expect_true(nzchar(RiboPilot:::SCAFFOLD_COL_HELP$chimera))
  js <- RiboPilot:::rt_chimera_badge()
  expect_s3_class(js, "JS_EVAL")
  expect_match(as.character(js), "severe")
})

test_that("synthetic fixtures get the expected label and a junction within tolerance", {
  skip_if_not_installed("withr")
  skip_if(!nzchar(Sys.which("blastn")) || !nzchar(Sys.which("makeblastdb")), "blast not installed")
  panel <- Sys.getenv("RIBOPILOT_TEST_PANEL", "")
  if (!nzchar(panel)) panel <- tryCatch(RiboPilot:::discover_panel(), error = function(e) NULL)
  skip_if(is.null(panel) || !nzchar(panel), "no local panel")
  # Channel 2's cross-taxon comparison needs a real pinned taxdump (see
  # RiboPilot:::have_taxdump); without one every window falls back to
  # identity-only buckets and host_fungal.fa cannot be told apart from
  # clean.fa. RIBOPILOT_TEST_TAXDUMP points at one; skip otherwise.
  taxdump <- Sys.getenv("RIBOPILOT_TEST_TAXDUMP", "")
  skip_if(!nzchar(taxdump) || !RiboPilot:::have_taxdump(taxdump), "no local taxdump")
  exp <- read.delim(test_path("fixtures/chimera/expected.tsv"), na.strings = "NA", stringsAsFactors = FALSE)
  for (i in seq_len(nrow(exp))) {
    seq <- paste(readLines(test_path("fixtures/chimera", exp$file[i]))[-1], collapse = "")
    con <- chimera_test_con("F")
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('F',1,1,?)", params = list(seq))
    r <- RiboPilot:::chimera_unit(con, "F", 1L, 1L, dbarg = panel, taxdump_dir = taxdump)
    DBI::dbDisconnect(con)
    expect_equal(r$label, exp$label[i], info = exp$file[i])
    if (!is.na(exp$junction[i])) {
      got <- as.integer(strsplit(r$junctions, ",", fixed = TRUE)[[1]])
      expect_true(length(got) > 0 && any(abs(got - exp$junction[i]) <= 300), info = exp$file[i])
    }
  }
})

test_that("chimera table column help covers every displayed column", {
  h <- RiboPilot:::CHIMERA_COL_HELP
  expect_true(all(c("start", "stop", "taxon", "species", "pident", "n_windows", "win_species", "lineage",
    "win_pident", "bitscore", "rank_distance", "depth_median", "clip_frac") %in% names(h)))
  expect_true(all(nzchar(unlist(h))))
})

test_that("the chimera panel shows the level and the reason in the open, and not run before a recheck", {
  panel <- as.character(RiboPilot:::chimera_panel())
  expect_match(panel, "chimera_result", fixed = TRUE)
  expect_false(grepl("chimera_plot", panel, fixed = TRUE))
  x <- as.character(RiboPilot:::chimera_result_text("minor", "18S resolves to Chordata: agree only to phylum"))
  expect_match(x, '<span class="mp-pill mp-pill-warning">minor</span>', fixed = TRUE)
  expect_match(x, "agree only to phylum", fixed = TRUE)
  expect_false(grepl("title=", x, fixed = TRUE))
  x <- as.character(RiboPilot:::chimera_result_text("clean", NA))
  expect_match(x, '<span class="mp-pill mp-pill-success">clean</span>', fixed = TRUE)
  x <- as.character(RiboPilot:::chimera_result_text(NA, NA, state = "not run"))
  expect_match(x, "not run", fixed = TRUE)
  expect_match(x, "press Update", fixed = TRUE)
})

test_that("the junction table names the checks behind each junction and why", {
  ch <- list(
    genes = list(score = 0, junctions = list()),
    windows = list(score = 0.5, junctions = 2400, reason = "1-2300 Actinopteri (Salmo salar, 99.0%); 2501-4500 Eurotiomycetes (Penicillium sp., 97.0%); different class"),
    ribotyper = list(score = 0, junctions = list(), reason = "ribotyper passed"),
    reads = list(score = 0.7, junctions = c(750, 2500), reason = "depth step at 750; clip pileup at 2500"))
  j <- RiboPilot:::chimera_junction_table("750,2400", as.character(jsonlite::toJSON(ch, auto_unbox = TRUE)))
  expect_equal(j$position, c(750L, 2400L))
  expect_equal(j$checks, c("reads", "windows, reads"))
  expect_match(j$why[1], "depth changes", fixed = TRUE)
  expect_match(j$why[1], "read evidence alone", fixed = TRUE)
  expect_match(j$why[2], "Eurotiomycetes", fixed = TRUE)
  expect_match(j$why[2], "clipped", fixed = TRUE)
  expect_equal(nrow(RiboPilot:::chimera_junction_table("", NA_character_)), 0L)
})

test_that("chimera_unit(review = TRUE) scores the copy and writes only its mirrors", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  orig <- withr::with_seed(3, paste(sample(c("A", "C", "G", "T"), 4000, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence) VALUES ('S1', 1, 1, 4000, ?),
    ('S1', 1, 2, 4000, ?)", params = list(orig, orig))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 2)")
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = '33208;7711;186623;8030;8016;8017;8018' WHERE ID = 'S1'")
  DBI::dbExecute(con, "UPDATE annotate SET chimeric = 0")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop)
    VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1',1,1,?,?)",
    params = list(paste(rep(20, 4000), collapse = ","), paste(c(rep(0, 2199), 9, rep(0, 1800)), collapse = ",")))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, reason) VALUES ('S1',1,1,0,'clean','original')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident)
    VALUES ('S1',1,1,1,500,99), ('S1',1,1,3501,4000,99)")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  # the copy: first 1000 bp trimmed off (18S cut to 1..800), chimeric per its own ITSx run, its own hits
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?, length = 3000, trim_start = 1001, trim_end = 4000,
    chimeric = 1", params = list(substr(orig, 1001, 4000)))
  DBI::dbExecute(con, "UPDATE review_annotations SET start = max(start - 1000, 1), stop = stop - 1000")
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  DBI::dbAppendTable(con, "review_blast_hits", data.frame(ID = "S1", path = 1L, scaffold = 1L,
    hit_rank = 1:6, hsp = 1L, accession = paste0("A", 1:6), species = "Salmo salar",
    lineage = "Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar",
    lineage_taxids = fish, pident = 99, bitscore = 900,
    qstart = rep(c(1L, 1500L), each = 3), qend = rep(c(800L, 3000L), each = 3)))
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    writeLines(sprintf("NC_1\tSalmo salar 18S\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t1800\t%s", ids), out)
    0L
  }
  fake_rt <- function(fasta, outdir, ...) {
    writeLines(c("#idx target class strnd p/f unexpected_features", "1 u18S SSU.Eukarya plus PASS -",
      "2 u28S LSU.Eukarya plus PASS -"), file.path(outdir, "rt.ribotyper.short.out")); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = 8030L,
    lineage = "Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar", lineage_taxids = fish)
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, chimera_run_ribotyper = fake_rt,
    tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  originals <- function() list(DBI::dbGetQuery(con, "SELECT ID, path, scaffold, score, label, reason FROM chimera"),
    DBI::dbGetQuery(con, "SELECT ID, path, scaffold, start, stop FROM chimera_windows ORDER BY start"))
  before <- originals()

  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel", review = TRUE)
  expect_identical(originals(), before)
  ch <- jsonlite::fromJSON(r$channels)
  expect_equal(ch$reads$junctions, 1200L)
  expect_equal(ch$ribotyper$score, 0.2)
  expect_equal(jsonlite::fromJSON(r$gene_taxa)[, c("element", "n_hits")],
    data.frame(element = c("18S", "28S"), n_hits = c(3L, 3L)))
  rc <- DBI::dbGetQuery(con, "SELECT junctions, reason FROM review_chimera")
  expect_equal(rc$junctions, "1200")
  expect_equal(rc$reason, r$reason)
  w <- DBI::dbGetQuery(con, "SELECT start, stop, clip_frac FROM review_chimera_windows ORDER BY start")
  expect_equal(w$start, RiboPilot:::chimera_tile(3000)$start)
  expect_equal(max(w$stop), 3000L)
  expect_equal(w$start[w$clip_frac > 0], c(1001L, 1101L))

  # no copy: nothing to score in review mode
  expect_null(RiboPilot:::chimera_unit(con, "S1", 1L, 2L, dbarg = "panel", review = TRUE))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_chimera WHERE scaffold = 2")$n, 0L)

  # original mode scores the original and leaves the copy alone
  RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel")
  expect_equal(DBI::dbGetQuery(con, "SELECT junctions, reason FROM review_chimera"), rc)
  expect_equal(DBI::dbGetQuery(con, "SELECT junctions FROM chimera")$junctions, "2200")
  expect_equal(DBI::dbGetQuery(con, "SELECT start FROM chimera_windows ORDER BY start")$start,
    RiboPilot:::chimera_tile(4000)$start)
  expect_equal(DBI::dbGetQuery(con, "SELECT start, stop, clip_frac FROM review_chimera_windows ORDER BY start"), w)
})

test_that("the chimera panel reads the review layer on Annotation Review and the originals on the Assemble step", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, reason, junctions, db_source)
    VALUES ('S1',1,1,'severe','original reason','3500','/x/orig.fa')")
  fixture_mark_done(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label, reason, junctions, db_source)
    VALUES ('S1',1,1,'minor','copy reason','1200','/x/copy.fa')")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    # not run: the row a discarded recheck left stays hidden
    expect_match(output$chimera_result$html, "not run", fixed = TRUE)
    expect_false(grepl("copy reason", output$chimera_result$html, fixed = TRUE))

    DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2")
    refresh(refresh() + 1)
    session$flushReact()
    expect_match(output$chimera_result$html, "minor", fixed = TRUE)
    expect_match(output$chimera_result$html, "copy reason", fixed = TRUE)
    expect_match(output$chimera_junctions$html, '"position":[1200]', fixed = TRUE)

    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()
    expect_match(output$chimera_result$html, "severe", fixed = TRUE)
    expect_match(output$chimera_junctions$html, '"position":[3500]', fixed = TRUE)
  })
})

test_that("clicking a chimera table row zooms the open modal's viewer", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  msgs <- list()
  local_mocked_bindings(seqview_goto = function(session, canvas, pos1, pos2)
    msgs[[length(msgs) + 1]] <<- list(canvas, pos1, pos2))
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(chim_zoom = list(pos1 = 100, pos2 = 300, nonce = 1))
    session$setInputs(mode = "Assemble")
    session$setInputs(chim_zoom = list(pos1 = 2625, pos2 = 2625, nonce = 2))
  })
  expect_equal(msgs[[1]], list("seqv_ann-canvas", 100, 300))
  expect_equal(msgs[[2]], list("seqv_asm-canvas", 2625, 2625))
  expect_match(as.character(RiboPilot:::rt_zoom_row("start", "stop")), "chim_zoom", fixed = TRUE)
})

test_that("the reason text leaves out the empty match of a no-hit segment", {
  seg <- data.frame(start = c(1, 501), stop = c(750, 2500), taxon = c("nohit", "Actinopterygii"),
    species = c(NA, "Salmo salar"), pident = c(NA, 99), stringsAsFactors = FALSE)
  r <- RiboPilot:::chimera_reason(seg, list(), integer(0))
  expect_equal(r, "1-750 nohit; 501-2500 Actinopterygii (Salmo salar, 99.0%); no junction")
})
