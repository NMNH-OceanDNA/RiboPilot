FB_DANIO <- c(lin = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
  tax = "33208;7711;186623;7952;2743709;7954;7955")
FB_DANIO2 <- c(lin = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio sp.",
  tax = "33208;7711;186623;7952;2743709;7954;1841699")
FB_GANNET <- c(lin = "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;Morus bassanus",
  tax = "33208;7711;8782;30410;30446;37577;37578")
FB_PLANT <- c(lin = "Viridiplantae;Streptophyta;Magnoliopsida;Rosales;Moraceae;Morus;Morus alba",
  tax = "33090;35493;3398;3744;3487;3497;3498")

fb_hit <- function(el, acc, title, bits, taxon, evalue = 1e-50) {
  data.frame(accession = acc, title = title, pident = 99, qcovs = 100, evalue = evalue,
    staxid = as.integer(sub(".*;", "", taxon[["tax"]])),
    bitscore = bits, qseqid = el, lineage = taxon[["lin"]], lineage_taxids = taxon[["tax"]],
    stringsAsFactors = FALSE)
}
fb_targets <- function(el, len) data.frame(element = el, start = 1L, stop = len, length = len,
  sequence = "", stringsAsFactors = FALSE)
fb_call <- function(hits, el = "ITS1", len = 300L, sample = FB_DANIO[["tax"]]) {
  RiboPilot:::feature_blast_calls(hits, fb_targets(el, len), sample)$calls
}

test_that("feature_blast_task() splits at 50 and 100 bp", {
  expect_equal(RiboPilot:::feature_blast_task(c(NA, 10, 49, 50, 99, 100, 5000)),
    c(NA, NA, NA, "blastn", "blastn", "megablast", "megablast"))
})

test_that("short, unmatched, and non-rDNA features get their status", {
  r <- RiboPilot:::feature_blast_calls(data.frame(), fb_targets(c("5.8S", "ITS1"), c(40L, 80L)),
    FB_DANIO[["tax"]])$calls
  expect_equal(r$status, c("too_small", "no_hit"))
  expect_equal(r$note[1], "no BLAST, feature too small")
  expect_equal(r$low_conf, c(0L, 1L))
  expect_equal(fb_call(fb_hit("ITS1", "A1", "Danio rerio ITS1", 500, FB_DANIO, evalue = 1))$status, "no_hit")
  expect_equal(fb_call(fb_hit("ITS1", "A1", "Danio rerio mitochondrion, complete genome", 500, FB_DANIO))$status,
    "non_rdna")
})

test_that("a distinct best hit reports the rank it shares with the sample", {
  h <- rbind(fb_hit("ITS1", "A1", "Morus bassanus internal transcribed spacer 1", 500, FB_GANNET),
    fb_hit("ITS1", "A2", "Danio rerio internal transcribed spacer 1", 300, FB_DANIO))
  r <- fb_call(h)
  expect_equal(r$status, "ok")
  expect_equal(r$match_rank, "phylum")
  expect_equal(r$species, "Morus bassanus")
  expect_true(is.na(r$n_tied))
  expect_equal(fb_call(fb_hit("ITS1", "A1", "Danio rerio 18S ribosomal RNA", 500, FB_DANIO))$match_rank, "species")
})

test_that("a hit outside the sample's phylum is off-target", {
  r <- fb_call(fb_hit("ITS2", "P1", "Morus alba internal transcribed spacer 2", 500, FB_PLANT), el = "ITS2")
  expect_equal(r$status, "off_target")
  expect_equal(r$match_rank, "none")
})

test_that("tied hits collapse to their shared clade", {
  h <- rbind(fb_hit("5.8S", "A1", "Danio rerio 5.8S ribosomal RNA", 300, FB_DANIO),
    fb_hit("5.8S", "A2", "Morus bassanus 5.8S ribosomal RNA", 300, FB_GANNET))
  r <- fb_call(h, el = "5.8S", len = 160L)
  expect_equal(r$status, "ok")
  expect_equal(r$species, "Chordata (2 tied hits)")
  expect_equal(r$n_tied, 2L)
  expect_equal(r$match_rank, "phylum")
  # same genus ties still name the genus
  r2 <- fb_call(rbind(fb_hit("18S", "A1", "Danio rerio 18S rRNA", 900, FB_DANIO),
    fb_hit("18S", "A2", "Danio sp. 18S rRNA", 900, FB_DANIO2)), el = "18S", len = 1800L)
  expect_equal(r2$species, "Danio (2 tied hits)")
  expect_equal(r2$match_rank, "genus")
})

test_that("a tie across kingdoms is not off-target when a tied hit shares the phylum", {
  h <- rbind(fb_hit("5.8S", "P1", "Morus alba 5.8S ribosomal RNA", 300, FB_PLANT),
    fb_hit("5.8S", "A1", "Danio rerio 5.8S ribosomal RNA", 300, FB_DANIO))
  r <- fb_call(h, el = "5.8S", len = 160L)
  expect_equal(r$status, "ok")
  expect_equal(r$n_tied, 2L)
  h2 <- rbind(fb_hit("5.8S", "P1", "Morus alba 5.8S ribosomal RNA", 300, FB_PLANT),
    fb_hit("5.8S", "P2", "Morus alba 5.8S ribosomal RNA", 299, FB_PLANT))
  expect_equal(fb_call(h2, el = "5.8S", len = 160L)$status, "off_target")
})

test_that("without a sample lineage nothing is off-target", {
  r <- fb_call(fb_hit("ITS2", "P1", "Morus alba ITS2", 500, FB_PLANT), el = "ITS2", sample = NA_character_)
  expect_equal(r$status, "ok")
  expect_equal(r$note, "sample lineage unknown")
  expect_true(is.na(r$match_rank))
})

test_that("only the top hits per feature are kept", {
  h <- do.call(rbind, lapply(1:15, function(i) fb_hit("18S", paste0("A", i), "Danio rerio 18S rRNA", 1000 - i,
    FB_DANIO)))
  k <- RiboPilot:::feature_blast_calls(h, fb_targets("18S", 1800L), FB_DANIO[["tax"]])$hits
  expect_equal(k$hit_rank, 1:10)
  expect_equal(k$accession[1], "A1")
})

test_that("feature_blast_warning_text() lists problems in element order", {
  fb <- data.frame(ID = "S1", path = 1L, scaffold = c(1L, 1L, 1L, 1L, 2L),
    element = c("ITS2", "5.8S", "ITS1", "18S", "18S"),
    status = c("off_target", "too_small", "no_hit", "ok", "ok"),
    lineage = c(FB_PLANT[["lin"]], NA, NA, NA, NA), stringsAsFactors = FALSE)
  w <- RiboPilot:::feature_blast_warning_text(fb)
  expect_equal(nrow(w), 1L)
  expect_equal(w$feature_blast_warnings,
    "ITS1: no hit; 5.8S: no BLAST, feature too small; ITS2: off-target (Streptophyta)")
  expect_equal(nrow(RiboPilot:::feature_blast_warning_text(fb[fb$status == "ok", ])), 0L)
})

fb_seed <- function(con) {
  seq <- withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 3000, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, ?)",
    params = list(seq))
  DBI::dbExecute(con, "UPDATE samples SET lineage_taxids = ? WHERE ID = 'S1'", params = list(FB_DANIO[["tax"]]))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES
    ('S1',1,1,'18S',1,1800), ('S1',1,1,'ITS1',1801,1880), ('S1',1,1,'5.8S',1881,1920),
    ('S1',1,1,'ITS2',1921,2200), ('S1',1,1,'28S',2201,3000), ('S1',1,1,'IGS',2990,3000)")
  seq
}

test_that("feature_blast_unit() searches by length and stores calls and hits", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  seq <- fb_seed(con)
  calls <- list()
  fake_run <- function(query, dbarg, task, out, ...) {
    q <- readLines(query)
    ids <- sub("^>", "", q[c(TRUE, FALSE)])
    calls[[task]] <<- ids
    if (task == "megablast") expect_equal(q[2], substring(seq, 1, 1800))
    rows <- list(
      "18S" = fb_hit("18S", "A1", "Danio rerio 18S ribosomal RNA", 3000, FB_DANIO),
      "ITS1" = fb_hit("ITS1", "P1", "Morus alba internal transcribed spacer 1", 100, FB_PLANT))
    do.call(rbind, rows[intersect(ids, names(rows))]) %||% data.frame()
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = c(7955L, 3498L),
    lineage = c(FB_DANIO[["lin"]], FB_PLANT[["lin"]]), lineage_taxids = c(FB_DANIO[["tax"]], FB_PLANT[["tax"]]))
  testthat::local_mocked_bindings(feature_blast_run = fake_run, tk_ranked_lineage = fake_lin,
    .package = "RiboPilot")
  RiboPilot:::feature_blast_unit(con, "S1", 1L, 1L, dbarg = "panel")
  expect_equal(calls$megablast, c("18S", "ITS2", "28S"))
  expect_equal(calls$`dc-megablast`, c("ITS2", "28S"))
  expect_equal(calls$blastn, "ITS1")
  fb <- DBI::dbGetQuery(con, "SELECT element, status, task, low_conf, db_source FROM feature_blast ORDER BY element")
  expect_equal(fb$element, c("18S", "28S", "5.8S", "ITS1", "ITS2"))
  expect_equal(fb$status[fb$element == "5.8S"], "too_small")
  expect_equal(fb$status[fb$element == "ITS2"], "no_hit")
  expect_equal(fb$low_conf[fb$element == "ITS1"], 1L)
  expect_equal(fb$status[fb$element == "ITS1"], "off_target")
  expect_equal(fb$status[fb$element == "18S"], "ok")
  expect_true(all(fb$db_source == "panel"))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM feature_blast_hits")$n, 2L)
})

test_that("feature_blast_unit() stores not_run without a panel or when blastn fails, and rewrites old rows", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  fb_seed(con)
  RiboPilot:::feature_blast_unit(con, "S1", 1L, 1L, dbarg = NULL)
  fb <- DBI::dbGetQuery(con, "SELECT element, status, note FROM feature_blast")
  expect_equal(sort(fb$status), c("not_run", "not_run", "not_run", "not_run", "too_small"))
  expect_true(all(fb$note[fb$status == "not_run"] == "no local BLAST panel"))
  testthat::local_mocked_bindings(feature_blast_run = function(...) stop("boom"), .package = "RiboPilot")
  RiboPilot:::feature_blast_unit(con, "S1", 1L, 1L, dbarg = "panel")
  fb <- DBI::dbGetQuery(con, "SELECT status, note FROM feature_blast WHERE status = 'not_run'")
  expect_equal(nrow(fb), 4L)
  expect_match(fb$note[1], "BLAST failed: boom")
})

test_that("feature_blast_unit(review = TRUE) reads and writes the review copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  fb_seed(con)
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 1)")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE element != 'ITS1'")
  RiboPilot:::feature_blast_unit(con, "S1", 1L, 1L, dbarg = NULL, review = TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT element FROM review_feature_blast")$element, "ITS1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM feature_blast")$n, 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT element FROM review_feature_blast_v")$element, "ITS1")
})
