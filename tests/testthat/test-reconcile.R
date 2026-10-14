# barrnap/ITSx parsing and reconciliation into a canonical feature model.

test_that("parse_barrnap maps elements and captures model coverage", {
  gff <- c(
    "s\tpybarrnap\trRNA\t8\t1820\t0\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA",
    "s\tpybarrnap\trRNA\t8\t1424\t0\t-\t.\tName=18S_rRNA;product=18S ribosomal RNA (partial);Note=aligned only 75.82 percent of the 18S ribosomal RNA",
    "s\tpybarrnap\trRNA\t1\t120\t0\t+\t.\tName=5_8S_rRNA;product=5.8S ribosomal RNA")
  d <- RiboPilot:::parse_barrnap(gff)
  expect_equal(nrow(d), 3)
  expect_equal(d$element, c("18S", "18S", "5.8S"))
  expect_equal(d$model_cov[1], 100)      # full hit
  expect_equal(d$model_cov[2], 75.82)    # partial
  expect_equal(d$strand[2], "-")
  expect_true(is.data.frame(RiboPilot:::parse_barrnap(character(0))))
})

test_that("parse_itsx reads regions and detects chimeras", {
  line <- "s\t3000 bp.\tSSU: 1-1800\tITS1: 1801-2000\t5.8S: 2001-2160\tITS2: 2161-2400\tLSU: 2401-2900\tChimeric! "
  d <- RiboPilot:::parse_itsx(line)
  expect_equal(nrow(d), 5)
  expect_true(all(c("18S", "ITS1", "5.8S", "ITS2", "28S") %in% d$element))
  expect_true(RiboPilot:::itsx_is_chimeric(line))
  expect_false(RiboPilot:::itsx_is_chimeric(sub("\tChimeric! ", "", line)))
})

test_that("reconcile preserves strand and fills consensus", {
  bar <- data.frame(element = c("18S", "5.8S", "28S"), start = c(1, 2201, 2401),
    stop = c(1800, 2360, 6000), strand = "-", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(element = c("ITS1", "ITS2"), start = c(1801, 2361),
    stop = c(2200, 2400), stringsAsFactors = FALSE)
  rec <- RiboPilot:::reconcile_features(bar, its)
  expect_true(all(rec$strand == "-"))                       # real strand kept, not forced +
  expect_equal(rec$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
})

test_that("reconcile clips a spacer that overlaps an anchor", {
  bar <- data.frame(element = c("18S", "5.8S", "28S"), start = c(935, 3153, 3953),
    stop = c(2749, 3308, 8027), strand = "+", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(element = c("ITS1", "ITS2"), start = c(2751, 3197),
    stop = c(3151, 3948), stringsAsFactors = FALSE)              # ITS2 overlaps 5.8S
  rec <- RiboPilot:::reconcile_features(bar, its)
  o <- rec[order(rec$start), ]
  expect_false(any(o$start[-1] <= o$stop[-nrow(o)]))            # no overlaps
  expect_equal(o$start[o$element == "ITS2"], 3309)             # clipped to 5.8S.stop+1
})

test_that("reconcile derives a missing spacer by subtraction", {
  bar <- data.frame(element = c("18S", "5.8S"), start = c(1, 2201),
    stop = c(1800, 2360), strand = "+", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame()   # no ITSx: ITS1 derived between 18S and 5.8S
  rec <- RiboPilot:::reconcile_features(bar, its)
  its1 <- rec[rec$element == "ITS1", ]
  expect_equal(its1$start, 1801)
  expect_equal(its1$stop, 2200)
  expect_equal(its1$tool, "derived")
})

test_that("parse_barrnap carries the GFF seqid", {
  gff <- c(
    "ACC1\tbarrnap\trRNA\t10\t1800\t.\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA",
    "ACC2\tbarrnap\trRNA\t20\t1900\t.\t-\t.\tName=18S_rRNA;product=18S ribosomal RNA")
  d <- RiboPilot:::parse_barrnap(gff)
  expect_equal(nrow(d), 2L)
  expect_equal(d$seqid, c("ACC1", "ACC2"))
  expect_equal(d$element, c("18S", "18S"))
  expect_equal(d$start, c(10L, 20L))
  expect_equal(d$strand, c("+", "-"))
})

test_that("parse_itsx reads every record, not just the first", {
  pos <- c(
    "ACC1\t7000 bp.\tSSU: 1-1800\tITS1: 1801-2100\t5.8S: 2101-2258\tITS2: 2259-2500\tLSU: 2501-6900",
    "ACC2\t4000 bp.\tSSU: Not found\tITS1: Not found\t5.8S: Not found\tITS2: Not found\tLSU: 1-3950")
  d <- RiboPilot:::parse_itsx(pos)
  expect_equal(sort(unique(d$seqid)), c("ACC1", "ACC2"))
  expect_equal(nrow(d[d$seqid == "ACC1", ]), 5L)
  expect_equal(nrow(d[d$seqid == "ACC2", ]), 1L)
  expect_equal(d$element[d$seqid == "ACC2"], "28S")
})

test_that("a single-record positions line still parses exactly as before", {
  pos <- "S1\t7000 bp.\tSSU: 1-1800\tITS1: 1801-2100\t5.8S: 2101-2258\tITS2: 2259-2500\tLSU: 2501-6900"
  d <- RiboPilot:::parse_itsx(pos)
  expect_equal(d$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
  expect_equal(d$start, c(1L, 1801L, 2101L, 2259L, 2501L))
})

test_that("reconcile moves a gene end up to 5 bp to meet its ITSx spacer", {
  bar <- data.frame(element = c("18S", "28S"), start = c(1, 2362), stop = c(1799, 6000),
    strand = "+", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(element = c("ITS1", "5.8S", "ITS2"), start = c(1801, 2001, 2161),
    stop = c(2000, 2160, 2357), stringsAsFactors = FALSE)
  rec <- RiboPilot:::reconcile_features(bar, its)
  expect_equal(rec$stop[rec$element == "18S"], 1800)          # 1 bp gap closed
  expect_equal(rec$start[rec$element == "28S"], 2358)         # 4 bp gap closed
  expect_equal(rec$tool[rec$element == "18S"], "barrnap+ITSx")
  expect_equal(rec$start[-1], rec$stop[-nrow(rec)] + 1)       # every feature abuts the next

  bar$stop[1] <- 1794                                         # 6 bp: left for curation to flag
  rec2 <- RiboPilot:::reconcile_features(bar, its)
  expect_equal(rec2$stop[rec2$element == "18S"], 1794)
  expect_equal(rec2$tool[rec2$element == "18S"], "barrnap")
})

test_that("reconcile anchors 18S and 28S on the copies nearest the ITS region", {
  # a 28S tail from the previous repeat opens the scaffold, and a second 28S
  # copy runs off its end
  bar <- data.frame(element = c("28S", "18S", "5.8S", "28S"),
    start = c(1L, 2001L, 4101L, 4501L), stop = c(900L, 3800L, 4258L, 6000L),
    strand = "+", model_cov = c(25, 100, 100, 40), stringsAsFactors = FALSE)
  its <- data.frame(element = c("ITS1", "5.8S", "ITS2"), start = c(3801L, 4101L, 4259L),
    stop = c(4100L, 4258L, 4500L), stringsAsFactors = FALSE)
  rec <- RiboPilot:::reconcile_features(bar, its)
  expect_equal(rec$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
  expect_equal(rec$start[rec$element == "28S"], 4501L)
  expect_equal(rec$stop[rec$element == "ITS2"], 4500L)
})

test_that("reconcile ignores a minority-strand anchor copy", {
  bar <- data.frame(element = c("18S", "5.8S", "28S", "28S"),
    start = c(1L, 2001L, 2300L, 2401L), stop = c(1800L, 2160L, 2350L, 6000L),
    strand = c("+", "+", "-", "+"), model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(element = c("ITS1", "5.8S", "ITS2"), start = c(1801L, 2001L, 2161L),
    stop = c(2000L, 2160L, 2400L), stringsAsFactors = FALSE)
  rec <- RiboPilot:::reconcile_features(bar, its)
  expect_equal(rec$start[rec$element == "28S"], 2401L)
  expect_true(all(rec$strand == "+"))
})

test_that("with no ITS region the longest copy of a gene is used", {
  bar <- data.frame(element = c("18S", "18S"), start = c(1L, 5001L), stop = c(300L, 6800L),
    strand = "+", model_cov = c(15, 100), stringsAsFactors = FALSE)
  rec <- RiboPilot:::reconcile_features(bar, data.frame())
  expect_equal(rec$start, 5001L)
})

# barrnap GFF and tuned-ITSx positions for unit S1.1.1, as annotate_tools writes them
reconcile_inputs <- function(dir) {
  gff <- file.path(dir, "unit.barrnap.gff")
  writeLines(c(
    "S1.1.1\tbarrnap:0.9\trRNA\t351\t2150\t0\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA",
    "S1.1.1\tbarrnap:0.9\trRNA\t2651\t2807\t0\t+\t.\tName=5_8S_rRNA;product=5.8S ribosomal RNA",
    "S1.1.1\tbarrnap:0.9\trRNA\t3991\t6310\t0\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA"), gff)
  pos <- file.path(dir, "unit.itsx.positions.txt")
  writeLines(paste0("S1.1.1\t6310 bp.\tSSU: 351-2150\tITS1: 2151-2650\t5.8S: 2651-2807\t",
    "ITS2: 2808-3990\tLSU: 3991-6310\tChimeric! "), pos)
  list(gff = gff, pos = pos)
}

test_that("reconcile_one writes chimeric and clears the stored ITSx check", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotate SET chimeric = 0, itsx_warnings = 'itsx_no_5.8S', itsx_origin = 'M',
    itsx_scores = 'M: 1 1e-10 40'")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, itsx_warnings, itsx_origin, itsx_scores)
    VALUES ('S1', 1, 2, '', 'M', 'M: 1 1e-10 40')")
  f <- reconcile_inputs(withr::local_tempdir())
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, gff = f$gff, positions = f$pos)
  an <- DBI::dbGetQuery(con, "SELECT annotate_switch, n_features, missing, chimeric, itsx_warnings,
    itsx_origin, itsx_scores FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")
  expect_equal(an$annotate_switch[1], 2L)
  expect_equal(an$n_features[1], 5L)
  expect_equal(an$missing[1], "")
  expect_equal(an$chimeric[1], 1L)
  expect_true(is.na(an$itsx_warnings[1]) && is.na(an$itsx_origin[1]) && is.na(an$itsx_scores[1]))
  expect_equal(an$itsx_origin[2], "M")
})

test_that("reconcile_one keeps the old features when the annotate update fails", {
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  f <- reconcile_inputs(withr::local_tempdir())
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, gff = f$gff, positions = f$pos)
  DBI::dbExecute(con,
    "CREATE TRIGGER block_annotate BEFORE UPDATE ON annotate BEGIN SELECT RAISE(ABORT, 'blocked'); END")
  expect_error(RiboPilot::reconcile_one(con, "S1", 1L, 1L), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations WHERE ID = 'S1'")$n, 5L)
})

test_that("rp_reconcile.R stores not_run when the ITSx check throws", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("pkgload")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!nzchar(src), "not a source checkout")
  con <- chimera_test_con()
  withr::defer(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  f <- reconcile_inputs(d)
  pre <- file.path(d, "unit.itsx_check")
  writeLines("ITSx run finished at Thu Oct  8 12:00:00 2026.", paste0(pre, ".summary.txt"))
  dir.create(paste0(pre, ".positions.txt"))   # unreadable: the check itself throws
  res <- processx::run("Rscript",
    c(file.path(system.file("nextflow", "bin", package = "RiboPilot"), "rp_reconcile.R"),
      DBI::dbGetInfo(con)$dbname, "S1", "1", "1", f$gff, f$pos, pre),
    env = c("current", RIBOPILOT_SRC = src), error_on_status = FALSE)
  expect_equal(res$status, 0L, info = res$stderr)
  an <- DBI::dbGetQuery(con, "SELECT n_features, itsx_warnings, itsx_origin FROM annotate
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  expect_equal(an$n_features, 5L)
  expect_match(an$itsx_warnings, "^not_run\\(ITSx check failed: .+\\)$")
  expect_true(is.na(an$itsx_origin))
})
