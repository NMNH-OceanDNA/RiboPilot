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
