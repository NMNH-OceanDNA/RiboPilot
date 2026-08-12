# Scientific defaults that decide what data reaches the curator: BLAST hit
# thresholds, read trimming, and which reads a coverage depth is allowed to
# count. Each is a default a user can retune, so pin the shipped value.

mk_default_db <- function() {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  fa <- file.path(dir, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  map <- file.path(dir, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  db <- file.path(dir, ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  withr::defer(DBI::dbDisconnect(con), envir = parent.frame())
  con
}

test_that("the default fastp options quality-trim and set a length floor", {
  skip_if_not_installed("withr")
  con <- mk_default_db()
  f <- DBI::dbGetQuery(con, "SELECT fastp FROM pre_opts WHERE pre_opts = 'default'")$fastp
  expect_length(f, 1)
  # 3' quality trimming: without it, degraded read tails reach the assembler and
  # become spurious variant groups at rDNA copy number
  expect_match(f, "--cut_tail", fixed = TRUE)
  q <- as.numeric(sub(".*--cut_mean_quality\\s+(\\d+).*", "\\1", f))
  expect_true(!is.na(q) && q >= 20)
  # fastp's own default length floor is 15 bp
  l <- as.numeric(sub(".*(?:^|\\s)-l\\s+(\\d+).*", "\\1", f))
  expect_true(!is.na(l) && l >= 50)
  # adapter handling is not dropped in the process
  expect_match(f, "--detect_adapter_for_pe", fixed = TRUE)
  # duplicate reads are real copy number for rDNA, never an artefact to remove
  expect_match(f, "--dont_eval_duplication", fixed = TRUE)
  expect_false(grepl("--dedup", f, fixed = TRUE))
})

test_that("the default remote BLAST carries a significance floor but no identity floor", {
  skip_if_not_installed("withr")
  con <- mk_default_db()
  x <- DBI::dbGetQuery(con,
    "SELECT extra_opts FROM blast_opts WHERE blast_opts = 'default'")$extra_opts
  expect_length(x, 1)
  # blastn's own default is -evalue 10, which keeps hits expected by chance
  e <- as.numeric(sub(".*-evalue\\s+(\\S+).*", "\\1", x))
  expect_true(!is.na(e) && e <= 1e-10)
  # One search feeds both the reference panel and the contamination screen, and
  # blast_screen is written from the unfiltered distribution. An identity floor
  # here would blind the screen to the distant cross-kingdom hits it exists to
  # catch, so the panel's floor lives in parse_blast_candidates() instead.
  expect_false(grepl("-perc_identity", x, fixed = TRUE))
})

test_that("the reference panel floor filters candidates without touching the screen", {
  hits <- data.frame(
    saccver = c("A1", "B1", "C1"),
    stitle = paste(c("Danio rerio", "Danio rerio", "Danio rerio"),
      "18S ribosomal RNA gene"),
    pident = c(99, 70, 99),          # B1 too divergent to align against
    qcovs = c(90, 90, 10),           # C1 covers a tenth of the query
    evalue = c(1e-99, 1e-30, 1e-30),
    staxids = "7955",
    bitscore = c(2000, 900, 400),
    length = c(1800, 1700, 200),
    qstart = 1, qend = c(1800, 1700, 200), sstart = 1, send = c(1800, 1700, 200),
    stringsAsFactors = FALSE)
  f <- withr::local_tempfile()
  utils::write.table(hits, f, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  parsed <- RiboPilot:::parse_blast_hits_file(f)
  keep <- !is.na(parsed$pident) & parsed$pident >= RiboPilot:::REF_CANDIDATE_MIN_PIDENT &
    !is.na(parsed$qcovs) & parsed$qcovs >= RiboPilot:::REF_CANDIDATE_MIN_QCOVS
  expect_equal(parsed$accession[keep], "A1")
  # all three still exist in the parsed distribution the screen reads
  expect_equal(nrow(parsed), 3)
})

nf_process_body <- function(file, name) {
  lines <- readLines(system.file("nextflow", "modules", file, package = "RiboPilot"))
  s <- grep(paste0("^process\\s+", name, "\\s*\\{"), lines)
  e <- s - 1 + which(lines[s:length(lines)] == "}")[1]
  lines[s:e]
}

test_that("blast_genbank passes the tunable options to both searches", {
  body <- nf_process_body("blast.nf", "blast_genbank")
  # thresholds live in blast_opts.extra_opts so a user can retune them; the
  # dc-megablast fallback must not silently search without them
  expect_length(grep("blastn -remote", body), 2)
  expect_length(grep("\\$\\{XO\\}", body[!grepl("^\\s*#", body)]), 2)
  expect_true(any(grepl("SELECT COALESCE\\(b.extra_opts", body)))
})

test_that("coverage drops multi-mapping reads before depth is counted", {
  body <- nf_process_body("coverage.nf", "coverage")
  txt <- paste(body, collapse = "\n")
  # the conserved 18S/28S core maps equally well to every scaffold and bowtie2
  # breaks those ties at random, so unfiltered depths are partly noise
  q <- as.numeric(sub("(?s).*samtools view[^|]*-q\\s+(\\d+).*", "\\1", txt, perl = TRUE))
  expect_true(!is.na(q) && q >= 10)
  # the filter has to sit inside the bowtie2 -> sort pipe, not after it
  expect_match(txt, "bowtie2 -x ref(?s).*\\|\\s*samtools view[^|]*\\|\\s*samtools sort", perl = TRUE)
  expect_match(txt, "errorStrategy 'ignore'", fixed = TRUE)
  expect_match(txt, "sqlite3 -cmd '.timeout", fixed = TRUE)
})

test_that("coverage maps one path at a time, so the MAPQ filter cannot empty the sample", {
  body <- nf_process_body("coverage.nf", "coverage")
  txt <- paste(body, collapse = "\n")
  # Paths are alternative resolutions of one molecule. Indexed together every
  # read is a multi-mapper, the MAPQ filter discards the lot, and the sample
  # ends up with no depth at all: measured at 2 surviving alignments out of
  # 1.65M read pairs on a 10-path sample.
  expect_match(txt, "SELECT DISTINCT path FROM assemblies", fixed = TRUE)
  expect_match(txt, "for P in \\$PATHS", fixed = TRUE)
  # the index must be built from ONE path's scaffolds, never the whole sample
  expect_match(txt, "AND path=\\$P", fixed = TRUE)
  expect_match(txt, "bowtie2-build --quiet path.fa ref", fixed = TRUE)
  # every path's depth accumulates into the one file the R step reads
  expect_match(txt, ">> depth.txt", fixed = TRUE)
  expect_false(grepl("bowtie2-build --quiet scaffolds.fa", txt, fixed = TRUE))
})

test_that("coverage records uncovered scaffolds as zero rather than omitting them", {
  body <- nf_process_body("coverage.nf", "coverage")
  txt <- paste(body, collapse = "\n")
  # `samtools depth -a` leaves out reference sequences with no alignments
  # entirely, so an emptied scaffold produced no depth rows, and
  # write_coverage_units() then wrote no coverage row for it at all. In the app
  # "covered by nothing" and "never measured" both render blank.
  expect_match(txt, "samtools depth -aa", fixed = TRUE)
  expect_false(grepl("samtools depth -a ", txt, fixed = TRUE))
})

test_that("an empty depth file still writes no coverage rows rather than erroring", {
  # samtools depth over a BAM the MAPQ filter emptied yields an empty file
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage (ID TEXT, path INT, scaffold INT,
    mean_depth REAL, median_depth REAL, min_depth INT, max_depth INT,
    depths TEXT, time_stamp INT)")
  f <- tempfile(); file.create(f); on.exit(unlink(f), add = TRUE)
  expect_equal(RiboPilot::write_coverage_units(con, "S1", f), 0L)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM coverage")), 0L)
})
