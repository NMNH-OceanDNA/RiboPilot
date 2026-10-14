# recheck_unit and its helpers: the per-contig recheck of a review copy
# (SPEC S7.3, S5.4).

RC_SEQ <- paste(readLines(test_path("fixtures", "chimera", "clean.fa"))[-1], collapse = "")
RC_UNIT <- "ID = 'S1' AND path = 1 AND scaffold = 1"

# GFF lines as (py)barrnap writes them; cov < 100 marks a partial hit.
rc_gff <- function(file, element, start, stop, cov = 100, strand = "+") {
  nm <- c("18S" = "18S_rRNA", "5.8S" = "5_8S_rRNA", "28S" = "28S_rRNA")[element]
  note <- ifelse(cov < 100, sprintf(" (partial);Note=aligned only %s percent of the model", cov), "")
  writeLines(paste("S1.1.1", "barrnap:0.9", "rRNA", start, stop, "0", strand, ".",
    paste0("Name=", nm, ";product=", element, " ribosomal RNA", note), sep = "\t"), file)
}

# An ITSx positions line for a contig of `len` bp whose 28S ends at `lsu_end`.
rc_positions <- function(file, len, lsu_end, chimeric = FALSE) {
  writeLines(paste0("S1.1.1\t", len, " bp.\tSSU: 351-2150\tITS1: 2151-2650\t5.8S: 2651-2807\t",
    "ITS2: 2808-3990\tLSU: 3991-", lsu_end, "\t", if (chimeric) "Chimeric! " else ""), file)
}

# A 6310 bp contig S1.1.1 annotated from barrnap and ITSx output, done, named
# by the seed database, with flat coverage and, unless hits = FALSE, two
# original subjects: A1.1 (HSPs 1-3000 and 5801-6310) and B2.1 (6101-6310).
rc_project <- function(hits = TRUE, env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  db <- file.path(d, ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE, seeds_db = fa, labels_db = fa)
  con <- RiboPilot:::open_project_db(db)
  withr::defer(DBI::dbDisconnect(con), envir = env)
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore,
      edited, blast_accession, blast_species, blast_pident, blast_source, time_stamp)
    VALUES ('S1', 1, 1, 'linear', ?, ?, 0, 0, 'SEED1', 'Danio rerio', 99.4, 'seeds', 1)",
    params = list(nchar(RC_SEQ), RC_SEQ))
  rc_gff(file.path(d, "o.gff"), c("18S", "5.8S", "28S"), c(351, 2651, 3991), c(2150, 2807, 6310))
  rc_positions(file.path(d, "o.pos"), 6310, 6310)
  RiboPilot::reconcile_one(con, "S1", 1L, 1L, file.path(d, "o.gff"), file.path(d, "o.pos"))
  DBI::dbExecute(con, paste("UPDATE annotate SET itsx_warnings = '' WHERE", RC_UNIT))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1', 1, 1, ?, ?)",
    params = list(paste(rep(30, 6310), collapse = ","), paste(rep(0, 6310), collapse = ",")))
  if (hits) {
    h <- data.frame(ID = "S1", path = 1L, scaffold = 1L, hit_rank = c(1L, 1L, 2L), hsp = c(1L, 2L, 1L),
      accession = c("A1.1", "A1.1", "B2.1"),
      title = c(rep("Danio rerio 18S ribosomal RNA gene, partial sequence", 2),
        "Danio aesculapii 28S ribosomal RNA gene, partial sequence"),
      species = c("Danio rerio", "Danio rerio", "Danio aesculapii"), staxid = c(7955L, 7955L, 1841699L),
      pident = c(99.5, 99.5, 98), qcovs = c(56, 56, 3), evalue = 0, bitscore = c(5000, 900, 380),
      aln_len = c(3000L, 510L, 210L), qstart = c(1L, 5801L, 6101L), qend = c(3000L, 6310L, 6310L),
      sstart = c(1L, 5801L, 1L), send = c(3000L, 6310L, 210L), slen = 7000L, is_nuclear_rrna = 1L,
      hit_class = "rdna", time_stamp = 1L)
    DBI::dbAppendTable(con, "blast_hits", h)
  }
  con
}

# Output of rp_rrna_tools.sh for a copy trimmed to 1..6000: barrnap finds the
# cut 28S at 70 percent of its model, and the tuned ITSx run calls it chimeric.
rc_tools <- function(status = c(barrnap = 0, itsx = 0, itsx_check = 0), finished = TRUE,
                     env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  if (!is.null(status)) writeLines(paste0(names(status), "=", status), file.path(d, "unit.tools_status"))
  rc_gff(file.path(d, "unit.barrnap.gff"), c("18S", "5.8S", "28S"), c(351, 2651, 3991),
    c(2150, 2807, 6000), cov = c(100, 100, 70))
  rc_positions(file.path(d, "unit.itsx.positions.txt"), 6000, 6000, chimeric = TRUE)
  rc_positions(file.path(d, "unit.itsx_check.positions.txt"), 6000, 6000)
  if (finished) writeLines("ITSx run finished at today.", file.path(d, "unit.itsx_check.summary.txt"))
  d
}

rc_rev <- function(con) {
  DBI::dbGetQuery(con, paste("SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE", RC_UNIT))$rev
}

rc_copy <- function(con) {
  DBI::dbGetQuery(con, paste("SELECT check_switch, check_note, checked_at, n_features, missing, chimeric,
    n_flags, curate_flags, itsx_warnings, blast_accession, blast_qcovs, blast_source, time_stamp
    FROM review_units WHERE", RC_UNIT))
}

rc_count <- function(con, tbl) {
  DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM", tbl, "WHERE", RC_UNIT))$n
}

rc_cov <- function(con, tbl = "review_annotations") {
  a <- DBI::dbGetQuery(con, paste("SELECT element, model_coverage FROM", tbl, "WHERE", RC_UNIT))
  stats::setNames(a$model_coverage, a$element)
}

# Every original table the recheck could touch, whole.
rc_originals <- function(con) {
  tb <- c("assemblies", "annotate", "annotations", "coverage", "screen", "chimera", "chimera_windows",
    "blast_hits", "blast_screen", "blast_ref_candidates", "blast_ref_alignment", "blast_searched", "path_scores")
  stats::setNames(lapply(tb, function(t) DBI::dbReadTable(con, t)), tb)
}

# A blastn stand-in: logs each call's arguments to RC_BLAST_LOG (one per line,
# then "--"), prints the file RC_MEGA or RC_DC for -task megablast or
# dc-megablast, and exits 1 when that variable is "fail".
rc_blastn <- function(env = parent.frame()) {
  skip_on_os("windows")
  d <- withr::local_tempdir(.local_envir = env)
  f <- file.path(d, "blastn")
  writeLines(c("#!/usr/bin/env bash",
    'prev=""; task=""',
    'for a in "$@"; do printf "%s\\n" "$a" >> "$RC_BLAST_LOG"; [ "$prev" = "-task" ] && task=$a; prev=$a; done',
    'echo "--" >> "$RC_BLAST_LOG"',
    'if [ "$task" = "dc-megablast" ]; then out=$RC_DC; else out=$RC_MEGA; fi',
    '[ "$out" = "fail" ] && exit 1',
    '[ -n "$out" ] && cat "$out"',
    "exit 0"), f)
  Sys.chmod(f, "0755")
  withr::local_envvar(RC_BLAST_LOG = file.path(d, "calls.log"), RC_MEGA = "", RC_DC = "", .local_envir = env)
  f
}

rc_calls <- function() {
  log <- Sys.getenv("RC_BLAST_LOG")
  if (!file.exists(log)) return(list())
  x <- readLines(log)
  unname(split(x[x != "--"], cumsum(x == "--")[x != "--"]))
}

rc_hit_lines <- function(file, qcovs, accession = "P1.1") {
  writeLines(paste(accession, "Danio rerio 18S ribosomal RNA gene, partial sequence", 99.1, qcovs,
    "0.0", 7955, 4000, 2000, 1, 2000, 1, 2000, 7000, "query", sep = "\t"), file)
  file
}

test_that("the recheck's BLAST settings match the local search in blast.nf", {
  nf <- paste(readLines(system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")),
    collapse = "\n")
  b <- RiboPilot:::RECHECK_BLAST
  for (s in c(sprintf('OUTFMT="%s"', b$outfmt), sprintf("-max_hsps %d -task %s", b$max_hsps, b$task),
              paste("-task", b$retry_task), sprintf("-lt %d ]", b$retry_below_qcovs),
              sprintf("COALESCE(b.max_target_seqs,%d)", b$default_max_target_seqs))) {
    expect_true(grepl(s, nf, fixed = TRUE), info = s)
  }
})

test_that("recheck_blast_local runs megablast with the sample's options and keeps a good result", {
  blastn <- rc_blastn()
  wd <- withr::local_tempdir()
  Sys.setenv(RC_MEGA = rc_hit_lines(file.path(wd, "mega.txt"), 80))
  out <- RiboPilot:::recheck_blast_local("ACGTACGT", "/db/panel /db/decoys", 7L,
    " -evalue 1e-20  -perc_identity 90 ", file.path(wd, "run"), threads = 2L, blastn = blastn)
  expect_equal(out, file.path(wd, "run", "hits.txt"))
  expect_equal(readLines(out), readLines(file.path(wd, "mega.txt")))
  expect_equal(readLines(file.path(wd, "run", "query.fa")), c(">query", "ACGTACGT"))
  calls <- rc_calls()
  expect_length(calls, 1)
  expect_equal(calls[[1]], c("-db", "/db/panel /db/decoys", "-query", file.path(wd, "run", "query.fa"),
    "-outfmt", RiboPilot:::RECHECK_BLAST$outfmt, "-max_target_seqs", "7", "-max_hsps", "5",
    "-task", "megablast", "-num_threads", "2", "-evalue", "1e-20", "-perc_identity", "90"))
})

test_that("recheck_blast_local retries with dc-megablast as blast.nf does", {
  blastn <- rc_blastn()
  wd <- withr::local_tempdir()
  mega <- rc_hit_lines(file.path(wd, "mega.txt"), 24, "MEGA.1")
  dc <- rc_hit_lines(file.path(wd, "dc.txt"), 90, "DC.1")
  empty <- file.path(wd, "empty.txt")
  file.create(empty)
  run <- function(m, d, sub) {
    Sys.setenv(RC_MEGA = m, RC_DC = d)
    unlink(Sys.getenv("RC_BLAST_LOG"))
    out <- RiboPilot:::recheck_blast_local("ACGT", "panel", 50L, "", file.path(wd, sub), blastn = blastn)
    list(lines = readLines(out), tasks = vapply(rc_calls(), function(a) a[which(a == "-task") + 1], ""))
  }
  # best qcovs under 25: the retry's hits replace megablast's
  r <- run(mega, dc, "a")
  expect_equal(r$tasks, c("megablast", "dc-megablast"))
  expect_equal(r$lines, readLines(dc))
  # no hits at all: the retry runs too
  r <- run(empty, dc, "b")
  expect_equal(r$tasks, c("megablast", "dc-megablast"))
  expect_equal(r$lines, readLines(dc))
  # a retry with no output keeps the sentinel, a failed retry keeps megablast's hits
  expect_equal(run(empty, empty, "c")$lines, "NO_SIGNIFICANT_HITS")
  expect_equal(run(mega, "fail", "d")$lines, readLines(mega))
  # 25 or more: no retry
  r <- run(rc_hit_lines(file.path(wd, "m25.txt"), 25), dc, "e")
  expect_equal(r$tasks, "megablast")
  # megablast itself failed: no hits file
  Sys.setenv(RC_MEGA = "fail")
  expect_null(RiboPilot:::recheck_blast_local("ACGT", "panel", 50L, "", file.path(wd, "f"), blastn = blastn))
})

test_that("map mode moves the original hits onto the copy's window", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 101L, 6000L)
  before <- rc_originals(con)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  h <- DBI::dbGetQuery(con, paste("SELECT accession, qstart, qend, aln_len, sstart, qcovs
    FROM review_blast_hits WHERE", RC_UNIT, "ORDER BY hit_rank, hsp"))
  expect_equal(h$accession, c("A1.1", "A1.1"))
  expect_equal(h$qstart, c(1L, 5701L))
  expect_equal(h$qend, c(2900L, 5900L))
  expect_equal(h$aln_len, c(2900L, 200L))
  expect_equal(h$sstart, c(NA_integer_, NA_integer_))
  expect_equal(h$qcovs, c(53, 53))
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "A1.1")
  expect_equal(u$blast_qcovs, 53)
  expect_equal(rc_count(con, "review_blast_screen"), 1)
  expect_equal(rc_originals(con), before)
})

test_that("map mode keeps the copied summary only when the original has no hits", {
  con <- rc_project(hits = FALSE)
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  # rows a discarded recheck left behind
  DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank, hsp, accession)
    VALUES ('S1', 1, 1, 1, 1, 'OLD.1')")
  DBI::dbExecute(con, "INSERT INTO review_blast_screen (ID, path, scaffold, accession) VALUES ('S1', 1, 1, 'OLD.1')")
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  expect_equal(rc_count(con, "review_blast_hits"), 0)
  expect_equal(rc_count(con, "review_blast_screen"), 0)
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "SEED1")
  expect_equal(u$blast_source, "seeds")

  # the original has hits, but none inside the window: the summary goes
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 3001L, 5800L)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "map"), "map")
  expect_equal(rc_count(con, "review_blast_hits"), 0)
  expect_true(is.na(rc_copy(con)$blast_accession))
})

test_that("local mode searches the panel with the sample's options and never NCBI", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  DBI::dbExecute(con, "UPDATE blast_opts SET extra_opts = '-evalue 1e-5', max_target_seqs = 9
    WHERE blast_opts = 'default'")
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('S1', 1, 1, 'x', 'local')")
  blastn <- rc_blastn()
  withr::local_envvar(PATH = paste(dirname(blastn), Sys.getenv("PATH"), sep = .Platform$path.sep))
  wd <- withr::local_tempdir()
  Sys.setenv(RC_MEGA = rc_hit_lines(file.path(wd, "mega.txt"), 95))
  before <- rc_originals(con)
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = "panel"), "local")
  a <- rc_calls()[[1]]
  expect_false("-remote" %in% a)
  expect_equal(a[which(a == "-max_target_seqs") + 1], "9")
  expect_equal(utils::tail(a, 2), c("-evalue", "1e-5"))
  expect_equal(DBI::dbGetQuery(con, paste("SELECT accession FROM review_blast_hits WHERE", RC_UNIT))$accession,
    "P1.1")
  u <- rc_copy(con)
  expect_equal(u$blast_accession, "P1.1")
  expect_equal(u$blast_source, "local_rrna")
  expect_equal(rc_originals(con), before)

  # a failed search falls back to the mapped original hits
  Sys.setenv(RC_MEGA = "fail")
  expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = "panel"), "map")
  expect_equal(rc_copy(con)$blast_accession, "A1.1")
})

test_that("local mode without a panel maps the hits and never calls blastn", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  blastn <- rc_blastn()
  withr::local_envvar(PATH = paste(dirname(blastn), Sys.getenv("PATH"), sep = .Platform$path.sep))
  for (dbarg in list(NULL, "", NA_character_)) {
    expect_equal(RiboPilot:::recheck_hits(con, "S1", 1L, 1L, "local", dbarg = dbarg), "map")
  }
  expect_length(rc_calls(), 0)
  expect_equal(rc_copy(con)$blast_accession, "A1.1")
})

test_that("recheck_model_coverage refreshes only barrnap features, by element, strand, and overlap", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  DBI::dbExecute(con, paste("UPDATE review_annotations SET tool = 'barrnap+ITSx+ITSx', model_coverage = 55
    WHERE element = '5.8S' AND", RC_UNIT))
  DBI::dbExecute(con, paste("UPDATE review_annotations SET model_coverage = 42 WHERE element = 'ITS1' AND", RC_UNIT))
  gff <- withr::local_tempfile()
  writeLines(c(
    paste("S1.1.1", "barrnap:0.9", "rRNA", 1, 400, "0", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA (partial);Note=aligned only 20 percent of the model", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 351, 2150, "0", "+", ".",
      "Name=18S_rRNA;product=18S ribosomal RNA (partial);Note=aligned only 91 percent of the model", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 2651, 2807, "0", "+", ".",
      "Name=5_8S_rRNA;product=5.8S ribosomal RNA", sep = "\t"),
    paste("S1.1.1", "barrnap:0.9", "rRNA", 3991, 6000, "0", "-", ".",
      "Name=28S_rRNA;product=28S ribosomal RNA", sep = "\t")), gff)
  orig <- rc_cov(con, "annotations")
  RiboPilot:::recheck_model_coverage(con, "S1", 1L, 1L, NULL)
  expect_equal(rc_cov(con)[c("18S", "28S")], orig[c("18S", "28S")])
  RiboPilot:::recheck_model_coverage(con, "S1", 1L, 1L, gff)
  cv <- rc_cov(con)
  expect_equal(unname(cv[c("18S", "5.8S", "28S", "ITS1")]), c(91, 100, 0, 42))
  expect_true(is.na(cv[["ITS2"]]))
  expect_equal(rc_cov(con, "annotations"), orig)
})

test_that("recheck_unit reruns every check on the copy and marks it done", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  # running, as query_recheck leaves it, with no counts yet
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 4, n_features = NULL, missing = NULL
    WHERE", RC_UNIT))
  # the original's best reference: its 18S differs at every tenth base
  ref <- strsplit(RC_SEQ, "")[[1]]
  i <- seq(360L, 2150L, by = 10L)
  ref[i] <- chartr("ACGT", "CGTA", ref[i])
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates (ID, path, scaffold, rank, accession)
    VALUES ('S1', 1, 1, 1, 'REF1.1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence, length) VALUES ('REF1.1', ?, ?)",
    params = list(paste(ref, collapse = ""), length(ref)))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations (accession, element, start, stop, strand)
    VALUES ('REF1.1', '18S', 351, 2150, '+')")
  before <- rc_originals(con)
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "done")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 2L)
  expect_true(is.na(u$check_note))
  expect_gte(u$checked_at, t0)
  expect_gte(u$time_stamp, t0)
  expect_equal(u$n_features, 5L)
  expect_equal(u$missing, "")
  expect_equal(u$chimeric, 1L)
  expect_false(is.na(u$n_flags))
  expect_equal(u$itsx_warnings, "")
  expect_equal(u$blast_accession, "A1.1")
  expect_equal(rc_cov(con)[c("18S", "28S")], c(`18S` = 100, `28S` = 70))
  ra <- DBI::dbGetQuery(con, paste("SELECT element, pid FROM review_ref_alignment WHERE", RC_UNIT))
  expect_equal(ra$element, "18S")
  expect_lt(ra$pid, 95)
  # curation saw the new model coverage and reference alignment, and the
  # chimera check the new ITSx flag
  fl <- DBI::dbGetQuery(con, paste("SELECT element, flags FROM review_annotations WHERE", RC_UNIT))
  expect_equal(fl$flags[fl$element == "28S"], "truncated")
  expect_match(fl$flags[fl$element == "18S"], "low_identity", fixed = TRUE)
  expect_match(DBI::dbGetQuery(con, paste("SELECT reason FROM review_chimera WHERE", RC_UNIT))$reason,
    "ITSx chimeric", fixed = TRUE)
  for (t in c("review_blast_hits", "review_blast_screen", "review_screen", "review_chimera")) {
    expect_gt(rc_count(con, t), 0, label = t)
  }
  # map mode has no panel: the per-feature BLAST is stored as not run
  fb <- DBI::dbGetQuery(con, paste("SELECT status, note FROM review_feature_blast WHERE", RC_UNIT))
  expect_equal(nrow(fb), 5L)
  expect_true(all(fb$status %in% c("not_run", "too_small")))
  expect_equal(rc_rev(con), rev)
  expect_equal(rc_originals(con), before)
})

test_that("failed tools leave their values as copied and the contig still finishes", {
  for (status in list(c(barrnap = 1, itsx = 1, itsx_check = 0), NULL)) {
    con <- rc_project()
    RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
    tools <- rc_tools(status = status, finished = FALSE)
    expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", tools), "done")
    u <- rc_copy(con)
    expect_equal(u$check_switch, 2L)
    expect_equal(u$chimeric, 0L)
    expect_equal(u$itsx_warnings, "not_run(ITSx did not finish)")
    expect_equal(rc_cov(con)[c("18S", "28S")], c(`18S` = 100, `28S` = 100))
  }
})

test_that("an ITSx check that throws is stored as not_run and the contig still finishes", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  local_mocked_bindings(itsx_check_one = function(...) stop("bad ITSx table"), .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", rc_tools()), "done")
  u <- DBI::dbGetQuery(con, paste("SELECT check_switch, itsx_warnings, itsx_origin, itsx_scores
    FROM review_units WHERE", RC_UNIT))
  expect_equal(u$check_switch, 2L)
  expect_equal(u$itsx_warnings, "not_run(ITSx check failed: bad ITSx table)")
  expect_true(is.na(u$itsx_origin))
  expect_true(is.na(u$itsx_scores))
})

test_that("an error in a step marks the copy failed, under the rev guard", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 4, time_stamp = 0 WHERE", RC_UNIT))
  local_mocked_bindings(screen_one = function(...) stop("screen broke"), .package = "RiboPilot")
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "failed")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 3L)
  expect_equal(u$check_note, "screen broke")
  expect_gte(u$time_stamp, t0)
  expect_equal(rc_count(con, "review_chimera"), 0)

  # an edit before the error: the failure is not recorded
  RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 1)
  rev <- rc_rev(con)
  local_mocked_bindings(screen_one = function(con, id, path, scaffold, ...) {
    RiboPilot:::review_trim(con, id, path, scaffold, 1L, 5900L)
    stop("screen broke")
  }, .package = "RiboPilot")
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  u <- rc_copy(con)
  expect_equal(u$check_switch, 1L)
  expect_true(is.na(u$check_note))
})

test_that("a copy that changed before the run starts is not rechecked", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  rev <- rc_rev(con)
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5900L)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_screen"), 0)
  expect_equal(rc_copy(con)$check_switch, 1L)

  rev <- rc_rev(con)
  RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 0)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_screen"), 0)
  RiboPilot:::review_drop(con, "S1", 1L, 1L)
  expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale")
  expect_equal(rc_count(con, "review_units"), 0)
  expect_equal(rc_count(con, "review_screen"), 0)
})

test_that("the rev guard discards a result after an edit, On hold, or Restore then a new trim", {
  edits <- list(
    edit = function(con) RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5900L),
    on_hold = function(con) RiboPilot:::review_set_state(con, data.frame(ID = "S1", path = 1L, scaffold = 1L), 0),
    restore_and_trim = function(con) {
      RiboPilot:::review_drop(con, "S1", 1L, 1L)
      Sys.sleep(0.01)
      RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 5800L)
    })
  for (nm in names(edits)) {
    con <- rc_project()
    RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
    rev <- rc_rev(con)
    # the edit lands while the last check runs
    local_mocked_bindings(chimera_unit = function(con, ...) edits[[nm]](con), .package = "RiboPilot")
    expect_equal(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rev, "map", rc_tools()), "stale", info = nm)
    u <- rc_copy(con)
    expect_equal(u$check_switch, if (nm == "on_hold") 0L else 1L, info = nm)
    expect_true(is.na(u$checked_at), info = nm)
    if (nm != "on_hold") expect_false(identical(rc_rev(con), rev), info = nm)
  }
})

test_that("recheck_release_running puts copies left running back to not run", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, check_switch,
        check_note, rev, time_stamp) VALUES ('S1', 1, ?, 'ACGT', 4, ?, ?, CAST(? AS INTEGER), 5)",
      params = list(s, c(4L, 4L, 1L, 2L)[s], c(NA, "old", NA, NA)[s], "1760000000000"))
  }
  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::recheck_release_running(con), 2L)
  r <- DBI::dbGetQuery(con, "SELECT check_switch, check_note, time_stamp FROM review_units ORDER BY scaffold")
  expect_equal(r$check_switch, c(1L, 1L, 1L, 2L))
  expect_equal(r$check_note, c("recheck did not finish", "recheck did not finish", NA, NA))
  expect_true(all(r$time_stamp[1:2] >= t0))
  expect_equal(r$time_stamp[3:4], c(5L, 5L))
  expect_equal(RiboPilot:::recheck_release_running(con), 0L)
})

test_that("NA dbarg reaches chimera_unit as NULL and a numeric rev errors", {
  con <- rc_project()
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1L, 6000L)
  seen <- "unset"
  local_mocked_bindings(chimera_unit = function(con, id, path, scaffold, dbarg = NULL, ...) seen <<- dbarg,
    .package = "RiboPilot")
  RiboPilot:::recheck_unit(con, "S1", 1L, 1L, rc_rev(con), "map", rc_tools(), dbarg = NA_character_)
  expect_null(seen)
  expect_error(RiboPilot:::recheck_unit(con, "S1", 1L, 1L, 1.76e12, "map", rc_tools()))
})
