# Each Update button shows one number, and it must agree with what its
# workflow's queries select, or the run window promises work the run will not
# do. The Assemble step's Update counts samples queued for assembly plus units
# queued for annotation only (main.nf); Annotation Review's counts the queued
# review copies (main_recheck.nf).

mk_pending_project <- function(dir) {
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, assemble_lock INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT,
    annotate_switch INT, annotate_lock INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxonomy_lock INT)")
  DBI::dbDisconnect(con)
  db
}

test_that("nf_pending_count adds assembly work to annotate-only work", {
  dir <- withr::local_tempdir()
  db <- mk_pending_project(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))

  expect_equal(RiboPilot:::nf_pending_count(dir), 0L)
  DBI::dbExecute(con, "INSERT INTO samples VALUES('A',1),('B',1)")

  # a sample queued for assembly
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('A',1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 1L)

  # a second sample, assembled already, with one unit queued for annotation
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('B',2,0)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('B',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('B',1,1,1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)

  # a queued unit on the sample that is about to be assembled is NOT extra work:
  # assembly re-seeds it, and the Nextflow query excludes it for that reason
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('A',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('A',1,1,1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)

  # a locked sample is not work
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'A'")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)  # A drops as assembly, returns as a unit
})

test_that("nf_pending_count skips samples whose taxonomy is not locked", {
  dir <- withr::local_tempdir()
  db <- mk_pending_project(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('A',1,0)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('A',0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 0L)
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  expect_equal(RiboPilot:::nf_pending_count(dir), 1L)
})

test_that("nf_pending_count takes no workflow argument", {
  # The main pipeline is one workflow; a caller passing "annotate" would
  # silently get a count for something that no longer exists.
  expect_named(formals(RiboPilot:::nf_pending_count), "path")
})

test_that("every nf_command() caller matches its signature", {
  # A real run caught run_assembly() still passing the removed workflow
  # argument, which the suite missed because nothing calls these drivers with a
  # project. Check the calls against the formals instead.
  fns <- list(run_assembly = RiboPilot::run_assembly,
    build_submit_script = RiboPilot:::build_submit_script,
    app_server = RiboPilot:::app_server)
  ok <- names(formals(RiboPilot:::nf_command))
  for (nm in names(fns)) {
    txt <- paste(deparse(fns[[nm]]), collapse = " ")
    calls <- regmatches(txt, gregexpr("nf_command\\([^)]*\\)", txt))[[1]]
    for (cl in calls) {
      args <- strsplit(sub("^nf_command\\(", "", sub("\\)$", "", cl)), ",")[[1]]
      args <- trimws(args)
      named <- grepl("=", args, fixed = TRUE)
      # only the project path may be positional. run_assembly() passed a second
      # positional argument, "assemble", which silently bound to `resume` once
      # the workflow parameter was removed, and the run died at launch.
      expect_lte(sum(!named), 1L, label = paste(nm, cl))
      for (a in sub("\\s*=.*", "", args[named])) {
        expect_true(a %in% ok, info = paste(nm, cl, a))
      }
    }
  }
})

# ---- Annotation Review's Update: the recheck queue ----

recheck_queue_ids <- c("QUEUED", "RUNNING", "LOCAL", "NO_BLAST", "DONE", "HELD", "FAILED",
  "LOCKED", "DROPPED", "ASM_OPEN", "TAX_OPEN", "NO_COPY")

# one contig per sample; only the first four are queued for a recheck
fill_recheck_queue <- function(con) {
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = CASE WHEN ID = 'TAX_OPEN' THEN 0 ELSE 1 END")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = CASE WHEN ID = 'ASM_OPEN' THEN 0 ELSE 1 END")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2,
    annotate_lock = CASE WHEN ID = 'LOCKED' THEN 1 ELSE 0 END")
  DBI::dbExecute(con, "UPDATE annotate_opts SET cpus = 3, memory = 5 WHERE annotate_opts = 'default'")
  DBI::dbExecute(con, "INSERT INTO blast_opts (blast_opts, run_blast) VALUES ('off', 0)")
  DBI::dbExecute(con, "UPDATE assemble SET blast_opts = 'off' WHERE ID = 'NO_BLAST'")
  for (id in recheck_queue_ids) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length, ignore)
      VALUES (?, 1, 1, 'ACGTACGT', 8, ?)", params = list(id, as.integer(id == "DROPPED")))
  }
  cs <- c(QUEUED = 1L, RUNNING = 4L, LOCAL = 1L, NO_BLAST = 1L, DONE = 2L, HELD = 0L,
    FAILED = 3L, LOCKED = 1L, DROPPED = 1L, ASM_OPEN = 1L, TAX_OPEN = 1L)
  for (i in seq_along(cs)) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length,
      trim_start, trim_end, edited, check_switch, rev, time_stamp)
      VALUES (?, 1, 1, 'ACGTACGT', 8, 1, 8, 1, ?, CAST(? AS INTEGER), 0)",
      params = list(names(cs)[i], cs[[i]], format(1791234567890 + i, scientific = FALSE)))
  }
  DBI::dbExecute(con, "INSERT INTO blast_searched (ID, path, scaffold, query_md5, source)
    VALUES ('QUEUED', 1, 1, 'x', 'remote'), ('LOCAL', 1, 1, 'x', 'local'),
           ('NO_BLAST', 1, 1, 'x', 'local')")
}

test_that("the recheck queue takes queued and running copies of locked, kept contigs", {
  con <- chimera_test_con(recheck_queue_ids)
  on.exit(DBI::dbDisconnect(con))
  fill_recheck_queue(con)
  got <- DBI::dbGetQuery(con, RiboPilot:::RECHECK_QUEUE_SQL)
  expect_named(got, c("ID", "path", "scaffold", "rev", "cpus", "mem", "mode"))
  got <- got[order(got$ID), ]
  expect_equal(got$ID, c("LOCAL", "NO_BLAST", "QUEUED", "RUNNING"))
  # rev as exact text: the workflow hands it back to the guard unchanged
  expect_identical(got$rev, c("1791234567893", "1791234567894", "1791234567891", "1791234567892"))
  # local only when the sample searches and the contig's original search was local
  expect_equal(got$mode, c("local", "map", "map", "map"))
  expect_equal(unique(got$cpus), 3L)
  expect_equal(unique(got$mem), 5L)
  expect_identical(RiboPilot:::nf_recheck_pending(dirname(DBI::dbGetInfo(con)$dbname)), 4L)
})

test_that("query_recheck runs on the sqlite3 CLI and marks what it takes as running", {
  skip_if(!nzchar(Sys.which("sqlite3")), "sqlite3 not installed")
  nf <- readLines(system.file("nextflow", "modules", "recheck.nf", package = "RiboPilot"))
  start <- grep("^process query_recheck \\{", nf)
  body <- nf[start:(start - 1 + which(nf[start:length(nf)] == "}")[1])]
  from <- grep("sqlite3 -cmd '.timeout 10000' -separator", body, fixed = TRUE)
  to <- grep("COMMIT;\" >> recheck_units.tsv", body, fixed = TRUE)
  expect_length(from, 1)
  expect_length(to, 1)
  sql <- sub("\" >> recheck_units.tsv$", "", paste(body[(from + 1):to], collapse = "\n"))

  con <- chimera_test_con(recheck_queue_ids)
  on.exit(DBI::dbDisconnect(con))
  fill_recheck_queue(con)
  db <- DBI::dbGetInfo(con)$dbname
  want <- DBI::dbGetQuery(con, RiboPilot:::RECHECK_QUEUE_SQL)
  want[] <- lapply(want, as.character)
  t0 <- as.integer(Sys.time())
  # the machine's sqlite3 CLI, which can be older than the SQLite inside RSQLite
  res <- processx::run("sqlite3", c("-cmd", ".timeout 10000", "-separator", "\t", db, sql))
  out <- utils::read.delim(text = res$stdout, header = FALSE, col.names = names(want),
    colClasses = "character")
  expect_equal(out[order(out$ID), ], want[order(want$ID), ], ignore_attr = TRUE)

  st <- DBI::dbGetQuery(con, "SELECT ID, check_switch, time_stamp FROM review_units ORDER BY ID")
  taken <- st$ID %in% want$ID
  expect_equal(st$check_switch[taken], rep(4L, 4))
  expect_true(all(st$time_stamp[taken] >= t0))
  expect_equal(stats::setNames(st$check_switch[!taken], st$ID[!taken]),
    c(ASM_OPEN = 1L, DONE = 2L, DROPPED = 1L, FAILED = 3L, HELD = 0L, LOCKED = 1L, TAX_OPEN = 1L))
  expect_true(all(st$time_stamp[!taken] == 0))
  # running copies stay in the queue, so the next Update takes them again
  expect_identical(RiboPilot:::nf_recheck_pending(dirname(db)), 4L)
})

test_that("nf_recheck_command and nf_recheck_pending take only the project", {
  expect_named(formals(RiboPilot:::nf_recheck_command), "path")
  expect_named(formals(RiboPilot:::nf_recheck_pending), "path")
  # the callers the app launch adds pass the project and nothing else
  ns <- asNamespace("RiboPilot")
  for (nm in c("app_server", "build_recheck_submit_script")) {
    fn <- get0(nm, envir = ns, inherits = FALSE)
    if (is.null(fn)) next
    txt <- paste(deparse(fn), collapse = " ")
    for (cl in regmatches(txt, gregexpr("nf_recheck_command\\([^)]*\\)", txt))[[1]]) {
      args <- trimws(strsplit(sub("^nf_recheck_command\\(", "", sub("\\)$", "", cl)), ",")[[1]])
      expect_length(args, 1)
      expect_false(grepl("=", sub("^path\\s*=", "", args), fixed = TRUE), label = paste(nm, cl))
    }
  }
})
