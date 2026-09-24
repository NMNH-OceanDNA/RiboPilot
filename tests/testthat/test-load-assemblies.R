# load_assemblies() replaces a sample's rows wholesale on every (re)assembly.
# The curator's keep/drop decisions must survive that, or changing assemble_opts
# (which re-queues the sample) silently un-ignores every dropped scaffold.

mk_asm_db <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT,
                       length INT, length_raw INT, sequence TEXT, source_contig TEXT, ignore INT, edited INT, time_stamp INT,
                       triage TEXT, triage_reason TEXT, triage_source TEXT, triage_rule TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, paths INT, scaffolds INT, length TEXT,
                       topology TEXT, assemble_switch INT, assemble_notes TEXT,
                       assemble_opts TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, max_paths INT,
                       max_scaffolds INT, min_assembly_length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
                       annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT,
                       partial TEXT, problematic TEXT, export_group TEXT,
                       PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assembly_backup(ID TEXT, path INT, scaffold INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default', 10, 10, NULL)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',NULL,NULL,NULL,NULL,1,NULL,'default',NULL)")
  con
}

# GetOrganelle-style output: one path fasta, two scaffolds.
write_asm_fasta <- function(dir, id = "S1") {
  ad <- file.path(dir, id, "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste0(">", id, ".1.1 circular scaffold_1--101+,102+"), strrep("ACGT", 200),
    paste0(">", id, ".1.2 linear scaffold_2--103+"), strrep("TGCA", 150)),
    file.path(ad, paste0(id, "_assembly_1.fasta")))
  invisible(ad)
}

test_that("a re-assembly preserves per-scaffold ignore flags", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  write_asm_fasta(dir)

  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assemblies")$n, 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate")$n, 2)

  # WF1 annotates every unit as it is assembled, so give scaffold 2 a feature
  # model before dropping it -- that is the work a destructive drop would lose.
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,2,'18S')")

  # curator drops scaffold 2
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE ID='S1' AND path=1 AND scaffold=2")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 1L)

  # HELD, not deleted: the unit is already invisible downstream (everything
  # filters on assemblies.ignore), so dropping its curation buys nothing and
  # makes "drop it, look again, put it back" cost a re-annotation.
  st <- DBI::dbGetQuery(con, "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(nrow(st), 1L)
  expect_equal(st$annotate_switch, 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotations WHERE scaffold=2")$n, 1)

  # re-run assembly over the same output
  RiboPilot:::load_assemblies(con, dir, "S1")
  ig <- DBI::dbGetQuery(con, "SELECT scaffold, ignore FROM assemblies WHERE ID='S1' ORDER BY scaffold")
  expect_equal(ig$ignore, c(0L, 1L))
  # a re-assembly must not treat a curator-dropped unit as an orphan
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate")$n, 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotations WHERE scaffold=2")$n, 1)
})

test_that("un-ignoring a held unit does not re-queue an already-annotated one", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  write_asm_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,2,'18S')")

  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 1L)   # drop
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 0L)   # put back
  st <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(st$annotate_switch, 2L)   # already annotated -> done, no WF2 round trip

  # a unit with no feature model still queues for annotation
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold=1")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 1L, 1L)
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 1L, 0L)
  st1 <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=1")
  expect_equal(st1$annotate_switch, 1L)
})

test_that("assemble.scaffolds records total units, not the largest path's count", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  ad <- write_asm_fasta(dir)
  # a second path with one scaffold: 3 units total, max 2 per path
  writeLines(c(">S1.2.1 circular scaffold_1--104+", strrep("ACGT", 180)),
    file.path(ad, "S1_assembly_2.fasta"))

  RiboPilot:::load_assemblies(con, dir, "S1")
  a <- DBI::dbGetQuery(con, "SELECT paths, scaffolds FROM assemble WHERE ID='S1'")
  expect_equal(a$paths, 2)
  expect_equal(a$scaffolds, 3)
})

# --- max_paths / max_scaffolds are gates, not caps -----------------------------
# Both were editable in the app and documented as thresholds "for a sample to
# continue", but max_scaffolds was read by nothing and max_paths quietly kept its
# top N instead of stopping. A sample that blew past both went on to cost one
# serial remote BLAST per scaffold.

write_many_scaffolds <- function(dir, n, id = "S1", path = 1) {
  ad <- file.path(dir, id, "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(unlist(lapply(seq_len(n), function(k) {
    c(sprintf(">%s.%d.%d linear scaffold_%d", id, path, k, k), strrep("ACGT", 100 + k))
  })), file.path(ad, sprintf("%s_assembly_%d.fasta", id, path)))
}

gate_state <- function(con) {
  DBI::dbGetQuery(con, "SELECT assemble_switch, assemble_notes, scaffolds FROM assemble WHERE ID='S1'")
}

test_that("too many scaffolds in a path fails the sample", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 3")
  write_many_scaffolds(dir, 5)

  RiboPilot:::load_assemblies(con, dir, "S1")

  st <- gate_state(con)
  expect_equal(st$assemble_switch, 3L)
  expect_match(st$assemble_notes, "max_scaffolds")
})

test_that("a failed gate stops the downstream work, but keeps the sequences", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 3")
  write_many_scaffolds(dir, 5)

  RiboPilot:::load_assemblies(con, dir, "S1")

  # kept, so raising the limit and re-running does not mean re-assembling
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assemblies")$n, 5)
  # every downstream step selects on ignore = 0, so this is what stops them
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemblies WHERE COALESCE(ignore,0) = 0")$n, 0)
  # and nothing is queued for annotation
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM annotate WHERE annotate_switch = 1")$n, 0)
})

test_that("too many paths fails the sample instead of keeping the longest", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_paths = 2, max_scaffolds = 10")
  for (p in 1:3) write_many_scaffolds(dir, 1, path = p)

  RiboPilot:::load_assemblies(con, dir, "S1")

  st <- gate_state(con)
  expect_equal(st$assemble_switch, 3L)
  expect_match(st$assemble_notes, "max_paths")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(DISTINCT path) n FROM assemblies")$n, 3)
})

test_that("a sample inside both limits is untouched by the gate", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_paths = 10, max_scaffolds = 10")
  write_many_scaffolds(dir, 4)

  RiboPilot:::load_assemblies(con, dir, "S1")

  st <- gate_state(con)
  expect_equal(st$assemble_switch, 4L)          # in progress, not failed
  expect_false(grepl("max_", st$assemble_notes))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM annotate WHERE annotate_switch = 1")$n, 4)
})

test_that("the limit counts scaffolds per path, not across the whole sample", {
  # 3 paths of 2 scaffolds each is 6 units but no path is over a limit of 3
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_paths = 5, max_scaffolds = 3")
  for (p in 1:3) write_many_scaffolds(dir, 2, path = p)

  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(gate_state(con)$assemble_switch, 4L)
})

test_that("an unset or zero limit means no gate", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_paths = 0, max_scaffolds = NULL")
  write_many_scaffolds(dir, 30)

  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(gate_state(con)$assemble_switch, 4L)
})

test_that("raising the limit and re-running recovers a gated sample", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 3")
  write_many_scaffolds(dir, 5)
  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(gate_state(con)$assemble_switch, 3L)

  # the curator raises the limit and re-queues; the blanket ignore the gate set
  # must not be mistaken for their own keep/drop decisions
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 10")
  RiboPilot:::load_assemblies(con, dir, "S1")

  expect_equal(gate_state(con)$assemble_switch, 4L)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemblies WHERE COALESCE(ignore,0) = 0")$n, 5)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM annotate WHERE annotate_switch = 1")$n, 5)
})

test_that("a curator's own drops still survive a re-assembly", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 10")
  write_many_scaffolds(dir, 4)
  RiboPilot:::load_assemblies(con, dir, "S1")

  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE scaffold = 2")
  RiboPilot:::load_assemblies(con, dir, "S1")

  ig <- DBI::dbGetQuery(con, "SELECT scaffold, ignore FROM assemblies ORDER BY scaffold")
  expect_equal(ig$ignore, c(0L, 1L, 0L, 0L))
})

# --- min_assembly_length is a filter, not a gate -------------------------------
# Short contigs are noise. They are parked rather than discarded (the sequence
# stays available if the curator lowers the limit) and counted out of the gates,
# so junk fragments cannot fail an otherwise good sample.
# Lengths written by write_many_scaffolds(): 4 * (100 + k), i.e. 404 .. 420.

test_that("scaffolds under min_assembly_length are parked, not dropped", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET min_assembly_length = 415")
  write_many_scaffolds(dir, 5)

  RiboPilot:::load_assemblies(con, dir, "S1")

  ig <- DBI::dbGetQuery(con,
    "SELECT scaffold, length, ignore, triage, triage_source, triage_rule, triage_reason
       FROM assemblies ORDER BY scaffold")
  expect_equal(nrow(ig), 5)                       # sequences kept
  expect_equal(ig$ignore, c(1L, 1L, 1L, 0L, 0L))  # three short ones parked
  # recorded as a triage decision, with a reason the curator can read, because
  # triage owns assemblies.ignore and would otherwise overturn a raw flag
  expect_equal(ig$triage, c("drop", "drop", "drop", NA, NA))
  expect_equal(ig$triage_rule[1], "too_short")
  expect_equal(ig$triage_source[1], "auto")
  expect_match(ig$triage_reason[1], "min_assembly_length")
  st <- gate_state(con)
  expect_equal(st$assemble_switch, 4L)
  expect_match(st$assemble_notes, "min_assembly_length")
  # the dropped units keep a held annotate row rather than a queued one, so they
  # are not pending work and the sample can settle
  sw <- DBI::dbGetQuery(con, "SELECT scaffold, annotate_switch FROM annotate ORDER BY scaffold")
  expect_equal(sw$annotate_switch, c(0L, 0L, 0L, 1L, 1L))
  # the summary row describes the scaffolds the sample is actually using, on one
  # basis: counts and the length list must not disagree about which those are
  expect_equal(st$scaffolds, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT length FROM assemble WHERE ID='S1'")$length,
    "420;416")
})

test_that("short scaffolds do not count towards max_scaffolds", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET max_scaffolds = 3, min_assembly_length = 415")
  write_many_scaffolds(dir, 5)   # 5 scaffolds would trip a limit of 3

  RiboPilot:::load_assemblies(con, dir, "S1")

  st <- gate_state(con)
  expect_equal(st$assemble_switch, 4L)
  expect_false(grepl("max_scaffolds", st$assemble_notes))
})

test_that("a sample whose every scaffold is too short fails", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble_opts SET min_assembly_length = 5000")
  write_many_scaffolds(dir, 3)

  RiboPilot:::load_assemblies(con, dir, "S1")

  st <- gate_state(con)
  expect_equal(st$assemble_switch, 3L)
  expect_match(st$assemble_notes, "min_assembly_length")
})

# --- header parsing ------------------------------------------------------------
# `>ID.path.scaffold` is read from the right: an ID with a dot in it would
# otherwise shift the keys, and a non-numeric field would key rows on NA, which
# every downstream join silently misses.

test_that("a dot in the sample ID does not shift the path and scaffold", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET ID = 'SRR123.1'")
  ad <- file.path(dir, "SRR123.1", "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(">SRR123.1.2.3 linear scaffold_1", strrep("ACGT", 100)),
    file.path(ad, "SRR123.1_assembly_1.fasta"))

  RiboPilot:::load_assemblies(con, dir, "SRR123.1")

  a <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM assemblies")
  expect_equal(a$ID, "SRR123.1")
  expect_equal(a$path, 2L)
  expect_equal(a$scaffold, 3L)
})

test_that("a header with a non-numeric path or scaffold is skipped with a warning", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  ad <- file.path(dir, "S1", "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    ">S1.v2.1 linear scaffold_bad", strrep("ACGT", 100),
    ">S1.1.x linear scaffold_bad", strrep("ACGT", 100),
    ">S1.1.1 linear scaffold_ok", strrep("ACGT", 100)),
    file.path(ad, "S1_assembly_1.fasta"))

  expect_warning(RiboPilot:::load_assemblies(con, dir, "S1"), "skipped 2")

  a <- DBI::dbGetQuery(con, "SELECT path, scaffold FROM assemblies")
  expect_equal(nrow(a), 1L)
  expect_equal(a$path, 1L)
  expect_equal(a$scaffold, 1L)
})

test_that("a fasta of nothing but unparsable headers fails the sample", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  ad <- file.path(dir, "S1", "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(">S1 linear scaffold_1", strrep("ACGT", 100)),
    file.path(ad, "S1_assembly_1.fasta"))

  expect_warning(RiboPilot:::load_assemblies(con, dir, "S1"), "skipped 1")

  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assemblies")$n, 0)
  expect_equal(gate_state(con)$assemble_switch, 3L)
})

test_that("the length filter wins over the previous run's ignore flags", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  write_many_scaffolds(dir, 5)
  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemblies WHERE COALESCE(ignore,0) = 0")$n, 5)

  # the curator raises the floor and re-runs: carrying the old ignore = 0 over
  # would leave the short scaffolds live
  DBI::dbExecute(con, "UPDATE assemble_opts SET min_assembly_length = 415")
  RiboPilot:::load_assemblies(con, dir, "S1")

  ig <- DBI::dbGetQuery(con, "SELECT ignore FROM assemblies ORDER BY scaffold")$ignore
  expect_equal(ig, c(1L, 1L, 1L, 0L, 0L))
})

test_that("the contig name from the source fasta is stored", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  write_asm_fasta(dir)

  RiboPilot:::load_assemblies(con, dir, "S1")
  src <- DBI::dbGetQuery(con,
    "SELECT source_contig FROM assemblies ORDER BY path, scaffold")$source_contig
  expect_equal(src, c("scaffold_1--101+,102+", "scaffold_2--103+"))
})

test_that("a header with no contig name leaves it empty", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  ad <- file.path(dir, "S1", "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(">S1.1.1 linear", strrep("ACGT", 200)),
    file.path(ad, "S1_assembly_1.fasta"))

  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT source_contig FROM assemblies")$source_contig))
})

test_that("a re-assembly that changes a unit's sequence clears its old evidence and lock", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  write_asm_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")

  # both units annotated, screened, locked, reviewed, and grouped for export
  for (s in 1:2) {
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO screen (ID,path,scaffold,screen_verdict) VALUES('S1',1,?,'on_target')",
      params = list(s))
  }
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, annotate_lock = 1, reviewed = 'yes',
    partial = 'yes', export_group = 'g1'")

  # scaffold 1 comes back as a different contig; scaffold 2 is unchanged
  fa <- file.path(dir, "S1", "assemble", "S1_assembly_1.fasta")
  writeLines(c(">S1.1.1 linear scaffold_1--101+", strrep("GATC", 120),
    ">S1.1.2 linear scaffold_2--103+", strrep("TGCA", 150)), fa)
  RiboPilot:::load_assemblies(con, dir, "S1")

  ev <- function(tbl, s) DBI::dbGetQuery(con,
    sprintf("SELECT COUNT(*) n FROM %s WHERE ID='S1' AND path=1 AND scaffold=%d", tbl, s))$n
  an <- DBI::dbGetQuery(con, "SELECT * FROM annotate WHERE ID='S1' ORDER BY scaffold")
  expect_equal(ev("annotations", 1), 0)
  expect_equal(ev("screen", 1), 0)
  expect_equal(an$annotate_lock, c(0L, 1L))
  expect_equal(an$annotate_switch, c(1L, 2L))
  expect_equal(an$reviewed, c("no", "yes"))
  expect_equal(an$export_group, c(NA, "g1"))
  # the unchanged unit keeps everything
  expect_equal(ev("annotations", 2), 1)
  expect_equal(ev("screen", 2), 1)
})

test_that("a failed assembly carries the assembler's reason into the note", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir); on.exit(DBI::dbDisconnect(con))
  ad <- file.path(dir, "S1", "assemble")
  dir.create(ad, recursive = TRUE)
  writeLines(">No assembly found: SPAdes exit 21: Invalid kmer coverage histogram",
    file.path(ad, "S1_assembly_0.fasta"))
  RiboPilot::load_assemblies(con, dir, "S1")
  r <- DBI::dbGetQuery(con, "SELECT assemble_switch, assemble_notes FROM assemble WHERE ID = 'S1'")
  expect_equal(r$assemble_switch, 3L)
  expect_equal(r$assemble_notes, "failed assembly: SPAdes exit 21: Invalid kmer coverage histogram")
  # the supplied-assembly marker says nothing new; its reason is in find_rrna_notes
  writeLines(">No rDNA contigs found", file.path(ad, "S1_assembly_0.fasta"))
  RiboPilot::load_assemblies(con, dir, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT assemble_notes FROM assemble")$assemble_notes,
    "failed assembly")
})
