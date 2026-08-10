# Schema migrations. v9 = retained hit distribution + rank-aligned lineages,
# v10 = median per-position depth.
# Migrations must be additive and idempotent, and must not disturb existing rows.

test_that("the schema version stamp is transactional", {
  # migrate_db() stamps user_version inside its transaction. That is only safe
  # because the PRAGMA rolls back with the transaction rather than being written
  # straight to the file header.
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "PRAGMA user_version = 7")
  DBI::dbExecute(con, "BEGIN IMMEDIATE")
  DBI::dbExecute(con, "PRAGMA user_version = 99")
  DBI::dbExecute(con, "ROLLBACK")
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 7)
})

# A minimal v8-shaped project DB, as a project created by the previous release
# would look.
mk_v8_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL,
    ignore INT, PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT,
    title TEXT, species TEXT, staxid INT, lineage TEXT, pident REAL, qcovs REAL, evalue REAL,
    is_nuclear_rrna INT, time_stamp INT, PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, time_stamp INT,
    PRIMARY KEY (ID, path, scaffold, rank))")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa; Chordata')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'ACGT','AA1','Danio rerio',99,80,1e-90,0)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',2,0,'default','default')")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,'AA1','t','Danio rerio',7955,'Metazoa',99,80,1e-90,1,1)")
  DBI::dbExecute(con, "PRAGMA user_version = 8")
  con
}

test_that("migrating a v8 project adds the v9 tables and columns", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("blast_hits" %in% DBI::dbListTables(con))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("blast_hits" %in% DBI::dbListTables(con))
  expect_true("blast_source" %in% DBI::dbListFields(con, "assemblies"))
  expect_true("lineage_taxids" %in% DBI::dbListFields(con, "samples"))
  for (f in c("lineage_taxids", "rank_distance", "bitscore")) {
    expect_true(f %in% DBI::dbListFields(con, "blast_screen"), info = f)
  }
  expect_true("bitscore" %in% DBI::dbListFields(con, "blast_ref_candidates"))
})

test_that("migration preserves existing rows and leaves new columns NULL", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  s <- DBI::dbGetQuery(con, "SELECT * FROM samples")
  expect_equal(s$ID, "S1")
  expect_equal(s$taxid, 7955L)
  expect_equal(s$lineage, "Metazoa; Chordata")
  expect_true(is.na(s$lineage_taxids))

  a <- DBI::dbGetQuery(con, "SELECT * FROM assemblies")
  expect_equal(a$blast_accession, "AA1")
  expect_equal(a$sequence, "ACGT")
  expect_true(is.na(a$blast_source))   # provenance unknown for pre-v9 rows

  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 1)
})

test_that("migrate_db is idempotent", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession)
                       VALUES('S1',1,1,1,1,'AA1')")
  expect_silent(migrate_db(con))
  expect_silent(migrate_db(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  # a second migrate must not wipe rows written after the first
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_hits")$n, 1)
})

test_that("a fresh new_db() matches a migrated one", {
  dir <- withr::local_tempdir()
  map <- file.path(dir, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "a_R1.fq", R2 = "a_R2.fq"),
    map, row.names = FALSE)
  db <- file.path(dir, ".sqlite")
  new_db(db_path = db, mapping_fn = map,
    seeds_db = map, labels_db = map)   # any existing file passes the DB check
  fresh <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(fresh), add = TRUE)

  expect_equal(DBI::dbGetQuery(fresh, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("blast_hits" %in% DBI::dbListTables(fresh))
  expect_true("blast_source" %in% DBI::dbListFields(fresh, "assemblies"))
  expect_true("lineage_taxids" %in% DBI::dbListFields(fresh, "samples"))

  mig <- mk_v8_db(); on.exit(DBI::dbDisconnect(mig), add = TRUE)
  migrate_db(mig)
  # the columns v9 introduces must exist in both, or code diverges by vintage
  for (tbl in c("blast_screen", "blast_ref_candidates", "blast_hits", "screen")) {
    expect_true(all(DBI::dbListFields(mig, tbl) %in% DBI::dbListFields(fresh, tbl)), info = tbl)
  }
})

test_that("v10 adds median_depth and back-fills it from the stored tracks", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    min_depth INT, max_depth INT, depths TEXT, time_stamp INT)")
  # a track with a heavy right tail: mean is pulled well above the typical base
  spiky <- paste(c(rep(10, 90), rep(1000, 10)), collapse = ",")
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,1,109,10,1000,?,1)", params = list(spiky))
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,2,0,0,0,'',1)")
  DBI::dbExecute(con, "PRAGMA user_version = 9")

  migrate_db(con)
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("median_depth" %in% DBI::dbListFields(con, "coverage"))
  # migration alone leaves it NULL: parsing every track would stall app startup
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT median_depth FROM coverage WHERE scaffold=1")$median_depth))

  expect_equal(backfill_coverage_median(con), 1L)   # the empty track is skipped
  got <- DBI::dbGetQuery(con, "SELECT mean_depth, median_depth FROM coverage WHERE scaffold=1")
  expect_equal(got$median_depth, 10)                # typical base, not 109
  expect_equal(got$mean_depth, 109)
  expect_equal(backfill_coverage_median(con), 0L)   # idempotent
})

test_that("write_coverage stores both statistics", {
  skip_if_not_installed("withr")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    median_depth REAL, min_depth INT, max_depth INT, depths TEXT, time_stamp INT)")
  f <- withr::local_tempfile()
  writeLines(sprintf("ref\t%d\t%d", 1:100, c(rep(10, 90), rep(1000, 10))), f)
  RiboPilot::write_coverage(con, "S1", 1, 1, f)
  d <- DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, max_depth FROM coverage")
  expect_equal(d$median_depth, 10)
  expect_equal(d$mean_depth, 109)
  expect_equal(d$max_depth, 1000L)
})

test_that("v12 adds triage_opts, its seed row, and the rule columns", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  expect_true("triage_opts" %in% DBI::dbListTables(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM triage_opts")$n, 1L)

  # every rule key has a column, and the seed carries the shipped defaults
  r <- DBI::dbGetQuery(con, "SELECT * FROM triage_opts WHERE triage_opts = 'default'")
  for (k in RiboPilot:::triage_rule_keys()) {
    expect_true(paste0("rule_", k) %in% names(r))
  }
  expect_equal(r$rule_no_feat_non_rdna, "drop")
  expect_equal(r$rule_odd_one_out, "flag")

  expect_true("triage_opts" %in% DBI::dbListFields(con, "assemble"))
  expect_true("triage_rule" %in% DBI::dbListFields(con, "assemblies"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
})

test_that("v12 backfills coverage.depths when a coverage table has path but not depths", {
  # Not a real migration path -- init_db.R and the pre-v8 coverage rebuild
  # always add path and depths together -- but unit_evidence() now hard-selects
  # depths, so its absence would break the whole scaffold query, not one cell.
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
  expect_false("depths" %in% DBI::dbListFields(con, "coverage"))

  migrate_db(con)

  expect_true("depths" %in% DBI::dbListFields(con, "coverage"))
})

test_that("v13 adds hit_class and backfills it from stored titles", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  # a v12-shaped blast_hits, i.e. with every column EXCEPT hit_class
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, hsp INT, accession TEXT, title TEXT, species TEXT, staxid INT,
    lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, aln_len INT, qstart INT, qend INT, sstart INT,
    send INT, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT,
    is_nuclear_rrna INT)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,2,0)")
  ins <- "INSERT INTO blast_hits (ID,path,scaffold,hit_rank,hsp,title,bitscore,is_nuclear_rrna)
          VALUES (?,?,?,?,1,?,?,?)"
  # unit 1: all chromosome records, the case that used to score 0
  for (i in 1:3) DBI::dbExecute(con, ins, params = list("S1", 1L, 1L, i,
    sprintf("Scomber japonicus DNA, primary haplotype, chromosome %d, sequence", i), 900 - i, 0L))
  # unit 2: genuine rDNA
  for (i in 1:3) DBI::dbExecute(con, ins, params = list("S1", 1L, 2L, i,
    "Danio rerio 18S ribosomal RNA gene, complete sequence", 900 - i, 1L))
  DBI::dbExecute(con, "PRAGMA user_version = 12")

  migrate_db(con)

  expect_true("hit_class" %in% DBI::dbListFields(con, "blast_hits"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  cls <- DBI::dbGetQuery(con, "SELECT scaffold, hit_class FROM blast_hits ORDER BY scaffold, hit_rank")
  expect_equal(unique(cls$hit_class[cls$scaffold == 1]), "uninformative")
  expect_equal(unique(cls$hit_class[cls$scaffold == 2]), "rdna")

  # the repair: unit 1 goes from a false 0 to NA, unit 2 to 1
  scr <- DBI::dbGetQuery(con, "SELECT scaffold, is_nuclear_rrna FROM blast_screen ORDER BY scaffold")
  expect_true(is.na(scr$is_nuclear_rrna[scr$scaffold == 1]))
  expect_equal(scr$is_nuclear_rrna[scr$scaffold == 2], 1L)
})

test_that("the backfill leaves units with no hits alone and is idempotent", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, title TEXT, bitscore REAL, hit_class TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT,
    is_nuclear_rrna INT)")
  # a screened unit with NO hits rows must keep whatever it had
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,9,1)")
  RiboPilot:::backfill_hit_class(con)
  expect_equal(DBI::dbGetQuery(con,
    "SELECT is_nuclear_rrna FROM blast_screen WHERE scaffold=9")$is_nuclear_rrna, 1L)

  DBI::dbExecute(con, "INSERT INTO blast_hits VALUES('S1',1,1,1,
    'Danio rerio 18S ribosomal RNA gene, complete sequence', 900, NULL)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,0)")
  n1 <- RiboPilot:::backfill_hit_class(con)
  expect_equal(n1, 1L)
  n2 <- RiboPilot:::backfill_hit_class(con)
  expect_equal(n2, 0L)   # nothing left to classify
  expect_equal(DBI::dbGetQuery(con,
    "SELECT is_nuclear_rrna FROM blast_screen WHERE scaffold=1")$is_nuclear_rrna, 1L)
})

test_that("migrate_db does not rescan blast_hits when already at v13", {
  # parse_blast_candidates() calls migrate_db() once per unit; an ungated
  # backfill would rescan the whole project every time, so a DB already at v13
  # must come out of migrate_db() with its NULL hit_class rows untouched.
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, title TEXT, bitscore REAL, hit_class TEXT)")
  DBI::dbExecute(con, "INSERT INTO blast_hits VALUES('S1',1,1,1,
    'Danio rerio 18S ribosomal RNA gene, complete sequence', 900, NULL)")
  DBI::dbExecute(con, "PRAGMA user_version = 13")

  migrate_db(con)

  cls <- DBI::dbGetQuery(con, "SELECT hit_class FROM blast_hits WHERE scaffold=1")$hit_class
  expect_true(is.na(cls))
})

test_that("v14 adds the screen table and its companions", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  expect_true("screen" %in% DBI::dbListTables(con))
  for (col in c("screen_verdict", "screen_class", "screen_confidence",
                "screen_reason", "screen_version", "lca_taxid", "lca_taxon",
                "lca_rank", "lca_n_hits", "host_concordant", "margin")) {
    expect_true(col %in% DBI::dbListFields(con, "screen"), info = col)
  }
  expect_true("reference_support" %in% DBI::dbListFields(con, "samples"))
  expect_true("rule_screen_off_target" %in% DBI::dbListFields(con, "triage_opts"))
  # migrate_db() always stamps the latest, so this reads the constant rather
  # than 14L; the v15 test below pins the literal for the current version
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)

  # the new rule ships as flag, per doc 08 section 5's validation gate
  expect_equal(DBI::dbGetQuery(con,
    "SELECT rule_screen_off_target v FROM triage_opts WHERE triage_opts='default'")$v,
    "flag")
})

test_that("the screen table is created empty rather than backfilled", {
  # scoring needs samples.lineage_taxids and a taxdump, neither of which a
  # migration can assume; an empty row must stay distinguishable from a computed one
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen")$n, 0L)
})

test_that("v15 adds the intragenomic variant columns to screen", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  for (col in c("variant_class", "variant_group", "variant_identity")) {
    expect_true(col %in% DBI::dbListFields(con, "screen"), info = col)
  }
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
})

test_that("a real v14 project gains the variant columns and keeps its verdict", {
  # The screen table must be built at the V14 SHAPE, not from the current
  # SCREEN_DDL: a v8 db migrates by creating the table from that DDL, which
  # already has the new columns, so add_col() is a no-op there and the path an
  # existing v14 project actually takes goes untested.
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE screen (
    ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
    screen_verdict TEXT, screen_class TEXT, screen_confidence REAL, screen_reason TEXT,
    screen_version INTEGER, lca_taxid INTEGER, lca_taxon TEXT, lca_rank TEXT,
    lca_n_hits INTEGER, host_concordant INTEGER, margin REAL, time_stamp INTEGER,
    PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con,
    "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version)
     VALUES ('S1',1,1,'on_target',2)")
  DBI::dbExecute(con, "PRAGMA user_version = 14")
  expect_false("variant_class" %in% DBI::dbListFields(con, "screen"))

  migrate_db(con)
  migrate_db(con)  # additive and idempotent: the second run is a no-op
  for (col in c("variant_class", "variant_group", "variant_identity")) {
    expect_true(col %in% DBI::dbListFields(con, "screen"), info = col)
  }
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen")
  expect_equal(nrow(r), 1L)
  expect_equal(r$screen_verdict, "on_target")
  expect_true(is.na(r$variant_class))
})

# --- v16: per-sample assembler ---

# assemble_opts at its v15 shape. mk_v8_db() has no assemble_opts at all, and
# add_col() guards on the table existing, so without this the migration path an
# existing project actually takes goes untested.
mk_v15_assemble_opts <- function(con) {
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (
    assemble_opts TEXT NOT NULL, cpus INTEGER, memory INTEGER, getOrganelle TEXT,
    seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INTEGER,
    max_scaffolds INTEGER, min_assembly_length INTEGER,
    PRIMARY KEY (assemble_opts))")
  DBI::dbExecute(con,
    "INSERT INTO assemble_opts (assemble_opts, cpus, getOrganelle)
     VALUES ('default', 4, '-F anonym')")
}

test_that("v16 adds assembler to assemble_opts and defaults existing presets", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  DBI::dbExecute(con, "PRAGMA user_version = 15")
  expect_false("assembler" %in% DBI::dbListFields(con, "assemble_opts"))

  migrate_db(con)
  migrate_db(con)  # additive and idempotent

  expect_true("assembler" %in% DBI::dbListFields(con, "assemble_opts"))
  r <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts")
  expect_equal(nrow(r), 1L)
  expect_equal(r$assembler, "getorganelle")   # backfilled by the column DEFAULT
  expect_equal(r$getOrganelle, "-F anonym")   # untouched
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
})

test_that("v16 returns fallback-stranded samples to the ordinary queue", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S5',5,0,'default','default')")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S4',4,0,'default','default')")
  DBI::dbExecute(con, "PRAGMA user_version = 15")

  migrate_db(con)

  sw <- DBI::dbGetQuery(con, "SELECT ID, assemble_switch FROM assemble ORDER BY ID")
  expect_equal(sw$assemble_switch[sw$ID == "S5"], 1L)  # was stranded at 5
  expect_equal(sw$assemble_switch[sw$ID == "S4"], 4L)  # other states untouched
  expect_equal(sw$assemble_switch[sw$ID == "S1"], 2L)
})

# --- v17: per-preset options for the fallback assembler ---

test_that("v17 adds the bowtie2 and spades option columns with working defaults", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  DBI::dbExecute(con, "PRAGMA user_version = 15")

  migrate_db(con)
  migrate_db(con)  # additive and idempotent

  r <- DBI::dbGetQuery(con, "SELECT * FROM assemble_opts")
  expect_equal(r$bowtie2, "--local")     # the flags the fallback used hardcoded
  expect_equal(r$spades, "--isolate")
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
})

test_that("concurrent migrations of the same project do not collide", {
  # Every pipeline process migrates on open, and WF2 fans out to one process per
  # unit, so a pending migration is attempted by a dozen processes at once. With
  # an unserialised migration two of them both see a column missing and both
  # ALTER, and the loser dies with "duplicate column name", failing the run.
  skip_on_os("windows")
  skip_on_cran()
  # the race needs real contention, and R CMD check caps a package at 2 cores
  skip_if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")), "R CMD check limits cores")
  db <- file.path(withr::local_tempdir(), "race.sqlite")
  seed <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(seed, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  mk_v15_assemble_opts(seed)
  DBI::dbExecute(seed, "PRAGMA user_version = 15")
  DBI::dbDisconnect(seed)

  start <- as.numeric(Sys.time()) + 2   # common start instant, so they collide
  res <- parallel::mclapply(1:12, function(i) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "PRAGMA busy_timeout=15000;")
    while (as.numeric(Sys.time()) < start) NULL
    tryCatch({ migrate_db(con); "ok" }, error = function(e) conditionMessage(e))
  }, mc.cores = 12)

  expect_equal(unique(unlist(res)), "ok")
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))
  expect_true(all(c("assembler", "bowtie2", "spades") %in%
    DBI::dbListFields(con, "assemble_opts")))
})

test_that("a fresh database and a migrated one agree on assemble_opts columns", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  DBI::dbExecute(con, "PRAGMA user_version = 15")
  migrate_db(con)

  dir <- withr::local_tempdir()
  map <- file.path(dir, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "a_R1.fq", R2 = "a_R2.fq"),
    map, row.names = FALSE)
  db <- file.path(dir, ".sqlite")
  new_db(db_path = db, mapping_fn = map, seeds_db = map, labels_db = map)
  fresh <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(fresh), add = TRUE)

  expect_true("assembler" %in% DBI::dbListFields(fresh, "assemble_opts"))
  expect_equal(DBI::dbGetQuery(fresh, "SELECT assembler FROM assemble_opts")$assembler,
    "getorganelle")
  expect_setequal(DBI::dbListFields(fresh, "assemble_opts"),
    DBI::dbListFields(con, "assemble_opts"))
})

test_that("v18 drops the stitched table, in a migrated and a fresh database", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE stitched(ID TEXT PRIMARY KEY, path INT, sequence TEXT,
    length INT, reference_accession TEXT, n_gaps INT, n_redundant INT, time_stamp INT)")
  DBI::dbExecute(con, "INSERT INTO stitched VALUES('S1',1,'ACGT',4,'AA1',0,0,1)")

  migrate_db(con)
  expect_false("stitched" %in% DBI::dbListTables(con))
  expect_silent(migrate_db(con))

  dir <- withr::local_tempdir()
  map <- file.path(dir, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "a_R1.fq", R2 = "a_R2.fq"),
    map, row.names = FALSE)
  db <- file.path(dir, ".sqlite")
  new_db(db_path = db, mapping_fn = map, seeds_db = map, labels_db = map)
  fresh <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(fresh), add = TRUE)
  expect_false("stitched" %in% DBI::dbListTables(fresh))
})

# --- serialization: busy is not the same thing as already-in-a-transaction ---

# A file DB is required for any of this: two connections to ":memory:" are two
# different databases, so they can never contend.
mk_locked_db <- function(env = parent.frame()) {
  db <- file.path(withr::local_tempdir(.local_envir = env), "lock.sqlite")
  seed <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(seed, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(seed, "PRAGMA user_version = 15")
  DBI::dbDisconnect(seed)
  db
}

test_that("db_in_transaction tells an open transaction from a busy database", {
  # The whole point of the guard: a write lock held by ANOTHER process must not
  # read as "the caller already gave me a transaction", or migrate_db() would
  # run every guarded ALTER in autocommit, unserialized.
  db <- mk_locked_db()
  holder <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(holder), add = TRUE)
  DBI::dbExecute(holder, "BEGIN IMMEDIATE")
  DBI::dbExecute(holder, "CREATE TABLE lock_probe(x)")

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA busy_timeout = 50;")
  expect_false(RiboPilot:::db_in_transaction(con))   # busy, not transacting

  DBI::dbExecute(con, "BEGIN DEFERRED")
  expect_true(RiboPilot:::db_in_transaction(con))
  DBI::dbExecute(con, "ROLLBACK")
  DBI::dbExecute(holder, "ROLLBACK")
  expect_false(RiboPilot:::db_in_transaction(con))
})

test_that("a lock that will not clear fails loudly instead of migrating unlocked", {
  db <- mk_locked_db()
  holder <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(holder), add = TRUE)
  DBI::dbExecute(holder, "BEGIN IMMEDIATE")
  DBI::dbExecute(holder, "CREATE TABLE lock_probe(x)")

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbExecute(con, "PRAGMA busy_timeout = 50;")
  expect_error(RiboPilot:::begin_migration_txn(con, tries = 2L), "Could not lock")
  # and it leaves no half-open transaction behind
  expect_false(RiboPilot:::db_in_transaction(con))
  DBI::dbExecute(holder, "ROLLBACK")
})

test_that("migrate_db waits for a concurrent write lock rather than working without it", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if_not_installed("parallel")
  db <- mk_locked_db()
  # A sentinel written just before the lock is released, rather than an elapsed
  # time: the holder's head start is not synchronised, so a wall-clock margin
  # fails under load on a busy machine while proving nothing extra.
  released <- file.path(dirname(db), "released")
  holder <- parallel::mcparallel({
    h <- DBI::dbConnect(RSQLite::SQLite(), db)
    DBI::dbExecute(h, "BEGIN IMMEDIATE")
    DBI::dbExecute(h, "CREATE TABLE lock_probe(x)")
    Sys.sleep(2)
    file.create(released)
    DBI::dbExecute(h, "ROLLBACK")
    DBI::dbDisconnect(h)
    TRUE
  })
  Sys.sleep(0.5)   # let the lock be taken before we try

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  migrate_db(con)
  # the lock cannot have been taken before the holder let go, so a migration
  # that returned without seeing the sentinel ran unserialized
  expect_true(file.exists(released))
  parallel::mccollect(holder)
  expect_true("triage_opts" %in% DBI::dbListTables(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
})

test_that("migrate_db sets its own busy timeout and never lowers the caller's", {
  # correctness must not depend on the caller having used open_project_db()
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA busy_timeout")[[1]], 0)
  migrate_db(con)
  expect_gte(DBI::dbGetQuery(con, "PRAGMA busy_timeout")[[1]],
    RiboPilot:::MIGRATE_BUSY_TIMEOUT_MS)

  other <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(other), add = TRUE)
  DBI::dbExecute(other, "PRAGMA busy_timeout = 60000;")
  migrate_db(other)
  expect_equal(DBI::dbGetQuery(other, "PRAGMA busy_timeout")[[1]], 60000)
})

test_that("migrate_db keeps a caller's own transaction and does not commit it", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "BEGIN IMMEDIATE")
  migrate_db(con)
  expect_true(RiboPilot:::db_in_transaction(con))   # still the caller's to end
  DBI::dbExecute(con, "ROLLBACK")
  # rolled back with the caller's transaction, schema stamp included
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 8)
  expect_false("blast_hits" %in% DBI::dbListTables(con))
})

test_that("the v13 backfill runs inside a caller's transaction without nesting", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE blast_hits(ID TEXT, path INT, scaffold INT,
    hit_rank INT, title TEXT, bitscore REAL, hit_class TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT,
    is_nuclear_rrna INT)")
  DBI::dbExecute(con, "INSERT INTO blast_hits VALUES('S1',1,1,1,
    'Danio rerio 18S ribosomal RNA gene, complete sequence', 900, NULL)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,0)")
  DBI::dbExecute(con, "PRAGMA user_version = 12")

  DBI::dbExecute(con, "BEGIN IMMEDIATE")
  expect_silent(migrate_db(con))
  DBI::dbExecute(con, "COMMIT")

  expect_equal(DBI::dbGetQuery(con, "SELECT hit_class FROM blast_hits")$hit_class, "rdna")
})

# --- value rewrites are gated on the version, like every other one ---

test_that("the assemble_switch = 5 rewrite does not run on a current database", {
  # It is a one-off v16 repair, but it ran on every call, i.e. once per unit of
  # a fanned-out run. A state 5 written after v16 means something else entirely.
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S5',5,0,'default','default')")
  DBI::dbExecute(con, sprintf("PRAGMA user_version = %d",
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con,
    "SELECT assemble_switch s FROM assemble WHERE ID='S5'")$s, 5L)
})

# --- the catalogue snapshot ---

test_that("the catalogue snapshot sees tables and columns created in the same call", {
  # migrate_db() answers ~40 guards from one snapshot of the catalogue. A stale
  # snapshot is silent: a table created at v12 and guarded at v14 would look
  # absent, and the v14 column would never be added.
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("triage_opts" %in% DBI::dbListTables(con))
  expect_false("reference" %in% DBI::dbListTables(con))

  migrate_db(con)

  # created at v12, then guarded at v14
  expect_true("rule_screen_off_target" %in% DBI::dbListFields(con, "triage_opts"))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT rule_screen_off_target v FROM triage_opts")$v, "flag")
  # created at v19, then guarded two lines later
  expect_true("preset_name" %in% DBI::dbListFields(con, "reference"))
  # created at v9, then read by the v9 cleanup
  expect_true("blast_hits" %in% DBI::dbListTables(con))
  # a column added at v11 and read by the v11 rewrite in the same call
  expect_true(all(c("triage", "triage_reason", "triage_source") %in%
    DBI::dbListFields(con, "assemblies")))
})

test_that("the snapshot tracks a table that is dropped and rebuilt mid-migration", {
  # v7 rebuilds annotate/annotations/coverage under the per-unit key by
  # DROP + CREATE (annotate via a rename), so the snapshot has to follow both.
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT,
    ignore INT, PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT PRIMARY KEY, annotate_switch INT,
    annotate_lock INT, reviewed TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, element TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, mean_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE stitched(ID TEXT PRIMARY KEY, sequence TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'ACGT',0)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',2,0,'yes')")
  DBI::dbExecute(con, "PRAGMA user_version = 6")

  migrate_db(con)

  # annotate came back under its new name with the v1/v2/v12 columns intact
  expect_true(all(c("path", "scaffold", "annotate_opts", "curate_opts", "n_flags") %in%
    DBI::dbListFields(con, "annotate")))
  expect_false("annotate_new" %in% DBI::dbListTables(con))
  expect_true("path" %in% DBI::dbListFields(con, "annotations"))
  # coverage is rebuilt at v7 and then gains median_depth at v10 and depths at v12
  expect_true(all(c("path", "median_depth", "depths") %in% DBI::dbListFields(con, "coverage")))
  # dropped tables leave the snapshot too
  expect_false("stitched" %in% DBI::dbListTables(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)

  expect_silent(migrate_db(con))
})

test_that("a second migrate_db changes nothing at all", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  mk_v15_assemble_opts(con)
  migrate_db(con)
  schema_of <- function(con) {
    DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY type, name")
  }
  before <- schema_of(con)
  rows <- lapply(c("samples", "assemblies", "assemble", "triage_opts", "curate_opts",
                   "annotate_opts", "assemble_opts", "reference_config"),
    function(t) DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", t)))

  expect_silent(migrate_db(con))

  expect_equal(schema_of(con), before)
  for (i in seq_along(rows)) {
    expect_equal(DBI::dbGetQuery(con, sprintf("SELECT * FROM %s",
      c("samples", "assemblies", "assemble", "triage_opts", "curate_opts",
        "annotate_opts", "assemble_opts", "reference_config")[i])), rows[[i]])
  }
})

# --- the pre-migration backup ---

test_that("the pre-migration backup is a real copy, WAL content included", {
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "PRAGMA journal_mode = wal;")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1')")
  DBI::dbExecute(con, "PRAGMA user_version = 15")
  DBI::dbDisconnect(con)

  b1 <- RiboPilot:::backup_project_db(dir)
  expect_true(file.exists(b1))
  # file.copy() of the .sqlite alone would miss everything still in the -wal
  bcon <- DBI::dbConnect(RSQLite::SQLite(), b1); on.exit(DBI::dbDisconnect(bcon))
  expect_equal(DBI::dbGetQuery(bcon, "SELECT ID FROM samples")$ID, "S1")
  expect_equal(DBI::dbGetQuery(bcon, "PRAGMA user_version")[[1]], 15)

  # two backups in the same second must both exist
  b2 <- RiboPilot:::backup_project_db(dir)
  expect_true(file.exists(b2))
  expect_false(identical(b1, b2))
})

test_that("backup_project_db says so when it cannot write a backup", {
  skip_on_os("windows")
  skip_if(as.integer(Sys.getenv("EUID", "1")) == 0 || unname(Sys.info()["user"]) == "root",
    "root ignores directory permissions")
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY)")
  DBI::dbDisconnect(con)
  bdir <- file.path(dir, ".old_sqlite_dbs")
  dir.create(bdir)
  Sys.chmod(bdir, "500")
  on.exit(Sys.chmod(bdir, "700"), add = TRUE)

  expect_warning(res <- RiboPilot:::backup_project_db(dir), "Could not back up")
  expect_null(res)
})
