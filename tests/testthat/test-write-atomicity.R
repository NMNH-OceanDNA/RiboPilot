# Every replace in the pipeline is a DELETE followed by an INSERT. Without a
# transaction the delete commits on its own, so an insert that fails takes the
# rows with it. Each test here blocks the insert with a trigger and asserts what
# was already stored is still there.

FISH_TAX <- "33208;7711;186623;7952;2743709;7954;7955"

block_inserts <- function(con, tbl) {
  DBI::dbExecute(con, sprintf(
    "CREATE TRIGGER block_%s BEFORE INSERT ON %s BEGIN SELECT RAISE(ABORT, 'blocked'); END",
    tbl, tbl))
}

test_that("a failed screen insert leaves the stored verdict in place", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage_taxids TEXT, reference_support INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT,
    sequence TEXT, length INT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',?,NULL)", params = list(FISH_TAX))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL)")
  for (i in 1:3) {
    DBI::dbExecute(con,
      "INSERT INTO blast_hits (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class)
       VALUES ('S1',1,1,?,1,?,?,'Eukaryota;Chordata;Actinopteri;o;f;g;s','rdna')",
      params = list(i, 900 - i, FISH_TAX))
  }

  screen_one(con, "S1", 1, 1)
  before <- DBI::dbGetQuery(con, "SELECT * FROM screen")
  expect_equal(nrow(before), 1L)

  block_inserts(con, "screen")
  expect_error(screen_one(con, "S1", 1, 1, force = TRUE), "blocked")

  after <- DBI::dbGetQuery(con, "SELECT * FROM screen")
  expect_equal(nrow(after), 1L)
  expect_equal(after$screen_verdict, before$screen_verdict)
})

test_that("a failed assemblies insert leaves the sample's scaffolds in place", {
  dir <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT,
    length INT, length_raw INT, sequence TEXT, source_contig TEXT, ignore INT, edited INT, time_stamp INT,
    triage TEXT, triage_reason TEXT, triage_source TEXT, triage_rule TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, paths INT, scaffolds INT, length TEXT,
    topology TEXT, assemble_switch INT, assemble_notes TEXT, assemble_opts TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, max_paths INT,
    max_scaffolds INT, min_assembly_length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT,
    PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default', 10, 10, NULL)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',NULL,NULL,NULL,NULL,1,NULL,'default',NULL)")

  ad <- file.path(dir, "S1", "assemble")
  dir.create(ad, recursive = TRUE)
  writeLines(c(
    ">S1.1.1 circular", strrep("ACGT", 200),
    ">S1.1.2 linear", strrep("TGCA", 150)),
    file.path(ad, "S1_assembly_1.fasta"))

  RiboPilot:::load_assemblies(con, dir, "S1")
  before <- DBI::dbGetQuery(con, "SELECT path, scaffold, length FROM assemblies ORDER BY path, scaffold")
  expect_equal(nrow(before), 2L)

  block_inserts(con, "assemblies")
  expect_error(RiboPilot:::load_assemblies(con, dir, "S1"), "blocked")

  after <- DBI::dbGetQuery(con, "SELECT path, scaffold, length FROM assemblies ORDER BY path, scaffold")
  expect_equal(after, before)
})

test_that("a failed reference insert leaves the groups and their membership in place", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT, ref_id TEXT, user_group TEXT, ref_source TEXT)")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_DDL)
  for (id in c("S1", "S2")) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,NULL,NULL,NULL,1,2,NULL,NULL,NULL)", params = list(id))
  }
  grp <- data.frame(ref_id = "asteroidea", clade_taxid = 7588L, clade_name = "Asteroidea",
    clade_rank = "class", n_records = 4102L, n_samples = 2L, stringsAsFactors = FALSE)
  mem <- data.frame(ID = c("S1", "S2"), ref_id = "asteroidea", stringsAsFactors = FALSE)

  RiboPilot:::write_reference_groups(con, grp, mem)
  # a built database somebody waited ten minutes for
  DBI::dbExecute(con, "UPDATE reference SET reference_switch = 2, source = 'built',
    seeds_db = '/refs/a_seeds.fa', labels_db = '/refs/a_labels.fa'")
  before <- DBI::dbGetQuery(con, "SELECT * FROM reference")
  members_before <- DBI::dbGetQuery(con, "SELECT ID, ref_id FROM samples ORDER BY ID")

  block_inserts(con, "reference")
  moved <- grp
  moved$clade_taxid <- 7586L
  expect_error(RiboPilot:::write_reference_groups(con, moved, mem), "blocked")

  expect_equal(DBI::dbGetQuery(con, "SELECT * FROM reference"), before)
  expect_equal(DBI::dbGetQuery(con, "SELECT ID, ref_id FROM samples ORDER BY ID"), members_before)
})

test_that("a failed path_scores insert leaves the stored scores in place", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE path_scores (ID TEXT NOT NULL, path INTEGER NOT NULL,
    length INTEGER, depth REAL, anchors TEXT, blast_pident REAL, hit_id TEXT, hit_taxid TEXT,
    hit_species TEXT, blast_lineage TEXT, s_tax REAL, s_identity REAL, s_anchors REAL,
    s_depth REAL, s_length REAL, score REAL, rank INTEGER, auto_select INTEGER,
    time_stamp INTEGER, PRIMARY KEY (ID, path))")
  scored <- data.frame(path = 1:2, length = c(5000L, 4000L), score = c(0.9, 0.8),
    rank = 1:2, auto_select = c(TRUE, FALSE), stringsAsFactors = FALSE)

  RiboPilot:::write_path_scores(con, "S1", scored)
  before <- DBI::dbGetQuery(con, "SELECT ID, path, score FROM path_scores ORDER BY path")
  expect_equal(nrow(before), 2L)

  block_inserts(con, "path_scores")
  expect_error(RiboPilot:::write_path_scores(con, "S1", scored), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT ID, path, score FROM path_scores ORDER BY path"), before)
})

test_that("a failed reference-record insert leaves the fetched sequence in place", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences (accession TEXT PRIMARY KEY,
    species TEXT, title TEXT, taxid TEXT, lineage TEXT, sequence TEXT, length INTEGER,
    time_stamp INTEGER)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations (accession TEXT, element TEXT,
    type TEXT, start INTEGER, stop INTEGER, strand TEXT, time_stamp INTEGER,
    PRIMARY KEY (accession, element))")

  d <- withr::local_tempdir()
  seq <- paste0(strrep("ACGT", 10), strrep("GT", 10))
  fa <- file.path(d, "MZ1.fa")
  writeLines(c(">MZ1.1 Fake species small subunit ribosomal RNA gene", seq), fa)
  ann <- data.frame(seqid = "MZ1.1", element = "18S", type = "rRNA", start = 1L,
    stop = nchar(seq), strand = "+", stringsAsFactors = FALSE)

  RiboPilot:::write_blast_refs(con, fa, ann)
  before <- DBI::dbGetQuery(con, "SELECT accession, sequence FROM blast_ref_sequences")
  expect_equal(nrow(before), 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_annotations")$n, 1L)

  # the features fail, and the record they belong to must not vanish with them
  block_inserts(con, "blast_ref_annotations")
  expect_error(RiboPilot:::write_blast_refs(con, fa, ann), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT accession, sequence FROM blast_ref_sequences"), before)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_annotations")$n, 1L)
})

test_that("a failed blast_hits insert leaves the unit's stored hits in place", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL,
    bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT,
    title TEXT, species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT,
    pident REAL, qcovs REAL, evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, blast_accession TEXT,
    blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT,
    ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,?)", params = list(FISH_TAX))

  f <- withr::local_tempfile()
  writeLines(paste("AA1", "Danio rerio 18S ribosomal RNA gene", 99, 80, "1e-90", 7955,
    1500, 1600, 1, 1600, 1, 1600, sep = "\t"), f)

  parse_blast_candidates(con, "S1", 1, 1, f)
  before <- DBI::dbGetQuery(con, "SELECT ID, accession FROM blast_hits")
  expect_equal(nrow(before), 1L)

  block_inserts(con, "blast_hits")
  expect_error(parse_blast_candidates(con, "S1", 1, 1, f), "blocked")
  expect_equal(DBI::dbGetQuery(con, "SELECT ID, accession FROM blast_hits"), before)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 1L)
})
