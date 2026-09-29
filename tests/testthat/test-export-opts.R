# Saved FASTA-header configurations, the token data templates resolve against,
# and the checks that stand between a bad template and a rejected submission.

mk_opts_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  # samples carries whatever columns the user's mapping file had: voucher_id is
  # one of those, and must be usable as a token with no extra plumbing
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT,
    specimen_voucher TEXT, voucher_id TEXT, R1 TEXT, R2 TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT,
    topology TEXT, sequence TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, n_features INT,
    missing TEXT, partial TEXT, reviewed TEXT, export_group TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT,
    start INT, stop INT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio',7955,
    'MCZ 12345','ABC-1','r1.fq','r2.fq')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,5000,'linear',
    ?, 0)", params = list(strrep("ACGT", 1250)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,5,'','no','yes','batch1')")
  for (f in list(c("18S", 1, 1800), c("ITS1", 1801, 2100), c("5.8S", 2101, 2260),
                 c("ITS2", 2261, 2600), c("28S", 2601, 5000))) {
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,?,?,?)",
      params = list(f[1], as.integer(f[2]), as.integer(f[3])))
  }
  con
}

one_unit <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)

test_that("the definition phrase groups runs the way GenBank writes them", {
  # ends run into the scaffold ends, middle is whole: three runs, last one "and"
  txt <- RiboPilot:::rrna_definition_text(
    c("18S", "ITS1", "5.8S", "ITS2", "28S"), c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(txt, paste0(
    "18S ribosomal RNA gene, partial sequence; ",
    "internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, ",
    "complete sequence; and 28S ribosomal RNA gene, partial sequence."))

  # a single complete feature is a one-clause sentence with no "and"
  expect_equal(RiboPilot:::rrna_definition_text("18S", FALSE),
    "18S ribosomal RNA gene, complete sequence.")
  # two of a kind read as "A and B"
  expect_equal(RiboPilot:::rrna_definition_text(c("ITS1", "ITS2"), c(FALSE, FALSE)),
    "internal transcribed spacer 1 and internal transcribed spacer 2, complete sequence.")
  # nothing annotated still yields usable prose
  expect_match(RiboPilot:::rrna_definition_text(character(0), logical(0)), "ribosomal")
})

test_that("token data carries the sample's own metadata columns", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::export_header_data(con, one_unit)

  expect_equal(nrow(d), 1L)
  expect_equal(d$voucher_id, "ABC-1")          # a mapping-file column, not a built-in
  expect_equal(d$seqid, "S1")
  expect_equal(d$elements, "18S-ITS1-5.8S-ITS2-28S")
  expect_equal(d$length, 5000L)
  expect_equal(d$export_group, "batch1")
  expect_match(d$rrna_definition, "^18S ribosomal RNA gene, partial sequence;")
  expect_equal(d$completeness, "partial sequence")   # both ends touch the scaffold end
  # read files are plumbing, never header material
  expect_false(any(c("R1", "R2") %in% names(d)))
})

test_that("the default templates resolve against real data", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::export_header_data(con, one_unit)
  h <- RiboPilot:::render_fasta_headers(RiboPilot:::DEFAULT_FASTA_HEADER, d)
  expect_equal(h, paste("S1 [organism=Danio rerio] [moltype=genomic DNA] Danio rerio",
    "18S ribosomal RNA gene, partial sequence;",
    "internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,",
    "complete sequence; and 28S ribosomal RNA gene, partial sequence."))

  d$element <- "ITS1"; d$product <- "internal transcribed spacer 1"
  expect_equal(RiboPilot:::render_fasta_headers(RiboPilot:::DEFAULT_FASTA_HEADER_MARKER, d),
    "S1_ITS1 [organism=Danio rerio] Danio rerio internal transcribed spacer 1")
})

test_that("a custom column from the mapping file is a usable token", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::export_header_data(con, one_unit)
  expect_equal(RiboPilot:::render_fasta_headers("{seqid} voucher {voucher_id}", d),
    "S1 voucher ABC-1")
})

test_that("validation catches the mistakes that get a submission rejected", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::export_header_data(con, one_unit)
  v <- function(...) RiboPilot:::validate_fasta_header(...)

  expect_equal(v(RiboPilot:::DEFAULT_FASTA_HEADER, d)$level, "ok")
  expect_equal(v("", d)$level, "error")
  expect_equal(v("   ", d)$level, "error")

  br <- v("{seqid} [organism={Taxon]", d)
  expect_equal(br$level, "error")
  expect_match(br$message, "Unmatched")

  nope <- v("{seqid} {not_a_column}", d)
  expect_equal(nope$level, "error")
  expect_match(nope$message, "not_a_column")

  # one unit per sample: a missing {seqid} is survivable, so warn
  expect_equal(v("[organism={Taxon}]", d)$level, "warn")
  # two units on one sample: the record and its feature table would disagree
  d2 <- rbind(d, d)
  expect_equal(v("[organism={Taxon}]", d2)$level, "error")
  # ...and the check can be waived where it does not apply
  expect_equal(v("[organism={Taxon}]", d2, require_seqid = FALSE)$level, "ok")
})

test_that("configurations are saved, listed and read back by name", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::ensure_export_opts(con)

  expect_equal(RiboPilot:::list_export_templates(con), "default")
  expect_equal(RiboPilot:::get_export_opts(con)$fasta_header, RiboPilot:::DEFAULT_FASTA_HEADER)

  RiboPilot:::set_export_opts(con, "{seqid} mine", "{seqid}_{element} mine", name = "voucher_first")
  expect_equal(RiboPilot:::list_export_templates(con), c("default", "voucher_first"))
  expect_equal(RiboPilot:::get_export_opts(con, "voucher_first")$fasta_header, "{seqid} mine")

  # saving the same name again updates rather than duplicating
  RiboPilot:::set_export_opts(con, "{seqid} mine v2", "{seqid}_{element}", name = "voucher_first")
  expect_equal(RiboPilot:::get_export_opts(con, "voucher_first")$fasta_header, "{seqid} mine v2")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM export_opts")$n, 2L)

  # a blank field falls back to the package default rather than exporting nothing
  RiboPilot:::set_export_opts(con, "", "", name = "empty")
  expect_equal(RiboPilot:::get_export_opts(con, "empty")$fasta_header, RiboPilot:::DEFAULT_FASTA_HEADER)
  # an unknown name does the same
  expect_equal(RiboPilot:::get_export_opts(con, "no_such_config")$fasta_header,
    RiboPilot:::DEFAULT_FASTA_HEADER)
})

test_that("a project made before this table existed still exports", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::list_export_templates(con), "default")
  expect_equal(RiboPilot:::get_export_opts(con)$fasta_header_marker,
    RiboPilot:::DEFAULT_FASTA_HEADER_MARKER)
})

test_that("[organism=] is the Taxonomy name only when it matched the whole Taxon", {
  eo <- RiboPilot:::export_organism
  expect_equal(eo("danio_rerio", "Danio rerio", "exact", "species"), "Danio rerio")
  expect_equal(eo("Upeneus sp. A", "Upeneus", "exact", "genus"), "Upeneus sp. A")
  expect_equal(eo("Upeneous parvus", "Upeneus parvus", "manual", "species"), "Upeneus parvus")
  expect_equal(eo("Gen. nov. 1", "Mullidae", "manual", "family"), "Gen. nov. 1")
  expect_equal(eo("Upeneus parvis", "Upeneus parvus", "fuzzy_species", "species"), "Upeneus parvis")
  expect_equal(eo("Danio rerio"), "Danio rerio")

  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN match_type TEXT")
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'danio rerio', match_type = 'exact'")
  d <- RiboPilot:::export_header_data(con, one_unit)
  expect_match(RiboPilot:::render_fasta_headers(RiboPilot:::DEFAULT_FASTA_HEADER, d),
    "^S1 \\[organism=Danio rerio\\] \\[moltype=genomic DNA\\] Danio rerio 18S")
})
