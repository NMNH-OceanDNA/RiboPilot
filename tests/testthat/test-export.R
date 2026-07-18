# GenBank export: source-modifier allowlist, N-run gaps, feature table.

test_that("source_modifiers allowlists real INSDC qualifiers only", {
  row <- data.frame(ID = "S1", Taxon = "Foo bar", taxonomy_switch = 1L,
    taxonomy_lock = 0L, rank = "species", specimen_voucher = "USNM:1",
    geo_loc_name = "USA", plate_well = "A1", stringsAsFactors = FALSE)
  mods <- RiboPilot:::source_modifiers(row)
  expect_true(any(grepl("specimen_voucher=USNM:1", mods)))
  expect_true(any(grepl("geo_loc_name=USA", mods)))
  expect_false(any(grepl("taxonomy_switch", mods)))   # internal columns dropped
  expect_false(any(grepl("taxonomy_lock", mods)))
  expect_false(any(grepl("plate_well", mods)))        # unrecognized column dropped
})

test_that("find_n_runs finds runs >= threshold only", {
  s <- paste0(strrep("A", 50), strrep("N", 30), strrep("C", 40), "NNN", strrep("G", 20))
  runs <- RiboPilot:::find_n_runs(s, min_len = 10)
  expect_equal(nrow(runs), 1)
  expect_equal(runs$start, 51)
  expect_equal(runs$stop, 80)
  expect_equal(nrow(RiboPilot:::find_n_runs(NA_character_)), 0)
})

test_that("feature_table_text_con emits partials, gaps, and ITS notes", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE stitched(ID TEXT, sequence TEXT, length INT)")
  seq <- paste0(strrep("A", 100), strrep("N", 20), strrep("C", 140))  # gap 101-120
  DBI::dbExecute(con, "INSERT INTO stitched VALUES('S1', ?, ?)", params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1','18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1','ITS1','misc_RNA',121,200)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1','28S','rRNA',201,260)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "<1\t100\trRNA")            # 5' partial
  expect_match(tbl, "201\t>260\trRNA")          # 3' partial
  expect_match(tbl, "101\t120\tassembly_gap")   # gap feature
  expect_match(tbl, "estimated_length\t20")
  expect_match(tbl, "internal transcribed spacer 1")
})

test_that("is_nuclear_rrna_title splits nuclear rDNA from organellar/other", {
  f <- RiboPilot:::is_nuclear_rrna_title
  expect_true(f("Symbiodinium sp. 18S ribosomal RNA gene, partial"))
  expect_true(f("Xenia sp. internal transcribed spacer 1, 5.8S ribosomal RNA"))
  expect_false(f("Homo sapiens mitochondrion, complete genome"))
  expect_false(f("Escherichia coli 16S ribosomal RNA gene"))
  expect_false(f("Some fish mitochondrial 12S ribosomal RNA"))
})
