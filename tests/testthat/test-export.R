# GenBank export: N-run gaps, feature table.

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
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  seq <- paste0(strrep("A", 100), strrep("N", 20), strrep("C", 140))  # gap 101-120
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, ?)", params = list(seq, nchar(seq)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'ITS1','misc_RNA',121,200)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',201,260)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "<1\t100\trRNA")            # 5' partial
  expect_match(tbl, "201\t>260\trRNA")          # 3' partial
  expect_match(tbl, "101\t120\tassembly_gap")   # gap feature
  expect_match(tbl, "estimated_length\t20")
  expect_match(tbl, "internal transcribed spacer 1")
})

test_that("force_partial marks a curator-flagged partial record at both ends", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  seq <- strrep("A", 500)
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 500)", params = list(seq))
  # 18S well inside the sequence (100-260), so no positional partial
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',100,260)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',300,460)")
  complete <- RiboPilot:::feature_table_text_con(con, "S1", force_partial = FALSE)
  expect_match(complete, "100\t260\trRNA")     # not marked partial
  partial <- RiboPilot:::feature_table_text_con(con, "S1", force_partial = TRUE)
  expect_match(partial, "<100\t260\trRNA")     # 5' partial on first feature
  expect_match(partial, "300\t>460\trRNA")     # 3' partial on last feature
})

test_that("submission_risks flags chimeric, in-span gaps, and misassembly", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, chimeric INT)")

  # clean, contiguous, no gaps
  clean <- paste0(strrep("A", 100), strrep("C", 100), strrep("G", 100))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('OK',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('OK',1,1, 0)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'5.8S','rRNA',101,200)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'28S','rRNA',201,300)")
  expect_equal(RiboPilot:::submission_risks(con, "OK")$note, "")

  # N-run inside the 18S..28S span -> gap_in_rDNA
  gapped <- paste0(strrep("A", 100), strrep("N", 20), strrep("C", 180))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('GAP',1,1, ?, ?)", params = list(gapped, nchar(gapped)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('GAP',1,1, 0)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('GAP',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('GAP',1,1,'28S','rRNA',121,300)")
  r <- RiboPilot:::submission_risks(con, "GAP")
  expect_true(r$gap_in_rdna)
  expect_match(r$note, "gap_in_rDNA")

  # chimeric flag from the annotate table
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('CHIM',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('CHIM',1,1, 1)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('CHIM',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('CHIM',1,1,'28S','rRNA',201,300)")
  expect_true(RiboPilot:::submission_risks(con, "CHIM")$chimeric)

  # overlapping / out-of-order features -> misassembly
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('BAD',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('BAD',1,1, 0)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('BAD',1,1,'18S','rRNA',1,150)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('BAD',1,1,'5.8S','rRNA',120,200)")
  expect_true(RiboPilot:::submission_risks(con, "BAD")$misassembly)
})

test_that("risk_flags reports suspect and likely chimera labels", {
  r <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "suspect", chimera_reason = "why")
  expect_true(r$chimeric)
  expect_match(r$note, "chimeric")
  expect_false(RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "clean")$chimeric)
})

test_that("risk_flags strips ';' out of the embedded chimera reason", {
  r <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "likely",
    chimera_reason = "18S resolves to Fish; 28S to Fungi")
  paren <- sub("^chimeric\\((.*)\\)$", "\\1", r$note)
  expect_false(grepl(";", paren, fixed = TRUE))
  expect_match(paren, "18S resolves to Fish, 28S to Fungi", fixed = TRUE)
})

test_that("export_summary_df reads the chimera table into submission_risk", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1',1,1)")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, reason) VALUES ('S1',1,1,'likely','18S vs 28S mismatch')")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)
  smry <- RiboPilot:::export_summary_df(con, units)
  expect_match(smry$submission_risk, "chimeric")
})

test_that("feature_map_svg_text draws one arrow per element, in order, in palette colours", {
  d <- data.frame(
    element = c("18S", "ITS1", "5.8S", "ITS2", "28S"),
    start = c(1L, 201L, 241L, 261L, 301L),
    stop = c(200L, 240L, 260L, 300L, 400L), stringsAsFactors = FALSE)
  svg <- RiboPilot:::feature_map_svg_text(d)
  expect_match(svg, "^<svg xmlns='http://www.w3.org/2000/svg'")
  expect_match(svg, "</svg>$")
  expect_equal(lengths(regmatches(svg, gregexpr("<polygon", svg, fixed = TRUE))), 5L)
  # labels in 5'->3' order
  labs <- regmatches(svg, gregexpr("(?<=font-size='11' fill='#000000'>)[^<]+", svg, perl = TRUE))[[1]]
  expect_equal(labs, d$element)
  # element fills come from the shared palette
  fills <- regmatches(svg, gregexpr("(?<=fill=')#[0-9A-F]{6}(?=' stroke='#333333')", svg, perl = TRUE))[[1]]
  expect_equal(fills, unname(RiboPilot:::RRNA_ELEMENT_COLORS[d$element]))
  # left-to-right layout matches the coordinate order
  x <- as.numeric(regmatches(svg, gregexpr("(?<=<polygon points=')[0-9.]+", svg, perl = TRUE))[[1]])
  expect_false(is.unsorted(x))
  # bp axis labels
  expect_match(svg, ">400</text>")
})

test_that("feature_map_svg falls back to a note when a unit has no annotations", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT, stop INT)")
  expect_equal(RiboPilot:::feature_map_svg(con, "S1"), "<em>no annotations</em>")
})

test_that("is_nuclear_rrna_title splits nuclear rDNA from organellar/other", {
  f <- RiboPilot:::is_nuclear_rrna_title
  expect_true(f("Symbiodinium sp. 18S ribosomal RNA gene, partial"))
  expect_true(f("Xenia sp. internal transcribed spacer 1, 5.8S ribosomal RNA"))
  expect_false(f("Homo sapiens mitochondrion, complete genome"))
  expect_false(f("Escherichia coli 16S ribosomal RNA gene"))
  expect_false(f("Some fish mitochondrial 12S ribosomal RNA"))
})
