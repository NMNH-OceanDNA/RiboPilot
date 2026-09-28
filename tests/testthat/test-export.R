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
  expect_match(tbl, "linkage_evidence\tunspecified")
  expect_match(tbl, "internal transcribed spacer 1")
})

test_that("partial ends are marked per feature, strand-aware, and minus is written stop first", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT,
    start INT, stop INT, strand TEXT, partial_start INT, partial_stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 500)", params = list(strrep("A", 500)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',100,260,'+',1,0)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',300,460,'+',0,1)")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "\n<100\t260\trRNA")
  expect_match(tbl, "\n300\t>460\trRNA")

  # minus strand: 5' end is the genomic right end, written first
  DBI::dbExecute(con, "UPDATE annotations SET strand = '-'")
  tbl <- RiboPilot:::feature_table_text_con(con, "S1")
  expect_match(tbl, "\n<260\t100\trRNA")
  expect_match(tbl, "\n460\t>300\trRNA")

  # a minus feature touching base 1 is 3' partial there, unless set otherwise
  DBI::dbExecute(con, "UPDATE annotations SET start = 1, partial_stop = NULL WHERE element = '18S'")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<260\t>1\trRNA")
  DBI::dbExecute(con, "UPDATE annotations SET partial_stop = 0 WHERE element = '18S'")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<260\t1\trRNA")
})

test_that("a scaffold end marks a gene partial only on an untrimmed unit", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT,
    start INT, stop INT, strand TEXT, partial_start INT, partial_stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, RiboPilot:::ASSEMBLY_BACKUP_DDL)
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 300)", params = list(strrep("A", 300)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',1,300,'+',NULL,NULL)")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t>300\trRNA")
  # trimmed: the snapshot holds a longer sequence, so the ends are cut, not assembly ends
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence) VALUES ('S1',1,1,?)",
    params = list(strrep("A", 400)))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n1\t300\trRNA")
  # restored: snapshot and sequence agree again
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = ?", params = list(strrep("A", 300)))
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n<1\t>300\trRNA")
})

test_that("a unit marked partial changes only the completeness wording", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, n_features INT, missing TEXT,
    partial TEXT, reviewed TEXT, export_group TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT,
    start INT, stop INT, partial_start INT, partial_stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Genus species')")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,1,'','yes','yes',NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 500, 'linear')", params = list(strrep("A", 500)))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S','rRNA',100,260,0,0)")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n100\t260\trRNA")
  h <- RiboPilot:::export_header_data(con, data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1"))
  expect_equal(h$completeness, "partial sequence")
  expect_equal(h$rrna_definition, "18S ribosomal RNA gene, complete sequence.")
})

test_that("a partial end is exported partial, in the feature table and the defline", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, n_features INT, missing TEXT,
    partial TEXT, reviewed TEXT, export_group TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT,
    start INT, stop INT, flags TEXT, partial_start INT, partial_stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT, topology TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Genus species')")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,1,'','no','yes',NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1, ?, 1694, 'linear')", params = list(strrep("A", 1694)))
  # 28S well inside the scaffold, 3' end not recovered
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'28S','rRNA',155,1488,'length,truncated',0,1)")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n155\t>1488\trRNA")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1")
  h <- RiboPilot:::export_header_data(con, units)
  expect_equal(h$completeness, "partial sequence")
  expect_match(h$rrna_definition, "28S ribosomal RNA gene, partial sequence")

  # the truncated flag alone does not mark ends: curation sets partial_start/stop
  DBI::dbExecute(con, "UPDATE annotations SET partial_stop = 0")
  expect_match(RiboPilot:::feature_table_text_con(con, "S1"), "\n155\t1488\trRNA")
  expect_equal(RiboPilot:::export_header_data(con, units)$completeness, "complete sequence")
})

test_that("export summary flags chimeric, in-span gaps, and misassembly", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, match_type TEXT, rank TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, length INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, chimeric INT, n_features INT, missing TEXT, reviewed TEXT, partial TEXT, problematic TEXT, curate_flags TEXT)")
  risk <- function(id) RiboPilot:::export_summary_df(con,
    data.frame(ID = id, path = 1L, scaffold = 1L, seqid = id))$submission_risk

  # clean, contiguous, no gaps
  clean <- paste0(strrep("A", 100), strrep("C", 100), strrep("G", 100))
  DBI::dbExecute(con, "INSERT INTO samples VALUES('OK','Danio rerio',NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('OK',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('OK',1,1, 0, 3, '', 'no', 'no', 'no', NULL)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'5.8S','rRNA',101,200)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('OK',1,1,'28S','rRNA',201,300)")
  expect_equal(risk("OK"), "")

  # N-run inside the 18S..28S span -> gap_in_rDNA
  gapped <- paste0(strrep("A", 100), strrep("N", 20), strrep("C", 180))
  DBI::dbExecute(con, "INSERT INTO samples VALUES('GAP','Danio rerio',NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('GAP',1,1, ?, ?)", params = list(gapped, nchar(gapped)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('GAP',1,1, 0, 2, '', 'no', 'no', 'no', NULL)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('GAP',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('GAP',1,1,'28S','rRNA',121,300)")
  expect_match(risk("GAP"), "gap_in_rDNA")

  # chimeric flag from the annotate table
  DBI::dbExecute(con, "INSERT INTO samples VALUES('CHIM','Danio rerio',NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('CHIM',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('CHIM',1,1, 1, 2, '', 'no', 'no', 'no', NULL)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('CHIM',1,1,'18S','rRNA',1,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('CHIM',1,1,'28S','rRNA',201,300)")
  expect_match(risk("CHIM"), "chimeric")

  # overlapping / out-of-order features -> misassembly
  DBI::dbExecute(con, "INSERT INTO samples VALUES('BAD','Danio rerio',NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('BAD',1,1, ?, ?)", params = list(clean, nchar(clean)))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('BAD',1,1, 0, 2, '', 'no', 'no', 'no', NULL)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('BAD',1,1,'18S','rRNA',1,150)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('BAD',1,1,'5.8S','rRNA',120,200)")
  expect_match(risk("BAD"), "misassembly")
})

test_that("risk_flags flags a severe chimera or ITSx's call, never a minor one", {
  r <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "severe", chimera_reason = "why")
  expect_true(r$chimeric)
  expect_equal(r$note, "chimeric(why)")
  m <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "minor", chimera_reason = "why")
  expect_false(m$chimeric)
  expect_equal(m$note, "")
  expect_false(RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "clean")$chimeric)
  expect_equal(RiboPilot:::risk_flags(1L, NULL, NA_character_, chimera_label = "minor", chimera_reason = "why")$note,
    "chimeric(ITSx)")
})

test_that("risk_flags strips ';' out of the embedded chimera reason", {
  r <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "severe",
    chimera_reason = "18S resolves to Fish; 28S to Fungi")
  paren <- sub("^chimeric\\((.*)\\)$", "\\1", r$note)
  expect_false(grepl(";", paren, fixed = TRUE))
  expect_match(paren, "18S resolves to Fish, 28S to Fungi", fixed = TRUE)
})

test_that("export_summary_df puts a severe chimera in submission_risk and a minor one in chimera_warning", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1',1,1)")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, reason) VALUES ('S1',1,1,'severe','18S vs 28S mismatch')")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)
  smry <- RiboPilot:::export_summary_df(con, units)
  expect_match(smry$submission_risk, "chimeric")
  expect_equal(smry$chimera_warning, "")
  DBI::dbExecute(con, "UPDATE chimera SET label = 'minor', reason = 'depth step at 3000'")
  smry <- RiboPilot:::export_summary_df(con, units)
  expect_false(grepl("chimeric", smry$submission_risk))
  expect_equal(smry$chimera_warning, "depth step at 3000")
})

test_that("export_summary_df puts off_target and contaminant into submission_risk", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1',1,1)")
  DBI::dbExecute(con, "UPDATE annotate SET curate_flags = 'missing_5.8S; contaminant' WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict, lca_taxon) VALUES ('S1',1,1,'off_target','Actiniaria')")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)
  risk <- RiboPilot:::export_summary_df(con, units)$submission_risk
  expect_match(risk, "off_target(Actiniaria)", fixed = TRUE)
  expect_match(risk, "(^|;)contaminant(;|$)")
  expect_true(RiboPilot:::has_screen_risk(risk))

  # an unchecked screen is not a contaminant call
  DBI::dbExecute(con, "UPDATE annotate SET curate_flags = 'contaminant_unchecked'")
  DBI::dbExecute(con, "UPDATE screen SET screen_verdict = 'on_target'")
  risk <- RiboPilot:::export_summary_df(con, units)$submission_risk
  expect_false(RiboPilot:::has_screen_risk(risk))
  expect_equal(RiboPilot:::has_screen_risk(c(NA, "", "chimeric(ITSx)", "gap_in_rDNA", "missing_organism")),
    c(FALSE, FALSE, TRUE, FALSE, FALSE))
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

test_that("a nuclear rDNA title classifies as rdna, organellar and bacterial ones do not", {
  f <- function(x) RiboPilot:::classify_hit_title(x) == "rdna"
  expect_true(f("Symbiodinium sp. 18S ribosomal RNA gene, partial"))
  expect_true(f("Xenia sp. internal transcribed spacer 1, 5.8S ribosomal RNA"))
  expect_false(f("Homo sapiens mitochondrion, complete genome"))
  expect_false(f("Escherichia coli 16S ribosomal RNA gene"))
  expect_false(f("Some fish mitochondrial 12S ribosomal RNA"))
})

test_that("export_summary_df flags a Taxon that differs from the Taxonomy name", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO annotate (ID, path, scaffold) VALUES ('S1',1,1)")
  units <- data.frame(ID = "S1", path = 1L, scaffold = 1L, seqid = "S1", stringsAsFactors = FALSE)
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'Upeneus sp. A', sci_name = 'Upeneus',
    match_type = 'exact', rank = 'genus'")
  s <- RiboPilot:::export_summary_df(con, units)
  expect_equal(s$organism, "Upeneus sp. A")
  expect_match(s$submission_risk, "taxon_differs(Upeneus)", fixed = TRUE)
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'Upeneus parvus', sci_name = 'Upeneus parvus', rank = 'species'")
  s <- RiboPilot:::export_summary_df(con, units)
  expect_false(grepl("taxon_differs", s$submission_risk))
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'upeneus_parvus '")
  s <- RiboPilot:::export_summary_df(con, units)
  expect_false(grepl("taxon_differs", s$submission_risk))
})

test_that("unit_checks names every finding, colours the worst, and agrees with the export prompt", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  expect_equal(nrow(RiboPilot:::unit_checks(con, "sample")), 0L)
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length, ignore) VALUES ('S1',1,?,'ACGT',4,0)",
      params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1',1,?)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1',1,?,'on_target')",
      params = list(s))
  }
  # s1 clean; s2 minor chimera + a truncated feature; s3 off target and contaminant; s4 dropped
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1',1,2,'minor')")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, status, flags) VALUES ('S1',1,2,'28S','flag','truncated')")
  DBI::dbExecute(con, "UPDATE screen SET screen_verdict = 'off_target', lca_taxon = 'Actiniaria' WHERE scaffold = 3")
  DBI::dbExecute(con, "UPDATE annotate SET curate_flags = 'contaminant' WHERE scaffold = 3")
  DBI::dbExecute(con, "UPDATE screen SET screen_verdict = 'off_target' WHERE scaffold = 4")
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE scaffold = 4")

  u <- RiboPilot:::unit_checks(con, "unit")
  expect_equal(u$scaffold, 1:3)
  expect_equal(u$checks, c("ok", "minor chimera, feature flags", "off target, contaminant"))
  expect_equal(u$checks_tone, c("success", "warning", "danger"))
  # red is exactly what the Export and Lock prompts ask about
  smry <- RiboPilot:::export_summary_df(con, cbind(u[1:3], seqid = "x"))
  expect_equal(u$checks_tone == "danger", RiboPilot:::has_screen_risk(smry$submission_risk))

  # one row per sample, over its kept units only, worst tone
  s <- RiboPilot:::unit_checks(con, "sample")
  expect_equal(s$ID, "S1")
  expect_equal(s$checks, "off target, contaminant, minor chimera, feature flags")
  expect_equal(s$checks_tone, "danger")

  # nothing flagged but no usable screen is grey, not green
  DBI::dbExecute(con, "DELETE FROM screen WHERE scaffold = 1")
  expect_equal(RiboPilot:::unit_checks(con, "unit")[1, c("checks", "checks_tone")],
    data.frame(checks = "no reference", checks_tone = "neutral"))
})

test_that("the export done dialog names flagged records and is green only when nothing is flagged", {
  smry <- data.frame(seqid = c("A", "B_p1_s1"), itsx_check = c("", ""),
    submission_risk = c("", "off_target(Actiniaria);contaminant"))
  html <- as.character(RiboPilot:::export_check_summary_ui(smry))
  expect_match(html, "B_p1_s1 (BLAST hits agree on Actiniaria; the best BLAST hit is from another phylum)",
    fixed = TRUE)
  expect_false(grepl("mp-fg-success", html, fixed = TRUE))
  expect_lt(regexpr("Submission-risk", html), regexpr("Final ITSx", html))
  smry$submission_risk <- ""
  expect_match(as.character(RiboPilot:::export_check_summary_ui(smry)), "mp-fg-success", fixed = TRUE)
})
