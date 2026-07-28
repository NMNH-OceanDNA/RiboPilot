# Tri-state triage. The rules are deliberately lopsided: a DROP needs positive
# evidence that the scaffold is not the target locus; anything merely suspicious
# becomes UNDECIDED, which still flows downstream.

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
# rank-aligned to FISH_T, in the 7-field form parse_rank_names() expects
FISH_L <- "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio"
COPE_T <- "33208;6656;6854;6843;1;2;3"
COPE_L <- "Metazoa;Arthropoda;Hexanauplia;Cyclopoida;Cyclopidae;Cyclops;Cyclops sp."

# ---- DB-backed behaviour -----------------------------------------------------
# Pure-rule behaviour (the rule registry itself) is owned by test-triage-rules.R.

mk_triage_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT,
    ignore INT, triage TEXT, triage_reason TEXT, triage_source TEXT, triage_rule TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, n_features INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL, median_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, species TEXT,
    lineage TEXT, lineage_taxids TEXT, rank_distance INT, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, sci_name TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, triage_opts TEXT)")
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio',?)", params = list(FISH_T))
  for (s in 1:3) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,?,5000,0,NULL,NULL,NULL,NULL)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,?,2,0,'default','default','no',1)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,?,100,100)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,?,'Danio rerio',?,?,0,1)",
      params = list(s, FISH_L, FISH_T))
  }
  con
}

# one more annotated, screened unit on S1
add_unit <- function(con, scaffold, sp, lin, tx, nuc = 1L) {
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,?,5000,0,NULL,NULL,NULL,NULL)",
    params = list(scaffold))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,?,2,0,'default','default','no',1)",
    params = list(scaffold))
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,?,100,100)", params = list(scaffold))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(scaffold))
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,?,?,?,?,0,?)",
    params = list(scaffold, sp, lin, tx, nuc))
}

test_that("unit_triage_state reports whether annotation actually ran", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::unit_triage_state(con, "S1")
  expect_true("annotated_n" %in% names(d))
  expect_false(any(is.na(d$annotated_n)))

  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL WHERE scaffold = 2")
  d2 <- RiboPilot:::unit_triage_state(con, "S1")
  expect_true(is.na(d2$annotated_n[2]))
})

test_that("a unit whose annotation never ran is never dropped", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # the exact shape that used to drop on a tool failure: no annotation rows and
  # a best hit that fails the nuclear-rDNA title regex
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule, ignore FROM assemblies WHERE scaffold = 1")
  expect_true(is.na(r$triage))
  expect_equal(r$triage_rule, "not_annotated")
  expect_equal(r$ignore, 0L)
})

test_that("the rule that fired is recorded alongside the reason", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0 WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule FROM assemblies WHERE scaffold = 1")
  expect_equal(r$triage, "drop")
  expect_equal(r$triage_rule, "no_feat_non_rdna")
})

test_that("auto-triage honours the sample's triage_opts preset", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0 WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")
  RiboPilot:::write_triage_opts(con, list(no_feat_non_rdna = "flag"))

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule FROM assemblies WHERE scaffold = 1")
  expect_true(is.na(r$triage))
  expect_equal(r$triage_rule, "no_feat_non_rdna")
})

test_that("preview reports per-rule counts alongside the totals", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens'")

  p <- triage_preview(con, "S1")
  expect_true(is.data.frame(p$by_rule))
  expect_equal(nrow(p$by_rule), 6L)
  expect_equal(sum(p$by_rule$n_decided), p$n_units)
  expect_true(all(p$by_rule$n_matching >= p$by_rule$n_decided))
  expect_equal(p$by_rule$n_decided[p$by_rule$key == "no_feat_non_rdna"], 3L)
})

test_that("a hand-decided unit is excluded from by_rule counts, but not from n_units", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens'")
  # scaffold 1 already has a curator's verdict; Apply will leave it alone
  set_unit_triage(con, "S1", 1, 1, "keep", source = "user")

  p <- triage_preview(con, "S1")
  expect_equal(p$n_units, 3L)              # still describes the whole selection
  expect_equal(p$n_user_skipped, 1L)
  # only the 2 units Apply would actually touch are counted per-rule
  expect_equal(p$by_rule$n_decided[p$by_rule$key == "no_feat_non_rdna"], 2L)
  expect_equal(sum(p$by_rule$n_decided), 2L)
})

test_that("preview and Apply agree on odd_one_out when a unit is skipped", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # sibling concordance is a property of the whole sample: three fish scaffolds
  # and two copepod ones have a strict species majority only while all five
  # vote. Counting on the unskipped subset alone loses it, and the preview then
  # reports a rule Apply goes on to fire.
  add_unit(con, 4, "Cyclops sp.", COPE_L, COPE_T)
  add_unit(con, 5, "Cyclops sp.", COPE_L, COPE_T)
  set_unit_triage(con, "S1", 1, 1, "keep", source = "user")

  p <- triage_preview(con, "S1")
  n_previewed <- p$by_rule$n_decided[p$by_rule$key == "odd_one_out"]
  expect_equal(n_previewed, 2L)

  auto_triage_sample(con, "S1")
  n_applied <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM assemblies WHERE triage_rule = 'odd_one_out'")$n
  expect_equal(n_applied, n_previewed)
})

test_that("the odd-one-out reason names the sibling taxon", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 4, "Cyclops sp.", COPE_L, COPE_T)

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con,
    "SELECT triage, triage_rule, triage_reason FROM assemblies WHERE scaffold = 4")
  expect_equal(r$triage_rule, "odd_one_out")
  expect_true(is.na(r$triage))
  expect_match(r$triage_reason, "Danio rerio at species", fixed = TRUE)
  expect_false(grepl("the majority", r$triage_reason, fixed = TRUE))
})

test_that("a unit with no annotate row at all is never dropped", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # not merely a NULL n_features: the whole row is absent
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "DELETE FROM annotate WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule, ignore FROM assemblies WHERE scaffold = 1")
  expect_true(is.na(r$triage))
  expect_equal(r$triage_rule, "not_annotated")
  expect_equal(r$ignore, 0L)
})

test_that("samples on different presets report a mixed action, not one sample's", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # a second sample, following a preset that disagrees with 'default' on one rule
  RiboPilot:::write_triage_opts(con, list(no_feat_non_rdna = "flag"), preset = "gentle")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S2','gentle')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Danio rerio',?)", params = list(FISH_T))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S2',1,1,5000,0,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S2',1,1,2,0,'default','default','no',1)")
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S2',1,1,100,100)")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S2',1,1,'18S')")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S2',1,1,'Danio rerio',?,?,0,1)",
    params = list(FISH_L, FISH_T))

  p <- triage_preview(con, c("S1", "S2"))
  # 'default' says drop, 'gentle' says flag: neither sample's action alone is honest
  expect_equal(p$by_rule$action[p$by_rule$key == "no_feat_non_rdna"], "mixed")
  # a rule neither preset customised must not be flattened to "mixed"
  expect_equal(p$by_rule$action[p$by_rule$key == "odd_one_out"], "flag")
})

test_that("set_unit_triage derives ignore and records who decided", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  set_unit_triage(con, "S1", 1, 2, "drop", reason = "test", source = "user")
  r <- DBI::dbGetQuery(con, "SELECT ignore, triage, triage_reason, triage_source FROM assemblies WHERE scaffold=2")
  expect_equal(r$ignore, 1L)
  expect_equal(r$triage, "drop")
  expect_equal(r$triage_source, "user")
  # held, not deleted
  expect_equal(DBI::dbGetQuery(con, "SELECT annotate_switch FROM annotate WHERE scaffold=2")$annotate_switch, 0L)

  set_unit_triage(con, "S1", 1, 2, "keep", source = "user")
  expect_equal(DBI::dbGetQuery(con, "SELECT ignore FROM assemblies WHERE scaffold=2")$ignore, 0L)

  # undecided is a real third state, and is INCLUDED downstream
  set_unit_triage(con, "S1", 1, 2, NA, reason = "unsure", source = "auto")
  r3 <- DBI::dbGetQuery(con, "SELECT ignore, triage FROM assemblies WHERE scaffold=2")
  expect_equal(r3$ignore, 0L)
  expect_true(is.na(r3$triage))

  expect_error(set_unit_triage(con, "S1", 1, 2, "banish"))
})

test_that("auto-triage never overwrites a human decision or a locked unit", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # make scaffold 1 droppable by the rules
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold=1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna=0, species='Homo sapiens' WHERE scaffold=1")
  # but the curator already said keep
  set_unit_triage(con, "S1", 1, 1, "keep", reason = "checked by hand", source = "user")
  # and scaffold 3 is annotation-locked
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock=1 WHERE scaffold=3")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT scaffold, triage, triage_source, ignore FROM assemblies ORDER BY scaffold")
  expect_equal(r$triage[1], "keep")                 # user decision survived
  expect_equal(r$triage_source[1], "user")
  expect_equal(r$ignore[1], 0L)
  expect_true(is.na(r$triage[3]) || r$triage_source[3] != "auto")  # locked, untouched

  # with overwrite_user the same run does re-decide it
  auto_triage_sample(con, "S1", overwrite_user = TRUE)
  r2 <- DBI::dbGetQuery(con, "SELECT triage, triage_source FROM assemblies WHERE scaffold=1")
  expect_equal(r2$triage, "drop")
  expect_equal(r2$triage_source, "auto")
})

test_that("auto-triage is idempotent", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  n1 <- auto_triage_sample(con, "S1")
  expect_true(n1 > 0)
  expect_equal(auto_triage_sample(con, "S1"), 0L)   # nothing left to change
})

test_that("preview counts the outcome and names samples that would be emptied", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # make every unit droppable
  DBI::dbExecute(con, "DELETE FROM annotations")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna=0, species='Homo sapiens'")
  p <- triage_preview(con, "S1")
  expect_equal(p$n_units, 3L)
  expect_equal(p$n_drop, 3L)
  expect_equal(p$emptied, "S1")      # the outcome that silently loses a sample
  expect_equal(p$n_samples, 1L)

  # a hand-kept unit is skipped and counted as such, and rescues the sample
  set_unit_triage(con, "S1", 1, 1, "keep", source = "user")
  p2 <- triage_preview(con, "S1")
  expect_equal(p2$n_user_skipped, 1L)
  expect_equal(p2$n_drop, 2L)
  expect_equal(p2$emptied, character(0))
})

test_that("the triage button emits full verdict words, not letters", {
  # A bare "k"/"d" used to resolve only through match.arg() partial matching,
  # and then surfaced to the user as "marked k".
  js <- as.character(RiboPilot:::rt_triage_bttn("scaf_triage"))
  expect_match(js, "'keep'", fixed = TRUE)
  expect_match(js, "'drop'", fixed = TRUE)
  expect_false(grepl("toLowerCase", js, fixed = TRUE))
})

test_that("set_unit_triage rejects abbreviations instead of expanding them", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  expect_error(set_unit_triage(con, "S1", 1, 1, "k"), "must be one of")
  expect_error(set_unit_triage(con, "S1", 1, 1, "dro"), "must be one of")
  expect_silent(set_unit_triage(con, "S1", 1, 1, "keep"))
  expect_silent(set_unit_triage(con, "S1", 1, 1, NA))
})
