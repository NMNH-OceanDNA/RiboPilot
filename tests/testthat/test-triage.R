# Tri-state triage. The rules are deliberately lopsided: a DROP needs positive
# evidence that the scaffold is not the target locus; anything merely suspicious
# becomes UNDECIDED, which still flows downstream.

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
COPE_T <- "33208;6656;6854;6843;1;2;3"

st <- function(...) {
  rows <- list(...)
  d <- do.call(rbind, lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    data.frame(path = 1L, scaffold = i, length = r$length %||% 5000L,
      ignore = 0L, triage = r$triage %||% NA_character_,
      triage_reason = NA_character_, triage_source = r$source %||% NA_character_,
      depth = r$depth %||% 100, screen_species = r$sp %||% NA_character_,
      screen_taxids = r$tx %||% NA_character_, rank_distance = r$rd %||% NA_integer_,
      is_nuclear_rrna = r$nuc %||% NA_integer_, sample_taxids = FISH_T,
      sci_name = "Danio rerio", n_features = r$feats %||% 5L,
      annotate_lock = r$lock %||% 0L, reviewed = "no", stringsAsFactors = FALSE)
  }))
  d
}

test_that("a drop needs positive evidence on BOTH available channels", {
  # no features AND the best core_nt match is not nuclear rDNA -> drop
  d <- RiboPilot:::triage_rules(st(list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T)))
  expect_equal(d$new_triage, "drop")
  expect_match(d$new_reason, "no rRNA or ITS features")
  expect_match(d$new_reason, "not a nuclear rDNA record")

  # no features alone is NOT enough: annotation may simply not have run
  d2 <- RiboPilot:::triage_rules(st(list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T)))
  expect_true(is.na(d2$new_triage))
  expect_match(d2$new_reason, "no rRNA or ITS features")

  # a non-rDNA top hit alone is NOT enough: the scaffold still has a feature model
  d3 <- RiboPilot:::triage_rules(st(list(feats = 5, nuc = 0, sp = "Homo sapiens", tx = FISH_T)))
  expect_true(is.na(d3$new_triage))
})

test_that("taxonomic distance alone never drops a scaffold", {
  # The failure mode that would make the tool useless for its target taxa: an
  # under-represented organism whose nearest reference is in another phylum.
  d <- RiboPilot:::triage_rules(st(
    list(feats = 5, nuc = 1, sp = "Distant relative", tx = COPE_T, rd = 6L)))
  expect_false(identical(d$new_triage, "drop"))
})

test_that("a clean scaffold is kept with no reason attached", {
  d <- RiboPilot:::triage_rules(st(list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T)))
  expect_equal(d$new_triage, "keep")
  expect_true(is.na(d$new_reason))
})

test_that("shallow and odd-one-out scaffolds are flagged, not dropped", {
  # 3 fish + 1 copepod, all annotated: the copepod is the odd one out
  d <- RiboPilot:::triage_rules(st(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Cyclops sp.", tx = COPE_T)))
  expect_equal(d$new_triage[1:3], rep("keep", 3))
  expect_true(is.na(d$new_triage[4]))
  expect_match(d$new_reason[4], "differs from this sample's other scaffolds")

  # depth well below the sample max is a flag with a legible reason
  s <- RiboPilot:::triage_rules(st(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T, depth = 200),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T, depth = 4)))
  expect_equal(s$new_triage[1], "keep")
  expect_true(is.na(s$new_triage[2]))
  expect_match(s$new_reason[2], "2% of this sample's deepest scaffold")
})

test_that("a missing GenBank screen is undecided, never a drop", {
  d <- RiboPilot:::triage_rules(st(list(feats = 5, nuc = NA, sp = NA, tx = NA)))
  expect_true(is.na(d$new_triage))
  expect_match(d$new_reason, "screen has not run")
})

# ---- DB-backed behaviour -----------------------------------------------------

mk_triage_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT,
    ignore INT, triage TEXT, triage_reason TEXT, triage_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL, median_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, species TEXT,
    lineage_taxids TEXT, rank_distance INT, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, sci_name TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio',?)", params = list(FISH_T))
  for (s in 1:3) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,?,5000,0,NULL,NULL,NULL)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,?,2,0,'default','default','no')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,?,100,100)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,?,'Danio rerio',?,0,1)",
      params = list(s, FISH_T))
  }
  con
}

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
