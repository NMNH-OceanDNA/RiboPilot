# Settings must resolve identically for the app and for the WF1 triage process,
# and a hand-edited database must never be able to arm an unsupported drop.

mk_opts_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, triage_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  con
}

test_that("the seeded default resolves to the shipped actions", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::triage_opts_for(con, "S1"),
    RiboPilot:::triage_actions_default())
})

test_that("a sample follows its own preset", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_triage_opts(con,
    list(no_feat_non_rdna = "flag", odd_one_out = "off"), preset = "cautious")
  DBI::dbExecute(con, "UPDATE assemble SET triage_opts = 'cautious' WHERE ID = 'S1'")
  a <- RiboPilot:::triage_opts_for(con, "S1")
  expect_equal(a$no_feat_non_rdna, "flag")
  expect_equal(a$odd_one_out, "off")
  expect_equal(a$no_screen, "flag")
})

test_that("an unsupported drop written into the DB is read back as flag", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE triage_opts SET rule_odd_one_out = 'drop' WHERE triage_opts = 'default'")
  expect_equal(RiboPilot:::triage_opts_for(con, "S1")$odd_one_out, "flag")
})

test_that("a missing table, preset or sample falls back to the defaults", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::triage_opts_row(con, "nope"),
    RiboPilot:::triage_actions_default())
  expect_equal(RiboPilot:::triage_opts_for(con, "NOT_A_SAMPLE"),
    RiboPilot:::triage_actions_default())

  bare <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(bare), add = TRUE)
  expect_equal(RiboPilot:::triage_opts_for(bare, "S1"),
    RiboPilot:::triage_actions_default())
})

test_that("write then read round-trips every rule", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  want <- RiboPilot:::triage_actions_default()
  want$not_annotated <- "off"
  want$no_feat_rdna_hit <- "drop"
  RiboPilot:::write_triage_opts(con, want)
  expect_equal(RiboPilot:::triage_opts_row(con, "default"), want)
})
