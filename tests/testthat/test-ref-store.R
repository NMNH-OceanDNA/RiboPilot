# The shared store is read by many projects and, on a cluster, by many users.
# Reading must not require write access.

mk_ro_db <- function(dir) {
  p <- file.path(dir, "ro.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "CREATE TABLE t (a TEXT)")
  DBI::dbExecute(con, "INSERT INTO t VALUES ('x')")
  DBI::dbDisconnect(con)
  Sys.chmod(p, "0444")
  p
}

test_that("open_sqlite_ro reads a write-protected database", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT a FROM t")$a, "x")
})

test_that("open_sqlite_ro refuses to write", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_error(DBI::dbExecute(con, "INSERT INTO t VALUES ('y')"))
})

test_that("open_sqlite_ro errors on a missing file rather than creating one", {
  d <- withr::local_tempdir()
  p <- file.path(d, "absent.sqlite")
  expect_error(RiboPilot:::open_sqlite_ro(p))
  expect_false(file.exists(p))
})
