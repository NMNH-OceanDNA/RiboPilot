geome_project <- function(bcids = c("ark:/21547/CYB2REEDY", ""), envir = parent.frame()) {
  m <- meta_reads(c("s1", "s2"))
  m$GEOME_BCID <- bcids
  meta_project(m, envir = envir)
}

geome_up <- function(p, up, ...) {
  f <- file.path(p$dir, "up.csv")
  up <- cbind(up, meta_reads(up$ID)[c("R1", "R2")])
  utils::write.csv(up, f, row.names = FALSE)
  suppressMessages(update_sample_metadata(p$proj, f, ...))
}

test_that("add_samples stores and fetches BCIDs for new samples only", {
  local_mocked_bindings(.geome_get = geome_fixture_get, url_status = function(u) "ok")
  p <- geome_project()
  up <- meta_reads("s3")
  up$MyBcid <- "https://n2t.net/ark:/21547/CYA2Reedy01"
  file.create(file.path(p$data_dir, c(up$R1, up$R2)))
  f <- file.path(p$dir, "add.csv")
  utils::write.csv(up, f, row.names = FALSE)
  suppressMessages(add_samples(p$proj, f, mapping_geome = "MyBcid"))
  expect_equal(meta_q(p$proj, "SELECT GEOME_BCID FROM samples WHERE ID='s3'")$GEOME_BCID,
    "ark:/21547/CYA2Reedy01")
  expect_false("MyBcid" %in% names(meta_q(p$proj, "SELECT * FROM samples")))
  expect_setequal(meta_q(p$proj, "SELECT ID FROM meta_status")$ID, c("s1", "s3"))
})

test_that("update_sample_metadata refetches changed, drops cleared, skips unchanged", {
  local_mocked_bindings(.geome_get = geome_fixture_get, url_status = function(u) "ok")
  p <- geome_project(bcids = c("ark:/21547/CYB2REEDY", "ark:/21547/CYA2Reedy01"))
  geome_reset_calls()
  up <- data.frame(ID = c("s1", "s2"), Taxon = "Fish a",
    GEOME_BCID = c("ark:/21547/CYB2REEDY", ""))
  geome_up(p, up)
  expect_equal(meta_q(p$proj, "SELECT ID FROM meta_status")$ID, "s1")
  expect_equal(meta_q(p$proj, "SELECT COUNT(*) n FROM meta_records WHERE ID='s2'")$n, 0L)
  expect_equal(length(ls(geome_calls)), 0L)
  up$GEOME_BCID <- c("ark:/21547/CYC2CMPI38181.1", "")
  geome_up(p, up)
  lv <- meta_q(p$proj, "SELECT DISTINCT level FROM meta_records WHERE ID='s1'")$level
  expect_true("Tissue" %in% lv)
})

test_that("update_sample_metadata with fetch_geome = FALSE drops stale records for a changed BCID", {
  local_mocked_bindings(.geome_get = geome_fixture_get, url_status = function(u) "ok")
  p <- geome_project()
  up <- data.frame(ID = c("s1", "s2"), Taxon = "Fish a",
    GEOME_BCID = c("ark:/21547/CYA2Reedy01", ""))
  geome_up(p, up, fetch_geome = FALSE)
  expect_equal(meta_q(p$proj, "SELECT COUNT(*) n FROM meta_records WHERE ID='s1'")$n, 0L)
  expect_equal(meta_q(p$proj, "SELECT COUNT(*) n FROM meta_status WHERE ID='s1'")$n, 0L)
})

test_that("update CSV without a BCID column leaves GEOME data alone", {
  local_mocked_bindings(.geome_get = geome_fixture_get, url_status = function(u) "ok")
  p <- geome_project()
  geome_up(p, data.frame(ID = "s1", Taxon = "Fish a"))
  expect_equal(meta_q(p$proj, "SELECT ID FROM meta_status")$ID, "s1")
  expect_equal(meta_q(p$proj, "SELECT GEOME_BCID FROM samples WHERE ID='s1'")$GEOME_BCID,
    "ark:/21547/CYB2REEDY")
})

test_that("a Taxon change still re-queues taxonomy while BCIDs are synced", {
  local_mocked_bindings(.geome_get = geome_fixture_get, url_status = function(u) "ok")
  p <- geome_project()
  geome_up(p, data.frame(ID = "s1", Taxon = "Fish b", GEOME_BCID = "ark:/21547/CYB2REEDY"))
  expect_equal(meta_q(p$proj, "SELECT taxonomy_switch, taxonomy_lock FROM samples WHERE ID='s1'")$taxonomy_switch, 1L)
  expect_equal(meta_q(p$proj, "SELECT ID FROM meta_status")$ID, "s1")
})
