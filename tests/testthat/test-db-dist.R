# Downloading a reference panel, and knowing which one it is.

test_that("rp_sha256 matches a known digest", {
  f <- tempfile(); on.exit(unlink(f))
  writeBin(charToRaw("abc"), f)
  # sha256("abc"), a fixed published vector
  expect_equal(rp_sha256(f),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
})

test_that("rp_sha256 returns NA for a missing file", {
  expect_true(is.na(rp_sha256(file.path(tempdir(), "definitely-absent-xyz"))))
})

test_that("read_panel_manifest reads db_version and db_name", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  writeLines('{"manifest_schema":1,"db_name":"ribopilot_rrna","db_version":"2026.08.1",
    "counts":{"records_after_dedup":597988}}', file.path(d, "manifest.json"))
  m <- read_panel_manifest(d)
  expect_equal(m$db_version, "2026.08.1")
  expect_equal(m$db_name, "ribopilot_rrna")
  expect_equal(m$counts$records_after_dedup, 597988)
})

test_that("read_panel_manifest returns NULL when absent or unusable", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  expect_null(read_panel_manifest(d))
  writeLines("not json at all {{{", file.path(d, "manifest.json"))
  expect_null(read_panel_manifest(d))
  writeLines('{"db_name":"x"}', file.path(d, "manifest.json"))
  expect_null(read_panel_manifest(d))  # no db_version
})
