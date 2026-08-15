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

test_that("ribopilot_db_dir honours RIBOPILOT_DB_DIR", {
  withr::with_envvar(c(RIBOPILOT_DB_DIR = "/custom/panels"), {
    expect_equal(ribopilot_db_dir(), "/custom/panels")
  })
})

test_that("ribopilot_db_dir falls back to blastdb under the ref dir", {
  withr::with_envvar(c(RIBOPILOT_DB_DIR = "", RIBOPILOT_REF_DIR = "/refs"), {
    expect_equal(ribopilot_db_dir(), file.path("/refs", "blastdb"))
  })
})

RP_DIGEST <- "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

test_that("parse_sha256_line reads sha256sum output", {
  expect_equal(parse_sha256_line(paste0(RP_DIGEST, "  /tmp/x")), RP_DIGEST)
})

test_that("parse_sha256_line reads shasum -a 256 output", {
  expect_equal(parse_sha256_line(paste0(RP_DIGEST, "  /tmp/x")), RP_DIGEST)
})

test_that("parse_sha256_line reads openssl 3.x output", {
  expect_equal(parse_sha256_line(paste0("SHA2-256(/tmp/x)= ", RP_DIGEST)), RP_DIGEST)
})

test_that("parse_sha256_line reads openssl 1.x output", {
  expect_equal(parse_sha256_line(paste0("SHA256(/tmp/x)= ", RP_DIGEST)), RP_DIGEST)
})

test_that("parse_sha256_line is not fooled by a hex-looking path component", {
  line <- paste0("SHA2-256(/tmp/", strrep("a", 64), "/x)= ", RP_DIGEST)
  expect_equal(parse_sha256_line(line), RP_DIGEST)
})

test_that("parse_sha256_line returns NA for unrecognised output", {
  expect_true(is.na(parse_sha256_line("not a digest line")))
})
