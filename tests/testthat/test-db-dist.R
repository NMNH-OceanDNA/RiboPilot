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

test_that("parse_sha256_line reads shasum -a 256 output (single-space separator)", {
  expect_equal(parse_sha256_line(paste0(RP_DIGEST, " /tmp/x")), RP_DIGEST)
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

test_that("parse_sha256_line rejects a 65-char hex run rather than truncating it", {
  expect_true(is.na(parse_sha256_line(paste0(strrep("a", 65), "  /tmp/x"))))
})

test_that("parse_sha256_line rejects a 63-char hex run", {
  expect_true(is.na(parse_sha256_line(paste0(strrep("a", 63), "  /tmp/x"))))
})

test_that("parse_sha256_line still accepts exactly 64 hex chars", {
  expect_equal(parse_sha256_line(paste0(RP_DIGEST, "  /tmp/x")), RP_DIGEST)
})

test_that("parse_sha256_line trims leading whitespace before matching", {
  expect_equal(parse_sha256_line(paste0("\t", RP_DIGEST, "  /tmp/x")), RP_DIGEST)
})

# One release asset: a .tar.gz holding `files` (name -> contents) plus its .sha256
make_fake_asset <- function(tagdir, name, files) {
  src <- tempfile(); dir.create(src); on.exit(unlink(src, recursive = TRUE))
  for (f in names(files)) writeLines(files[[f]], file.path(src, f))
  tgz <- file.path(tagdir, name)
  withr::with_dir(src, utils::tar(tgz, files = list.files("."), compression = "gzip"))
  writeLines(paste0(rp_sha256(tgz), "  ", basename(tgz)), paste0(tgz, ".sha256"))
  tgz
}

make_fake_release <- function(root, version = "9.9.9", decoys = FALSE) {
  # assets live under a tag dir, mirroring the real releases/download/<tag>/<asset> layout
  rel <- file.path(root, "rel")
  tagdir <- file.path(rel, paste0("db-", version)); dir.create(tagdir, recursive = TRUE)
  make_fake_asset(tagdir, sprintf("ribopilot_rrna-%s.tar.gz", version), list(
    manifest.json = sprintf(
      '{"manifest_schema":1,"db_name":"ribopilot_rrna","db_version":"%s"}', version),
    ribopilot_rrna.nin = "x", ribopilot_rrna.nsq = "x", ribopilot_rrna.nhr = "x",
    ribopilot_rrna.ndb = "x"))
  if (isTRUE(decoys)) {
    make_fake_asset(tagdir, sprintf("decoys-%s.tar.gz", version), list(decoy.nin = "x"))
  }
  paste0("file://", rel)
}

# End users do not build or install panels: the shipped container carries one.
# Kept internal on purpose, so a stray devtools::document() must not re-export it.
test_that("install_ribopilot_db is internal, not part of the public surface", {
  expect_false("install_ribopilot_db" %in% getNamespaceExports("RiboPilot"))
  expect_true(is.function(RiboPilot:::install_ribopilot_db))
})

test_that("install_ribopilot_db extracts a verified panel and returns its prefix", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs")
  p <- install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)
  expect_equal(p, file.path(dest, "9.9.9", "ribopilot_rrna"))
  expect_true(blast_db_exists(p))
  expect_equal(read_panel_manifest(dirname(p))$db_version, "9.9.9")
})

test_that("install_ribopilot_db refuses a tarball whose checksum does not match", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  sums <- list.files(file.path(root, "rel"), pattern = "\\.sha256$", full.names = TRUE,
    recursive = TRUE)
  writeLines(paste0(strrep("0", 64), "  ", sub("\\.sha256$", "", basename(sums[1]))), sums[1])
  dest <- file.path(root, "dbs")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base),
    "[Cc]hecksum")
  # a failed verification must leave nothing behind to be mistaken for a panel
  expect_false(blast_db_exists(file.path(dest, "9.9.9", "ribopilot_rrna")))
})

test_that("install_ribopilot_db is a no-op when the panel is already present", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs")
  p <- install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)
  unlink(file.path(root, "rel"), recursive = TRUE)   # make any refetch impossible
  expect_equal(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base), p)
})

test_that("install_ribopilot_db reads an openssl-shaped checksum file", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  tgz <- file.path(root, "rel", "db-9.9.9", "ribopilot_rrna-9.9.9.tar.gz")
  # the path holds a 64-hex component, which an unanchored search would take
  writeLines(sprintf("SHA2-256(/b/%s/%s)= %s", strrep("a", 64), basename(tgz), rp_sha256(tgz)),
    paste0(tgz, ".sha256"))
  dest <- file.path(root, "dbs")
  expect_true(blast_db_exists(
    install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)))
})

test_that("install_ribopilot_db leaves no panel when a later asset fails to verify", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root, decoys = TRUE)
  # the panel tarball verifies and unpacks; the decoys tarball then fails
  sum_file <- file.path(root, "rel", "db-9.9.9", "decoys-9.9.9.tar.gz.sha256")
  writeLines(paste0(strrep("0", 64), "  ", "decoys-9.9.9.tar.gz"), sum_file)
  dest <- file.path(root, "dbs")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, base_url = base), "[Cc]hecksum")
  expect_false(blast_db_exists(file.path(dest, "9.9.9", "ribopilot_rrna")))
  expect_false(dir.exists(file.path(dest, "9.9.9")))
})

test_that("install_ribopilot_db leaves no panel when a later asset is missing", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)   # the panel is published, the decoys are not
  dest <- file.path(root, "dbs")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, base_url = base),
    "Could not download")
  # the underlying reason has to survive, or a TLS or proxy failure is unreadable
  expect_error(install_ribopilot_db("9.9.9", dest = dest, base_url = base), "cannot open")
  expect_false(blast_db_exists(file.path(dest, "9.9.9", "ribopilot_rrna")))
  expect_false(dir.exists(file.path(dest, "9.9.9")))
})

test_that("install_ribopilot_db rejects an archive holding only the index file", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  tagdir <- file.path(root, "rel", "db-9.9.9"); dir.create(tagdir, recursive = TRUE)
  make_fake_asset(tagdir, "ribopilot_rrna-9.9.9.tar.gz", list(ribopilot_rrna.nin = "x"))
  dest <- file.path(root, "dbs")
  base <- paste0("file://", file.path(root, "rel"))
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base),
    "did not contain a usable panel")
  # naming the absentees is what makes a future .nal-volume panel diagnosable
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base),
    "ribopilot_rrna\\.nsq, ribopilot_rrna\\.nhr, ribopilot_rrna\\.ndb, manifest\\.json")
  expect_false(blast_db_exists(file.path(dest, "9.9.9", "ribopilot_rrna")))
})

test_that("install_ribopilot_db refuses an archive it cannot unpack", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  tagdir <- file.path(root, "rel", "db-9.9.9"); dir.create(tagdir, recursive = TRUE)
  # a gzip header with no body: the checksum verifies, the unpack cannot
  tgz <- file.path(tagdir, "ribopilot_rrna-9.9.9.tar.gz")
  writeBin(as.raw(c(0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00)), tgz)
  writeLines(paste0(rp_sha256(tgz), "  ", basename(tgz)), paste0(tgz, ".sha256"))
  dest <- file.path(root, "dbs")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE,
    base_url = paste0("file://", file.path(root, "rel"))), "Could not unpack")
  expect_false(blast_db_exists(file.path(dest, "9.9.9", "ribopilot_rrna")))
})

test_that("a failed install leaves no staging directory behind", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)   # the panel is published, the decoys are not
  dest <- file.path(root, "dbs")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, base_url = base),
    "Could not download")
  # a stage can hold a few hundred MB, so a failure must not orphan one
  left <- list.files(dest, all.files = TRUE, no.. = TRUE)
  expect_false(any(startsWith(left, ".stage-")))
})

test_that("a failed reinstall leaves the installed panel whole, not re-extracted", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root, decoys = TRUE)
  dest <- file.path(root, "dbs")
  p <- install_ribopilot_db("9.9.9", dest = dest, base_url = base)
  # republish the panel asset with a different manifest, then break the decoys:
  # the first asset verifies and unpacks, so only staging keeps it out of the panel
  tagdir <- file.path(root, "rel", "db-9.9.9")
  make_fake_asset(tagdir, "ribopilot_rrna-9.9.9.tar.gz", list(
    manifest.json =
      '{"manifest_schema":1,"db_name":"ribopilot_rrna","db_version":"9.9.9-new"}',
    ribopilot_rrna.nin = "y", ribopilot_rrna.nsq = "y", ribopilot_rrna.nhr = "y"))
  writeLines(paste0(strrep("0", 64), "  decoys-9.9.9.tar.gz"),
    file.path(tagdir, "decoys-9.9.9.tar.gz.sha256"))
  expect_error(install_ribopilot_db("9.9.9", dest = dest, base_url = base, force = TRUE),
    "[Cc]hecksum")
  expect_true(blast_db_exists(p))
  expect_equal(read_panel_manifest(dirname(p))$db_version, "9.9.9")
})

test_that("a set-aside directory from a crash is never the target of a new install", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs"); dir.create(dest, recursive = TRUE)
  # a survivor of an earlier crash, named for a pid the OS has since recycled
  survivor <- file.path(dest, paste0("9.9.9.old-", Sys.getpid()))
  dir.create(survivor); writeLines("only copy", file.path(survivor, "ribopilot_rrna.nsq"))
  install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)
  expect_true(file.exists(file.path(survivor, "ribopilot_rrna.nsq")))
})

test_that("install sweeps leftovers from interrupted runs but spares fresh ones", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs"); dir.create(dest, recursive = TRUE)
  stale <- file.path(dest, c(".stage-9.9.9-4242", "9.9.9.old-4242"))
  fresh <- file.path(dest, c(".stage-9.9.9-4243", ".stage-8.8.8-4242"))
  for (d in c(stale, fresh)) dir.create(d)
  for (d in stale) Sys.setFileTime(d, Sys.time() - 3 * 86400)
  expect_message(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base),
    "Reclaiming space")
  expect_false(any(dir.exists(stale)))
  expect_true(all(dir.exists(fresh)))   # a live concurrent install, and another version
})

test_that("install_ribopilot_db with force=TRUE replaces what is there", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs")
  p <- install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)
  writeLines("stale", file.path(dirname(p), "leftover.txt"))
  expect_equal(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base,
    force = TRUE), p)
  expect_true(blast_db_exists(p))
  expect_false(file.exists(file.path(dirname(p), "leftover.txt")))
})

test_that("a failed re-install leaves the installed panel intact", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  base <- make_fake_release(root)
  dest <- file.path(root, "dbs")
  p <- install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base)
  unlink(file.path(root, "rel"), recursive = TRUE)
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE, base_url = base,
    force = TRUE), "Could not download")
  expect_true(blast_db_exists(p))
  expect_equal(read_panel_manifest(dirname(p))$db_version, "9.9.9")
})

test_that("install_ribopilot_db refuses a destination holding a space", {
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  dest <- file.path(root, "panel dir")
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE,
    base_url = "file:///nowhere"), "whitespace")
  expect_false(dir.exists(dest))   # refused before anything is created or fetched
})

test_that("install_ribopilot_db refuses a destination it cannot write", {
  skip_on_os("windows")
  skip_if(unname(Sys.info()[["user"]]) == "root", "root ignores directory permissions")
  root <- tempfile(); dir.create(root); on.exit(unlink(root, recursive = TRUE))
  dest <- file.path(root, "readonly"); dir.create(dest)
  Sys.chmod(dest, "0555")
  on.exit(Sys.chmod(dest, "0755"), add = TRUE, after = FALSE)
  expect_error(install_ribopilot_db("9.9.9", dest = dest, decoys = FALSE,
    base_url = "file:///nowhere"), "not writable")
})

mk_prov_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  con
}

test_that("record_blastdb_provenance stores the manifest version beside the panel", {
  con <- mk_prov_db(); on.exit(DBI::dbDisconnect(con))
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(d, "manifest.json"))
  file.create(file.path(d, "ribopilot_rrna.nin"))
  expect_equal(record_blastdb_provenance(con, file.path(d, "ribopilot_rrna")), "2026.08.1")
  expect_equal(blastdb_provenance(con), "ribopilot_rrna 2026.08.1")
})

test_that("an unset or unmanifested panel records nothing and reports empty", {
  con <- mk_prov_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(record_blastdb_provenance(con, ""), "")
  expect_equal(blastdb_provenance(con), "")
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  file.create(file.path(d, "ribopilot_rrna.nin"))   # panel, no manifest
  expect_equal(record_blastdb_provenance(con, file.path(d, "ribopilot_rrna")), "")
})

test_that("provenance can be recorded from a discovered panel", {
  con <- mk_prov_db(); on.exit(DBI::dbDisconnect(con))
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE), add = TRUE)
  d <- file.path(r, "2026.08.1"); dir.create(d)
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", file.path(d, paste0("ribopilot_rrna.", e)))
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}', file.path(d, "manifest.json"))
  p <- discover_panel(r)
  expect_equal(record_blastdb_provenance(con, p), "2026.08.1")
  expect_equal(blastdb_provenance(con), "ribopilot_rrna 2026.08.1")
})

# rp_searched.R's DB fixture: everything derive_blast_products() and
# mark_unit_searched() touch, plus reference_config for the provenance write.
mk_searched_script_db <- function(db_path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, blast_accession TEXT,
    blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT, ignore INT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,NULL,NULL,NULL,NULL,NULL,NULL,0)")
  DBI::dbDisconnect(con)
}

run_searched_script <- function(script, src, db_path, panel) {
  err <- withr::local_tempfile()
  status <- system2("Rscript",
    c(shQuote(script), shQuote(db_path), "S1", "1", "1", "abc", "0", "local", shQuote(panel)),
    env = paste0("RIBOPILOT_SRC=", src), stdout = FALSE, stderr = err)
  list(status = status, stderr = paste(readLines(err, warn = FALSE), collapse = "\n"))
}

test_that("rp_searched.R records nothing for a panel prefix with no manifest, without erroring", {
  skip_on_cran()
  skip_if_not_installed("pkgload")
  skip_if_not_installed("withr")
  src <- RiboPilot:::ribopilot_src_default()
  # An INSTALLED package directory also has a DESCRIPTION and an R/, so under
  # R CMD check this ran with a src that rp_env.R rightly refuses, and the
  # refusal notice on stderr failed the assertion below. Same predicate as
  # rp_env.R's rp_is_source_dir(): .R files under R/ are what tell them apart.
  skip_if(!file.exists(file.path(src, "DESCRIPTION")) ||
      length(list.files(file.path(src, "R"), pattern = "\\.[Rr]$")) == 0,
    "not a source checkout")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src, "inst", "nextflow", "bin")
  script <- file.path(bin_dir, "rp_searched.R")

  panel_dir <- withr::local_tempdir()  # no manifest.json beside this panel
  panel <- file.path(panel_dir, "ribopilot_rrna")

  db_path <- withr::local_tempfile(fileext = ".sqlite")
  mk_searched_script_db(db_path)

  res <- run_searched_script(script, src, db_path, panel)
  expect_equal(res$status, 0, info = res$stderr)
  expect_equal(res$stderr, "")

  con2 <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con2), add = TRUE)
  expect_equal(blastdb_provenance(con2), "")
})

test_that("provenance reads degrade to empty on a project predating the column", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE reference_config (rowid_key INTEGER PRIMARY KEY,
    ref_dir TEXT, min_records INTEGER, max_records INTEGER, blastdb TEXT)")
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  expect_equal(blastdb_provenance(con), "")
  expect_silent(record_blastdb_provenance(con, ""))

  # a genuine prefix with a readable manifest reaches the UPDATE itself, which
  # is the write-side half of the same guard
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(d, "manifest.json"))
  file.create(file.path(d, "ribopilot_rrna.nin"))
  expect_silent(record_blastdb_provenance(con, file.path(d, "ribopilot_rrna")))
  expect_equal(blastdb_provenance(con), "")
})
