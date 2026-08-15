# Stage 4b: Database Distribution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the local rDNA panel obtainable with one command and configurable from the app, so a user never hand-places a BLAST database or hand-edits a config row.

**Architecture:** A versioned directory holds the panel volumes and the NCBI decoys side by side, exactly where Stage 4a's `rp_dbarg.R` already looks for them. A download helper fetches a checksummed tarball from a GitHub release and extracts it there. The Reference settings modal gains a panel-path control that validates read-only, because the panel may live on a read-only shared mount. Build-side shell scripts fetch the decoys and publish the release.

**Tech Stack:** R (>= 4.4.0), jsonlite, DBI/RSQLite, shiny, testthat 3e, BLAST+ 2.17.0, `update_blastdb.pl`, GitHub CLI.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md` (Stage 4 section)

## Global Constraints

- **ASCII only in R code.** No non-ASCII characters anywhere in `R/`, `tests/`, or `inst/`.
- **Minimal comments.** Comment why, not what. Never narrate a bugfix in a comment.
- **No schema migrations.** This is pre-release. Never write a migration for existing projects. New columns go in the DDL; every read of a possibly-missing column is wrapped so an older project degrades to a default instead of erroring. Follow the existing `nf_blastdb()` pattern (`R/nf_command.R:85-93`).
- **R floor is 4.4.0.** `tools::sha256sum()` was introduced in R 4.5.0 and MUST NOT be called unconditionally.
- **Do not add package dependencies.** `jsonlite` is already in Imports; `digest` and `openssl` are not, and must not be added.
- **Never push, never create a release from a task.** Task 6 writes scripts; it does not run them against the real repo.
- **The panel may be read-only.** Never call `validate_ref_dir()` on a panel path: it requires write access and creates directories.
- **Never name anything `manifest.tsv` or reuse `read_ref_manifest()`.** That is the clade-database manifest (`R/reference_dir.R:88`), a different file for a different purpose. The panel manifest is `manifest.json`.
- Tests that shell out to network or BLAST binaries must be guarded with `testthat::skip_on_cran()` and a `Sys.which()` guard.

## Layout This Plan Establishes

Stage 4a's `inst/nextflow/bin/rp_dbarg.R` defaults the decoy directory to `dirname(panel_prefix)`. Decoys therefore live **beside** the panel, in the same version directory. Do not invent a separate decoy folder.

```
<db_dir>/<db_version>/
  manifest.json
  ribopilot_rrna.{ndb,nhr,nin,njs,nog,nos,not,nsq,ntf,nto}
  16S_ribosomal_RNA.*
  SSU_eukaryote_rRNA.*
  LSU_eukaryote_rRNA.*
  taxdb.bt{d,i}
```

The configured `reference_config.blastdb` is the prefix `<db_dir>/<db_version>/ribopilot_rrna`.

Measured sizes (real, from the built 2026.08.1 panel): volumes + manifest are 332 MB unpacked, **75.7 MB gzipped**. The `.fasta` (770 MB), `.taxid.map` (11.6 MB), and `.dropped.tsv` (15.6 MB) are build artifacts and are NOT shipped.

## File Structure

- Create `R/db_dist.R` — checksum, panel manifest, panel location, install helper. One responsibility: getting a panel onto disk and knowing what it is.
- Modify `R/reference_dir.R` — `set_ref_config()` gains `blastdb`; add `validate_blastdb_prefix()`.
- Modify `R/app_reference.R` — panel control in the settings modal and its observer.
- Modify `R/db_schema.R` — one provenance column.
- Modify `R/nf_command.R` — record provenance when a run starts.
- Create `dev/refdb/fetch_decoys.sh`, `dev/refdb/package_release.sh`.
- Create `tests/testthat/test-db-dist.R`; extend `test-reference-dir.R`, `test-reference-app.R`.

---

### Task 1: Checksum and panel manifest

**Files:**
- Create: `R/db_dist.R`
- Test: `tests/testthat/test-db-dist.R`

**Interfaces:**
- Consumes: nothing.
- Produces: `rp_sha256(path)` returning a lowercase 64-char hex string or `NA_character_`; `read_panel_manifest(dir)` returning a named list or `NULL`.

- [ ] **Step 1: Write the failing tests**

```r
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
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-db-dist.R")'`
Expected: FAIL, "could not find function rp_sha256".

- [ ] **Step 3: Implement**

```r
# Getting a reference panel onto disk, and knowing which one it is.

#' SHA-256 of a file, as lowercase hex
#'
#' tools::sha256sum() arrived in R 4.5.0 and the package floor is 4.4.0, so it
#' is used only when present and a command-line digest is the fallback.
#' @return 64-char lowercase hex, or NA_character_ when unobtainable.
#' @noRd
rp_sha256 <- function(path) {
  if (length(path) != 1 || is.na(path) || !file.exists(path)) return(NA_character_)
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    v <- tryCatch(unname(get("sha256sum", envir = asNamespace("tools"))(path)),
      error = function(e) NA_character_)
    if (length(v) == 1 && !is.na(v)) return(tolower(v))
  }
  cand <- list(c("sha256sum", path), c("shasum", "-a", "256", path),
    c("openssl", "dgst", "-sha256", path))
  for (cc in cand) {
    if (!nzchar(Sys.which(cc[1]))) next
    out <- tryCatch(suppressWarnings(
      system2(cc[1], shQuote(cc[-1]), stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0))
    if (length(out) < 1) next
    hit <- regmatches(out[1], regexpr("[0-9a-fA-F]{64}", out[1]))
    if (length(hit) == 1) return(tolower(hit))
  }
  NA_character_
}

#' Read a panel directory's manifest.json
#'
#' Distinct from read_ref_manifest(), which indexes built clade databases.
#' @return Named list, or NULL when there is no usable manifest.
#' @noRd
read_panel_manifest <- function(dir) {
  if (length(dir) != 1 || is.na(dir)) return(NULL)
  p <- file.path(path.expand(dir), "manifest.json")
  if (!file.exists(p)) return(NULL)
  m <- tryCatch(jsonlite::fromJSON(p, simplifyVector = TRUE), error = function(e) NULL)
  if (!is.list(m)) return(NULL)
  ok <- function(x) length(x) == 1 && !is.na(x) && nzchar(as.character(x))
  if (!ok(m$db_version) || !ok(m$db_name)) return(NULL)
  m
}
```

- [ ] **Step 4: Run to verify pass**

Run the same command. Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/db_dist.R tests/testthat/test-db-dist.R
git commit -m "feat: add file checksum and panel manifest reader"
```

---

### Task 2: Panel location and read-only validation

**Files:**
- Modify: `R/db_dist.R`
- Modify: `R/reference_dir.R` (add `validate_blastdb_prefix()` near `validate_ref_dir()`)
- Test: `tests/testthat/test-db-dist.R`, `tests/testthat/test-reference-dir.R`

**Interfaces:**
- Consumes: `blast_db_exists(prefix)` from `R/blast_local.R`; `ribopilot_ref_dir()` from `R/build_reference_db.R`.
- Produces: `ribopilot_db_dir()` returning a path; `validate_blastdb_prefix(prefix)` returning `list(ok, msg)` where `ok=TRUE, msg=""` means "unset, use remote".

- [ ] **Step 1: Write the failing tests**

```r
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

test_that("an empty panel prefix is valid and means remote", {
  r <- validate_blastdb_prefix("")
  expect_true(r$ok); expect_equal(r$msg, "")
  expect_true(validate_blastdb_prefix(NULL)$ok)
})

test_that("a prefix with no volume files is rejected", {
  r <- validate_blastdb_prefix(file.path(tempdir(), "nope", "ribopilot_rrna"))
  expect_false(r$ok)
  expect_match(r$msg, "No BLAST database")
})

test_that("a prefix with volume files is accepted without needing write access", {
  d <- tempfile(); dir.create(d)
  on.exit({Sys.chmod(d, "0755"); unlink(d, recursive = TRUE)})
  file.create(file.path(d, "ribopilot_rrna.nin"))
  Sys.chmod(d, "0555")  # read-only, as on a shared HPC mount
  skip_if(unname(file.access(d, 2)) == 0, "cannot drop write permission here")
  r <- validate_blastdb_prefix(file.path(d, "ribopilot_rrna"))
  expect_true(r$ok)
  expect_equal(r$msg, file.path(d, "ribopilot_rrna"))
})

test_that("a whitespace prefix is rejected, since blastn cannot express it", {
  d <- tempfile("has space"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  file.create(file.path(d, "ribopilot_rrna.nin"))
  r <- validate_blastdb_prefix(file.path(d, "ribopilot_rrna"))
  expect_false(r$ok)
  expect_match(r$msg, "space")
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-db-dist.R")'`
Expected: FAIL, "could not find function ribopilot_db_dir".

- [ ] **Step 3: Implement**

Append to `R/db_dist.R`:

```r
#' Where downloaded BLAST panels live
#'
#' Honours RIBOPILOT_DB_DIR so a cluster user can put a multi-hundred-MB panel
#' off their home quota, then falls back to a subfolder of the reference dir.
#' @noRd
ribopilot_db_dir <- function() {
  d <- Sys.getenv("RIBOPILOT_DB_DIR", "")
  if (!nzchar(d)) d <- file.path(ribopilot_ref_dir(), "blastdb")
  path.expand(d)
}
```

Append to `R/reference_dir.R`:

```r
#' Is `prefix` usable as a local BLAST panel?
#'
#' Read-only on purpose: a panel is often on a shared read-only mount, so this
#' must not reuse validate_ref_dir(), which requires write access and creates
#' the directory. An empty prefix is valid and selects the remote search.
#' @return list(ok, msg); on success msg is the normalised prefix.
#' @noRd
validate_blastdb_prefix <- function(prefix) {
  prefix <- trimws(prefix %||% "")
  if (!nzchar(prefix)) return(list(ok = TRUE, msg = ""))
  prefix <- path.expand(prefix)
  if (grepl("\\s", prefix)) {
    return(list(ok = FALSE, msg = "Panel path cannot contain spaces; blastn cannot read it."))
  }
  if (!blast_db_exists(prefix)) {
    return(list(ok = FALSE, msg = paste0(
      "No BLAST database at ", prefix,
      ". Give the prefix without a file extension, for example ",
      "/data/blastdb/2026.08.1/ribopilot_rrna")))
  }
  list(ok = TRUE, msg = prefix)
}
```

- [ ] **Step 4: Run to verify pass**

Run both test files. Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/db_dist.R R/reference_dir.R tests/testthat/test-db-dist.R tests/testthat/test-reference-dir.R
git commit -m "feat: resolve panel location and validate a panel prefix read-only"
```

---

### Task 3: The download helper

**Files:**
- Modify: `R/db_dist.R`
- Test: `tests/testthat/test-db-dist.R`

**Interfaces:**
- Consumes: `rp_sha256()`, `read_panel_manifest()`, `ribopilot_db_dir()` from Tasks 1-2.
- Produces: `install_ribopilot_db(version, dest, decoys, base_url, force)` returning the panel prefix invisibly.

Release assets, at tag `db-<version>`:
- `ribopilot_rrna-<version>.tar.gz` and `ribopilot_rrna-<version>.tar.gz.sha256`
- `decoys-<version>.tar.gz` and `decoys-<version>.tar.gz.sha256`

Both tarballs extract members with no leading directory, straight into `<dest>/<version>/`.

The `.sha256` file holds one line in `sha256sum` format: `<64 hex><two spaces><filename>`.

- [ ] **Step 1: Write the failing tests**

These use a `file://` base URL, so they need no network. Build a real tarball in `tempdir()`.

```r
make_fake_release <- function(root, version = "9.9.9") {
  stage <- file.path(root, "stage"); dir.create(stage, recursive = TRUE)
  writeLines(sprintf('{"manifest_schema":1,"db_name":"ribopilot_rrna","db_version":"%s"}', version),
    file.path(stage, "manifest.json"))
  writeLines("x", file.path(stage, "ribopilot_rrna.nin"))
  writeLines("x", file.path(stage, "ribopilot_rrna.nsq"))
  rel <- file.path(root, "rel"); dir.create(rel, recursive = TRUE)
  tgz <- file.path(rel, sprintf("ribopilot_rrna-%s.tar.gz", version))
  withr::with_dir(stage, utils::tar(tgz, files = list.files("."), compression = "gzip"))
  writeLines(paste0(rp_sha256(tgz), "  ", basename(tgz)), paste0(tgz, ".sha256"))
  paste0("file://", rel)
}

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
  sums <- list.files(file.path(root, "rel"), pattern = "\\.sha256$", full.names = TRUE)
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
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, "could not find function install_ribopilot_db".

- [ ] **Step 3: Implement**

Append to `R/db_dist.R`:

```r
#' Base URL of the panel release assets
#' @noRd
RP_DB_RELEASE_BASE <- "https://github.com/dmacguigan/RiboPilot/releases/download"

#' Default panel version
#' @noRd
RP_DB_VERSION <- "2026.08.1"

#' Download one asset and check it against its published .sha256
#' @noRd
fetch_verified <- function(base_url, tag, asset, dir) {
  url <- paste(base_url, tag, asset, sep = "/")
  dest <- file.path(dir, asset)
  if (utils::download.file(url, dest, mode = "wb", quiet = TRUE) != 0 || !file.exists(dest)) {
    stop("Could not download ", url, call. = FALSE)
  }
  sum_dest <- paste0(dest, ".sha256")
  if (utils::download.file(paste0(url, ".sha256"), sum_dest, mode = "wb", quiet = TRUE) != 0) {
    stop("Could not download the checksum for ", asset, call. = FALSE)
  }
  line <- readLines(sum_dest, warn = FALSE)[1]
  want <- regmatches(line, regexpr("[0-9a-fA-F]{64}", line))
  got <- rp_sha256(dest)
  if (length(want) != 1 || is.na(got)) {
    stop("Could not verify the checksum of ", asset, call. = FALSE)
  }
  if (!identical(tolower(want), got)) {
    stop("Checksum mismatch for ", asset, ": the download is corrupt or truncated.",
      call. = FALSE)
  }
  dest
}

#' Download and install a RiboPilot BLAST panel
#'
#' Extraction is staged in a sibling directory and moved into place only after
#' every asset verifies, so an interrupted or corrupt download never leaves a
#' half-panel that blast_db_exists() would call usable.
#'
#' @param version Panel version, matching release tag `db-<version>`.
#' @param dest Parent directory; the panel lands in `<dest>/<version>`.
#' @param decoys Also fetch the NCBI decoy databases (adds roughly 110 MB).
#' @param base_url Release download base; overridable for testing.
#' @param force Re-download even when the panel is already present.
#' @return (invisibly) the panel prefix.
#' @export
install_ribopilot_db <- function(version = RP_DB_VERSION,
                                 dest = ribopilot_db_dir(),
                                 decoys = TRUE,
                                 base_url = RP_DB_RELEASE_BASE,
                                 force = FALSE) {
  dest <- path.expand(dest)
  final <- file.path(dest, version)
  prefix <- file.path(final, "ribopilot_rrna")
  if (!force && blast_db_exists(prefix)) {
    message("Panel ", version, " is already installed at ", final)
    return(invisible(prefix))
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (unname(file.access(dest, 2)) != 0) {
    stop(dest, " is not writable. Set RIBOPILOT_DB_DIR to a writable location.",
      call. = FALSE)
  }
  stage <- file.path(dest, paste0(".stage-", version))
  unlink(stage, recursive = TRUE)
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  tag <- paste0("db-", version)
  assets <- sprintf("ribopilot_rrna-%s.tar.gz", version)
  if (isTRUE(decoys)) assets <- c(assets, sprintf("decoys-%s.tar.gz", version))
  for (a in assets) {
    message("Downloading ", a, " ...")
    tgz <- fetch_verified(base_url, tag, a, stage)
    utils::untar(tgz, exdir = stage)
    unlink(c(tgz, paste0(tgz, ".sha256")))
  }
  if (!blast_db_exists(file.path(stage, "ribopilot_rrna"))) {
    stop("The downloaded archive did not contain a usable panel.", call. = FALSE)
  }
  unlink(final, recursive = TRUE)
  if (!file.rename(stage, final)) {
    stop("Could not move the panel into ", final, call. = FALSE)
  }
  message("Panel ", version, " installed at ", final)
  invisible(prefix)
}
```

- [ ] **Step 4: Run to verify pass**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-db-dist.R")'`
Expected: all PASS.

- [ ] **Step 5: Document the export**

Run: `Rscript -e 'roxygen2::roxygenise(".")'` and confirm `NAMESPACE` gains `export(install_ribopilot_db)`.

- [ ] **Step 6: Commit**

```bash
git add R/db_dist.R tests/testthat/test-db-dist.R NAMESPACE man/
git commit -m "feat: add install_ribopilot_db to fetch and verify a panel"
```

---

### Task 4: Persist and configure the panel path

**Files:**
- Modify: `R/reference_dir.R:32-46` (`set_ref_config`), `R/reference_dir.R:15-25` (`ref_config`)
- Modify: `R/app_reference.R:138-161` (modal), `R/app_reference.R:769-776` (observer)
- Test: `tests/testthat/test-reference-dir.R`, `tests/testthat/test-reference-app.R`

**Interfaces:**
- Consumes: `validate_blastdb_prefix()` from Task 2.
- Produces: `set_ref_config(con, ..., blastdb=)`; `ref_config(con)$blastdb`.

The `blastdb` column already exists in `REFERENCE_CONFIG_DDL` (`R/db_schema.R:151`) and has no writer. `ref_config()` currently does not return it. Reuse the existing test helper `mk_cfg_db()` at `tests/testthat/test-reference-dir.R:3`.

- [ ] **Step 1: Write the failing tests**

```r
# in test-reference-dir.R
test_that("ref_config reports the stored panel prefix, empty when unset", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(ref_config(con)$blastdb, "")
  set_ref_config(con, blastdb = "/data/blastdb/2026.08.1/ribopilot_rrna")
  expect_equal(ref_config(con)$blastdb, "/data/blastdb/2026.08.1/ribopilot_rrna")
})

test_that("set_ref_config clears the panel prefix when given an empty string", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  set_ref_config(con, blastdb = "/data/p/ribopilot_rrna")
  set_ref_config(con, blastdb = "")
  expect_equal(ref_config(con)$blastdb, "")
  expect_true(is.na(DBI::dbGetQuery(con,
    "SELECT blastdb FROM reference_config WHERE rowid_key = 1")$blastdb[1]))
})

test_that("set_ref_config leaves the panel prefix alone when not supplied", {
  con <- mk_cfg_db(); on.exit(DBI::dbDisconnect(con))
  set_ref_config(con, blastdb = "/data/p/ribopilot_rrna")
  set_ref_config(con, min_records = 500)
  expect_equal(ref_config(con)$blastdb, "/data/p/ribopilot_rrna")
})

# in test-reference-app.R
test_that("the settings modal exposes a panel path input carrying the current value", {
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000,
    blastdb = "/data/blastdb/2026.08.1/ribopilot_rrna")
  html <- as.character(RiboPilot:::reference_settings_modal(cfg))
  expect_match(html, "ref_blastdb_path")
  expect_match(html, "/data/blastdb/2026.08.1/ribopilot_rrna", fixed = TRUE)
})

test_that("the settings modal renders when no panel is configured", {
  cfg <- list(ref_dir = "/refs", min_records = 1000, max_records = 25000, blastdb = "")
  expect_match(as.character(RiboPilot:::reference_settings_modal(cfg)), "ref_blastdb_path")
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-reference-dir.R")'`
Expected: FAIL, `ref_config(con)$blastdb` is NULL.

- [ ] **Step 3: Implement the config layer**

In `R/reference_dir.R`, add to the list `ref_config()` returns:

```r
    max_records = if (nrow(d) > 0) as.integer(d$max_records[1] %|NA|% 25000L) else 25000L,
    blastdb = {
      b <- if (nrow(d) > 0 && "blastdb" %in% names(d)) d$blastdb[1] else NA_character_
      if (!is.na(b) && nzchar(trimws(b))) path.expand(trimws(b)) else ""
    }
```

Add the parameter to `set_ref_config()`, mirroring the `ref_dir` branch:

```r
set_ref_config <- function(con, ref_dir = NULL, min_records = NULL, max_records = NULL,
                           blastdb = NULL) {
```

and inside, after the `max_records` branch:

```r
  if (!is.null(blastdb)) {
    DBI::dbExecute(con, "UPDATE reference_config SET blastdb = ? WHERE rowid_key = 1",
      params = list(if (nzchar(trimws(blastdb))) path.expand(trimws(blastdb)) else NA_character_))
  }
```

- [ ] **Step 4: Implement the UI**

In `R/app_reference.R`, inside `reference_settings_modal()`, after the reference-folder `textInput` and before the records heading:

```r
    shiny::h5("Local BLAST panel", style = SETTINGS_HEAD),
    shiny::div(
      style = "font-size:0.85em; color:#666; margin-bottom:0.4em;",
      "Path prefix of a downloaded panel, without a file extension. ",
      "Leave empty to search NCBI over the network."),
    shiny::textInput("ref_blastdb_path", NULL, width = "100%",
      value = cfg$blastdb %||% "",
      placeholder = "/data/blastdb/2026.08.1/ribopilot_rrna"),
```

In the `ref_dir_apply` observer, validate the panel before saving anything, so one bad field saves nothing:

```r
shiny::observeEvent(input$ref_dir_apply, {
  v <- validate_ref_dir(input$ref_dir_path %||% "")
  if (!isTRUE(v$ok)) { notify(v$msg, "error"); return() }
  b <- validate_blastdb_prefix(input$ref_blastdb_path %||% "")
  if (!isTRUE(b$ok)) { notify(b$msg, "error"); return() }
  with_con(function(con) set_ref_config(con, ref_dir = v$msg,
    min_records = input$ref_min, max_records = input$ref_max, blastdb = b$msg))
  shiny::removeModal(); bump()
  notify(paste("Reference folder set to", v$msg))
})
```

- [ ] **Step 5: Run to verify pass**

Run both test files. Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add R/reference_dir.R R/app_reference.R tests/testthat/test-reference-dir.R tests/testthat/test-reference-app.R
git commit -m "feat: configure the local BLAST panel from reference settings"
```

---

### Task 5: Record which panel produced a project's results

**Files:**
- Modify: `R/db_schema.R:146-153` (`REFERENCE_CONFIG_DDL`)
- Modify: `R/db_dist.R`
- Modify: `R/nf_command.R:39` (record at launch)
- Test: `tests/testthat/test-db-dist.R`

**Interfaces:**
- Consumes: `read_panel_manifest()` from Task 1; `nf_blastdb(path)` from `R/nf_command.R:85`.
- Produces: `record_blastdb_provenance(con, prefix)` returning the version invisibly; `blastdb_provenance(con)` returning a string for display.

Per Global Constraints there is no migration: the column is added to the DDL and every read is guarded, so an older project reports `""`.

- [ ] **Step 1: Write the failing tests**

```r
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

test_that("provenance reads degrade to empty on a project predating the column", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE reference_config (rowid_key INTEGER PRIMARY KEY,
    ref_dir TEXT, min_records INTEGER, max_records INTEGER, blastdb TEXT)")
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  expect_equal(blastdb_provenance(con), "")
  expect_silent(record_blastdb_provenance(con, ""))
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, "could not find function record_blastdb_provenance".

- [ ] **Step 3: Add the column**

In `R/db_schema.R`, in `REFERENCE_CONFIG_DDL`, after `blastdb TEXT,`:

```sql
  blastdb_version TEXT,
```

- [ ] **Step 4: Implement**

Append to `R/db_dist.R`:

```r
#' Record which panel a run searched against
#'
#' Written when a run launches so a result can be traced to the database that
#' produced it. Silent and non-fatal on a project predating the column.
#' @return (invisibly) the recorded version, or "" when there is nothing to record.
#' @noRd
record_blastdb_provenance <- function(con, prefix) {
  prefix <- trimws(prefix %||% "")
  if (!nzchar(prefix)) return(invisible(""))
  m <- read_panel_manifest(dirname(path.expand(prefix)))
  if (is.null(m)) return(invisible(""))
  v <- paste(m$db_name, m$db_version)
  ok <- tryCatch({
    DBI::dbExecute(con, "UPDATE reference_config SET blastdb_version = ? WHERE rowid_key = 1",
      params = list(v))
    TRUE
  }, error = function(e) FALSE)
  invisible(if (ok) as.character(m$db_version) else "")
}

#' The recorded panel, for display; "" when none
#' @noRd
blastdb_provenance <- function(con) {
  v <- tryCatch(DBI::dbGetQuery(con,
    "SELECT blastdb_version FROM reference_config WHERE rowid_key = 1")$blastdb_version,
    error = function(e) NULL)
  if (length(v) < 1 || is.na(v[1]) || !nzchar(trimws(v[1]))) "" else trimws(v[1])
}
```

- [ ] **Step 5: Call it at launch**

In `R/nf_command.R`, immediately after `blastdb <- nf_blastdb(path)` (line 39):

```r
  if (nzchar(blastdb)) {
    # local() so the connection closes here rather than on nf_command's exit,
    # and so `con` cannot shadow anything in the enclosing frame
    local({
      pcon <- tryCatch(open_project_db(file.path(path, ".sqlite")), error = function(e) NULL)
      if (is.null(pcon)) return(invisible(NULL))
      on.exit(DBI::dbDisconnect(pcon), add = TRUE)
      tryCatch(record_blastdb_provenance(pcon, blastdb), error = function(e) NULL)
    })
  }
```

Note the surrounding function already opens its own connection for `nf_blastdb()`; this is a second, short-lived one. Do not merge them, and do not name this one `con`.

- [ ] **Step 6: Run to verify pass**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-db-dist.R")'`
Then the nextflow tests: `testthat::test_file("tests/testthat/test-nf-command-quoting.R")`.
Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add R/db_schema.R R/db_dist.R R/nf_command.R tests/testthat/test-db-dist.R
git commit -m "feat: record the panel version a project searched against"
```

---

### Task 6: Build-side decoy fetch and release packaging

**Files:**
- Create: `dev/refdb/fetch_decoys.sh`
- Create: `dev/refdb/package_release.sh`
- Test: `tests/testthat/test-db-dist.R` (a shell-syntax check only)

**Interfaces:**
- Consumes: the built panel directory at `<build>/blastdb/<version>/`.
- Produces: release assets consumed by `install_ribopilot_db()` (Task 3), named exactly `ribopilot_rrna-<version>.tar.gz`, `decoys-<version>.tar.gz`, and a `.sha256` beside each.

These scripts are run by hand by the maintainer. **Neither is executed during this task**, and `package_release.sh` must never push without `--publish` being passed explicitly.

- [ ] **Step 1: Write `dev/refdb/fetch_decoys.sh`**

```bash
#!/usr/bin/env bash
# Fetch the NCBI decoy databases into a panel directory.
#
# The decoys sit beside the panel because rp_dbarg.R defaults the decoy
# directory to dirname(panel_prefix). taxdb comes along so staxids resolve.
set -ueo pipefail

DEST="${1:?usage: fetch_decoys.sh <panel_dir>}"
DBS="16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb"

command -v update_blastdb.pl >/dev/null || { echo "update_blastdb.pl not on PATH" >&2; exit 1; }
mkdir -p "$DEST"
cd "$DEST"

for d in $DBS; do
  echo "fetching $d ..."
  update_blastdb.pl --decompress "$d"
done

# update_blastdb.pl verifies each archive's md5 as it goes; record what landed
: > DECOY_VERSION.txt
for d in $DBS; do
  printf '%s\t%s\n' "$d" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> DECOY_VERSION.txt
done
echo "decoys installed in $DEST"
```

- [ ] **Step 2: Write `dev/refdb/package_release.sh`**

```bash
#!/usr/bin/env bash
# Package a built panel (and its decoys) as GitHub release assets.
#
# Ships the BLAST volume set and the manifest only. The .fasta, .taxid.map and
# .dropped.tsv in the build directory are build artifacts, not query-time files.
set -ueo pipefail

SRC="${1:?usage: package_release.sh <panel_dir> <version> [--publish]}"
VER="${2:?usage: package_release.sh <panel_dir> <version> [--publish]}"
PUBLISH="${3:-}"
OUT="${SRC%/}/../release-${VER}"
REPO="dmacguigan/RiboPilot"

mkdir -p "$OUT"
cd "$SRC"

PANEL_FILES="manifest.json"
for e in ndb nhr nin njs nog nos not nsq ntf nto; do
  [ -f "ribopilot_rrna.$e" ] && PANEL_FILES="$PANEL_FILES ribopilot_rrna.$e"
done

echo "packaging panel ..."
tar -czf "$OUT/ribopilot_rrna-${VER}.tar.gz" $PANEL_FILES

DECOYS=""
for n in 16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb; do
  for f in "$n".*; do [ -e "$f" ] && DECOYS="$DECOYS $f"; done
done
if [ -n "$DECOYS" ]; then
  echo "packaging decoys ..."
  tar -czf "$OUT/decoys-${VER}.tar.gz" $DECOYS
fi

cd "$OUT"
for f in *.tar.gz; do
  sha256sum "$f" > "$f.sha256"
done
ls -lh "$OUT"

if [ "$PUBLISH" != "--publish" ]; then
  echo "assets staged in $OUT; re-run with --publish to create the release"
  exit 0
fi

command -v gh >/dev/null || { echo "gh not on PATH" >&2; exit 1; }
# --clobber is deliberately absent: a published database version is immutable
gh release create "db-${VER}" ./*.tar.gz ./*.sha256 \
  --repo "$REPO" --title "rDNA panel ${VER}" --latest=false \
  --notes "Metazoan nuclear rDNA BLAST panel ${VER}. Install with RiboPilot::install_ribopilot_db(\"${VER}\")."
```

- [ ] **Step 3: Make them executable and syntax-check both**

```bash
chmod +x dev/refdb/fetch_decoys.sh dev/refdb/package_release.sh
bash -n dev/refdb/fetch_decoys.sh && bash -n dev/refdb/package_release.sh && echo "syntax OK"
```

Expected: `syntax OK`.

- [ ] **Step 4: Add a guard test that the asset names agree with the downloader**

The downloader and the packager must agree on names; a mismatch is only discoverable in production otherwise.

```r
test_that("the packaging script emits the asset names the installer requests", {
  skip_on_cran()
  s <- readLines("../../dev/refdb/package_release.sh", warn = FALSE)
  txt <- paste(s, collapse = "\n")
  expect_match(txt, 'ribopilot_rrna-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'decoys-${VER}.tar.gz', fixed = TRUE)
  expect_match(txt, 'db-${VER}', fixed = TRUE)
  # the installer builds exactly these
  expect_equal(sprintf("ribopilot_rrna-%s.tar.gz", "1.2.3"), "ribopilot_rrna-1.2.3.tar.gz")
  expect_equal(paste0("db-", "1.2.3"), "db-1.2.3")
})
```

- [ ] **Step 5: Run the full suite both ways**

```bash
Rscript -e 'pkgload::load_all("."); r <- testthat::test_dir("tests/testthat", reporter="silent"); d <- as.data.frame(r); cat("failed:", sum(d$failed), "error:", sum(d$error), "skipped:", sum(d$skipped), "passed:", sum(d$passed), "\n")'
NOT_CRAN=true Rscript -e '...same...'
```

Expected: 0 failed, 0 errors in both.

- [ ] **Step 6: Commit**

```bash
git add dev/refdb/fetch_decoys.sh dev/refdb/package_release.sh tests/testthat/test-db-dist.R
git commit -m "feat: add decoy fetch and release packaging scripts"
```

---

## Deferred, with reasons

- **`:ro` on the panel bind mount.** The panel directory is mounted read-write today. BLAST never writes to a database, so the exposure is an accidental write by unrelated code, not a live bug. Making all mounts read-only is wrong (the reference directory is written during clade builds), so this needs a per-mount flag in `container_bind_opts()`, which is a wider change than this stage warrants.
- **`.dropped.tsv` as a release asset.** 15.6 MB of build diagnostics answering "why is my taxon missing". Useful, but not needed to run a search; publish separately if users ask.
- **Docker image layer for the panel.** Deliberately not baked in: the image is already 5.98 GB and versions on every package bump, while the panel versions quarterly. Runtime mounting is already wired through `reference_mount_dirs()`.

## Self-Review

**Spec coverage.** The Stage 4a scope-out note lists six items for 4b: decoy fetch at build time (Task 6), Docker layer strategy (deferred above, with the reason the spec itself gives), the GitHub release (Task 6), the manifest (Tasks 1 and 5), the download helper (Task 3), and recording the database version per project (Task 5). The panel-path UI setter, flagged as missing after Stage 4a, is Task 4.

**Placeholders.** None. Every code step carries the code to write.

**Type consistency.** `validate_blastdb_prefix()` returns `list(ok, msg)`, matching `validate_ref_dir()`, and Task 4 consumes `b$msg` as the value to store. `install_ribopilot_db()` returns the prefix; `record_blastdb_provenance()` takes that same prefix and calls `dirname()` on it to find the manifest, which is where Task 3 puts it. Asset names appear in Tasks 3 and 6 and are pinned together by the Task 6 guard test.

**Known coupling.** Task 5 adds a column that Task 4's `ref_config()` does not read; they touch the same table but different columns, so order between them does not matter. Task 3's tests depend on Tasks 1 and 2, so the order 1 -> 2 -> 3 is required. Tasks 4, 5, 6 are independent of each other.
