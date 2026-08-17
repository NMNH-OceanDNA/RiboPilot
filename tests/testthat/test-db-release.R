# The release packaging is the one piece of the refdb rebuild CI can afford: a
# mock panel directory of correctly-named files, packaged for real, in about a
# second. Text-grepping the scripts proved only that a string was present.

pkg_sh <- "../../dev/refdb/package_release.sh"
decoys_sh <- "../../dev/refdb/fetch_decoys.sh"
refdb_sh <- "../../dev/refdb/build_refdb.sh"

# dev/ is in .Rbuildignore, so the guard belongs inside each test rather than at
# file level, where one skip would drop the docker tests too.
skip_without_packager <- function() {
  skip_on_cran()
  skip_if_not(file.exists(pkg_sh),
    "dev/ is Rbuildignored; this runs from the source tree only")
  for (tool in c("bash", "jq", "tar", "sha256sum")) {
    skip_if(Sys.which(tool) == "", paste(tool, "is not on PATH"))
  }
}

RP_TAXONOMY_FILES <- c("taxdb.btd", "taxdb.bti", "taxonomy4blast.sqlite3")
RP_PANEL_EXT <- c("ndb", "nhr", "nin", "njs", "nog", "nos", "not", "nsq", "ntf", "nto")
# Exactly the volume sets in the read-only fixture at ~/.cache/ribopilot-blastdb:
# 16S carries .nnd/.nni, the two eukaryote sets do not. A mock with four
# extensions apiece would let a packaging bug through that a real panel hits.
RP_DECOY_EXT <- list(
  "16S_ribosomal_RNA" = c("ndb", "nhr", "nin", "nnd", "nni", "nog", "nos", "not",
    "nsq", "ntf", "nto"),
  "SSU_eukaryote_rRNA" = c("ndb", "nhr", "nin", "nog", "nos", "not", "nsq", "ntf", "nto"),
  "LSU_eukaryote_rRNA" = c("ndb", "nhr", "nin", "nog", "nos", "not", "nsq", "ntf", "nto"))
RP_DECOY_VOLUMES <- unlist(mapply(function(n, e) paste0(n, ".", e),
  names(RP_DECOY_EXT), RP_DECOY_EXT, USE.NAMES = FALSE))
# update_blastdb.pl preserves these deliberately, as the cache key for its next
# run. They are in the panel directory after every real fetch and must not ship.
RP_DECOY_MARKERS <- paste0(c(names(RP_DECOY_EXT), "taxdb"), ".tar.gz.md5")

# A panel directory holding what the build and decoys stages really leave: the
# panel volume set, a manifest, the three decoy databases at their real volume
# counts, the cache markers, the taxonomy files and the version record. Files
# are named, not built; nothing here reads them.
mock_panel_dir <- function(ver = "9.9.9", drop = character(0), extra = character(0),
                           manifest_ver = ver, .envir = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .envir)
  d <- file.path(root, "blastdb", ver)
  dir.create(d, recursive = TRUE)
  files <- c(paste0("ribopilot_rrna.", RP_PANEL_EXT), RP_DECOY_VOLUMES,
    RP_DECOY_MARKERS, RP_TAXONOMY_FILES, "DECOY_VERSION.txt")
  for (f in c(setdiff(files, drop), extra)) writeLines(f, file.path(d, f))
  # the drop ledger, which the build stage leaves in the panel directory: it is
  # not in either tarball but it does ship, as its own release asset
  if (!"ribopilot_rrna.dropped.tsv" %in% drop) {
    writeLines(c("AAAA00000000\twgs_master_unavailable\t",
      "MZ000001\ttaxid_unresolved\t"), file.path(d, "ribopilot_rrna.dropped.tsv"))
  }
  writeLines(sprintf('{"db_name": "ribopilot_rrna", "db_version": "%s"}', manifest_ver),
    file.path(d, "manifest.json"))
  d
}

run_packager <- function(src, ver = basename(src)) {
  out <- suppressWarnings(system2("bash",
    c(shQuote(normalizePath(pkg_sh)), shQuote(src), shQuote(ver)),
    stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out,
    rel = file.path(dirname(src), paste0("release-", ver)))
}

tar_contents <- function(path) {
  sort(suppressWarnings(system2("tar", c("-tzf", shQuote(path)), stdout = TRUE)))
}

rp_db_version_string <- function() {
  sub('.*"(.*)".*', "\\1",
    grep("^RP_DB_VERSION", readLines("../../R/db_dist.R", warn = FALSE), value = TRUE)[1])
}

test_that("a packaged release installs whole and is discoverable", {
  skip_without_packager()
  src <- mock_panel_dir()
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))

  # install_ribopilot_db() reads <base_url>/<tag>/<asset>, so the staged assets
  # are served from a file:// URL exactly as a release would serve them.
  serve <- withr::local_tempdir()
  dir.create(file.path(serve, "db-9.9.9"))
  file.copy(list.files(r$rel, full.names = TRUE), file.path(serve, "db-9.9.9"))
  dest <- withr::local_tempdir()
  suppressMessages(install_ribopilot_db("9.9.9", dest = dest,
    base_url = paste0("file://", serve)))

  final <- file.path(dest, "9.9.9")
  expect_equal(RiboPilot:::panel_missing_parts(final), character(0))
  # the whole point of packaging the decoys: without these two the installed
  # panel searches without taxonomy and nobody can say which decoys it holds
  expect_true(file.exists(file.path(final, "taxonomy4blast.sqlite3")))
  expect_true(file.exists(file.path(final, "DECOY_VERSION.txt")))
  expect_equal(RiboPilot:::discover_panel(dest), file.path(final, "ribopilot_rrna"))
})

test_that("the decoy tarball carries every taxonomy file, including the one taxdb.* misses", {
  skip_without_packager()
  src <- mock_panel_dir()
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  got <- tar_contents(file.path(r$rel, "decoys-9.9.9.tar.gz"))
  expect_setequal(got, c(RP_DECOY_VOLUMES, RP_TAXONOMY_FILES, "DECOY_VERSION.txt"))
  # taxonomy4blast.sqlite3 is not a taxdb.* file, and shipping without it makes
  # -taxids return unfiltered rows while blastn still exits 0
  expect_true("taxonomy4blast.sqlite3" %in% got)
  # the cache markers were in the directory and must not be in the tarball:
  # shipped, they become the key update_blastdb.pl reads in a user's decoy dir
  expect_true(all(file.exists(file.path(src, RP_DECOY_MARKERS))))
  expect_false(any(RP_DECOY_MARKERS %in% got))
  expect_false(any(grepl("[.]md5$", got)))
  # flat, never nested under a directory: the installer untars into the panel dir
  expect_false(any(grepl("/", got, fixed = TRUE)))
})

test_that("packaging refuses when any taxonomy file or the version record is absent", {
  skip_without_packager()
  for (f in c(RP_TAXONOMY_FILES, "DECOY_VERSION.txt")) {
    src <- mock_panel_dir(drop = f)
    r <- run_packager(src)
    expect_gt(r$status, 0L)
    expect_true(any(grepl(paste("no", f), r$out, fixed = TRUE)), label = f)
    # nothing at all: a lone panel tarball with no checksum and no decoys, or
    # even an empty release-<VER>, is a version-carrying artifact of a run that
    # refused, and --check-version reports the directory as a mismatch
    expect_false(dir.exists(r$rel), label = f)
  }
})

test_that("packaging refuses a version that disagrees with the manifest or the directory", {
  skip_without_packager()
  src <- mock_panel_dir(manifest_ver = "8.8.8")
  r <- run_packager(src)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("db_version is '8.8.8'", r$out, fixed = TRUE)))

  src <- mock_panel_dir()
  r <- run_packager(src, ver = "7.7.7")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("is named '9.9.9', not '7.7.7'", r$out, fixed = TRUE)))
  expect_false(dir.exists(r$rel))
})

test_that("packaging refuses while a decoy archive is still sitting in the panel", {
  skip_without_packager()
  src <- mock_panel_dir(extra = "16S_ribosomal_RNA.tar.gz")
  r <- run_packager(src)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("the decoy fetch did not finish", r$out, fixed = TRUE)))
  expect_false(dir.exists(r$rel))
})

test_that("the panel tarball ships the volume set and the manifest and nothing else", {
  skip_without_packager()
  src <- mock_panel_dir(extra = c("ribopilot_rrna.fasta", "ribopilot_rrna.taxid.map"))
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  expect_setequal(tar_contents(file.path(r$rel, "ribopilot_rrna-9.9.9.tar.gz")),
    c("manifest.json", paste0("ribopilot_rrna.", RP_PANEL_EXT)))
})

# Spec section 6: the drop ledger is "the only record of the discarded records",
# and `clean` deletes it from the panel directory once the release assets verify.
# If packaging does not ship it, it never leaves the build machine.
test_that("the drop ledger ships as its own gzipped asset, with a checksum, outside both tarballs", {
  skip_without_packager()
  src <- mock_panel_dir()
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))

  asset <- file.path(r$rel, "dropped-9.9.9.tsv.gz")
  expect_true(file.exists(asset))
  expect_equal(readLines(gzfile(asset), warn = FALSE),
    readLines(file.path(src, "ribopilot_rrna.dropped.tsv"), warn = FALSE))

  # same sidecar shape as the tarballs, and it actually verifies
  line <- readLines(paste0(asset, ".sha256"), warn = FALSE)
  expect_length(line, 1L)
  expect_match(line, "^[0-9a-f]{64}  dropped-9[.]9[.]9[.]tsv[.]gz$")
  expect_equal(system2("bash", c("-c", shQuote(sprintf(
    "cd %s && sha256sum -c --status dropped-9.9.9.tsv.gz.sha256", shQuote(r$rel)))),
    stdout = FALSE, stderr = FALSE), 0L)

  # an asset in its own right: nobody should unpack an 82 MB index to read it
  expect_false("ribopilot_rrna.dropped.tsv" %in%
    tar_contents(file.path(r$rel, "ribopilot_rrna-9.9.9.tar.gz")))
  expect_false("ribopilot_rrna.dropped.tsv" %in%
    tar_contents(file.path(r$rel, "decoys-9.9.9.tar.gz")))
})

test_that("packaging refuses a panel with no drop ledger, and the publish list names the asset", {
  skip_without_packager()
  src <- mock_panel_dir(drop = "ribopilot_rrna.dropped.tsv")
  r <- run_packager(src)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("no ribopilot_rrna.dropped.tsv", r$out, fixed = TRUE)))
  expect_false(dir.exists(r$rel))

  # the asset has to be in the gh upload list too, or it is staged and never
  # published, which is the same as not shipping it
  s <- readLines(pkg_sh, warn = FALSE)
  gh <- grep("^gh release create", s, value = TRUE)
  expect_length(gh, 1L)
  expect_match(gh, "dropped-${VER}.tsv.gz", fixed = TRUE)
})

test_that("packaging refuses a panel missing any of the required volume files", {
  skip_without_packager()
  for (e in RiboPilot:::PANEL_REQUIRED_EXT) {
    src <- mock_panel_dir(drop = paste0("ribopilot_rrna.", e))
    r <- run_packager(src)
    expect_gt(r$status, 0L)
    expect_true(any(grepl("the panel is incomplete", r$out, fixed = TRUE)), label = e)
  }
})

test_that("the checksum sidecars are in the shape sha256sum -c and the installer both parse", {
  skip_without_packager()
  src <- mock_panel_dir()
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  for (a in c("ribopilot_rrna-9.9.9.tar.gz", "decoys-9.9.9.tar.gz")) {
    line <- readLines(file.path(r$rel, paste0(a, ".sha256")), warn = FALSE)
    expect_length(line, 1L)
    expect_match(line, sprintf("^[0-9a-f]{64}  %s$", a), label = a)
    expect_equal(RiboPilot:::parse_sha256_line(line),
      RiboPilot:::rp_sha256(file.path(r$rel, a)), label = a)
  }
})

test_that("the shell decoy lists are RP_DECOY_DBS plus exactly taxdb", {
  skip_on_cran()
  skip_if_not(file.exists(refdb_sh),
    "dev/ is Rbuildignored; this runs from the source tree only")
  # taxdb is a taxonomy sidecar, not a searchable database, so it must stay out
  # of RP_DECOY_DBS: rp_dbarg composes that vector into the blastn -db list.
  shell_list <- function(path, name) {
    line <- grep(paste0("^", name, "="), readLines(path, warn = FALSE), value = TRUE)[1]
    strsplit(trimws(gsub('"', "", sub(paste0("^", name, "="), "", line))), " +")[[1]]
  }
  expect_setequal(shell_list(refdb_sh, "DECOY_DBS"),
    c(RiboPilot:::RP_DECOY_DBS, "taxdb"))
  expect_setequal(shell_list(decoys_sh, "SEARCHABLE_DBS"), RiboPilot:::RP_DECOY_DBS)
  expect_setequal(shell_list(decoys_sh, "TAXONOMY_FILES"), RP_TAXONOMY_FILES)
  expect_setequal(shell_list(decoys_sh, "REQUIRED_EXT"), RiboPilot:::PANEL_REQUIRED_EXT)
  expect_setequal(names(RP_DECOY_EXT), RiboPilot:::RP_DECOY_DBS)
})

# --- fetch_decoys.sh ---------------------------------------------------------

# update_blastdb.pl and the three NCBI endpoints are stubbed; makeblastdb and
# blastdbcmd are the real ones, because the verification this stage exists for
# is "does blastdbcmd accept what landed", which a stubbed blastdbcmd cannot
# answer. RP_STUB_BREAK damages exactly one thing so one guard is under test.
RP_DECOY_SEQS <- c(`16S_ribosomal_RNA` = 27644L, SSU_eukaryote_rRNA = 8784L,
  LSU_eukaryote_rRNA = 6575L)

decoy_stub_bin <- function(meta_null = character(0), .envir = parent.frame()) {
  bindir <- withr::local_tempdir(.local_envir = .envir)
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -e",
    'db="${!#}"',
    'brk="${RP_STUB_BREAK:-}"',
    'here="$(cd "$(dirname "$0")" && pwd)"',
    "# preserved by the real tool as the cache key for its next run",
    'echo "0000  $db.tar.gz" > "$db.tar.gz.md5"',
    '[ "$brk" != uptodate ] || { echo "$db is up to date"; exit 0; }',
    'if [ "$db" = taxdb ]; then',
    '  for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do echo stub > "$f"; done',
    '  [ "$brk" != no_taxonomy ] || rm -f taxonomy4blast.sqlite3',
    "  exit 0",
    "fi",
    'n="$(awk -v d="$db" \'$1 == d {print $2}\' "$here/counts.txt")"',
    "# the read-only fixture's own staleness: 27,649 installed against 27,644 live",
    '[ "$brk" != stale ] || [ "$db" != 16S_ribosomal_RNA ] || n=$((n + 5))',
    'awk -v n="$n" \'BEGIN {for (i = 1; i <= n; i++) printf(">s%d stub\\nACGTACGTACGTACGTACGT\\n", i)}\' > "$db.stub.fasta"',
    'makeblastdb -in "$db.stub.fasta" -dbtype nucl -blastdb_version 5 -out "$db" > /dev/null',
    'rm -f "$db.stub.fasta"',
    'case "$brk" in',
    "  missing_volume) rm -f 16S_ribosomal_RNA.nsq ;;",
    "  empty_volume) : > 16S_ribosomal_RNA.nsq ;;",
    '  bad_index) printf "not a blast index" > 16S_ribosomal_RNA.nin ;;',
    "esac"
  ), file.path(bindir, "update_blastdb.pl"))
  writeLines(sprintf("%s %d", names(RP_DECOY_SEQS), RP_DECOY_SEQS),
    file.path(bindir, "counts.txt"))

  meta <- data.frame(
    dbname = c(RiboPilot:::RP_DECOY_DBS, "taxdb"),
    `last-updated` = c("2026-08-15T00:00:00", "2022-12-05T00:00:00",
      "2022-12-05T00:00:00", "2026-08-13T16:23:29.394961"),
    `number-of-sequences` = c(unname(RP_DECOY_SEQS[RiboPilot:::RP_DECOY_DBS]), 0L),
    check.names = FALSE, stringsAsFactors = FALSE)
  # a present entry with a null field, which jq -r renders as the string "null"
  meta$`last-updated`[meta$dbname %in% meta_null] <- NA_character_
  writeLines(jsonlite::toJSON(meta, auto_unbox = FALSE, na = "null"),
    file.path(bindir, "meta.json"))
  md5 <- c("8ea926b50cb87831abd58e97758b37ff", "f47c1bad5051bda9c5c00a236087c62b",
    "4f8310a6e3e344017594928ce0b05b9a", "5fc70fa3d152fee1608f134ca4cb1507")
  writeLines(sprintf("%s  %s.tar.gz", md5, c(RiboPilot:::RP_DECOY_DBS, "taxdb")),
    file.path(bindir, "md5.txt"))

  writeLines(c(
    "#!/usr/bin/env bash",
    "set -e",
    'here="$(cd "$(dirname "$0")" && pwd)"',
    'url=""; out=""; head=0',
    'while [ $# -gt 0 ]; do',
    '  case "$1" in',
    '    -o) out="$2"; shift 2 ;;',
    "    -sI | -Is | -I) head=1; shift ;;",
    "    -*) shift ;;",
    '    *) url="$1"; shift ;;',
    "  esac",
    "done",
    '[ "$head" != 1 ] || { printf "Last-Modified: Sat, 15 Aug 2026 09:36:25 GMT\\r\\n"; exit 0; }',
    'case "$url" in',
    '  *blastdb-metadata-1-1.json) body="$(cat "$here/meta.json")" ;;',
    '  *.tar.gz.md5) name="${url##*/}"; body="$(grep -F "  ${name%.md5}" "$here/md5.txt")" ;;',
    '  *) echo "stub curl: unexpected url $url" >&2; exit 22 ;;',
    "esac",
    '[ -n "$out" ] || { printf "%s\\n" "$body"; exit 0; }',
    'printf "%s\\n" "$body" > "$out"'
  ), file.path(bindir, "curl"))

  Sys.chmod(list.files(bindir, full.names = TRUE, pattern = "^(curl|update_blastdb.pl)$"),
    "0755")
  bindir
}

run_fetch_decoys <- function(dest, bindir, brk = "") {
  out <- suppressWarnings(system2("bash",
    c(shQuote(normalizePath(decoys_sh)), shQuote(dest)),
    stdout = TRUE, stderr = TRUE,
    env = c(paste0("PATH=", bindir, ":", Sys.getenv("PATH")),
      paste0("RP_STUB_BREAK=", brk))))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out)
}

skip_without_decoy_tools <- function() {
  skip_on_cran()
  skip_if_not(file.exists(decoys_sh),
    "dev/ is Rbuildignored; this runs from the source tree only")
  for (tool in c("bash", "jq", "makeblastdb", "blastdbcmd")) {
    skip_if(Sys.which(tool) == "", paste(tool, "is not on PATH"))
  }
}

test_that("fetch_decoys.sh verifies what landed and records the source's own versions", {
  skip_without_decoy_tools()
  dest <- withr::local_tempdir()
  r <- run_fetch_decoys(dest, decoy_stub_bin())
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  for (f in RP_TAXONOMY_FILES) expect_true(file.exists(file.path(dest, f)), label = f)

  v <- readLines(file.path(dest, "DECOY_VERSION.txt"), warn = FALSE)
  expect_equal(sub("\t.*", "", v),
    c("# database", RiboPilot:::RP_DECOY_DBS, "taxdb"))
  # every field comes from the source, none from the marker left on disk
  expect_equal(v[2], paste("16S_ribosomal_RNA", "8ea926b50cb87831abd58e97758b37ff",
    "Sat, 15 Aug 2026 09:36:25 GMT", "2026-08-15T00:00:00", "27644", sep = "\t"))
  expect_equal(v[5], paste("taxdb", "5fc70fa3d152fee1608f134ca4cb1507",
    "Sat, 15 Aug 2026 09:36:25 GMT", "2026-08-13T16:23:29.394961", "0", sep = "\t"))
})

test_that("fetch_decoys.sh refuses whatever the fetch failed to leave behind", {
  skip_without_decoy_tools()
  cases <- list(
    # update_blastdb.pl prints "up to date" and installs nothing whenever the
    # .tar.gz.md5 marker is there, so the marker can never stand in for a check
    uptodate = "is missing or empty after fetching",
    missing_volume = "16S_ribosomal_RNA.nsq is missing or empty",
    empty_volume = "16S_ribosomal_RNA.nsq is missing or empty",
    bad_index = "blastdbcmd -info failed",
    # complete, openable and one release out of date: every other check passes
    # and the live version would otherwise be written beside the old database
    stale = "holds 27649 sequence(s) but the source publishes 27644",
    no_taxonomy = "taxonomy4blast.sqlite3 is missing or empty")
  for (brk in names(cases)) {
    dest <- withr::local_tempdir()
    # a previous run's record must survive a failed one rather than be replaced
    # by a half-written file the packager would then ship
    writeLines("sentinel", file.path(dest, "DECOY_VERSION.txt"))
    r <- run_fetch_decoys(dest, decoy_stub_bin(), brk)
    expect_gt(r$status, 0L)
    expect_true(any(grepl(cases[[brk]], r$out, fixed = TRUE)), label = brk)
    expect_equal(readLines(file.path(dest, "DECOY_VERSION.txt"), warn = FALSE),
      "sentinel", label = brk)
  }
})

test_that("fetch_decoys.sh refuses a null metadata field rather than recording it as a version", {
  skip_without_decoy_tools()
  dest <- withr::local_tempdir()
  r <- run_fetch_decoys(dest, decoy_stub_bin(meta_null = "16S_ribosomal_RNA"))
  expect_gt(r$status, 0L)
  expect_true(any(grepl("has no last-updated", r$out, fixed = TRUE)))
  expect_false(file.exists(file.path(dest, "DECOY_VERSION.txt")))
})

# --- the two stages in build_refdb.sh ----------------------------------------

# Sourcing build_refdb.sh leaves BUILD_DIR, RP_STATE_FILE, VER and OUT_DIR
# empty, so the driver sets all four, the same way the build stage's harness in
# test-refdb-build.R does.
run_stage_sourced <- function(stage, build, ver = "9.9.9", path_prefix = NULL,
                              env_extra = character(0)) {
  driver <- file.path(build, "drive_stage.sh")
  writeLines(c(
    "set -uo pipefail",
    'refdb="$1"; b="$2"',
    "shift 2",
    "set --",
    'source "$refdb"',
    "set +Ee",
    "trap - ERR",
    'BUILD_DIR="$b"',
    'RP_STATE_FILE="$BUILD_DIR/state.txt"',
    sprintf('VER="%s"', ver),
    'OUT_DIR="$BUILD_DIR/blastdb/$VER"',
    stage,
    'echo "RC=$?"'
  ), driver)
  env <- c(env_extra, if (is.null(path_prefix)) character(0) else {
    paste0("PATH=", path_prefix, ":", Sys.getenv("PATH"))
  })
  out <- suppressWarnings(system2("timeout", c("120", "bash", shQuote(driver),
    shQuote(normalizePath(refdb_sh)), shQuote(build)),
    stdout = TRUE, stderr = TRUE, env = env))
  st <- attr(out, "status")
  status <- if (is.null(st)) 0L else as.integer(st)
  expect_false(status == 124L, label = paste(stage, "timed out"))
  list(status = status, out = out, rc = grep("^RC=", out, value = TRUE))
}

test_that("the package stage stages both assets, verifies them and records its marker", {
  skip_without_packager()
  skip_if(Sys.which("timeout") == "", "timeout is not on PATH")
  src <- mock_panel_dir()
  build <- dirname(dirname(src))
  r <- run_stage_sourced("stage_package", build)
  expect_equal(r$rc, "RC=0", label = paste(r$out, collapse = "\n"))

  rel <- file.path(build, "blastdb", "release-9.9.9")
  for (a in c("ribopilot_rrna-9.9.9.tar.gz", "decoys-9.9.9.tar.gz")) {
    expect_true(file.exists(file.path(rel, a)), label = a)
    expect_true(file.exists(file.path(rel, paste0(a, ".sha256"))), label = a)
  }
  state <- readLines(file.path(build, "state.txt"), warn = FALSE)
  expect_match(state, "^package\t[0-9a-f]{64}\t")
})

test_that("the package stage refuses when the build stage has not run", {
  skip_without_packager()
  skip_if(Sys.which("timeout") == "", "timeout is not on PATH")
  build <- withr::local_tempdir()
  r <- run_stage_sourced("stage_package", build)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("run the build stage first", r$out, fixed = TRUE)))
})

test_that("the decoys stage fetches into the panel directory and records what landed", {
  skip_without_decoy_tools()
  skip_if(Sys.which("timeout") == "", "timeout is not on PATH")
  src <- mock_panel_dir(drop = c(RP_TAXONOMY_FILES, "DECOY_VERSION.txt"))
  build <- dirname(dirname(src))
  r <- run_stage_sourced("stage_decoys", build, path_prefix = decoy_stub_bin())
  expect_equal(r$rc, "RC=0", label = paste(r$out, collapse = "\n"))
  for (f in c(RP_TAXONOMY_FILES, "DECOY_VERSION.txt")) {
    expect_true(file.exists(file.path(src, f)), label = f)
  }
  expect_match(readLines(file.path(build, "state.txt"), warn = FALSE),
    "^decoys\t[0-9a-f]{64}\t")
})

test_that("the decoys stage fails, with no marker, when verification fails", {
  skip_without_decoy_tools()
  skip_if(Sys.which("timeout") == "", "timeout is not on PATH")
  src <- mock_panel_dir(drop = c(RP_TAXONOMY_FILES, "DECOY_VERSION.txt"))
  build <- dirname(dirname(src))
  r <- run_stage_sourced("stage_decoys", build, path_prefix = decoy_stub_bin(),
    env_extra = "RP_STUB_BREAK=no_taxonomy")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("fetch_decoys.sh failed", r$out, fixed = TRUE)))
  # a green marker on an unverifiable decoy set is the whole failure mode
  expect_false(file.exists(file.path(build, "state.txt")))
})

# Every REQUIRED_TOOLS entry preflight looks for that neither of these two
# stages actually runs. curl is left alone: preflight demands that it resolve
# inside the repo shim directory.
write_preflight_stubs <- function(.envir = parent.frame()) {
  bindir <- withr::local_tempdir(.local_envir = .envir)
  for (t in c("esearch", "efetch", "epost", "esummary", "xtract",
              "update_blastdb.pl", "makeblastdb", "blastdbcmd", "seqkit")) {
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bindir, t))
  }
  Sys.chmod(list.files(bindir, full.names = TRUE), "0755")
  bindir
}

run_refdb_cli <- function(build, args, bindir) {
  out <- suppressWarnings(system2("timeout",
    c("120", "bash", shQuote(normalizePath(refdb_sh)),
      "--build-dir", shQuote(build), args),
    stdout = TRUE, stderr = TRUE,
    env = c(paste0("PATH=", bindir, ":", Sys.getenv("PATH")), "NCBI_API_KEY=stub")))
  st <- attr(out, "status")
  status <- if (is.null(st)) 0L else as.integer(st)
  expect_false(status == 124L, label = "the CLI run timed out instead of finishing")
  list(status = status, out = out)
}

test_that("the package stage runs on a built panel without --force, and the panel stages still refuse", {
  skip_without_packager()
  skip_if(Sys.which("timeout") == "", "timeout is not on PATH")
  ver <- rp_db_version_string()
  src <- mock_panel_dir(ver = ver)
  build <- dirname(dirname(src))
  bindir <- write_preflight_stubs()

  # the normal case for both new stages is a panel that is already built, so
  # requiring --force there would mean saying "rebuild the panel in place"
  r <- run_refdb_cli(build, c("--stage", "package"), bindir)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  expect_false(any(grepl("already exists", r$out, fixed = TRUE)))
  expect_true(file.exists(file.path(build, "blastdb", paste0("release-", ver),
    paste0("ribopilot_rrna-", ver, ".tar.gz"))))

  # the refusal itself is untouched for anything that does rebuild the panel
  for (a in list(c("--stage", "build"), character(0))) {
    r <- run_refdb_cli(build, a, bindir)
    expect_gt(r$status, 0L)
    expect_true(any(grepl("already exists", r$out, fixed = TRUE)),
      label = paste(a, collapse = " "))
  }
})

# --- docker and the installer ------------------------------------------------

test_that("the Dockerfile and the deploy scripts name the assets the packager produces", {
  skip_without_packager()
  dockerfile <- "../../docker/Dockerfile"
  skip_if_not(file.exists(dockerfile),
    "docker/ is Rbuildignored; this runs from the source tree only")
  ver <- rp_db_version_string()
  src <- mock_panel_dir(ver = ver)
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  produced <- list.files(r$rel, pattern = "[.]tar[.]gz$")
  expect_length(produced, 2L)

  for (f in c(dockerfile, "../../docker/deploy-dockerhub.sh", "../../docker/deploy-local.sh")) {
    txt <- gsub("${RP_DB_VERSION}", ver,
      paste(readLines(f, warn = FALSE), collapse = "\n"), fixed = TRUE)
    for (a in produced) expect_match(txt, a, fixed = TRUE, label = paste(f, a))
  }
  df <- paste(readLines(dockerfile, warn = FALSE), collapse = "\n")
  expect_match(df, "ENV RIBOPILOT_DB_DIR=/opt/ribopilot/blastdb", fixed = TRUE)
})

test_that("the deploy scripts can build an image with no panel at all", {
  # the release the scripts stage from may not exist yet; that must not block
  # every image build, and a panel-less image must not advertise a panel
  skip_on_cran()
  dockerfile <- "../../docker/Dockerfile"
  skip_if_not(file.exists(dockerfile), "docker/ is Rbuildignored; this runs from the source tree only")
  df <- readLines(dockerfile, warn = FALSE)
  expect_true(any(grepl("^FROM db-\\$\\{RP_DB\\} AS final", df)))
  # RIBOPILOT_DB_DIR is set only in the stage that actually carries the panel
  db_env <- grep("^ENV RIBOPILOT_DB_DIR=", df)
  expect_length(db_env, 1)
  expect_gt(db_env, max(grep("^FROM base AS db-1", df)))

  for (s in c("../../docker/deploy-local.sh", "../../docker/deploy-dockerhub.sh")) {
    ds <- paste(readLines(s, warn = FALSE), collapse = "\n")
    expect_match(ds, "--no-db", fixed = TRUE, info = s)
    expect_match(ds, "RP_DB=${RP_DB}", fixed = TRUE, info = s)
    # gh, not curl: the repo is private, so an unauthenticated fetch 404s
    expect_match(ds, "gh release download", fixed = TRUE, info = s)
    expect_false(grepl("curl -fsSL", ds, fixed = TRUE), info = s)
  }
})

test_that("the deploy script checksum-verifies the db panel tarballs before building", {
  skip_on_cran()
  deploy <- "../../docker/deploy-dockerhub.sh"
  skip_if_not(file.exists(deploy), "docker/ is Rbuildignored; this runs from the source tree only")
  ds <- paste(readLines(deploy, warn = FALSE), collapse = "\n")
  # fetches the published checksum alongside the tarball
  expect_match(ds, ".sha256", fixed = TRUE)
  # and actually runs the verification, not just downloads it
  expect_match(ds, '-c "$sumfile"', fixed = TRUE)
})

# --- docker build wiring -----------------------------------------------------

RP_DOCKERFILE <- "../../docker/Dockerfile"
RP_DEPLOY_SH <- c("../../docker/deploy-local.sh", "../../docker/deploy-dockerhub.sh")

skip_without_docker_dir <- function() {
  skip_on_cran()
  skip_if_not(file.exists(RP_DOCKERFILE),
    "docker/ is Rbuildignored; this runs from the source tree only")
  skip_if(Sys.which("bash") == "", "bash is not on PATH")
}

# The sweep is lifted out of the deploy script and run for real against a staged
# docker/: a grep for the rm would pass just as well on a loop that deleted the
# current version's cache, which is the failure this loop exists to avoid.
# Lifted from the empty-version guard, not from the loop: with an empty
# RP_DB_VERSION every case pattern becomes *"-.tar.gz", matches nothing, and the
# loop deletes the entire 160 MB cache, so the guard is part of the sweep.
sweep_block <- function(path) {
  x <- readLines(path, warn = FALSE)
  i <- grep('^\\[ -n "\\$RP_DB_VERSION" \\]', x)[1]
  if (is.na(i)) i <- grep("^for f in docker/ribopilot_rrna-", x)[1]
  if (is.na(i)) return(character(0))
  ends <- grep("^done$", x)
  x[i:ends[ends > i][1]]
}

run_sweep <- function(path, root, ver) {
  script <- file.path(root, "sweep.sh")
  writeLines(c("set -euo pipefail", sprintf('cd "%s"', root),
    sprintf('RP_DB_VERSION="%s"', ver), sweep_block(path)), script)
  out <- suppressWarnings(system2("bash", shQuote(script),
    stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out)
}

test_that("the deploy scripts sweep staged assets from other panel versions only", {
  skip_without_docker_dir()
  # the current version's two tarballs are a 160 MB cache the scripts deliberately
  # reuse and checksum-verify, so sweeping them would make every build re-download
  keep <- c("ribopilot_rrna-9.9.9.tar.gz", "ribopilot_rrna-9.9.9.tar.gz.sha256",
    "decoys-9.9.9.tar.gz", "decoys-9.9.9.tar.gz.sha256",
    "RiboPilot_0.0.0.tar.gz", "Dockerfile")
  drop <- c("ribopilot_rrna-8.8.8.tar.gz", "ribopilot_rrna-8.8.8.tar.gz.sha256",
    "decoys-8.8.8.tar.gz",
    # a sidecar whose tarball is already gone still carries a stale version
    "decoys-7.7.7.tar.gz.sha256")
  for (s in RP_DEPLOY_SH) {
    skip_if_not(file.exists(s), paste(s, "is Rbuildignored"))
    root <- withr::local_tempdir()
    dir.create(file.path(root, "docker"))
    for (f in c(keep, drop)) writeLines(f, file.path(root, "docker", f))
    r <- run_sweep(s, root, "9.9.9")
    expect_equal(r$status, 0L, label = paste(s, paste(r$out, collapse = "\n")))
    expect_setequal(list.files(file.path(root, "docker")), keep)
    expect_setequal(r$out, paste("removing stale staged asset", drop))

    # nothing staged at all: an unmatched glob must not become a filename, and
    # must not be announced as a removal either
    bare <- withr::local_tempdir()
    dir.create(file.path(bare, "docker"))
    writeLines("x", file.path(bare, "docker", "Dockerfile"))
    r <- run_sweep(s, bare, "9.9.9")
    expect_equal(r$status, 0L, label = paste(s, paste(r$out, collapse = "\n")))
    expect_equal(list.files(file.path(bare, "docker")), "Dockerfile")
    expect_length(r$out, 0L)

    # an empty RP_DB_VERSION (constant missing or renamed in R/db_dist.R) must
    # abort, not sweep: every case pattern would match nothing and the whole
    # staged cache would go, including on a --no-db build
    empty <- withr::local_tempdir()
    dir.create(file.path(empty, "docker"))
    for (f in c(keep, drop)) writeLines(f, file.path(empty, "docker", f))
    r <- run_sweep(s, empty, "")
    expect_gt(r$status, 0L, label = s)
    expect_true(any(grepl("RP_DB_VERSION not found", r$out, fixed = TRUE)), label = s)
    expect_setequal(list.files(file.path(empty, "docker")), c(keep, drop))

    # after RP_DB_VERSION is known and before the fetch, so the sweep can tell
    # which version is current and cannot delete what the fetch just verified
    x <- readLines(s, warn = FALSE)
    at <- grep("^for f in docker/ribopilot_rrna-", x)
    expect_length(at, 1L)
    expect_gt(at, grep("^RP_DB_VERSION=[$][(]awk", x))
    expect_lt(at, grep('^  for A in "ribopilot_rrna-', x))
  }
})

# The verification RUN as the Dockerfile actually carries it, pointed at a local
# directory instead of the image path.
dockerfile_db_check <- function(dir) {
  x <- readLines(RP_DOCKERFILE, warn = FALSE)
  i <- grep("^RUN blastdbcmd ", x)[1]
  if (is.na(i)) return("false")
  j <- i
  while (grepl("\\\\$", x[j])) j <- j + 1L
  gsub("/opt/ribopilot/blastdb/${RP_DB_VERSION}", dir,
    sub("^RUN +", "", paste(x[i:j], collapse = "\n")), fixed = TRUE)
}

run_db_check <- function(dir) {
  out <- suppressWarnings(system2("bash", c("-c", shQuote(dockerfile_db_check(dir))),
    stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  if (is.null(st)) 0L else as.integer(st)
}

test_that("the panel verification sits in db-1, after both ADDs, with no pipe", {
  skip_without_docker_dir()
  x <- readLines(RP_DOCKERFILE, warn = FALSE)
  at <- grep("^RUN blastdbcmd ", x)
  expect_length(at, 1L)
  adds <- grep("^ADD docker/", x)
  expect_length(adds, 2L)
  expect_gt(at, max(adds))
  # nothing but comments between the ADDs and the check: another instruction
  # there would be a layer that got built on an unverified database
  expect_true(all(grepl("^\\s*(#|$)", x[(max(adds) + 1L):(at - 1L)])))
  # and still inside db-1, above the ENV that advertises the panel
  expect_lt(at, grep("^ENV RIBOPILOT_DB_DIR=", x))
  # a piped command returns 141 on this machine even when it succeeds, which
  # would pass the check for the wrong reason
  expect_false(grepl("|", dockerfile_db_check("/nowhere"), fixed = TRUE))
})

test_that("the Dockerfile's panel check accepts a real database and refuses a broken one", {
  skip_without_docker_dir()
  for (tool in c("makeblastdb", "blastdbcmd")) {
    skip_if(Sys.which(tool) == "", paste(tool, "is not on PATH"))
  }
  d <- withr::local_tempdir()
  fa <- file.path(d, "in.fasta")
  writeLines(c(">s1", "ACGTACGTACGTACGTACGTACGT", ">s2", "TTTTGGGGCCCCAAAATTTTGGGG"), fa)
  expect_equal(system2("makeblastdb",
    c("-in", shQuote(fa), "-dbtype", "nucl", "-blastdb_version", "5",
      "-out", shQuote(file.path(d, "ribopilot_rrna"))),
    stdout = FALSE, stderr = FALSE), 0L)
  writeLines('{"db_name": "ribopilot_rrna"}', file.path(d, "manifest.json"))
  expect_equal(run_db_check(d), 0L)

  # ADD plain-copies a corrupt-but-gzip-shaped archive instead of failing, so the
  # destination holds the tarball and no volumes at all
  vols <- list.files(d, pattern = "^ribopilot_rrna[.]", full.names = TRUE)
  stash <- withr::local_tempdir()
  expect_true(all(file.rename(vols, file.path(stash, basename(vols)))))
  writeLines("gzip-shaped junk", file.path(d, "ribopilot_rrna-9.9.9.tar.gz"))
  expect_gt(run_db_check(d), 0L)
  expect_true(all(file.rename(file.path(stash, basename(vols)), vols)))
  expect_equal(run_db_check(d), 0L)

  # a volume set that is present but unreadable: blastdbcmd -info is what tells
  # these apart, and it is the half of the check a file test cannot do
  idx <- file.path(d, "ribopilot_rrna.nin")
  intact <- readBin(idx, "raw", file.size(idx))
  writeBin(charToRaw("not a blast index"), idx)
  expect_gt(run_db_check(d), 0L)
  writeBin(intact, idx)
  expect_equal(run_db_check(d), 0L)

  # an openable database whose manifest never landed: -s, not -e, because an
  # empty file is what a truncated extraction leaves
  writeLines(character(0), file.path(d, "manifest.json"))
  expect_gt(run_db_check(d), 0L)
  expect_true(file.remove(file.path(d, "manifest.json")))
  expect_gt(run_db_check(d), 0L)
})
