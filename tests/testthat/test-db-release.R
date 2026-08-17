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

# A panel directory holding exactly what the build and decoys stages leave: the
# panel volume set, a manifest, the three decoy databases, the taxonomy files
# and the version record. Files are named, not built; nothing here reads them.
mock_panel_dir <- function(ver = "9.9.9", drop = character(0), extra = character(0),
                           manifest_ver = ver, .envir = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = .envir)
  d <- file.path(root, "blastdb", ver)
  dir.create(d, recursive = TRUE)
  files <- c(
    paste0("ribopilot_rrna.", RP_PANEL_EXT),
    paste0(rep(RiboPilot:::RP_DECOY_DBS, each = 4), ".", c("ndb", "nhr", "nin", "nsq")),
    RP_TAXONOMY_FILES, "DECOY_VERSION.txt")
  for (f in c(setdiff(files, drop), extra)) writeLines(f, file.path(d, f))
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
  expect_setequal(got, c(
    paste0(rep(RiboPilot:::RP_DECOY_DBS, each = 4), ".", c("ndb", "nhr", "nin", "nsq")),
    RP_TAXONOMY_FILES, "DECOY_VERSION.txt"))
  # taxonomy4blast.sqlite3 is not a taxdb.* file, and shipping without it makes
  # -taxids return unfiltered rows while blastn still exits 0
  expect_true("taxonomy4blast.sqlite3" %in% got)
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
    expect_false(file.exists(file.path(r$rel, "decoys-9.9.9.tar.gz")), label = f)
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
  expect_false(file.exists(file.path(r$rel, "ribopilot_rrna-7.7.7.tar.gz")))
})

test_that("packaging refuses while a decoy archive is still sitting in the panel", {
  skip_without_packager()
  src <- mock_panel_dir(extra = "16S_ribosomal_RNA.tar.gz")
  r <- run_packager(src)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("the decoy fetch did not finish", r$out, fixed = TRUE)))
})

test_that("the panel tarball ships the volume set and the manifest and nothing else", {
  skip_without_packager()
  src <- mock_panel_dir(extra = c("ribopilot_rrna.fasta", "ribopilot_rrna.taxid.map",
    "ribopilot_rrna.dropped.tsv"))
  r <- run_packager(src)
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))
  expect_setequal(tar_contents(file.path(r$rel, "ribopilot_rrna-9.9.9.tar.gz")),
    c("manifest.json", paste0("ribopilot_rrna.", RP_PANEL_EXT)))
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
})

# --- fetch_decoys.sh ---------------------------------------------------------

# update_blastdb.pl and the three NCBI endpoints are stubbed; makeblastdb and
# blastdbcmd are the real ones, because the verification this stage exists for
# is "does blastdbcmd accept what landed", which a stubbed blastdbcmd cannot
# answer. RP_STUB_BREAK damages exactly one thing so one guard is under test.
decoy_stub_bin <- function(.envir = parent.frame()) {
  bindir <- withr::local_tempdir(.local_envir = .envir)
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -e",
    'db="${!#}"',
    'brk="${RP_STUB_BREAK:-}"',
    '[ "$brk" != uptodate ] || { echo "$db is up to date"; exit 0; }',
    'if [ "$db" = taxdb ]; then',
    '  for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do echo stub > "$f"; done',
    '  [ "$brk" != no_taxonomy ] || rm -f taxonomy4blast.sqlite3',
    "  exit 0",
    "fi",
    'printf ">seq1 stub\\nACGTACGTACGTACGTACGTACGTACGT\\n" > "$db.stub.fasta"',
    'makeblastdb -in "$db.stub.fasta" -dbtype nucl -blastdb_version 5 -out "$db" > /dev/null',
    'rm -f "$db.stub.fasta"',
    'case "$brk" in',
    "  missing_volume) rm -f 16S_ribosomal_RNA.nsq ;;",
    "  empty_volume) : > 16S_ribosomal_RNA.nsq ;;",
    '  bad_index) printf "not a blast index" > 16S_ribosomal_RNA.nin ;;',
    "esac"
  ), file.path(bindir, "update_blastdb.pl"))

  meta <- data.frame(
    dbname = c(RiboPilot:::RP_DECOY_DBS, "taxdb"),
    `last-updated` = c("2026-08-15T00:00:00", "2022-12-05T00:00:00",
      "2022-12-05T00:00:00", "2026-08-13T16:23:29.394961"),
    `number-of-sequences` = c(27644L, 8784L, 6575L, 0L),
    check.names = FALSE, stringsAsFactors = FALSE)
  writeLines(jsonlite::toJSON(meta, auto_unbox = FALSE), file.path(bindir, "meta.json"))
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
  # publishing is manual: the stage must never reach gh release create
  expect_false(any(grepl("gh release", r$out, fixed = TRUE)))
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
