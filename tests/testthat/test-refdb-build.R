refdb_dir <- "../../dev/refdb"
refdb_sh <- file.path(refdb_dir, "build_refdb.sh")
refdb_lib <- file.path(refdb_dir, "lib.sh")

# dev/ is in .Rbuildignore, so these guards belong inside the tests that need
# them: at file level a single skip drops every test in the file, including
# ones that never touch dev/.

test_that("build_refdb.sh declares the exact Entrez query the panel is defined by", {
  skip_if_not(file.exists(refdb_sh))
  s <- paste(readLines(refdb_sh, warn = FALSE), collapse = "\n")
  expect_match(s, "txid33208[Organism:exp]", fixed = TRUE)
  expect_match(s, "NOT mitochondrion[Filter]", fixed = TRUE)
  expect_match(s, "NOT plastid[Filter]", fixed = TRUE)
  # the filter that silently removed every RefSeq record must never reappear
  expect_false(grepl("efilter", s, fixed = TRUE))
})

test_that("every refdb script sets the strict shell options and C locale", {
  skip_if_not(dir.exists(refdb_dir))
  files <- list.files(refdb_dir, pattern = "[.]sh$", full.names = TRUE)
  expect_true(length(files) > 0)
  for (f in files) {
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_match(s, "set -Eeuo pipefail", fixed = TRUE, label = basename(f))
    expect_match(s, "LC_ALL=C", fixed = TRUE, label = basename(f))
  }
})

test_that("the docker ARG default agrees with RP_DB_VERSION", {
  skip_if_not(file.exists("../../docker/Dockerfile"))
  v <- sub('.*"(.*)".*', "\\1",
    grep("^RP_DB_VERSION", readLines("../../R/db_dist.R", warn = FALSE), value = TRUE)[1])
  arg <- grep("^ARG RP_DB_VERSION", readLines("../../docker/Dockerfile", warn = FALSE), value = TRUE)[1]
  expect_equal(sub(".*=", "", arg), v)
})

# --- lock -------------------------------------------------------------------

test_that("two concurrent rp_lock runs never both acquire the lock", {
  skip_if_not(file.exists(refdb_lib))
  skip_if(Sys.which("bash") == "")
  work <- withr::local_tempdir()
  harness <- file.path(work, "race.sh")
  writeLines(c(
    "set -uo pipefail",
    'source "$1"',
    "set +Ee",
    "trap - ERR",
    "bad=0",
    'for i in $(seq 1 "$3"); do',
    '  d="$2/p$i"',
    '  mkdir -p "$d"',
    "  for r in a b; do",
    "    (",
    '      while [ ! -e "$d/go" ]; do :; done',
    '      rp_lock "$d/lock" >/dev/null 2>&1 && : > "$d/won.$r"',
    "    ) &",
    "  done",
    "  sleep 0.05",
    '  : > "$d/go"',
    "  wait",
    "  n=0",
    "  for r in a b; do",
    '    if [ -e "$d/won.$r" ]; then n=$((n + 1)); fi',
    "  done",
    '  if [ "$n" -ne 1 ]; then bad=$((bad + 1)); fi',
    "done",
    'echo "bad=$bad"'
  ), harness)
  out <- system2("bash", c(
    shQuote(harness), shQuote(normalizePath(refdb_lib)), shQuote(work), "40"
  ), stdout = TRUE, stderr = TRUE)
  expect_equal(grep("^bad=", out, value = TRUE), "bad=0")
})

test_that("rp_lock recovers from a dead recorded pid and from an empty lockfile", {
  skip_if_not(file.exists(refdb_lib))
  skip_if(Sys.which("bash") == "")
  work <- withr::local_tempdir()
  dead <- system2("bash", c("-c", shQuote("echo $$")), stdout = TRUE)
  try_lock <- function(contents, name) {
    lock <- file.path(work, name)
    writeLines(contents, lock)
    script <- sprintf(
      'source %s; rp_lock %s >/dev/null 2>&1; echo "rc=$?"',
      shQuote(normalizePath(refdb_lib)), shQuote(lock))
    system2("bash", c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
  }
  expect_true(any(grepl("rc=0", try_lock(dead, "dead.pid"), fixed = TRUE)))
  expect_true(any(grepl("rc=0", try_lock(character(0), "empty.pid"), fixed = TRUE)))
})

# --- secret scrubbing -------------------------------------------------------

test_that("the api key never reaches stderr, in any position it really appears in", {
  skip_if_not(file.exists(refdb_lib))
  skip_if(Sys.which("bash") == "")
  secret <- "d3adb33fcafe0123456789ab"
  work <- withr::local_tempdir()
  driver <- file.path(work, "emit.sh")
  writeLines(c(
    'source "$1"',
    "rp_scrub_output",
    "trap rp_scrub_flush EXIT",
    'K="$2"',
    'echo "api_key=$K" >&2',
    'echo "fetching with api_key=$K now" >&2',
    'echo "https://eutils.ncbi.nlm.nih.gov/x.fcgi?api_key=$K" >&2',
    'echo "https://eutils.ncbi.nlm.nih.gov/x.fcgi?db=nuccore&api_key=$K" >&2',
    'echo "two api_key=$K and api_key=$K on one line" >&2',
    'echo "https://x?api_key=$K&retmax=500" >&2',
    'echo "api_key=$K-suffix" >&2',
    'echo "NCBI_API_KEY=$K" >&2',
    'export NCBI_API_KEY="$K"',
    'export -p | grep NCBI_API_KEY >&2 || true',
    'rp_log "esearch failed with api_key=$K"',
    "set -x",
    ': "url?api_key=$K"'
  ), driver)
  out <- system2("bash", c(shQuote(driver), shQuote(normalizePath(refdb_lib)), secret),
    stdout = TRUE, stderr = TRUE)
  expect_gt(length(out), 10)
  expect_false(any(grepl(secret, out, fixed = TRUE)))
  expect_true(any(grepl("api_key=REDACTED", out, fixed = TRUE)))
})

test_that("build_refdb.sh scrubs both spellings of the key on its own stderr", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  secret <- "d3adb33fcafe0123456789ab"
  for (arg in c(paste0("api_key=", secret), paste0("NCBI_API_KEY=", secret))) {
    out <- suppressWarnings(system2("bash", c(shQuote(refdb_sh), shQuote(arg)),
      stdout = TRUE, stderr = TRUE))
    expect_false(any(grepl(secret, out, fixed = TRUE)), label = arg)
    expect_true(any(grepl("REDACTED", out, fixed = TRUE)), label = arg)
  }
})

# --- check-version / clean --------------------------------------------------

rp_db_version <- function() {
  sub('.*"(.*)".*', "\\1",
    grep("^RP_DB_VERSION", readLines("../../R/db_dist.R", warn = FALSE), value = TRUE)[1])
}

run_refdb <- function(...) {
  out <- suppressWarnings(system2("bash", c(shQuote(refdb_sh), ...),
    stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out)
}

test_that("--check-version reports a version directory that exists on disk", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  build <- withr::local_tempdir()
  stale <- file.path(build, "blastdb", "2026.07.9")
  dir.create(stale, recursive = TRUE)
  writeLines('{"db_version": "2026.07.9"}', file.path(stale, "manifest.json"))
  r <- run_refdb("--build-dir", shQuote(build), "--check-version")
  expect_gt(r$status, 0)
  expect_true(any(grepl("MISMATCH", r$out, fixed = TRUE)))
  expect_true(any(grepl("2026.07.9", r$out, fixed = TRUE)))
})

# Build tree with everything clean is meant to look at.
make_build_tree <- function(ver, tarballs = c("good", "corrupt", "none")) {
  tarballs <- match.arg(tarballs)
  build <- withr::local_tempdir(.local_envir = parent.frame())
  panel <- file.path(build, "blastdb", ver)
  rel <- file.path(build, "blastdb", paste0("release-", ver))
  dir.create(file.path(build, "chunks"), recursive = TRUE)
  dir.create(panel, recursive = TRUE)
  writeLines("fa", file.path(build, "chunks", "fa_0001"))
  writeLines("raw", file.path(build, "metazoa_rrna.raw.fasta"))
  for (f in c("ribopilot_rrna.nin", "manifest.json", "ribopilot_rrna.fasta",
              "ribopilot_rrna.taxid.map", "dropped.tsv")) {
    writeLines("x", file.path(panel, f))
  }
  writeLines(paste(
    c("accessions", "taxidmap", "fetch", "build", "decoys", "package"),
    "hash", "2026-08-15T00:00:00Z", sep = "\t"),
    file.path(build, "state.txt"))
  if (tarballs != "none") {
    dir.create(rel, recursive = TRUE)
    for (f in paste0(c("ribopilot_rrna-", "decoys-"), ver, ".tar.gz")) {
      writeLines("tarball", file.path(rel, f))
    }
    system2("bash", c("-c", shQuote(sprintf(
      "cd %s && for f in *.tar.gz; do sha256sum \"$f\" > \"$f.sha256\"; done", shQuote(rel)))))
    if (tarballs == "corrupt") {
      writeLines("tampered", file.path(rel, paste0("decoys-", ver, ".tar.gz")))
    }
  }
  build
}

test_that("clean refuses until the release tarballs verify", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  ver <- rp_db_version()
  for (state in c("none", "corrupt")) {
    build <- make_build_tree(ver, state)
    r <- run_refdb("--build-dir", shQuote(build), "clean")
    expect_gt(r$status, 0)
    expect_true(any(grepl("refusing to clean", r$out, fixed = TRUE)), label = state)
    expect_true(file.exists(file.path(build, "chunks", "fa_0001")), label = state)
  }
})

test_that("clean drops the unshipped artifacts and keeps the surviving stages' markers", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  ver <- rp_db_version()
  build <- make_build_tree(ver, "good")
  panel <- file.path(build, "blastdb", ver)
  r <- run_refdb("--build-dir", shQuote(build), "clean")
  expect_equal(r$status, 0L)
  expect_false(dir.exists(file.path(build, "chunks")))
  expect_false(file.exists(file.path(build, "metazoa_rrna.raw.fasta")))
  expect_false(file.exists(file.path(panel, "ribopilot_rrna.fasta")))
  expect_false(file.exists(file.path(panel, "ribopilot_rrna.taxid.map")))
  expect_false(file.exists(file.path(panel, "dropped.tsv")))
  expect_true(file.exists(file.path(panel, "ribopilot_rrna.nin")))
  expect_true(file.exists(file.path(panel, "manifest.json")))
  expect_true(file.exists(file.path(build, "blastdb", paste0("release-", ver),
    paste0("ribopilot_rrna-", ver, ".tar.gz"))))
  stages <- sub("\t.*", "", readLines(file.path(build, "state.txt"), warn = FALSE))
  expect_setequal(stages, c("accessions", "build", "decoys", "package"))
  expect_false(file.exists(file.path(build, "build.pid")))
})

test_that("clean refuses while a build holds the lock", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  ver <- rp_db_version()
  build <- make_build_tree(ver, "good")
  lock <- file.path(build, "build.pid")
  stop_file <- file.path(build, "stop")
  # holder releases on the stop file, and after 10s regardless, so a failure
  # here cannot leave a spinning process behind
  system2("bash", c("-c", shQuote(sprintf(
    "echo $$ > %s; for i in $(seq 1 200); do if [ -e %s ]; then break; fi; sleep 0.05; done",
    shQuote(lock), shQuote(stop_file)))), wait = FALSE)
  for (i in 1:100) {
    if (file.exists(lock) && length(readLines(lock, warn = FALSE))) break
    Sys.sleep(0.05)
  }
  r <- run_refdb("--build-dir", shQuote(build), "clean")
  file.create(stop_file)
  expect_gt(r$status, 0)
  expect_true(any(grepl("already running", r$out, fixed = TRUE)))
  expect_true(file.exists(file.path(build, "chunks", "fa_0001")))
})
