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
  # BUILD_DIR-root leftovers of the accessions and taxidmap stages
  dir.create(file.path(build, "taxidmap.tmp"))
  for (f in c("accessions.raw.txt", "accessions.esearch.xml", "taxidmap.tsv",
              "taxidmap.sources.tsv", "taxidmap.sha256", "dropped.taxidmap.tsv",
              "accessions.txt", "accessions.sha256")) {
    writeLines("x", file.path(build, f))
  }
  writeLines("x", file.path(build, "taxidmap.tmp", "combined.tsv"))
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
  # both harvest stages' BUILD_DIR-root leftovers go too
  for (f in c("accessions.raw.txt", "accessions.esearch.xml", "taxidmap.tsv",
              "taxidmap.sources.tsv", "taxidmap.sha256", "dropped.taxidmap.tsv")) {
    expect_false(file.exists(file.path(build, f)), label = f)
  }
  expect_false(dir.exists(file.path(build, "taxidmap.tmp")))
  # accessions.txt stays: it is the fingerprint fetch refuses to resume against
  expect_true(file.exists(file.path(build, "accessions.txt")))
  expect_true(file.exists(file.path(build, "accessions.sha256")))
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

# --- accessions stage --------------------------------------------------------

test_that("build_refdb.sh commits to the spec's EXPECTED_COUNT", {
  skip_if_not(file.exists(refdb_sh))
  s <- paste(readLines(refdb_sh, warn = FALSE), collapse = "\n")
  expect_match(s, "EXPECTED_COUNT=892123", fixed = TRUE)
})

test_that("build_refdb.sh filters WGS master accessions with the exact spec pattern", {
  skip_if_not(file.exists(refdb_sh))
  s <- paste(readLines(refdb_sh, warn = FALSE), collapse = "\n")
  expect_match(s, "[A-Z]{4}0{6}", fixed = TRUE)
})

test_that("build_refdb.sh never uses the broken grep -c-with-fallback line-count idiom", {
  skip_if_not(file.exists(refdb_sh))
  s <- paste(readLines(refdb_sh, warn = FALSE), collapse = "\n")
  # on an existing EMPTY file this idiom prints "0" twice, and `[ -eq ]` on
  # the two-line result errors rc=2, silently skipping the guard
  expect_false(grepl("grep -c . ", s, fixed = TRUE))
  expect_true(grepl("wc -l", s, fixed = TRUE))
})

# Sources build_refdb.sh with no positional args (so its own arg-parser sees
# $# == 0 and never runs main_build/dispatch), the same way lib.sh is sourced
# elsewhere in this file, so guard functions can be called directly. `$dir`
# in body_lines is the driver's own directory, a fresh tempdir per call.
run_sourced <- function(body_lines) {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  work <- withr::local_tempdir(.local_envir = parent.frame())
  driver <- file.path(work, "drive.sh")
  writeLines(c(
    "set -uo pipefail",
    'dir="$(cd "$(dirname "$0")" && pwd)"',
    'refdb="$1"',
    "shift",
    "set --",
    'source "$refdb"',
    "set +Ee",
    "trap - ERR",
    body_lines
  ), driver)
  out <- suppressWarnings(system2("bash", c(shQuote(driver), shQuote(normalizePath(refdb_sh))),
    stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out, work = work)
}

test_that("rp_line_count returns a clean single 0 on an empty file", {
  r <- run_sourced(c(
    ': > "$dir/empty.txt"',
    'n=$(rp_line_count "$dir/empty.txt")',
    'echo "N=$n"'
  ))
  expect_equal(r$status, 0L)
  expect_equal(r$out, "N=0")
})

test_that("rp_assert_count aborts on a mismatch derived from an empty file, passes on a match", {
  r_fail <- run_sourced(c(
    ': > "$dir/empty.txt"',
    'n=$(rp_line_count "$dir/empty.txt")',
    'rp_assert_count "test label" "$n" "5"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r_fail$status, 0L)
  expect_false(any(r_fail$out == "UNREACHABLE"))
  expect_true(any(grepl("test label: expected 5, got 0", r_fail$out, fixed = TRUE)))

  r_pass <- run_sourced(c(
    ': > "$dir/empty.txt"',
    'n=$(rp_line_count "$dir/empty.txt")',
    'rp_assert_count "test label" "$n" "0"',
    'echo "REACHED"'
  ))
  expect_equal(r_pass$status, 0L)
  expect_true(any(r_pass$out == "REACHED"))
})

test_that("rp_assert_unique aborts on a duplicate accession, passes on empty input", {
  r_dup <- run_sourced(c(
    'printf "AB123456\nAB123456\n" > "$dir/dup.txt"',
    'rp_assert_unique "test label" "$dir/dup.txt"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r_dup$status, 0L)
  expect_false(any(r_dup$out == "UNREACHABLE"))

  r_empty <- run_sourced(c(
    ': > "$dir/empty.txt"',
    'rp_assert_unique "test label" "$dir/empty.txt"',
    'echo "REACHED"'
  ))
  expect_equal(r_empty$status, 0L)
  expect_true(any(r_empty$out == "REACHED"))
})

test_that("rp_filter_wgs_masters drops matching accessions as a recorded drop, no-op on empty input", {
  r <- run_sourced(c(
    'printf "AAAA00000000\nMZ000001\nBBBB00000000\nAB123456\n" > "$dir/in.txt"',
    'rp_filter_wgs_masters "$dir/in.txt" "$dir/out.txt"',
    'cat "$dir/out.txt"'
  ))
  expect_equal(r$status, 0L)
  kept <- r$out[!grepl("^\\[", r$out)]
  expect_setequal(kept, c("MZ000001", "AB123456"))
  expect_true(any(grepl("dropped 2 WGS master accession", r$out, fixed = TRUE)))

  r_empty <- run_sourced(c(
    ': > "$dir/in.txt"',
    'rp_filter_wgs_masters "$dir/in.txt" "$dir/out.txt"',
    'wc -l < "$dir/out.txt" | tr -d " "'
  ))
  expect_equal(r_empty$status, 0L)
  expect_true(any(grepl("dropped 0 WGS master accession", r_empty$out, fixed = TRUE)))
  expect_true(any(r_empty$out == "0"))
})

test_that("rp_check_drift refuses beyond tolerance without --accept-drift, accepts (and logs) with it", {
  r_refuse <- run_sourced(c(
    'rp_check_drift "800000" "0"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r_refuse$status, 0L)
  expect_false(any(r_refuse$out == "UNREACHABLE"))
  expect_true(any(grepl("accept-drift", r_refuse$out, fixed = TRUE)))

  r_accept <- run_sourced(c(
    'rp_check_drift "800000" "1"',
    'echo "REACHED"'
  ))
  expect_equal(r_accept$status, 0L)
  expect_true(any(r_accept$out == "REACHED"))
  expect_true(any(grepl("accepted via --accept-drift", r_accept$out, fixed = TRUE)))

  r_within <- run_sourced(c(
    'rp_check_drift "902123" "0"',
    'echo "REACHED"'
  ))
  expect_equal(r_within$status, 0L)
  expect_true(any(r_within$out == "REACHED"))
  expect_false(any(grepl("drift", r_within$out, fixed = TRUE)))
})

# Static stub tools for the full-CLI accessions-stage tests: esearch/xtract/
# efetch behave per-test via companion data files; the rest are no-ops that
# only need to exist so run_preflight's tool check passes.
write_accessions_stubs <- function(bindir) {
  dir.create(bindir, showWarnings = FALSE, recursive = TRUE)
  write_stub <- function(name, body) {
    p <- file.path(bindir, name)
    writeLines(c("#!/usr/bin/env bash", "set -eu", body), p)
    Sys.chmod(p, "0755")
  }
  # ENTREZ_DIRECT is esearch's real default chaining wrapper (verified
  # against a live esearch, see task-2-report.md), not raw eSearchResult XML.
  write_stub("esearch", c(
    'self_dir="$(cd "$(dirname "$0")" && pwd)"',
    'printf x >> "$self_dir/esearch_calls.log"',
    'n="$(cat "$self_dir/esearch_count.txt")"',
    "cat <<XML",
    '<ENTREZ_DIRECT><Db>nuccore</Db><WebEnv>NCID_test</WebEnv><QueryKey>1</QueryKey><Count>$n</Count></ENTREZ_DIRECT>',
    "XML"
  ))
  write_stub("xtract", c(
    'file="" pattern=""',
    'while [ $# -gt 0 ]; do case "$1" in -input) file="$2"; shift 2 ;; -pattern) pattern="$2"; shift 2 ;; *) shift ;; esac; done',
    # a -pattern that does not actually wrap the file is exactly the bug this
    # stub exists to catch, so fail instead of parsing anyway
    'grep -q "<$pattern>" "$file"',
    'grep -o "<Count>[0-9]*</Count>" "$file" | head -1 | grep -oE "[0-9]+"'
  ))
  write_stub("efetch", c(
    'self_dir="$(cd "$(dirname "$0")" && pwd)"',
    'cat "$self_dir/efetch_output.txt"'
  ))
  for (noop in c("epost", "esummary", "update_blastdb.pl", "makeblastdb", "blastdbcmd", "seqkit")) {
    write_stub(noop, "exit 0")
  }
  writeLines("0", file.path(bindir, "esearch_count.txt"))
  writeLines(character(0), file.path(bindir, "efetch_output.txt"))
  file.create(file.path(bindir, "esearch_calls.log"))
  bindir
}

run_accessions_stage <- function(count, efetch_lines, extra_args = character(0),
                                 min_accessions = 1L) {
  build <- withr::local_tempdir(.local_envir = parent.frame())
  bindir <- file.path(build, "stubbin")
  write_accessions_stubs(bindir)
  writeLines(as.character(count), file.path(bindir, "esearch_count.txt"))
  writeLines(efetch_lines, file.path(bindir, "efetch_output.txt"))
  env <- c(
    paste0("PATH=", bindir, ":", Sys.getenv("PATH")),
    "NCBI_API_KEY=test-key-not-real"
  )
  # RP_MIN_ACCESSIONS is the test seam for the collapse floor: these fixtures
  # are a handful of records, not ~892,000. min_accessions = NULL leaves it
  # unset, so the stage runs against the real production floor.
  if (!is.null(min_accessions)) {
    env <- c(env, paste0("RP_MIN_ACCESSIONS=", min_accessions))
  }
  out <- suppressWarnings(system2("bash",
    c(shQuote(refdb_sh), "--build-dir", shQuote(build), "--stage", "accessions", extra_args),
    stdout = TRUE, stderr = TRUE, env = env))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out, build = build, bindir = bindir)
}

test_that("accessions stage writes accessions.txt/.sha256 and two distinctly named counts, drops WGS masters, calls esearch once", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  lines <- c("AAAA00000000", "BBBB00000000", "MZ000001", "AB123456",
             "KX999999", "OP111111", "QR222222", "ST333333")
  r <- run_accessions_stage(count = length(lines), efetch_lines = lines,
    extra_args = "--accept-drift")
  expect_equal(r$status, 0L)

  acc_file <- file.path(r$build, "accessions.txt")
  expect_true(file.exists(acc_file))
  expect_setequal(readLines(acc_file), setdiff(lines, c("AAAA00000000", "BBBB00000000")))

  sha_file <- file.path(r$build, "accessions.sha256")
  expect_true(file.exists(sha_file))
  want <- system2("sha256sum", shQuote(acc_file), stdout = TRUE)
  want_hash <- sub("\\s.*", "", want)
  got <- readLines(sha_file)
  expect_true(any(startsWith(got, want_hash)))

  # The pre-filter esearch <Count> and the length of the list are DIFFERENT
  # numbers by construction (the WGS masters are dropped between them), so
  # each gets its own name. A single ambiguous "accessions.count" is what
  # would let a chunking stage request records the list does not contain.
  esearch_n <- as.integer(readLines(file.path(r$build, "accessions.esearch_count")))
  final_n <- as.integer(readLines(file.path(r$build, "accessions.final_count")))
  expect_equal(esearch_n, length(lines))
  expect_equal(final_n, length(readLines(acc_file)))
  expect_equal(esearch_n - final_n, 2L)
  expect_false(file.exists(file.path(r$build, "accessions.count")))

  calls <- readLines(file.path(r$bindir, "esearch_calls.log"), warn = FALSE)
  expect_equal(paste(calls, collapse = ""), "x")
})

# THE REQUIRED TEST: --accept-drift becomes routine (about five weeks of the
# source's measured drift exceeds the tolerance band), so it must not be able
# to promote a collapsed result set. Every other guard passes honestly here --
# 0 == 0 lines, 0 unique, 0 dropped -- which is exactly why the floor has to
# be its own unconditional check.
test_that("an empty harvest is refused even with --accept-drift, and promotes nothing", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  r <- run_accessions_stage(count = 0, efetch_lines = character(0),
    extra_args = "--accept-drift", min_accessions = NULL)
  expect_gt(r$status, 0L)
  expect_false(file.exists(file.path(r$build, "accessions.txt")))
  expect_false(file.exists(file.path(r$build, "accessions.sha256")))
  expect_false(file.exists(file.path(r$build, "accessions.final_count")))
  # no stage marker: state.txt is only ever created by rp_marker_set
  expect_false(file.exists(file.path(r$build, "state.txt")))
  expect_true(any(grepl("collapse floor", r$out, fixed = TRUE)))
})

test_that("the collapse floor also refuses a harvest that is merely tiny, not just empty", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  lines <- c("MZ000001", "AB123456", "KX999999")
  r <- run_accessions_stage(count = length(lines), efetch_lines = lines,
    extra_args = "--accept-drift", min_accessions = NULL)
  expect_gt(r$status, 0L)
  expect_false(file.exists(file.path(r$build, "accessions.txt")))
  expect_true(any(grepl("collapse floor", r$out, fixed = TRUE)))
})

test_that("the accession fingerprint is over the sorted list, so the same set in a different order hashes the same", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  lines <- c("MZ000001", "AB123456", "KX999999", "OP111111")
  hash_of <- function(build) {
    sub("\\s.*", "", readLines(file.path(build, "accessions.sha256"))[1])
  }
  r1 <- run_accessions_stage(count = length(lines), efetch_lines = lines,
    extra_args = "--accept-drift")
  r2 <- run_accessions_stage(count = length(lines), efetch_lines = rev(lines),
    extra_args = "--accept-drift")
  expect_equal(r1$status, 0L)
  expect_equal(r2$status, 0L)
  # efetch history-server order is not contractually stable; the fetch stage
  # refuses to resume on a hash mismatch, so an order-sensitive fingerprint
  # would discard a valid 2 h 20 m fetch
  expect_equal(hash_of(r1$build), hash_of(r2$build))
  expect_equal(readLines(file.path(r1$build, "accessions.txt")),
    sort(lines, method = "radix"))
  expect_equal(readLines(file.path(r2$build, "accessions.txt")),
    sort(lines, method = "radix"))
})

test_that("accessions stage refuses an efetch error line that every other guard lets through", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  # efetch exits 0 with QUERY FAILURE on total failure, and such a line is
  # unique, non-WGS and counts as a line
  lines <- c("MZ000001", "QUERY FAILURE", "AB123456")
  r <- run_accessions_stage(count = length(lines), efetch_lines = lines,
    extra_args = "--accept-drift")
  expect_gt(r$status, 0L)
  expect_false(file.exists(file.path(r$build, "accessions.txt")))
  expect_true(any(grepl("QUERY FAILURE", r$out, fixed = TRUE)))
})

test_that("sourcing build_refdb.sh defines functions and constants and nothing else", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  work <- withr::local_tempdir()
  driver <- file.path(work, "src.sh")
  writeLines(c(
    "set -uo pipefail",
    'dir="$(cd "$(dirname "$0")" && pwd)"',
    # a positional argument must not make the sourced script parse it and exit
    'source "$1" some-positional-argument',
    "set +Ee",
    "trap - ERR",
    'trap -p EXIT > "$dir/exit_traps.txt"',
    'echo "BUILD_DIR=[$BUILD_DIR]"',
    'echo "STATE=[$RP_STATE_FILE]"',
    'echo "SCRUB=[${RP_SCRUB_PIDS:-}]"',
    'echo "EXITTRAPS=$(wc -l < "$dir/exit_traps.txt" | tr -d " ")"',
    'echo "SURVIVED"'
  ), driver)
  out <- suppressWarnings(system2("bash",
    c(shQuote(driver), shQuote(normalizePath(refdb_sh))),
    stdout = TRUE, stderr = TRUE))
  expect_true(any(out == "SURVIVED"))
  # no real-tree default: a sourced stage function cannot write into dev/refdb/build
  expect_true(any(out == "BUILD_DIR=[]"))
  expect_true(any(out == "STATE=[]"))
  # no scrub coprocesses forked into the caller, no EXIT trap installed in it
  expect_true(any(out == "SCRUB=[]"))
  expect_true(any(out == "EXITTRAPS=0"))
})

test_that("a sourced stage function refuses to run rather than writing into a defaulted build tree", {
  r <- run_sourced(c(
    "stage_accessions",
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
  expect_true(any(grepl("BUILD_DIR is not set", r$out, fixed = TRUE)))
})

test_that("accessions stage aborts on a truncated harvest and does not promote a partial accessions.txt", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  full <- c("MZ000001", "AB123456", "KX999999", "OP111111",
            "QR222222", "ST333333", "UV444444", "WX555555")
  r <- run_accessions_stage(count = length(full), efetch_lines = full[-length(full)],
    extra_args = "--accept-drift")
  expect_gt(r$status, 0L)
  expect_false(file.exists(file.path(r$build, "accessions.txt")))
  expect_true(any(grepl("expected 8, got 7", r$out, fixed = TRUE)))
})

test_that("accessions stage refuses a drifted count without --accept-drift, proceeds with it", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  lines <- c("MZ000001", "AB123456")

  r_refuse <- run_accessions_stage(count = length(lines), efetch_lines = lines)
  expect_gt(r_refuse$status, 0L)
  expect_false(file.exists(file.path(r_refuse$build, "accessions.txt")))
  expect_true(any(grepl("accept-drift", r_refuse$out, fixed = TRUE)))

  r_accept <- run_accessions_stage(count = length(lines), efetch_lines = lines,
    extra_args = "--accept-drift")
  expect_equal(r_accept$status, 0L)
  expect_true(file.exists(file.path(r_accept$build, "accessions.txt")))
  expect_true(any(grepl("accepted via --accept-drift", r_accept$out, fixed = TRUE)))
})

# --- taxidmap stage ----------------------------------------------------------

test_that("build_refdb.sh streams the exact three accession2taxid dumps, never landing them", {
  skip_if_not(file.exists(refdb_sh))
  lines <- readLines(refdb_sh, warn = FALSE)
  s <- paste(lines, collapse = "\n")
  expect_match(s, "nucl_gb", fixed = TRUE)
  expect_match(s, "nucl_wgs", fixed = TRUE)
  expect_match(s, "pdb", fixed = TRUE)
  expect_match(s, "accession2taxid", fixed = TRUE)
  # the actual scan code, not just a comment mentioning it
  start <- grep("^rp_taxidmap_scan_dump\\(\\)", lines)
  expect_equal(length(start), 1L)
  end <- start + which(lines[(start + 1):length(lines)] == "}")[1]
  code_lines <- lines[start:end]
  code_lines <- code_lines[!grepl("^[[:space:]]*#", code_lines)]
  expect_true(any(grepl("curl --fail", code_lines, fixed = TRUE)))
  expect_true(any(grepl("zcat", code_lines, fixed = TRUE)))
  expect_true(any(grepl("awk", code_lines, fixed = TRUE)))
})

test_that("build_refdb.sh never suppresses zcat's stderr in the taxidmap scan", {
  skip_if_not(file.exists(refdb_sh))
  lines <- readLines(refdb_sh, warn = FALSE)
  zcat_lines <- grep("zcat", lines, value = TRUE)
  # gzip's own CRC is the truncation detector; suppressing its stderr was the
  # historical defect (356,851 taxid-0 records from one 40% truncation)
  expect_true(length(zcat_lines) > 0)
  expect_false(any(grepl("2>/dev/null", zcat_lines, fixed = TRUE)))
})

test_that("build_refdb.sh's taxid map validator checks the subset before the uniqueness check", {
  skip_if_not(file.exists(refdb_sh))
  lines <- readLines(refdb_sh, warn = FALSE)
  start <- grep("^rp_taxidmap_validate\\(\\)", lines)
  expect_equal(length(start), 1L)
  end <- start + which(lines[(start + 1):length(lines)] == "}")[1]
  body <- lines[start:end]
  pos_subset <- grep("not in the harvest list", body, fixed = TRUE)
  pos_unique <- grep("duplicate accession row", body, fixed = TRUE)
  expect_equal(length(pos_subset), 1L)
  expect_equal(length(pos_unique), 1L)
  expect_true(pos_subset < pos_unique)
})

test_that("stage_taxidmap's own status messages go through rp_log (stderr), never a bare echo that a caller could capture as map content", {
  skip_if_not(file.exists(refdb_sh))
  s <- readLines(refdb_sh, warn = FALSE)
  start <- grep("^stage_taxidmap\\(\\)", s)
  end <- grep("^stage_fetch\\(\\)", s)
  expect_equal(length(start), 1L)
  expect_true(length(end) >= 1)
  body <- s[start:(end[1] - 1)]
  # the historical bug was a helper writing a human-readable status line to a
  # stdout that was being captured into the taxid map
  suspect <- grep('^[[:space:]]*(echo|printf) ["\x27][A-Za-z]', body, value = TRUE)
  expect_equal(suspect, character(0))
})

# rp_taxidmap_validate -------------------------------------------------------

test_that("taxid map validator's line-shape check rejects a malformed row", {
  r <- run_sourced(c(
    'printf "AB123456.1\n" > "$dir/harvest.txt"',
    'printf "AB123456.1\t9606\nbad row with spaces\t123\n" > "$dir/map.tsv"',
    'rp_taxidmap_validate "$dir/map.tsv" "$dir/harvest.txt"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
  expect_true(any(grepl("not matching accession", r$out, fixed = TRUE)))
})

# THE REQUIRED TEST: reproduces the historical corruption exactly (two
# phantom log lines injected into the map as fake accession rows) and proves
# the guard ORDER is load-bearing -- the subset check catches it, the weaker
# "unique == row count" check alone does not.
test_that("taxid map validator's subset check catches injected log lines that the weaker unique-count check would miss", {
  r <- run_sourced(c(
    'printf "AB123456.1\nCD789012.1\n" > "$dir/harvest.txt"',
    # two real accession rows plus two phantom "log lines" shaped exactly
    # like the historical corruption: syntactically valid map rows
    # (alnum/dot/underscore, tab, digits) that are not accessions at all
    'printf "AB123456.1\t9606\nCD789012.1\t10090\ntaxidmap_scan_nucl_gb_complete\t1\ntaxidmap_scan_nucl_wgs_complete\t1\n" > "$dir/corrupt_map.tsv"',
    # measure the weaker check directly: unique accessions (column 1) == row count
    'total=$(wc -l < "$dir/corrupt_map.tsv" | tr -d " ")',
    'uniq=$(cut -f1 "$dir/corrupt_map.tsv" | sort -u | wc -l | tr -d " ")',
    'if [ "$total" = "$uniq" ]; then echo "WEAK_CHECK_WOULD_PASS total=$total uniq=$uniq"; else echo "WEAK_CHECK_WOULD_REJECT total=$total uniq=$uniq"; fi',
    # now the real, ordered validator
    'rp_taxidmap_validate "$dir/corrupt_map.tsv" "$dir/harvest.txt"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
  # measured result 1: the weak check alone would have passed the corrupt map
  expect_true(any(grepl("WEAK_CHECK_WOULD_PASS total=4 uniq=4", r$out, fixed = TRUE)),
    label = paste(r$out, collapse = "\n"))
  # measured result 2: the real validator rejects it via the subset check,
  # naming the offending phantom line
  expect_true(any(grepl("not in the harvest list", r$out, fixed = TRUE)))
  expect_true(any(grepl("taxidmap_scan_nucl_gb_complete", r$out, fixed = TRUE)))
})

test_that("taxid map validator's uniqueness check catches a duplicated accession row that passes format and subset", {
  r <- run_sourced(c(
    'printf "AB123456.1\nCD789012.1\n" > "$dir/harvest.txt"',
    'printf "AB123456.1\t9606\nAB123456.1\t9606\nCD789012.1\t10090\n" > "$dir/map.tsv"',
    'rp_taxidmap_validate "$dir/map.tsv" "$dir/harvest.txt"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
  expect_true(any(grepl("duplicate accession row", r$out, fixed = TRUE)))
})

test_that("taxid map validator passes a clean, valid map", {
  r <- run_sourced(c(
    'printf "AB123456.1\nCD789012.1\n" > "$dir/harvest.txt"',
    'printf "AB123456.1\t9606\nCD789012.1\t10090\n" > "$dir/map.tsv"',
    'rp_taxidmap_validate "$dir/map.tsv" "$dir/harvest.txt"',
    'echo "REACHED"'
  ))
  expect_equal(r$status, 0L)
  expect_true(any(r$out == "REACHED"))
})

# rp_assert_coverage ---------------------------------------------------------

test_that("rp_assert_coverage enforces the 99.5% floor with integer boundary math, passes trivially on empty input", {
  r_pass_boundary <- run_sourced(c(
    'rp_assert_coverage "test" "199" "200"',
    'echo "REACHED"'
  ))
  expect_equal(r_pass_boundary$status, 0L)
  expect_true(any(r_pass_boundary$out == "REACHED"))

  r_fail_just_below <- run_sourced(c(
    'rp_assert_coverage "test" "198" "200"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r_fail_just_below$status, 0L)
  expect_false(any(r_fail_just_below$out == "UNREACHABLE"))
  expect_true(any(grepl("198 of 200", r_fail_just_below$out, fixed = TRUE)))

  r_empty <- run_sourced(c(
    'rp_assert_coverage "test" "0" "0"',
    'echo "REACHED"'
  ))
  expect_equal(r_empty$status, 0L)
  expect_true(any(r_empty$out == "REACHED"))
})

# rp_taxidmap_scan_dump -------------------------------------------------------

make_gz_fixture <- function(path, lines) {
  con <- gzfile(path, "w")
  writeLines(lines, con)
  close(con)
}

test_that("rp_taxidmap_scan_dump matches wanted accessions from a streamed gz dump, skips the header and unwanted rows, tags the source", {
  fixdir <- withr::local_tempdir()
  dump_path <- file.path(fixdir, "nucl_gb.accession2taxid.gz")
  make_gz_fixture(dump_path, c(
    "accession\taccession.version\ttaxid\tgi",
    "AB111111\tAB111111.1\t100\t1",
    "AB222222\tAB222222.1\t200\t2",
    "ZZ999999\tZZ999999.1\t999\t9"
  ))
  r <- run_sourced(c(
    'printf "AB111111.1\nAB222222.1\n" > "$dir/wanted.txt"',
    sprintf('rp_taxidmap_scan_dump "$dir/wanted.txt" "file://%s" ftp "$dir/out.tsv"', dump_path),
    'cat "$dir/out.tsv"'
  ))
  expect_equal(r$status, 0L)
  rows <- r$out[!grepl("^\\[", r$out)]
  expect_setequal(rows, c("AB111111.1\t100\tftp", "AB222222.1\t200\tftp"))
})

test_that("rp_taxidmap_scan_dump aborts on a truncated gz dump, exit status checked explicitly", {
  fixdir <- withr::local_tempdir()
  good_path <- file.path(fixdir, "good.gz")
  make_gz_fixture(good_path, c(
    "accession\taccession.version\ttaxid\tgi",
    "AB111111\tAB111111.1\t100\t1",
    "AB222222\tAB222222.1\t200\t2"
  ))
  raw <- readBin(good_path, "raw", n = file.info(good_path)$size)
  truncated_path <- file.path(fixdir, "truncated.gz")
  writeBin(raw[seq_len(length(raw) - 5L)], truncated_path)

  r <- run_sourced(c(
    'printf "AB111111.1\nAB222222.1\n" > "$dir/wanted.txt"',
    sprintf('rp_taxidmap_scan_dump "$dir/wanted.txt" "file://%s" ftp "$dir/out.tsv"', truncated_path),
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
  expect_true(any(grepl("scan of", r$out, fixed = TRUE)))
})

test_that("rp_taxidmap_scan_dump aborts when the dump URL does not resolve (curl --fail), not a silent empty match", {
  r <- run_sourced(c(
    'printf "AB111111.1\n" > "$dir/wanted.txt"',
    'rp_taxidmap_scan_dump "$dir/wanted.txt" "file://$dir/does_not_exist.gz" ftp "$dir/out.tsv"',
    'echo "UNREACHABLE"'
  ))
  expect_gt(r$status, 0L)
  expect_false(any(r$out == "UNREACHABLE"))
})

# Full CLI stage -- epost/esummary/xtract stubbed (per the brief: real
# network calls are limited to fetching one dump's HTTP headers, done once,
# separately, during development). The three accession2taxid dumps are real
# gz files served via RP_A2T_BASE_URL=file://..., exercised through the real
# curl | zcat | awk pipeline, not stubbed.
write_taxidmap_stubs <- function(bindir) {
  dir.create(bindir, showWarnings = FALSE, recursive = TRUE)
  write_stub <- function(name, body) {
    p <- file.path(bindir, name)
    writeLines(c("#!/usr/bin/env bash", "set -eu", body), p)
    Sys.chmod(p, "0755")
  }
  write_stub("epost", c(
    'self_dir="$(cd "$(dirname "$0")" && pwd)"',
    'printf x >> "$self_dir/epost_calls.log"',
    'exit 0'
  ))
  write_stub("esummary", "exit 0")
  # Verified against real EDirect (test-edirect, esample, ecommon.sh on this
  # machine; no live epost/esummary call made): esummary is `efetch -docsum`,
  # whose root pattern is DocumentSummary (not DocSum/eSummaryResult), and
  # AccessionVersion/TaxId are real element names on it. Failing on a wrong
  # -pattern/-element here, instead of silently returning nothing, is what
  # would catch this exact class of bug (Task 2's ENTREZ_DIRECT lesson).
  write_stub("xtract", c(
    'pattern="" have_accver=0 have_taxid=0',
    'while [ $# -gt 0 ]; do',
    '  case "$1" in',
    '    -pattern) pattern="$2"; shift 2 ;;',
    '    -element)',
    '      shift',
    '      while [ $# -gt 0 ] && [[ "$1" != -* ]]; do',
    '        [ "$1" = "AccessionVersion" ] && have_accver=1',
    '        [ "$1" = "TaxId" ] && have_taxid=1',
    '        shift',
    '      done',
    '      ;;',
    '    *) shift ;;',
    '  esac',
    'done',
    '[ "$pattern" = "DocumentSummary" ] || { echo "stub xtract: unexpected -pattern $pattern" >&2; exit 1; }',
    '[ "$have_accver" = 1 ] && [ "$have_taxid" = 1 ] || { echo "stub xtract: unexpected -element" >&2; exit 1; }',
    'self_dir="$(cd "$(dirname "$0")" && pwd)"',
    'cat "$self_dir/xtract_output.tsv"'
  ))
  for (noop in c("esearch", "efetch", "update_blastdb.pl", "makeblastdb", "blastdbcmd", "seqkit")) {
    write_stub(noop, "exit 0")
  }
  writeLines(character(0), file.path(bindir, "xtract_output.tsv"))
  file.create(file.path(bindir, "epost_calls.log"))
  bindir
}

run_taxidmap_stage <- function(accessions, dumps = list(), xtract_output = character(0),
                                extra_args = character(0)) {
  build <- withr::local_tempdir(.local_envir = parent.frame())
  bindir <- file.path(build, "stubbin")
  write_taxidmap_stubs(bindir)
  writeLines(xtract_output, file.path(bindir, "xtract_output.tsv"))
  writeLines(accessions, file.path(build, "accessions.txt"))

  fixdir <- file.path(build, "a2t_fixtures")
  dir.create(fixdir)
  hdr <- "accession\taccession.version\ttaxid\tgi"
  all_dumps <- c("nucl_gb", "nucl_wgs", "pdb")
  for (nm in all_dumps) {
    lines <- if (!is.null(dumps[[nm]])) dumps[[nm]] else hdr
    make_gz_fixture(file.path(fixdir, paste0(nm, ".accession2taxid.gz")), lines)
  }

  env <- c(
    paste0("PATH=", bindir, ":", Sys.getenv("PATH")),
    "NCBI_API_KEY=test-key-not-real",
    paste0("RP_A2T_BASE_URL=file://", fixdir)
  )
  out <- suppressWarnings(system2("bash",
    c(shQuote(refdb_sh), "--build-dir", shQuote(build), "--stage", "taxidmap", extra_args),
    stdout = TRUE, stderr = TRUE, env = env))
  st <- attr(out, "status")
  list(status = if (is.null(st)) 0L else as.integer(st), out = out, build = build, bindir = bindir)
}

test_that("taxidmap stage maps accessions across all three ftp dumps plus the esummary fallback, tags each row's source", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  hdr <- "accession\taccession.version\ttaxid\tgi"
  r <- run_taxidmap_stage(
    accessions = c("AB111111.1", "AB222222.1", "CD333333.1", "CD444444.1", "EF555555.1", "GH666666.1"),
    dumps = list(
      nucl_gb = c(hdr, "AB111111\tAB111111.1\t100\t1", "AB222222\tAB222222.1\t200\t2", "ZZ000000\tZZ000000.1\t999\t9"),
      nucl_wgs = c(hdr, "CD333333\tCD333333.1\t300\t3", "CD444444\tCD444444.1\t400\t4"),
      pdb = c(hdr, "EF555555\tEF555555.1\t500\t5")
    ),
    xtract_output = "GH666666.1\t600"
  )
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))

  map <- readLines(file.path(r$build, "taxidmap.tsv"))
  expect_setequal(map, c(
    "AB111111.1\t100", "AB222222.1\t200", "CD333333.1\t300",
    "CD444444.1\t400", "EF555555.1\t500", "GH666666.1\t600"
  ))

  sources <- readLines(file.path(r$build, "taxidmap.sources.tsv"))
  expect_true(any(sources == "GH666666.1\t600\tesummary"))
  expect_true(any(sources == "AB111111.1\t100\tftp"))

  sha_file <- file.path(r$build, "taxidmap.sha256")
  expect_true(file.exists(sha_file))
  want <- system2("sha256sum", shQuote(file.path(r$build, "taxidmap.tsv")), stdout = TRUE)
  want_hash <- sub("\\s.*", "", want)
  expect_true(any(startsWith(readLines(sha_file), want_hash)))

  # residual (one accession) went through epost exactly once, not once per accession
  calls <- readLines(file.path(r$bindir, "epost_calls.log"), warn = FALSE)
  expect_equal(paste(calls, collapse = ""), "x")

  # one dropped file per stage, never a shared dropped.tsv that a later stage
  # would have to append to (see the contract comment in stage_taxidmap)
  expect_true(file.exists(file.path(r$build, "dropped.taxidmap.tsv")))
  expect_equal(length(readLines(file.path(r$build, "dropped.taxidmap.tsv"))), 0L)
  expect_false(file.exists(file.path(r$build, "dropped.tsv")))
})

test_that("taxidmap stage drops distinct unresolved accessions individually into dropped.tsv, never a shared key", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  hdr <- "accession\taccession.version\ttaxid\tgi"
  n <- 1000L
  accs <- sprintf("AB%06d.1", seq_len(n))
  bare <- sprintf("AB%06d", 1:998)
  dump_lines <- c(hdr, sprintf("%s\t%s.1\t%d\t%d", bare, bare, 1:998, 1:998))
  unresolved <- accs[999:1000]

  r <- run_taxidmap_stage(
    accessions = accs,
    dumps = list(nucl_gb = dump_lines),
    xtract_output = character(0)
  )
  expect_equal(r$status, 0L, label = paste(r$out, collapse = "\n"))

  dropped <- read.table(file.path(r$build, "dropped.taxidmap.tsv"), sep = "\t", header = FALSE,
    col.names = c("accession", "reason", "detail"), stringsAsFactors = FALSE, fill = TRUE)
  expect_equal(nrow(dropped), 2L)
  expect_setequal(dropped$accession, unresolved)
  expect_true(all(dropped$reason == "taxid_unresolved"))
  # each dropped record keeps its own accession as the key -- nothing here
  # collapses two different unresolved accessions into one shared row
  expect_equal(length(unique(dropped$accession)), 2L)
})

test_that("taxidmap stage aborts below the 99.5% coverage floor and promotes neither the map nor dropped.tsv", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  hdr <- "accession\taccession.version\ttaxid\tgi"
  r <- run_taxidmap_stage(
    accessions = c("AB111111.1", "AB222222.1", "AB333333.1", "AB444444.1", "AB555555.1"),
    dumps = list(nucl_gb = c(hdr, "AB111111\tAB111111.1\t100\t1", "AB222222\tAB222222.1\t200\t2", "AB333333\tAB333333.1\t300\t3")),
    xtract_output = character(0)
  )
  expect_gt(r$status, 0L)
  expect_false(file.exists(file.path(r$build, "taxidmap.tsv")))
  expect_false(file.exists(file.path(r$build, "dropped.taxidmap.tsv")))
  expect_true(any(grepl("3 of 5 accessions mapped, below the 99.5% floor", r$out, fixed = TRUE)))
})

test_that("taxidmap stage refuses to run before the accessions stage has produced accessions.txt", {
  skip_if_not(file.exists(refdb_sh))
  skip_if(Sys.which("bash") == "")
  build <- withr::local_tempdir()
  bindir <- file.path(build, "stubbin")
  write_taxidmap_stubs(bindir)
  env <- c(paste0("PATH=", bindir, ":", Sys.getenv("PATH")), "NCBI_API_KEY=test-key-not-real")
  out <- suppressWarnings(system2("bash",
    c(shQuote(refdb_sh), "--build-dir", shQuote(build), "--stage", "taxidmap"),
    stdout = TRUE, stderr = TRUE, env = env))
  st <- attr(out, "status")
  expect_gt(if (is.null(st)) 0L else as.integer(st), 0L)
  expect_true(any(grepl("run the accessions stage first", out, fixed = TRUE)))
})
