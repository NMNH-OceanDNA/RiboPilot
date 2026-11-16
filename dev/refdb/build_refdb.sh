#!/usr/bin/env bash
# Orchestrator for the RiboPilot metazoan rDNA BLAST panel build.
#
# Usage: build_refdb.sh [--version V] [--build-dir D] [--stage NAME]
#                        [--jobs N] [--accept-drift] [--force] [--check-version]
#                        [status|clean]
#
# See dev/design/18_REFDB_REBUILD.md for the spec this implements: argument
# parsing, preflight, version sync, and the six build stages (accessions,
# taxidmap, fetch, build, decoys, package).
set -Eeuo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=dev/refdb/lib.sh
source "$SCRIPT_DIR/lib.sh"

# Sourcing must define functions and constants and do nothing else. Tests
# source this file to drive guards directly, and a side effect here would
# rewire the caller's streams, install traps in it, or leave BUILD_DIR
# pointing at the developer's real build tree.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then RP_EXECUTED=1; else RP_EXECUTED=0; fi

RP_LOCK_HELD=""

# ---------------------------------------------------------------------------
# Chunk-worker bookkeeping. Defined here, far from the fetch stage that uses
# it, because rp_cleanup below calls rp_fetch_reap and the EXIT trap that runs
# rp_cleanup is installed a few lines further down: defined any later, every
# exit in between -- --help, an unknown argument -- would run the trap against
# a function that does not exist yet, and command-not-found under errexit would
# abort the trap before it could flush the scrubbed output.
# ---------------------------------------------------------------------------
RP_FETCH_PIDS=""
RP_FETCH_DONE_PID=""

# rp_fetch_reap
# Kill every chunk worker's process group and wait for it. Job control gives
# each worker its own group, so this takes its efetch child down too. It has to
# happen before any exit, for two different reasons: an orphan still holding
# the scrubbed stderr open makes rp_scrub_flush's wait block forever, so the
# script hangs instead of failing; and an orphan that outlives the run can
# publish a chunk into a build tree the next run is already using, which a pid
# lockfile cannot see. The wait is per pid and never bare: a bare wait would
# also wait on the scrub filters, which do not exit until the streams close.
rp_fetch_reap() {
  local p
  # Nothing may interrupt a reap in progress, and job control goes off before
  # the kill rather than after: a signal-killed job makes a shell still in
  # monitor mode print "[1]- Terminated" per worker, interleaved with the real
  # error. Turning monitor off is the whole fix -- redirecting the wait's
  # stderr instead would be far worse, because bash parks the real stderr on a
  # spare descriptor for the duration of a redirected command, and exiting
  # from a trap that fired inside one leaves it parked, so rp_scrub_flush's
  # `exec 2>&-` closes the wrong descriptor and its wait never returns.
  # Process groups are assigned at fork time, so the kill still finds them.
  trap - INT TERM
  set +m
  for p in ${RP_FETCH_PIDS:-}; do kill -- "-$p" 2>/dev/null || true; done
  for p in ${RP_FETCH_PIDS:-}; do wait "$p" || true; done
  RP_FETCH_PIDS=""
}

# rp_fetch_wait
# wait -n, then drop the worker that finished. Without the pruning the list
# would still name every worker of the run -- 179 of them in production,
# nearly all long dead -- and the kill above would then signal whatever
# process group has since been handed a recycled pid.
rp_fetch_wait() {
  local status=0 p keep=""
  wait -n -p RP_FETCH_DONE_PID || status=$?
  for p in ${RP_FETCH_PIDS:-}; do
    [ "$p" = "${RP_FETCH_DONE_PID:-}" ] || keep="$keep $p"
  done
  RP_FETCH_PIDS="$keep"
  return "$status"
}

# rp_fetch_interrupt
# The maintainer runs this stage by hand for over two hours, so Ctrl-C is a
# likely action, and the workers are in their own process groups precisely so
# the terminal's SIGINT does not reach them.
rp_fetch_interrupt() {
  trap - INT TERM
  rp_fetch_reap
  rp_die "fetch: interrupted; stopped the chunk workers"
}

rp_cleanup() {
  # Workers first: rp_scrub_flush waits on the scrub filters, and a live child
  # still holding the scrubbed stderr open never lets that wait return.
  rp_fetch_reap
  [ -z "$RP_LOCK_HELD" ] || rm -f "$RP_LOCK_HELD"
  rp_scrub_flush
}

if [ "$RP_EXECUTED" -eq 1 ]; then
  # Before anything can print: the canonical invocation is
  # `./build_refdb.sh 2>&1 | tee build.log`, and the key must not reach it.
  rp_scrub_output
  trap rp_cleanup EXIT
  trap 'rp_die "unexpected error at ${BASH_SOURCE[0]}:${LINENO}"' ERR
fi

# ---------------------------------------------------------------------------
# Verbatim definitional constants (dev/design/18_REFDB_REBUILD.md section 3).
# Changing any of these is a new panel version, never a silent fix.
# ---------------------------------------------------------------------------
TAXID=33208
RRNA_TERMS='18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal RNA[Title] OR internal transcribed spacer[Title]'
# txid33208 kept literal (not ${TAXID}) so the panel-definition test can match
# the query text directly.
QUERY="(txid33208[Organism:exp]) AND (${RRNA_TERMS}) NOT mitochondrion[Filter] NOT plastid[Filter]"
# Do not add any molecule-type filter here. Measured live: query alone
# 892,123 records; through such a filter, 741,526 -- it kills 100% of RefSeq
# plus 2,232 non-RefSeq records. Off-target records (EST/cDNA "similar to
# rRNA", PDB) are IN the panel by definition; excluding them is a new panel
# version, not a bugfix.
CHUNK_SIZE=5000
DECOY_DBS="16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA ITS_eukaryote_sequences taxdb"
MAKEBLASTDB_FLAGS=(-dbtype nucl -parse_seqids -blastdb_version 5)

# accession->taxid source (dev/design/18_REFDB_REBUILD.md section 1). ~9.4 GB
# total across the three dumps, not the 2.5 GB an old comment once claimed;
# streamed, never landed. Overridable so tests can point it at a file://
# fixture directory instead of the network.
RP_A2T_BASE_URL="${RP_A2T_BASE_URL:-https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid}"
RP_A2T_DUMPS="nucl_gb nucl_wgs pdb"

# Committed count as of 2026-08-15 (dev/design/18_REFDB_REBUILD.md section 4).
# Tolerance covers about 5 weeks of the source's measured drift of ~370
# records/day (891,382 on 2026-08-13 -> 892,123 on 2026-08-15), so a rebuild
# run every few weeks does not require --accept-drift, while a change big
# enough to matter -- a molecule-type filter alone swings the count by
# ~150,000 -- still trips it.
EXPECTED_COUNT=892123
EXPECTED_COUNT_TOLERANCE=15000
# Half the expected count. Far below any plausible real shrink of a database
# that only grows, and far above the collapse cases (0, or a handful of rows
# from an error response) this exists to catch.
EXPECTED_COUNT_FLOOR="${RP_MIN_ACCESSIONS:-446000}"

# rp_panel_key ARTIFACT_HASH
# The key a stage marker records: the panel's definitional identity folded
# together with the artifact hash the stage would otherwise record alone.
#
# A marker keyed on the artifact alone matches across a change of panel
# DEFINITION, which is the thing spec section 7 explicitly plans for ("query
# filters for EST/GSS/PDB ... each is a new panel version"). Bump
# RP_DB_VERSION, edit QUERY, re-run in the same build dir: RP_STATE_FILE is
# $BUILD_DIR/state.txt, not version-scoped, and run_preflight's existing-OUT_DIR
# refusal does not fire because blastdb/<NEW_VER> is a new path. Every stage
# would then skip on a stale marker and blastdb/<NEW_VER> would ship the
# PREVIOUS version's records under a manifest quoting the NEW query, exit 0.
# Folding QUERY and EXPECTED_COUNT in makes that marker unmatchable.
rp_panel_key() {
  printf '%s\n%s\n%s\n' "$QUERY" "$EXPECTED_COUNT" "$1" | sha256sum | cut -d' ' -f1
}

# Fetch stage. Rounds are per chunk and each one requests only what is still
# missing, so a response that dies at 4,900 of 5,000 records costs 100
# accessions on the next round rather than the whole chunk. The ceiling is on
# what is permanently unavailable across the entire build: the WGS masters are
# already filtered by the accessions stage, so the measured expectation is 0
# and this is a ceiling on surprises, not a budget.
FETCH_ROUNDS=3
MAX_UNAVAILABLE="${RP_MAX_UNAVAILABLE:-100}"

# esearch's <Count> and the number of accessions the history actually yields are
# NOT equal, and nothing in E-utilities promises they will be. Measured twice
# against the live panel query: 892,364 reported, 891,825 and 891,876 returned,
# a shortfall of ~0.06% that does not track network failures (the run with
# hundreds of curl errors lost FEWER records than the run with three). The
# historical build that produced the shipped panel never compared them: it took
# its total from `wc -l` of the list it got (01_harvest.sh:53) and guarded only
# against gross truncation. So does this one. The ceiling is on how far the two
# may diverge before the harvest is considered broken rather than merely
# incomplete; 5,000 is ten times the measured shortfall and far below any
# truncation the floor would not already catch.
MAX_UNRETURNED="${RP_MAX_UNRETURNED:-5000}"

REQUIRED_TOOLS="esearch efetch epost esummary xtract update_blastdb.pl makeblastdb blastdbcmd seqkit jq curl sha256sum tar gzip zcat comm"
ALL_STAGES="accessions taxidmap fetch build decoys package"
MIN_FREE_GB=8

# nquire (EDirect) does `command -v curl`, so the shim has to win the PATH
# race before any EDirect call, not just inside run_preflight.
SHIM_DIR="$REPO_ROOT/inst/ref_dbs/getOrganelle/GenBankDownload/bin"
PATH="$SHIM_DIR:$PATH"
export PATH

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
  cat <<'EOF'
Usage: build_refdb.sh [--version V] [--build-dir D] [--stage NAME]
                       [--jobs N] [--accept-drift] [--force] [--check-version]
                       [status|clean]

  --version V       Refuse to run unless R/db_dist.R's RP_DB_VERSION == V.
  --build-dir D     Build tree root (default: dev/refdb/build).
  --stage NAME      Run only one stage (accessions|taxidmap|fetch|build|decoys|package).
  --jobs N          Parallelism for the fetch stage (default: 4).
  --accept-drift    Allow the harvested count to differ from EXPECTED_COUNT.
  --force           Allow rebuilding an existing blastdb/<VER> in place. With
                    --stage NAME it also re-runs that stage even when the state
                    file records it complete.
  --check-version   Diff RP_DB_VERSION across the build's version-carrying files, then exit.
  status            Print the last recorded stage state and exit.
  clean             Remove intermediate build artifacts and exit.
EOF
}

VERSION_ARG=""
# No real-tree default when sourced: a sourced stage function must refuse
# rather than write into the developer's build directory.
if [ "$RP_EXECUTED" -eq 1 ]; then BUILD_DIR="$SCRIPT_DIR/build"; else BUILD_DIR=""; fi
STAGE=""
# Four concurrent fetches take the chunk fetch from 2 h 20 m to roughly 40 min.
# run_preflight already refuses to run without NCBI_API_KEY, so the 10 req/s
# budget this spends against is guaranteed to exist.
JOBS=4
ACCEPT_DRIFT=0
FORCE=0
CHECK_VERSION=0
COMMAND=""

while [ "$RP_EXECUTED" -eq 1 ] && [ $# -gt 0 ]; do
  case "$1" in
    --version)
      VERSION_ARG="${2:?--version requires an argument}"
      shift 2
      ;;
    --build-dir)
      BUILD_DIR="${2:?--build-dir requires an argument}"
      shift 2
      ;;
    --stage)
      STAGE="${2:?--stage requires an argument}"
      shift 2
      ;;
    --jobs)
      JOBS="${2:?--jobs requires an argument}"
      shift 2
      ;;
    --accept-drift)
      ACCEPT_DRIFT=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --check-version)
      CHECK_VERSION=1
      shift
      ;;
    status | clean)
      COMMAND="$1"
      shift
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      rp_die "unknown argument: $1"
      ;;
  esac
done

if [ -n "$STAGE" ]; then
  case " $ALL_STAGES " in
    *" $STAGE "*) : ;;
    *) rp_die "unknown --stage: $STAGE (expected one of: $ALL_STAGES)" ;;
  esac
fi

RP_STATE_FILE="${BUILD_DIR:+$BUILD_DIR/state.txt}"

# ---------------------------------------------------------------------------
# Version helpers (dev/design/18_REFDB_REBUILD.md section 5)
# ---------------------------------------------------------------------------

# rp_read_version
# RP_DB_VERSION in R/db_dist.R is the single source of truth.
rp_read_version() {
  local line ver
  line="$(grep '^RP_DB_VERSION' "$REPO_ROOT/R/db_dist.R" || true)"
  [ -n "$line" ] || rp_die "RP_DB_VERSION not found in R/db_dist.R"
  ver="$(printf '%s\n' "$line" | sed -E 's/^RP_DB_VERSION *<- *"([^"]*)".*/\1/')"
  [ -n "$ver" ] || rp_die "could not parse RP_DB_VERSION from R/db_dist.R: $line"
  printf '%s\n' "$ver"
}

# rp_valid_version VERSION
# Dotted numeric only. numeric_version("v2026.08.1") and "2026.08.1-rc1" are
# both NA, which makes the panel invisible to discover_panel() in R/blast_local.R.
rp_valid_version() {
  [[ "$1" =~ ^[0-9]+(\.[0-9]+)*$ ]]
}

# cmd_check_version
# Diffs RP_DB_VERSION against every place the version string must agree
# (spec section 5). The output-tree, manifest, tarball and release-tag places
# only exist after a build has run; those are checked when present, skipped
# otherwise. The release tag (db-<VER>) is deliberately not checked here: it
# is a network-authoritative fact and --publish is manual by design, so no
# local artifact represents it before a release exists.
cmd_check_version() {
  local ver docker_ver ok=1 d dver manifest_ver f base pat fver sumfile

  ver="$(rp_read_version)"
  rp_valid_version "$ver" || rp_die "RP_DB_VERSION '$ver' is not dotted-numeric"

  docker_ver="$(grep '^ARG RP_DB_VERSION=' "$REPO_ROOT/docker/Dockerfile" || true)"
  docker_ver="${docker_ver#ARG RP_DB_VERSION=}"
  if [ "$docker_ver" != "$ver" ]; then
    rp_log "MISMATCH docker/Dockerfile ARG RP_DB_VERSION=$docker_ver != $ver"
    ok=0
  fi

  # Every version-carrying directory that actually exists on disk, not the one
  # derived from $ver: a derived path can only ever compare $ver to itself, so
  # a stale blastdb/<old>/ next to it was reported OK.
  for d in "$BUILD_DIR"/blastdb/*/; do
    [ -d "$d" ] || continue
    base="$(basename "$d")"
    dver="${base#release-}"
    if [ "$dver" != "$ver" ]; then
      rp_log "MISMATCH $d carries version $dver != $ver"
      ok=0
    fi
    if [ -f "$d/manifest.json" ]; then
      manifest_ver="$(jq -r '.db_version // empty' "$d/manifest.json" 2>/dev/null || true)"
      if [ "$manifest_ver" != "$ver" ]; then
        rp_log "MISMATCH ${d}manifest.json db_version=$manifest_ver != $ver"
        ok=0
      fi
    fi
  done

  for pat in ribopilot_rrna decoys; do
    for f in "$REPO_ROOT"/docker/"${pat}"-*.tar.gz; do
      [ -e "$f" ] || continue
      base="$(basename "$f")"
      fver="${base#"${pat}"-}"
      fver="${fver%.tar.gz}"
      if [ "$fver" != "$ver" ]; then
        rp_log "MISMATCH $base carries version $fver != $ver"
        ok=0
      fi
      sumfile="$f.sha256"
      if [ -f "$sumfile" ] && ! grep -Eq "^[0-9a-f]{64}  ${base}\$" "$sumfile"; then
        rp_log "MISMATCH $sumfile is not in '<hex>  <basename>' form"
        ok=0
      fi
    done
  done

  [ "$ok" -eq 1 ] || rp_die "version check failed: see MISMATCH lines above"
  rp_log "check-version OK: $ver"
}

# ---------------------------------------------------------------------------
# status / clean
# ---------------------------------------------------------------------------

# cmd_status
# Reads RP_STATE_FILE only: one small line per completed stage. Never
# rescans the (up to 1.2 GB) build tree, never globs blastdb/*/manifest.json.
cmd_status() {
  if [ ! -f "$RP_STATE_FILE" ]; then
    echo "no build state at $RP_STATE_FILE (no stage has completed)"
    return 0
  fi
  echo "build state ($RP_STATE_FILE):"
  awk -F'\t' '{printf "  %-12s hash=%s completed=%s\n", $1, $2, $3}' "$RP_STATE_FILE"
}

# cmd_clean
# Spec section 6: drops chunks/, the raw FASTA, TMPDIR's sort spill, and the
# .fasta, .taxid.map and .dropped.tsv -- but only once the release assets verify,
# because until they do these intermediates are the only thing standing between
# a mistake and a 2 h 41 m refetch. Never touches blastdb/<VER> output;
# overwriting that is gated separately by --force.
cmd_clean() {
  local ver panel_dir rel f stage
  rp_need sha256sum
  ver="$(rp_read_version)"
  rp_valid_version "$ver" || rp_die "RP_DB_VERSION '$ver' is not dotted-numeric"
  panel_dir="$BUILD_DIR/blastdb/$ver"
  rel="$BUILD_DIR/blastdb/release-$ver"

  mkdir -p "$BUILD_DIR"
  # The same lock a build takes: unlocked, clean deletes chunks/ out from under
  # a running fetch.
  rp_lock "$BUILD_DIR/build.pid"
  RP_LOCK_HELD="$BUILD_DIR/build.pid"

  # The drop ledger is in this list because the very next thing clean does is
  # delete it from the panel directory, and it is the only record of what the
  # build discarded: it has to have been packaged for release first.
  for f in "ribopilot_rrna-$ver.tar.gz" "decoys-$ver.tar.gz" "dropped-$ver.tsv.gz"; do
    [ -f "$rel/$f" ] \
      || rp_die "refusing to clean: no $rel/$f (run package_release.sh first)"
    [ -f "$rel/$f.sha256" ] \
      || rp_die "refusing to clean: no $rel/$f.sha256"
    (cd "$rel" && sha256sum -c --status "$f.sha256") \
      || rp_die "refusing to clean: $rel/$f does not match $rel/$f.sha256"
  done
  rp_log "release assets verify; cleaning intermediates under $BUILD_DIR"

  # $BUILD_DIR/tmp is run_preflight's TMPDIR default and holds up to 2 GB of
  # sort spill. Only that path, never $TMPDIR: an operator who pointed TMPDIR
  # at a shared scratch directory must not have it deleted from under them.
  rm -rf "$BUILD_DIR/chunks" "$BUILD_DIR/chunks.missing" "$BUILD_DIR/chunks.tmp" \
    "$BUILD_DIR/taxidmap.tmp" "$BUILD_DIR/build.tmp" "$BUILD_DIR/tmp"
  rm -f "$BUILD_DIR"/*.raw.fasta
  # The resume key goes with the chunks it keys: left behind, it would tell the
  # next fetch that a chunk dir which no longer exists was proven.
  rm -f "$BUILD_DIR/chunks.accessions.sha256" "$BUILD_DIR/dropped.fetch.tsv"
  # Harvest-stage leftovers. accessions.txt and accessions.sha256 stay: they
  # are the fingerprint the fetch stage refuses to resume against.
  rm -f "$BUILD_DIR/accessions.raw.txt" "$BUILD_DIR/accessions.esearch.xml" \
    "$BUILD_DIR/taxidmap.tsv" "$BUILD_DIR/taxidmap.sources.tsv" \
    "$BUILD_DIR/taxidmap.sha256" "$BUILD_DIR/dropped.taxidmap.tsv" \
    "$BUILD_DIR/dropped.accessions.tsv"
  rm -f "$panel_dir"/*.fasta "$panel_dir"/*.taxid.map "$panel_dir"/*dropped.tsv

  # Only the markers whose outputs this removed. Deleting the whole state file
  # would re-run stages whose artifacts are still there, fetch included.
  for stage in taxidmap fetch; do
    rp_marker_clear "$stage"
  done
  rp_log "clean complete (blastdb/ output and release tarballs left in place)"
}

# ---------------------------------------------------------------------------
# Preflight (dev/design/18_REFDB_REBUILD.md section 4, guard list)
# ---------------------------------------------------------------------------

# rp_require_build_dir
# A stage reached by sourcing has no BUILD_DIR, and must refuse rather than
# fall back to a default that would write into the developer's real tree.
rp_require_build_dir() {
  [ -n "$BUILD_DIR" ] || rp_die "BUILD_DIR is not set; run this script rather than sourcing it, or pass --build-dir"
}

rp_require_free_space() {
  local dir="$1" min_gb="$2" avail_kb
  mkdir -p "$dir"
  avail_kb="$(df -Pk "$dir" | awk 'NR==2 {print $4}')"
  [ -n "$avail_kb" ] || rp_die "could not determine free space for $dir"
  if [ "$avail_kb" -lt $((min_gb * 1024 * 1024)) ]; then
    rp_die "need at least ${min_gb}GB free on the filesystem holding $dir (have $((avail_kb / 1024 / 1024))GB)"
  fi
}

run_preflight() {
  local tool resolved_curl ver

  # rp_fetch_wait uses `wait -n -p`, which is bash 5.1. On 4.4 the first
  # saturated reap prints "wait: -p: invalid option" and the run aborts with
  # the misleading "fetch: chunk worker(s) failed" -- two hours in.
  if [ "${BASH_VERSINFO[0]}" -lt 5 ] \
    || { [ "${BASH_VERSINFO[0]}" -eq 5 ] && [ "${BASH_VERSINFO[1]}" -lt 1 ]; }; then
    rp_die "bash 5.1 or newer is required (this is ${BASH_VERSION}); the fetch stage uses 'wait -n -p'"
  fi

  for tool in $REQUIRED_TOOLS; do
    rp_need "$tool"
  done

  resolved_curl="$(command -v curl || true)"
  [ "$resolved_curl" = "$SHIM_DIR/curl" ] \
    || rp_die "curl resolves to '$resolved_curl', not the repo shim at $SHIM_DIR/curl; check PATH"

  rp_require_free_space "$BUILD_DIR" "$MIN_FREE_GB"

  # sort -S / TMPDIR must never be hardcoded (e.g. -S 4G): wrong on a shared node.
  RP_SORT_MEM="${RP_SORT_MEM:-2G}"
  TMPDIR="${TMPDIR:-$BUILD_DIR/tmp}"
  mkdir -p "$TMPDIR"
  export RP_SORT_MEM TMPDIR

  [ -n "${NCBI_API_KEY:-}" ] \
    || rp_die "NCBI_API_KEY is not set; without it EDirect rate-limits to 3 req/s instead of 10 and the fetch stage takes hours longer"

  ver="$(rp_read_version)"
  rp_valid_version "$ver" || rp_die "RP_DB_VERSION '$ver' is not dotted-numeric"
  if [ -n "$VERSION_ARG" ] && [ "$VERSION_ARG" != "$ver" ]; then
    rp_die "--version $VERSION_ARG does not match RP_DB_VERSION $ver in R/db_dist.R"
  fi
  VER="$ver"

  OUT_DIR="$BUILD_DIR/blastdb/$VER"
  # The decoys and package stages read the built panel and write beside it;
  # neither rebuilds it, and both are idempotent by spec. Applying the overwrite
  # refusal to them would make their normal case -- re-fetch or re-package a
  # panel that is already built -- reachable only by passing a flag that means
  # "rebuild the panel in place", which is not what the operator wants to say.
  case "$STAGE" in
    decoys | package) ;;
    *)
      if [ -d "$OUT_DIR" ] && [ "$FORCE" -ne 1 ]; then
        rp_die "$OUT_DIR already exists; pass --force to rebuild it in place, or finish the run without rebuilding it by running the remaining stages one at a time (--stage decoys, then --stage package), which do not need --force"
      fi
      ;;
  esac
}

# ---------------------------------------------------------------------------
# Accessions stage guards (dev/design/18_REFDB_REBUILD.md section 4).
# ---------------------------------------------------------------------------

# rp_line_count FILE
# wc -l, never a grep-match-count with an `|| echo 0` fallback: on an
# existing EMPTY file, grep's own count and the fallback both print, so that
# idiom emits two lines of "0"; `[ "$n" -eq ... ]` then sees a two-line
# argument and errors with rc 2, which an unguarded caller treats as the
# guard being skipped, promoting an empty file as success. wc -l has none of
# that: one line of output, "0" on an empty file.
rp_line_count() {
  wc -l < "$1" | tr -d ' '
}

# rp_assert_count LABEL ACTUAL EXPECTED
# Exact-match guard: replaces the historical `-lt 500000` floor, which sat
# 44% below the true count and passed a stream that died at 600k.
rp_assert_count() {
  local label="$1" actual="$2" expected="$3"
  [ "$actual" = "$expected" ] \
    || rp_die "$label: expected $expected, got $actual"
}

# rp_assert_unique LABEL FILE
# Count == lines == unique was verified on three clades; a duplicate means
# the harvest looped or the history-server replay landed the wrong page.
rp_assert_unique() {
  local label="$1" file="$2" total uniq
  total="$(rp_line_count "$file")"
  uniq="$(sort -u "$file" | wc -l | tr -d ' ')"
  [ "$total" = "$uniq" ] \
    || rp_die "$label: $file has $((total - uniq)) duplicate line(s) among $total"
}

# rp_filter_wgs_masters IN OUT DROPPED
# WGS master accessions (4 letters + 6 zeros) are permanently unavailable via
# efetch by design, so dropping them here is a recorded expected drop, not a
# silent one. Recorded in DROPPED as well as logged: the build stage's ledger
# identity is computed from files, so a drop that exists only in the log makes
# it uncomputable.
rp_filter_wgs_masters() {
  local in="$1" out="$2" dropped="$3" total kept
  grep -Ev '^[A-Z]{4}0{6}' "$in" > "$out" || true
  grep -E '^[A-Z]{4}0{6}' "$in" \
    | awk -v OFS='\t' '{print $1, "wgs_master_unavailable", ""}' > "$dropped" || true
  total="$(rp_line_count "$in")"
  kept="$(rp_line_count "$out")"
  rp_log "accessions: dropped $((total - kept)) WGS master accession(s) (expected, permanently unavailable)"
}

# rp_assert_shape LABEL FILE
# EDirect reports total failure as a line of prose with exit status 0, so
# "QUERY FAILURE" otherwise passes counting, uniqueness and WGS filtering and
# lands in the accession list as though it were an accession.
rp_assert_shape() {
  local label="$1" file="$2" bad
  bad="$(grep -vE '^[A-Za-z0-9._]+$' "$file" | head -1 || true)"
  [ -z "$bad" ] || rp_die "$label: not an accession: $bad"
}

# rp_check_drift ACTUAL ACCEPT_DRIFT
# Beyond EXPECTED_COUNT_TOLERANCE, refuse unless ACCEPT_DRIFT=1, which logs
# the new count instead of aborting. Within tolerance, silent: the source
# moves ~370 records/day and that is not news.
# rp_check_collapse ACTUAL
# A floor --accept-drift cannot waive. Drift is gradual (~370 records/day), so
# --accept-drift becomes routine within weeks and would otherwise also wave
# through a transiently rejected query returning Count=0: every later guard
# then passes honestly on an empty set and an empty panel ships.
rp_check_collapse() {
  local actual="$1"
  [ "$actual" -lt "$EXPECTED_COUNT_FLOOR" ] || return 0
  rp_die "accessions: count $actual is below the collapse floor of $EXPECTED_COUNT_FLOOR; --accept-drift does not waive this. The query or the source is broken, not merely drifting."
}

rp_check_drift() {
  local actual="$1" accept_drift="$2" delta
  if [ "$actual" -ge "$EXPECTED_COUNT" ]; then
    delta=$((actual - EXPECTED_COUNT))
  else
    delta=$((EXPECTED_COUNT - actual))
  fi
  [ "$delta" -gt "$EXPECTED_COUNT_TOLERANCE" ] || return 0
  if [ "$accept_drift" -eq 1 ]; then
    rp_log "accessions: count drifted by $delta (EXPECTED_COUNT=$EXPECTED_COUNT, got $actual); accepted via --accept-drift"
  else
    rp_die "accessions: count $actual is $delta from EXPECTED_COUNT=$EXPECTED_COUNT (tolerance $EXPECTED_COUNT_TOLERANCE); pass --accept-drift to proceed and record the new count"
  fi
}

# ---------------------------------------------------------------------------
# Taxidmap stage guards (dev/design/18_REFDB_REBUILD.md section 4, Taxonomy).
# ---------------------------------------------------------------------------

# RP_TAXIDMAP_LINE_RE
# accession<TAB>taxid, nothing else -- the exact shape makeblastdb -taxid_map
# requires. Built with an actual embedded tab, not the two characters "\t":
# POSIX extended regex does not guarantee that escape means tab.
RP_TAXIDMAP_LINE_RE='^[A-Za-z0-9._]+'$'\t''[0-9]+$'

# rp_taxidmap_scan_dump WANTED_FILE URL SOURCE_TAG OUT_FILE
# Streams URL, never landed: curl --fail | zcat | awk. Matches column 2
# (accession.version) of the dump against WANTED_FILE (one accession per
# line) and writes accession<TAB>taxid<TAB>SOURCE_TAG per match to OUT_FILE.
#
# zcat's stderr is left connected here, never routed to the bit bucket:
# gzip's own CRC catches truncation and corruption for free, so a
# stderr-visible zcat failure plus `pipefail` (lib.sh) is enough to abort.
# `curl --fail` catches the other failure mode, an HTTP error page fed to
# zcat as if it were gzip data. The exit status is still checked explicitly
# (the `if !` below), not left to implicit propagation, so the abort message
# can name the URL: a 40% truncation of one of these dumps once silently
# produced 356,851 taxid-0 records.
rp_taxidmap_scan_dump() {
  local wanted="$1" url="$2" source_tag="$3" out="$4"
  if ! curl --fail -sL "$url" \
      | zcat \
      | awk -F'\t' -v OFS='\t' -v src="$source_tag" '
          NR==FNR { want[$1] = 1; next }
          $1 == "accession" && $2 == "accession.version" { next }
          ($2 in want) { print $2, $3, src; delete want[$2] }
        ' "$wanted" - \
      > "$out"; then
    rp_die "taxidmap: scan of $url failed (curl/zcat/awk pipeline exit nonzero)"
  fi
}

# rp_taxidmap_validate MAP HARVEST
# Three checks, in this exact order -- the order is load-bearing, not
# stylistic:
#   1. every line matches accession<TAB>taxid.
#   2. the map's accession set is a SUBSET of the harvest list.
#   3. unique accessions in the map == row count.
# Step 2 is the one that catches the historical corruption: a helper wrote
# two progress lines to a stdout that was being captured as the map, and
# each injected line -- e.g. "taxidmap_scan_nucl_gb_complete<TAB>1" -- is a
# syntactically valid, individually unique row. Step 3 run alone PASSES a map
# corrupted this way (it only compares the map to itself); step 2 fails it,
# because neither injected line is a harvested accession. Hence step 2 must
# run before step 3, and it must run at all.
rp_taxidmap_validate() {
  local map="$1" harvest="$2" bad total uniq extra map_ids harvest_ids

  bad="$(grep -Ev "$RP_TAXIDMAP_LINE_RE" "$map" || true)"
  [ -z "$bad" ] \
    || rp_die "taxidmap: map has line(s) not matching accession<TAB>taxid: $(printf '%s\n' "$bad" | head -1)"

  map_ids="$(mktemp)"
  harvest_ids="$(mktemp)"
  cut -f1 "$map" | sort -u > "$map_ids"
  sort -u "$harvest" > "$harvest_ids"
  extra="$(comm -23 "$map_ids" "$harvest_ids" || true)"
  rm -f "$map_ids" "$harvest_ids"
  [ -z "$extra" ] \
    || rp_die "taxidmap: map has accession(s) not in the harvest list: $(printf '%s\n' "$extra" | head -1)"

  total="$(rp_line_count "$map")"
  uniq="$(cut -f1 "$map" | sort -u | wc -l | tr -d ' ')"
  [ "$total" = "$uniq" ] \
    || rp_die "taxidmap: map has $((total - uniq)) duplicate accession row(s) among $total"
}

# rp_assert_coverage LABEL MAPPED REQUESTED
# Mapped rows must be >= 99.5% of requested accessions (historical coverage
# was 100%). Integer comparison (mapped*1000 vs requested*995) avoids a
# floating-point percentage; on REQUESTED=0 this is 0 >= 0, which passes
# rather than dividing by zero.
rp_assert_coverage() {
  local label="$1" mapped="$2" requested="$3"
  [ "$((mapped * 1000))" -ge "$((requested * 995))" ] \
    || rp_die "$label: $mapped of $requested accessions mapped, below the 99.5% floor"
}

# ---------------------------------------------------------------------------
# Fetch stage helpers (dev/design/18_REFDB_REBUILD.md section 4, Fetch).
# ---------------------------------------------------------------------------

# rp_fasta_complete FILE
# Pass through only whole records: a header, at least one non-empty sequence
# line, and -- for the final record -- a stream that actually ended with a
# newline. A response cut after its last header is the case every count-based
# check waves through, and appending the next round onto its unterminated line
# would splice two records into one.
rp_fasta_complete() {
  local f="$1" ended=0
  if [ -s "$f" ] && [ "$(tail -c 1 "$f" | wc -l)" -eq 1 ]; then ended=1; fi
  awk -v ended="$ended" '
    /^>/ { if (n > 0) printf "%s", rec; rec = $0 "\n"; n = 0; next }
    rec != "" { rec = rec $0 "\n"; if ($0 != "") n++ }
    END { if (n > 0 && ended == 1) printf "%s", rec }
  ' "$f"
}

# rp_fasta_accessions FILE...
# The accession set a fasta actually carries, sorted and unique. efetch returns
# a PDB record as ">pdb|28ZU|S2" for an accession the list requests as
# "28ZU_S2", so without this normalisation the panel's 298 PDB records read as
# 298 missing accessions and 298 unrequested ones at the same time.
#
# This is the ONLY place that normalisation is spelled. The build stage feeds
# its PDB deflines back through this function rather than repeating the rule in
# awk: two spellings that drift would put an accession in the drop ledger that
# the harvest never held, at the right row count, which is exactly the shape of
# error the ledger identity was blind to before it compared sets.
rp_fasta_accessions() {
  sed -n -E 's/^>([^ ]+).*/\1/p' "$@" \
    | sed -E 's/^pdb\|([^|]+)\|(.+)$/\1_\2/' \
    | sort -u
}

# rp_fetch_chunk IDX FIRST_LINE
# One chunk end to end, run as a background worker. Cuts its own request list
# from the sorted, hash-pinned accession list by line range (positional cutting
# is safe only because of that pinning), then fetches and retries ONLY what is
# still missing. An existing chunk file seeds the accumulator, so resuming a
# chunk that was interrupted before its missing list landed re-requests a
# handful of accessions rather than all 5000.
rp_fetch_chunk() {
  local idx="$1" first="$2"
  local tmp="$BUILD_DIR/chunks.tmp"
  local want="$tmp/want_$idx" need="$tmp/need_$idx" got="$tmp/got_$idx"
  local extra="$tmp/extra_$idx" resp="$tmp/resp_$idx" part="$tmp/fa_$idx.tmp"
  local last=$((first + CHUNK_SIZE - 1)) round n_got n_missing

  sed -n "${first},${last}p;$((last + 1))q" "$BUILD_DIR/accessions.txt" \
    | sort -u > "$want"
  if [ -f "$BUILD_DIR/chunks/fa_$idx" ]; then
    rp_fasta_complete "$BUILD_DIR/chunks/fa_$idx" > "$part"
  else
    : > "$part"
  fi
  rp_fasta_accessions "$part" > "$got"
  comm -23 "$want" "$got" > "$need"

  for ((round = 1; round <= FETCH_ROUNDS; round++)); do
    [ -s "$need" ] || break
    # Never branch on efetch's exit status. It exits 0 with zero bytes and a
    # QUERY FAILURE line on total failure, which is what made the historical
    # repair script's `&& break` retry loop inert. Only the content decides.
    efetch -db nuccore -input "$need" -format fasta > "$resp" || true
    rp_fasta_complete "$resp" >> "$part"
    rp_fasta_accessions "$part" > "$got"
    comm -13 "$want" "$got" > "$extra"
    [ ! -s "$extra" ] \
      || rp_die "fetch: chunk $idx returned $(rp_line_count "$extra") record(s) that were never requested, e.g. $(head -1 "$extra")"
    comm -23 "$want" "$got" > "$need"
  done

  n_got="$(rp_line_count "$got")"
  n_missing="$(rp_line_count "$need")"
  [ "$n_got" -gt 0 ] \
    || rp_die "fetch: chunk $idx returned no usable records in $FETCH_ROUNDS round(s); an empty body with exit status 0 is exactly how efetch reports total failure"

  # Both files or neither: the missing list lands first, then the chunk, each
  # by atomic rename. A worker killed between the two leaves a chunk whose
  # completeness resume must not assume, and re-verifying one is cheap.
  mv "$need" "$BUILD_DIR/chunks.missing/fa_$idx.txt"
  mv "$part" "$BUILD_DIR/chunks/fa_$idx"
  rm -f "$want" "$got" "$extra" "$resp"
  rp_log "fetch: chunk $idx done: $n_got record(s), $n_missing unavailable"
}

# ---------------------------------------------------------------------------
# Build stage helpers (dev/design/18_REFDB_REBUILD.md section 4, Build/package).
# ---------------------------------------------------------------------------

# RP_PANEL_FILES
# Exactly what package_release.sh:21-31 puts in the panel tarball.
# unpacked_bytes is measured over this list and nothing else: `du -sb` over the
# panel directory reported 1.14 GB for a 331 MiB database, because it also
# counted the .fasta, .taxid.map and .dropped.tsv, which are build artifacts and
# are never shipped.
RP_PANEL_FILES="manifest.json ribopilot_rrna.ndb ribopilot_rrna.nhr ribopilot_rrna.nin ribopilot_rrna.njs ribopilot_rrna.nog ribopilot_rrna.nos ribopilot_rrna.not ribopilot_rrna.nsq ribopilot_rrna.ntf ribopilot_rrna.nto"

# RP_DROP_PRECEDENCE
# Which reason wins when an accession is recorded by more than one stage. The
# per-stage ledgers do not partition and are positively correlated: a record
# suppressed at NCBI between the harvest and the fetch fails both the
# accession2taxid dump and efetch's live lookup, so it lands in
# dropped.taxidmap.tsv and dropped.fetch.tsv alike. The order is how far the
# record got, earliest failure winning, so the reason names where it stopped.
RP_DROP_PRECEDENCE="wgs_master_unavailable efetch_unavailable taxid_unresolved curated_exclusion pdb_structure_record_nonstandard_defline byte_identical_to_kept_record_same_taxid"

# Curated exclusions: GenBank records known to be wrong (contaminated genome
# annotations, misidentified vouchers). One row per accession, first column the
# versioned accession; the rest is the evidence trail. Applied after the dedup:
# a byte-identical same-taxid duplicate of an excluded record is the same
# sequence from the same genome, so it goes too, naming the listed record.
RP_EXCLUSIONS="${RP_EXCLUSIONS:-$SCRIPT_DIR/exclusions.csv}"

# rp_shipped_bytes DIR
# Total size of the shipped files present in DIR.
rp_shipped_bytes() {
  local dir="$1" f total=0
  for f in $RP_PANEL_FILES; do
    [ -f "$dir/$f" ] || continue
    total=$((total + $(stat -c %s "$dir/$f")))
  done
  printf '%s\n' "$total"
}

# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

# stage_accessions
# One esearch against nuccore with the verbatim QUERY. The XML is captured
# once and replayed into efetch on STDIN (not a second esearch), so the
# accession list and the <Count> it is checked against come from the same
# history-server result set. Stdin, because that is how EDirect chains: -input
# reads IDENTIFIERS from a file, so handing it this envelope leaves efetch with
# no database at all.
stage_accessions() {
  local xml raw filtered count harvested unreturned
  rp_require_build_dir

  mkdir -p "$BUILD_DIR"
  xml="$BUILD_DIR/accessions.esearch.xml"
  raw="$BUILD_DIR/accessions.raw.txt"
  filtered="$BUILD_DIR/accessions.txt"

  # Marker + input-hash gate, spec section 6: "re-run the identical command;
  # completed stages skip". Without it a plain restart re-harvests, and the
  # source moves ~370 records/day, so any drift inside the tolerance rewrites
  # accessions.txt with a new hash -- the very key stage_fetch refuses to
  # resume against, destroying the 2 h 20 m fetch the restart existed to
  # resume. --accept-drift does NOT re-harvest: it relaxes a guard, it does not
  # mean "redo", and making it re-harvest would put the destruction straight
  # back.
  #
  # --force re-harvests only when the operator NAMED this stage. A bare --force
  # on a full run means "rebuild blastdb/<VER> in place" -- run_preflight
  # demands it once that directory exists -- and re-harvesting on the strength
  # of that would put the operator back on the destructive path for every
  # restart after the build stage completed, which is the only route the gate
  # above does not already cover. `--stage accessions --force` is the precise
  # way to say re-harvest, and is the one thing that does it.
  #
  # The marker key is rp_panel_key, not the bare list hash: see there. A gate
  # keyed on the artifact alone skips across a QUERY change and ships the
  # previous panel's records under the new version.
  if { [ "$FORCE" -eq 0 ] || [ "$STAGE" != "accessions" ]; } \
      && [ -f "$filtered" ] && [ -f "$BUILD_DIR/accessions.sha256" ] \
      && (cd "$BUILD_DIR" && sha256sum -c --status accessions.sha256) \
      && rp_marker_ok accessions "$(rp_panel_key "$(cut -d' ' -f1 < "$BUILD_DIR/accessions.sha256")")"; then
    rp_log "accessions: already harvested for this panel definition ($(rp_line_count "$filtered") accession(s)); skipping (--stage accessions --force re-harvests)"
    return 0
  fi

  rp_log "accessions: esearch (nuccore) for the panel query"
  esearch -db nuccore -query "$QUERY" > "$xml"

  # esearch's default output is its own ENTREZ_DIRECT chaining wrapper
  # (WebEnv/QueryKey/Count), not raw eSearchResult XML -- verified against a
  # live esearch, not assumed.
  count="$(xtract -input "$xml" -pattern ENTREZ_DIRECT -element Count)"
  [ -n "$count" ] && [[ "$count" =~ ^[0-9]+$ ]] \
    || rp_die "accessions: could not parse a numeric <Count> from esearch output"

  rp_check_collapse "$count"
  rp_check_drift "$count" "$ACCEPT_DRIFT"

  rp_log "accessions: replaying the same WebEnv/QueryKey into efetch -format acc ($count expected)"
  efetch -format acc < "$xml" > "$raw"

  rp_assert_shape "accessions: harvested list" "$raw"
  rp_assert_unique "accessions: harvested accession list" "$raw"

  # The list that came back is the requested set, not <Count>. Guarded the way
  # the historical build guarded it -- a floor against truncation -- plus a
  # bound on how far the two may disagree, because a LARGE divergence means the
  # replay died partway rather than that the source is merely inexact.
  harvested="$(rp_line_count "$raw")"
  rp_check_collapse "$harvested"
  [ "$harvested" -le "$count" ] \
    || rp_die "accessions: the history returned $harvested accession(s) but esearch reported only $count; the replayed result set is not the one that was searched"
  unreturned=$((count - harvested))
  [ "$unreturned" -le "$MAX_UNRETURNED" ] \
    || rp_die "accessions: esearch reported $count record(s) but the history returned only $harvested, $unreturned short, above the ceiling of $MAX_UNRETURNED; the replay was truncated"
  [ "$unreturned" -eq 0 ] \
    || rp_log "accessions: esearch reported $count, the history returned $harvested ($unreturned not returned, within the ceiling of $MAX_UNRETURNED)"

  rp_filter_wgs_masters "$raw" "$BUILD_DIR/accessions.wgsfiltered.txt" \
    "$BUILD_DIR/dropped.accessions.tsv"
  # Sorted before hashing: efetch history-server order is not contractually
  # stable, and the fetch stage refuses to resume on a fingerprint mismatch,
  # so an unsorted hash would discard a valid multi-hour fetch.
  sort "$BUILD_DIR/accessions.wgsfiltered.txt" > "$filtered"
  rm -f "$BUILD_DIR/accessions.wgsfiltered.txt"

  (cd "$BUILD_DIR" && sha256sum "$(basename "$filtered")" > accessions.sha256)
  # Three distinctly named counts, because conflating any two of them lets a
  # later stage reconcile against records it never held. esearch_count is what
  # the SOURCE claimed and drives the drift guards only; requested_count is what
  # was actually harvested and is the ledger identity's `requested`;
  # final_count is post-WGS-filter and is what the chunking divides.
  printf '%s\n' "$count" > "$BUILD_DIR/accessions.esearch_count"
  printf '%s\n' "$harvested" > "$BUILD_DIR/accessions.requested_count"
  rp_line_count "$filtered" > "$BUILD_DIR/accessions.final_count"
  rp_marker_set accessions "$(rp_panel_key "$(cut -d' ' -f1 < "$BUILD_DIR/accessions.sha256")")"
  rp_log "accessions: wrote $(rp_line_count "$filtered") accession(s) to $filtered"
}

# stage_taxidmap
# Depends only on accessions.txt, so it is safe to run concurrently with
# fetch (spec: both depend on the same accession list, not on each other).
# Three sources, tried in order, each accession claimed by at most one:
#   1-3. the three accession2taxid dumps (rp_taxidmap_scan_dump), each dump
#        searched only for accessions the prior sources did not resolve.
#   4. epost | esummary | xtract for whatever is left (newer than the dumps).
# Everything is assembled in a scratch dir and validated before anything is
# written under $BUILD_DIR, so a failed run (bad coverage, a corrupt map)
# never promotes a partial taxidmap.tsv or dropped.tsv -- same reasoning as
# stage_accessions never promoting a truncated accessions.txt.
stage_taxidmap() {
  local harvest requested tmp combined dump url out residual candidate_map
  rp_require_build_dir
  local candidate_sources mapped unresolved_n dropped_tmp ftp_n esummary_n

  mkdir -p "$BUILD_DIR"
  harvest="$BUILD_DIR/accessions.txt"
  [ -f "$harvest" ] || rp_die "taxidmap: $harvest not found; run the accessions stage first"
  requested="$(rp_line_count "$harvest")"

  tmp="$BUILD_DIR/taxidmap.tmp"
  rm -rf "$tmp"
  mkdir -p "$tmp"

  sort -u "$harvest" > "$tmp/still_wanted.txt"
  combined="$tmp/combined.tsv"
  : > "$combined"

  for dump in $RP_A2T_DUMPS; do
    url="$RP_A2T_BASE_URL/$dump.accession2taxid.gz"
    out="$tmp/matched.$dump.tsv"
    rp_log "taxidmap: scanning $dump dump ($(rp_line_count "$tmp/still_wanted.txt") accession(s) still wanted)"
    rp_taxidmap_scan_dump "$tmp/still_wanted.txt" "$url" ftp "$out"
    cat "$out" >> "$combined"
    cut -f1 "$out" | sort -u > "$tmp/found.$dump.txt"
    comm -23 "$tmp/still_wanted.txt" "$tmp/found.$dump.txt" > "$tmp/still_wanted.next.txt"
    mv "$tmp/still_wanted.next.txt" "$tmp/still_wanted.txt"
  done

  residual="$tmp/still_wanted.txt"
  if [ "$(rp_line_count "$residual")" -gt 0 ]; then
    rp_log "taxidmap: $(rp_line_count "$residual") accession(s) not in the ftp dumps; falling back to epost | esummary"
    if ! epost -db nuccore -format acc -input "$residual" \
        | esummary \
        | xtract -pattern DocumentSummary -element AccessionVersion TaxId \
        > "$tmp/esummary.raw.tsv"; then
      rp_die "taxidmap: epost | esummary | xtract pipeline failed for $(rp_line_count "$residual") residual accession(s)"
    fi
    # xtract prints an empty second column for a DocumentSummary with no
    # TaxId (withdrawn/obsolete record); keep only rows that actually
    # resolved. Never `grep -c ... || echo 0` here either -- rp_line_count
    # below is the only count taken from this.
    grep -E $'^[^\t]+\t[0-9]+$' "$tmp/esummary.raw.tsv" > "$tmp/esummary.matched.tsv" || true
    awk -F'\t' -v OFS='\t' '{print $1, $2, "esummary"}' "$tmp/esummary.matched.tsv" >> "$combined"
    cut -f1 "$tmp/esummary.matched.tsv" | sort -u > "$tmp/found.esummary.txt"
    comm -23 "$residual" "$tmp/found.esummary.txt" > "$tmp/still_wanted.next.txt"
    mv "$tmp/still_wanted.next.txt" "$residual"
  fi

  candidate_map="$tmp/taxidmap.tsv"
  candidate_sources="$tmp/taxidmap.sources.tsv"
  cut -f1,2 "$combined" > "$candidate_map"
  cp "$combined" "$candidate_sources"

  rp_taxidmap_validate "$candidate_map" "$harvest"

  mapped="$(rp_line_count "$candidate_map")"
  rp_assert_coverage "taxidmap: mapped rows vs requested accessions" "$mapped" "$requested"

  # Only now, with both guards passed, commit. Each stage owns its own
  # dropped-record file and rewrites it wholesale; the build stage merges them
  # into the shipped dropped.tsv. A single shared file would make this stage's
  # idempotent rewrite silently discard another stage's rows.
  unresolved_n="$(rp_line_count "$residual")"
  dropped_tmp="$tmp/dropped.tsv"
  if [ "$unresolved_n" -gt 0 ]; then
    # One row per accession, never a shared "NA" placeholder: a shared key
    # would let unrelated species that happen to share byte-identical rRNA
    # collide in the build stage's (taxid, sequence) dedup and be mislabeled
    # as duplicates of each other.
    awk -v OFS='\t' '{print $1, "taxid_unresolved", ""}' "$residual" > "$dropped_tmp"
  else
    : > "$dropped_tmp"
  fi
  cp "$dropped_tmp" "$BUILD_DIR/dropped.taxidmap.tsv"

  cp "$candidate_map" "$BUILD_DIR/taxidmap.tsv"
  cp "$candidate_sources" "$BUILD_DIR/taxidmap.sources.tsv"
  (cd "$BUILD_DIR" && sha256sum taxidmap.tsv > taxidmap.sha256)
  rp_marker_set taxidmap "$(cut -d' ' -f1 < "$BUILD_DIR/taxidmap.sha256")"

  ftp_n="$(awk -F'\t' '$3=="ftp"' "$combined" | wc -l | tr -d ' ')"
  esummary_n="$(awk -F'\t' '$3=="esummary"' "$combined" | wc -l | tr -d ' ')"
  rp_log "taxidmap: wrote $mapped mapped accession(s) ($ftp_n ftp, $esummary_n esummary), $unresolved_n unresolved (dropped) of $requested requested"

  rm -rf "$tmp"
}
# stage_fetch
# The panel FASTA, 5000 accessions per chunk, JOBS chunks at a time. Resumable
# at chunk granularity and keyed to the sha256 of the accession list, because
# re-cutting a list that moved would positionally keep withdrawn records and
# miss new ones without a single chunk name changing.
stage_fetch() {
  rp_require_build_dir
  local acc sidecar hex recorded total n_chunks chunks missing tmp key f
  local i idx running failed
  local allwant allgot notgot notwant n_unavailable n_got

  acc="$BUILD_DIR/accessions.txt"
  sidecar="$BUILD_DIR/accessions.sha256"
  [ -f "$acc" ] || rp_die "fetch: $acc not found; run the accessions stage first"
  [ -f "$sidecar" ] || rp_die "fetch: $sidecar not found; run the accessions stage first"
  # Never zero chunks and a green marker: every later guard passes honestly on
  # an empty set, and an empty panel ships.
  total="$(rp_line_count "$acc")"
  [ "$total" -gt 0 ] || rp_die "fetch: $acc is empty; there is nothing to fetch"
  # Otherwise a non-numeric --jobs makes the pool's `-ge` test error rc 2,
  # which reads as "pool not full" and launches every chunk at once.
  [[ "$JOBS" =~ ^[1-9][0-9]*$ ]] \
    || rp_die "fetch: --jobs must be a positive integer, got '$JOBS'"
  # Checked here, with everything else, rather than where it is used: that
  # comparison is the last thing the stage does, so a bad value would otherwise
  # abort with a bare "unexpected error" after the whole 2 h 20 m fetch.
  [[ "$MAX_UNAVAILABLE" =~ ^[0-9]+$ ]] \
    || rp_die "fetch: RP_MAX_UNAVAILABLE must be a non-negative integer, got '$MAX_UNAVAILABLE'"

  # Recompute rather than trust the sidecar: this hash is both the resume key
  # and the stage marker, so it has to describe the file that is about to be
  # chunked, not what a previous stage said about it.
  hex="$(sha256sum "$acc" | cut -d' ' -f1)"
  recorded="$(cut -d' ' -f1 < "$sidecar")"
  [ "$hex" = "$recorded" ] \
    || rp_die "fetch: $acc hashes to $hex but $sidecar records $recorded; the accession list changed after the accessions stage wrote it"

  # The list's PROVENANCE, which no hash of the list can supply. rp_panel_key
  # folds the CURRENT QUERY with the list hash, so recomputing it here says
  # nothing about which panel definition the list was HARVESTED under --
  # stage_accessions is the only stage that knows, and it records that in its
  # own marker. Without this, editing QUERY and running --stage fetch directly
  # re-stamps the PREVIOUS panel's chunks with the new key, stage_build's
  # precondition then passes, and blastdb/<NEW_VER> ships the old record set
  # under a manifest quoting the new query, exit 0. A full run is saved only
  # incidentally, because stage_accessions re-harvests before this runs.
  # rp_marker_present, not rp_marker_ok: a tree with no accessions marker has no
  # provenance to contradict, and must still be fetchable.
  if rp_marker_present accessions \
      && ! rp_marker_ok accessions "$(rp_panel_key "$hex")"; then
    rp_die "fetch: $acc was harvested under a different panel definition than the one in effect now, so the chunks beside it hold the previous panel's record set. Re-harvest under the current definition with --stage accessions --force, then --stage fetch."
  fi

  # rp_panel_key, not the bare list hash: otherwise this marker inherits the
  # same staleness the accessions gate now refuses. An operator who edits QUERY
  # and runs --stage fetch would find a marker recorded under the previous panel
  # definition still matching, and skip.
  if [ "$FORCE" -eq 0 ] && rp_marker_ok fetch "$(rp_panel_key "$hex")"; then
    rp_log "fetch: already complete for this accession list and panel definition; skipping (--stage fetch --force re-runs it)"
    return 0
  fi

  chunks="$BUILD_DIR/chunks"
  missing="$BUILD_DIR/chunks.missing"
  tmp="$BUILD_DIR/chunks.tmp"
  key="$BUILD_DIR/chunks.accessions.sha256"

  if [ -f "$key" ]; then
    recorded="$(tr -d '[:space:]' < "$key")"
    [ "$recorded" = "$hex" ] \
      || rp_die "fetch: refusing to resume: $chunks was cut from an accession list hashing $recorded, but $acc now hashes $hex. The source moved. The chunks are reusable only with the exact list they were cut from, and both files have to agree: restore a copy hashing $recorded as $acc, rewrite the sidecar from it (cd $BUILD_DIR && sha256sum accessions.txt > accessions.sha256), then --stage fetch resumes without re-harvesting. The accessions stage overwrites $acc in place and keeps no backup, so if no copy was kept the chunks cannot be matched again: remove $chunks, $missing and $key and fetch fresh."
  else
    for f in "$chunks"/fa_[0-9][0-9][0-9][0-9]; do
      [ -e "$f" ] || break
      rp_die "fetch: refusing to resume: $chunks holds chunk(s) but $key is absent, so the accession list they were cut from cannot be proven. Remove $chunks to fetch it fresh."
    done
  fi

  mkdir -p "$chunks" "$missing"
  # All scratch lives here, never in $chunks: the historical build kept request
  # lists beside the chunks, where a loose fa_* glob would have eaten them.
  rm -rf "$tmp"
  mkdir -p "$tmp"
  # Written before the first chunk exists. A fetch interrupted at chunk 3 of
  # 179 has to be resumable, and chunks with no key beside them are refused.
  printf '%s\n' "$hex" > "$key"
  # Emptied before the fetch, rewritten after it: an abort must not leave an
  # earlier run's ledger sitting there claiming records were unavailable in a
  # fetch that never finished.
  : > "$BUILD_DIR/dropped.fetch.tsv"

  n_chunks=$(( (total + CHUNK_SIZE - 1) / CHUNK_SIZE ))
  rp_log "fetch: $total accession(s) in $n_chunks chunk(s) of $CHUNK_SIZE, $JOBS job(s)"

  running=0
  failed=0
  RP_FETCH_PIDS=""
  # Job control so each worker leads its own process group, which is what lets
  # an abort take a worker's efetch child down with it rather than orphan it.
  set -m
  trap rp_fetch_interrupt INT TERM
  for ((i = 0; i < n_chunks; i++)); do
    [ "$failed" -eq 0 ] || break
    idx="$(printf '%04d' "$i")"
    # Complete means BOTH files. A chunk without its missing list was written
    # by a worker that did not live to publish the other half.
    if [ -f "$chunks/fa_$idx" ] && [ -f "$missing/fa_$idx.txt" ]; then
      continue
    fi
    if [ "$running" -ge "$JOBS" ]; then
      rp_fetch_wait || failed=1
      running=$((running - 1))
    fi
    # `set -m` drops bash's implicit </dev/null for async commands, so without
    # this a worker inherits the script's stdin, and a background process group
    # that reads a terminal takes SIGTTIN.
    rp_fetch_chunk "$idx" $((i * CHUNK_SIZE + 1)) < /dev/null &
    RP_FETCH_PIDS="$RP_FETCH_PIDS $!"
    running=$((running + 1))
  done
  # On failure, kill and drain in one step rather than sitting through the rest
  # of an in-flight chunk. Nothing may exit while a worker is alive.
  while [ "$running" -gt 0 ]; do
    if [ "$failed" -ne 0 ]; then
      rp_fetch_reap
      break
    fi
    rp_fetch_wait || failed=1
    running=$((running - 1))
  done
  trap - INT TERM
  set +m
  RP_FETCH_PIDS=""
  [ "$failed" -eq 0 ] || rp_die "fetch: chunk worker(s) failed; see the errors above"

  allwant="$tmp/want.all"
  allgot="$tmp/got.all"
  notgot="$tmp/want_not_got.all"
  notwant="$tmp/got_not_want.all"
  sort -u "$acc" > "$allwant"
  # The strict glob, so a stray file in the chunk dir cannot join the panel.
  rp_fasta_accessions "$chunks"/fa_[0-9][0-9][0-9][0-9] > "$allgot"
  comm -23 "$allwant" "$allgot" > "$notgot"
  comm -13 "$allwant" "$allgot" > "$notwant"

  # The reverse reconcile, which the historical build never computed: a shrunk
  # accession list re-cut positionally shipped 12 withdrawn records in
  # simulation. Measured over the real 179-chunk build this is exactly 0 once
  # the PDB headers are normalised, so it gets no tolerance at all.
  [ ! -s "$notwant" ] \
    || rp_die "fetch: $(rp_line_count "$notwant") record(s) under $chunks were never requested, e.g. $(head -1 "$notwant")"

  # The per-chunk missing lists have to describe the chunks that are on disk
  # now, not the ones a previous run published. Resume trusts a chunk that has
  # both its files, and the dropped ledger is built from these lists, so a
  # chunk that changed underneath us would otherwise ship a panel whose funnel
  # silently does not add up.
  # -u because an accession duplicated in the list lands in two chunks, so it
  # would be recorded missing twice here but appear once in the reconcile, and
  # this guard would then abort a legitimate build.
  sort -u "$missing"/fa_[0-9][0-9][0-9][0-9].txt > "$tmp/missing.all"
  comm -3 "$notgot" "$tmp/missing.all" > "$tmp/missing.disagree"
  [ ! -s "$tmp/missing.disagree" ] \
    || rp_die "fetch: the chunk missing lists record $(rp_line_count "$tmp/missing.all") unavailable accession(s) but $(rp_line_count "$notgot") are actually absent from $chunks; a published chunk changed. Remove $chunks and $missing and fetch again."

  n_unavailable="$(rp_line_count "$notgot")"
  [ "$n_unavailable" -le "$MAX_UNAVAILABLE" ] \
    || rp_die "fetch: $n_unavailable accession(s) unavailable after $FETCH_ROUNDS round(s) per chunk, above the ceiling of $MAX_UNAVAILABLE; e.g. $(head -1 "$notgot")"

  # Rewritten wholesale from chunks.missing/ every run. Each stage owns its own
  # dropped-record file and the build stage merges them, so appending to a
  # shared dropped.tsv would let one stage's rewrite discard another's rows.
  awk -v OFS='\t' '{print $1, "efetch_unavailable", ""}' \
    "$missing"/fa_[0-9][0-9][0-9][0-9].txt > "$BUILD_DIR/dropped.fetch.tsv"

  n_got="$(rp_line_count "$allgot")"
  rp_marker_set fetch "$(rp_panel_key "$hex")"
  rp_log "fetch: $n_got record(s) in $n_chunks chunk(s) of $total requested, $n_unavailable permanently unavailable"
  rm -rf "$tmp"
}
# stage_build
# The panel: concatenate the verified chunks, dereplicate byte-exactly on
# (taxid, uppercased sequence), run makeblastdb, and write the shipped drop
# ledger and the manifest. Everything is assembled in scratch and cross-checked
# against the artifacts makeblastdb actually produced before anything is
# promoted into $OUT_DIR, the same way the fetch stage publishes a chunk only
# once it has verified it.
#
# The raw concatenated fasta is never landed: the historical 1.2 GB
# metazoa_rrna.raw.fasta was regenerated every run and is pure peak-disk waste,
# so the chunk stream is piped into each consumer instead.
stage_build() {
  rp_require_build_dir
  local v f tmp acc chunks hex requested final_count want_chunks n_chunks
  local stats n_kept sum_len n_taxids n_dropped n_unavailable n_downloaded
  local mean_defline reasons njs n_sample n_resolved has_taxids
  local built_utc blast_ver seqkit_ver pkg_ver git_sha mkcmd_json unpacked total
  local round settled n_excluded
  local -a mkcmd chunk_files

  # A stage reached by sourcing has none of these, and must refuse rather than
  # build into a path it invented.
  for v in VER OUT_DIR RP_SORT_MEM TMPDIR; do
    [ -n "${!v:-}" ] \
      || rp_die "build: $v is not set; run this script rather than sourcing it"
  done

  acc="$BUILD_DIR/accessions.txt"
  chunks="$BUILD_DIR/chunks"
  for f in "$acc" "$BUILD_DIR/accessions.requested_count" \
      "$BUILD_DIR/accessions.final_count" "$BUILD_DIR/taxidmap.tsv" \
      "$BUILD_DIR/dropped.accessions.tsv" "$BUILD_DIR/dropped.taxidmap.tsv" \
      "$BUILD_DIR/dropped.fetch.tsv"; do
    [ -f "$f" ] || rp_die "build: $f not found; run the earlier stages first"
  done

  # The pre-filter esearch <Count> is what the ledger reconciles against; the
  # post-filter count is what the chunks were cut from. They differ by the WGS
  # masters, so using one where the other belongs silently shifts the identity.
  requested="$(tr -d '[:space:]' < "$BUILD_DIR/accessions.requested_count")"
  final_count="$(tr -d '[:space:]' < "$BUILD_DIR/accessions.final_count")"
  [[ "$requested" =~ ^[0-9]+$ ]] && [[ "$final_count" =~ ^[0-9]+$ ]] \
    || rp_die "build: accessions.requested_count ('$requested') and accessions.final_count ('$final_count') must both be numeric"

  hex="$(sha256sum "$acc" | cut -d' ' -f1)"
  rp_marker_ok fetch "$(rp_panel_key "$hex")" \
    || rp_die "build: the fetch stage is not recorded complete for $acc (sha256 $hex) under this panel definition; run --stage fetch first. A marker recorded under a different QUERY or EXPECTED_COUNT does not match, deliberately: those chunks hold a different record set."

  # The strict glob, and the exact chunk count the accession list implies. An
  # empty directory passes the historical `0 -ne 0`, and a stray non-empty
  # fa_retry on its own builds and ships a 3-sequence database with a full
  # manifest.
  chunk_files=("$chunks"/fa_[0-9][0-9][0-9][0-9])
  [ -e "${chunk_files[0]}" ] \
    || rp_die "build: no chunk files matching $chunks/fa_[0-9][0-9][0-9][0-9]; run --stage fetch first"
  n_chunks="${#chunk_files[@]}"
  want_chunks=$(( (final_count + CHUNK_SIZE - 1) / CHUNK_SIZE ))
  rp_assert_count "build: chunk files under $chunks" "$n_chunks" "$want_chunks"

  tmp="$BUILD_DIR/build.tmp"
  rm -rf "$tmp"
  mkdir -p "$tmp"
  : > "$tmp/dup.acc"
  : > "$tmp/unusable.acc"

  # The PDB records' accessions, normalised by the one function that owns that
  # rule rather than by a second spelling of it here. `|| true` because zero
  # matching records is a legitimate state that bare errexit turned into a
  # failed build at the historical 03_build.sh:101.
  grep -h '^>pdb|' "${chunk_files[@]}" > "$tmp/pdb.deflines" || true
  rp_fasta_accessions "$tmp/pdb.deflines" > "$tmp/pdb.acc"

  # Accessions only; CR and blank lines tolerated, header skipped by name.
  : > "$tmp/exclude.acc"
  if [ -f "$RP_EXCLUSIONS" ]; then
    awk -F',' 'NR > 1 { gsub(/[ \r]/, "", $1); if ($1 != "" && $1 !~ /^#/) print $1 }' \
      "$RP_EXCLUSIONS" | sort -u > "$tmp/exclude.acc"
  fi

  rp_log "build: dereplicating $n_chunks chunk(s) on (taxid, uppercased sequence)"
  cat "${chunk_files[@]}" \
    | seqkit fx2tab -i \
    | awk -F'\t' -v OFS='\t' -v MAP="$BUILD_DIR/taxidmap.tsv" \
        -v GONE="$BUILD_DIR/dropped.taxidmap.tsv" -v BAD="$tmp/unusable.acc" \
        '
        BEGIN {
          while ((getline line < MAP) > 0) { split(line, a, "\t"); tax[a[1]] = a[2] }
          while ((getline line < GONE) > 0) { split(line, a, "\t"); gone[a[1]] = 1 }
        }
        {
          acc = $1
          # A pipe-form PDB defline never reaches the dedup, and is dropped
          # before it rather than after: 28ZU_S2 sorts ahead of any
          # letter-initial accession under LC_ALL=C, so a PDB record left in
          # would win the tie-break, be dropped anyway, and take the real
          # GenBank record it duplicated out of the panel with it. Their ledger
          # rows come from pdb.acc above.
          if (acc ~ /^pdb\|/) next
          # Already recorded as taxid_unresolved by the taxonomy stage, which
          # does not remove them from the accession list the fetch works from.
          if (acc in gone) next
          if ($2 == "") { print acc, "empty sequence" > BAD; next }
          if (!(acc in tax)) { print acc, "no taxid in taxidmap.tsv" > BAD; next }
          print tax[acc], toupper($2), acc
        }' \
    | sort -t$'\t' -k1,1 -k2,2 -k3,3 -S "$RP_SORT_MEM" -T "$TMPDIR" \
    | awk -F'\t' -v OFS='\t' -v DUP="$tmp/dup.acc" '
        { key = $1 FS $2
          if (key != prev) { prev = key; keeper = $3; print $3 }
          else { print $3, keeper > DUP } }' \
    > "$tmp/keep.all.acc"

  # keeper first, then its duplicates by the keeper column; both leave the
  # panel, the duplicates ledgered against the listed accession
  awk -F'\t' -v OFS='\t' -v EXCL="$tmp/exclude.acc" '
      BEGIN { while ((getline line < EXCL) > 0) excl[line] = 1 }
      ($1 in excl) { print $1, $1 } !($1 in excl) { print $1 > "/dev/stderr" }' \
    "$tmp/keep.all.acc" 2> "$tmp/keep.acc" > "$tmp/excluded.acc"
  awk -F'\t' -v OFS='\t' -v EXCL="$tmp/exclude.acc" '
      BEGIN { while ((getline line < EXCL) > 0) excl[line] = 1 }
      ($2 in excl) { print $1, $2 } !($2 in excl) { print $1, $2 > "/dev/stderr" }' \
    "$tmp/dup.acc" 2> "$tmp/dup.kept.acc" >> "$tmp/excluded.acc"
  mv "$tmp/dup.kept.acc" "$tmp/dup.acc"

  # Never a shared "NA" taxid for an unmapped record: that key collapses
  # unrelated species with byte-identical rRNA and mislabels them as duplicates
  # of each other.
  [ ! -s "$tmp/unusable.acc" ] \
    || rp_die "build: $(rp_line_count "$tmp/unusable.acc") fetched record(s) cannot enter the panel, e.g. $(head -1 "$tmp/unusable.acc" | tr '\t' ' ')"

  rp_log "build: writing the panel fasta and its restricted taxid map"
  cat "${chunk_files[@]}" \
    | awk -v KEEP="$tmp/keep.acc" '
        BEGIN { while ((getline line < KEEP) > 0) keep[line] = 1 }
        /^>/ { k = (substr($1, 2) in keep) }
        k' \
    > "$tmp/ribopilot_rrna.fasta"
  # Restricted to survivors so makeblastdb does not warn on absent ids.
  awk -F'\t' 'NR==FNR { k[$1]; next } ($1 in k)' \
    "$tmp/keep.acc" "$BUILD_DIR/taxidmap.tsv" > "$tmp/ribopilot_rrna.taxid.map"

  # Every count below is read back off a final artifact rather than carried
  # forward from an earlier step, so a write that lost records fails the ledger
  # instead of being reported as a success.
  stats="$(seqkit stats -T "$tmp/ribopilot_rrna.fasta" | awk -F'\t' 'NR==2 {print $4, $5}')"
  n_kept="${stats%% *}"
  sum_len="${stats##* }"
  n_taxids="$(cut -f2 "$tmp/ribopilot_rrna.taxid.map" | sort -u | wc -l | tr -d ' ')"
  mean_defline="$(awk '/^>/ {s += length($0); n++} END {printf "%d", n ? s / n : 0}' \
    "$tmp/ribopilot_rrna.fasta")"

  # No header row. The historical one made 293,693 lines read as 293,693
  # dropped records, 110 more than the manifest claimed.
  cat "$BUILD_DIR/dropped.accessions.tsv" "$BUILD_DIR/dropped.taxidmap.tsv" \
    "$BUILD_DIR/dropped.fetch.tsv" > "$tmp/dropped.all.tsv"
  awk -v OFS='\t' '{print $1, "pdb_structure_record_nonstandard_defline", ""}' \
    "$tmp/pdb.acc" >> "$tmp/dropped.all.tsv"
  awk -F'\t' -v OFS='\t' '{print $1, "byte_identical_to_kept_record_same_taxid", $2}' \
    "$tmp/dup.acc" >> "$tmp/dropped.all.tsv"
  # Only records the fetch actually held: a listed accession the harvest never
  # returned has no place in a ledger that must reconcile against the request.
  n_excluded="$(rp_line_count "$tmp/excluded.acc")"
  awk -F'\t' -v OFS='\t' '{print $1, "curated_exclusion", $2}' "$tmp/excluded.acc" >> "$tmp/dropped.all.tsv"
  rp_log "build: $n_excluded of $(rp_line_count "$tmp/exclude.acc") curated exclusion(s) found in the fetched records and dropped"
  # One row per accession, RP_DROP_PRECEDENCE deciding which reason survives.
  # The stage ledgers overlap, so concatenating them records the same record
  # twice -- the historical bug's exact class, and a build that aborts on
  # arithmetic after the 2 h 20 m fetch rather than naming the accession.
  awk -F'\t' -v OFS='\t' -v ORDER="$RP_DROP_PRECEDENCE" '
      BEGIN { n = split(ORDER, r, " "); for (i = 1; i <= n; i++) rank[r[i]] = i }
      { print (($2 in rank) ? rank[$2] : n + 1), $0 }' "$tmp/dropped.all.tsv" \
    | sort -t$'\t' -k2,2 -k1,1n \
    | awk -F'\t' -v OFS='\t' '$2 "" != prev "" { prev = $2; print $2, $3, $4 }' \
    > "$tmp/ribopilot_rrna.dropped.tsv"

  n_dropped="$(rp_line_count "$tmp/ribopilot_rrna.dropped.tsv")"
  reasons="$(awk -F'\t' '{c[$2]++} END {for (r in c) printf "%s=%d ", r, c[r]}' \
    "$tmp/ribopilot_rrna.dropped.tsv")"

  # Sets, not counts. The historical failure was one set of records counted
  # twice while another was missed, and a count identity balances silently
  # whenever those two happen to coincide. Comparing the sorted sets also names
  # the offending accession instead of reporting arithmetic.
  cut -f1 "$tmp/ribopilot_rrna.dropped.tsv" | sort > "$tmp/dropped.acc"
  [ -z "$(uniq -d "$tmp/dropped.acc" | head -1)" ] \
    || rp_die "build: the drop ledger records $(uniq -d "$tmp/dropped.acc" | head -1) more than once after collapsing by reason"
  # The requested set is the list the fetch worked from plus the masters
  # filtered out before it -- i.e. what was actually harvested, which is what
  # accessions.requested_count records. NOT esearch's <Count>: the source
  # reports more matches than its history returns, so reconciling against it
  # would demand rows for accessions that were never in any list.
  cut -f1 "$BUILD_DIR/dropped.accessions.tsv" > "$tmp/wgs.acc"
  sort -u "$tmp/wgs.acc" "$acc" > "$tmp/requested.acc"
  rp_assert_count "build: reconstructed requested set vs accessions.requested_count" \
    "$(rp_line_count "$tmp/requested.acc")" "$requested"
  rp_fasta_accessions "$tmp/ribopilot_rrna.fasta" > "$tmp/kept.acc"
  sort -u "$tmp/kept.acc" "$tmp/dropped.acc" > "$tmp/accounted.acc"
  comm -23 "$tmp/requested.acc" "$tmp/accounted.acc" > "$tmp/unaccounted.acc"
  comm -13 "$tmp/requested.acc" "$tmp/accounted.acc" > "$tmp/unrequested.acc"
  [ ! -s "$tmp/unaccounted.acc" ] \
    || rp_die "build: $(rp_line_count "$tmp/unaccounted.acc") requested accession(s) are neither in the panel nor in the drop ledger, e.g. $(head -1 "$tmp/unaccounted.acc")"
  [ ! -s "$tmp/unrequested.acc" ] \
    || rp_die "build: $(rp_line_count "$tmp/unrequested.acc") accession(s) in the panel or the ledger were never requested, e.g. $(head -1 "$tmp/unrequested.acc")"
  # With both sets equal and each unique, the counts can still disagree only if
  # an accession is on both sides at once, which is the double-count itself.
  [ "$((n_kept + n_dropped))" -eq "$requested" ] \
    || rp_die "build: the drop ledger does not balance: $requested requested != $n_kept kept + $n_dropped dropped ($reasons)"
  rp_log "build: $n_kept kept, $n_dropped dropped ($reasons)"

  # Relative paths, run from the panel directory, so the command recorded in the
  # manifest is re-runnable wherever the panel ends up.
  mkcmd=(makeblastdb -in ribopilot_rrna.fasta "${MAKEBLASTDB_FLAGS[@]}"
    -out ribopilot_rrna -title "RiboPilot Metazoa nrDNA $VER"
    -taxid_map ribopilot_rrna.taxid.map)
  rp_log "build: makeblastdb"
  ( cd "$tmp" && "${mkcmd[@]}" ) >&2 || rp_die "build: makeblastdb failed"

  # makeblastdb's exit status proves nothing: a truncated FASTA gives "added 2
  # sequences" and exit 0, and a taxid map matching nothing prints Error: and
  # still exits 0 with every staxid 0. The .njs it generates is the only
  # statement of what the database actually holds, so it is what gets checked.
  njs="$tmp/ribopilot_rrna.njs"
  [ -f "$njs" ] || rp_die "build: makeblastdb produced no $njs"
  rp_assert_count "build: .njs number-of-sequences vs panel fasta records" \
    "$(jq -r '.["number-of-sequences"]' "$njs")" "$n_kept"
  rp_assert_count "build: .njs number-of-letters vs seqkit sum_len" \
    "$(jq -r '.["number-of-letters"]' "$njs")" "$sum_len"
  rp_assert_count "build: .njs number-of-taxids vs the restricted taxid map" \
    "$(jq -r '.["number-of-taxids"]' "$njs")" "$n_taxids"

  # Sampled, never hardcoded: the historical manifest asserted has_taxids true
  # unconditionally, which is exactly what a taxid map matching nothing would
  # have shipped as.
  head -n 100 "$tmp/keep.acc" > "$tmp/sample.acc"
  n_sample="$(rp_line_count "$tmp/sample.acc")"
  ( cd "$tmp" && blastdbcmd -db ribopilot_rrna -entry_batch sample.acc -outfmt '%T' ) \
    > "$tmp/sample.taxid" || true
  n_resolved="$(rp_line_count "$tmp/sample.taxid")"
  has_taxids=false
  if [ "$n_sample" -gt 0 ] && [ "$n_resolved" = "$n_sample" ] \
      && ! grep -qvE '^[1-9][0-9]*$' "$tmp/sample.taxid"; then
    has_taxids=true
  fi
  [ "$has_taxids" = true ] \
    || rp_log "build: WARNING has_taxids is false: blastdbcmd returned $n_resolved line(s) for $n_sample sampled accession(s), not all of them a nonzero taxid; this panel cannot filter by taxonomy"

  n_unavailable="$(rp_line_count "$BUILD_DIR/dropped.fetch.tsv")"
  n_downloaded=$((final_count - n_unavailable))
  built_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  blast_ver="$(makeblastdb -version | awk 'NR==1 {print $2}')"
  seqkit_ver="$(seqkit version | awk '{print $2}')"
  pkg_ver="$(awk '/^Version:/ {print $2}' "$REPO_ROOT/DESCRIPTION")"
  git_sha="$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || true)"
  # A JSON array, not a string: the historical makeblastdb_cmd interpolated
  # -title RiboPilot Metazoa nrDNA 2026.08.1 unquoted and cannot be re-run.
  mkcmd_json="$(printf '%s\n' "${mkcmd[@]}" | jq -Rs 'split("\n")[:-1]')"

  # manifest.json is itself one of the shipped files, so its own size belongs in
  # unpacked_bytes: render it, re-measure, and settle once the number stops
  # moving.
  unpacked="$(rp_shipped_bytes "$tmp")"
  settled=0
  for round in 1 2 3; do
    jq -n \
      --arg db_version "$VER" \
      --arg built_utc "$built_utc" \
      --arg ribopilot_version "$pkg_ver" \
      --arg git_sha "$git_sha" \
      --arg platform "$(uname -sm)" \
      --arg query "$QUERY" \
      --arg taxid_source "$RP_A2T_BASE_URL" \
      --arg blast "$blast_ver" \
      --arg seqkit "$seqkit_ver" \
      --argjson requested "$requested" \
      --argjson downloaded "$n_downloaded" \
      --argjson unreconciled "$n_unavailable" \
      --argjson kept "$n_kept" \
      --argjson dropped "$n_dropped" \
      --argjson taxids "$n_taxids" \
      --argjson bases "$sum_len" \
      --argjson defline "$mean_defline" \
      --argjson unpacked "$unpacked" \
      --argjson has_taxids "$has_taxids" \
      --argjson makeblastdb_cmd "$mkcmd_json" \
      --argjson excluded "$n_excluded" \
      '{
        manifest_schema: 1,
        db_name: "ribopilot_rrna",
        db_version: $db_version,
        built_utc: $built_utc,
        builder: {
          ribopilot_version: $ribopilot_version,
          git_sha: $git_sha,
          platform: $platform
        },
        source: {
          provider: "NCBI Nucleotide (Entrez)",
          query: $query,
          molecule_filter: "none (RefSeq XR_/NR_ rRNA records are MolType=rna and would be dropped)",
          taxid_source: $taxid_source,
          requested_accessions: $requested,
          records_downloaded: $downloaded,
          unreconciled_accessions: $unreconciled
        },
        pipeline: {
          tier: 1,
          dedup: { mode: "byte_exact", key: ["taxid", "uppercased_sequence"] },
          annotated: false,
          species_selection: false,
          misid_screen: false,
          curated_exclusions: $excluded
        },
        counts: {
          records_after_dedup: $kept,
          records_dropped: $dropped,
          distinct_taxids: $taxids,
          total_bases: $bases,
          mean_defline_chars: $defline,
          unpacked_bytes: $unpacked
        },
        tools: { blast: $blast, seqkit: $seqkit },
        blastdb: {
          version: 5,
          prefix: "ribopilot_rrna",
          makeblastdb_cmd: $makeblastdb_cmd,
          defline_format: "verbatim GenBank: >ACC.V description",
          has_taxids: $has_taxids
        },
        decoy_panel: {
          mode: "ncbi_prebuilt",
          databases: ["SSU_eukaryote_rRNA", "LSU_eukaryote_rRNA", "16S_ribosomal_RNA", "ITS_eukaryote_sequences"],
          note: "composed at query time via a space-separated -db list; not merged into this build. ITS_eukaryote_sequences serves the GenBank BLAST check only and is never searched by the screen"
        }
      }' > "$tmp/manifest.json"
    total="$(rp_shipped_bytes "$tmp")"
    if [ "$total" = "$unpacked" ]; then settled=1; break; fi
    unpacked="$total"
  done
  # Falling out of the loop would ship a byte count describing a previous
  # render, which is the one thing every other check in this stage exists to
  # prevent.
  [ "$settled" -eq 1 ] \
    || rp_die "build: unpacked_bytes did not settle after 3 renders (last $unpacked, measured $total)"

  mkdir -p "$OUT_DIR"
  mv "$tmp/manifest.json" "$tmp"/ribopilot_rrna.* "$OUT_DIR/"
  rp_marker_set build "$hex"
  rp_log "build: $n_kept record(s), $n_taxids taxid(s), $sum_len base(s), $unpacked shipped byte(s) in $OUT_DIR"
  rm -rf "$tmp"
}
# stage_decoys
# The NCBI decoys are fetched into the panel directory itself, because
# local_blast_db_arg() (R/blast_local.R) defaults the decoy directory to the
# panel's own directory. fetch_decoys.sh owns the fetch, the post-fetch
# verification and the version capture; this stage only points it at $OUT_DIR.
#
# Deliberately not marker-gated: the decoys are a ~291 MB idempotent download
# that update_blastdb.pl skips when it is already current, and re-running the
# verification against the volumes on disk is the point of the stage.
stage_decoys() {
  rp_require_build_dir
  local v
  for v in VER OUT_DIR; do
    [ -n "${!v:-}" ] \
      || rp_die "decoys: $v is not set; run this script rather than sourcing it"
  done
  [ -d "$OUT_DIR" ] || rp_die "decoys: $OUT_DIR not found; run the build stage first"

  rp_log "decoys: fetching and verifying $DECOY_DBS in $OUT_DIR"
  bash "$SCRIPT_DIR/fetch_decoys.sh" "$OUT_DIR" >&2 \
    || rp_die "decoys: fetch_decoys.sh failed; see the errors above"
  # The marker records what landed, not what was asked for: DECOY_VERSION.txt
  # is written only after every volume, taxonomy file and blastdbcmd -info
  # check has passed.
  rp_marker_set decoys "$(sha256sum "$OUT_DIR/DECOY_VERSION.txt" | cut -d' ' -f1)"
  rp_log "decoys: verified, versions recorded in $OUT_DIR/DECOY_VERSION.txt"
}

# stage_package
# Stages the two release tarballs and their checksums. --publish is never
# passed: publishing is manual by spec section 7, so this stage stops at
# staged assets.
stage_package() {
  rp_require_build_dir
  local v rel f
  for v in VER OUT_DIR; do
    [ -n "${!v:-}" ] \
      || rp_die "package: $v is not set; run this script rather than sourcing it"
  done
  [ -d "$OUT_DIR" ] || rp_die "package: $OUT_DIR not found; run the build stage first"

  rel="$BUILD_DIR/blastdb/release-$VER"
  bash "$SCRIPT_DIR/package_release.sh" "$OUT_DIR" "$VER" >&2 \
    || rp_die "package: package_release.sh failed; see the errors above"

  # Named rather than left to the marker's own read below, so a missing asset
  # aborts with the asset's name instead of a bare unexpected-error. The
  # tarballs are not re-verified here: clean checksums them again before it
  # deletes anything on their strength, which is where that check belongs.
  for f in "ribopilot_rrna-$VER.tar.gz" "decoys-$VER.tar.gz" "dropped-$VER.tsv.gz"; do
    [ -f "$rel/$f" ] || rp_die "package: $rel/$f was not produced"
    [ -f "$rel/$f.sha256" ] || rp_die "package: $rel/$f.sha256 was not produced"
  done
  rp_marker_set package "$(cut -d' ' -f1 < "$rel/ribopilot_rrna-$VER.tar.gz.sha256")"
  rp_log "package: release assets staged in $rel (publish manually with package_release.sh --publish)"
}

run_stage() {
  case "$1" in
    accessions) stage_accessions ;;
    taxidmap) stage_taxidmap ;;
    fetch) stage_fetch ;;
    build) stage_build ;;
    decoys) stage_decoys ;;
    package) stage_package ;;
    *) rp_die "unknown stage: $1" ;;
  esac
}

main_build() {
  local s
  run_preflight
  mkdir -p "$BUILD_DIR"
  rp_lock "$BUILD_DIR/build.pid"
  RP_LOCK_HELD="$BUILD_DIR/build.pid"

  if [ -n "$STAGE" ]; then
    run_stage "$STAGE"
  else
    for s in $ALL_STAGES; do
      run_stage "$s"
    done
  fi
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
# Guarded so tests can `source` this file (functions and constants only, no
# side effects) to exercise guard functions directly, the same way lib.sh is
# sourced in isolation elsewhere in the test suite.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  if [ "$CHECK_VERSION" -eq 1 ]; then
    cmd_check_version
    exit $?
  fi

  case "$COMMAND" in
    status)
      cmd_status
      exit 0
      ;;
    clean)
      cmd_clean
      exit 0
      ;;
  esac

  main_build
fi
