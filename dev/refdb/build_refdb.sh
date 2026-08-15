#!/usr/bin/env bash
# Orchestrator for the RiboPilot metazoan rDNA BLAST panel build.
#
# Usage: build_refdb.sh [--version V] [--build-dir D] [--stage NAME]
#                        [--jobs N] [--accept-drift] [--force] [--check-version]
#                        [status|clean]
#
# See dev/design/18_REFDB_REBUILD.md for the spec this implements. This file
# is the skeleton: argument parsing, preflight, and version sync. The build
# stages themselves (accessions, taxidmap, fetch, build, decoys, package) land
# in later work and are stubs here.
set -Eeuo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=dev/refdb/lib.sh
source "$SCRIPT_DIR/lib.sh"

# Before anything can print: the canonical invocation is
# `./build_refdb.sh 2>&1 | tee build.log`, and the key must not reach that log.
rp_scrub_output

RP_LOCK_HELD=""
rp_cleanup() {
  [ -z "$RP_LOCK_HELD" ] || rm -f "$RP_LOCK_HELD"
  rp_scrub_flush
}
trap rp_cleanup EXIT
trap 'rp_die "unexpected error at ${BASH_SOURCE[0]}:${LINENO}"' ERR

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
DECOY_DBS="16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb"
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

REQUIRED_TOOLS="esearch efetch epost esummary xtract update_blastdb.pl makeblastdb blastdbcmd seqkit jq curl sha256sum tar zcat comm"
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
  --jobs N          Parallelism for the fetch stage (default: 1).
  --accept-drift    Allow the harvested count to differ from EXPECTED_COUNT.
  --force           Allow rebuilding an existing blastdb/<VER> in place.
  --check-version   Diff RP_DB_VERSION across the build's version-carrying files, then exit.
  status            Print the last recorded stage state and exit.
  clean             Remove intermediate build artifacts and exit.
EOF
}

VERSION_ARG=""
BUILD_DIR="$SCRIPT_DIR/build"
STAGE=""
JOBS=1
ACCEPT_DRIFT=0
FORCE=0
CHECK_VERSION=0
COMMAND=""

while [ $# -gt 0 ]; do
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

RP_STATE_FILE="$BUILD_DIR/state.txt"

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
# Spec section 6: drops chunks/, the raw FASTA, and the unshipped .fasta,
# .taxid.map and .dropped.tsv -- but only once the release tarballs verify,
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

  for f in "ribopilot_rrna-$ver.tar.gz" "decoys-$ver.tar.gz"; do
    [ -f "$rel/$f" ] \
      || rp_die "refusing to clean: no $rel/$f (run package_release.sh first)"
    [ -f "$rel/$f.sha256" ] \
      || rp_die "refusing to clean: no $rel/$f.sha256"
    (cd "$rel" && sha256sum -c --status "$f.sha256") \
      || rp_die "refusing to clean: $rel/$f does not match $rel/$f.sha256"
  done
  rp_log "release tarballs verify; cleaning intermediates under $BUILD_DIR"

  rm -rf "$BUILD_DIR/chunks"
  rm -f "$BUILD_DIR"/*.raw.fasta
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
  if [ -d "$OUT_DIR" ] && [ "$FORCE" -ne 1 ]; then
    rp_die "$OUT_DIR already exists; pass --force to rebuild it in place"
  fi
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

# rp_filter_wgs_masters IN OUT
# WGS master accessions (4 letters + 6 zeros) are permanently unavailable via
# efetch by design, so dropping them here is a recorded expected drop, not a
# silent one.
rp_filter_wgs_masters() {
  local in="$1" out="$2" total kept
  grep -Ev '^[A-Z]{4}0{6}' "$in" > "$out" || true
  total="$(rp_line_count "$in")"
  kept="$(rp_line_count "$out")"
  rp_log "accessions: dropped $((total - kept)) WGS master accession(s) (expected, permanently unavailable)"
}

# rp_check_drift ACTUAL ACCEPT_DRIFT
# Beyond EXPECTED_COUNT_TOLERANCE, refuse unless ACCEPT_DRIFT=1, which logs
# the new count instead of aborting. Within tolerance, silent: the source
# moves ~370 records/day and that is not news.
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
# Stages -- not implemented yet, land in later work. Each is a hard stub so a
# partial run cannot be mistaken for a real one.
# ---------------------------------------------------------------------------

# stage_accessions
# One esearch against nuccore with the verbatim QUERY. The XML is captured
# once and replayed into efetch (-input, not a second esearch), so the
# accession list and the <Count> it is checked against come from the same
# history-server result set.
stage_accessions() {
  local xml raw filtered count

  mkdir -p "$BUILD_DIR"
  xml="$BUILD_DIR/accessions.esearch.xml"
  raw="$BUILD_DIR/accessions.raw.txt"
  filtered="$BUILD_DIR/accessions.txt"

  rp_log "accessions: esearch (nuccore) for the panel query"
  esearch -db nuccore -query "$QUERY" > "$xml"

  # esearch's default output is its own ENTREZ_DIRECT chaining wrapper
  # (WebEnv/QueryKey/Count), not raw eSearchResult XML -- verified against a
  # live esearch, not assumed.
  count="$(xtract -input "$xml" -pattern ENTREZ_DIRECT -element Count)"
  [ -n "$count" ] && [[ "$count" =~ ^[0-9]+$ ]] \
    || rp_die "accessions: could not parse a numeric <Count> from esearch output"

  rp_check_drift "$count" "$ACCEPT_DRIFT"

  rp_log "accessions: replaying the same WebEnv/QueryKey into efetch -format acc ($count expected)"
  efetch -input "$xml" -format acc > "$raw"

  rp_assert_count "accessions: harvested line count vs esearch <Count>" "$(rp_line_count "$raw")" "$count"
  rp_assert_unique "accessions: harvested accession list" "$raw"

  rp_filter_wgs_masters "$raw" "$filtered"

  (cd "$BUILD_DIR" && sha256sum "$(basename "$filtered")" > accessions.sha256)
  printf '%s\n' "$count" > "$BUILD_DIR/accessions.count"
  rp_marker_set accessions "$(cut -d' ' -f1 < "$BUILD_DIR/accessions.sha256")"
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

  # Only now, with both guards passed, commit. dropped.tsv is written fresh
  # each run (never appended across retries): the stage is spec'd idempotent
  # ("re-run from scratch, cheap"), and appending would let a retried run
  # accumulate stale unresolved rows from a prior attempt.
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
  cp "$dropped_tmp" "$BUILD_DIR/dropped.tsv"

  cp "$candidate_map" "$BUILD_DIR/taxidmap.tsv"
  cp "$candidate_sources" "$BUILD_DIR/taxidmap.sources.tsv"
  (cd "$BUILD_DIR" && sha256sum taxidmap.tsv > taxidmap.sha256)
  rp_marker_set taxidmap "$(cut -d' ' -f1 < "$BUILD_DIR/taxidmap.sha256")"

  ftp_n="$(awk -F'\t' '$3=="ftp"' "$combined" | wc -l | tr -d ' ')"
  esummary_n="$(awk -F'\t' '$3=="esummary"' "$combined" | wc -l | tr -d ' ')"
  rp_log "taxidmap: wrote $mapped mapped accession(s) ($ftp_n ftp, $esummary_n esummary), $unresolved_n unresolved (dropped) of $requested requested"

  rm -rf "$tmp"
}
stage_fetch() { rp_die "stage 'fetch' is not implemented yet"; }
stage_build() { rp_die "stage 'build' is not implemented yet"; }
stage_decoys() { rp_die "stage 'decoys' is not implemented yet"; }
stage_package() { rp_die "stage 'package' is not implemented yet"; }

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
