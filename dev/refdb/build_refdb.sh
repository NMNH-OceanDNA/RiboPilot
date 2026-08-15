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

REQUIRED_TOOLS="esearch efetch epost esummary xtract update_blastdb.pl makeblastdb blastdbcmd seqkit jq curl sha256sum tar"
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
  local ver docker_ver ok=1 out_dir manifest_ver f base pat fver sumfile

  ver="$(rp_read_version)"
  rp_valid_version "$ver" || rp_die "RP_DB_VERSION '$ver' is not dotted-numeric"

  docker_ver="$(grep '^ARG RP_DB_VERSION=' "$REPO_ROOT/docker/Dockerfile" || true)"
  docker_ver="${docker_ver#ARG RP_DB_VERSION=}"
  if [ "$docker_ver" != "$ver" ]; then
    rp_log "MISMATCH docker/Dockerfile ARG RP_DB_VERSION=$docker_ver != $ver"
    ok=0
  fi

  out_dir="$BUILD_DIR/blastdb/$ver"
  if [ -f "$out_dir/manifest.json" ]; then
    manifest_ver="$(jq -r '.db_version // empty' "$out_dir/manifest.json" 2>/dev/null || true)"
    if [ "$manifest_ver" != "$ver" ]; then
      rp_log "MISMATCH $out_dir/manifest.json db_version=$manifest_ver != $ver"
      ok=0
    fi
  fi

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
# Removes cheap-to-regenerate intermediates only. Never touches blastdb/<VER>
# output; overwriting that is gated separately by --force.
cmd_clean() {
  rp_log "cleaning intermediate build artifacts under $BUILD_DIR"
  rm -rf "$BUILD_DIR/chunks"
  rm -f "$BUILD_DIR"/*.raw.fasta
  rm -f "$RP_STATE_FILE"
  rp_log "clean complete (blastdb/ output, if any, is left in place)"
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
# Stages -- not implemented yet, land in later work. Each is a hard stub so a
# partial run cannot be mistaken for a real one.
# ---------------------------------------------------------------------------

stage_accessions() { rp_die "stage 'accessions' is not implemented yet"; }
stage_taxidmap() { rp_die "stage 'taxidmap' is not implemented yet"; }
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
  trap 'rm -f "$BUILD_DIR/build.pid"' EXIT

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
