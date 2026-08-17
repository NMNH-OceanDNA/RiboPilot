#!/bin/bash
# Build the RiboPilot Docker image locally.
# Usage: docker/deploy-local.sh [tag] [--no-db]
# Run from the repository root. Use docker/deploy-dockerhub.sh to publish.
#
# The tag defaults to the DESCRIPTION Version rather than "latest": a floating
# tag makes "which code is in this image" unanswerable, and it is the question
# you ask first when a run misbehaves.
#
# --no-db builds without the rDNA panel, for when the panel release has not been
# published yet. That image searches NCBI remotely and says so.

set -euo pipefail

tag=""
no_db=0
for arg in "$@"; do
  case "$arg" in
    --no-db) no_db=1 ;;
    -*)      echo "Unknown option: $arg" >&2; exit 2 ;;
    *)       tag="$arg" ;;
  esac
done
tag="${tag:-$(awk '/^Version:/{print $2}' DESCRIPTION)}"
repo="ribopilot"

if [ ! -f DESCRIPTION ] || [ ! -f docker/Dockerfile ]; then
  echo "Run this from the repository root (no DESCRIPTION / docker/Dockerfile here)." >&2
  exit 1
fi

# Build the R package tarball into docker/ for the Dockerfile to COPY. Clear old
# ones first: the Dockerfile COPYs docker/RiboPilot_*.tar.gz to a single
# destination, so a tarball left over from another Version breaks the build.
# NOTE: no devtools::document() here -- RiboPilot's NAMESPACE is hand-maintained.
rm -f docker/RiboPilot_*.tar.gz
Rscript -e 'devtools::build(path = "docker", vignettes = FALSE)'

# db panel tarballs: fetched once and reused, these are 82 MB and 77 MB, so a
# per-build download would make every build slow for no reason. Every asset,
# freshly fetched or already on disk, is checksum-verified before the build:
# Docker's ADD does not reliably fail on a corrupt-but-gzip-shaped tarball, so
# an unverified stale or interrupted download would get baked in silently.
# gh, not curl: the repo is private, so an unauthenticated fetch 404s.
RP_DB_VERSION=$(awk -F'"' '/^RP_DB_VERSION/{print $2}' R/db_dist.R)
RP_DB_REPO="dmacguigan/RiboPilot"
RP_DB_TAG="db-${RP_DB_VERSION}"
RP_DB=1

# Sweep staged assets from OTHER panel versions, so no artifact carrying a stale
# version is left in the build context. The current version's tarballs stay:
# they are the cache described above, and the checksum below is a stronger
# guarantee against a corrupt one than deleting and re-downloading 160 MB.
for f in docker/ribopilot_rrna-*.tar.gz docker/ribopilot_rrna-*.tar.gz.sha256 \
         docker/decoys-*.tar.gz docker/decoys-*.tar.gz.sha256; do
  [ -e "$f" ] || continue
  case "$f" in
    *"-${RP_DB_VERSION}.tar.gz" | *"-${RP_DB_VERSION}.tar.gz.sha256") continue ;;
  esac
  echo "removing stale staged asset ${f#docker/}"
  rm -f "$f"
done

rp_db_unavailable() {
  cat >&2 <<MSG
Could not get $1 from release ${RP_DB_TAG} of ${RP_DB_REPO}.
The panel release probably has not been published yet. Do one of:
  dev/refdb/package_release.sh <panel_dir> ${RP_DB_VERSION} --publish
  # or place ribopilot_rrna-${RP_DB_VERSION}.tar.gz, decoys-${RP_DB_VERSION}.tar.gz
  # and their .sha256 files in docker/ by hand
  # or build an image with no panel in it:
  $0 --no-db
MSG
  exit 1
}

if [ "$no_db" -eq 1 ]; then
  RP_DB=0
  echo "Building WITHOUT the rDNA panel (--no-db): this image searches NCBI remotely."
else
  if command -v sha256sum >/dev/null 2>&1; then
    RP_SHA_CMD=(sha256sum)
  elif command -v shasum >/dev/null 2>&1; then
    RP_SHA_CMD=(shasum -a 256)
  else
    echo "Need sha256sum or shasum on PATH to verify the db panel tarballs." >&2
    exit 1
  fi

  for A in "ribopilot_rrna-${RP_DB_VERSION}.tar.gz" "decoys-${RP_DB_VERSION}.tar.gz"; do
    if [ ! -f "docker/$A" ]; then
      command -v gh >/dev/null 2>&1 || rp_db_unavailable "$A"
      echo "fetching $A ..."
      gh release download "$RP_DB_TAG" --repo "$RP_DB_REPO" --pattern "$A" \
        --output "docker/$A" --clobber || rp_db_unavailable "$A"
    fi
    sumfile="$(mktemp)"
    # a hand-placed tarball brings its own .sha256; only reach for the release
    # when there is none beside it
    if [ -f "docker/$A.sha256" ]; then
      cp "docker/$A.sha256" "$sumfile"
    else
      command -v gh >/dev/null 2>&1 || rp_db_unavailable "$A.sha256"
      gh release download "$RP_DB_TAG" --repo "$RP_DB_REPO" --pattern "$A.sha256" \
        --output "$sumfile" --clobber || rp_db_unavailable "$A.sha256"
    fi
    if ! ( cd docker && "${RP_SHA_CMD[@]}" -c "$sumfile" ); then
      echo "Checksum mismatch for $A: the file is corrupt or truncated. Removing it so the next run re-fetches it." >&2
      rm -f "docker/$A" "$sumfile"
      exit 1
    fi
    rm -f "$sumfile"
  done
fi

docker build --platform linux/amd64 -f docker/Dockerfile --progress=plain \
  --build-arg "RP_DB=${RP_DB}" \
  --build-arg "RP_DB_VERSION=${RP_DB_VERSION}" \
  -t "${repo}:${tag}" .

if [ "$RP_DB" -eq 1 ]; then
  echo "Local image build successful: ${repo}:${tag} (rDNA panel ${RP_DB_VERSION} baked in)"
else
  echo "Local image build successful: ${repo}:${tag} (NO rDNA panel; searches NCBI remotely)"
fi
