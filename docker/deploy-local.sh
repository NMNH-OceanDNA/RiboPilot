#!/bin/bash
# Build the RiboPilot Docker image locally. Usage: docker/deploy-local.sh [tag]
# Run from the repository root. Use docker/deploy-dockerhub.sh to publish.
#
# The tag defaults to the DESCRIPTION Version rather than "latest": a floating
# tag makes "which code is in this image" unanswerable, and it is the question
# you ask first when a run misbehaves.

set -euo pipefail

tag="${1:-$(awk '/^Version:/{print $2}' DESCRIPTION)}"
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
RP_DB_VERSION=$(awk -F'"' '/^RP_DB_VERSION/{print $2}' R/db_dist.R)
RP_DB_BASE="https://github.com/dmacguigan/RiboPilot/releases/download/db-${RP_DB_VERSION}"

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
    echo "fetching $A ..."
    curl -fsSL -o "docker/$A" "${RP_DB_BASE}/$A"
  fi
  sumfile="$(mktemp)"
  curl -fsSL -o "$sumfile" "${RP_DB_BASE}/$A.sha256"
  if ! ( cd docker && "${RP_SHA_CMD[@]}" -c "$sumfile" ); then
    echo "Checksum mismatch for $A: the file is corrupt or truncated. Removing it so the next run re-fetches it." >&2
    rm -f "docker/$A" "$sumfile"
    exit 1
  fi
  rm -f "$sumfile"
done

docker build --platform linux/amd64 -f docker/Dockerfile --progress=plain \
  --build-arg "RP_DB_VERSION=${RP_DB_VERSION}" \
  -t "${repo}:${tag}" .

echo "Local image build successful: ${repo}:${tag}"
