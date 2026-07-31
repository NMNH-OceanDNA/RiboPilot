#!/bin/bash
# Build the RiboPilot Docker image locally. Usage: docker/deploy-local.sh [tag]
# Run from the repository root. Use docker/deploy-dockerhub.sh to publish.

set -euo pipefail

tag="${1:-latest}"
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

docker build --platform linux/amd64 -f docker/Dockerfile --progress=plain -t "${repo}:${tag}" .

echo "Local image build successful: ${repo}:${tag}"
