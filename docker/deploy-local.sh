#!/bin/bash
# Build the RiboPilot Docker image locally. Usage: docker/deploy-local.sh [tag]
# Run from the repository root.

tag="${1:-latest}"
repo="ribopilot"

# Build the R package tarball into docker/ for the Dockerfile to COPY.
# NOTE: no devtools::document() here -- RiboPilot's NAMESPACE is hand-maintained.
Rscript -e 'devtools::build(path = "docker", vignettes = FALSE)'
if [ $? -ne 0 ]; then
  echo "Failed to build the R package tarball"
  exit 1
fi

docker build -f docker/Dockerfile --progress=plain -t ${repo}:${tag} .
if [ $? -ne 0 ]; then
  echo "Failed to build the Docker image"
  exit 1
fi

echo "Local image build successful: ${repo}:${tag}"
