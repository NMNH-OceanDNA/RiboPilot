#!/bin/bash
# Build the RiboPilot image and push it to Docker Hub as macguigand/ribopilot.
# Usage, from the repository root:
#
#   docker/deploy-dockerhub.sh                 # tag = DESCRIPTION Version
#   docker/deploy-dockerhub.sh v0.1.0          # explicit tag
#   docker/deploy-dockerhub.sh v0.1.0 --latest # also push :latest
#   docker/deploy-dockerhub.sh --dry-run       # build + test, push nothing
#
# Environment:
#   DOCKERHUB_USER   push under a different account (default macguigand)
#   DOCKERHUB_TOKEN  access token, so the script never prompts
#
# The tag defaults to the DESCRIPTION Version because new_project() derives its
# default `container` from the same place. Let those two drift and every project
# asks for a tag that was never published.

set -euo pipefail

NAMESPACE="${DOCKERHUB_USER:-macguigand}"
REPO='ribopilot'
TAG=''
ALSO_LATEST=0
DRY_RUN=0

for arg in "$@"; do
  case "$arg" in
    --latest)  ALSO_LATEST=1 ;;
    --dry-run) DRY_RUN=1 ;;
    -*)        echo "Unknown option: $arg" >&2; exit 2 ;;
    *)         TAG="$arg" ;;
  esac
done

if [ ! -f DESCRIPTION ] || [ ! -f docker/Dockerfile ]; then
  echo "Run this from the repository root (no DESCRIPTION / docker/Dockerfile here)." >&2
  exit 1
fi

TAG="${TAG:-$(awk '/^Version:/{print $2}' DESCRIPTION)}"
if [ -z "$TAG" ]; then
  echo "Could not determine a tag: no argument, and no Version in DESCRIPTION." >&2
  exit 1
fi

LOCAL="${REPO}:${TAG}"
REMOTE="${NAMESPACE}/${REPO}:${TAG}"

# A published image should be traceable to a commit. A dirty tree is not fatal
# (you may be publishing a deliberate one-off) but it must not pass silently.
GIT_SHA="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
  GIT_SHA="${GIT_SHA}-dirty"
  echo "WARNING: the working tree is dirty. This image will not match any commit."
fi

echo "Image:  ${REMOTE}"
echo "Commit: ${GIT_SHA}"
[ "$ALSO_LATEST" -eq 1 ] && echo "Also:   ${NAMESPACE}/${REPO}:latest"
[ "$DRY_RUN" -eq 1 ] && echo "Mode:   dry run, nothing will be pushed"

# --- package tarball -------------------------------------------------------
# The Dockerfile COPYs docker/RiboPilot_*.tar.gz to a single destination, so a
# leftover tarball from an older Version breaks the build (two sources, one
# non-directory destination) or, worse, installs the wrong one. Clear them.
# NOTE: no devtools::document() -- RiboPilot's NAMESPACE is hand-maintained.
rm -f docker/RiboPilot_*.tar.gz
Rscript -e 'devtools::build(path = "docker", vignettes = FALSE)'

count=$(ls -1 docker/RiboPilot_*.tar.gz 2>/dev/null | wc -l)
if [ "$count" -ne 1 ]; then
  echo "Expected exactly one package tarball in docker/, found ${count}." >&2
  exit 1
fi

# --- build -----------------------------------------------------------------
# Pinned to linux/amd64: the conda channels this image installs from do not
# carry every tool for arm64, and an image that silently builds for the host
# architecture is an image that fails on somebody else's cluster.
docker build \
  --platform linux/amd64 \
  -f docker/Dockerfile \
  --progress=plain \
  --build-arg "GIT_VERSION_TAG=${GIT_SHA}" \
  --label "org.opencontainers.image.revision=${GIT_SHA}" \
  --label "org.opencontainers.image.version=${TAG}" \
  --label "org.opencontainers.image.source=https://github.com/dmacguigan/RiboPilot" \
  -t "${LOCAL}" .

# --- smoke test ------------------------------------------------------------
# Publishing is the one step that is awkward to undo, so check the image can do
# what the pipeline immediately asks of it. Several processes degrade to a
# SILENT no-op when a tool is missing, which is exactly the failure worth
# catching here rather than three hours into somebody else's run.
echo "Testing ${LOCAL}..."
docker run --rm --platform linux/amd64 "${LOCAL}" bash -lc '
  set -e
  for t in fastp get_organelle_from_reads.py spades.py blastn makeblastdb \
           bowtie2 bowtie2-build samtools seqkit barrnap ITSx hmmscan \
           efetch taxonkit sqlite3 Rscript; do
    command -v "$t" >/dev/null || { echo "MISSING TOOL: $t" >&2; exit 1; }
  done
  Rscript -e "suppressMessages(library(RiboPilot)); cat(\"RiboPilot\", as.character(packageVersion(\"RiboPilot\")), \"loads OK\n\")"
'
echo "Image looks sane."

if [ "$DRY_RUN" -eq 1 ]; then
  echo "Dry run: built and tested ${LOCAL}, pushed nothing."
  exit 0
fi

# --- login -----------------------------------------------------------------
# Reuse an existing session where there is one; only ask when there is not.
if ! docker system info 2>/dev/null | grep -q "Username: ${NAMESPACE}"; then
  if [ -n "${DOCKERHUB_TOKEN:-}" ]; then
    echo "Logging in as ${NAMESPACE} with DOCKERHUB_TOKEN..."
    echo "${DOCKERHUB_TOKEN}" | docker login -u "${NAMESPACE}" --password-stdin
  else
    echo "Logging in as ${NAMESPACE}. Use an access token, not your password:"
    echo "  https://app.docker.com/settings/personal-access-tokens"
    docker login -u "${NAMESPACE}"
  fi
fi

# --- push ------------------------------------------------------------------
read -r -p "Push ${REMOTE} to Docker Hub? [y/N] " reply
case "$reply" in
  [yY]|[yY][eE][sS]) ;;
  *) echo "Aborted. The local image ${LOCAL} is still here."; exit 0 ;;
esac

docker tag "${LOCAL}" "${REMOTE}"
docker push "${REMOTE}"
echo "Pushed ${REMOTE}"

if [ "$ALSO_LATEST" -eq 1 ]; then
  docker tag "${LOCAL}" "${NAMESPACE}/${REPO}:latest"
  docker push "${NAMESPACE}/${REPO}:latest"
  echo "Pushed ${NAMESPACE}/${REPO}:latest"
fi

echo
echo "Projects can use it with:"
echo "  new_project(..., executor = \"docker\", container = \"${REMOTE}\")"
