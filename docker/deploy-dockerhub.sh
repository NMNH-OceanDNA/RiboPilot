#!/bin/bash
# Build the RiboPilot image and push it to Docker Hub as macguigand/ribopilot.
# Usage, from the repository root:
#
#   docker/deploy-dockerhub.sh                 # tag = DESCRIPTION Version
#   docker/deploy-dockerhub.sh v0.1.0          # explicit tag
#   docker/deploy-dockerhub.sh --dry-run       # build + test, push nothing
#   docker/deploy-dockerhub.sh --no-db         # build without the rDNA panel
#
# There is deliberately no :latest. A floating tag means two machines can run
# different code while naming the same image, and every project pins an exact
# version anyway (new_project() reads it from DESCRIPTION).
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
DRY_RUN=0
NO_DB=0

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --no-db)   NO_DB=1 ;;
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

# --- db panel tarballs -------------------------------------------------------
# Fetched once and reused: these are 82 MB and 77 MB, so a per-build download
# would make every build slow for no reason. Not covered by the rm -f above.
# Every asset, freshly fetched or already on disk, is checksum-verified before
# the build: Docker's ADD does not reliably fail on a corrupt-but-gzip-shaped
# tarball, so an unverified stale or interrupted download would get baked into
# the image silently, and every user of it would then search a broken database.
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

if [ "$NO_DB" -eq 1 ]; then
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

# --- build -----------------------------------------------------------------
# Pinned to linux/amd64: the conda channels this image installs from do not
# carry every tool for arm64, and an image that silently builds for the host
# architecture is an image that fails on somebody else's cluster.
docker build \
  --platform linux/amd64 \
  -f docker/Dockerfile \
  --progress=plain \
  --build-arg "GIT_VERSION_TAG=${GIT_SHA}" \
  --build-arg "RP_DB=${RP_DB}" \
  --build-arg "RP_DB_VERSION=${RP_DB_VERSION}" \
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
           efetch sqlite3 Rscript; do
    command -v "$t" >/dev/null || { echo "MISSING TOOL: $t" >&2; exit 1; }
  done
  Rscript -e "suppressMessages(library(RiboPilot)); cat(\"RiboPilot\", as.character(packageVersion(\"RiboPilot\")), \"loads OK\n\")"
'
# No pipe here: a piped command returns 141 on this machine even when it
# succeeds, which would pass this check for the wrong reason.
if [ "$RP_DB" -eq 1 ]; then
  docker run --rm --platform linux/amd64 "${LOCAL}" bash -lc \
    'blastdbcmd -db "$RIBOPILOT_DB_DIR/'"${RP_DB_VERSION}"'/ribopilot_rrna" -info > /dev/null'
  echo "Image looks sane (rDNA panel ${RP_DB_VERSION} baked in)."
else
  # a --no-db image must not advertise a panel it does not carry
  docker run --rm --platform linux/amd64 "${LOCAL}" bash -lc \
    '[ -z "${RIBOPILOT_DB_DIR:-}" ]'
  echo "Image looks sane (NO rDNA panel; searches NCBI remotely)."
fi

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

echo
echo "Projects can use it with:"
echo "  new_project(..., executor = \"docker\", container = \"${REMOTE}\")"
