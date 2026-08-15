#!/usr/bin/env bash
# Package a built panel (and its decoys) as GitHub release assets.
#
# Ships the BLAST volume set and the manifest only. The .fasta, .taxid.map and
# .dropped.tsv in the build directory are build artifacts, not query-time files.
set -Eeuo pipefail
export LC_ALL=C
trap 'echo "package_release.sh: error at line $LINENO" >&2' ERR

SRC="${1:?usage: package_release.sh <panel_dir> <version> [--publish]}"
VER="${2:?usage: package_release.sh <panel_dir> <version> [--publish]}"
PUBLISH="${3:-}"
[ -d "$SRC" ] || { echo "no such panel directory: $SRC" >&2; exit 1; }
# Absolute, because the script cd's into SRC before writing anything
OUT="$(cd "${SRC%/}/.." && pwd)/release-${VER}"
REPO="dmacguigan/RiboPilot"

mkdir -p "$OUT"
cd "$SRC"

PANEL_FILES="manifest.json"
for e in ndb nhr nin njs nog nos not nsq ntf nto; do
  [ -f "ribopilot_rrna.$e" ] && PANEL_FILES="$PANEL_FILES ribopilot_rrna.$e"
done
# PANEL_REQUIRED_EXT in R/blast_local.R. Without these the installer rejects the
# panel and discovery cannot see it, so a release missing any is a release
# nobody can use.
for e in nin nsq nhr ndb; do
  [ -f "ribopilot_rrna.$e" ] || {
    echo "no ribopilot_rrna.$e in $SRC; the panel is incomplete, rebuild it" >&2; exit 1; }
done

echo "packaging panel ..."
tar -czf "$OUT/ribopilot_rrna-${VER}.tar.gz" $PANEL_FILES

DECOYS=""
HAVE_TAXDB=""
for n in 16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb; do
  for f in "$n".*; do
    if [ -e "$f" ]; then
      DECOYS="$DECOYS $f"
      if [ "$n" = taxdb ]; then HAVE_TAXDB=1; fi
    fi
  done
done
# install_ribopilot_db() fetches the decoys by default, so a panel-only release
# 404s for every default install after the panel has already been downloaded.
[ -n "$DECOYS" ] || { echo "no decoy databases in $SRC; run fetch_decoys.sh first" >&2; exit 1; }
[ -n "$HAVE_TAXDB" ] || {
  echo "no taxdb in $SRC; staxids would stop resolving, run fetch_decoys.sh" >&2; exit 1; }

echo "packaging decoys ..."
tar -czf "$OUT/decoys-${VER}.tar.gz" $DECOYS

cd "$OUT"
for f in *.tar.gz; do
  sha256sum "$f" > "$f.sha256"
done
ls -lh "$OUT"

if [ "$PUBLISH" != "--publish" ]; then
  echo "assets staged in $OUT; re-run with --publish to create the release"
  exit 0
fi

command -v gh >/dev/null || { echo "gh not on PATH" >&2; exit 1; }
# The db-${VER} tag must not already exist: gh release create cannot overwrite a
# release, so a bad one has to be deleted before this can be re-run.
gh release create "db-${VER}" ./*.tar.gz ./*.sha256 \
  --repo "$REPO" --title "rDNA panel ${VER}" --latest=false \
  --notes "Metazoan nuclear rDNA BLAST panel ${VER}. Install with RiboPilot::install_ribopilot_db(\"${VER}\")."
