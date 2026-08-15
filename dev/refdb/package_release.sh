#!/usr/bin/env bash
# Package a built panel (and its decoys) as GitHub release assets.
#
# Ships the BLAST volume set and the manifest only. The .fasta, .taxid.map and
# .dropped.tsv in the build directory are build artifacts, not query-time files.
set -ueo pipefail

SRC="${1:?usage: package_release.sh <panel_dir> <version> [--publish]}"
VER="${2:?usage: package_release.sh <panel_dir> <version> [--publish]}"
PUBLISH="${3:-}"
OUT="${SRC%/}/../release-${VER}"
REPO="dmacguigan/RiboPilot"

mkdir -p "$OUT"
cd "$SRC"

PANEL_FILES="manifest.json"
for e in ndb nhr nin njs nog nos not nsq ntf nto; do
  [ -f "ribopilot_rrna.$e" ] && PANEL_FILES="$PANEL_FILES ribopilot_rrna.$e"
done

echo "packaging panel ..."
tar -czf "$OUT/ribopilot_rrna-${VER}.tar.gz" $PANEL_FILES

DECOYS=""
for n in 16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb; do
  for f in "$n".*; do [ -e "$f" ] && DECOYS="$DECOYS $f"; done
done
if [ -n "$DECOYS" ]; then
  echo "packaging decoys ..."
  tar -czf "$OUT/decoys-${VER}.tar.gz" $DECOYS
fi

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
# --clobber is deliberately absent: a published database version is immutable
gh release create "db-${VER}" ./*.tar.gz ./*.sha256 \
  --repo "$REPO" --title "rDNA panel ${VER}" --latest=false \
  --notes "Metazoan nuclear rDNA BLAST panel ${VER}. Install with RiboPilot::install_ribopilot_db(\"${VER}\")."
