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
command -v jq >/dev/null || { echo "jq not on PATH" >&2; exit 1; }
# Absolute, because the script cd's into SRC before writing anything
OUT="$(cd "${SRC%/}/.." && pwd)/release-${VER}"
REPO="dmacguigan/RiboPilot"

cd "$SRC"

# The version has to be the same string in all three places or the release is
# internally inconsistent: discover_panel() reads the directory basename, the
# installer reads the tarball name, and the app reads manifest.json.
BASENAME="${PWD##*/}"
[ "$BASENAME" = "$VER" ] \
  || { echo "$SRC is named '$BASENAME', not '$VER'; packaging it would ship a tarball whose name disagrees with the panel it holds" >&2; exit 1; }
[ -f manifest.json ] || { echo "no manifest.json in $SRC" >&2; exit 1; }
MANIFEST_VER="$(jq -r '.db_version // empty' manifest.json)"
[ "$MANIFEST_VER" = "$VER" ] \
  || { echo "manifest.json db_version is '$MANIFEST_VER', not '$VER'" >&2; exit 1; }

# update_blastdb.pl keeps the .tar.gz.md5 markers but never the archive itself,
# so a decoy archive still sitting here means a fetch was interrupted
# mid-decompress; packaging it would put a 72 MB tarball inside a tarball.
for f in *.tar.gz; do
  [ -e "$f" ] || break
  echo "$SRC still holds $f; the decoy fetch did not finish, re-run fetch_decoys.sh" >&2
  exit 1
done

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

DECOYS=""
# .n* rather than .*: every nucleotide volume extension begins with n, while
# update_blastdb.pl deliberately preserves a <db>.tar.gz.md5 beside them as the
# cache key for its next run. Shipping that marker would plant a stale key in
# every user's decoy directory, and this task's whole premise is that the marker
# is never evidence.
for n in 16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA; do
  for f in "$n".n*; do
    if [ -e "$f" ]; then DECOYS="$DECOYS $f"; fi
  done
done
# install_ribopilot_db() fetches the decoys by default, so a panel-only release
# 404s for every default install after the panel has already been downloaded.
[ -n "$DECOYS" ] || { echo "no decoy databases in $SRC; run fetch_decoys.sh first" >&2; exit 1; }

# Named one by one, never a taxdb.* glob: taxonomy4blast.sqlite3 does not match
# that glob, so it was silently left out of every release, and without it
# -taxids returns unfiltered rows and blastn still exits 0.
for f in taxdb.btd taxdb.bti taxonomy4blast.sqlite3 DECOY_VERSION.txt; do
  [ -f "$f" ] \
    || { echo "no $f in $SRC; run fetch_decoys.sh" >&2; exit 1; }
  DECOYS="$DECOYS $f"
done

# Only once every refusal above is past: a refused run must leave nothing
# behind, neither a lone panel tarball with no checksum nor an empty
# release-<VER>, which is a version-carrying directory --check-version reports.
mkdir -p "$OUT"
echo "packaging panel ..."
tar -czf "$OUT/ribopilot_rrna-${VER}.tar.gz" $PANEL_FILES
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
  --notes "Metazoan nuclear rDNA BLAST panel ${VER}. Maintainer assets: the panel ships baked into the RiboPilot container image, so end users do not install it themselves."
