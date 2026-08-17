#!/usr/bin/env bash
# Fetch, verify and version the NCBI decoy databases into a panel directory.
#
# The decoys sit beside the panel because local_blast_db_arg() (R/blast_local.R)
# defaults the decoy directory to the panel's own directory. taxdb comes along
# so staxids resolve.
set -Eeuo pipefail
export LC_ALL=C
trap 'echo "fetch_decoys.sh: error at line $LINENO" >&2' ERR

DEST="${1:?usage: fetch_decoys.sh <panel_dir>}"
# taxdb LAST: every decoy tarball bundles its own 2024 taxdb, so fetching taxdb
# first would have it overwritten by the next decoy and the panel would ship a
# two-year-old taxonomy.
SEARCHABLE_DBS="16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA"
DBS="$SEARCHABLE_DBS taxdb"
# PANEL_REQUIRED_EXT in R/blast_local.R. A decoy the installer would reject must
# not pass this stage either.
REQUIRED_EXT="nin nsq nhr ndb"
# taxonomy4blast.sqlite3 is not a taxdb.* file and is easy to lose. Without it
# -taxids returns unfiltered rows and still exits 0: wrong answers, no error.
TAXONOMY_FILES="taxdb.btd taxdb.bti taxonomy4blast.sqlite3"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/blast/db"
METADATA_URL="$BASE_URL/blastdb-metadata-1-1.json"

for t in update_blastdb.pl blastdbcmd curl jq; do
  command -v "$t" >/dev/null || { echo "$t not on PATH" >&2; exit 1; }
done
mkdir -p "$DEST"
cd "$DEST"

for d in $DBS; do
  echo "fetching $d ..."
  # --source ncbi skips ~30 s of AWS/GCP metadata probes per database and stops
  # update_blastdb.pl silently switching mirrors between runs. --decompress is
  # documented as applying only when the source is ncbi.
  update_blastdb.pl --source ncbi --decompress "$d"
done

# Everything below runs AFTER the fetch and looks only at the installed volumes.
# The <db>.tar.gz.md5 left on disk is update_blastdb.pl's cache key, not
# evidence: with that marker present and the volumes deleted it prints "up to
# date", installs nothing and exits 0.
for d in $SEARCHABLE_DBS; do
  for e in $REQUIRED_EXT; do
    [ -s "$d.$e" ] \
      || { echo "$DEST/$d.$e is missing or empty after fetching $d" >&2; exit 1; }
  done
  blastdbcmd -db "$d" -info >/dev/null \
    || { echo "blastdbcmd -info failed for $DEST/$d" >&2; exit 1; }
done
for f in $TAXONOMY_FILES; do
  [ -s "$f" ] \
    || { echo "$DEST/$f is missing or empty; without it -taxids returns unfiltered rows and still exits 0" >&2; exit 1; }
done

# Real versions, all three from the source rather than from disk. Assembled in a
# temp file and moved into place, so an interrupted capture cannot leave a
# half-written DECOY_VERSION.txt for package_release.sh to ship.
echo "recording decoy versions ..."
meta="$(mktemp)"
tmp="$(mktemp "DECOY_VERSION.txt.XXXXXX")"
trap 'rm -f "$meta" "$tmp"' EXIT
curl --fail -sL "$METADATA_URL" -o "$meta" \
  || { echo "could not fetch $METADATA_URL" >&2; exit 1; }
printf '# database\tsource_md5\tlast_modified\tlast_updated\tnumber_of_sequences\n' > "$tmp"
for d in $DBS; do
  # Each capture falls back to an empty string rather than letting errexit take
  # the script out, so the abort below names the source that failed.
  md5="$(curl --fail -sL "$BASE_URL/$d.tar.gz.md5" | awk 'NR==1 {print $1}')" || md5=""
  [ -n "$md5" ] || { echo "could not fetch $BASE_URL/$d.tar.gz.md5" >&2; exit 1; }
  modified="$(curl --fail -sI "$BASE_URL/$d.tar.gz" \
    | tr -d '\r' | sed -n 's/^[Ll]ast-[Mm]odified: //p')" || modified=""
  [ -n "$modified" ] || { echo "no Last-Modified for $BASE_URL/$d.tar.gz" >&2; exit 1; }
  updated="$(jq -r --arg d "$d" '.[] | select(.dbname == $d) | .["last-updated"]' "$meta")" || updated=""
  nseq="$(jq -r --arg d "$d" '.[] | select(.dbname == $d) | .["number-of-sequences"]' "$meta")" || nseq=""
  [ -n "$updated" ] && [ -n "$nseq" ] \
    || { echo "$d has no entry in $METADATA_URL" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\t%s\n' "$d" "$md5" "$modified" "$updated" "$nseq" >> "$tmp"
done
# mktemp makes it 0600, and this one ships inside the decoy tarball.
chmod 0644 "$tmp"
mv "$tmp" DECOY_VERSION.txt
echo "decoys installed and verified in $DEST"
