#!/usr/bin/env bash
# Fetch the NCBI decoy databases into a panel directory.
#
# The decoys sit beside the panel because rp_dbarg.R defaults the decoy
# directory to dirname(panel_prefix). taxdb comes along so staxids resolve.
set -ueo pipefail

DEST="${1:?usage: fetch_decoys.sh <panel_dir>}"
DBS="16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb"

command -v update_blastdb.pl >/dev/null || { echo "update_blastdb.pl not on PATH" >&2; exit 1; }
mkdir -p "$DEST"
cd "$DEST"

for d in $DBS; do
  echo "fetching $d ..."
  update_blastdb.pl --decompress "$d"
done

# update_blastdb.pl verifies each archive's md5 as it goes; record what landed
: > DECOY_VERSION.txt
for d in $DBS; do
  printf '%s\t%s\n' "$d" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> DECOY_VERSION.txt
done
echo "decoys installed in $DEST"
