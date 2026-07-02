#!/bin/bash
# Format GetOrganelle path fastas -> <id>_assembly_<n>.fasta with headers
# ">id.<path>.<scaffold> <topology>" (topology inferred from the filename).
set -euo pipefail
id=$1
godir=$2
shopt -s nullglob
files=("$godir"/*path_sequence.fasta)
if [ ${#files[@]} -eq 0 ]; then
    echo ">No assembly found" > "${id}_assembly_0.fasta"
    exit 0
fi
n=0
for f in "${files[@]}"; do
    n=$((n + 1))
    topo="linear"
    if grep -qi "circular\|complete" <<< "$f"; then topo="circular"; fi
    awk -v topo="$topo" -v id="$id" -v p="$n" \
        '/^>/ {print ">" id "." p "." ++c " " topo} !/^>/ {print}' \
        "$f" > "${id}_assembly_${n}.fasta"
done
