#!/bin/bash
# Format GetOrganelle path fastas -> <id>_assembly_<n>.fasta with headers
# ">id.<path>.<scaffold> <topology> <original_header>". The original GetOrganelle
# header (scaffold_N--<edge ids>) is kept so downstream scoring can look up
# per-contig graph depth; loaders parse only the leading id.path.scaffold token.
#
# n is a primary key: the curator's keep/drop decisions are stored against it,
# so it must come from the graph/path numbers GetOrganelle puts in the filename
# (...graph<G>.<P>.path_sequence.fasta), not from glob order, which reorders
# once G or P reaches 10 (graph1.10 sorts before graph1.2).
set -euo pipefail
id=$1
godir=$2
shopt -s nullglob
mapfile -t files < <(
    for f in "$godir"/*path_sequence.fasta; do
        key=$(sed -n 's/.*graph\([0-9][0-9]*\)\.\([0-9][0-9]*\)\.path_sequence\.fasta$/\1 \2/p' <<< "${f##*/}")
        printf '%s %s\n' "${key:-999999 999999}" "$f"
    done | sort -k1,1n -k2,2n -k3 | cut -d' ' -f3-
)
if [ ${#files[@]} -eq 0 ]; then
    why=$(grep -hE ' - (ERROR|WARNING): ' "$godir/get_org.log.txt" 2>/dev/null \
        | grep -v 'still trying' \
        | sed -E 's/^.* - (ERROR|WARNING): //; s/system call for: "[^"]*" //' \
        | tail -n 3 | cut -c1-200 | awk 'NR > 1 { printf "; " } { printf "%s", $0 }' || true)
    echo ">No assembly found${why:+: $why}" > "${id}_assembly_0.fasta"
    exit 0
fi
n=0
for f in "${files[@]}"; do
    n=$((n + 1))
    topo="linear"
    if grep -qi "circular\|complete" <<< "$f"; then topo="circular"; fi
    awk -v topo="$topo" -v id="$id" -v p="$n" \
        '/^>/ {orig=substr($0,2); print ">" id "." p "." ++c " " topo " " orig} !/^>/ {print}' \
        "$f" > "${id}_assembly_${n}.fasta"
done
