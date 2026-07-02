#!/bin/bash
# Run the RiboPilot assembly pipeline (fastp + GetOrganelle nr) on test samples.
#
#   conda activate ribopilot
#   bash dev/run_test_assembly.sh [DATADIR] [OUTDIR]
#
# DATADIR  directory of paired reads named <ID>_R1.fastq.gz / <ID>_R2.fastq.gz
#          (default: ./testdata_full). Get the fish test reads with, in R:
#              RiboPilot::new_test_project("some_dir", n = 2, full_size = TRUE)
#          or download from ENA; then point DATADIR at that data folder.
# OUTDIR   where the samplesheet, work dir, and out/ go (default: ./test_run).
#
# Uses the shipped Actinopterygii seed/label DBs (local paths, so no network).
# Results: <OUTDIR>/out/<ID>/assemble/<ID>_assembly_<n>.fasta (+ .gfa, log).
set -euo pipefail

RP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATADIR="${1:-$RP/testdata_full}"
OUTDIR="${2:-$RP/test_run}"

SEEDS="$RP/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta"
LABELS="$RP/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta"
FASTP="--trim_poly_g --correction --detect_adapter_for_pe --dont_eval_duplication"
GETORG="-F anonym -R 10 -k 21,45,65,85,105,115 --larger-auto-ws"
CPUS=8
MEM=24

for t in nextflow fastp get_organelle_from_reads.py; do
    command -v "$t" >/dev/null || { echo "ERROR: '$t' not on PATH (activate the ribopilot conda env)." >&2; exit 1; }
done

mkdir -p "$OUTDIR"
SS="$OUTDIR/samplesheet.tsv"
# TAB-separated so the comma-containing -k list needs no quoting
printf 'ID\tR1\tR2\tcpus\tmemory\tfastp\tgetOrganelle\tseeds_db\tlabels_db\n' > "$SS"
shopt -s nullglob
n=0
for r1 in "$DATADIR"/*_R1.fastq.gz; do
    id=$(basename "$r1" _R1.fastq.gz)
    r2="$DATADIR/${id}_R2.fastq.gz"
    [ -f "$r2" ] || continue
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$id" "$(basename "$r1")" "$(basename "$r2")" "$CPUS" "$MEM" \
        "$FASTP" "$GETORG" "$SEEDS" "$LABELS" >> "$SS"
    n=$((n + 1))
done

if [ "$n" -eq 0 ]; then
    echo "No <ID>_R1.fastq.gz / <ID>_R2.fastq.gz pairs found in: $DATADIR" >&2
    echo "Fetch reads first (see the header of this script)." >&2
    exit 1
fi
echo "Samplesheet: $SS  (${n} sample(s))"

# Minimal token-free config (rawDir/minDepth come from the CLI here; the tokens
# in inst/config.conda are for the new_project()-generated project .config).
CFG="$OUTDIR/nextflow.config"
cat > "$CFG" <<'EOF'
process { executor = 'local' }
params { publishDir = 'out' }
EOF

cd "$OUTDIR"
nextflow run "$RP/inst/nextflow/main.nf" -c "$CFG" \
    -work-dir "$OUTDIR/work" -resume \
    --samplesheet "$SS" --rawDir "$DATADIR" --minDepth 500 --publishDir out

echo
echo "Done. Assemblies: $OUTDIR/out/<ID>/assemble/<ID>_assembly_<n>.fasta"
echo "Load them into a project database from R with:"
echo "  con <- DBI::dbConnect(RSQLite::SQLite(), '<project>/.sqlite')"
echo "  RiboPilot::load_assemblies(con, '$OUTDIR/out', '<ID>')"
