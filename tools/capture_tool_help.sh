#!/bin/bash
# Regenerate inst/tool_help/<tool>.txt from the RiboPilot Docker image, so the
# app's CLI-help popovers match the tool versions the pipeline actually runs.
# Usage: tools/capture_tool_help.sh [image]   (default: ribopilot:<DESCRIPTION Version>)
set -euo pipefail

image="${1:-ribopilot:$(awk '/^Version:/{print $2}' DESCRIPTION)}"
here="$(cd "$(dirname "$0")/.." && pwd)"
out="$here/inst/tool_help"
mkdir -p "$out"

# tool name -> help command inside the image
declare -A CMD=(
  [fastp]="fastp --help"
  [getOrganelle]="get_organelle_from_reads.py --help"
  [blastn]="blastn -help"
  [pybarrnap]="pybarrnap --help"
  [barrnap]="barrnap --help"
  [ITSx]="ITSx --help"
  [bowtie2]="bowtie2 --help"
  [spades]="spades.py --help"
)

for t in "${!CMD[@]}"; do
  echo "capturing $t ..."
  docker run --rm "$image" bash -lc "${CMD[$t]}" > "$out/$t.txt" 2>&1 || true
done
echo "wrote help text to $out"
