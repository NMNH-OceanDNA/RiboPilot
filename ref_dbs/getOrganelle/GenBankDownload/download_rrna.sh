#!/bin/bash
# Build a GetOrganelle nuclear ribosomal (nrDNA) reference database for a clade.
#
# Downloads nuclear ribosomal RNA records (18S / 5.8S / 28S / ITS) from GenBank
# for a taxonomic query and builds two FASTA files for GetOrganelle:
#   <prefix>_seeds.fasta   -- read-recruiting sequences (18S/28S/whole cistron),
#                             headers ">ACCESSION_Genus_species"
#   <prefix>_labels.fasta  -- per-gene sequences for contig classification,
#                             headers ">18S" / ">5.8S" / ">28S"
#
# Adapted from GenBankMitoReview (github.com/dmacguigan/GenBankMitoReview),
# retargeted from mitochondrial genomes to the nuclear ribosomal cistron.
#
# author: Dan MacGuigan
#
# Run inside a conda env with Entrez Direct + biopython (e.g. 'mitoreview' or
# the RiboPilot env). The bundled bin/curl shim forces HTTP/1.1 to avoid the
# "curl (56) unexpected eof" failure Entrez Direct hits on OpenSSL 3.x.
#
# The query is an advanced GenBank search term, e.g. '"Actinopterygii"[Organism]'.
#
# Usage:
#   bash download_rrna.sh [--outdir DIR] [--api-key KEY] [--refseq] '<query>' <prefix>
#
# Example:
#   bash download_rrna.sh --outdir out '"Actinopterygii"[Organism]' actinopterygii_rrna

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH="${SCRIPT_DIR}/bin:${PATH}"

print_help() {
  cat <<'EOF'
Usage:
  bash download_rrna.sh [OPTIONS] '<query>' <prefix>

Description:
  Downloads nuclear ribosomal RNA records (18S / 5.8S / 28S / ITS) from NCBI
  GenBank for a taxonomic query and builds GetOrganelle nr seed + label FASTAs.

Arguments:
  '<query>'   Advanced GenBank search term (single-quoted), e.g.
                '"Actinopterygii"[Organism]'
  <prefix>    Output file prefix (e.g. 'actinopterygii_rrna')

Options:
  --outdir <dir>   Output directory (default: current directory).
  --api-key <key>  NCBI API key (raises rate limit 3 -> 10 req/s).
  --refseq         Restrict to RefSeq records (NR_*/NG_* etc).
  --count-only     Print the matching record count and exit (no download).
  -h, --help       Show this help and exit.

Output:
  <outdir>/<prefix>.raw.fasta     all matching records (fasta)
  <outdir>/<prefix>_seeds.fasta   read-recruiting sequences for GetOrganelle
  <outdir>/<prefix>_labels.fasta  per-gene label sequences for GetOrganelle

Dependencies (conda): Entrez Direct (esearch, efilter, efetch), python, biopython
EOF
}

# nuclear rRNA element terms (nuclear-specific: excludes mito 12S/16S)
RRNA_TERMS='18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal RNA[Title] OR internal transcribed spacer[Title]'

COMPLETE=0
REFSEQ=0
COUNT_ONLY=0
OUTDIR="."
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)    print_help; exit 0 ;;
    --refseq)     REFSEQ=1; shift ;;
    --count-only) COUNT_ONLY=1; shift ;;
    --outdir)     [[ $# -ge 2 ]] || { echo "ERROR: --outdir requires a value" >&2; exit 1; }
                  OUTDIR="$2"; shift 2 ;;
    --api-key)    [[ $# -ge 2 ]] || { echo "ERROR: --api-key requires a value" >&2; exit 1; }
                  export NCBI_API_KEY="$2"; shift 2 ;;
    --) shift; break ;;
    -*) echo "ERROR: unknown option '$1'" >&2; echo "Use -h for help." >&2; exit 1 ;;
    *)  break ;;
  esac
done

if [[ $# -ne 2 ]]; then
  echo "Usage: bash download_rrna.sh [--outdir DIR] [--api-key KEY] [--refseq] [--count-only] '<query>' <prefix>" >&2
  echo "Use -h for full help." >&2
  exit 1
fi

USER_QUERY="$1"
PREFIX="$2"

mkdir -p "${OUTDIR}"
OUTDIR="${OUTDIR%/}"
FULL_PREFIX="${OUTDIR}/${PREFIX}"

# Build the GenBank query: taxon AND nuclear rRNA titles, exclude organelle rRNA.
QUERY="(${USER_QUERY}) AND (${RRNA_TERMS}) NOT mitochondrion[Filter] NOT plastid[Filter]"
if [[ "${REFSEQ}" -eq 1 ]]; then
  QUERY="${QUERY} AND refseq[filter]"
fi

echo "GenBank query:"
echo "  ${QUERY}"
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  echo "  (NCBI API key set -- 10 req/s)"
fi
echo ""

N_REC=$( esearch -db nuccore -query "${QUERY}" \
  | efilter -molecule genomic \
  | grep "Count" | head -1 | cut -f2 -d">" | cut -f1 -d"<" )

echo "Found ${N_REC} matching nuclear rRNA record(s)."
if [[ "${COUNT_ONLY}" -eq 1 ]]; then
  exit 0
fi
echo "Downloading (fasta)..."
echo ""

esearch -db nuccore -query "${QUERY}" \
  | efilter -molecule genomic \
  | efetch -format fasta > "${FULL_PREFIX}.raw.fasta"

echo "Raw records written to ${FULL_PREFIX}.raw.fasta"
echo "Building GetOrganelle seed + label FASTAs..."
echo ""

python "${SCRIPT_DIR}/parse_rrna_fasta.py" \
  "${FULL_PREFIX}.raw.fasta" \
  --seeds  "${FULL_PREFIX}_seeds.fasta" \
  --labels "${FULL_PREFIX}_labels.fasta"

echo ""
echo "Done."
echo "  seeds:  ${FULL_PREFIX}_seeds.fasta"
echo "  labels: ${FULL_PREFIX}_labels.fasta"
