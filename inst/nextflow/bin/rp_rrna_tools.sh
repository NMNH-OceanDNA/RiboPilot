#!/bin/bash
# rRNA anchors (barrnap or pybarrnap), the tuned ITSx run, and the stock ITSx
# check on one contig. Writes unit.barrnap.gff, unit.itsx.*, unit.itsx_check.*
# (the stock run's stderr in unit.itsx_check.err), and unit.tools_status (each
# tool's exit status). Runs all three tools, then exits with the barrnap or
# pybarrnap status; ITSx never changes it.
set -euo pipefail
if [ $# -ne 7 ]; then
    echo "usage: rp_rrna_tools.sh <fasta> <tool> <kingdom> <itsx_taxon> <itsx_opts> <rrna_opts> <cpus>" >&2
    exit 2
fi
fa=$1; tool=${2:-pybarrnap}; kingdom=${3:-euk}; taxon=${4:-all}
# $itsxo and $rrnao stay unquoted below: each holds several flags
itsxo=$5; rrnao=$6; cpus=$7
b=0; i=0; c=0

# Same kingdom translation as R's rrna_kingdom_for(): the Perl barrnap fork
# calls its eukaryote Rfam models "fun" and has no euk/all; pybarrnap has no
# "fun". Passing the stored kingdom through untranslated makes the tool exit
# non-zero and the unit gets no anchors at all.
if [ "$tool" = "barrnap" ]; then
    case "$kingdom" in euk|all) kingdom=fun ;; esac
    barrnap --kingdom "$kingdom" --quiet $rrnao "$fa" > unit.barrnap.gff || b=$?
else
    case "$kingdom" in fun) kingdom=euk ;; esac
    pybarrnap -k "$kingdom" -t "$cpus" -q $rrnao "$fa" > unit.barrnap.gff || b=$?
fi

ITSx -i "$fa" -o unit.itsx --save_regions none --graphical F \
    --cpu "$cpus" -t "$taxon" --preserve T $itsxo || i=$?
[ -f unit.itsx.positions.txt ] || : > unit.itsx.positions.txt

# second, stock run (defaults, no organism group): what NCBI's Submission
# Portal will see. Its output becomes itsx_warnings.
ITSx -i "$fa" -o unit.itsx_check --preserve T --detailed_results T \
    --save_regions none --graphical F --cpu "$cpus" 2> unit.itsx_check.err || c=$?
[ -f unit.itsx_check.positions.txt ] || : > unit.itsx_check.positions.txt

printf 'barrnap=%s\nitsx=%s\nitsx_check=%s\n' "$b" "$i" "$c" > unit.tools_status
exit "$b"
