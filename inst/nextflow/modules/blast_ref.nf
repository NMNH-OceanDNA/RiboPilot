// WF1 reference fetch + align. blast_ref_fetch pulls each candidate accession's
// sequence + feature table from NCBI (edirect efetch, skipping already-fetched
// accessions) into blast_ref_sequences / blast_ref_annotations. blast_ref_align
// then computes per-feature sample-vs-reference nucleotide alignments (pure R).

process blast_ref_fetch {
    cache false
    tag "${id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    val id

    output:
    val id

    script:
    """
    ACCS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT DISTINCT accession FROM blast_ref_candidates WHERE ID = '${id}';")
    if [ -z "\$ACCS" ]; then echo "no candidates for ${id}"; exit 0; fi
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 30 )); fi
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'
    for acc in \$ACCS; do
        have=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT 1 FROM blast_ref_sequences WHERE accession = '\$acc';")
        if [ -n "\$have" ]; then continue; fi
        efetch -db nuccore -id "\$acc" -format fasta > "\$acc.fa" 2>/dev/null || true
        efetch -db nuccore -id "\$acc" -format ft    > "\$acc.ft" 2>/dev/null || true
        if [ ! -s "\$acc.fa" ]; then echo "efetch failed for \$acc"; continue; fi
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_fetch.R '${params.db}' "\$acc" "\$acc.fa" "\$acc.ft"
    done
    """
}

process blast_ref_align {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'

    input:
    tuple val(id), val(upath), val(uscaf)

    output:
    tuple val(id), val(upath), val(uscaf)

    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_align.R '${params.db}' '${id}' ${upath} ${uscaf}
    """
}
