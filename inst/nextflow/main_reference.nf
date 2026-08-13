// Reference clade databases: count GenBank records, or build GetOrganelle seed
// and label FASTAs. Driven by a work-list TSV written by the app
// (write_clade_manifest): ref_id, taxid, clade_name, query.
//
// Two entry points rather than two files: the processes share their inputs,
// their container, and their retry policy, and a second file drifts out of step
// (see dev/design/15_WF2_DECOMMISSION.md).

params.clades = ''
params.refDir = ''
params.outDir = ''
params.ribopilot_src = ''
params.ncbi_api_key = ''
params.refseq = false

process count_clade {
    cache false
    tag "${clade_name}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir params.outDir, mode: 'copy'

    input:
    tuple val(taxid), val(clade_name), val(query)

    output:
    path "${taxid}.count"

    script:
    """
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 20 )); fi
    export RP_REFSEQ='${params.refseq}'
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_db.R count '${taxid}' '${query}' '${taxid}.count'
    """
}

process build_ref_db {
    cache false
    tag "${ref_id}"
    errorStrategy 'ignore'
    // the wrapper writes the seeds/ and labels/ layout, so publishing the
    // relative paths unchanged lands them where reference_db_paths() looks
    publishDir params.refDir, mode: 'copy', pattern: 'seeds/*_seeds.fasta'
    publishDir params.refDir, mode: 'copy', pattern: 'labels/*_labels.fasta'

    input:
    tuple val(ref_id), val(query)

    output:
    path "seeds/*_seeds.fasta", optional: true
    path "labels/*_labels.fasta", optional: true

    script:
    """
    export RP_REFSEQ='${params.refseq}'
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_db.R build '${ref_id}' '${query}' .
    """
}

workflow COUNT {
    Channel.fromPath(params.clades)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.taxid, r.clade_name, r.query) }
        | count_clade
}

workflow BUILD {
    Channel.fromPath(params.clades)
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.ref_id, r.query) }
        | build_ref_db
}
