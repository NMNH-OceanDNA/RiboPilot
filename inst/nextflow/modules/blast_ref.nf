// WF1 reference fetch + align. blast_ref_fetch batch-fetches sequence only
// (edirect efetch, skipping already-fetched accessions) and stores every
// record with locally-derived annotation into blast_ref_sequences /
// blast_ref_annotations. blast_ref_align then computes per-feature
// sample-vs-reference nucleotide alignments (pure R).

// Declared here, not only in main.nf: includes are resolved before the
// including script's body runs, so a main.nf default is still null when this
// module reads it, and a null interpolates into the shell as the word 'null'.
params.ref_dir = ''

process blast_ref_fetch {
    cache false
    tag "${id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    // talks to NCBI; unthrottled, the host gets rate-limited
    maxForks 3
    cpus 4

    input:
    val id

    output:
    val id

    script:
    """
    if [ -n '${params.ref_dir}' ]; then
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_prime.R '${params.db}' '${params.ref_dir}' '${id}' || true
    fi
    ACCS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT DISTINCT c.accession FROM blast_ref_candidates c LEFT JOIN blast_ref_sequences s ON s.accession = c.accession WHERE c.ID = '${id}' AND s.accession IS NULL;")
    if [ -z "\$ACCS" ]; then echo "nothing to fetch for ${id}"; exit 0; fi
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 30 )); fi
    # edirect throttles within one invocation, not across a loop of them; the
    # loop below makes one call per chunk of 200 accessions. NCBI allows 10
    # requests/s with a key and 3 without, so pace the loop to stay under the
    # lower bound either way.
    if [ -n "\${NCBI_API_KEY:-}" ]; then PACE=0.2; else PACE=0.7; fi

    printf '%s\\n' \$ACCS > wanted.txt
    WANTED=\$(grep -c . wanted.txt)

    # one request per chunk instead of two per accession
    : > refs.fa
    split -l 200 -d -a 4 wanted.txt chunk_
    for f in chunk_*; do
        IDS=\$(paste -sd, "\$f")
        if ! efetch -db nuccore -id "\$IDS" -format fasta >> refs.fa; then
            echo "efetch failed for chunk \$f" >&2
        fi
        sleep \$PACE
    done

    if [ ! -s refs.fa ]; then
        echo "no reference sequence returned for ${id}" >&2
        exit 1
    fi
    grep -q -e '<ERROR>' -e '<eFetchResult>' -e 'F a i l e d' -e '^Error' refs.fa && { echo "efetch error body for ${id}" >&2; exit 1; }

    RIBOPILOT_SRC='${params.ribopilot_src}' rp_ref_store.R '${params.db}' refs.fa ${task.cpus} '${params.ref_dir}'

    GOT=\$(grep -c '^>' refs.fa || true)
    echo "reference fetch for ${id}: returned \${GOT} of \${WANTED} requested accession(s)"
    if [ "\$GOT" -lt "\$WANTED" ]; then
        echo "reference fetch for ${id}: got fewer records (\${GOT}) than requested (\${WANTED})" >&2
    fi
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
