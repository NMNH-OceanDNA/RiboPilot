// WF1 remote GenBank BLAST: megablast the stitched nrDNA sequence against
// core_nt, then parse the hits into (a) ranked nuclear-rDNA reference
// candidates and (b) a contamination screen. Best-effort: retries on
// connection / rate-limit failures, then 'ignore' so other samples keep
// running.
//
// The search is deliberately UNRESTRICTED by default. Restricting it to
// nuclear-rRNA titles (the old default entrez_query) also restricted what a
// contaminant could look like: a symbiont, bacterial 16S, or mitochondrial
// rRNA assembly could only ever match some distant nuclear rDNA record and so
// looked plausible. The nuclear-rDNA filter now runs client-side over the hit
// list (parse_blast_candidates), which keeps the reference panel identical
// while letting the screen see the whole database.

process blast_genbank {
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
    # Skip samples with run_blast = 0
    RB=\$(sqlite3 '${params.db}' "SELECT b.run_blast FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    if [ "\${RB:-1}" != "1" ]; then echo "blast disabled for ${id}"; exit 0; fi

    EQ=\$(sqlite3 '${params.db}' "SELECT COALESCE(b.entrez_query,'') FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    XO=\$(sqlite3 '${params.db}' "SELECT COALESCE(b.extra_opts,'') FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    MTS=\$(sqlite3 '${params.db}' "SELECT COALESCE(b.max_target_seqs,50) FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    EQ_OPT=""
    if [ -n "\${EQ}" ]; then EQ_OPT="-entrez_query"; fi

    # Back off on retries to let NCBI recover from rate limits
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 60 )); fi
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'

    # BLAST each non-ignored scaffold (unit) separately -> per-unit candidates.
    # `for` (not a pipe) so a blastn failure can `exit 1` and trigger a retry.
    UNITS=\$(sqlite3 -separator '.' '${params.db}' "SELECT path||'.'||scaffold FROM assemblies WHERE ID='${id}' AND COALESCE(ignore,0)=0;")
    for U in \$UNITS; do
        P=\${U%.*}; S=\${U#*.}
        sqlite3 '${params.db}' "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM assemblies WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" > query.fa
        [ -s query.fa ] || continue
        if blastn -remote -db core_nt -query query.fa \
            -outfmt "6 saccver stitle pident qcovs evalue staxids" \
            -max_target_seqs "\${MTS:-50}" -max_hsps 1 -task megablast \
            \${EQ_OPT:+\$EQ_OPT "\$EQ"} \${XO} > hits.txt; then
            [ -s hits.txt ] || echo "NO_SIGNIFICANT_HITS" > hits.txt
        else
            exit 1
        fi
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_blast.R '${params.db}' '${id}' "\$P" "\$S" hits.txt '${params.taxdump_dir}'
    done
    """
}
