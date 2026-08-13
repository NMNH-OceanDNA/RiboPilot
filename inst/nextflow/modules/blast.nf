// WF1 remote GenBank BLAST: megablast each nrDNA scaffold against
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
    RB=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT b.run_blast FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    if [ "\${RB:-1}" != "1" ]; then echo "blast disabled for ${id}"; exit 0; fi

    EQ=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT COALESCE(b.entrez_query,'') FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    XO=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT COALESCE(b.extra_opts,'') FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    MTS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT COALESCE(b.max_target_seqs,50) FROM assemble a JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = '${id}';")
    EQ_OPT=""
    if [ -n "\${EQ}" ]; then EQ_OPT="-entrez_query"; fi

    # Back off on retries to let NCBI recover from rate limits
    if [ "${task.attempt}" -gt 1 ]; then sleep \$(( (${task.attempt} - 1) * 60 )); fi
    export NCBI_API_KEY='${params.ncbi_api_key ?: ""}'

    # BLAST each non-ignored scaffold (unit) separately -> per-unit candidates.
    # `for` (not a pipe) so a blastn failure can `exit 1` and trigger a retry.
    UNITS=\$(sqlite3 -cmd '.timeout 10000' -separator '.' '${params.db}' "SELECT path||'.'||scaffold FROM assemblies WHERE ID='${id}' AND COALESCE(ignore,0)=0;")
    for U in \$UNITS; do
        P=\${U%.*}; S=\${U#*.}
        # Already screened on an earlier attempt: do not re-submit to NCBI. The
        # loop exits 1 on any unit's failure under a shared retry budget, so
        # without this a 30-unit sample re-runs every completed search.
        DONE=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT COUNT(*) FROM blast_hits WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" 2>/dev/null || echo 0)
        if [ "\${DONE:-0}" -gt 0 ]; then continue; fi

        sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM assemblies WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" > query.fa
        [ -s query.fa ] || continue

        # \${XO} carries the significance and identity floors (blast_opts.extra_opts,
        # "Extra blastn options" in the app); blastn's own defaults keep chance hits.
        # Keep the hit distribution, not just one row per subject: bitscore (the
        # statistic that actually orders alignments), alignment length, and the
        # query coordinates that make a chimera visible as two disjoint HSP sets
        # with different subjects. Same single remote round trip as before.
        # slen enables a subject-coverage filter without a second remote re-run;
        # qseqid is captured for a future batched query. Appended last so the
        # positional parser in parse_blast_hits_file() needs no change.
        OUTFMT="6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"
        if blastn -remote -db core_nt -query query.fa \
            -outfmt "\$OUTFMT" \
            -max_target_seqs "\${MTS:-50}" -max_hsps 5 -task megablast \
            \${EQ_OPT:+\$EQ_OPT "\$EQ"} \${XO} > hits.txt; then
            [ -s hits.txt ] || echo "NO_SIGNIFICANT_HITS" > hits.txt
        else
            exit 1
        fi

        # megablast needs a 28 bp exact match, which a divergent fungal or
        # bacterial contaminant against a metazoan query often cannot provide,
        # so "no hits" reads as unchecked rather than as a signal. One bounded
        # discontiguous-megablast retry for those units only.
        BESTCOV=\$(awk -F'\\t' 'BEGIN{m=0} \$4+0>m{m=\$4+0} END{print m+0}' hits.txt 2>/dev/null || echo 0)
        if grep -q NO_SIGNIFICANT_HITS hits.txt || [ "\${BESTCOV%%.*}" -lt 25 ]; then
            if blastn -remote -db core_nt -query query.fa \
                -outfmt "\$OUTFMT" \
                -max_target_seqs "\${MTS:-50}" -max_hsps 5 -task dc-megablast \
                \${EQ_OPT:+\$EQ_OPT "\$EQ"} \${XO} > hits_dc.txt 2>/dev/null; then
                if [ -s hits_dc.txt ]; then
                    echo "dc-megablast fallback used for ${id} \$P.\$S" >&2
                    mv hits_dc.txt hits.txt
                fi
            fi
        fi
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_blast.R '${params.db}' '${id}' "\$P" "\$S" hits.txt '${params.taxdump_dir}'
    done
    """
}
