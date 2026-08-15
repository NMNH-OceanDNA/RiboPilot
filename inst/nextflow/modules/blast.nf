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
    cpus 4

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

    # BLAST each non-ignored scaffold (unit) separately -> per-unit candidates.
    # `for` (not a pipe) so a blastn failure can `exit 1` and trigger a retry.
    UNITS=\$(sqlite3 -cmd '.timeout 10000' -separator '.' '${params.db}' "SELECT path||'.'||scaffold FROM assemblies WHERE ID='${id}' AND COALESCE(ignore,0)=0;")

    DBARG=""
    SRC="remote"
    if [ -n '${params.blastdb}' ]; then
        if DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}'); then
            SRC="local"
        else
            DBARG=""
            echo "configured local panel '${params.blastdb}' could not be used, falling back to remote for ${id}" >&2
        fi
    else
        # empty blastdb is the normal unconfigured case, not an error: call
        # rp_dbarg.R anyway (a later task uses this for panel discovery) but
        # discard its usage-error noise for this argument-less invocation
        if DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '' 2>/dev/null); then
            SRC="local"
        else
            DBARG=""
        fi
    fi
    if [ "\$SRC" = "local" ] && [ -n "\${EQ}" ]; then
        echo "entrez_query is remote-only, ignored for local search on ${id}" >&2
    fi
    if [ "\$SRC" = "local" ] && echo "\$XO" | grep -q -- '-evalue'; then
        EVAL_VAL=\$(echo "\$XO" | awk '{for (i=1;i<=NF;i++) if (\$i=="-evalue") print \$(i+1)}')
        echo "WARNING: extra_opts sets -evalue \${EVAL_VAL} for the local search on ${id}; the local panel (~7x10^8 bp) is a much smaller search space than core_nt (~10^12 bp), so this threshold does not mean the same thing here" >&2
    fi

    for U in \$UNITS; do
        P=\${U%.*}; S=\${U#*.}
        # per-unit copy: a local failure on one unit must not leak "remote"
        # into the plain shell loop's next iteration. SRC stays the task-level
        # resolution of whether a usable panel exists at all.
        UNIT_SRC="\$SRC"
        sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM assemblies WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" > query.fa
        [ -s query.fa ] || continue

        # Already searched with this exact query on an earlier attempt or run: do
        # not re-submit to NCBI. Keyed on a hash of the query sequence rather than
        # a blast_hits row count, so a unit that found nothing is not re-searched
        # forever, and a re-assembled or trimmed unit is searched again instead of
        # inheriting the previous contig's hits.
        QMD5=\$( (grep -v '^>' query.fa || true) | tr -d '[:space:]' | (md5sum 2>/dev/null || md5 -q) | awk '{print \$1}')
        PREV=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "SELECT query_md5 FROM blast_searched WHERE ID='${id}' AND path=\$P AND scaffold=\$S;" 2>/dev/null || echo "")
        if [ "\${PREV}" = "\${QMD5}" ]; then echo "already searched ${id} \$P.\$S"; continue; fi

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
        if [ "\$UNIT_SRC" = "local" ]; then
            if blastn -db "\$DBARG" -query query.fa \
                -outfmt "\$OUTFMT" \
                -max_target_seqs "\${MTS:-50}" -max_hsps 5 -task megablast \
                -num_threads ${task.cpus} \
                \${XO} > hits.txt; then
                [ -s hits.txt ] || echo "NO_SIGNIFICANT_HITS" > hits.txt
            else
                # a panel that exists and passed rp_dbarg.R but was built without
                # -parse_seqids (or mixed with one that was) fails inside blastn
                # itself, not before it. That is a config error, not "no answer":
                # fall back to remote for this unit rather than losing it entirely.
                echo "local BLAST failed for ${id} \$P.\$S (misconfigured panel?), falling back to remote" >&2
                UNIT_SRC="remote"
            fi
        fi
        if [ "\$UNIT_SRC" = "remote" ]; then
            if blastn -remote -db core_nt -query query.fa \
                -outfmt "\$OUTFMT" \
                -max_target_seqs "\${MTS:-50}" -max_hsps 5 -task megablast \
                \${EQ_OPT:+\$EQ_OPT "\$EQ"} \${XO} > hits.txt; then
                [ -s hits.txt ] || echo "NO_SIGNIFICANT_HITS" > hits.txt
            else
                exit 1
            fi
        fi

        # megablast needs a 28 bp exact match, which a divergent fungal or
        # bacterial contaminant against a metazoan query often cannot provide,
        # so "no hits" reads as unchecked rather than as a signal. One bounded
        # discontiguous-megablast retry for those units only.
        BESTCOV=\$(awk -F'\\t' 'BEGIN{m=0} \$4+0>m{m=\$4+0} END{print m+0}' hits.txt 2>/dev/null || echo 0)
        if grep -q NO_SIGNIFICANT_HITS hits.txt || [ "\${BESTCOV%%.*}" -lt 25 ]; then
            if [ "\$UNIT_SRC" = "local" ]; then
                if blastn -db "\$DBARG" -query query.fa \
                    -outfmt "\$OUTFMT" \
                    -max_target_seqs "\${MTS:-50}" -max_hsps 5 -task dc-megablast \
                    -num_threads ${task.cpus} \
                    \${XO} > hits_dc.txt 2>/dev/null; then
                    if [ -s hits_dc.txt ]; then
                        echo "dc-megablast fallback used for ${id} \$P.\$S" >&2
                        mv hits_dc.txt hits.txt
                    fi
                fi
            else
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
        fi
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_blast.R '${params.db}' '${id}' "\$P" "\$S" hits.txt '${params.taxdump_dir}'

        NHITS=\$(grep -c . hits.txt 2>/dev/null || echo 0)
        grep -q NO_SIGNIFICANT_HITS hits.txt && NHITS=0
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_searched.R '${params.db}' '${id}' "\$P" "\$S" "\$QMD5" "\$NHITS" "\$UNIT_SRC"
    done
    """
}
