// Processes of the recheck workflow (main_recheck.nf): rerun every offline
// check on the review copies queued on Annotation Review. All are uncached and,
// after the query, run per contig. dev/specs/2026-10-08-contig-recheck-design.md,
// sections 5.4 and 7.

// Declared here, not only in main_recheck.nf: includes are resolved before the
// including script's body runs, so a main_recheck.nf default is still null when
// this module reads it, and a null interpolates into the shell as 'null'.
params.blastdb = ''
params.taxdump_dir = ''

// Take the queued copies and mark them running (4) in one transaction, with
// the rev each copy had at that moment. An edit after this raises the rev, so
// recheck_tools and the recheck discard the result. The selection is
// RECHECK_QUEUE_SQL in R/recheck.R, word for word (test-nextflow.R checks it);
// nf_recheck_pending() counts the same rows.
process query_recheck {
    cache false
    executor 'local'
    // A read that loses a lock race must not terminate the run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3
    output:
    path 'recheck_units.tsv'
    script:
    """
    printf 'ID\tpath\tscaffold\trev\tcpus\tmem\tmode\n' > recheck_units.tsv
    sqlite3 -cmd '.timeout 10000' -separator '\t' '${params.db}' "
      BEGIN IMMEDIATE;
      CREATE TEMP TABLE q AS
      SELECT r.ID, r.path, r.scaffold, CAST(r.rev AS TEXT) AS rev,
          COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem,
          CASE WHEN COALESCE(bo.run_blast,1) = 1 AND EXISTS (SELECT 1 FROM blast_searched bs
            WHERE bs.ID = r.ID AND bs.path = r.path AND bs.scaffold = r.scaffold AND bs.source = 'local')
          THEN 'local' ELSE 'map' END AS mode
        FROM review_units r
        JOIN annotate an ON an.ID = r.ID AND an.path = r.path AND an.scaffold = r.scaffold
        JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
        JOIN assemblies asm ON asm.ID = r.ID AND asm.path = r.path AND asm.scaffold = r.scaffold
        JOIN assemble a ON a.ID = r.ID
        JOIN samples s ON s.ID = r.ID
        LEFT JOIN blast_opts bo ON bo.blast_opts = a.blast_opts
        WHERE r.check_switch IN (1,4) AND COALESCE(an.annotate_lock,0) = 0
          AND COALESCE(asm.ignore,0) = 0 AND a.assemble_lock = 1 AND s.taxonomy_lock = 1;
      SELECT * FROM q;
      UPDATE review_units SET check_switch = 4, time_stamp = CAST(strftime('%s','now') AS INTEGER)
        WHERE EXISTS (SELECT 1 FROM q WHERE q.ID = review_units.ID AND q.path = review_units.path
          AND q.scaffold = review_units.scaffold);
      COMMIT;" >> recheck_units.tsv
    """
}

// barrnap/pybarrnap, tuned ITSx, and the stock ITSx check on one copy, through
// the script annotate_tools uses, with the same options read. The fasta is
// written only while the copy still has the queued rev; otherwise nothing runs
// and nothing is emitted. A tool failure is recorded in unit.tools_status for
// the recheck; only an out-of-memory exit fails the task, so nextflow.config
// retries it with more memory.
process recheck_tools {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    cpus { nthreads as int }
    // scales on retry; see the OOM retry policy in nextflow.config
    memory { "${(mem as int) * task.attempt} GB" }
    input:
    tuple val(id), val(upath), val(uscaf), val(rev), val(nthreads), val(mem), val(mode)
    output:
    tuple val(id), val(upath), val(uscaf), val(rev), val(mode), path('unit.tools_status'),
          path('unit.barrnap.gff'), path('unit.itsx.positions.txt'), path('unit.itsx_check*'), val(nthreads), optional: true
    script:
    """
    sqlite3 -cmd '.timeout 10000' '${params.db}' \
      "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM review_units
       WHERE ID='${id}' AND path=${upath} AND scaffold=${uscaf} AND rev=${rev}
         AND check_switch IN (1,4) AND sequence IS NOT NULL;" > unit.fasta

    read TOOL KINGDOM ITSXT ITSXO <<< "\$(sqlite3 -cmd '.timeout 10000' -separator ' ' '${params.db}' "
      SELECT COALESCE(ao.rrna_tool,'pybarrnap'), COALESCE(ao.rrna_kingdom,'euk'),
             COALESCE(ao.itsx_taxon,'all'), COALESCE(ao.itsx_opts,'')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"
    TOOL="\${TOOL:-pybarrnap}"; KINGDOM="\${KINGDOM:-euk}"; ITSXT="\${ITSXT:-all}"
    # own query: an empty itsx_opts would otherwise shift the columns under read
    RRNAO="\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "
      SELECT COALESCE(CASE WHEN ao.rrna_tool = 'barrnap' THEN ao.barrnap_opts ELSE ao.pybarrnap_opts END, '')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"

    rc=0
    if [ -s unit.fasta ]; then
        rp_rrna_tools.sh unit.fasta "\$TOOL" "\$KINGDOM" "\$ITSXT" "\$ITSXO" "\$RRNAO" ${task.cpus} || rc=\$?
    fi
    case "\$rc" in 104|134|137|139|140|143|247) exit "\$rc" ;; esac
    """
}

// Score the copy: recheck_unit() reruns the checks on the review tables and
// marks the copy done only if its rev is unchanged. The local panel and the
// labels database are resolved as the main pipeline's chimera process does.
// rp_recheck.R exits 0 for done, stale, and failed; the database holds which.
process recheck {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    cpus { nthreads as int }
    input:
    tuple val(id), val(upath), val(uscaf), val(rev), val(mode), path(status), path(gff), path(pos), path(chk), val(nthreads)
    script:
    """
    DBARG=""
    if [ -z "\${RIBOPILOT_NO_LOCAL_BLAST:-}" ]; then
        DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}' 2>/dev/null || echo "")
    fi
    LABELS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "
      SELECT COALESCE(ao.labels_db,'') FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      WHERE a.ID='${id}' LIMIT 1;" 2>/dev/null || echo "")
    RP_THREADS=${task.cpus} RIBOPILOT_SRC='${params.ribopilot_src}' rp_recheck.R \
      '${params.db}' '${id}' ${upath} ${uscaf} ${rev} ${mode} . "\$DBARG" '${params.taxdump_dir}' "\$LABELS"
    """
}
