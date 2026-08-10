// WF1 downstream processes: read/write the project SQLite DB directly (sqlite3
// CLI for the sample query; R-script processes for scoring and annotation
// reconciliation). No R around Nextflow, no nf-sqldb plugin.

// Query the samples ready to assemble straight from the DB (replaces the R
// samplesheet prep). Header + tab separated for splitCsv.
process query_samples {
    cache false
    executor 'local'
    // A read that loses a lock race must not terminate the run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3
    output:
    path 'samplesheet.tsv'
    script:
    """
    # Header is written explicitly (sqlite3 -header emits nothing for 0 rows, and
    # splitCsv(header:true) then fails); an empty result becomes a header-only
    # file so the workflow simply does nothing instead of crashing.
    printf 'ID\tR1\tR2\tpre_cpus\tpre_mem\tasm_cpus\tasm_mem\tfastp\tgetOrganelle\tseeds_db\tlabels_db\tassembler\n' > samplesheet.tsv
    sqlite3 -cmd '.timeout 10000' -separator '\t' '${params.db}' "
      SELECT p.ID, p.R1, p.R2,
             COALESCE(po.cpus,4) AS pre_cpus, COALESCE(po.memory,8) AS pre_mem,
             COALESCE(ao.cpus,4) AS asm_cpus, COALESCE(ao.memory,16) AS asm_mem,
             po.fastp AS fastp,
             ao.getOrganelle || ' --max-paths-num ' || COALESCE(ao.max_paths, 10) AS getOrganelle,
             ao.seeds_db, ao.labels_db,
             COALESCE(ao.assembler, 'getorganelle') AS assembler
      FROM assemble a
      JOIN preprocess p     ON a.ID = p.ID
      JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      JOIN pre_opts po      ON p.pre_opts = po.pre_opts
      WHERE a.assemble_switch IN (1,4) AND a.assemble_lock = 0;" >> samplesheet.tsv
    """
}

// WF2: assembled samples queued for (re-)annotation, with their annotate_opts
// compute settings.
// Units queued for annotation ONLY: a trim, an annotate_opts or curate_opts
// change, a restored scaffold, an edited sequence. WF2 used to own this; there
// is one workflow now (dev/design/15_WF2_DECOMMISSION.md).
//
// Samples this run is about to assemble are EXCLUDED. Assembly re-seeds a
// sample's annotate units itself, so without the exclusion a sample queued for
// both would be annotated twice in one run, racing on the same rows.
//
// WF2's assemble_lock = 1 gate is deliberately NOT carried over: it existed so
// re-annotation could not run before the assembly had been reviewed, and with
// one workflow a unit is annotate-only work because a user queued it.
process query_annotate_units {
    cache false
    executor 'local'
    // A read that loses a lock race must not terminate the run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3
    output:
    path 'annotate_units.tsv'
    script:
    """
    printf 'ID\tpath\tscaffold\tcpus\tmem\n' > annotate_units.tsv
    sqlite3 -cmd '.timeout 10000' -separator '\t' '${params.db}' "
      SELECT an.ID, an.path, an.scaffold, COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem
      FROM annotate an
      JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
      JOIN assemble a ON a.ID = an.ID
      WHERE an.annotate_switch IN (1,4) AND an.annotate_lock = 0
        AND COALESCE(asm.ignore,0) = 0
        AND NOT (a.assemble_switch IN (1,4) AND COALESCE(a.assemble_lock,0) = 0);" >> annotate_units.tsv
    """
}

// WF1: the assembly units of ONE sample that has just been assembled, so
// annotation can run over every scaffold as it is produced rather than waiting
// for the user to review and lock. Unlike query_annotate_units this does not
// gate on annotate_switch: nothing has been reviewed yet.
process query_sample_units {
    cache false
    executor 'local'
    // A read that loses a lock race must not terminate the run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3
    tag "${id}"
    input:
    val id
    output:
    tuple val(id), path('units.tsv')
    script:
    """
    printf 'ID\tpath\tscaffold\tcpus\tmem\n' > units.tsv
    sqlite3 -cmd '.timeout 10000' -separator '\t' '${params.db}' "
      SELECT an.ID, an.path, an.scaffold, COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem
      FROM annotate an
      JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
      WHERE an.ID = '${id}' AND an.annotate_lock = 0 AND COALESCE(asm.ignore,0) = 0;" >> units.tsv
    """
}

// Load assemblies + score candidate paths (auto-select). Reads the published
// assembly fastas from the project out/ dir; the staged inputs force ordering
// after assemble.
process score {
    cache false
    tag "${id}"
    // Retry first: the failure this is most likely to hit is DB lock contention
    // with another task, which clears on its own.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    input:
    tuple val(id), path(fastas)
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_score.R '${params.db}' '${id}' \
        '${params.projectDir}/${params.publishDir}' '${params.taxdump_dir}'
    """
}

// rRNA anchors (pybarrnap or barrnap) + ITSx (spacers) on ONE assembly unit's
// scaffold sequence (pulled from the assemblies table), with tool/kingdom/taxon
// per-unit from annotate_opts. Annotation now runs per (ID, path, scaffold).
process annotate_tools {
    cache false
    tag "${id}.${upath}.${uscaf}"
    // One scaffold's barrnap/ITSx failure must not sink the whole sample: WF1
    // now annotates every unit, including short fragments and contaminants.
    errorStrategy 'ignore'
    cpus { nthreads as int }
    // scales on retry; see the OOM retry policy in nextflow.config
    memory { "${(mem as int) * task.attempt} GB" }
    input:
    tuple val(id), val(upath), val(uscaf), val(nthreads), val(mem)
    output:
    tuple val(id), val(upath), val(uscaf), path("unit.barrnap.gff"), path("unit.itsx.positions.txt")
    script:
    """
    # write this unit's scaffold sequence to a fasta straight from the DB
    sqlite3 -cmd '.timeout 10000' '${params.db}' \
      "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM assemblies
       WHERE ID='${id}' AND path=${upath} AND scaffold=${uscaf};" > unit.fasta

    read TOOL KINGDOM ITSXT ITSXO <<< "\$(sqlite3 -cmd '.timeout 10000' -separator ' ' '${params.db}' "
      SELECT COALESCE(ao.rrna_tool,'pybarrnap'), COALESCE(ao.rrna_kingdom,'euk'),
             COALESCE(ao.itsx_taxon,'all'), COALESCE(ao.itsx_opts,'')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"
    TOOL="\${TOOL:-pybarrnap}"; KINGDOM="\${KINGDOM:-euk}"; ITSXT="\${ITSXT:-all}"

    # Same kingdom translation as R's rrna_kingdom_for(): the Perl barrnap fork
    # calls its eukaryote Rfam models "fun" and has no euk/all; pybarrnap has no
    # "fun". Passing the stored kingdom through untranslated makes the tool exit
    # non-zero and the unit gets no anchors at all.
    if [ "\$TOOL" = "barrnap" ]; then
        case "\$KINGDOM" in euk|all) KINGDOM=fun ;; esac
        barrnap --kingdom "\$KINGDOM" --quiet unit.fasta > unit.barrnap.gff
    else
        case "\$KINGDOM" in fun) KINGDOM=euk ;; esac
        pybarrnap -k "\$KINGDOM" -t ${task.cpus} -q unit.fasta > unit.barrnap.gff
    fi

    ITSx -i unit.fasta -o unit.itsx --save_regions none --graphical F \
        --cpu ${task.cpus} -t "\$ITSXT" --preserve T \$ITSXO || true
    [ -f unit.itsx.positions.txt ] || : > unit.itsx.positions.txt
    """
}

// Reconcile barrnap + ITSx into the feature model -> annotations table (per unit).
process reconcile {
    cache false
    tag "${id}.${upath}.${uscaf}"
    // One unit's malformed GFF must not terminate the run, same as annotate_tools.
    errorStrategy 'ignore'
    input:
    tuple val(id), val(upath), val(uscaf), path(gff), path(pos)
    output:
    tuple val(id), val(upath), val(uscaf)
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_reconcile.R '${params.db}' '${id}' ${upath} ${uscaf} ${gff} ${pos}
    """
}

// Curate the reconciled feature model per unit: structural + reference flags ->
// per-feature status/flags + per-unit n_flags/notes.
// Propose keep/drop/undecided per unit once a sample's evidence is complete
// (features from reconcile, screen from blast_genbank, depth from coverage).
// Advisory only: a drop needs positive evidence, and a user decision is sticky.
process triage {
    cache false
    tag "${id}"
    errorStrategy 'ignore'
    input:
    val id
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_triage.R '${params.db}' '${id}'
    """
}

// Machine contamination verdict per unit, from the retained BLAST hit
// distribution. Best-effort: a taxonomy failure must not sink a sample.
process screen {
    cache false
    tag "${id}"
    errorStrategy 'ignore'
    input:
    val id
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_screen.R '${params.db}' '${id}'
    """
}

process curate {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    input:
    tuple val(id), val(upath), val(uscaf)
    output:
    tuple val(id), val(upath), val(uscaf)
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_curate.R '${params.db}' '${id}' ${upath} ${uscaf}
    """
}
