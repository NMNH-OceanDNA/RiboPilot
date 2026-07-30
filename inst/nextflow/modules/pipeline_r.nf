// WF1 downstream processes: read/write the project SQLite DB directly (sqlite3
// CLI for the sample query; R-script processes for scoring / stitching /
// annotation reconciliation). No R around Nextflow, no nf-sqldb plugin.

// Query the samples ready to assemble straight from the DB (replaces the R
// samplesheet prep). Header + tab separated for splitCsv.
process query_samples {
    cache false
    executor 'local'
    output:
    path 'samplesheet.tsv'
    script:
    """
    # Header is written explicitly (sqlite3 -header emits nothing for 0 rows, and
    # splitCsv(header:true) then fails); an empty result becomes a header-only
    # file so the workflow simply does nothing instead of crashing.
    printf 'ID\tR1\tR2\tpre_cpus\tpre_mem\tasm_cpus\tasm_mem\tfastp\tgetOrganelle\tseeds_db\tlabels_db\tassembler\n' > samplesheet.tsv
    sqlite3 -separator '\t' '${params.db}' "
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
process query_annotate_samples {
    cache false
    executor 'local'
    output:
    path 'samplesheet.tsv'
    script:
    """
    printf 'ID\tpath\tscaffold\tcpus\tmem\n' > samplesheet.tsv
    sqlite3 -separator '\t' '${params.db}' "
      SELECT an.ID, an.path, an.scaffold, COALESCE(ao.cpus,4) AS cpus, COALESCE(ao.memory,8) AS mem
      FROM annotate an
      JOIN assemble a ON a.ID = an.ID
      JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
      WHERE an.annotate_switch IN (1,4) AND an.annotate_lock = 0 AND a.assemble_lock = 1
        AND COALESCE(asm.ignore,0) = 0;" >> samplesheet.tsv
    """
}

// WF1: the assembly units of ONE sample that has just been assembled, so
// annotation can run over every scaffold as it is produced rather than waiting
// for the user to review and lock. Unlike query_annotate_samples this does not
// gate on assemble_lock or annotate_switch: nothing has been reviewed yet.
process query_sample_units {
    cache false
    executor 'local'
    tag "${id}"
    input:
    val id
    output:
    tuple val(id), path('units.tsv')
    script:
    """
    printf 'ID\tpath\tscaffold\tcpus\tmem\n' > units.tsv
    sqlite3 -separator '\t' '${params.db}' "
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

// Stitch the selected path into one sequence -> stitched table + fasta.
// A sample with no usable assembly has nothing to stitch; that is a recorded
// outcome (load_assemblies marks it failed), not an error, so the fasta is
// optional and the sample simply drops out of the downstream channel instead
// of aborting every other sample in the run.
process stitch {
    cache false
    tag "${id}"
    input:
    val id
    output:
    tuple val(id), path("${id}.stitched.fasta"), optional: true
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_stitch.R '${params.db}' '${id}' ${id}.stitched.fasta
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
    memory { "${mem} GB" }
    input:
    tuple val(id), val(upath), val(uscaf), val(nthreads), val(mem)
    output:
    tuple val(id), val(upath), val(uscaf), path("unit.barrnap.gff"), path("unit.itsx.positions.txt")
    script:
    """
    # write this unit's scaffold sequence to a fasta straight from the DB
    sqlite3 '${params.db}' \
      "SELECT '>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM assemblies
       WHERE ID='${id}' AND path=${upath} AND scaffold=${uscaf};" > unit.fasta

    read TOOL KINGDOM ITSXT ITSXO <<< "\$(sqlite3 -separator ' ' '${params.db}' "
      SELECT COALESCE(ao.rrna_tool,'pybarrnap'), COALESCE(ao.rrna_kingdom,'euk'),
             COALESCE(ao.itsx_taxon,'all'), COALESCE(ao.itsx_opts,'')
      FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts
      WHERE an.ID = '${id}' AND an.path = ${upath} AND an.scaffold = ${uscaf};")"
    TOOL="\${TOOL:-pybarrnap}"; KINGDOM="\${KINGDOM:-euk}"; ITSXT="\${ITSXT:-all}"

    if [ "\$TOOL" = "barrnap" ]; then
        barrnap --kingdom "\$KINGDOM" --quiet unit.fasta > unit.barrnap.gff
    else
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
    input:
    tuple val(id), val(upath), val(uscaf)
    output:
    tuple val(id), val(upath), val(uscaf)
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_curate.R '${params.db}' '${id}' ${upath} ${uscaf}
    """
}
