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
    sqlite3 -header -separator '\t' '${params.db}' "
      SELECT p.ID, p.R1, p.R2, ao.cpus, po.fastp AS fastp, ao.getOrganelle AS getOrganelle,
             ao.seeds_db, ao.labels_db
      FROM assemble a
      JOIN preprocess p     ON a.ID = p.ID
      JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      JOIN pre_opts po      ON p.pre_opts = po.pre_opts
      WHERE a.assemble_switch IN (1,4) AND a.assemble_lock = 0;" > samplesheet.tsv
    """
}

// Fallback: samples queued for an assembly fallback (assemble_switch = 5).
process query_fallback_samples {
    cache false
    executor 'local'
    output:
    path 'samplesheet.tsv'
    script:
    """
    sqlite3 -header -separator '\t' '${params.db}' "
      SELECT p.ID, p.R1, p.R2, ao.cpus, po.fastp AS fastp
      FROM assemble a
      JOIN preprocess p     ON a.ID = p.ID
      JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      JOIN pre_opts po      ON p.pre_opts = po.pre_opts
      WHERE a.assemble_switch = 5 AND a.assemble_lock = 0;" > samplesheet.tsv
    """
}

// WF2: IDs of assembled samples queued for (re-)annotation.
process query_annotate_samples {
    cache false
    executor 'local'
    output:
    path 'ids.txt'
    script:
    """
    sqlite3 '${params.db}' "
      SELECT an.ID FROM annotate an
      JOIN assemble a ON a.ID = an.ID
      WHERE an.annotate_switch IN (1,4) AND an.annotate_lock = 0 AND a.assemble_lock = 1
        AND EXISTS (SELECT 1 FROM assemblies x WHERE x.ID = an.ID);" > ids.txt
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
process stitch {
    cache false
    tag "${id}"
    input:
    val id
    output:
    tuple val(id), path("${id}.stitched.fasta")
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_stitch.R '${params.db}' '${id}' ${id}.stitched.fasta
    """
}

// barrnap (Rfam anchors) + ITSx (spacers) on the stitched sequence.
process annotate_tools {
    cache false
    tag "${id}"
    cpus 4
    input:
    tuple val(id), path(stitched)
    output:
    tuple val(id), path("${id}.barrnap.gff"), path("${id}.itsx.positions.txt")
    script:
    """
    barrnap --kingdom fun --quiet ${stitched} > ${id}.barrnap.gff
    ITSx -i ${stitched} -o ${id}.itsx --save_regions none --graphical F \
        --cpu ${task.cpus} -t all --preserve T || true
    [ -f ${id}.itsx.positions.txt ] || : > ${id}.itsx.positions.txt
    """
}

// Reconcile barrnap + ITSx into the feature model -> annotations table.
process reconcile {
    cache false
    tag "${id}"
    input:
    tuple val(id), path(gff), path(pos)
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_reconcile.R '${params.db}' '${id}' ${gff} ${pos}
    """
}

// Curate the reconciled feature model: structural + reference flags -> per-
// feature status/flags + per-sample n_flags/notes. Reads any reference identity
// already in the DB (blast_ref_alignment); never blocks on remote BLAST.
process curate {
    cache false
    tag "${id}"
    input:
    val id
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_curate.R '${params.db}' '${id}'
    """
}
