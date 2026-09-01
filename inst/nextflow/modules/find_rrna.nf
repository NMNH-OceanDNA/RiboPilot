// Find the nuclear rDNA contigs inside a user-supplied assembly, so a genome
// (or a single contig pulled out of one) enters WF1 where an assembly would.
// See dev/design/20_USER_ASSEMBLIES.md.
//
// The output is the same fasta contract the assembler emits
// (<id>_assembly_1.fasta, headers ">id.path.scaffold topology"), published to
// the same out/<id>/assemble folder, so score() and the whole tail of WF1 need
// no knowledge of where the sequences came from.

// Contigs per chunk of a user-supplied assembly sent to one rRNA scan task.
// Declared HERE, not in main.nf: includes are resolved before the including
// script's body runs, so a params default set in main.nf is still null when
// this module reads it. A config or --find_rrna_chunk_size value still wins.
params.find_rrna_chunk_size = 5000

// Scan one chunk of the assembly for rRNA genes. Chunked so a draft genome is
// searched in parallel and never read whole.
process find_rrna_scan {

    tag "${id}"
    cpus { nthreads as int }
    memory { "${(mem as int) * task.attempt} GB" }
    errorStrategy 'ignore'

    input:
    tuple val(id), path(chunk), val(nthreads), val(mem), val(tool), val(kingdom)

    output:
    tuple val(id), path('chunk.gff')

    script:
    """
    # Same kingdom translation as annotate_tools / R's rrna_kingdom_for(): the
    # Perl barrnap fork calls its eukaryote models "fun" and has no euk/all;
    # pybarrnap has no "fun".
    TOOL='${tool}'; KINGDOM='${kingdom}'
    TOOL="\${TOOL:-pybarrnap}"; KINGDOM="\${KINGDOM:-euk}"
    if [ "\$TOOL" = "barrnap" ]; then
        case "\$KINGDOM" in euk|all) KINGDOM=fun ;; esac
        barrnap --kingdom "\$KINGDOM" --threads ${task.cpus} --quiet ${chunk} > chunk.gff
    else
        case "\$KINGDOM" in fun) KINGDOM=euk ;; esac
        pybarrnap -k "\$KINGDOM" -t ${task.cpus} -q ${chunk} > chunk.gff
    fi
    """
}

// Pick the rDNA contigs out of the assembly using the merged scan, write them
// in the assembler's fasta format, and record the evidence.
process find_rrna_pick {

    tag "${id}"
    cpus 1
    memory { "${4 * task.attempt} GB" }
    publishDir { "${launchDir}/${params.publishDir}/${id}/assemble" }, overwrite: true, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(assembly), path(gffs)

    output:
    tuple val(id), path("${id}_assembly_*.fasta"), path('status.txt'), emit: picked
    path 'rrna_scan.gff'

    script:
    """
    cat ${gffs} > rrna_scan.gff
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_find_rrna.R \\
        '${params.db}' '${id}' ${assembly} rrna_scan.gff '${id}'
    # publishDir needs a file to match even when nothing was found; the sentinel
    # name is the one load_assemblies() already skips.
    [ -f ${id}_assembly_1.fasta ] || echo ">No rDNA contigs found" > ${id}_assembly_0.fasta
    """
}

// User-supplied assemblies in, confirmed rDNA contigs out. A sample where
// nothing is confirmed is failed by rp_find_rrna.R and emits nothing here, so
// it never reaches score().
workflow find_rrna {
    take:
        find_in     // tuple(id, assembly, cpus, mem, tool, kingdom)

    main:
        find_in
            .splitFasta(by: params.find_rrna_chunk_size, file: true, elem: 1)
            .set { chunks }

        find_rrna_scan(chunks)
            .groupTuple()
            .set { gffs_by_id }

        find_in
            .map { id, assembly, cpus, mem, tool, kingdom -> tuple(id, assembly) }
            .join(gffs_by_id)
            .set { pick_in }

        find_rrna_pick(pick_in)

        find_rrna_pick.out.picked
            .filter { id, fastas, status -> status.text.trim() == 'ok' }
            .map { id, fastas, status -> tuple(id, fastas) }
            .set { found }

    emit:
        paths = found
}

// A sample with the rDNA search switched off: take the supplied file as given.
// Only the header contract changes, so the rest of WF1 cannot tell the
// difference. Intended for a single supplied rDNA fragment; a whole genome sent
// through here trips the max_scaffolds gate, which is the correct answer.
process stage_assembly {

    tag "${id}"
    cpus 1
    memory 2.GB
    publishDir { "${launchDir}/${params.publishDir}/${id}/assemble" }, overwrite: true, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(assembly)

    output:
    tuple val(id), path("${id}_assembly_1.fasta"), emit: paths

    script:
    """
    awk -v id='${id}' \
        '/^>/ {orig=substr(\$0,2); print ">" id ".1." ++c " linear " orig; next} {print}' \
        ${assembly} > ${id}_assembly_1.fasta

    sqlite3 -cmd '.timeout 10000' '${params.db}' \\
      "UPDATE assemble SET find_rrna_notes =
         'using the supplied assembly as given (rDNA search off)' WHERE ID = '${id}';"
    """
}
