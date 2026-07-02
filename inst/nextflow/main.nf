// RiboPilot assembly pipeline (WF1): preprocess (fastp) + targeted nr assembly
// (GetOrganelle). Samplesheet-driven so it runs on any executor without the
// nf-sqldb plugin; the R driver builds the samplesheet from the project database
// and loads the published outputs back into it (assemblies + scoring).

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'

// params: samplesheet, rawDir, publishDir, minDepth
params.publishDir = 'out'
params.minDepth = 0

workflow {
    // samplesheet (TAB-separated so comma-containing options like GetOrganelle's
    // -k list need no quoting): ID,R1,R2,cpus,memory,fastp,getOrganelle,seeds_db,labels_db
    rows = channel.fromPath(params.samplesheet).splitCsv(header: true, sep: '\t')

    pre_in = rows.map { r ->
        tuple(r.ID,
              file(params.rawDir + '/' + r.R1),
              file(params.rawDir + '/' + r.R2),
              r.cpus, r.fastp)
    }
    preprocess(pre_in)

    // assembly options keyed by ID
    asm_opts = rows.map { r -> tuple(r.ID, r.cpus, r.getOrganelle, file(r.seeds_db), file(r.labels_db)) }

    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(asm_opts)
        .map { id, reads, cpus, go, seeds, labels -> tuple(id, reads, cpus, go, seeds, labels) }

    assemble(asm_in)
}
