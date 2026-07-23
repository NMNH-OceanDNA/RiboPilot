// RiboPilot WF1 (self-contained): read samples from the project DB via sqlite3,
// preprocess (fastp) + assemble (GetOrganelle), then score the candidate paths.
// score() also loads the per-scaffold `assemblies` rows and seeds one annotate
// unit per non-ignored scaffold. Assemblies are NOT stitched; each scaffold is
// annotated independently in WF2 after the user reviews and locks the assembly.
// (Per-scaffold GenBank BLAST + coverage are added in a later step.)

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { query_samples; score } from './modules/pipeline_r.nf'
include { blast_genbank } from './modules/blast.nf'
include { blast_ref_fetch } from './modules/blast_ref.nf'
include { coverage } from './modules/coverage.nf'

params.publishDir = 'out'
params.minDepth = 0
params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.ncbi_api_key = ''

workflow {
    // samples straight from the DB (tab-separated, with header)
    rows = query_samples().splitCsv(header: true, sep: '\t')

    pre_in = rows.map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.pre_cpus, r.pre_mem, r.fastp) }
    preprocess(pre_in)

    asm_opts = rows.map { r -> tuple(r.ID, r.asm_cpus, r.asm_mem, r.getOrganelle, file(r.seeds_db), file(r.labels_db)) }
    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(asm_opts)
        .map { id, reads, cpus, mem, go, seeds, labels -> tuple(id, reads, cpus, mem, go, seeds, labels) }
    assemble(asm_in)

    // Load assemblies + score paths (also seeds annotate units + per-scaffold
    // BLAST). Then map reads to the scaffolds for per-unit coverage. Stops at the
    // assembly-review stage; annotation/curation happen per unit in WF2.
    score(assemble.out.paths)
    // per-scaffold remote GenBank BLAST -> ranked reference candidates + screen,
    // then fetch each candidate reference (sequence + features) from NCBI
    blast_genbank(score.out)
    blast_ref_fetch(blast_genbank.out)
    // per-unit coverage (independent branch off score)
    cov_in = score.out.join(preprocess.out.reads.map { id, reads, after -> tuple(id, reads) })
    coverage(cov_in)
}
