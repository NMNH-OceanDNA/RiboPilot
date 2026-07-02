// RiboPilot annotation workflow (WF2): run barrnap (18S/28S anchors) + ITSx
// (ITS1/5.8S/ITS2 boundaries) on each sample's stitched selected path. The R
// driver prepares the stitched sequences (WF1 output + selected path) and loads
// + reconciles the tool outputs into the feature model.

include { annotate } from './modules/annotate.nf'

params.publishDir = 'out'

workflow {
    // samplesheet (TAB): ID, stitched (fasta path), cpus
    rows = channel.fromPath(params.samplesheet).splitCsv(header: true, sep: '\t')
    ann_in = rows.map { r -> tuple(r.ID, file(r.stitched), r.cpus) }
    annotate(ann_in)
}
