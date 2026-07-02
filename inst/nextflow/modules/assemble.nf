process assemble {

    executor params.assemble.executor
    container params.assemble.container
    maxRetries params.assemble.maxRetries

    publishDir "${launchDir}/${params.publishDir}", overwrite: true, pattern: "${id}/assemble/**", mode: 'copy'

    errorStrategy 'ignore'

    tag "${id}"

    input:
    tuple val(id), path(reads), val(opts), path(dbs)

    output:
    // formatted per-path assembly fastas for DB loading
    tuple val("${id}"), path("${id}/assemble/${opts_id}/${id}_assembly_*.fasta"), val("${opts_id}")
    // assembler summary + logs published
    path "${id}/assemble/${opts_id}/*"

    shell:
    opts_id = "default"
    workingDir = "${id}/assemble"
    outDir = "${workingDir}/${opts_id}"
    '''
    mkdir -p !{workingDir}
    get_organelle_from_reads.py \
        -1 !{reads[0]} \
        -2 !{reads[1]} \
        -o !{workingDir}/ --overwrite \
        -s !{dbs[0]} \
        --genes !{dbs[1]} \
        -t !{opts.cpus} \
        !{opts.getOrganelle}

    mkdir -p !{outDir}
    ### LOGS ###
    cp !{workingDir}/get_org.log.txt !{outDir}/assembler.log.txt 2>/dev/null || true
    echo "!{opts.getOrganelle}" > !{outDir}/opts.txt

    ### FORMAT ASSEMBLIES ###
    # rename each GetOrganelle path fasta to <id>_assembly_<n>.fasta and tag
    # each header ">id.<path>.<scaffold> <topology>"
    shopt -s nullglob
    files=(!{workingDir}/*.fasta)
    if [ ${#files[@]} -eq 0 ]; then
        echo ">No assembly found" > !{outDir}/!{id}_assembly_0.fasta
    else
        n=0
        for f in "${files[@]}"; do
            n=$((n+1))
            # infer topology from filename (GetOrganelle marks circular/complete)
            topo="linear"
            if echo "$f" | grep -qi "circular\|complete"; then topo="circular"; fi
            awk -v topo="$topo" -v id="!{id}" -v p="$n" \
                '/^>/ {print ">" id "." p "." ++c " " topo} !/^>/ {print}' \
                "$f" > !{outDir}/!{id}_assembly_${n}.fasta
        done
    fi
    '''
}
