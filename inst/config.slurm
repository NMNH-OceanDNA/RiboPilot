// RiboPilot SLURM config template. Turn it into a reusable profile with
//   generate_config(name = "my_cluster", scheduler = "slurm", queue = ..., ...)
// which fills the container-engine, queue, and cluster-options settings and
// saves it; then run new_project(executor = "my_cluster").
//
// RiboPilot does not use nf-sqldb: each process opens the project SQLite DB
// directly by path and pulls its own CPU count from the DB. The container must
// therefore see the project files, so the generated container-engine block
// bind-mounts the project directory, and new_project() adds the read,
// assembly, and reference folders. Requires a shared filesystem visible to
// all compute nodes.

// Throttle job submission to be a good cluster citizen (adjust as needed)
executor {
    queueSize = 50
    pollInterval = '30 sec'
    submitRateLimit = '10/1min'
}

<<CONTAINER_ENGINE>>

process {
    executor = 'slurm'
    container = '<<CONTAINER_ID>>'
    queue = '<<QUEUE>>'                         // SLURM partition (maps to -p)
    time = '<<WALLTIME>>'                       // per-job walltime (maps to -t); many sites reject a job without one
    clusterOptions = '<<CLUSTER_OPTIONS>>'      // extra sbatch directives, e.g. --account=...
    memory = { 16.GB * task.attempt }           // CPUs come from the pipeline; tune memory as needed
}

// Pipeline parameters
params {
    publishDir = 'out'                          // Directory to publish pipeline outputs
    rawDir = '<<RAW_DIR>>'                       // Directory containing raw read data
    asmbDir = '<<ASMB_DIR>>'                     // Directory containing user-supplied assemblies
    minDepth = <<MIN_DEPTH>>                     // Skip samples with fewer than this many reads
    ncbi_api_key = '<<NCBI_API_KEY>>'            // Optional NCBI API key for BLAST + GenBank fetch
}
