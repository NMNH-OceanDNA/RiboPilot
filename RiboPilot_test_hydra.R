# RiboPilot test run on the NMNH Hydra cluster (SGE + Singularity)
# ---------------------------------------------------------------
# Run this from an INTERACTIVE job, not the login node:
#
#   qrsh -pe mthread 4 -l mres=16G,h_data=4G,h_vmem=4G
#   module load tools/R            # check `module avail` for exact names
#   module load bio/nextflow
#   module load bio/sratoolkit     # needed: ENA is blocked on SI network
#   module load tools/singularity  # or `tools/apptainer`
#   R
#
# One-time install of RiboPilot from the private repo (see notes at bottom):
#   Sys.setenv(GITHUB_PAT = "github_pat_xxx")   # better: put in ~/.Renviron
#   install.packages(c("remotes", "BiocManager"))
#   BiocManager::install(c("Biostrings", "DECIPHER", "pwalign"))
#   remotes::install_github("dmacguigan/RiboPilot")

library(RiboPilot)

# ---- paths -------------------------------------------------------------
# Use /scratch02 or /pool for project + work dirs. Home has a small quota and
# the read downloads alone are several GB.
proj <- "/scratch02/<YOUR_DIR>/ribopilot_test"

# Nextflow scratch. Keep off home.
Sys.setenv(NXF_HOME = "/scratch02/<YOUR_DIR>/.nextflow")
Sys.setenv(NXF_SINGULARITY_CACHEDIR = "/scratch02/<YOUR_DIR>/singularity_cache")
Sys.setenv(NXF_OPTS = "-Xms500M -Xmx4G")

# Optional but recommended: speeds up BLAST + GenBank lookups.
# Sys.setenv(NCBI_API_KEY = "your_ncbi_key")

# ---- sanity check ------------------------------------------------------
Sys.which(c("nextflow", "singularity", "fasterq-dump"))

# ---- create the test project ------------------------------------------
# Downloads reads, writes mapping.csv, initializes the DB and .config.
# n = 2 keeps the download small for a first smoke test; drop it for all 10.
new_test_project(
  path      = proj,
  n         = 2,
  executor  = "NMNH_Hydra",
  container = paste0("macguigand/ribopilot:", packageVersion("RiboPilot")),
  min_depth = 500,
  # ncbi_api_key = Sys.getenv("NCBI_API_KEY"),
  # fish rRNA seeds/labels; drop these to use the defaults
  custom_seeds_db  = system.file("ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
                                 package = "RiboPilot"),
  custom_labels_db = system.file("ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
                                 package = "RiboPilot"),
  Rproj = FALSE
)

# ---- REVIEW .config BEFORE RUNNING ------------------------------------
# file.edit(file.path(proj, ".config"))
#
# Two things to fix by hand in the singularity runOptions line:
#   1. "-B /scratch02/<<EXTRA_MOUNTS>>"  -> replace with your real mounts,
#      e.g. "-B /scratch02/<YOUR_DIR> -B /pool/public"
#   2. confirm the container tag exists on Docker Hub; if not, point it at a
#      real tag or a local .sif path.
# Also tune queueSize / memory in the process block if your queue is busy.

# ---- run the pipeline --------------------------------------------------
# Assembly + scoring + BLAST + annotation, submitted to SGE by Nextflow.
# Keep work dir on scratch.
run_assembly(proj, work = "/scratch02/<YOUR_DIR>/ribopilot_test_work")

# Re-run annotation only, after changing selected paths in the app.
# run_annotate(proj)

# ---- curate ------------------------------------------------------------
# Hydra has no browser. Either
#   (a) rsync the project dir back to your laptop and run the app locally, or
#   (b) SSH tunnel:  ssh -L 8080:localhost:8080 user@hydra
#       then in R on Hydra:
# RiboPilot(proj, options = list(port = 8080, host = "127.0.0.1", launch.browser = FALSE))

# ---- notes on the private repo ----------------------------------------
# Make a fine-grained token at github.com/settings/tokens with access limited
# to dmacguigan/RiboPilot and permission "Contents: read". Then on Hydra:
#
#   echo 'GITHUB_PAT=github_pat_xxx' >> ~/.Renviron
#   chmod 600 ~/.Renviron
#
# remotes::install_github() picks GITHUB_PAT up automatically. Login nodes have
# internet; compute nodes may not, so install and download reads from an
# interactive session that can reach the network.
