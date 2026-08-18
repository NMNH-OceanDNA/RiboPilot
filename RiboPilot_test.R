#library(RiboPilot)   # or
pkgload::load_all("~/Documents/GitHub/RiboPilot")

Sys.setenv(PATH = paste("<PATH_TO_CONDA>/envs/ribopilot/bin",
                        Sys.getenv("PATH"), sep = ":"))
# verify:
Sys.which(c("nextflow", "fastp", "get_organelle_from_reads.py"))

# 1-sample project pointing at the reads you already have
new_project("<YOUR_DIR>/RiboPilot_test/run_03",
            mapping_fn  = "<YOUR_DIR>/RiboPilot_test/mapping.csv",
            data_path   = "~/Documents/GitHub/RiboPilot/testdata_full",
            custom_seeds_db  = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
            custom_labels_db = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
            executor = "conda", min_depth = 500,
            Rproj = FALSE)

run_assembly("<YOUR_DIR>/RiboPilot_test/run_03", conda_env = "ribopilot")
# review in RiboPilot::run_app(...), adjust selected paths, then re-run
# run_assembly() (or the Assemble tab's Update button) to re-annotate what changed
# then Export mode / export_bundle()

RiboPilot::RiboPilot("<YOUR_DIR>/RiboPilot_test/run_03")
#
