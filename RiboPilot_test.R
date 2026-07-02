#library(RiboPilot)   # or
pkgload::load_all("~/Documents/GitHub/RiboPilot")

Sys.setenv(PATH = paste("/home/dmacguig/miniconda3/envs/ribopilot/bin",
                        Sys.getenv("PATH"), sep = ":"))
# verify:
Sys.which(c("nextflow", "fastp", "get_organelle_from_reads.py"))

# 1-sample project pointing at the reads you already have
new_project("/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/run_02",
            mapping_fn  = "/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/mapping.csv",
            data_path   = "~/Documents/GitHub/RiboPilot/testdata_full",
            custom_seeds_db  = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
            custom_labels_db = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
            executor = "conda", min_depth = 500,
            Rproj = FALSE)

run_assembly("~/Documents/DJMCompBio/RiboPilot_test/run_02", conda_env = "ribopilot")
# review in RiboPilot::run_app(...), adjust selected paths
run_annotate("~/Documents/DJMCompBio/RiboPilot_test", conda_env = "ribopilot")
# then Export mode / export_bundle()

RiboPilot::run_app("~/Documents/DJMCompBio/RiboPilot_test/run_02")
