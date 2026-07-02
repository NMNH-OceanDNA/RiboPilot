#library(RiboPilot)   # or
pkgload::load_all("~/Documents/GitHub/RiboPilot")

Sys.setenv(PATH = paste("/home/dmacguig/miniconda3/envs/ribopilot/bin",
                        Sys.getenv("PATH"), sep = ":"))
# verify:
Sys.which(c("nextflow", "fastp", "get_organelle_from_reads.py"))

# 1-sample project pointing at the reads you already have
new_project("/home/dmacguig/Documents/DJMCompBio/RiboPilot_test",
            mapping_fn  = "/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/mapping.csv",
            data_path   = "~/Documents/GitHub/RiboPilot/testdata_full",
            custom_seeds_db  = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
            custom_labels_db = "~/Documents/GitHub/RiboPilot/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
            executor = "conda", min_depth = 500,
            Rproj = FALSE)

run_assembly("/home/dmacguig/Documents/DJMCompBio/RiboPilot_test")   # assemble -> load -> score

# inspect
con <- DBI::dbConnect(RSQLite::SQLite(), "/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/.sqlite")
DBI::dbGetQuery(con, "SELECT path,depth,anchors,blast_pident,score,rank,auto_select FROM path_scores")
