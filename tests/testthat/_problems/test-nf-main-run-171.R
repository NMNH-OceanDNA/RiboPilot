# Extracted from test-nf-main-run.R:171

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RiboPilot", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_on_cran()
skip_on_os("windows")
skip_if(is.na(RiboPilot:::nf_installed_version()), "Nextflow not installed or not runnable")
if (RiboPilot:::nf_version_status() != "supported") testthat::skip_if_offline()
for (b in c("sqlite3", "jq")) skip_if(!nzchar(Sys.which(b)), paste(b, "not installed"))
tools <- local_fake_pipeline_tools()
proj_dir <- withr::local_tempdir()
p <- file.path(proj_dir, "proj")
asmb_dir <- file.path(proj_dir, "asmb")
dir.create(asmb_dir)
file.copy(tools$scaffold, file.path(asmb_dir, "SUPPLIED_assembly.fa"))
map_fn <- file.path(proj_dir, "mapping.csv")
writeLines(c("ID,Taxon,Assembly", "SUPPLIED,Upeneus parvus,SUPPLIED_assembly.fa"), map_fn)
suppressMessages(suppressWarnings(RiboPilot::new_project(path = p, mapping_fn = map_fn,
    assembly_path = asmb_dir, executor = "conda", custom_seeds_db = tools$scaffold,
    custom_labels_db = tools$scaffold, Rproj = FALSE, fetch_metadata = FALSE)))
