#library(RiboPilot)   # or
pkgload::load_all("~/Documents/GitHub/RiboPilot")

proj <- "/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/2026_09_02"

# 1-sample project pointing at the reads you already have
new_test_project(
  path      = proj,
  #n         = 0,
  executor  = "local",
  container = "ribopilot:userAsmb",
  min_depth = 500,
  Rproj     = FALSE,
  force     = TRUE,
  sra_bin = "~/miniconda3/envs/SRA_toolkit/bin"
)

RiboPilot::RiboPilot(proj)

