#library(RiboPilot)   # or
pkgload::load_all("~/Documents/GitHub/RiboPilot")

proj <- "/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/2026_10_01"

Sys.setenv(GITHUB_PAT = "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")   # better: put in ~/.Renviron
Sys.setenv(NCBI_API_KEY = "a114385052cd7266dd59adff6814d68fa608")

# 1-sample project pointing at the reads you already have
RiboPilot::new_test_project(
  path      = proj,
  #n         = 0,
  full_size = FALSE,
  executor  = "local",
  container = "ribopilot:devel",
  min_depth = 500,
  Rproj     = FALSE,
  force     = TRUE,
  sra_bin = "~/miniconda3/envs/SRA_toolkit/bin"
)

RiboPilot::RiboPilot(proj)

