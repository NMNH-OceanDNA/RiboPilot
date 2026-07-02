# Shared bootstrap sourced by the WF1/WF2 R-script processes. Loads RiboPilot
# (installed package if present, else the source dir given by RIBOPILOT_SRC).
suppressMessages({
  if (requireNamespace("RiboPilot", quietly = TRUE)) {
    library(RiboPilot)
  } else {
    src <- Sys.getenv("RIBOPILOT_SRC")
    if (!nzchar(src)) stop("RiboPilot not installed and RIBOPILOT_SRC not set")
    pkgload::load_all(src, quiet = TRUE, attach_testthat = FALSE)
  }
})
