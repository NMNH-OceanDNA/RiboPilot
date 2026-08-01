#!/usr/bin/env Rscript
# Run the bundled GenBank rRNA script for one clade.
# Args: count <taxid> <query> <outfile> | build <ref_id> <query> <outdir>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))

api_key <- Sys.getenv("NCBI_API_KEY", "")
refseq <- identical(Sys.getenv("RP_REFSEQ", ""), "true")

if (identical(a[1], "count")) {
  script <- RiboPilot:::app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
  if (!nzchar(script) || !file.exists(script)) stop("download_rrna.sh not found")
  args <- c(shQuote(script))
  if (refseq) args <- c(args, "--refseq")
  if (nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))
  args <- c(args, "--count-only", shQuote(a[3]), "count_probe")
  # stdin must be detached: Entrez Direct reads stdin when it is open and blocks
  out <- suppressWarnings(system2("bash", args, stdout = TRUE, stderr = TRUE,
    stdin = "/dev/null"))
  n <- RiboPilot:::parse_record_count(out)
  if (is.na(n)) { cat(paste(out, collapse = "\n"), "\n"); stop("count failed") }
  writeLines(paste(a[2], n, sep = "\t"), a[4])
} else {
  # naming scheme and the seeds/labels layout stay defined in one place
  RiboPilot::build_reference_db(query = a[3], clade = a[2], out_dir = a[4],
    api_key = if (nzchar(api_key)) api_key else NULL, refseq = refseq,
    overwrite = TRUE)
}
