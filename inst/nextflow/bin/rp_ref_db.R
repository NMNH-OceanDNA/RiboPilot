#!/usr/bin/env Rscript
# Run the bundled GenBank rRNA script for one clade.
# Args: count <taxid> <query> <outfile> | build <prefix> <query> <outdir>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))

script <- app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
if (!nzchar(script) || !file.exists(script)) stop("download_rrna.sh not found")
api_key <- Sys.getenv("NCBI_API_KEY", "")
refseq <- identical(Sys.getenv("RP_REFSEQ", ""), "true")

args <- c(shQuote(script))
if (refseq) args <- c(args, "--refseq")
if (nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))

# stdin must be detached: Entrez Direct reads stdin when it is open and blocks
if (identical(a[1], "count")) {
  args <- c(args, "--count-only", shQuote(a[3]), "count_probe")
  out <- suppressWarnings(system2("bash", args, stdout = TRUE, stderr = TRUE,
    stdin = "/dev/null"))
  n <- parse_record_count(out)
  if (is.na(n)) { cat(paste(out, collapse = "\n"), "\n"); stop("count failed") }
  writeLines(paste(a[2], n, sep = "\t"), a[4])
} else {
  args <- c(args, "--outdir", shQuote(a[4]), shQuote(a[3]), shQuote(a[2]))
  status <- system2("bash", args, stdin = "/dev/null")
  if (status != 0) stop("build failed for ", a[2])
}
