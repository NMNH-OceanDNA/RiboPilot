# Builds the synthetic chimera-detector fixtures. Run once; outputs are
# committed, this script is not re-run by tests.
#
# Raw inputs (not committed, regenerate before re-running this script):
#   sqlite3 /home/dmacguig/Documents/DJMCompBio/RiboPilot_test/run_02/.sqlite \
#     "SELECT sequence FROM assemblies WHERE ID='SRR22396758' AND path=1 AND scaffold=3;" \
#     > _host_raw.txt
#   blastdbcmd -db /opt/ribopilot/blastdb/2026.08.18/LSU_eukaryote_rRNA \
#     -entry AF003355.1 -outfmt "%f" > fungal_28S_raw.fa
#
# Host scaffold: SRR22396758 path 1 scaffold 3 (Upeneus parvus, 6310 bp) from
# the local test project run_02 DB. That project predates per-scaffold
# annotation, so no stored 18S/ITS1/5.8S/ITS2/28S coordinates exist; boundaries
# below were estimated with a local blastn against the shipped rDNA panel and
# the NCBI SSU/LSU decoy sets (best hit was a "SSU partial; ITS1, 5.8S, ITS2
# complete; LSU partial" GenBank record aligning host 2151-3994 at ~99.5%
# identity). Boundaries are approximate, not annotated ground truth; they only
# need to be biologically plausible splice points for a synthetic test.
#
# Fungal 28S: AF003355.1, Penicillium aurantiogriseum, fetched from the NCBI
# LSU_eukaryote_rRNA decoy db via blastdbcmd. No fungal ITS record was at
# hand, so host_fungal.fa splices host 18S+ITS1 directly onto the fungal 28S
# (skipping 5.8S/ITS2), noted in expected.tsv.

set.seed(42)

host <- paste(readLines("_host_raw.txt"), collapse = "")
fungal <- paste(readLines("fungal_28S_raw.fa")[-1], collapse = "")

stopifnot(nchar(host) == 6310)

HOST_18S_END <- 2150L
HOST_ITS1_END <- 2650L
HOST_28S_START <- 3991L

write_fasta <- function(id, seq, file) {
  wrapped <- paste(strwrap(seq, width = 70), collapse = "\n")
  writeLines(c(paste0(">", id), wrapped), file)
}

# clean.fa: unmodified host scaffold
write_fasta("clean", host, "clean.fa")

# host_fungal.fa: host 18S+ITS1 (1..2650) joined to fungal 28S
host_fungal <- paste0(substr(host, 1, HOST_ITS1_END), fungal)
write_fasta("host_fungal", host_fungal, "host_fungal.fa")

# host_variant.fa: 28S positions 500..1500 (relative to 28S start) mutated at
# 6% of sites
bases <- c("A", "C", "G", "T")
variant <- host
lo <- HOST_28S_START + 500L - 1L
hi <- HOST_28S_START + 1500L - 1L
window <- lo:hi
n_mut <- round(length(window) * 0.06)
mut_pos <- sort(sample(window, n_mut))
chars <- strsplit(variant, "")[[1]]
for (p in mut_pos) {
  orig <- chars[p]
  alt <- setdiff(bases, toupper(orig))
  chars[p] <- sample(alt, 1)
}
variant <- paste(chars, collapse = "")
write_fasta("host_variant", variant, "host_variant.fa")
cat("host_variant mutated positions (host coords):", mut_pos[1], "..", mut_pos[length(mut_pos)],
  "(", n_mut, "sites )\n")

# tandem_18S.fa: host 18S duplicated in tandem before ITS1
tandem <- paste0(substr(host, 1, HOST_18S_END), substr(host, 1, HOST_18S_END), substr(host, HOST_18S_END + 1, nchar(host)))
write_fasta("tandem_18S", tandem, "tandem_18S.fa")

cat("host_fungal.fa junction (host/fungal):", HOST_ITS1_END, "\n")
cat("tandem_18S.fa junctions (dup1/dup2, dup2/ITS1):", HOST_18S_END, HOST_18S_END * 2, "\n")
