# Try out user-supplied assemblies end to end.
#
# Builds a small synthetic assembly dataset, sets up a project from it, and
# opens the app. Nothing is downloaded: the contigs are built from the real fish
# nuclear rDNA repeat units bundled in ref_dbs/, padded with random sequence to
# look like genome scaffolds. See dev/design/20_USER_ASSEMBLIES.md.
#
# Run the pipeline from the app (Assemble tab -> Update), not from here.

pkgload::load_all("~/Documents/GitHub/RiboPilot")

Sys.setenv(PATH = paste("<PATH_TO_CONDA>/envs/ribopilot/bin",
                        Sys.getenv("PATH"), sep = ":"))
# verify: pybarrnap is what searches the supplied assemblies
Sys.which(c("nextflow", "pybarrnap", "ITSx"))

REPO     <- "~/Documents/GitHub/RiboPilot"
PROJECT  <- "~/RiboPilot_test/userAsmb_01"   # project folder to create
DATA     <- "~/RiboPilot_test/userAsmb_data" # where the test assemblies are written

SEEDS  <- file.path(REPO, "ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta")
LABELS <- file.path(REPO, "ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta")

# ---------------------------------------------------------------- test data --

# Records taken from the bundled seed database. Each is one complete nuclear
# rDNA repeat unit. The 5.8S window is where pybarrnap calls 5.8S in that
# record; it is used to build a contig that carries 5.8S and nothing else, which
# is what a contig looks like when the search decides not to carry it forward.
UNITS <- list(
  furzeri  = list(acc = "EU780557.1_Nothobranchius_furzeri",
                  taxon = "Nothobranchius furzeri", s58 = c(5419, 5571)),
  scitulus = list(acc = "OP151354.1_Prionotus_scitulus",
                  taxon = "Prionotus scitulus", s58 = c(1, 152)),
  rubio    = list(acc = "OP151435.1_Prionotus_rubio",
                  taxon = "Prionotus rubio", s58 = c(2354, 2506))
)

seeds <- Biostrings::readDNAStringSet(SEEDS)
names(seeds) <- sub("\\s.*$", "", names(seeds))
missing <- setdiff(vapply(UNITS, `[[`, "", "acc"), names(seeds))
if (length(missing)) {
  stop("the bundled seed database no longer has: ", paste(missing, collapse = ", "))
}

set.seed(20260901)
rand_contig <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
unit_seq <- function(key) as.character(seeds[[UNITS[[key]]$acc]])

# one rDNA unit dropped into the middle of a longer scaffold, the way a repeat
# sits inside a real genome contig
rdna_scaffold <- function(key, flank = 4000) {
  paste0(rand_contig(flank), unit_seq(key), rand_contig(flank))
}

# 5.8S and its immediate surroundings only: no 18S or 28S anchor to find
s58_only <- function(key, pad = 300) {
  w <- UNITS[[key]]$s58
  u <- unit_seq(key)
  substr(u, max(1, w[1] - pad), min(nchar(u), w[2] + pad))
}

write_asmb <- function(dir, id, seqs) {
  s <- Biostrings::DNAStringSet(unlist(seqs, use.names = FALSE))
  names(s) <- names(seqs)
  Biostrings::writeXStringSet(s, file.path(dir, paste0(id, ".fasta")))
}

dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
asmb_dir <- file.path(DATA, "assemblies")
dir.create(asmb_dir, showWarnings = FALSE)

# 1. a whole genome assembly: decoy scaffolds, one carrying the rDNA repeat, and
#    one carrying 5.8S on its own
write_asmb(asmb_dir, "ASM_genome", c(
  stats::setNames(
    lapply(c(31000, 18500, 12000, 9000, 7500, 5200, 4100), rand_contig),
    sprintf("scaffold_%02d", 1:7)),
  list(scaffold_08_rdna = rdna_scaffold("furzeri"),
       scaffold_09_partial = s58_only("furzeri"),
       scaffold_10 = rand_contig(3000))))

# 2. a single rDNA contig, the way it comes out of a targeted assembly
write_asmb(asmb_dir, "ASM_contig",
  list(Pscitulus_rDNA_unit = unit_seq("scitulus")))

# 3. a genome carrying two separate rDNA-bearing scaffolds
write_asmb(asmb_dir, "ASM_two_units", c(
  stats::setNames(lapply(c(22000, 14000, 6000), rand_contig),
    sprintf("ctg_%03d", 1:3)),
  list(ctg_004_rdna = rdna_scaffold("rubio", flank = 2500),
       ctg_005_rdna = rdna_scaffold("rubio", flank = 6000))))

# 4. a genome with no rDNA at all: the sample should fail the search and say so
write_asmb(asmb_dir, "ASM_no_rdna",
  stats::setNames(lapply(c(26000, 15000, 8000, 4000), rand_contig),
    sprintf("contig%d", 1:4)))

mapping <- data.frame(
  ID       = c("ASM_genome", "ASM_contig", "ASM_two_units", "ASM_no_rdna"),
  Taxon    = c(UNITS$furzeri$taxon, UNITS$scitulus$taxon,
               UNITS$rubio$taxon, "Danio rerio"),
  R1       = "",
  R2       = "",
  Assembly = c("ASM_genome.fasta", "ASM_contig.fasta",
               "ASM_two_units.fasta", "ASM_no_rdna.fasta"),
  Note     = c("genome, one rDNA scaffold + a 5.8S-only decoy",
               "single rDNA contig",
               "genome, two rDNA scaffolds",
               "genome with no rDNA (expected to fail the search)"),
  stringsAsFactors = FALSE
)
mapping_fn <- file.path(DATA, "mapping.csv")
utils::write.csv(mapping, mapping_fn, row.names = FALSE)

message("test assemblies written to ", asmb_dir)
print(mapping[, c("ID", "Taxon", "Assembly", "Note")])

# ------------------------------------------------------------------ project --

# No data_path: none of these samples has reads. Coverage is the only thing they
# lose by that; everything else runs as usual.
new_project(
  PROJECT,
  mapping_fn       = mapping_fn,
  assembly_path    = asmb_dir,
  custom_seeds_db  = SEEDS,
  custom_labels_db = LABELS,
  executor         = "conda",
  Rproj            = FALSE
)

# Now run it from the app: Assemble tab, tick the samples, press Update.
#
# What to look for:
#  - ASM_genome    one scaffold carried forward, and a "Contigs not carried
#                  forward" section in its Details modal listing the 5.8S-only
#                  scaffold with the reason
#  - ASM_contig    one scaffold, the whole supplied contig
#  - ASM_two_units two scaffolds, one per rDNA-bearing contig
#  - ASM_no_rdna   fails with "no rRNA genes found in the supplied assembly"
#  - every sample  a "Your contig" column naming the contig it came from, and
#                  empty coverage, since there are no reads
#
# To try the opt-out: Assemble tab -> assembly options -> untick "Search the
# supplied assembly for nuclear rDNA contigs". The file is then used exactly as
# supplied, which is what a short ITS-only fragment needs. Sending a whole
# genome through that way trips the max scaffolds limit, which is the right
# answer.
RiboPilot::RiboPilot(PROJECT)
