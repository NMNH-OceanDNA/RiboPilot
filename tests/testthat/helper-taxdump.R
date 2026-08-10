# A tiny real taxdump (dmp files) so the taxonomy tests exercise the actual
# index build and SQL, not a mock of it. Two genus *Morus* (mulberries and
# gannets) carry the cross-kingdom homonym case; Danio carries the ordinary one.

TEST_NODES <- list(
  # taxid, parent, rank
  c(1, 1, "no rank"),
  c(131567, 1, "no rank"),
  c(2759, 131567, "domain"),
  # plants
  c(33090, 2759, "kingdom"),
  c(35493, 33090, "phylum"),
  c(3398, 35493, "class"),
  c(3744, 3398, "order"),
  c(3487, 3744, "family"),
  c(3497, 3487, "genus"),
  c(3498, 3497, "species"),
  # animals
  c(33208, 2759, "kingdom"),
  c(7711, 33208, "phylum"),
  c(8782, 7711, "class"),
  c(30410, 8782, "order"),
  c(30446, 30410, "family"),
  c(37577, 30446, "genus"),
  c(186623, 7711, "class"),
  c(7952, 186623, "order"),
  c(2743709, 7952, "family"),
  c(7954, 2743709, "genus"),
  c(7955, 7954, "species"),
  c(1841699, 7954, "species")
)

TEST_NAMES <- list(
  # taxid, name, class
  c(1, "root", "scientific name"),
  c(131567, "cellular organisms", "scientific name"),
  c(2759, "Eukaryota", "scientific name"),
  c(33090, "Viridiplantae", "scientific name"),
  c(35493, "Streptophyta", "scientific name"),
  c(3398, "Magnoliopsida", "scientific name"),
  c(3744, "Rosales", "scientific name"),
  c(3487, "Moraceae", "scientific name"),
  c(3497, "Morus", "scientific name"),
  c(3498, "Morus alba", "scientific name"),
  c(33208, "Metazoa", "scientific name"),
  c(33208, "Animalia", "synonym"),
  c(7711, "Chordata", "scientific name"),
  c(8782, "Aves", "scientific name"),
  c(30410, "Suliformes", "scientific name"),
  c(30446, "Sulidae", "scientific name"),
  c(37577, "Morus", "scientific name"),
  c(186623, "Actinopteri", "scientific name"),
  c(7952, "Cypriniformes", "scientific name"),
  c(2743709, "Danionidae", "scientific name"),
  c(7954, "Danio", "scientific name"),
  c(7955, "Danio rerio", "scientific name"),
  c(7955, "zebrafish", "genbank common name"),
  c(7955, "Brachydanio rerio", "synonym"),
  c(7955, "Danio rerio Hamilton, 1822", "authority"),
  c(1841699, "Danio aesculapii", "scientific name")
)

# old taxid -> current taxid
TEST_MERGED <- list(c(29144, 7955))

#' Write the fixture dmp files into a fresh directory and return its path
local_taxdump <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  dmp <- function(rows) {
    paste0(vapply(rows, function(r) paste(r, collapse = "\t|\t"), character(1)), "\t|")
  }
  writeLines(dmp(TEST_NODES), file.path(dir, "nodes.dmp"))
  # names.dmp carries a "unique name" column between the name and its class
  writeLines(dmp(lapply(TEST_NAMES, function(r) c(r[1], r[2], "", r[3]))),
    file.path(dir, "names.dmp"))
  writeLines(dmp(TEST_MERGED), file.path(dir, "merged.dmp"))
  dir
}
