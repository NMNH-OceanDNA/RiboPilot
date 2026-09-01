# A project may supply its own assemblies instead of (or alongside) reads.
# See dev/design/20_USER_ASSEMBLIES.md.

mk_asmb_project <- function(mapping, assemblies = character(0), reads = character(0)) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  data_dir <- file.path(dir, "data")
  asmb_dir <- file.path(dir, "assemblies")
  dir.create(data_dir)
  dir.create(asmb_dir)
  file.create(file.path(data_dir, reads))
  file.create(file.path(asmb_dir, assemblies))
  map_fn <- file.path(dir, "mapping.csv")
  utils::write.csv(mapping, map_fn, row.names = FALSE)
  list(dir = dir, data_dir = data_dir, asmb_dir = asmb_dir, mapping_fn = map_fn,
    db = file.path(dir, ".sqlite"))
}

mk_asmb_db <- function(p, ...) {
  new_db(db_path = p$db, mapping_fn = p$mapping_fn, data_path = p$data_dir,
    assembly_path = p$asmb_dir, seeds_db = p$mapping_fn, labels_db = p$mapping_fn, ...)
}

test_that("a sample with an assembly and no reads is accepted", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "", R2 = "",
      Assembly = "S1.fa", stringsAsFactors = FALSE),
    assemblies = "S1.fa")
  expect_no_error(mk_asmb_db(p))

  con <- DBI::dbConnect(RSQLite::SQLite(), p$db)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT Assembly FROM samples")$Assembly, "S1.fa")
  expect_true("rrna_candidates" %in% DBI::dbListTables(con))
  expect_true("find_rrna_notes" %in% DBI::dbListFields(con, "assemble"))
  expect_equal(DBI::dbGetQuery(con, "SELECT find_rrna FROM assemble_opts")$find_rrna, 1L)
})

test_that("read and assembly samples mix in one project", {
  p <- mk_asmb_project(
    data.frame(ID = c("S1", "S2"), Taxon = "Danio rerio",
      R1 = c("S1_R1.fq.gz", ""), R2 = c("S1_R2.fq.gz", ""),
      Assembly = c("", "S2.fa"), stringsAsFactors = FALSE),
    assemblies = "S2.fa", reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  expect_no_error(mk_asmb_db(p))
})

test_that("a sample with neither reads nor an assembly errors", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "", R2 = "",
      Assembly = "", stringsAsFactors = FALSE))
  expect_error(mk_asmb_db(p), "either paired reads or an assembly")
})

test_that("half a read pair still errors, assembly or not", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz", R2 = "",
      Assembly = "S1.fa", stringsAsFactors = FALSE),
    assemblies = "S1.fa", reads = "S1_R1.fq.gz")
  expect_error(mk_asmb_db(p), "Single-end reads are not supported")
})

test_that("a missing assembly file errors and names the sample", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "", R2 = "",
      Assembly = "S1.fa", stringsAsFactors = FALSE))
  expect_error(mk_asmb_db(p), "S1: .*S1\\.fa")
  expect_error(mk_asmb_db(p), "not found")
})

test_that("check_assemblies = FALSE skips the file check", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "", R2 = "",
      Assembly = "S1.fa", stringsAsFactors = FALSE))
  expect_no_error(mk_asmb_db(p, check_assemblies = FALSE))
})

test_that("naming an assembly without an assembly_path errors", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "", R2 = "",
      Assembly = "S1.fa", stringsAsFactors = FALSE),
    assemblies = "S1.fa")
  expect_error(
    new_db(db_path = p$db, mapping_fn = p$mapping_fn, data_path = p$data_dir,
      seeds_db = p$mapping_fn, labels_db = p$mapping_fn),
    "no assembly_path was given")
})

test_that("a mapping with no Assembly column still gets the column", {
  p <- mk_asmb_project(
    data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
      R2 = "S1_R2.fq.gz", stringsAsFactors = FALSE),
    reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  expect_no_error(mk_asmb_db(p))
  con <- DBI::dbConnect(RSQLite::SQLite(), p$db)
  on.exit(DBI::dbDisconnect(con))
  expect_true("Assembly" %in% DBI::dbListFields(con, "samples"))
})
