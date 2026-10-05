# Mapping-file and project-path checks at setup: common spreadsheet mistakes
# must stop early with a message that names the problem and the fix.

# Writes `txt` byte for byte (so CRLF and a BOM survive) and returns the path.
write_map <- function(txt, name = "mapping.csv", bom = FALSE, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  fn <- file.path(dir, name)
  con <- file(fn, "wb")
  if (bom) writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeBin(charToRaw(txt), con)
  close(con)
  fn
}

test_that("file-level mistakes are refused with the fix named", {
  expect_error(RiboPilot:::read_mapping(file.path(tempdir(), "nope.csv")), "Mapping file not found")
  expect_error(RiboPilot:::read_mapping(write_map("")), "is empty")
  expect_error(RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2\n")), "no samples")
  expect_error(RiboPilot:::read_mapping(write_map("ID;Taxon;R1;R2\nS1;G;a;b\n")), "semicolons")
  expect_error(RiboPilot:::read_mapping(write_map("ID\tTaxon\tR1\tR2\nS1\tG\ta\tb\n")), "tabs")
  xlsx <- write_map("", name = "mapping.xlsx")
  writeBin(as.raw(c(0x50, 0x4b, 0x03, 0x04, 0x14, 0x00, 0x06, 0x00)), xlsx)
  expect_error(RiboPilot:::read_mapping(xlsx), "not a plain-text CSV")
})

test_that("values are kept as typed, trimmed, and blank rows dropped", {
  m <- RiboPilot:::read_mapping(write_map(
    "ID,Taxon,R1,R2\n001,Genus sp,a.fq,b.fq\n\" 002 \", Genus sp , c.fq ,d.fq\n,,,\n"))
  expect_equal(m$ID, c("001", "002"))
  expect_equal(m$Taxon, c("Genus sp", "Genus sp"))
  expect_equal(m$R1, c("a.fq", "c.fq"))
  # Excel's CSV UTF-8 adds a BOM, and Windows line ends are CRLF
  m <- RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2\r\nS1,G,a,b\r\n", bom = TRUE))
  expect_equal(names(m), c("ID", "Taxon", "R1", "R2"))
  expect_equal(m$R2, "b")
})

test_that("column names match ignoring case, and a missing one lists what was found", {
  map <- write_map("id,taxon,r1,r2\nS1,G,a,b\n")
  msgs <- testthat::capture_messages(m <- RiboPilot:::read_mapping(map))
  expect_true(any(grepl("Using mapping column 'taxon' as 'Taxon'", msgs, fixed = TRUE)))
  expect_equal(names(m), c("ID", "Taxon", "R1", "R2"))
  err <- expect_error(RiboPilot:::read_mapping(write_map("Sample,Species,R1,R2\nS1,G,a,b\n")))
  expect_match(conditionMessage(err), "missing column(s): ID, Taxon", fixed = TRUE)
  expect_match(conditionMessage(err), "Columns found: Sample, Species, R1, R2", fixed = TRUE)
  # an assembly-only mapping needs no read columns
  m <- RiboPilot:::read_mapping(write_map("ID,Taxon,Assembly\nA1,G,a.fa\n"))
  expect_true(all(is.na(m$R1)) && all(is.na(m$R2)))
})

test_that("columns RiboPilot fills in, or that differ only by case, are refused", {
  expect_error(RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2,TaxID,lineage\nS1,G,a,b,1,x\n")),
    "fills in itself: TaxID, lineage")
  expect_error(RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2,Site,site\nS1,G,a,b,x,y\n")),
    "more than upper/lower case")
})

test_that("the reserved list covers every sample column new_db adds", {
  skip_if_not_installed("withr")
  map <- write_map("ID,Taxon,R1,R2\nS1,Danio rerio,a.fq,b.fq\n")
  fa <- file.path(dirname(map), "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  db <- file.path(dirname(map), ".sqlite")
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE, seeds_db = fa, labels_db = fa)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  added <- setdiff(DBI::dbListFields(con, "samples"), c("ID", "Taxon", "R1", "R2", "Assembly", RiboPilot:::META_SAMPLE_COLS))
  expect_setequal(added, RiboPilot:::SAMPLE_BUILTIN_COLS)
})

test_that("the vignette's reserved names are the reserved list", {
  rmd <- test_path("..", "..", "vignettes", "Your-Own-Project.Rmd")
  skip_if_not(file.exists(rmd), "no vignette source")
  txt <- paste(readLines(rmd, warn = FALSE), collapse = " ")
  txt <- sub("\\.\\s.*$", "", sub("^.*use these names:", "", txt))
  expect_setequal(trimws(strsplit(gsub("\\band\\b", "", txt), ",")[[1]]),
    RiboPilot:::SAMPLE_BUILTIN_COLS)
})

test_that("unsafe, repeated, or missing IDs are refused and shown", {
  expect_error(RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2\nCHIM00,G,a,b\nchim00,G,c,d\n")),
    "Repeated: CHIM00, chim00")
  err <- expect_error(RiboPilot:::read_mapping(write_map(paste0(
    "ID,Taxon,R1,R2\nCHIM 00,G,a,b\nUSNM.1,G,c,d\nA/B,G,e,f\nCHIM\u00e900,G,g,h\n"))))
  for (id in c("'CHIM 00'", "'USNM.1'", "'A/B'", "'CHIM\u00e900'")) {
    expect_match(conditionMessage(err), id, fixed = TRUE)
  }
  expect_error(RiboPilot:::read_mapping(write_map("ID,Taxon,R1,R2\nS1,G,a,b\n,G,c,d\nNA,G,e,f\n")),
    "line\\(s\\) 3, 4")
  expect_error(RiboPilot:::read_mapping(write_map(paste0("ID,Taxon,R1,R2\n", strrep("A", 41), ",G,a,b\n"))),
    "40 characters")
})

test_that("unsafe characters in paths are named, one line per path", {
  err <- expect_error(RiboPilot:::stop_unsafe_paths(c(data_path = "/data/My Reads (1)",
    assembly_path = "/data/ok", ref = "https://example.org/a b")))
  expect_match(conditionMessage(err), "data_path: '/data/My Reads (1)' (remove: space, '(', ')')", fixed = TRUE)
  expect_false(grepl("assembly_path|example.org", conditionMessage(err)))
  # a colon is fine in a file name, not in a folder that gets mounted
  expect_no_error(RiboPilot:::stop_unsafe_paths(c("S1 R1" = "S1:x_R1.fq.gz"), mount = FALSE))
  expect_error(RiboPilot:::stop_unsafe_paths(c(data_path = "/data/a:b")), "':'")
})

# a project whose mapping names `reads` under data/
mk_setup <- function(mapping, reads = character(0), env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  data_dir <- file.path(dir, "data")
  dir.create(data_dir)
  for (r in reads) {
    dir.create(dirname(file.path(data_dir, r)), showWarnings = FALSE, recursive = TRUE)
    file.create(file.path(data_dir, r))
  }
  map_fn <- file.path(dir, "mapping.csv")
  utils::write.csv(mapping, map_fn, row.names = FALSE)
  list(dir = dir, data_dir = data_dir, mapping_fn = map_fn)
}

np <- function(p, proj, ...) {
  suppressWarnings(suppressMessages(new_project(path = proj, mapping_fn = p$mapping_fn,
    executor = "conda", Rproj = FALSE, custom_seeds_db = p$mapping_fn,
    custom_labels_db = p$mapping_fn, ...)))
}

test_that("a failed check leaves nothing behind, and a rerun can succeed", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Species = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  expect_error(np(p, proj, data_path = p$data_dir), "missing column\\(s\\): Taxon")
  expect_false(dir.exists(proj))
  p2 <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  expect_error(suppressMessages(new_project(path = proj, mapping_fn = p2$mapping_fn,
    data_path = p2$data_dir, executor = "dokcer", Rproj = FALSE)), "No config found for executor 'dokcer'")
  expect_false(dir.exists(proj))
  np(p2, proj, data_path = p2$data_dir)
  expect_true(file.exists(file.path(proj, ".sqlite")))
})

test_that("an existing project is refused and its mapping.csv left alone", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)
  before <- tools::md5sum(file.path(proj, "mapping.csv"))
  p2 <- mk_setup(data.frame(ID = "S2", Taxon = "Danio rerio", R1 = "S2_R1.fq.gz",
    R2 = "S2_R2.fq.gz"), reads = c("S2_R1.fq.gz", "S2_R2.fq.gz"))
  expect_error(np(p2, proj, data_path = p2$data_dir), "already exists")
  expect_equal(tools::md5sum(file.path(proj, "mapping.csv")), before)
})

test_that("setup refuses paths with spaces and reads with no data_path", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  spaced <- file.path(p$dir, "my reads")
  dir.create(spaced)
  file.copy(file.path(p$data_dir, c("S1_R1.fq.gz", "S1_R2.fq.gz")), spaced)
  expect_error(np(p, file.path(p$dir, "proj"), data_path = spaced), "data_path: .*remove: space")
  expect_error(np(p, file.path(p$dir, "my proj"), data_path = p$data_dir), "path: .*remove: space")
  expect_error(np(p, file.path(p$dir, "proj")), "no data_path was given")
})

test_that("read pairs, full paths, and a missing data_path folder are caught", {
  p <- mk_setup(data.frame(ID = c("S1", "S2"), Taxon = "Danio rerio",
    R1 = c("S1_R1.fq.gz", "a/x.fq.gz"), R2 = c("S1_R1.fq.gz", "b/x.fq.gz")),
    reads = c("S1_R1.fq.gz", "a/x.fq.gz", "b/x.fq.gz"))
  m <- RiboPilot:::read_mapping(p$mapping_fn)
  err <- expect_error(RiboPilot:::check_mapping(m, data_path = p$data_dir))
  expect_match(conditionMessage(err), "two different read files")
  expect_match(conditionMessage(err), "S1: S1_R1.fq.gz", fixed = TRUE)
  expect_match(conditionMessage(err), "S2: a/x.fq.gz", fixed = TRUE)
  m$R2 <- c("S1_R2.fq.gz", "b/y.fq.gz")
  m$R1[1] <- file.path(p$data_dir, "S1_R1.fq.gz")
  expect_error(RiboPilot:::check_mapping(m, data_path = p$data_dir), "not full paths")
  expect_error(RiboPilot:::check_mapping(m, data_path = file.path(p$dir, "nope")),
    "data_path does not exist or is not a folder")
  m$R1[1] <- "S1 R1.fq.gz"
  expect_error(RiboPilot:::check_mapping(m, data_path = p$data_dir), "S1 R1: 'S1 R1.fq.gz'")
})

test_that("assemblies must be FASTA and all gzipped or all plain", {
  dir <- withr::local_tempdir()
  writeLines(c("@r1", "ACGT", "+", "IIII"), file.path(dir, "A1.fq"))
  writeLines(c(">c1", "ACGT"), file.path(dir, "A2.fa"))
  gz <- gzfile(file.path(dir, "A3.fa.gz"), "w"); writeLines(c(">c1", "ACGT"), gz); close(gz)
  m <- data.frame(ID = c("A1", "A2"), Taxon = "G", R1 = NA_character_, R2 = NA_character_,
    Assembly = c("A1.fq", "A2.fa"))
  expect_error(RiboPilot:::check_mapping(m, assembly_path = dir), "A1: .*A1\\.fq")
  m$Assembly <- c("A3.fa.gz", "A2.fa")
  expect_error(RiboPilot:::check_mapping(m, assembly_path = dir), "all gzipped")
  m <- m[1, ]
  expect_no_error(RiboPilot:::check_mapping(m, assembly_path = dir))
})

# add_samples() / update_sample_metadata() ----

test_that("add_samples appends new samples, seeded the same as new_project", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)

  file.create(file.path(p$data_dir, c("S2_R1.fq.gz", "S2_R2.fq.gz")))
  map2 <- file.path(p$dir, "mapping2.csv")
  utils::write.csv(data.frame(ID = "S2", Taxon = "Danio rerio", site = "Reef A",
    R1 = "S2_R1.fq.gz", R2 = "S2_R2.fq.gz"), map2, row.names = FALSE)
  add_samples(project = proj, mapping_fn = map2, data_path = p$data_dir)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  samples <- DBI::dbReadTable(con, "samples")
  preprocess <- DBI::dbReadTable(con, "preprocess")
  assemble <- DBI::dbReadTable(con, "assemble")
  annotate <- DBI::dbReadTable(con, "annotate")

  expect_setequal(samples$ID, c("S1", "S2"))
  expect_equal(samples$site[samples$ID == "S2"], "Reef A")
  expect_true(is.na(samples$site[samples$ID == "S1"]))
  expect_setequal(preprocess$ID, c("S1", "S2"))
  expect_setequal(assemble$ID, c("S1", "S2"))
  expect_setequal(annotate$ID, c("S1", "S2"))
  # S2 is Ready the same way S1 (created by new_project) is Ready
  s2_asmb <- assemble[assemble$ID == "S2", ]
  expect_equal(s2_asmb$assemble_switch, 1)
  expect_equal(s2_asmb$assemble_lock, 0)
  s2_ann <- annotate[annotate$ID == "S2", ]
  expect_equal(s2_ann$annotate_switch, 1)
  expect_equal(s2_ann$annotate_lock, 0)
  expect_equal(s2_ann$reviewed, "no")
})

test_that("add_samples only takes files from the folders the project's config uses", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)
  other <- file.path(p$dir, "data2")
  dir.create(other)
  file.create(file.path(other, c("S2_R1.fq.gz", "S2_R2.fq.gz")))
  writeLines(c(">c1", "ACGT"), file.path(other, "S3.fa"))
  file.create(file.path(p$data_dir, c("S2_R1.fq.gz", "S2_R2.fq.gz")))
  map2 <- file.path(p$dir, "mapping2.csv")
  utils::write.csv(data.frame(ID = "S2", Taxon = "Danio rerio", R1 = "S2_R1.fq.gz",
    R2 = "S2_R2.fq.gz"), map2, row.names = FALSE)
  expect_error(add_samples(project = proj, mapping_fn = map2, data_path = other),
    "must be in this project's read folder")
  map3 <- file.path(p$dir, "mapping3.csv")
  utils::write.csv(data.frame(ID = "S3", Taxon = "Danio rerio", Assembly = "S3.fa"),
    map3, row.names = FALSE)
  expect_error(add_samples(project = proj, mapping_fn = map3, assembly_path = other),
    "no assembly folder")
  # the project's read folder is the default
  expect_message(add_samples(project = proj, mapping_fn = map2), "Added 1 sample")
})

test_that("add_samples refuses IDs already present, and bad IDs", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)

  dup_map <- file.path(p$dir, "dup.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), dup_map, row.names = FALSE)
  expect_error(add_samples(project = proj, mapping_fn = dup_map, data_path = p$data_dir),
    "Already in the project.*S1")

  file.create(file.path(p$data_dir, c("bad_R1.fq.gz", "bad_R2.fq.gz")))
  bad_map <- file.path(p$dir, "bad.csv")
  utils::write.csv(data.frame(ID = "S 2", Taxon = "Danio rerio", R1 = "bad_R1.fq.gz",
    R2 = "bad_R2.fq.gz"), bad_map, row.names = FALSE)
  expect_error(add_samples(project = proj, mapping_fn = bad_map, data_path = p$data_dir),
    "IDs may use only letters")

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbReadTable(con, "samples")$ID, "S1")
})

test_that("update_sample_metadata updates metadata but not pipeline state", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 0, assemble_lock = 1 WHERE ID = 'S1'")
  DBI::dbDisconnect(con)

  meta_map <- file.path(p$dir, "meta.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio nov", site = "Reef B",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), meta_map, row.names = FALSE)
  update_sample_metadata(project = proj, mapping_fn = meta_map, data_path = p$data_dir)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  samples <- DBI::dbReadTable(con, "samples")
  assemble <- DBI::dbReadTable(con, "assemble")
  expect_equal(samples$Taxon[samples$ID == "S1"], "Danio rerio nov")
  expect_equal(samples$site[samples$ID == "S1"], "Reef B")
  expect_equal(samples$R1[samples$ID == "S1"], "S1_R1.fq.gz")
  # pipeline state and locks are untouched
  expect_equal(assemble$assemble_switch[assemble$ID == "S1"], 0)
  expect_equal(assemble$assemble_lock[assemble$ID == "S1"], 1)
  bk <- file.path(proj, ".old_sqlite_dbs", ".sqlite.1")
  expect_true(file.exists(bk))
  old <- DBI::dbConnect(RSQLite::SQLite(), bk)
  expect_equal(DBI::dbReadTable(old, "samples")$Taxon, "Danio rerio")
  DBI::dbDisconnect(old)
})

test_that("update_sample_metadata writes to the stored ID when the case differs", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)
  meta_map <- file.path(p$dir, "meta.csv")
  utils::write.csv(data.frame(ID = "s1", Taxon = "Danio rerio", site = "Reef B",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), meta_map, row.names = FALSE)
  suppressMessages(update_sample_metadata(project = proj, mapping_fn = meta_map, data_path = p$data_dir))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  samples <- DBI::dbReadTable(con, "samples")
  expect_equal(samples$ID, "S1")
  expect_equal(samples$site, "Reef B")
})

test_that("a changed Taxon clears, unlocks, and re-queues that sample's taxonomy", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = c("S1", "S2"), Taxon = "Danio rerio", R1 = c("S1_R1.fq.gz", "S2_R1.fq.gz"),
    R2 = c("S1_R2.fq.gz", "S2_R2.fq.gz")),
    reads = c("S1_R1.fq.gz", "S1_R2.fq.gz", "S2_R1.fq.gz", "S2_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET taxid = 7955, sci_name = 'Danio rerio', rank = 'species',
    match_type = 'exact', lineage = 'x', lineage_taxids = '1', taxonomy_switch = 2, taxonomy_lock = 1")
  meta_map <- file.path(p$dir, "meta.csv")
  utils::write.csv(data.frame(ID = c("S1", "S2"), Taxon = c("Danio aesculapii", "Danio rerio"),
    R1 = c("S1_R1.fq.gz", "S2_R1.fq.gz"), R2 = c("S1_R2.fq.gz", "S2_R2.fq.gz")), meta_map, row.names = FALSE)
  expect_message(update_sample_metadata(project = proj, mapping_fn = meta_map, data_path = p$data_dir),
    "Taxon changed for S1")
  s <- DBI::dbGetQuery(con, "SELECT * FROM samples ORDER BY ID")
  expect_true(is.na(s$sci_name[1]) && is.na(s$taxid[1]) && is.na(s$lineage[1]))
  expect_equal(c(s$taxonomy_switch[1], s$taxonomy_lock[1]), c(1, 0))
  expect_equal(s$sci_name[2], "Danio rerio")
  expect_equal(c(s$taxonomy_switch[2], s$taxonomy_lock[2]), c(2, 1))
})

test_that("update_sample_metadata refuses an ID not in the project", {
  skip_on_os("windows")
  p <- mk_setup(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), reads = c("S1_R1.fq.gz", "S1_R2.fq.gz"))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir)

  meta_map <- file.path(p$dir, "meta.csv")
  utils::write.csv(data.frame(ID = "S9", Taxon = "Danio rerio", R1 = "S1_R1.fq.gz",
    R2 = "S1_R2.fq.gz"), meta_map, row.names = FALSE)
  expect_error(update_sample_metadata(project = proj, mapping_fn = meta_map, data_path = p$data_dir),
    "Not in the project.*S9")
})

test_that("GEOME_BCID, GBIF_ID, and BioSample columns are accepted, matched ignoring case, and stored normalized", {
  skip_on_os("windows")
  local_mocked_bindings(.geome_get = geome_fixture_get, .ncbi_get = function(...) stop("x"),
    url_status = function(u) "ok")
  p <- mk_setup(data.frame(ID = c("S1", "S2"), Taxon = "Danio rerio", R1 = c("a_1.fq", "b_1.fq"),
    R2 = c("a_2.fq", "b_2.fq"), geome_bcid = c("https://n2t.net/ark:/21547/CYB2REEDY", NA),
    Biosample = c(NA, "samn00000001")), reads = c("a_1.fq", "b_1.fq", "a_2.fq", "b_2.fq"))
  m <- RiboPilot:::read_mapping(p$mapping_fn)
  expect_true(all(c("GEOME_BCID", "BioSample") %in% names(m)))
  expect_false(any(c("geome_bcid", "Biosample") %in% names(m)))
  proj <- file.path(p$dir, "proj")
  np(p, proj, data_path = p$data_dir, fetch_biosample = FALSE)
  s <- meta_q(proj, "SELECT ID, GEOME_BCID, BioSample FROM samples ORDER BY ID")
  expect_equal(s$GEOME_BCID, c("ark:/21547/CYB2REEDY", NA))
  expect_equal(s$BioSample, c(NA, "SAMN00000001"))
})
