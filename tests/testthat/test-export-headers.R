# Templated FASTA headers at export time: defaults, custom mapping-file columns,
# the checks that block a write, and the new folder layout.

# S1 has two exportable units, S2 has one; `collection_site` stands in for an
# extra column of the user's mapping file.
mk_hdr_project <- function() {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT,
                       collection_site TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT,
                       length INT, topology TEXT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT,
                       type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_lock INT,
                       problematic TEXT, partial TEXT, reviewed TEXT, chimeric INT,
                       n_features INT, missing TEXT, export_group TEXT)")
  seq <- paste0(strrep("A", 200), strrep("C", 100), strrep("G", 100))
  for (u in list(c("S1", 1, 1), c("S1", 1, 2), c("S2", 1, 1))) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES(?,?,?,?,?,'linear',0)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3]), seq, nchar(seq)))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES(?,?,?,1,'no','no','yes',0,3,'',NULL)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'18S','rRNA',1,150)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'ITS1','misc_RNA',151,200)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'5.8S','rRNA',201,250)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'ITS2','misc_RNA',251,300)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES(?,?,?,'28S','rRNA',301,400)",
      params = list(u[1], as.integer(u[2]), as.integer(u[3])))
  }
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Genus speciesA',NULL,NULL,'Reef A')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S2','Genus speciesB',NULL,NULL,'Reef B')")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1),('S2',1)")
  dir
}

hdr_bundle <- function(dir, ...) {
  suppressWarnings(suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"),
    report = FALSE, zip = FALSE, ...)))
}

deflines <- function(gdir, label) {
  x <- readLines(file.path(gdir, paste0(label, ".fasta")))
  x[grepl("^>", x)]
}

test_that("the default templates resolve for every record", {
  dir <- mk_hdr_project()
  g <- hdr_bundle(dir, label = "dflt")
  d <- deflines(g, "dflt")
  expect_length(d, 3)
  expect_true(all(grepl("[organism=Genus species", d, fixed = TRUE)))
  expect_true(all(grepl("[moltype=genomic DNA]", d, fixed = TRUE)))
  expect_true(all(grepl("ribosomal RNA gene", d)))
  expect_true(any(grepl("^>S2 ", d)))
  expect_true(any(grepl("^>S1_p1_s1 ", d)))
  # marker deflines come from the marker template
  its <- readLines(file.path(g, "markers", "ITS", "all_ITS.fasta"))
  expect_true(any(grepl("_ITS \\[organism=", its)))
  expect_true(any(grepl("internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2",
    its, fixed = TRUE)))
})

test_that("an extra mapping-file column is usable as a header token", {
  dir <- mk_hdr_project()
  g <- hdr_bundle(dir, label = "site",
    fasta_header = "{seqid} site={collection_site} {Taxon}")
  d <- deflines(g, "site")
  expect_true(any(grepl(">S2 site=Reef B Genus speciesB", d, fixed = TRUE)))
  expect_true(all(grepl("site=Reef ", d, fixed = TRUE)))
})

test_that("the template is the whole defline: no metadata is appended", {
  dir <- mk_hdr_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN specimen_voucher TEXT")
  DBI::dbExecute(con, "UPDATE samples SET specimen_voucher = 'USNM:1' WHERE ID = 'S2'")
  DBI::dbDisconnect(con)
  g <- hdr_bundle(dir, label = "mods", fasta_header = "{seqid}")
  d <- deflines(g, "mods")
  expect_true(">S2" %in% d)
  expect_false(any(grepl("specimen_voucher", d, fixed = TRUE)))
})

test_that("a template without {seqid} is rejected before anything is written", {
  dir <- mk_hdr_project()
  expect_error(hdr_bundle(dir, label = "noid", fasta_header = "{Taxon}"),
    "seqid", ignore.case = TRUE)
  expect_false(dir.exists(file.path(dir, "out", "export", "noid")))
})

test_that("an unknown column and an unmatched brace are rejected", {
  dir <- mk_hdr_project()
  expect_error(hdr_bundle(dir, label = "badcol", fasta_header = "{seqid} {not_a_column}"),
    "not_a_column")
  expect_error(hdr_bundle(dir, label = "brace", fasta_header = "{seqid} {Taxon"),
    "Unmatched")
})

test_that("duplicate headers are rejected", {
  # two one-unit samples of the same taxon collide once {seqid} is dropped
  dir <- mk_hdr_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE ID='S1' AND scaffold=2")
  DBI::dbExecute(con, "UPDATE samples SET Taxon = 'Genus speciesA'")
  DBI::dbDisconnect(con)
  expect_error(hdr_bundle(dir, label = "dups", fasta_header = "{Taxon}"),
    "duplicated header")
  expect_false(dir.exists(file.path(dir, "out", "export", "dups")))
})

test_that("the bundle uses the samples/ layout and README.txt", {
  dir <- mk_hdr_project()
  g <- hdr_bundle(dir, label = "layout")
  expect_setequal(basename(Sys.glob(file.path(g, "samples", "*.fasta"))),
    c("S1_p1_s1.fasta", "S1_p1_s2.fasta", "S2.fasta"))
  expect_setequal(basename(Sys.glob(file.path(g, "samples", "*.tbl"))),
    c("S1_p1_s1.tbl", "S1_p1_s2.tbl", "S2.tbl"))
  # per-unit files no longer sit in the group root
  expect_length(Sys.glob(file.path(g, "*.fsa")), 0)
  expect_setequal(basename(Sys.glob(file.path(g, "*.tbl"))), "layout.tbl")
  expect_true(file.exists(file.path(g, "README.txt")))
  expect_false(file.exists(file.path(g, "README_SUBMISSION.txt")))
  expect_true(file.exists(file.path(g, "layout.fasta")))
  expect_true(file.exists(file.path(g, "layout_summary.csv")))
  expect_match(paste(readLines(file.path(g, "README.txt")), collapse = "\n"), "samples/")
})

test_that("a saved configuration is used when named", {
  dir <- mk_hdr_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  RiboPilot:::set_export_opts(con, "{seqid} custom {Taxon}", "{seqid}_{element} custom",
    name = "mine")
  DBI::dbDisconnect(con)
  g <- hdr_bundle(dir, label = "named", export_opts = "mine")
  expect_true(all(grepl(" custom ", deflines(g, "named"), fixed = TRUE)))
  # an explicit argument still wins over the saved configuration
  g2 <- hdr_bundle(dir, label = "override", export_opts = "mine",
    fasta_header = "{seqid} explicit")
  expect_true(all(grepl(" explicit", deflines(g2, "override"), fixed = TRUE)))
})

test_that("multi-unit samples get an automatic [isolate=seqid]; single-unit ones do not", {
  dir <- mk_hdr_project()
  d <- deflines(hdr_bundle(dir, label = "iso"), "iso")
  expect_true(any(grepl("^>S1_p1_s1 .* \\[isolate=S1_p1_s1\\]$", d)))
  expect_true(any(grepl("^>S1_p1_s2 .* \\[isolate=S1_p1_s2\\]$", d)))
  expect_false(any(grepl("^>S2 .*isolate", d)))
  # the marker files carry it too: same organism + voucher, several 18S records
  m <- readLines(file.path(dir, "out", "export", "iso", "markers", "18S", "all_18S.fasta"))
  expect_true(any(grepl("^>S1_p1_s1_18S .*\\[isolate=S1_p1_s1\\]$", m)))
})

test_that("[isolate=] in a template is refused only where RiboPilot sets the isolate", {
  dir <- mk_hdr_project()
  expect_error(hdr_bundle(dir, label = "iso2", fasta_header = "{seqid} [isolate={ID}] {Taxon}"),
    "isolate")
  multi <- data.frame(ID = c("S1", "S1"), seqid = c("S1_p1_s1", "S1_p1_s2"))
  single <- data.frame(ID = "S2", seqid = "S2")
  expect_equal(validate_fasta_header("{seqid} [Isolate=x]", multi)$level, "error")
  expect_equal(validate_fasta_header("{seqid} [isolate={ID}]", single)$level, "ok")
  expect_equal(validate_fasta_header("{seqid} [note=isolated from gut]", multi)$level, "ok")
})

test_that("an empty source modifier on any record is a warning, not only on the first", {
  dir <- mk_hdr_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN voucher TEXT")
  DBI::dbExecute(con, "UPDATE samples SET voucher = 'USNM:1' WHERE ID = 'S1'")
  DBI::dbDisconnect(con)
  expect_warning(suppressMessages(export_bundle(dir, out_dir = file.path(dir, "out"), report = FALSE,
    zip = FALSE, itsx_check = FALSE, label = "vou", fasta_header = "{seqid} [specimen-voucher={voucher}]")),
    "Empty \\[specimen-voucher=\\] on 1 of 3 record\\(s\\) \\(S2\\)")
  d <- data.frame(ID = "S1", seqid = "S1", t = NA)
  expect_equal(validate_fasta_header("{seqid} [note={t}]", d)$level, "warn")
})

test_that("a header template is filled by lookup and never runs code", {
  dir <- mk_hdr_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  u <- data.frame(ID = "S2", path = 1L, scaffold = 1L, seqid = "S2")
  d <- RiboPilot:::export_header_data(con, u)
  f <- tempfile("rp_hdr_")
  tpl <- sprintf('{seqid} {file.create("%s")}', f)
  expect_equal(RiboPilot:::validate_fasta_header(tpl, d)$level, "error")
  expect_error(RiboPilot:::render_fasta_headers(tpl, d), "is not in the database")
  expect_false(file.exists(f))
  expect_equal(RiboPilot:::render_fasta_headers("{seqid} [organism={Taxon}] site={collection_site}", d),
    "S2 [organism=Genus speciesB] site=Reef B")
})
