test_that("write_clade_manifest writes the four columns Nextflow reads", {
  d <- withr::local_tempdir()
  p <- RiboPilot:::write_clade_manifest(file.path(d, "clades.tsv"),
    data.frame(ref_id = "actinopteri", taxid = 186623L, clade_name = "Actinopteri",
      query = '"Actinopteri"[Organism]', stringsAsFactors = FALSE))
  x <- utils::read.delim(p, stringsAsFactors = FALSE)
  expect_equal(names(x), c("ref_id", "taxid", "clade_name", "query"))
  expect_equal(x$taxid, 186623L)
  lines <- readLines(p)
  expect_equal(lines[1], "ref_id\ttaxid\tclade_name\tquery")
  expect_equal(lines[2], 'actinopteri\t186623\tActinopteri\t"Actinopteri"[Organism]')
})

test_that("read_count_files reads every count file and skips junk", {
  d <- withr::local_tempdir()
  writeLines("7898\t11750", file.path(d, "7898.count"))
  writeLines("186623\t9000", file.path(d, "186623.count"))
  writeLines("", file.path(d, "6656.count"))
  x <- RiboPilot:::read_count_files(d)
  expect_equal(nrow(x), 2L)
  expect_equal(x$n_records[x$taxid == 7898L], 11750L)
})

test_that("merge_counts replaces stale rows and keeps the rest", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = c(7898L, 33208L), n_records = c(10L, 500L)),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out <- RiboPilot:::merge_counts(d, data.frame(taxid = 7898L, n_records = 11750L))
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 500L)
  expect_equal(nrow(RiboPilot:::read_counts_cache(d)), 2L)
})

test_that("nf_reference_command names the entry point and every param", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  a <- RiboPilot:::nf_reference_command(d, "COUNT", file.path(d, "clades.tsv"),
    "/refs", file.path(d, "out"))
  expect_equal(a[1], "run")
  expect_true("-entry" %in% a)
  expect_equal(a[which(a == "-entry") + 1], "COUNT")
  expect_equal(a[which(a == "--refDir") + 1], "/refs")
  expect_true("--clades" %in% a && "--outDir" %in% a)
})

test_that("the reference folder is bind-mounted for docker and for singularity", {
  docker <- c("docker {", "  enabled = true", "}")
  sing <- c("singularity {", "  enabled = true", "}")
  expect_equal(RiboPilot:::ref_bind_option(docker, "/refs"),
    c("-process.containerOptions", "-v /refs:/refs"))
  expect_equal(RiboPilot:::ref_bind_option(sing, "/refs"),
    c("-process.containerOptions", "--bind /refs"))
  expect_length(RiboPilot:::ref_bind_option("process { executor = 'local' }", "/refs"), 0L)
})

test_that("refseq and the api key only appear when asked for", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  plain <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o")
  expect_false("--refseq" %in% plain)
  expect_false("--ncbi_api_key" %in% plain)
  keyed <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o",
    api_key = "abc", refseq = TRUE)
  expect_equal(keyed[which(keyed == "--ncbi_api_key") + 1], "abc")
  expect_equal(keyed[which(keyed == "--refseq") + 1], "true")
})

test_that("nf_reference_command bind-mounts refDir when the project config has a real container block", {
  d <- withr::local_tempdir()
  ref_dir <- file.path(d, "refs")
  docker_cfg <- RiboPilot:::container_engine_block("docker")
  writeLines(docker_cfg, file.path(d, ".config"))
  a <- RiboPilot:::nf_reference_command(d, "COUNT", "c.tsv", ref_dir, "o")
  idx <- which(a == "-process.containerOptions")
  expect_length(idx, 1L)
  expect_equal(a[idx + 1], paste0("-v ", ref_dir, ":", ref_dir))

  sing_cfg <- RiboPilot:::container_engine_block("singularity")
  writeLines(sing_cfg, file.path(d, ".config"))
  b <- RiboPilot:::nf_reference_command(d, "COUNT", "c.tsv", ref_dir, "o")
  idx2 <- which(b == "-process.containerOptions")
  expect_length(idx2, 1L)
  expect_equal(b[idx2 + 1], paste0("--bind ", ref_dir))
})
