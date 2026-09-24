# Cross-kingdom homonyms and failed lineage lookups, against a real (tiny)
# taxdump built by helper-taxdump.R. No network, no external binaries.

test_that("a cross-kingdom homonym is flagged, not silently resolved", {
  res <- RiboPilot:::resolve_taxonomy("Morus", local_taxdump())
  expect_equal(res$match_type, "ambiguous")
  expect_false(is.na(res$alternatives))
  expect_true(grepl("3497", res$alternatives, fixed = TRUE))
  expect_true(grepl("37577", res$alternatives, fixed = TRUE))
  # the two kingdoms are what tells the user which Morus they were handed
  expect_true(grepl("Viridiplantae", res$alternatives, fixed = TRUE))
  expect_true(grepl("Metazoa", res$alternatives, fixed = TRUE))
  expect_true(is.na(res$taxid))
  expect_true(is.na(res$lineage_taxids))
})

test_that("an unambiguous name is still an exact match", {
  res <- RiboPilot:::resolve_taxonomy("Danio rerio", local_taxdump())
  expect_equal(res$match_type, "exact")
  expect_equal(res$taxid, 7955L)
  expect_true(is.na(res$alternatives))
  expect_equal(res$lineage, "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio")
  expect_equal(res$lineage_taxids, "33208;7711;186623;7952;2743709;7954;7955")
})

test_that("td_name2taxid keeps every scientific-name candidate and drops other taxa's synonyms", {
  dir <- local_taxdump()
  cand <- RiboPilot:::with_tax_index(dir, function(con)
    RiboPilot:::td_name2taxid(con, c("Morus", "Danio rerio", "danio RERIO")))
  expect_equal(cand$Morus, c(3497L, 37577L))
  expect_equal(cand[["Danio rerio"]], 7955L)
  expect_equal(cand[["danio RERIO"]], 7955L)
})

test_that("a mistyped species in a homonym genus searches every genus", {
  res <- RiboPilot:::resolve_taxonomy("Morus bassanu", local_taxdump())
  expect_equal(res$match_type, "fuzzy_species")
  expect_equal(res$taxid, 37578L)
})

test_that("a homonym genus with no close species stays unresolved", {
  res <- RiboPilot:::resolve_taxonomy("Morus zzzzqq", local_taxdump())
  expect_equal(res$match_type, "ambiguous")
  expect_true(is.na(res$taxid))
  expect_true(grepl("Metazoa", res$alternatives, fixed = TRUE))
})

test_that("search results name the kingdom that tells two homonyms apart", {
  dir <- local_taxdump()
  r <- RiboPilot:::tax_search(dir, "Morus", 10)
  expect_setequal(r$kingdom[r$name == "Morus"], c("Viridiplantae", "Metazoa"))
  expect_true("Morus   [genus, Metazoa]   (taxid 37577)" %in% names(RiboPilot:::tax_search_choices(r)))
})

test_that("synonyms and common names resolve, case-insensitively", {
  dir <- local_taxdump()
  res <- RiboPilot:::resolve_taxonomy(c("Brachydanio rerio", "danio rerio"), dir)
  expect_equal(res$taxid, c(7955L, 7955L))
  expect_equal(res$match_type, c("exact", "exact"))
})

test_that("a taxid that falls out of a batch is reported, not dropped in silence", {
  dir <- local_taxdump()
  expect_warning(
    info <- RiboPilot:::tk_taxid_info(c(3497L, 999999L), dir),
    "999999")
  expect_equal(info$taxid, 3497L)
})

test_that("a batch that returns nothing warns and yields an empty frame", {
  dir <- local_taxdump()
  expect_warning(info <- RiboPilot:::tk_taxid_info(c(999998L, 999999L), dir), "99999")
  expect_equal(nrow(info), 0L)
  expect_equal(names(info), c("taxid", "sci_name", "rank", "lineage", "lineage_taxids"))
})

test_that("a retired taxid resolves through merged.dmp", {
  info <- RiboPilot:::tk_taxid_info(29144L, local_taxdump())
  expect_equal(info$sci_name, "Danio rerio")
})

test_that("a missing rank is filled, keeping the seven slots aligned", {
  # a genus has no species rank: the slot must be filled, not collapsed away
  info <- RiboPilot:::tk_taxid_info(7954L, local_taxdump())
  expect_equal(RiboPilot:::parse_rank_names(info$lineage)[["genus"]], "Danio")
  expect_true(is.na(RiboPilot:::parse_rank_names(info$lineage)[["species"]]))
  expect_true(is.na(RiboPilot:::parse_rank_taxids(info$lineage_taxids)[["species"]]))
  expect_equal(length(strsplit(info$lineage, ";")[[1]]), 7L)
})

test_that("a mistyped epithet falls back to a fuzzy species match", {
  res <- RiboPilot:::resolve_taxonomy("Danio rerlo", local_taxdump())
  expect_equal(res$match_type, "fuzzy_species")
  expect_equal(res$taxid, 7955L)
})

test_that("an unrecognisable name stays unresolved", {
  res <- RiboPilot:::resolve_taxonomy("Zzzzyx qqqqq", local_taxdump())
  expect_equal(res$match_type, "unresolved")
  expect_true(is.na(res$taxid))
})

test_that("the rank-aligned lineage and full ancestor path agree", {
  dir <- local_taxdump()
  lin <- RiboPilot:::tk_ranked_lineage(7955L, dir)
  expect_equal(lin$lineage_taxids, "33208;7711;186623;7952;2743709;7954;7955")
  path <- RiboPilot:::tk_full_lineage_taxids(7955L, dir)
  expect_equal(path[["7955"]], c(1L, 131567L, 2759L, 33208L, 7711L, 186623L, 7952L,
    2743709L, 7954L, 7955L))
})

test_that("the name search index finds a genus and its species", {
  r <- RiboPilot:::tax_search(local_taxdump(), "Danio", 10)
  expect_true("Danio" %in% r$name)
  expect_true("Danio rerio" %in% r$name)
  expect_equal(r$rank[r$name == "Danio"], "genus")
})

test_that("an index with too few taxa is not trusted", {
  withr::local_options(ribopilot.tax_min_nodes = 1000L)
  expect_warning(idx <- RiboPilot:::build_tax_index(local_taxdump()), "incomplete")
  expect_null(idx)
})

test_that("a truncated taxdump tarball is refused, not half-indexed", {
  dir <- withr::local_tempdir()
  src <- local_taxdump()
  withr::with_dir(src, utils::tar(file.path(dir, "full.tar.gz"),
    files = c("names.dmp", "nodes.dmp", "merged.dmp"), compression = "gzip"))
  raw <- readBin(file.path(dir, "full.tar.gz"), "raw", file.info(file.path(dir, "full.tar.gz"))$size)
  writeBin(raw[seq_len(length(raw) %/% 2)], RiboPilot:::taxdump_tarball(dir))
  expect_error(RiboPilot:::extract_taxdump(dir), "incomplete or corrupt")
  expect_error(RiboPilot:::resolve_taxonomy("Danio rerio", dir), "Could not read the NCBI taxonomy")
})

test_that("a taxdump download that fails NCBI's checksum is not kept", {
  dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(download_or_stop = function(url, dest, what, ...) {
    writeLines(if (endsWith(url, ".md5")) paste0(strrep("0", 32), "  taxdump.tar.gz") else "partial", dest)
    dest
  }, .package = "RiboPilot")
  expect_error(RiboPilot:::ensure_taxdump(dir), "checksum")
  expect_false(file.exists(RiboPilot:::taxdump_tarball(dir)))
  expect_length(list.files(dir, pattern = "part"), 0L)
})

test_that("a verified taxdump replaces the old index", {
  dir <- withr::local_tempdir()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  testthat::local_mocked_bindings(download_or_stop = function(url, dest, what, ...) {
    if (endsWith(url, ".md5")) {
      writeLines(paste0(unname(tools::md5sum(sub("\\.md5$", "", dest))), "  taxdump.tar.gz"), dest)
    } else {
      writeLines("whole", dest)
    }
    dest
  }, .package = "RiboPilot")
  file.create(RiboPilot:::tax_index_path(dir))
  RiboPilot:::ensure_taxdump(dir, force = TRUE)
  expect_true(file.exists(RiboPilot:::taxdump_tarball(dir)))
  expect_false(file.exists(RiboPilot:::tax_index_path(dir)))
  expect_true(any(grepl("^md5: ", readLines(file.path(dir, "TAXDUMP_VERSION.txt")))))
})
