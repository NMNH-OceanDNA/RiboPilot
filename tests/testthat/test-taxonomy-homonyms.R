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
})

test_that("an unambiguous name is still an exact match", {
  res <- RiboPilot:::resolve_taxonomy("Danio rerio", local_taxdump())
  expect_equal(res$match_type, "exact")
  expect_equal(res$taxid, 7955L)
  expect_true(is.na(res$alternatives))
  expect_equal(res$lineage, "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio")
  expect_equal(res$lineage_taxids, "33208;7711;186623;7952;2743709;7954;7955")
})

test_that("tk_name2taxid keeps every candidate taxid for a name", {
  dir <- local_taxdump()
  n2t <- RiboPilot:::with_tax_index(dir, function(con)
    RiboPilot:::tk_name2taxid(c("Morus", "Danio rerio"), con))
  expect_equal(unname(n2t[["Morus"]]), 3497L)
  expect_equal(attr(n2t, "candidates")$Morus, c(3497L, 37577L))
  expect_equal(attr(n2t, "candidates")[["Danio rerio"]], 7955L)
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
