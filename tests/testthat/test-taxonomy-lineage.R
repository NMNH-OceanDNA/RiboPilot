# Lineage strings are read positionally by parse_rank_names/parse_rank_taxids, so
# every producer must emit the same seven rank-aligned slots. A producer that
# closes up an absent rank shifts every deeper rank one slot shallower, which is
# silent: the string still looks like a lineage.

test_that("rank_aligned_lineage leaves absent ranks empty rather than closing up", {
  # Chondromyces apiculatus has no assigned class. Collapsing that slot would put
  # the order in the class slot, the family in the order slot, and so on.
  lin <- data.frame(
    rank = c("kingdom", "phylum", "order", "family", "genus", "species"),
    id = c(2L, 29, 1, 49, 3, 51),
    name = c("Pseudomonadati", "Myxococcota", "Polyangiales", "Polyangiaceae",
             "Chondromyces", "Chondromyces apiculatus"),
    stringsAsFactors = FALSE)
  out <- RiboPilot:::rank_aligned_lineage(lin)
  expect_equal(out$names,
    "Pseudomonadati;Myxococcota;;Polyangiales;Polyangiaceae;Chondromyces;Chondromyces apiculatus")
  expect_equal(out$taxids, "2;29;;1;49;3;51")

  nm <- RiboPilot:::parse_rank_names(out$names)
  expect_true(is.na(nm[["class"]]))
  expect_equal(unname(nm[["order"]]), "Polyangiales")
  expect_equal(unname(RiboPilot:::parse_rank_taxids(out$taxids)[["species"]]), 51L)
})

test_that("rank_aligned_lineage drops domain and always emits seven fields", {
  lin <- data.frame(rank = c("domain", "kingdom", "species"),
    id = c(2L, 33208L, 7955L),
    name = c("Bacteria", "Metazoa", "Danio rerio"), stringsAsFactors = FALSE)
  out <- RiboPilot:::rank_aligned_lineage(lin)
  # domain is not one of the seven ranks, so it must not occupy the kingdom slot
  expect_equal(out$names, "Metazoa;;;;;;Danio rerio")
  expect_length(strsplit(out$names, ";")[[1]], 7)
  expect_equal(unname(RiboPilot:::parse_rank_names(out$names)[["kingdom"]]), "Metazoa")
})

test_that("rank_aligned_lineage handles an empty classification", {
  out <- RiboPilot:::rank_aligned_lineage(data.frame())
  expect_true(is.na(out$names))
  expect_true(is.na(out$taxids))
})
