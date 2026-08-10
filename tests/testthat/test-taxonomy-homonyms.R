# Cross-kingdom homonyms and failed lineage lookups. taxonkit is mocked: no
# taxdump, no network. Output shapes match real `taxonkit lineage`/`reformat`.

TK_DB <- list(
  "3497" = list(name = "Morus", rank = "genus",
    lin = "Viridiplantae;Streptophyta;Magnoliopsida;Rosales;Moraceae;Morus;",
    tax = "33090;35493;3398;3744;3487;3497;"),
  "37577" = list(name = "Morus", rank = "genus",
    lin = "Metazoa;Chordata;Aves;Suliformes;Sulidae;Morus;",
    tax = "33208;7711;8782;30410;30446;37577;"),
  "7955" = list(name = "Danio rerio", rank = "species",
    lin = "Metazoa;Chordata;Actinopteri;Cypriniformes;Danionidae;Danio;Danio rerio",
    tax = "33208;7711;186623;7952;2743709;7954;7955")
)

TK_NAMES <- list(Morus = c(3497L, 37577L), "Danio rerio" = 7955L)

# args/input/dir mirror tk_run(); `drop` lists taxids taxonkit fails to report.
fake_tk_run <- function(drop = integer(0)) {
  function(args, input, dir) {
    lines <- strsplit(input, "\n")[[1]]
    lines <- lines[nzchar(lines)]
    if (identical(args[1], "name2taxid")) {
      out <- character(0)
      for (nm in lines) {
        for (tx in TK_NAMES[[nm]] %||% integer(0)) out <- c(out, paste0(nm, "\t", tx))
      }
      return(out)
    }
    ids <- vapply(strsplit(lines, "\t"), `[`, character(1), 1)
    ids <- setdiff(ids, as.character(drop))
    rec <- TK_DB[ids]
    if (identical(args[1], "reformat")) {
      return(paste(ids, "full", vapply(rec, `[[`, character(1), "lin"),
        vapply(rec, `[[`, character(1), "tax"), sep = "\t"))
    }
    if (length(args) > 1) { # lineage -n -r
      return(paste(ids, "full", vapply(rec, `[[`, character(1), "name"),
        vapply(rec, `[[`, character(1), "rank"), sep = "\t"))
    }
    paste(ids, "full", sep = "\t")
  }
}

test_that("a cross-kingdom homonym is flagged, not silently resolved", {
  testthat::local_mocked_bindings(tk_run = fake_tk_run(), .package = "RiboPilot")
  res <- RiboPilot:::resolve_taxonomy_local("Morus", "dir")
  expect_equal(res$match_type, "ambiguous")
  expect_false(is.na(res$alternatives))
  expect_true(grepl("3497", res$alternatives, fixed = TRUE))
  expect_true(grepl("37577", res$alternatives, fixed = TRUE))
  # the two kingdoms are what tells the user which Morus they were handed
  expect_true(grepl("Viridiplantae", res$alternatives, fixed = TRUE))
  expect_true(grepl("Metazoa", res$alternatives, fixed = TRUE))
})

test_that("an unambiguous name is still an exact match", {
  testthat::local_mocked_bindings(tk_run = fake_tk_run(), .package = "RiboPilot")
  res <- RiboPilot:::resolve_taxonomy_local("Danio rerio", "dir")
  expect_equal(res$match_type, "exact")
  expect_equal(res$taxid, 7955L)
  expect_true(is.na(res$alternatives))
})

test_that("tk_name2taxid keeps every candidate taxid for a name", {
  testthat::local_mocked_bindings(tk_run = fake_tk_run(), .package = "RiboPilot")
  n2t <- RiboPilot:::tk_name2taxid(c("Morus", "Danio rerio"), "dir")
  expect_equal(unname(n2t[["Morus"]]), 3497L)
  expect_equal(attr(n2t, "candidates")$Morus, c(3497L, 37577L))
  expect_equal(attr(n2t, "candidates")[["Danio rerio"]], 7955L)
})

test_that("a taxid that falls out of a batch is reported, not dropped in silence", {
  testthat::local_mocked_bindings(tk_run = fake_tk_run(drop = 37577L), .package = "RiboPilot")
  expect_warning(
    info <- RiboPilot:::tk_taxid_info(c(3497L, 37577L), "dir"),
    "37577")
  expect_equal(info$taxid, 3497L)
})

test_that("a batch that returns nothing warns and yields an empty frame", {
  testthat::local_mocked_bindings(tk_run = fake_tk_run(drop = c(3497L, 37577L)),
    .package = "RiboPilot")
  expect_warning(info <- RiboPilot:::tk_taxid_info(c(3497L, 37577L), "dir"), "3497")
  expect_equal(nrow(info), 0L)
  expect_equal(names(info), c("taxid", "sci_name", "rank", "lineage", "lineage_taxids"))
})
