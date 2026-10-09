# Export reads per-contig data only through the review views (S9). A static
# check over the SQL strings of the export code: the original per-contig tables
# and the check columns of annotate are named only where the plan allows it.

export_view_root <- function() {
  root <- testthat::test_path("../..")
  if (!file.exists(file.path(root, "R", "export.R"))) skip("package sources not available")
  root
}

EXPORT_VIEW_FILES <- c("R/export.R", "R/export_opts.R", "R/export_report.R")
EXPORT_ORIG_TABLES <- c("assemblies", "annotations", "screen", "chimera", "chimera_windows",
  "coverage", "blast_hits", "blast_screen", "blast_ref_alignment")
EXPORT_CHECK_COLS <- c("n_features", "missing", "chimeric", "n_flags", "curate_flags",
  "itsx_warnings", "itsx_origin", "itsx_scores")
# Top-level object -> original tables its SQL may name. unit_fasta_text is Copy
# FASTA on the Assemble step; the others read assemblies only for ignore; the
# review = FALSE branches of unit_checks and export_summary_df serve the
# Assemble step.
EXPORT_VIEW_ALLOW <- list(
  unit_fasta_text = "assemblies", UNIT_MULTI_SQL = "assemblies",
  exportable_units = "assemblies", export_left_out = "assemblies",
  unit_checks = c(EXPORT_ORIG_TABLES, "annotate"),
  export_summary_df = c(EXPORT_ORIG_TABLES, "annotate"))

# Every string constant of a file, with the top-level object it sits in
export_strings <- function(f) {
  exprs <- parse(f, keep.source = TRUE)
  pd <- utils::getParseData(exprs)
  out <- lapply(seq_along(exprs), function(i) {
    e <- exprs[[i]]
    nm <- if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=")) as.character(e[[2]]) else ""
    sr <- attr(exprs, "srcref")[[i]]
    s <- pd[pd$token == "STR_CONST" & pd$line1 >= sr[1] & pd$line2 <= sr[3], ]
    if (nrow(s) == 0) return(NULL)
    data.frame(fn = nm, file = basename(f),
      txt = vapply(s$text, function(x) eval(parse(text = x)), character(1)),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

sql_tables <- function(txt) {
  m <- regmatches(txt, gregexpr("\\b(FROM|JOIN)\\s+[A-Za-z_]+", txt, perl = TRUE))[[1]]
  unique(sub("^(FROM|JOIN)\\s+", "", m))
}

test_that("export code names the original per-contig tables only where allowed", {
  root <- export_view_root()
  s <- do.call(rbind, lapply(file.path(root, EXPORT_VIEW_FILES), export_strings))
  expect_gt(nrow(s), 50)
  bad <- character(0)
  for (i in seq_len(nrow(s))) {
    tb <- sql_tables(s$txt[i])
    allow <- EXPORT_VIEW_ALLOW[[s$fn[i]]] %||% character(0)
    orig <- setdiff(intersect(tb, EXPORT_ORIG_TABLES), allow)
    if (length(orig)) bad <- c(bad, sprintf("%s %s reads %s", s$file[i], s$fn[i], paste(orig, collapse = ",")))
    if (s$fn[i] %in% c("unit_checks", "export_summary_df")) next
    cols <- paste(EXPORT_CHECK_COLS, collapse = "|")
    if (grepl(sprintf("\\ban\\.(%s)\\b", cols), s$txt[i], perl = TRUE) ||
        ("annotate" %in% tb && grepl(sprintf("(?<![.\\w])(%s)\\b", cols), s$txt[i], perl = TRUE))) {
      bad <- c(bad, sprintf("%s %s reads a check column of annotate", s$file[i], s$fn[i]))
    }
  }
  expect_equal(bad, character(0))
  # the functions allowed to read originals also read the views
  for (fn in c("unit_checks", "export_summary_df")) {
    txt <- paste(s$txt[s$fn == fn], collapse = "\n")
    for (v in c("review_units_v", "review_annotations_v", "review_screen_v", "review_chimera_v")) {
      expect_match(txt, v, fixed = TRUE, info = paste(fn, v))
    }
  }
})

test_that("every export summary in the export code reads the review layer", {
  root <- export_view_root()
  n <- 0
  for (f in file.path(root, EXPORT_VIEW_FILES)) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    fc <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "export_summary_df"]
    for (id in pd$parent[match(fc, pd$id)]) {
      expect_match(utils::getParseText(pd, id), "review = TRUE", fixed = TRUE, info = basename(f))
      n <- n + 1
    }
  }
  expect_gte(n, 2)
  rmd <- readLines(file.path(root, "inst", "export_report.Rmd"), warn = FALSE)
  calls <- grep("export_summary_df(", rmd, fixed = TRUE, value = TRUE)
  expect_length(calls, 1)
  expect_match(calls, "review = TRUE", fixed = TRUE)
  expect_false(any(grepl("dbGetQuery|dbExecute|dbReadTable|dbSendQuery", rmd)))
})
