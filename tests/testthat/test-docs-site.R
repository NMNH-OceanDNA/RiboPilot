test_that("site pages are ASCII and follow the app route", {
  root <- testthat::test_path("../..")
  vig <- file.path(root, "vignettes")
  if (!dir.exists(vig)) skip("package sources not available")
  pages <- c(list.files(vig, pattern = "[.]Rmd$", full.names = TRUE),
    file.path(root, "README.Rmd"))
  for (f in pages) {
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    expect_false(any(grepl("[^\\x01-\\x7F]", txt, perl = TRUE)), info = basename(f))
    expect_false(any(grepl("^\\s*RiboPilot\\(", txt)), info = basename(f))
    expect_false(any(grepl("PASTE_YOUR_KEY_HERE", txt, fixed = TRUE)), info = basename(f))
  }
  gs <- readLines(file.path(vig, "RiboPilot.Rmd"), warn = FALSE)
  expect_false(any(grepl("run_assembly\\(|run_fallback\\(|export_project\\(|conda_env", gs)))
})

test_that("the reference index names every Rd topic and only those", {
  root <- testthat::test_path("../..")
  skip_if_not(file.exists(file.path(root, "_pkgdown.yml")), "package sources not available")
  skip_if_not_installed("yaml")
  listed <- unlist(lapply(yaml::read_yaml(file.path(root, "_pkgdown.yml"))$reference, `[[`, "contents"))
  rd <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  aliases <- unlist(lapply(rd, function(f) sub("^\\\\alias\\{(.*)\\}$", "\\1",
    grep("^\\\\alias\\{", readLines(f, warn = FALSE), value = TRUE))))
  expect_equal(setdiff(listed, aliases), character(0))
  expect_equal(setdiff(sub("[.]Rd$", "", basename(rd)), listed), character(0))
})

test_that("docs and app use the same names", {
  root <- testthat::test_path("../..")
  skip_if_not(dir.exists(file.path(root, "vignettes")), "package sources not available")
  read <- function(fs) unlist(lapply(fs, readLines, warn = FALSE))
  docs <- read(c(list.files(file.path(root, "vignettes"), "[.]Rmd$", full.names = TRUE),
    file.path(root, "README.Rmd"), list.files(file.path(root, "man"), "[.]Rd$", full.names = TRUE)))
  app <- read(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE))
  app <- app[!grepl("^\\s*#", app)]
  for (b in c("GenBank screen", "all of GenBank", "anywhere in core_nt", "found remotely",
              "Taxonomy mode", "Assemble mode", "local rDNA panel", "local reference panel",
              ">local panel<", "built-in fish database", "GenBank evidence", "Top GenBank hit",
              "No curation warnings for this sample", "sample warning(s)",
              "Remote BLAST in the BLAST options")) {
    expect_false(any(grepl(b, docs, fixed = TRUE)), info = paste("docs:", b))
    expect_false(any(grepl(b, app, fixed = TRUE)), info = paste("app:", b))
  }
  tab <- "\\b(Taxonomy|Reference|Assemble|Annotate|Annotation Review|Export) tabs?\\b"
  expect_false(any(grepl(tab, docs, perl = TRUE)))
  expect_false(any(grepl(tab, app, perl = TRUE)))
})

test_that("help pages carry no raw markdown", {
  rd <- list.files(file.path(testthat::test_path("../.."), "man"), "[.]Rd$", full.names = TRUE)
  skip_if(length(rd) == 0, "package sources not available")
  for (f in rd) {
    txt <- readLines(f, warn = FALSE)
    expect_false(any(grepl("`", txt, fixed = TRUE)), info = basename(f))
    expect_false(any(grepl("\\[[A-Za-z_.][A-Za-z0-9_.]*\\(\\)\\]", txt)), info = basename(f))
  }
})

test_that("the tool-version table matches the image pins, and the citation parses", {
  root <- testthat::test_path("../..")
  df <- file.path(root, "docker", "Dockerfile")
  vig <- file.path(root, "vignettes", "Citing.Rmd")
  skip_if_not(file.exists(df) && file.exists(vig), "package sources not available")
  pins <- paste(readLines(df, warn = FALSE), collapse = "\n")
  rows <- grep("^\\| .* \\| `[A-Za-z0-9.-]+` \\| [0-9]", readLines(vig, warn = FALSE), value = TRUE)
  expect_gte(length(rows), 10)
  for (r in rows) {
    cell <- trimws(strsplit(r, "|", fixed = TRUE)[[1]])
    pkg <- gsub("`", "", cell[3])
    expect_true(grepl(paste0("'", pkg, "=", cell[4]), pins, fixed = TRUE) ||
      grepl(paste0(pkg, "-", cell[4]), pins, fixed = TRUE), info = r)
  }
  cit <- utils::readCitationFile(file.path(root, "inst", "CITATION"),
    meta = list(Version = "9.9.9", Encoding = "UTF-8"))
  expect_match(format(cit, style = "text"), "9.9.9", fixed = TRUE)
})
