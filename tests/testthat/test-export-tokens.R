tok_data <- function() {
  data.frame(
    ID = c("s1", "s2"), seqid = c("s1", "s2"), path = 1L, scaffold = 1L,
    Taxon = c("Fish a", "Fish b"), organism = c("Fish a", "Fish b"), sci_name = "Fish",
    topology = "linear", completeness = "complete sequence", site = c("Pond", "Lake"),
    GEOME_BCID = c("ark:/1/A", NA), GBIF_ID = c(NA, "6186461308"),
    geome_lat_lon = c("17.5 S 149.8 W", ""), gbif_Occurrence_countryCode = c("", "US"),
    length = 1500L, reviewed = "yes", element = "18S", product = "18S ribosomal RNA",
    export_group = "g1", specimen = "ok", specimen_message = "GEOME: fetched", specimen_icons = "GEOME:ok"
  )
}
tok_cols <- c("ID", "Taxon", "R1", "R2", "site", "GEOME_BCID", "GBIF_ID", "taxid", "lineage")

test_that("export_token_groups sorts columns into the six groups", {
  g <- export_token_groups(tok_data(), tok_cols,
                           c("geome:combo:lat_lon", "gbif:raw:Occurrence:countryCode"))
  expect_equal(names(g), c("Basics", "Your mapping columns", "GEOME", "GBIF", "NCBI", "RiboPilot columns"))
  expect_equal(g$Basics$tokens$token,
               c("seqid", "ID", "Taxon", "organism", "sci_name", "topology", "completeness", "path", "scaffold"))
  expect_equal(g$Basics$tokens$insert[3], "{Taxon}")
  expect_equal(g$`Your mapping columns`$tokens$token, "site")
  expect_equal(g$GEOME$tokens$token, c("GEOME_BCID", "geome_lat_lon"))
  expect_equal(g$GEOME$tokens$insert, c("{GEOME_BCID}", "[lat_lon={geome_lat_lon}]"))
  expect_equal(g$GEOME$tokens$example, c("ark:/1/A", "17.5 S 149.8 W"))
  expect_equal(g$GBIF$tokens$insert, c("{GBIF_ID}", "{gbif_Occurrence_countryCode}"))
  expect_equal(g$GBIF$tokens$example, c("", ""))
  expect_equal(g$`RiboPilot columns`$tokens$token, c("length", "reviewed", "element", "product"))
  expect_equal(unname(vapply(g, function(x) x$open, logical(1))),
               c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  listed <- unlist(lapply(g, function(x) x$tokens$token))
  expect_false(any(c("export_group", "specimen", "specimen_message", "specimen_icons", "R1", "R2",
                     "taxid", "lineage") %in% listed))
  expect_null(g$Basics$hint)
})

test_that("groups with nothing ticked or no extra columns say how to add some", {
  d <- tok_data()[, c("ID", "seqid", "Taxon", "GEOME_BCID", "GBIF_ID")]
  g <- export_token_groups(d, c("ID", "Taxon", "GEOME_BCID", "GBIF_ID"), character())
  expect_false(g$GEOME$open)
  expect_false(g$GBIF$open)
  expect_match(g$GEOME$hint, "Metadata button", fixed = TRUE)
  expect_match(g$GBIF$hint, "Metadata button", fixed = TRUE)
  expect_true(g$`Your mapping columns`$open)
  expect_null(g$`Your mapping columns`$hint)
  expect_equal(nrow(g$`Your mapping columns`$tokens), 2L)
  expect_equal(g$GEOME$tokens$token, character(0))
})

test_that("export_token_groups copes with no rows", {
  g <- export_token_groups(tok_data()[0, ], tok_cols, "geome:combo:lat_lon")
  expect_equal(g$GEOME$tokens$example, c("", ""))
  g0 <- export_token_groups(data.frame(), tok_cols, character())
  expect_equal(nrow(g0$Basics$tokens), 0L)
})

test_that("export_token_ui renders chips with insert text, a filter box, and picker links", {
  g <- export_token_groups(tok_data(), tok_cols, "geome:combo:lat_lon")
  html <- as.character(export_token_ui(g, target_id = "fasta_header", ns = NS(NULL)))
  expect_match(html, "data-target=\"fasta_header\"", fixed = TRUE)
  expect_match(html, "data-insert=\"[lat_lon={geome_lat_lon}]\"", fixed = TRUE)
  expect_match(html, "data-insert=\"{Taxon}\"", fixed = TRUE)
  expect_match(html, "mp-token-filter", fixed = TRUE)
  expect_match(html, "id=\"token_fields_geome\"", fixed = TRUE)
  expect_match(html, "id=\"token_fields_gbif\"", fixed = TRUE)
  expect_match(html, "First row: Pond", fixed = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr("<details open", html))), 3L)
})

test_that("token chips carry missing-value counts per export group", {
  d <- data.frame(ID = c("a", "b", "c"), Taxon = "x", export_group = c("g1", "g1", "g2"),
                  site = c("reef", NA, ""))
  g <- export_token_groups(d, c("ID", "Taxon", "site"), character())
  site <- g$`Your mapping columns`$tokens
  expect_equal(jsonlite::fromJSON(site$missing[site$token == "site"]), list(g1 = 1L, g2 = 1L))
  ids <- g$Basics$tokens
  expect_equal(jsonlite::fromJSON(ids$missing[ids$token == "ID"]), list(g1 = 0L, g2 = 0L))
  html <- as.character(export_token_ui(g, target_id = "x", ns = NS(NULL), totals = '{"g1":2,"g2":1}'))
  expect_match(html, "data-totals=\"{&quot;g1&quot;:2,&quot;g2&quot;:1}\"", fixed = TRUE)
  expect_match(html, "Columns by type", fixed = TRUE)
})

test_that("NCBI biosample and bioproject combos insert BioSample/BioProject modifiers", {
  expect_equal(meta_export_token(c("ncbi:combo:biosample", "ncbi:combo:bioproject",
                                   "ncbi:combo:biosample_sra", "ncbi:combo:bioproject_sra",
                                   "ncbi:combo:lat_lon")),
               c("[BioSample={ncbi_biosample}]", "[BioProject={ncbi_bioproject}]",
                 "[BioSample={ncbi_biosample_sra}]", "[BioProject={ncbi_bioproject_sra}]",
                 "[lat_lon={ncbi_lat_lon}]"))
})
