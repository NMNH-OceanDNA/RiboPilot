test_that("export_metadata_cols keeps only the user's mapping-file columns", {
  sample_cols <- c("ID", "Taxon", "R1", "R2", "Assembly", "taxid", "lineage",
                   "taxonomy_lock", "Donors", "Expected")
  declared <- c(".selection", "ID", "path", "Taxon", "curate_opts", "export_group")
  expect_equal(export_metadata_cols(sample_cols, declared), c("Donors", "Expected"))
})

test_that("export_metadata_cols is empty when the mapping file had no extras", {
  expect_equal(
    export_metadata_cols(c("ID", "Taxon", "R1", "R2"), c("ID", "Taxon")),
    character(0)
  )
})

test_that("export_metadata_cols treats the specimen status columns as owned", {
  expect_equal(export_metadata_cols(c("ID", "specimen", "specimen_message", "specimen_icons", "site"),
                                    character()), "site")
})
