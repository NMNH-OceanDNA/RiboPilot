mem_samples <- function(samples, envir = parent.frame()) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  withr::defer(DBI::dbDisconnect(con), envir = envir)
  DBI::dbWriteTable(con, "samples", samples)
  .meta_ensure_tables(con)
  con
}

test_that(".meta_save_fields replaces the selection", {
  con <- mem_samples(data.frame(ID = "s1", Taxon = "x"))
  .meta_save_fields(con, c("geome:combo:lat_lon", "geome:raw:Event:country"))
  .meta_save_fields(con, "geome:combo:lat_lon")
  expect_equal(DBI::dbGetQuery(con, "SELECT key FROM meta_export_fields")$key, "geome:combo:lat_lon")
  .meta_save_fields(con, character())
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM meta_export_fields")$n, 0L)
})

test_that("specimen_fields_modal has a section per source", {
  con <- mem_samples(data.frame(ID = "s1", Taxon = "x"))
  DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "s1", source = c("GEOME", "GBIF"), level = c("Event", "Occurrence"), depth = 0L,
    ref = "r", field = c("country", "countryCode"), value = c("Peru", "PE")))
  .meta_save_fields(con, c("geome:combo:lat_lon", "gbif:combo:sex"))
  sm <- lapply(stats::setNames(nm = names(META_SOURCES)), function(s) meta_field_summary(con, s))
  html <- as.character(specimen_fields_modal(NS("exp"), sm))
  expect_match(html, "Metadata fields for export", fixed = TRUE)
  for (s in c("GEOME", "GBIF", "NCBI")) expect_match(html, paste0("<h4>", s, "</h4>"), fixed = TRUE)
  for (k in c("geome", "gbif", "ncbi")) {
    expect_match(html, paste0("exp-", k, "_combos"), fixed = TRUE)
    expect_match(html, paste0("exp-", k, "_raw"), fixed = TRUE)
  }
  expect_match(html, "{gbif_specimen_voucher}", fixed = TRUE)
  expect_match(html, "value=\"gbif:combo:sex\" checked", fixed = TRUE)
  expect_lt(regexpr("exp-geome_combos", html), regexpr("exp-gbif_combos", html))
})

status_db <- function(envir = parent.frame()) {
  con <- mem_samples(data.frame(
    ID = c("a", "b", "c", "d"), Taxon = "x",
    GEOME_BCID = c("ark:/1/A", NA, "ark:/1/C", NA), GBIF_ID = c(NA, "2", "3", NA)), envir)
  DBI::dbAppendTable(con, "meta_status", data.frame(
    ID = c("a", "b", "c", "c"), source = c("GEOME", "GBIF", "GEOME", "GBIF"),
    ref = c("ark:/1/A", "2", "ark:/1/C", "3"), status = c("ok", "failed", "ok", "ok"),
    message = c(NA, "GBIF returned HTTP 503", NA, NA), fetched_at = 1L))
  DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "c", source = c("GEOME", "GBIF", "GEOME", "GBIF"),
    level = c("Event", "Occurrence", "Event", "Occurrence"), depth = 0L, ref = "r",
    field = c("country", "countryCode", "collectorList", "recordedBy"),
    value = c("Peru", "CL", "A. B", "Ann B")))
  con
}

test_that("specimen_status: failed beats ok beats none", {
  con <- status_db()
  s <- specimen_status(con)
  s <- s[order(s$ID), ]
  expect_equal(s$specimen[c(1, 2, 4)], c("ok", "failed", "none"))
  expect_equal(s$specimen_message[1], "GEOME: fetched")
  expect_equal(s$specimen_message[2], "GBIF: failed (GBIF returned HTTP 503)")
  expect_equal(s$specimen_message[4], "No GEOME, GBIF, or NCBI ID")
  expect_equal(s$specimen_icons[c(1, 2, 4)], c("GEOME:ok", "GBIF:failed", ""))
})

test_that("a conflict outranks ok once the comparison rules are present", {
  con <- status_db()
  s <- specimen_status(con)
  s <- s[order(s$ID), ]
  expect_equal(s$specimen[3], "conflict")
  expect_equal(s$specimen_message[3],
               "GEOME: fetched\nGBIF: fetched\nConflicts: country\nNot checked: collector")
  expect_equal(s$specimen_icons[3], "GEOME:ok GBIF:ok conflict")
})

test_that("a set but never fetched ID says so", {
  con <- mem_samples(data.frame(ID = "a", Taxon = "x", GBIF_ID = "5"))
  s <- specimen_status(con)
  expect_equal(s$specimen, "none")
  expect_equal(s$specimen_message, "GBIF: not fetched yet")
  expect_equal(s$specimen_icons, "GBIF:pending")
})

test_that(".specimen_status_join adds specimen columns to a collected table", {
  con <- status_db()
  out <- .specimen_status_join(data.frame(ID = c("b", "a", "zz")), con)
  expect_equal(out$specimen, c("failed", "ok", "none"))
  expect_true(is.na(out$specimen_message[3]))
})

test_that("rt_specimen draws one logo per source, a conflict flag, or a plus", {
  js <- as.character(rt_specimen("specimen_open"))
  expect_match(js, "specimen_open", fixed = TRUE)
  expect_match(js, "setInputValue", fixed = TRUE)
  expect_match(js, "dataset.id", fixed = TRUE)
  expect_match(js, "specimen_message", fixed = TRUE)
  expect_match(js, "specimen_icons", fixed = TRUE)
  for (x in c("www/specimen/geome_g.png", "www/specimen/gbif_leaf.png", "www/specimen/ncbi_helix.png",
              "fa-flag", "fa-square-plus", "fa-triangle-exclamation", "mp-spec-faded")) {
    expect_match(js, x, fixed = TRUE)
  }
})

test_that("specimen logo files and the countries list ship with the app", {
  for (f in c("geome_g.png", "gbif_leaf.png", "ncbi_helix.png")) {
    expect_true(file.exists(app_sys("app", "www", "specimen", f)))
  }
  expect_true(file.exists(app_sys("extdata", "countries.csv")))
})

test_that("specimen_col_def is a Metadata column", {
  cd <- specimen_col_def("specimen_open")
  expect_null(cd$class)
  expect_equal(cd$name, "Metadata")
  expect_equal(specimen_col_def("specimen_open", sticky = "left")$sticky, "left")
})

test_that("meta_record_view orders GEOME root first and links BCIDs", {
  recs <- data.frame(
    level = c("Tissue", "Event", "Project"), depth = c(0L, 2L, 3L),
    ref = c("ark:/1/T", "ark:/1/E", NA), field = c("tissueID", "country", "projectTitle"),
    value = c("T1", "Peru", "My proj"))
  html <- as.character(meta_record_view(recs, "GEOME", "m-geome_records"))
  expect_lt(regexpr("Project", html), regexpr("Event", html))
  expect_lt(regexpr("Event", html), regexpr("Tissue", html))
  expect_match(html, "https://geome-db.org/record/ark:/1/E", fixed = TRUE)
  expect_match(html, "id=\"m-geome_records\"[^>]*overflow-y: auto")
  expect_match(html, "#m-geome_records details", fixed = TRUE)
  expect_match(html, "Expand all", fixed = TRUE)
  expect_false(grepl("Expand all", as.character(meta_record_view(recs[0, ], "GEOME", "x"))))
})

test_that("meta_record_view shows GBIF occurrence first, issues as badges, and the citation", {
  recs <- data.frame(
    level = c("Occurrence", "Occurrence", "Dataset", "Dataset", "Organization"),
    depth = c(0L, 0L, 1L, 1L, 2L), ref = c("61", "61", "ds1", "ds1", "org1"),
    field = c("issues", "country", "title", "citation", "title"),
    value = c("ZERO_COORDINATE,COORDINATE_ROUNDED", "Peru", "UF Fish",
              "Robins R (2026). UF Fish.", "Florida Museum"))
  html <- as.character(meta_record_view(recs, "GBIF", "b2"))
  expect_lt(regexpr(">Occurrence<", html), regexpr(">Dataset<", html))
  expect_lt(regexpr(">Dataset<", html), regexpr(">Organization<", html))
  expect_match(html, "https://www.gbif.org/occurrence/61", fixed = TRUE)
  expect_match(html, "https://www.gbif.org/dataset/ds1", fixed = TRUE)
  expect_match(html, "https://www.gbif.org/publisher/org1", fixed = TRUE)
  expect_match(html, "mp-pill mp-pill-warning\">ZERO_COORDINATE<", fixed = TRUE)
  expect_match(html, "mp-pill mp-pill-warning\">COORDINATE_ROUNDED<", fixed = TRUE)
  expect_match(html, "class=\"mp-meta-citation\"", fixed = TRUE)
})

test_that("meta_record_view links NCBI levels in SRA, BioSample, BioProject order", {
  recs <- data.frame(level = c("BioProject", "SRA", "BioSample"), depth = c(2L, 0L, 1L),
                     ref = c("PRJNA1", "SRR1", "SAMN1"), field = "accession",
                     value = c("PRJNA1", "SRR1", "SAMN1"))
  html <- as.character(meta_record_view(recs, "NCBI", "m-ncbi_records"))
  expect_lt(regexpr("sra/SRR1", html), regexpr("biosample/SAMN1", html))
  expect_lt(regexpr("biosample/SAMN1", html), regexpr("bioproject/PRJNA1", html))
})

test_that("specimen_compare_view marks conflicts and notes and names the CSV column", {
  cf <- data.frame(ID = "s1", concept = c("country", "collector", "sex", "voucher"),
                   csv_column = c("geo_loc_name", NA, NA, NA),
                   csv_value = c("USA: Florida", NA, NA, NA),
                   geome_value = c("Canada", "A. B", NA, NA),
                   gbif_value = c(NA, "Ann B", "male", NA),
                   ncbi_value = c(NA, NA, NA, "Canada"),
                   status = c("conflict", "not checked", "single", NA))
  html <- as.character(specimen_compare_view(cf))
  expect_match(html, "<tr class=\"mp-spec-conflict\">", fixed = TRUE)
  expect_match(html, "<tr class=\"text-muted\">", fixed = TRUE)
  expect_match(html, "(geo_loc_name)", fixed = TRUE)
  expect_match(html, "<th>NCBI</th>", fixed = TRUE)
  expect_match(html, "Mapping file (column)", fixed = TRUE)
})

test_that("specimen_csv_map_ui offers every concept but taxon, with auto and none", {
  current <- list(coordinates = c("Latitude", "Longitude"), collection_date = character(),
                  country = "Where", locality = character(), voucher = character(),
                  collector = character(), sex = character(), dev_stage = character(),
                  taxon = "Taxon")
  html <- as.character(specimen_csv_map_ui(NS("v"), current, c(country = "Where", sex = ""),
                                           c("Latitude", "Longitude", "Where")))
  expect_match(html, "Mapping file columns...", fixed = TRUE)
  expect_match(html, "v-map_coordinates", fixed = TRUE)
  expect_false(grepl("v-map_taxon", html, fixed = TRUE))
  expect_match(html, "auto: Latitude + Longitude", fixed = TRUE)
  expect_match(html, "__none__", fixed = TRUE)
  expect_match(html, "\"maxItems\":2", fixed = TRUE)
  expect_match(html, "v-map_save", fixed = TRUE)
})

viewer_proj <- function(envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  path <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  DBI::dbWriteTable(con, "samples", data.frame(ID = c("s1", "s2"), Taxon = "x", country = "USA"))
  .meta_ensure_tables(con)
  DBI::dbExecute(con, "UPDATE samples SET BioSample = 'SAMN1', GBIF_ID = '77' WHERE ID = 's1'")
  DBI::dbExecute(con, "INSERT INTO meta_links VALUES ('s1', 'GBIF', '77', 'NCBI BioSample voucherURI', NULL)")
  DBI::dbExecute(con, "INSERT INTO meta_links VALUES ('s1', 'GEOME', NULL, NULL, '2 possible GBIF matches; not linked')")
  for (src in c("NCBI", "GBIF")) {
    DBI::dbAppendTable(con, "meta_records", data.frame(ID = "s1", source = src, level = "L", depth = 0L,
                                                       ref = "r", field = "f", value = "v"))
    DBI::dbAppendTable(con, "meta_status", data.frame(ID = "s1", source = src, ref = "r", status = "ok",
                                                      message = NA_character_, fetched_at = 1L))
  }
  DBI::dbDisconnect(con)
  path
}
open_path <- function(path) function() DBI::dbConnect(RSQLite::SQLite(), path)
q_path <- function(path, sql) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, sql)
}

test_that("the viewer NCBI tab has an ID box", {
  path <- viewer_proj()
  open <- shiny::reactiveVal(NULL)
  shiny::testServer(specimen_viewer_server, args = list(open = open, con = open_path(path)), {
    open("s2")
    session$flushReact()
    html <- output$ncbi_detail$html
    expect_match(html, "SAMN29555051 or SRR21844202", fixed = TRUE)
    expect_match(html, "This sample has no BioSample", fixed = TRUE)
  })
})

test_that("the viewer shows where a linked ID came from and link notes", {
  path <- viewer_proj()
  open <- shiny::reactiveVal(NULL)
  shiny::testServer(specimen_viewer_server, args = list(open = open, con = open_path(path)), {
    open("s1")
    session$flushReact()
    expect_match(output$gbif_detail$html, "Found through NCBI BioSample voucherURI", fixed = TRUE)
    expect_match(output$geome_detail$html, "2 possible GBIF matches; not linked", fixed = TRUE)
  })
  con <- open_path(path)()
  on.exit(DBI::dbDisconnect(con))
  st <- specimen_status(con)
  expect_match(st$specimen_message[st$ID == "s1"], "GBIF: fetched (linked from NCBI BioSample voucherURI)", fixed = TRUE)
})

test_that("the viewer's remove panel deletes only the chosen fetched data and calls on_change", {
  path <- viewer_proj()
  open <- shiny::reactiveVal(NULL)
  changed <- 0L
  shiny::testServer(specimen_viewer_server,
    args = list(open = open, con = open_path(path), on_change = function() changed <<- changed + 1L), {
    open("s1")
    session$flushReact()
    session$setInputs(remove_sources = "GBIF", remove_scope = "sample")
    session$setInputs(remove_go = 1)
  })
  expect_equal(changed, 1L)
  expect_equal(q_path(path, "SELECT source FROM meta_records")$source, "NCBI")
  s <- q_path(path, "SELECT * FROM samples WHERE ID = 's1'")
  expect_true(is.na(s$GBIF_ID))
  expect_equal(s$BioSample, "SAMN1")
  expect_equal(s$country, "USA")
})

test_that("the link checkbox saves the project switch", {
  path <- viewer_proj()
  open <- shiny::reactiveVal(NULL)
  shiny::testServer(specimen_viewer_server, args = list(open = open, con = open_path(path)), {
    session$setInputs(gbif_link = TRUE)
  })
  con <- open_path(path)()
  on.exit(DBI::dbDisconnect(con))
  expect_true(.meta_link_enabled(con))
})

test_that("the app adds a Metadata column after Taxon in Taxonomy, Assemble, Annotate, and Export", {
  skip_if_not_installed("withr")
  withr::local_dir(mk_app_project())
  shiny::testServer(RiboPilot:::app_server, {
    for (m in c("Taxonomy", "Assemble", "Annotate")) {
      session$setInputs(mode = m)
      d <- all_rows()
      i <- match("Taxon", names(d))
      expect_equal(names(d)[i + 1:3], c("specimen", "specimen_message", "specimen_icons"), info = m)
      expect_equal(d$specimen[1], "none")
    }
  })
})

test_that("Set Export Metadata saves ticked fields and the Export window offers them as chips", {
  skip_if_not_installed("withr")
  withr::local_dir(mk_app_project())
  con <- RiboPilot:::app_con()
  DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "S1", source = "GEOME", level = "Event", depth = 0L, ref = NA_character_,
    field = "country", value = "Peru"))
  DBI::dbDisconnect(con)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(specimen_fields = 1)
    expect_null(export_snap())
    session$setInputs(geome_combos = "geome:combo:lat_lon")
    session$setInputs(specimen_fields_save = 1)
    ticked <- RiboPilot:::app_query("SELECT key FROM meta_export_fields")$key
    expect_true("geome:combo:lat_lon" %in% ticked)
    expect_equal(meta_fields_ver(), 1L)
  })
})

test_that("opening a sample from the Metadata column runs the viewer without error", {
  skip_if_not_installed("withr")
  withr::local_dir(mk_app_project())
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Taxonomy")
    expect_no_error({
      session$setInputs(specimen_open = "S1")
      session$flushReact()
    })
  })
})

fake_fetch <- function(status, calls) {
  function(con, source, ids, refs, cache = new.env(), link = FALSE) {
    calls$n <- calls$n + length(ids)
    calls$ids <- c(calls$ids, ids)
    data.frame(ID = ids, status = status, message = if (status == "failed") "boom" else NA_character_)
  }
}

test_that("Fetch stores the typed ID, calls the fetcher, and fires on_change once", {
  path <- viewer_proj()
  calls <- new.env(); calls$n <- 0L
  local_mocked_bindings(.meta_fetch_into = fake_fetch("ok", calls), .package = "RiboPilot")
  open <- shiny::reactiveVal(NULL)
  changed <- 0L
  shiny::testServer(specimen_viewer_server,
    args = list(open = open, con = open_path(path), on_change = function() changed <<- changed + 1L), {
    open("s2")
    session$flushReact()
    session$setInputs(ncbi_ref = "SAMN29555051")
    session$setInputs(ncbi_fetch = 1)
  })
  expect_equal(q_path(path, "SELECT BioSample FROM samples WHERE ID = 's2'")$BioSample, "SAMN29555051")
  expect_equal(calls$n, 1L)
  expect_equal(changed, 1L)
})

test_that("a failed fetch does not throw", {
  path <- viewer_proj()
  calls <- new.env(); calls$n <- 0L
  local_mocked_bindings(.meta_fetch_into = fake_fetch("failed", calls), .package = "RiboPilot")
  open <- shiny::reactiveVal(NULL)
  changed <- 0L
  shiny::testServer(specimen_viewer_server,
    args = list(open = open, con = open_path(path), on_change = function() changed <<- changed + 1L), {
    open("s2")
    session$flushReact()
    session$setInputs(ncbi_ref = "SAMN29555051")
    expect_no_error(session$setInputs(ncbi_fetch = 1))
  })
  expect_equal(calls$n, 1L)
  expect_equal(changed, 1L)
})

test_that("Refresh all with no IDs only notifies, and with IDs fetches each job", {
  path <- viewer_proj()
  con <- open_path(path)()
  DBI::dbExecute(con, "UPDATE samples SET BioSample = NULL, GBIF_ID = NULL")
  DBI::dbDisconnect(con)
  calls <- new.env(); calls$n <- 0L
  local_mocked_bindings(.meta_fetch_into = fake_fetch("ok", calls), .package = "RiboPilot")
  open <- shiny::reactiveVal(NULL)
  changed <- 0L
  shiny::testServer(specimen_viewer_server,
    args = list(open = open, con = open_path(path), on_change = function() changed <<- changed + 1L), {
    expect_no_error(session$setInputs(refresh_all = 1))
  })
  expect_equal(calls$n, 0L)
  expect_equal(changed, 0L)

  con <- open_path(path)()
  DBI::dbExecute(con, "UPDATE samples SET BioSample = 'SAMN1', GBIF_ID = '77' WHERE ID = 's1'")
  DBI::dbDisconnect(con)
  shiny::testServer(specimen_viewer_server,
    args = list(open = open, con = open_path(path), on_change = function() changed <<- changed + 1L), {
    session$setInputs(refresh_all = 1)
  })
  expect_equal(calls$n, 2L)
  expect_equal(changed, 1L)
})

test_that("metadata joins accept a zero-row table whose ID column is logical", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  d0 <- data.frame(ID = logical(0), Taxon = character(0))
  expect_no_error(out <- .specimen_status_join(d0, con))
  expect_equal(nrow(out), 0)
  expect_no_error(meta_view_join(d0, con))
  DBI::dbExecute(con, "INSERT OR IGNORE INTO meta_export_fields VALUES ('geome:combo:lat_lon')")
  expect_no_error(.meta_join(d0, con))
})
