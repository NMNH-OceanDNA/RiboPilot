bc_project <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  RiboPilot::new_db(db_path = file.path(d, ".sqlite"), mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  d
}

bc_history <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "SELECT * FROM schema_history ORDER BY rowid")
}

bc_schema <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  tbls <- sort(setdiff(DBI::dbListTables(con), "sqlite_sequence"))
  stats::setNames(lapply(tbls, function(t) sort(DBI::dbListFields(con, t))), tbls)
}

bc_views <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  v <- DBI::dbGetQuery(con, "SELECT name, sql FROM sqlite_master WHERE type = 'view' ORDER BY name")
  stats::setNames(gsub("\\s+", " ", v$sql), v$name)
}

fake_step <- function(note = "add samples.extra") {
  list(note = note, fn = function(con) DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN extra TEXT"))
}

test_that("new_db records a created row at the current schema version", {
  h <- bc_history(bc_project())
  expect_equal(nrow(h), 1)
  expect_equal(h$event, "created")
  expect_equal(h$schema_version, RiboPilot:::current_schema_version())
  expect_equal(h$ribopilot_version, as.character(utils::packageVersion("RiboPilot")))
})

test_that("an up-to-date project writes nothing", {
  dir <- bc_project()
  expect_message(backwards_compatibility(dir), "up to date")
  expect_equal(nrow(bc_history(dir)), 1)
  expect_length(list.files(dir, pattern = "\\.bak\\."), 0)
})

test_that("a pending step upgrades, logs history, and backs up first", {
  dir <- bc_project()
  v <- RiboPilot:::current_schema_version()
  steps <- RiboPilot:::db_migrations()
  local_mocked_bindings(db_migrations = function() c(steps, list(fake_step())))
  suppressMessages(backwards_compatibility(dir))
  h <- bc_history(dir)
  expect_equal(h$event, c("created", "upgraded"))
  expect_equal(h$schema_version, c(v, v + 1L))
  expect_equal(h$note[2], "add samples.extra")
  expect_equal(RiboPilot:::project_schema_version(file.path(dir, ".sqlite")), v + 1L)
  bak <- list.files(dir, pattern = "^\\.sqlite\\.bak\\.", all.files = TRUE, full.names = TRUE)
  expect_length(bak, 1)
  expect_false("extra" %in% bc_schema(bak)$samples)
  expect_message(backwards_compatibility(dir), "up to date")
})

test_that("a failing step rolls back and keeps earlier steps", {
  dir <- bc_project()
  v <- RiboPilot:::current_schema_version()
  bad <- list(note = "boom", fn = function(con) {
    DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN half TEXT")
    stop("step failed")
  })
  steps <- RiboPilot:::db_migrations()
  local_mocked_bindings(db_migrations = function() c(steps, list(fake_step(), bad)))
  expect_error(suppressMessages(backwards_compatibility(dir)), "original database")
  h <- bc_history(dir)
  expect_equal(h$schema_version, c(v, v + 1L))
  s <- bc_schema(file.path(dir, ".sqlite"))$samples
  expect_true("extra" %in% s)
  expect_false("half" %in% s)
})

test_that("dry run reports pending steps and writes nothing", {
  dir <- bc_project()
  v <- RiboPilot:::current_schema_version()
  steps <- RiboPilot:::db_migrations()
  local_mocked_bindings(db_migrations = function() c(steps, list(fake_step())))
  expect_message(backwards_compatibility(dir, dry_run = TRUE), "add samples.extra")
  expect_equal(nrow(bc_history(dir)), 1)
  expect_length(list.files(dir, pattern = "\\.bak\\.", all.files = TRUE), 0)
})

test_that("unversioned and too-new databases stop with the fix", {
  dir <- bc_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE schema_history SET schema_version = 99")
  DBI::dbDisconnect(con)
  expect_error(backwards_compatibility(dir), "newer RiboPilot")

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "DROP TABLE schema_history")
  DBI::dbDisconnect(con)
  expect_error(backwards_compatibility(dir), "new_project")
})

test_that("the official container line is bumped and logged; custom images are not", {
  dir <- bc_project()
  cfg <- file.path(dir, ".config")
  writeLines(c("process {", "    container = 'macguigand/ribopilot:0.0.0.9001'", "}"), cfg)
  suppressMessages(backwards_compatibility(dir))
  want <- paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot"))
  expect_equal(readLines(cfg)[2], paste0("    container = '", want, "'"))
  expect_length(list.files(dir, pattern = "^\\.config\\.bak\\.", all.files = TRUE), 1)
  h <- bc_history(dir)
  expect_equal(h$event, c("created", "config"))

  writeLines(c("process {", "    container = 'me/custom:1'", "}"), cfg)
  expect_message(backwards_compatibility(dir), "up to date")
  expect_equal(readLines(cfg)[2], "    container = 'me/custom:1'")
})

test_that("a v1 database run through every step matches a fresh new_db()", {
  dir <- bc_project()
  old <- file.path(withr::local_tempdir(), ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), old)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  DBI::dbDisconnect(con)
  suppressMessages(backwards_compatibility(dirname(old)))
  expect_equal(bc_schema(old), bc_schema(file.path(dir, ".sqlite")))
  expect_equal(bc_views(old), bc_views(file.path(dir, ".sqlite")))
  expect_setequal(names(bc_views(old)), names(RiboPilot:::REVIEW_VIEWS_DDL))
})

test_that("run_app stops on a project that needs an update", {
  dir <- bc_project()
  v <- RiboPilot:::current_schema_version()
  steps <- RiboPilot:::db_migrations()
  local_mocked_bindings(db_migrations = function() c(steps, list(fake_step())))
  expect_error(run_app(dir, launch.browser = FALSE), "backwards_compatibility")
})

test_that("a v1 database with a lower-case biosample column upgrades cleanly", {
  old <- file.path(withr::local_tempdir(), ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), old)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "ALTER TABLE samples ADD COLUMN biosample TEXT")
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  DBI::dbDisconnect(con)
  suppressMessages(backwards_compatibility(dirname(old)))
  con <- DBI::dbConnect(RSQLite::SQLite(), old)
  on.exit(DBI::dbDisconnect(con))
  f <- DBI::dbListFields(con, "samples")
  expect_true("BioSample" %in% f)
  expect_false("biosample" %in% f)
  expect_equal(sum(tolower(f) == "biosample"), 1)
})

test_that("meta_apply_schema twice is a no-op", {
  dir <- bc_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  before <- DBI::dbListFields(con, "samples")
  expect_no_error(meta_apply_schema(con))
  expect_no_error(meta_apply_schema(con))
  expect_equal(DBI::dbListFields(con, "samples"), before)
})

