test_that("the viewer payload carries genes, depth, and clipped-read share", {
  f <- data.frame(element = c("18S", "ITS1"), start = c(1, 11), stop = c(10, 18),
    strand = c("+", "-"), partial_start = c(1, 0), partial_stop = c(0, NA))
  p <- RiboPilot:::seqview_payload("acgtacgtacgtacgtacgt", f, depths = c(rep(20, 19), 5),
    clips = c(10, rep(0, 19)), selected = "ITS1")
  expect_equal(p$len, 20L)
  expect_equal(p$seq, "ACGTACGTACGTACGTACGT")
  expect_equal(p$selected, "ITS1")
  expect_equal(vapply(p$features, `[[`, "", "row"), c("18S", "ITS1"))
  expect_true(p$features[[1]]$partial5)
  expect_false(p$features[[2]]$partial3)
  expect_equal(p$features[[2]]$dir, "-")
  expect_equal(p$features[[1]]$color, "#0072B2")
  expect_equal(p$clip[1], 0.5)
  expect_true(is.na(p$clip[20]))
  expect_equal(p$cutoff, RiboPilot:::CHIMERA_CLIP_FRAC)
})

test_that("the viewer payload leaves out tracks it has no data for", {
  f <- data.frame(element = character(0), start = numeric(0), stop = numeric(0), strand = character(0))
  p <- RiboPilot:::seqview_payload("ACGT", f)
  expect_null(p$depth); expect_null(p$clip); expect_null(p$selected)
  expect_length(p$features, 0L)
})

test_that("both details windows show the viewer instead of the coverage plot", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  html <- paste(as.character(RiboPilot:::annotate_modal("S1", 1L, 1L)), collapse = " ")
  expect_match(html, "seqv_ann-canvas", fixed = TRUE)
  expect_false(grepl("coverage_plot", html, fixed = TRUE))
})

test_that("the Annotation Review viewer shows the copy and the Assemble viewer the original", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  sq <- strrep("ACGTACGTAC", 10)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 100 WHERE path = 1 AND scaffold = 1",
    params = list(sq))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,40,'+'), ('S1',1,1,'28S',61,90,'+')")
  DBI::dbExecute(con, "UPDATE coverage SET depths = ?, clips = ? WHERE path = 1 AND scaffold = 1",
    params = list(paste(1:100 * 10, collapse = ","), paste(c(rep(0, 20), 105, rep(0, 79)), collapse = ",")))
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 90)
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  feat <- function(el) Filter(function(f) f$row == el, sent$features)[[1]]

  RiboPilot:::seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", "S1", 1L, 1L, review = TRUE)
  expect_equal(sent$len, 70L)
  expect_equal(sent$seq, substr(sq, 21, 90))
  expect_equal(c(feat("18S")$pos1, feat("18S")$pos2), c(1L, 20L))
  expect_true(feat("18S")$partial5)    # cut by the trim
  expect_equal(c(feat("28S")$pos1, feat("28S")$pos2), c(41L, 70L))
  expect_false(feat("28S")$partial3)   # trimmed exactly at its end: not an assembly end
  expect_equal(sent$depth, as.integer(21:90 * 10))
  expect_equal(sent$clip[1], 0.5)
  expect_equal(sent$input, "seqv_ann_pick")

  RiboPilot:::seqview_send(session, "seqv_asm-canvas", NULL, "S1", 1L, 1L)
  expect_equal(sent$len, 100L)
  expect_equal(sent$seq, sq)
  expect_equal(c(feat("18S")$pos1, feat("18S")$pos2), c(1L, 40L))
  expect_true(feat("18S")$partial5)    # starts at the contig's first base
  expect_false(feat("28S")$partial3)
  expect_equal(sent$depth, as.integer(1:100 * 10))
  expect_equal(sent$clip[21], 0.5)
})

test_that("Annotation Review shows curation warnings only once the contig's checks are done", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand, status, flags)
    VALUES ('S1',1,1,'18S',1,40,'+','flag','original feature flag'), ('S1',1,1,'28S',61,90,'+','ok',NULL)")
  DBI::dbExecute(con, "UPDATE annotate SET n_flags = 1, curate_flags = 'original contig flag'
    WHERE path = 1 AND scaffold = 1")
  fixture_mark_done(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    expect_match(output$curate_notes$html, "1 contig warning(s): ", fixed = TRUE)
    expect_match(output$curate_notes$html, "original contig flag", fixed = TRUE)
    expect_equal(flag_rows()$flags, "original feature flag")

    # an edit makes a copy and clears its checks
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_match(output$curate_notes$html, "Checks not run. Curation warnings appear once they are done.", fixed = TRUE)
    expect_error(flag_rows())

    # the recheck finished
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, n_flags = 2, curate_flags = 'copy contig flag'")
    DBI::dbExecute(con, "UPDATE review_annotations SET status = 'flag', flags = 'copy feature flag' WHERE element = '28S'")
    refresh(refresh() + 1)
    session$flushReact()
    expect_match(output$curate_notes$html, "2 contig warning(s): ", fixed = TRUE)
    expect_match(output$curate_notes$html, "copy contig flag", fixed = TRUE)
    expect_equal(flag_rows()$element, "28S")
    expect_equal(flag_rows()$flags, "copy feature flag")
  })
})

test_that("each details window sends its own layer to the sequence viewer", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  sent <- list()
  local_mocked_bindings(seqview_send = function(session, canvas, input, id, path, scaffold, selected = NULL,
    review = FALSE) sent[[canvas]] <<- review)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()
  })
  expect_true(sent[["seqv_ann-canvas"]])
  expect_false(sent[["seqv_asm-canvas"]])
})

test_that("chimera_origin de-overlaps segments for drawing and keeps the host grey", {
  seg <- data.frame(start = c(1, 501, 2251), stop = c(750, 2500, 4500),
    taxon = c("nohit", "Actinopterygii", "Eurotiomycetes"), species = c(NA, "Salmo salar", "Penicillium sp."),
    pident = c(NA, 99, 97), n_windows = c(3, 10, 9),
    lineage_taxids = c(NA, "33208;7711;7898", "4751;4890;147545"), stringsAsFactors = FALSE)
  win <- data.frame(start = c(1, 2251), stop = c(500, 2750), pident = c(NA, 97))
  o <- RiboPilot:::chimera_origin(seg, win, "2375", "severe", NA, "33208;7711;7898", 4500L)
  expect_equal(vapply(o$origin, `[[`, 0, "pos2"), c(500, 2250, 4500))
  expect_true(o$origin[[2]]$host); expect_false(o$origin[[3]]$host)
  expect_equal(o$origin[[2]]$color, "#B3B3B3")
  expect_equal(o$origin[[1]]$color, "#FFFFFF")
  expect_equal(o$origin[[3]]$color, RiboPilot:::CHIMERA_LINEAGE_COLORS[1])
  expect_equal(o$junctions[[1]]$pos, 2375L); expect_true(o$junctions[[1]]$severe)
  expect_length(o$ident, 4500L)
  expect_true(is.na(o$ident[1])); expect_equal(o$ident[2300], 97)
})

test_that("chimera_origin with no segments or junctions sends nothing", {
  o <- RiboPilot:::chimera_origin(data.frame(), data.frame(), NA, "clean", NA, NA_character_, 100L)
  expect_length(o$origin, 0L); expect_length(o$junctions, 0L); expect_null(o$ident)
})

test_that("the lineage palette has three colours", {
  expect_equal(RiboPilot:::CHIMERA_LINEAGE_COLORS, c("#EEDD88", "#77AADD", "#AAAA00"))
})
