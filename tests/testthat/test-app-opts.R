# the modal's tool-help icons need a reactive domain to namespace their inputs
ao_modal_html <- function(cur) {
  shiny::withReactiveDomain(shiny::MockShinySession$new(),
    as.character(assemble_opts_modal(cur, choices = c("default"), n = 1)))
}

# a full assemble_opts row, as the modal gets it from the database
ao_cur <- function(...) {
  utils::modifyList(list(assemble_opts = "default", cpus = 4L, memory = 16L,
    getOrganelle = "", seeds_db = "s.fa", labels_db = "l.fa", target = NA_character_,
    max_paths = 0L, max_scaffolds = 0L, min_assembly_length = 0L,
    assembler = "getorganelle", bowtie2 = "--local", spades = "--isolate"), list(...))
}

test_that("the assembly modal offers an assembler picker with the three assemblers", {
  html <- ao_modal_html(ao_cur(assembler = "bait"))
  expect_true(grepl("ao_assembler", html, fixed = TRUE))
  for (a in c("getorganelle", "bait", "denovo")) {
    expect_true(grepl(sprintf('value="%s"', a), html, fixed = TRUE))
  }
  # the sample's current assembler is the selected option
  expect_true(grepl('value="bait" selected', html, fixed = TRUE))
})

test_that("the assembly modal defaults to getorganelle when the value is missing", {
  cur <- ao_cur(); cur$assembler <- NULL
  expect_true(grepl('value="getorganelle" selected', ao_modal_html(cur), fixed = TRUE))
})

test_that("the picker comes before the tool options it switches between", {
  html <- ao_modal_html(ao_cur())
  expect_lt(regexpr("ao_assembler", html, fixed = TRUE),
    regexpr("ao_getorg", html, fixed = TRUE))
})

test_that("the fallback assemblers get their own option boxes", {
  html <- ao_modal_html(ao_cur(assembler = "bait", bowtie2 = "--local", spades = "--isolate"))
  expect_true(grepl("ao_bt2", html, fixed = TRUE))
  expect_true(grepl("ao_spades", html, fixed = TRUE))
})

test_that("each tool option's help sits inside its own input container", {
  # so that hiding the input hides its help with it; a help paragraph outside
  # the container would be left behind when the assembler is switched
  # walk the tag tree for the input's own container, so a help paragraph that
  # merely sits next to the input in the rendered html does not pass
  container_of <- function(tag, id) {
    if (inherits(tag, "list")) {
      for (ch in tag) {
        got <- container_of(ch, id)
        if (!is.null(got)) return(got)
      }
      return(NULL)
    }
    if (!inherits(tag, "shiny.tag")) return(NULL)
    if (grepl("shiny-input-container", tag$attribs$class %||% "") &&
        grepl(sprintf('id="%s"', id), as.character(tag), fixed = TRUE)) {
      return(tag)
    }
    container_of(tag$children, id)
  }

  ui <- shiny::withReactiveDomain(shiny::MockShinySession$new(),
    assemble_opts_modal(ao_cur(), choices = c("default"), n = 1))
  for (id in c("ao_getorg", "ao_bt2", "ao_spades")) {
    box <- container_of(ui, id)
    expect_false(is.null(box), info = id)
    expect_true(grepl(sprintf('id="help_%s"', id), as.character(box), fixed = TRUE), info = id)
  }
})

test_that("assembler_fields names the inputs to show for each assembler", {
  expect_equal(assembler_fields("getorganelle"),
    list(show = "ao_getorg", hide = c("ao_bt2", "ao_spades")))
  expect_equal(assembler_fields("bait"),
    list(show = c("ao_bt2", "ao_spades"), hide = "ao_getorg"))
  # de-novo assembles all reads, so there is nothing to bait with
  expect_equal(assembler_fields("denovo"),
    list(show = "ao_spades", hide = c("ao_getorg", "ao_bt2")))
})

test_that("opts_row carries the picked assembler into the assemble_opts row", {
  input <- list(ao_cpus = 8, ao_mem = 32, ao_getorg = "-R 10", ao_minlen = 100,
    ao_maxpaths = 2, ao_maxscaf = 3, ao_assembler = "denovo",
    ao_bt2 = "--local", ao_spades = "--careful")
  row <- opts_row("ao", "custom", input, base = list(seeds_db = "s.fa", labels_db = "l.fa"))
  expect_equal(row$assemble_opts, "custom")
  expect_equal(row$assembler, "denovo")
  expect_equal(row$cpus, 8L)
  expect_equal(row$seeds_db, "s.fa")
})

test_that("opts_row falls back to getorganelle if the picker is missing", {
  input <- list(ao_cpus = 4, ao_mem = 16, ao_getorg = "", ao_minlen = 0,
    ao_maxpaths = 0, ao_maxscaf = 0)
  expect_equal(opts_row("ao", "custom", input, base = list())$assembler, "getorganelle")
})

test_that("opts_row still builds the other opts tables", {
  expect_equal(opts_row("po", "p", list(po_cpus = 2, po_mem = 4, po_fastp = "-x"),
    base = list())$fastp, "-x")
  expect_equal(opts_row("co", "c", list(co_ruleset = "strict", co_pid = 0.9,
    co_contam = TRUE), base = list())$check_contaminant, 1L)
})

test_that("annotation and curation presets are set from the Assemble surface", {
  # One workflow runs assembly AND annotation, so every option a run consumes
  # belongs on the surface the run is launched from (doc 15). The Annotate tab
  # is a review surface now and must not offer them.
  srv <- paste(deparse(RiboPilot:::app_server), collapse = "\n")

  # the Assemble-mode query supplies both presets...
  expect_true(grepl("AS annotate_opts", srv, fixed = TRUE))
  expect_true(grepl("AS curate_opts", srv, fixed = TRUE))
  # ...and they are still editable cells wired to the opts handlers
  expect_true(grepl("set_annotate_opts", srv, fixed = TRUE))
  expect_true(grepl("set_curate_opts", srv, fixed = TRUE))

  # ...and the Annotate-mode query no longer carries them. Anchored on that
  # query specifically: the Assemble one names the same columns in subselects,
  # so a whole-file search would pass while the tab still showed them.
  i <- regexpr("an.annotate_lock, an.annotate_switch", srv, fixed = TRUE)
  expect_true(i > 0)
  ann_q <- substr(srv, i, i + 700)
  expect_false(grepl("annotate_opts", ann_q, fixed = TRUE))
  expect_false(grepl("curate_opts", ann_q, fixed = TRUE))
})

test_that("the Annotate tab has no Update button", {
  ui <- paste(deparse(RiboPilot:::app_ui), collapse = "\n")
  i <- regexpr("input.mode == 'Annotate'", ui, fixed = TRUE)
  expect_true(i > 0)
  block <- substr(ui, i, i + 400)
  expect_false(grepl("run_btn", block, fixed = TRUE))
})
