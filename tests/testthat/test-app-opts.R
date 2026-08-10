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

test_that("the mode picker shows Annotation Review but keys on Annotate", {
  # Renaming the tab must not rename the mode: every branch in app_server keys
  # on the VALUE, so a label-only change keeps them working.
  ui <- paste(deparse(RiboPilot:::app_ui), collapse = " ")
  # deparse backtick-quotes the name, so match on the pairing, not the source text
  expect_true(grepl("`Annotation Review` = \"Annotate\"", ui, fixed = TRUE))
  srv <- paste(deparse(RiboPilot:::app_server), collapse = " ")
  expect_true(grepl('mode == "Annotate"', srv, fixed = TRUE))
})

# ---- the assembly modal exposes the reference database ----

test_that("the assembly modal shows the preset's seed and label databases", {
  cur <- list(assemble_opts = "actinopterygii", cpus = 4, memory = 16,
    getOrganelle = "", assembler = "getorganelle", bowtie2 = "--local",
    spades = "--isolate", min_assembly_length = 0, max_paths = 0, max_scaffolds = 0,
    seeds_db = "/refs/seeds/actinopterygii_rrna_seeds.fasta",
    labels_db = "/refs/labels/actinopterygii_rrna_labels.fasta",
    clade_name = "Actinopterygii")
  html <- shiny::withReactiveDomain(shiny::MockShinySession$new(),
    as.character(RiboPilot:::assemble_opts_modal(cur, "actinopterygii", 10)))
  expect_true(grepl("ao_seeds", html, fixed = TRUE))
  expect_true(grepl("ao_labels", html, fixed = TRUE))
  expect_true(grepl("/refs/seeds/actinopterygii_rrna_seeds.fasta", html, fixed = TRUE))
  expect_true(grepl("/refs/labels/actinopterygii_rrna_labels.fasta", html, fixed = TRUE))
  # says where the paths came from
  expect_true(grepl("Reference tab for clade Actinopterygii", html, fixed = TRUE))
  # the stale claim that databases are set at project creation is gone
  expect_false(grepl("databases are set at project creation", html, fixed = TRUE))
})

test_that("a preset no reference group owns does not claim the Reference tab set it", {
  cur <- list(assemble_opts = "default", cpus = 4, memory = 16, getOrganelle = "",
    assembler = "getorganelle", bowtie2 = "--local", spades = "--isolate",
    min_assembly_length = 0, max_paths = 0, max_scaffolds = 0,
    seeds_db = NA_character_, labels_db = NA_character_)
  html <- shiny::withReactiveDomain(shiny::MockShinySession$new(),
    as.character(RiboPilot:::assemble_opts_modal(cur, "default", 1)))
  expect_false(grepl("Reference tab for clade", html, fixed = TRUE))
  expect_true(grepl("Set at project creation", html, fixed = TRUE))
})

test_that("opts_row takes the databases from the inputs, not the base row", {
  base <- list(seeds_db = "/old/seeds.fasta", labels_db = "/old/labels.fasta",
    target = NA_character_)
  input <- list(ao_cpus = 4, ao_mem = 16, ao_getorg = "", ao_assembler = "getorganelle",
    ao_bt2 = "--local", ao_spades = "--isolate", ao_minlen = 0, ao_maxpaths = 0,
    ao_maxscaf = 0, ao_seeds = "/new/seeds.fasta", ao_labels = "/new/labels.fasta")
  row <- RiboPilot:::opts_row("ao", "actinopterygii", input, base)
  expect_equal(row$seeds_db, "/new/seeds.fasta")
  expect_equal(row$labels_db, "/new/labels.fasta")
  # an emptied box clears the path rather than storing ""
  input$ao_seeds <- "  "
  expect_true(is.na(RiboPilot:::opts_row("ao", "x", input, base)$seeds_db))
})

test_that("saving a database that is missing or empty warns, naming the file", {
  d <- withr::local_tempdir()
  good <- file.path(d, "good.fasta"); writeLines(c(">18S", "ACGT"), good)
  empty <- file.path(d, "empty.fasta"); file.create(empty)
  w <- RiboPilot:::unusable_db_warnings(
    data.frame(seeds_db = file.path(d, "gone.fasta"), labels_db = empty,
      stringsAsFactors = FALSE))
  expect_length(w, 2L)
  expect_true(any(grepl("gone.fasta", w, fixed = TRUE)))
  expect_true(any(grepl("Label database", w, fixed = TRUE)))
  expect_length(RiboPilot:::unusable_db_warnings(
    data.frame(seeds_db = good, labels_db = good, stringsAsFactors = FALSE)), 0L)
  # nothing set is not a complaint
  expect_length(RiboPilot:::unusable_db_warnings(
    data.frame(seeds_db = NA_character_, labels_db = NA_character_,
      stringsAsFactors = FALSE)), 0L)
})
