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
    assembler = "getorganelle"), list(...))
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

test_that("opts_row carries the picked assembler into the assemble_opts row", {
  input <- list(ao_cpus = 8, ao_mem = 32, ao_getorg = "-R 10", ao_minlen = 100,
    ao_maxpaths = 2, ao_maxscaf = 3, ao_assembler = "denovo")
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
