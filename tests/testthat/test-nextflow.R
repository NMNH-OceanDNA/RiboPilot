# Nextflow version classification + NXF_VER pin selection (pure helpers).

test_that("nf_version_status classifies against the supported range", {
  f <- RiboPilot:::nf_version_status
  expect_equal(f(NA_character_), "missing")
  expect_equal(f("23.10.0"), "too_old")
  expect_equal(f("24.04.4"), "too_old")
  expect_equal(f("24.10.0"), "supported")   # lower bound inclusive
  expect_equal(f("25.04.8"), "supported")
  expect_equal(f("25.10.6"), "supported")
  expect_equal(f("26.0.0"), "too_new")      # broken-from inclusive
  expect_equal(f("26.1.0"), "too_new")
})

test_that("nf_pin_for pins in-range to itself and out-of-range to the max", {
  f <- RiboPilot:::nf_pin_for
  expect_true(is.na(f(NA_character_)))
  # supported -> preserve the original string, including a leading zero
  expect_equal(f("25.04.8"), "25.04.8")
  expect_equal(f("24.10.0"), "24.10.0")
  # out of range -> latest verified-good
  expect_equal(f("23.10.0"), RiboPilot:::NF_MAX_SUPPORTED)
  expect_equal(f("26.1.0"), RiboPilot:::NF_MAX_SUPPORTED)
})

test_that("nxf_ver_prefix is empty when Nextflow is missing", {
  # nf_pin_version(<bogus binary>) -> NA -> empty prefix, so nothing is pinned.
  expect_equal(RiboPilot:::nxf_ver_prefix("definitely-not-a-real-nextflow-bin"), "")
})

# WF1 pipeline wiring: the screen must run per sample and triage must wait for
# it, so a sample is triaged only once both its curation and its screen exist.

main_nf <- readLines(system.file("nextflow", "main.nf", package = "RiboPilot"))
pipeline_r_nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))

test_that("main.nf includes screen from pipeline_r.nf", {
  include_line <- grep("^include \\{.*\\} from './modules/pipeline_r.nf'$", main_nf, value = TRUE)
  expect_length(include_line, 1)
  expect_match(include_line, "\\bscreen\\b")
})

test_that("screen process exists in pipeline_r.nf with errorStrategy 'ignore'", {
  start <- grep("^process screen \\{", pipeline_r_nf)
  expect_length(start, 1)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]
  expect_true(any(grepl("errorStrategy 'ignore'", body)))
  expect_true(any(grepl("cache false", body)))
})

test_that("triage's input joins on screen.out", {
  start <- grep("^\\s*triage\\(", main_nf)
  expect_length(start, 1)
  end <- start - 1 + which(grepl("^\\s*\\}\\s*$", main_nf[start:length(main_nf)]) |
    grepl("\\)\\)$", main_nf[start:length(main_nf)]))[1]
  block <- paste(main_nf[start:end], collapse = "\n")
  expect_match(block, "\\.join\\(screen\\.out")
  # screen carries errorStrategy 'ignore', so an inner join would silently drop
  # a fully curated sample whenever its screen (or its upstream BLAST) failed.
  # remainder: true is what keeps that sample in triage; guard it specifically,
  # not just the presence of a join.
  expect_match(block, "\\.join\\(screen\\.out\\.map \\{ id -> tuple\\(id, true\\) \\}, remainder: true\\)")
})

# The assembly fallback is a per-sample assembler inside WF1, not a second
# workflow. It used to be a second workflow, and the copy drifted out of step
# with main.nf until it no longer ran. These assertions read the workflow files
# as text because the defect that motivated them was a channel arity mismatch,
# which no R-level test can see.

test_that("query_samples selects the assembler and declares it in the header", {
  start <- grep("^process query_samples \\{", pipeline_r_nf)
  expect_length(start, 1)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]

  hdr <- grep("printf 'ID", body, value = TRUE)
  expect_length(hdr, 1)
  expect_match(hdr, "assembler")
  expect_true(any(grepl("AS assembler", body)))

  # The header is hand-written because sqlite3 -header emits nothing for zero
  # rows, so it and the select list must agree on column count or splitCsv
  # misreads every field.
  n_hdr <- length(strsplit(sub(".*printf '", "", sub("\\\\n' > .*", "", hdr)), "\\\\t")[[1]])
  expect_equal(n_hdr, 12L)
})

test_that("main.nf branches on the assembler and merges both arms into score", {
  br <- grep("asm_branched = asm_in\\.branch", main_nf)
  expect_length(br, 1)
  block <- paste(main_nf[br:(br + 4)], collapse = "\n")
  expect_match(block, "go: it\\[7\\] == 'getorganelle'")
  # `fb: true`, not an equality test: branch silently DROPS an item matching no
  # condition, so an unrecognised assembler would vanish from the run.
  expect_match(block, "fb: true")

  expect_true(any(grepl("^\\s*assemble\\(asm_branched\\.go", main_nf)))
  expect_true(any(grepl("^\\s*fallback_assemble\\(asm_branched\\.fb", main_nf)))
  expect_true(any(grepl("score\\(assemble\\.out\\.paths\\.mix\\(fallback_assemble\\.out\\.paths\\)\\)", main_nf)))
})

test_that("main.nf includes fallback_assemble and still passes seeds to GetOrganelle", {
  expect_true(any(grepl("^include \\{ fallback_assemble \\} from './modules/fallback.nf'$", main_nf)))
  go_call <- grep("assemble\\(asm_branched\\.go", main_nf)
  expect_length(go_call, 1)
  block <- paste(main_nf[go_call:(go_call + 1)], collapse = " ")
  expect_match(block, "tuple\\(id, reads, cpus, mem, go, seeds, labels\\)")
})

test_that("there is ONE workflow: no second file left to drift", {
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  # main_reference.nf is the reference-database workflow, not a second WF1
  wf <- setdiff(list.files(nf_dir, pattern = "\\.nf$"),
    c("nextflow.config", "main_reference.nf"))
  expect_setequal(wf, "main.nf")
  expect_false(file.exists(file.path(nf_dir, "main_fallback.nf")))
  expect_false(file.exists(file.path(nf_dir, "main_annotate.nf")))

  # each entry point went with the workflow that called it
  expect_false(any(grepl("query_fallback_samples", pipeline_r_nf)))
  expect_false(any(grepl("query_annotate_samples", pipeline_r_nf)))

  # doc 08 open question 2 is answered: the stitch machinery is gone, not merely
  # unreachable. Both the process and the script it called must stay gone.
  for (nm in list(main = main_nf, modules = pipeline_r_nf)) {
    expect_false(any(grepl("\\bstitch\\b", nm[!grepl("^\\s*//", nm)])))
  }
  expect_false(file.exists(system.file("nextflow", "bin", "rp_stitch.R", package = "RiboPilot")))
})

test_that("WF1 mixes re-queued units into the annotation chain", {
  expect_true(any(grepl("query_annotate_units", pipeline_r_nf)))
  code <- main_nf[!grepl("^\\s*//", main_nf)]
  expect_true(any(grepl("query_annotate_units", code)))
  # the two arms must merge, or re-queued units get a second chain of their own
  # and the drift this sub-project removed comes straight back
  expect_true(any(grepl("annotate_tools(units.mix(requeued))", code, fixed = TRUE)))
})

test_that("the wiring the annotate-only arm depends on is still there", {
  # An annotate-only run fetches no references and runs no screen. Both of these
  # keep it from stalling or being dropped, and neither is obviously load-bearing
  # when read in isolation, so pin them: see dev/design/15_WF2_DECOMMISSION.md.
  code <- paste(main_nf[!grepl("^\\s*//", main_nf)], collapse = " ")
  expect_true(grepl("ifEmpty('ready')", code, fixed = TRUE))
  expect_true(grepl("remainder: true", code, fixed = TRUE))
})

# The annotate-only arm's selection is SQL inside a Nextflow process, so the
# suite would normally only see it as text. Extract that SQL and run it against a
# fixture database: the race it prevents (a sample queued for BOTH assembly and
# annotation being annotated twice in one run) is invisible to a text assertion.

annotate_units_sql <- function() {
  start <- grep("^process query_annotate_units \\{", pipeline_r_nf)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]
  sql <- body[grep("SELECT an.ID", body):grep(";\" >> annotate_units.tsv", body, fixed = TRUE)]
  sql <- sub(";\" >> annotate_units.tsv", "", paste(sql, collapse = "\n"), fixed = TRUE)
  gsub("\\$\\{params\\.db\\}", "", sql)
}

mk_units_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT,
    annotate_switch INT, annotate_lock INT, annotate_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate_opts(annotate_opts TEXT, cpus INT, memory INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, ignore INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, assemble_lock INT)")
  DBI::dbExecute(con, "INSERT INTO annotate_opts VALUES('default',4,8)")
  con
}

add_unit <- function(con, id, scaffold = 1, sw = 1, lock = 0, ignore = 0,
                     asm_sw = 2, asm_lock = 0) {
  DBI::dbExecute(con, "INSERT OR REPLACE INTO assemble VALUES(?,?,?)",
    params = list(id, asm_sw, asm_lock))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES(?,1,?,?)",
    params = list(id, scaffold, ignore))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES(?,1,?,?,?,'default')",
    params = list(id, scaffold, sw, lock))
}

test_that("query_annotate_units selects re-queued units and skips the rest", {
  con <- mk_units_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, "QUEUED")                                  # queued, assembly done
  add_unit(con, "DONE", sw = 2)                            # nothing outstanding
  add_unit(con, "HELD", sw = 0)                            # deliberately held
  add_unit(con, "LOCKED", lock = 1)                        # annotation locked
  add_unit(con, "DROPPED", ignore = 1)                     # parked scaffold

  got <- DBI::dbGetQuery(con, annotate_units_sql())$ID
  expect_equal(got, "QUEUED")
})

test_that("query_annotate_units excludes samples this run will assemble", {
  # Assembly re-seeds a sample's annotate units itself, so a sample queued for
  # both would otherwise be annotated twice in one run, racing on the same rows.
  con <- mk_units_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, "REASSEMBLING", asm_sw = 1)                # queued for assembly
  add_unit(con, "ANNOTATE_ONLY", asm_sw = 2)

  got <- DBI::dbGetQuery(con, annotate_units_sql())$ID
  expect_equal(got, "ANNOTATE_ONLY")

  # ...but a sample queued for assembly whose assembly is LOCKED is not going to
  # be assembled, so its queued units are ours to annotate
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'REASSEMBLING'")
  expect_setequal(DBI::dbGetQuery(con, annotate_units_sql())$ID,
    c("REASSEMBLING", "ANNOTATE_ONLY"))
})

test_that("query_annotate_units does not require a locked assembly", {
  # WF2 gated on assemble_lock = 1 so it could not run before review. With one
  # workflow a unit is annotate-only work because a user queued it.
  con <- mk_units_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, "UNREVIEWED", asm_sw = 2, asm_lock = 0)
  expect_equal(DBI::dbGetQuery(con, annotate_units_sql())$ID, "UNREVIEWED")
})
