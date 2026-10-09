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

test_that("nf_parse_version reads the version line and drops edge suffixes", {
  expect_equal(RiboPilot:::nf_parse_version(c("  N E X T F L O W", "  version 25.10.6 build 10000")), "25.10.6")
  expect_equal(RiboPilot:::nf_parse_version(c("version 24.10.0-edge build 1")), "24.10.0")
  expect_true(is.na(RiboPilot:::nf_parse_version("no version here")))
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
  expect_equal(n_hdr, 19L)
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
  # every assembly source merges into one score() call
  sc <- grep("score\\(assemble\\.out\\.paths", main_nf)
  expect_length(sc, 1)
  block <- paste(main_nf[sc:(sc + 3)], collapse = " ")
  expect_match(block, "mix\\(fallback_assemble\\.out\\.paths\\)")
  expect_match(block, "mix\\(find_rrna\\.out\\.paths\\)")
  expect_match(block, "mix\\(stage_assembly\\.out\\.paths\\)")
})

test_that("main.nf includes fallback_assemble and still passes seeds to GetOrganelle", {
  expect_true(any(grepl("^include \\{ fallback_assemble \\} from './modules/fallback.nf'$", main_nf)))
  go_call <- grep("assemble\\(asm_branched\\.go", main_nf)
  expect_length(go_call, 1)
  block <- paste(main_nf[go_call:(go_call + 1)], collapse = " ")
  expect_match(block, "tuple\\(id, reads, cpus, mem, go, file\\(seeds\\), file\\(labels\\)\\)")
})

test_that("there is ONE workflow: no second file left to drift", {
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  # main_reference.nf builds reference databases and main_recheck.nf rechecks
  # review copies; neither is a second WF1
  wf <- setdiff(list.files(nf_dir, pattern = "\\.nf$"),
    c("nextflow.config", "main_reference.nf", "main_recheck.nf"))
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
  DBI::dbExecute(con, "CREATE TABLE blast_searched(ID TEXT, path INT, scaffold INT)")
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

# ---- reference databases are mounted into the container ----

mk_mount_project <- function(dir, seeds, labels, engine = "docker") {
  writeLines(c(RiboPilot:::container_engine_block(engine),
    "params {", "  rawDir = '/data/raw'", "}"), file.path(dir, ".config"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (assemble_opts TEXT, seeds_db TEXT, labels_db TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('actinopterygii', ?, ?)",
    params = list(seeds, labels))
  DBI::dbDisconnect(con)
  dir
}

test_that("a database outside the project is bind-mounted for the pipeline run", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  dir.create(file.path(refs, "seeds")); dir.create(file.path(refs, "labels"))
  s <- file.path(refs, "seeds", "x_seeds.fasta"); l <- file.path(refs, "labels", "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l)
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  i <- which(a == "-process.containerOptions")
  expect_length(i, 1L)
  expect_true(grepl(file.path(refs, "seeds"), a[i + 1], fixed = TRUE))
  expect_true(grepl(file.path(refs, "labels"), a[i + 1], fixed = TRUE))
})

test_that("a database inside the project needs no extra mount", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  dir.create(file.path(proj, "refs"))
  s <- file.path(proj, "refs", "x_seeds.fasta"); l <- file.path(proj, "refs", "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l)
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_false(any(grepl(normalizePath(proj), a[which(a == "-process.containerOptions") + 1] %||% "",
    fixed = TRUE)))
})

test_that("a config with no container engine gets no mount option", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  s <- file.path(refs, "x_seeds.fasta"); l <- file.path(refs, "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  writeLines(c("process { executor = 'local' }", "params {", "  rawDir = '/data/raw'", "}"),
    file.path(proj, ".config"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (assemble_opts TEXT, seeds_db TEXT, labels_db TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('a', ?, ?)", params = list(s, l))
  DBI::dbDisconnect(con)
  expect_false("-process.containerOptions" %in% RiboPilot:::nf_command(proj, resume = FALSE))
})

test_that("singularity gets --bind, not -v", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  s <- file.path(refs, "x_seeds.fasta"); l <- file.path(refs, "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l, engine = "singularity")
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  v <- a[which(a == "-process.containerOptions") + 1]
  expect_true(grepl("--bind", v, fixed = TRUE))
  expect_false(grepl("-v ", v, fixed = TRUE))
})

test_that("directories the config already mounts are not mounted twice", {
  skip_if_not_installed("withr")
  # a fixture project has no reference_config row, so ref_config() would
  # otherwise fall back to this machine's real, populated default ref dir
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  s <- file.path(refs, "x_seeds.fasta"); l <- file.path(refs, "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  # the mount the config carries names refs explicitly, as new_project() writes it
  writeLines(c("docker {", "  enabled = true",
    sprintf('  runOptions = "-v ${params.projectDir}:${params.projectDir} -v %s:%s"', refs, refs),
    "}", "params {", "  rawDir = '/data/raw'", "}"), file.path(proj, ".config"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (assemble_opts TEXT, seeds_db TEXT, labels_db TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('a', ?, ?)", params = list(s, l))
  DBI::dbDisconnect(con)
  expect_false("-process.containerOptions" %in% RiboPilot:::nf_command(proj, resume = FALSE))
})

test_that("a config mount is parsed whole, not truncated at a letter s", {
  # the default regex engine reads [^\\s:] as "not backslash, s or colon", which
  # silently cut /home/u/Documents/... down to /home/u/Document
  cfg <- c("docker {", "  enabled = true",
    '  runOptions = "-v /home/u/Documents/refs:/home/u/Documents/refs"', "}")
  expect_equal(RiboPilot:::config_mounted_dirs(cfg), "/home/u/Documents/refs")
})

test_that("seed and label subfolders are each mounted, and neither is dropped", {
  skip_if_not_installed("withr")
  # a fixture project has no reference_config row, so ref_config() would
  # otherwise fall back to this machine's real, populated default ref dir
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  dir.create(file.path(refs, "seeds")); dir.create(file.path(refs, "labels"))
  s <- file.path(refs, "seeds", "x_seeds.fasta"); l <- file.path(refs, "labels", "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l)
  v <- RiboPilot:::nf_command(proj, resume = FALSE)
  v <- v[which(v == "-process.containerOptions") + 1]
  expect_true(grepl(file.path(refs, "seeds"), v, fixed = TRUE))
  expect_true(grepl(file.path(refs, "labels"), v, fixed = TRUE))
  expect_equal(lengths(regmatches(v, gregexpr("-v ", v, fixed = TRUE))), 2L)
})

test_that("a folder collapses its own subfolders when both are named", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  dir.create(file.path(d, "seeds"))
  expect_equal(RiboPilot:::drop_covered(c(d, file.path(d, "seeds")), c(d, file.path(d, "seeds"))), d)
})

# nf_ref_dir must use the reference folder only when reference_config.ref_dir is
# explicitly set, or an unconfigured project silently picks up whatever this
# machine's default reference folder resolves to.

mk_ref_config_project <- function(dir, ref_dir_value, blastdb = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records, blastdb)
    VALUES (1, ?, 1000, 25000, ?)", params = list(ref_dir_value, blastdb))
  DBI::dbDisconnect(con)
  dir
}

test_that("nf_ref_dir is blank for a project with no explicit ref_dir", {
  skip_if_not_installed("withr")
  # unset so a real machine default (which ribopilot_ref_dir() always resolves
  # to something for) cannot mask the fix
  withr::local_envvar(RIBOPILOT_REF_DIR = "")
  d <- withr::local_tempdir()
  mk_ref_config_project(d, NA_character_)
  expect_equal(RiboPilot:::nf_ref_dir(d), "")
})

test_that("nf_ref_dir returns the configured folder when explicitly set", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  mk_ref_config_project(d, refs)
  expect_equal(RiboPilot:::nf_ref_dir(d), path.expand(refs))
})

# Failure handling. A pipeline that reports success while quietly dropping work
# is worse than one that stops, so pin the three guarantees that prevent it.

test_that("every DB-touching process states an error strategy", {
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  mods <- list.files(file.path(nf_dir, "modules"), pattern = "\\.nf$", full.names = TRUE)
  for (f in c(mods, file.path(nf_dir, "main_reference.nf"))) {
    lines <- readLines(f)
    starts <- grep("^process\\s+\\S+\\s*\\{", lines)
    for (s in starts) {
      nm <- sub("^process\\s+(\\S+)\\s*\\{.*$", "\\1", lines[s])
      body <- lines[s:min(length(lines), s + 20)]
      body <- body[seq_len(which(c(grepl("^\\s*(input|script|output):", body), TRUE))[1] - 1)]
      expect_true(any(grepl("errorStrategy", body)),
        info = paste("process", nm, "in", basename(f), "has no errorStrategy"))
    }
  }
})

test_that("sqlite3 CLI calls set a busy timeout", {
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  mods <- list.files(file.path(nf_dir, "modules"), pattern = "\\.nf$", full.names = TRUE)
  for (f in mods) {
    lines <- readLines(f)
    calls <- grep("sqlite3\\s+(-|')", lines[!grepl("^\\s*(//|#)", lines)], value = TRUE)
    for (cl in calls) {
      # the CLI default busy timeout is 0, so a lock race fails the task outright
      expect_match(cl, "sqlite3 -cmd '\\.timeout", info = paste(basename(f), cl))
    }
  }
})

test_that("the shipped config sets pipefail and throttles the NCBI processes", {
  cfg <- readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot"))
  txt <- paste(cfg, collapse = " ")
  # without pipefail, `bowtie2 | samtools sort` reports samtools' exit status
  expect_match(txt, "pipefail")
  expect_match(txt, "blast_genbank\\|blast_ref_fetch")
  expect_match(txt, "maxForks")
})

test_that("an out-of-memory kill retries with more memory instead of dropping", {
  txt <- paste(readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot")),
    collapse = " ")
  expect_match(txt, "assemble\\|fallback_assemble\\|preprocess\\|annotate_tools")
  for (code in c(104, 134, 137, 139, 140, 143, 247)) {
    expect_match(txt, paste0("\\b", code, "\\b"))
  }
  expect_match(txt, "'retry' : 'ignore'", fixed = TRUE)
  expect_match(txt, "maxRetries = 3", fixed = TRUE)
  # without task.attempt in the memory closure the retry asks for the same
  # allocation that just failed. Derive the process list from the selector
  # rather than hardcoding it, so adding a process to the retry policy cannot
  # quietly skip this check.
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  pipeline_cfg <- readLines(file.path(nf_dir, "nextflow.config"))
  sel <- sub(".*withName: *'([^']+)'.*", "\\1",
    grep("withName.*assemble", pipeline_cfg, value = TRUE)[1])
  procs <- strsplit(sel, "|", fixed = TRUE)[[1]]
  expect_gt(length(procs), 1)
  mods <- list.files(file.path(nf_dir, "modules"), pattern = "\\.nf$", full.names = TRUE)
  for (p in procs) {
    mem <- NULL
    for (f in mods) {
      lines <- readLines(f)
      s <- grep(paste0("^process\\s+", p, "\\s*\\{"), lines)
      if (!length(s)) next
      body <- lines[s:min(length(lines), s + 20)]
      body <- body[seq_len(which(c(grepl("^\\s*(input|script|output):", body), TRUE))[1] - 1)]
      mem <- grep("^\\s*memory \\{", body, value = TRUE)
    }
    expect_length(mem, 1)
    expect_match(mem, "mem", info = p)          # DB value is still the starting point
    expect_match(mem, "task.attempt", fixed = TRUE, info = p)
  }
})

test_that("a sample with no seeds database goes to the fallback, not to file('')", {
  # sqlite3 renders a SQL NULL as '', and file('') throws during channel
  # construction, killing the run before any task is submitted
  txt <- paste(main_nf, collapse = "\n")
  expect_match(txt, "asm_opts = rows.map \\{ r -> tuple\\(r.ID, r.asm_cpus, r.asm_mem, r.getOrganelle, live_db\\(r.seeds_db\\), live_db\\(r.labels_db\\), r.assembler, r.bowtie2, r.spades\\)")
  expect_match(txt, "go: it[7] == 'getorganelle' && it[5] && it[6]", fixed = TRUE)
  expect_match(txt, "fb: true", fixed = TRUE)
})

test_that("a seeds path that no longer exists goes to the fallback too", {
  # a dead path used to be staged as a dangling symlink; GetOrganelle then exits
  # 0 after "seeds not found" and the sample reports success with no assembly
  txt <- paste(main_nf, collapse = "\n")
  expect_match(txt, "def live_db(p)", fixed = TRUE)
  expect_match(txt, "file(p).exists()", fixed = TRUE)
})

test_that("the Nextflow version probe runs once per session", {
  # a cold launcher downloads its dependencies before answering, so attaching
  # the package must not pay for the probe more than once
  rm(list = ls(RiboPilot:::.nf_probe_cache), envir = RiboPilot:::.nf_probe_cache)
  on.exit(rm(list = ls(RiboPilot:::.nf_probe_cache), envir = RiboPilot:::.nf_probe_cache))
  calls <- 0L
  local_mocked_bindings(
    nf_installed_version = function(...) { calls <<- calls + 1L; NA_character_ },
    .package = "RiboPilot")
  expect_true(is.na(RiboPilot:::nf_version_session()))
  expect_true(is.na(RiboPilot:::nf_version_session()))
  expect_equal(calls, 1L)
})

test_that("the API key is not interpolated into task scripts", {
  for (f in c("modules/blast.nf", "modules/blast_ref.nf", "main_reference.nf")) {
    txt <- readLines(system.file("nextflow", f, package = "RiboPilot"), warn = FALSE)
    expect_false(any(grepl("export NCBI_API_KEY='\\$\\{params", txt, perl = TRUE)),
      info = paste(f, "interpolates the API key into the task script"))
  }
})

test_that("efetch failures are not swallowed", {
  txt <- readLines(system.file("nextflow", "modules/blast_ref.nf", package = "RiboPilot"),
    warn = FALSE)
  efetch <- grep("efetch -db nuccore", txt, value = TRUE)
  expect_true(length(efetch) > 0)
  expect_false(any(grepl("\\|\\| true", efetch)),
    info = "efetch still ends with || true")
})

test_that("nextflow.config still supplies NCBI_API_KEY from params", {
  # a static grep for the removed exports gives false confidence on its own;
  # this pins the replacement so the key cannot be silently dropped again
  txt <- readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot"),
    warn = FALSE)
  line <- grep("NCBI_API_KEY\\s*[:=]", txt, value = TRUE)
  expect_length(line, 1)
  expect_match(line, "ncbi_api_key", fixed = TRUE)
})

test_that("with no NCBI key the variable is left out, and tasks read it unset-safe", {
  txt <- readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot"), warn = FALSE)
  k <- grep("NCBI_API_KEY\\s*:", txt)
  expect_false(grepl("?: ''", txt[k], fixed = TRUE))
  close <- grep("^\\].findAll \\{ it.value \\}", txt)
  expect_length(close, 1)
  expect_gt(close, k)
  nf <- readLines(system.file("nextflow", "modules/blast_ref.nf", package = "RiboPilot"), warn = FALSE)
  expect_false(any(grepl("\\$NCBI_API_KEY\"", nf, fixed = TRUE)))
})

test_that("nextflow.config forwards RIBOPILOT_NO_LOCAL_BLAST to tasks", {
  # docker tasks are launched with only NXF_TASK_WORKDIR forwarded, so a
  # shell-set opt-out reaches the container only if it is named in env{}
  txt <- readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot"),
    warn = FALSE)
  line <- grep("RIBOPILOT_NO_LOCAL_BLAST\\s*:", txt, value = TRUE)
  expect_length(line, 1)
  expect_match(line, "System.getenv('RIBOPILOT_NO_LOCAL_BLAST')", fixed = TRUE)
  # part of the env map the task is launched with, not stranded elsewhere. The
  # map form (env = [...]) is what lets the key be omitted entirely when unset:
  # an empty or null value makes Nextflow warn once per task.
  env_open <- grep("^env\\s*=\\s*\\[", txt)
  expect_length(env_open, 1)
  expect_gt(grep("RIBOPILOT_NO_LOCAL_BLAST\\s*:", txt), env_open)
  expect_match(paste(txt, collapse = "\n"), "\\+ \\(System.getenv\\('RIBOPILOT_NO_LOCAL_BLAST'\\)")
})

test_that("blast.nf gates re-search on blast_searched.query_md5, not a blast_hits COUNT", {
  # pins the fix: re-searching used to be gated on a COUNT(*) of blast_hits,
  # which re-submitted a unit forever once it found nothing. A revert to that
  # must not pass silently.
  path <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(path)) path <- "inst/nextflow/modules/blast.nf"
  txt <- readLines(path, warn = FALSE)
  code <- paste(txt, collapse = "\n")
  expect_match(code, "SELECT query_md5 FROM blast_searched", fixed = TRUE)
  expect_false(grepl("COUNT\\(\\*\\)[^\\n]*blast_hits|blast_hits[^\\n]*COUNT\\(\\*\\)", code),
    info = "blast.nf still gates re-search on a COUNT of blast_hits")
})

test_that("main.nf declares a ref_dir param", {
  p <- system.file("nextflow", "main.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/main.nf"
  expect_true(any(grepl("params.ref_dir", readLines(p, warn = FALSE), fixed = TRUE)))
})

test_that("the fetch module primes from the store before fetching", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- readLines(p, warn = FALSE)
  prime <- grep("rp_ref_prime.R", txt)
  accs <- grep("SELECT DISTINCT c.accession", txt)
  expect_length(prime, 1L)
  expect_true(prime < accs[1],
    info = "priming must run before the un-fetched accession list is computed")
})

test_that("main.nf declares a blastdb param", {
  p <- system.file("nextflow", "main.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/main.nf"
  expect_true(any(grepl("params.blastdb", readLines(p, warn = FALSE), fixed = TRUE)))
})

test_that("an empty ribopilot_src argument does not become the string true", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(proj, ".config"))
  a <- RiboPilot:::nf_command(proj, resume = FALSE, ribopilot_src = "")
  v <- a[which(a == "--ribopilot_src") + 1]
  expect_true(all(nzchar(v)))
})

test_that("an empty minDepth in the config falls back to 0", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  writeLines(c("params {", "  rawDir = '/data/raw'", "  minDepth = ''", "}"),
    file.path(proj, ".config"))
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_equal(a[which(a == "--minDepth") + 1], "0")
})

test_that("nf_blastdb is blank when the project configures no panel", {
  d <- withr::local_tempdir()
  proj <- mk_ref_config_project(d, NA_character_)
  expect_equal(RiboPilot:::nf_blastdb(proj), "")
})

test_that("nf_blastdb returns the configured panel prefix", {
  d <- withr::local_tempdir()
  proj <- mk_ref_config_project(d, NA_character_, blastdb = "/opt/rp/ribopilot_rrna")
  expect_equal(RiboPilot:::nf_blastdb(proj), "/opt/rp/ribopilot_rrna")
})

test_that("building the command preview does not record provenance; launching does", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(d, ".config"))
  panel_dir <- withr::local_tempdir()
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(panel_dir, "manifest.json"))
  prefix <- file.path(panel_dir, "ribopilot_rrna")
  file.create(paste0(prefix, ".nin"))
  mk_ref_config_project(d, NA_character_, blastdb = prefix)

  RiboPilot:::nf_command(d, resume = FALSE)   # the preview path: must record nothing
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  expect_equal(RiboPilot:::blastdb_provenance(con), "")
  DBI::dbDisconnect(con)

  RiboPilot:::nf_record_blastdb_launch(d)   # the launch path: must record it
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  expect_equal(RiboPilot:::blastdb_provenance(con), "ribopilot_rrna 2026.08.1")
  DBI::dbDisconnect(con)
})

test_that("blast.nf searches locally when a panel is configured", {
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("params.blastdb", txt, fixed = TRUE))
  expect_true(grepl("rp_dbarg.R", txt, fixed = TRUE))
  expect_true(grepl("-remote", txt, fixed = TRUE),
    info = "the remote path must remain for projects with no panel")
})

test_that("blast.nf records the search source", {
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  lines <- readLines(p, warn = FALSE)
  call_line <- grep("rp_searched.R", lines, fixed = TRUE, value = TRUE)
  expect_length(call_line, 1)
  # UNIT_SRC (not the task-level SRC) so one unit's local failure can't leak
  # into the source recorded for the next; UNIT_PANEL follows it so a unit
  # that fell back to remote does not get credited with a local panel
  expect_true(endsWith(trimws(call_line), '"\\$UNIT_SRC" "\\$UNIT_PANEL"'))
})

test_that("blast.nf claims a fixed 4 cpus, not a params.blastdb-gated closure", {
  # the cpus{} closure evaluates before the script block runs, so it cannot
  # see the local/remote decision rp_dbarg.R makes at task time; params.blastdb
  # also no longer reliably reflects "no panel configured" (nextflow coerces
  # an empty --blastdb '' CLI value to the string "true"), so a conditional
  # here was never correct
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  lines <- readLines(p, warn = FALSE)
  cpus_line <- grep("^\\s*cpus\\b", lines, value = TRUE)
  expect_length(cpus_line, 1)
  expect_match(cpus_line, "^\\s*cpus\\s+4\\s*$")
})

test_that("rp_dbarg.R does not suppress the decoy warnings", {
  p <- system.file("nextflow", "bin", "rp_dbarg.R", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/bin/rp_dbarg.R"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_false(grepl("suppressWarnings", txt, fixed = TRUE),
    info = "decoy not-found/whitespace warnings must reach the task log")
})

test_that("rp_dbarg.R falls back to discovery when given no panel", {
  txt <- readLines(system.file("nextflow/bin/rp_dbarg.R", package = "RiboPilot"), warn = FALSE)
  s <- paste(txt, collapse = "\n")
  expect_match(s, "discover_panel", fixed = TRUE)
  # an empty first argument must not be treated as a usable prefix
  expect_no_match(s, "rp_args(a, 1", fixed = TRUE)
})

test_that("rp_dbarg.R discovers a panel from RIBOPILOT_DB_DIR with no argument", {
  skip_on_cran()
  skip_if(!nzchar(Sys.which("Rscript")))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  d <- file.path(root, "2026.08.1"); dir.create(d)
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", file.path(d, paste0("ribopilot_rrna.", e)))
  script <- system.file("nextflow/bin/rp_dbarg.R", package = "RiboPilot")
  # NOTE: not getwd() -- testthat runs tests with the working directory set to
  # tests/testthat, not the package root, so getwd() here would send rp_env.R
  # a non-source RIBOPILOT_SRC and it would silently fall back to whatever
  # RiboPilot happens to be installed (which may predate discover_panel).
  out <- suppressWarnings(system2("Rscript", c(script, "''"),
    env = c(paste0("RIBOPILOT_DB_DIR=", root), paste0("RIBOPILOT_SRC=", RiboPilot:::ribopilot_src_default())),
    stdout = TRUE, stderr = FALSE))
  expect_match(paste(out, collapse = " "), "2026.08.1", fixed = TRUE)
})

test_that("blast.nf falls back to remote when a configured local panel fails inside blastn", {
  # A panel that exists and passes rp_dbarg.R's checks (no missing/whitespace
  # path) can still fail INSIDE blastn, e.g. built without -parse_seqids. That
  # used to `exit 1` and lose the unit entirely (no hits, no blast_searched
  # row). Runs the real blast.nf process end to end with a stub blastn on
  # PATH, so this is a genuine behavioural check, not a text grep.
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  blast_nf <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(blast_nf)) blast_nf <- file.path(src, "inst", "nextflow", "modules", "blast.nf")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src, "inst", "nextflow", "bin")

  root <- withr::local_tempdir()
  fake_bin <- file.path(root, "fakebin"); dir.create(fake_bin)
  panel_dir <- file.path(root, "panel"); dir.create(panel_dir)
  wf_dir <- file.path(root, "wf"); dir.create(wf_dir)
  work_dir <- file.path(root, "work")

  # stub blastn: fails any local (-db, no -remote) call like a panel built
  # without -parse_seqids; succeeds on -remote with one fixed rDNA hit.
  stub <- file.path(fake_bin, "blastn")
  writeLines(c(
    "#!/usr/bin/env bash",
    "case \" $* \" in",
    "  *\" -remote \"*)",
    "    printf 'REMOTE1\\tDanio rerio 18S ribosomal RNA gene\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'",
    "    exit 0",
    "    ;;",
    "  *)",
    "    echo 'BLAST Database error: No alias or index file found [TEST STUB]' >&2",
    "    exit 1",
    "    ;;",
    "esac"
  ), stub)
  Sys.chmod(stub, "0755")

  # a whole panel WITH a manifest beside it, so a unit wrongly credited to it
  # would leave a recorded provenance string behind: that is what makes the
  # UNIT_PANEL reset in blast.nf's blastn-failure branch checkable as behaviour
  panel <- file.path(panel_dir, "ribopilot_rrna")
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", paste0(panel, ".", e))
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(panel_dir, "manifest.json"))

  db_path <- file.path(root, "test.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER,
    entrez_query TEXT, extra_opts TEXT, max_target_seqs INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',1,'','',50)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, ignore INT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',0,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,NULL)")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  DBI::dbDisconnect(con)

  writeLines(c(
    sprintf("include { blast_genbank } from '%s'", blast_nf),
    "workflow {",
    "    blast_genbank(Channel.of('S1'))",
    "}"
  ), file.path(wf_dir, "main.nf"))
  writeLines(c(
    sprintf("params.db = '%s'", db_path),
    sprintf("params.blastdb = '%s'", panel),
    "params.taxdump_dir = ''",
    sprintf("params.ribopilot_src = '%s'", src),
    "process {",
    "    resourceLimits = [cpus: 1]",
    "    withName: 'blast_genbank' {",
    "        errorStrategy = 'terminate'",
    "        maxRetries = 0",
    "    }",
    "}"
  ), file.path(wf_dir, "nextflow.config"))

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_bin, bin_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  withr::local_dir(wf_dir)
  status <- system2("nextflow", c("run", "main.nf", "-w", shQuote(work_dir)),
    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0,
    info = "the run must succeed by falling back to remote, not exit 1 on the local failure")

  con2 <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con2), add = TRUE)
  searched <- DBI::dbGetQuery(con2, "SELECT source FROM blast_searched WHERE ID='S1'")
  expect_equal(nrow(searched), 1L,
    info = "the unit must end up recorded as searched, not silently dropped")
  expect_equal(searched$source[1], "remote",
    info = "the source must be recorded truthfully as remote, not the originally configured local")
  # fails if blast.nf's blastn-failure branch stops clearing UNIT_PANEL
  expect_equal(RiboPilot:::blastdb_provenance(con2), "",
    info = "a unit that fell back to remote must not be credited with the local panel")
})

# Shared fixture for the two rp_dbarg.R end-to-end cases below: an unconfigured
# project (blastdb "") must fall back to remote in silence, while a configured
# panel that rp_dbarg.R cannot use must fall back to remote AND say so. Both
# stub blastn (a real -remote call needs network and a real install) and read
# .command.err from the real work dir, since nextflow does not stream a task's
# own stderr to its own console log.
mk_dbarg_project <- function(root, blastdb) {
  fake_bin <- file.path(root, "fakebin"); dir.create(fake_bin)
  stub <- file.path(fake_bin, "blastn")
  writeLines(c(
    "#!/usr/bin/env bash",
    "case \" $* \" in",
    "  *\" -remote \"*)",
    "    printf 'REMOTE1\\tDanio rerio 18S ribosomal RNA gene\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'",
    "    exit 0",
    "    ;;",
    "  *)",
    "    echo 'unexpected local blastn call [TEST STUB]' >&2",
    "    exit 1",
    "    ;;",
    "esac"
  ), stub)
  Sys.chmod(stub, "0755")

  db_path <- file.path(root, "test.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER,
    entrez_query TEXT, extra_opts TEXT, max_target_seqs INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',1,'','',50)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, ignore INT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',0,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,NULL)")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbDisconnect(con)

  src_root <- RiboPilot:::ribopilot_src_default()
  blast_nf <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(blast_nf)) blast_nf <- file.path(src_root, "inst", "nextflow", "modules", "blast.nf")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src_root, "inst", "nextflow", "bin")

  wf_dir <- file.path(root, "wf"); dir.create(wf_dir)
  writeLines(c(
    sprintf("include { blast_genbank } from '%s'", blast_nf),
    "workflow {",
    "    blast_genbank(Channel.of('S1'))",
    "}"
  ), file.path(wf_dir, "main.nf"))
  writeLines(c(
    sprintf("params.db = '%s'", db_path),
    sprintf("params.blastdb = '%s'", blastdb),
    "params.taxdump_dir = ''",
    sprintf("params.ribopilot_src = '%s'", src_root),
    "process {",
    "    resourceLimits = [cpus: 1]",
    "    withName: 'blast_genbank' {",
    "        errorStrategy = 'terminate'",
    "        maxRetries = 0",
    "    }",
    "}"
  ), file.path(wf_dir, "nextflow.config"))
  list(wf_dir = wf_dir, db_path = db_path, fake_bin = fake_bin, bin_dir = bin_dir)
}

run_dbarg_project <- function(root, proj) {
  work_dir <- file.path(root, "work")
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(proj$fake_bin, proj$bin_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE, after = FALSE)
  withr::local_dir(proj$wf_dir, .local_envir = parent.frame())
  status <- system2("nextflow", c("run", "main.nf", "-w", shQuote(work_dir)),
    stdout = FALSE, stderr = FALSE)
  task_err <- list.files(work_dir, pattern = "^\\.command\\.err$",
    recursive = TRUE, full.names = TRUE, all.files = TRUE)
  list(status = status, stderr = paste(unlist(lapply(task_err, readLines, warn = FALSE)), collapse = "\n"))
}

test_that("an unconfigured project falls back to remote without usage noise", {
  # empty params.blastdb still calls rp_dbarg.R (a later task needs that for
  # discovery), but rp_args() rejecting the empty argument must not leak
  # "Error: usage: ..." / "Execution halted" into the task log
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")

  root <- withr::local_tempdir()
  # pinned at an empty directory, or a panel installed on the host makes this
  # pass through the "local blastn failed, fall back to remote" path instead of
  # the "nothing discovered" path it exists to test
  empty_db_dir <- file.path(root, "no_panels"); dir.create(empty_db_dir)
  withr::local_envvar(RIBOPILOT_DB_DIR = empty_db_dir)
  proj <- mk_dbarg_project(root, "")
  res <- run_dbarg_project(root, proj)
  expect_equal(res$status, 0)
  expect_false(grepl("usage:", res$stderr, fixed = TRUE))
  expect_false(grepl("Execution halted", res$stderr, fixed = TRUE))
  expect_false(grepl("configured local panel", res$stderr, fixed = TRUE),
    info = "no panel configured is the normal case, not something to report")
  expect_false(grepl("unexpected local blastn call", res$stderr, fixed = TRUE),
    info = "nothing should have been discovered, so blastn must never be called locally")

  con <- DBI::dbConnect(RSQLite::SQLite(), proj$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  searched <- DBI::dbGetQuery(con, "SELECT source FROM blast_searched WHERE ID='S1'")
  expect_equal(searched$source[1], "remote")
})

test_that("a changed BLAST setting searches a unit again, an unchanged one does not", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  root <- withr::local_tempdir()
  empty_db_dir <- file.path(root, "no_panels"); dir.create(empty_db_dir)
  withr::local_envvar(RIBOPILOT_DB_DIR = empty_db_dir)
  proj <- mk_dbarg_project(root, "")
  calls <- file.path(root, "blastn_calls.txt")
  writeLines(c("#!/usr/bin/env bash", sprintf("echo x >> '%s'", calls),
    "printf 'REMOTE1\\tDanio rerio 18S ribosomal RNA gene\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'"),
    file.path(proj$fake_bin, "blastn"))
  n_calls <- function() if (file.exists(calls)) length(readLines(calls)) else 0L
  expect_equal(run_dbarg_project(root, proj)$status, 0)
  expect_equal(n_calls(), 1L)
  run_dbarg_project(root, proj)
  expect_equal(n_calls(), 1L)
  con <- DBI::dbConnect(RSQLite::SQLite(), proj$db_path)
  DBI::dbExecute(con, "UPDATE blast_opts SET max_target_seqs = 20")
  DBI::dbDisconnect(con)
  run_dbarg_project(root, proj)
  expect_equal(n_calls(), 2L)
})

test_that("a failed escalation to NCBI is tried again on the next run", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  root <- withr::local_tempdir()
  panel_dir <- file.path(root, "panel"); dir.create(panel_dir)
  panel <- file.path(panel_dir, "ribopilot_rrna")
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", paste0(panel, ".", e))
  proj <- mk_dbarg_project(root, panel)
  calls <- file.path(root, "blastn_calls.txt")
  # the panel names nothing and NCBI cannot be reached
  writeLines(c("#!/usr/bin/env bash", sprintf("echo \"$*\" >> '%s'", calls),
    "case \" $* \" in *\" -remote \"*) exit 1 ;; *) exit 0 ;; esac"),
    file.path(proj$fake_bin, "blastn"))
  remote_calls <- function() sum(grepl("-remote", readLines(calls), fixed = TRUE))
  expect_equal(run_dbarg_project(root, proj)$status, 0)
  expect_equal(remote_calls(), 1L)
  con <- DBI::dbConnect(RSQLite::SQLite(), proj$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  expect_equal(DBI::dbGetQuery(con, "SELECT query_md5 FROM blast_searched")$query_md5,
    "remote_fallback_failed")
  run_dbarg_project(root, proj)
  expect_equal(remote_calls(), 2L)
})

test_that("a configured panel that rp_dbarg.R cannot use is reported, not silently dropped", {
  # local_blast_db_arg()'s commonest failure (blast_db_exists() FALSE, e.g. the
  # panel path was never bind-mounted into the container) returns NULL with no
  # R warning of its own, so blast.nf itself must say something or the run
  # silently spends hours on a remote search nobody asked for
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")

  root <- withr::local_tempdir()
  missing_panel <- file.path(root, "not_mounted", "ribopilot_rrna")  # no .nin ever created
  proj <- mk_dbarg_project(root, missing_panel)
  res <- run_dbarg_project(root, proj)
  expect_equal(res$status, 0)
  expect_true(grepl("configured local panel", res$stderr, fixed = TRUE))
  expect_true(grepl(missing_panel, res$stderr, fixed = TRUE),
    info = "must name the panel that could not be used")

  con <- DBI::dbConnect(RSQLite::SQLite(), proj$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  searched <- DBI::dbGetQuery(con, "SELECT source FROM blast_searched WHERE ID='S1'")
  expect_equal(searched$source[1], "remote")
})

test_that("blast_genbank records provenance for a discovered, unconfigured panel", {
  # The Docker-image scenario: no params.blastdb, a panel found only via
  # RIBOPILOT_DB_DIR (discover_panel()). Provenance has to be written from
  # inside the task, since the host never learns which prefix was used.
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  blast_nf <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(blast_nf)) blast_nf <- file.path(src, "inst", "nextflow", "modules", "blast.nf")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src, "inst", "nextflow", "bin")

  root <- withr::local_tempdir()
  fake_bin <- file.path(root, "fakebin"); dir.create(fake_bin)
  db_dir <- file.path(root, "panels"); dir.create(db_dir)
  panel_dir <- file.path(db_dir, "2026.08.1"); dir.create(panel_dir)
  wf_dir <- file.path(root, "wf"); dir.create(wf_dir)
  work_dir <- file.path(root, "work")

  for (e in c("nin", "nsq", "nhr", "ndb")) {
    writeLines("x", file.path(panel_dir, paste0("ribopilot_rrna.", e)))
  }
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(panel_dir, "manifest.json"))

  # stub blastn: any local (-db, no -remote) call succeeds with one plausible
  # hit, so the search stays local -- this test has no network for a real -remote
  stub <- file.path(fake_bin, "blastn")
  writeLines(c(
    "#!/usr/bin/env bash",
    "case \" $* \" in",
    "  *\" -remote \"*)",
    "    echo 'discovered-panel test must not go remote [TEST STUB]' >&2",
    "    exit 1",
    "    ;;",
    "  *)",
    "    printf 'LOCAL1\\tSome nrRNA\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'",
    "    exit 0",
    "    ;;",
    "esac"
  ), stub)
  Sys.chmod(stub, "0755")

  db_path <- file.path(root, "test.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER,
    entrez_query TEXT, extra_opts TEXT, max_target_seqs INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',1,'','',50)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, ignore INT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',0,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,NULL)")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  DBI::dbDisconnect(con)

  writeLines(c(
    sprintf("include { blast_genbank } from '%s'", blast_nf),
    "workflow {",
    "    blast_genbank(Channel.of('S1'))",
    "}"
  ), file.path(wf_dir, "main.nf"))
  writeLines(c(
    sprintf("params.db = '%s'", db_path),
    "params.blastdb = ''",
    "params.taxdump_dir = ''",
    sprintf("params.ribopilot_src = '%s'", src),
    "process {",
    "    resourceLimits = [cpus: 1]",
    "    withName: 'blast_genbank' {",
    "        errorStrategy = 'terminate'",
    "        maxRetries = 0",
    "    }",
    "}"
  ), file.path(wf_dir, "nextflow.config"))

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_bin, bin_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  withr::local_envvar(RIBOPILOT_DB_DIR = db_dir)

  withr::local_dir(wf_dir)
  status <- system2("nextflow", c("run", "main.nf", "-w", shQuote(work_dir)),
    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0)

  con2 <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con2), add = TRUE)
  searched <- DBI::dbGetQuery(con2, "SELECT source, db FROM blast_searched WHERE ID='S1'")
  expect_equal(searched$source[1], "local",
    info = "the stub only answers non-remote calls; a wrong source means the panel was not used")
  expect_equal(searched$db[1], "ribopilot_rrna 2026.08.1")
  expect_equal(RiboPilot:::blastdb_provenance(con2), "ribopilot_rrna 2026.08.1")

  # discovery must never leak the container-side path into the host-config column
  cfg <- DBI::dbGetQuery(con2, "SELECT blastdb FROM reference_config WHERE rowid_key=1")$blastdb
  expect_true(length(cfg) == 0 || is.na(cfg[1]) || !nzchar(cfg[1]))
})

test_that("blast.nf warns when extra_opts sets -evalue for a local search", {
  # core_nt is ~10^12 bp and the local panel ~7x10^8 bp, so the shipped
  # default -evalue 1e-20 admits far weaker alignments locally. No units are
  # ignore=0, so the per-unit loop never runs and no blastn call is needed;
  # this only checks the warning that fires right after SRC/DBARG resolve.
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  blast_nf <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(blast_nf)) blast_nf <- file.path(src, "inst", "nextflow", "modules", "blast.nf")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src, "inst", "nextflow", "bin")

  run_once <- function(extra_opts) {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    panel_dir <- file.path(root, "panel"); dir.create(panel_dir)
    wf_dir <- file.path(root, "wf"); dir.create(wf_dir)
    work_dir <- file.path(root, "work")
    panel <- file.path(panel_dir, "ribopilot_rrna")
    file.create(paste0(panel, ".nin"))

    db_path <- file.path(root, "test.sqlite")
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
    DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER,
      entrez_query TEXT, extra_opts TEXT, max_target_seqs INTEGER)")
    DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
    DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',1,'',?,50)", params = list(extra_opts))
    # no assemblies rows: UNITS resolves empty, so the per-unit loop never runs
    DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, ignore INT)")
    DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
    DBI::dbDisconnect(con)

    writeLines(c(
      sprintf("include { blast_genbank } from '%s'", blast_nf),
      "workflow {",
      "    blast_genbank(Channel.of('S1'))",
      "}"
    ), file.path(wf_dir, "main.nf"))
    writeLines(c(
      sprintf("params.db = '%s'", db_path),
      sprintf("params.blastdb = '%s'", panel),
      "params.taxdump_dir = ''",
      sprintf("params.ribopilot_src = '%s'", src),
      "process {",
      "    resourceLimits = [cpus: 1]",
      "    withName: 'blast_genbank' {",
      "        errorStrategy = 'terminate'",
      "        maxRetries = 0",
      "    }",
      "}"
    ), file.path(wf_dir, "nextflow.config"))

    withr::local_dir(wf_dir, .local_envir = parent.frame())
    status <- system2("nextflow", c("run", "main.nf", "-w", shQuote(work_dir)),
      stdout = FALSE, stderr = FALSE)
    # nextflow does not stream a task's own stderr to its own console log;
    # read it from the task work dir instead.
    task_err <- list.files(work_dir, pattern = "^\\.command\\.err$",
      recursive = TRUE, full.names = TRUE, all.files = TRUE)
    task_stderr <- paste(unlist(lapply(task_err, readLines, warn = FALSE)), collapse = "\n")
    list(status = status, stderr = task_stderr)
  }

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  with_evalue <- run_once("-evalue 1e-20")
  expect_equal(with_evalue$status, 0)
  expect_match(with_evalue$stderr, "WARNING:.*-evalue 1e-20",
    info = "must name the value and warn the local search space differs from core_nt")

  without_evalue <- run_once("")
  expect_equal(without_evalue$status, 0)
  expect_false(grepl("-evalue", without_evalue$stderr, fixed = TRUE),
    info = "no -evalue in extra_opts means no warning")
})

test_that("a local failure on one unit does not disable local search for the rest of the sample", {
  # blast.nf's per-unit loop is a plain shell `for`, not a subshell, so a
  # mutation of SRC inside the loop body used to leak into every later
  # iteration. Two units, only the first fails locally: the second must
  # still search locally, not silently fall back because of the first.
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  skip_if(Sys.which("sqlite3") == "", "sqlite3 not installed")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!file.exists(file.path(src, "DESCRIPTION")), "not a source checkout")
  blast_nf <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(blast_nf)) blast_nf <- file.path(src, "inst", "nextflow", "modules", "blast.nf")
  bin_dir <- system.file("nextflow", "bin", package = "RiboPilot")
  if (!nzchar(bin_dir)) bin_dir <- file.path(src, "inst", "nextflow", "bin")

  root <- withr::local_tempdir()
  fake_bin <- file.path(root, "fakebin"); dir.create(fake_bin)
  panel_dir <- file.path(root, "panel"); dir.create(panel_dir)
  wf_dir <- file.path(root, "wf"); dir.create(wf_dir)
  work_dir <- file.path(root, "work")

  # stub blastn: fails the LOCAL call only for scaffold 1 (header "S1.1.1"),
  # succeeds locally for scaffold 2, and always succeeds on -remote.
  stub <- file.path(fake_bin, "blastn")
  writeLines(c(
    "#!/usr/bin/env bash",
    "case \" $* \" in",
    "  *\" -remote \"*)",
    "    printf 'REMOTE1\\tDanio rerio 18S ribosomal RNA gene\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'",
    "    exit 0",
    "    ;;",
    "esac",
    "Q=\"\"; prev=\"\"",
    "for a in \"$@\"; do",
    "  if [ \"$prev\" = \"-query\" ]; then Q=\"$a\"; fi",
    "  prev=\"$a\"",
    "done",
    "HDR=$(head -n1 \"$Q\" 2>/dev/null || echo '')",
    "case \"$HDR\" in",
    "  *'.1.1'*)",
    "    echo 'BLAST Database error: No alias or index file found [TEST STUB]' >&2",
    "    exit 1",
    "    ;;",
    "  *)",
    "    printf 'LOCAL1\\tDanio rerio 18S ribosomal RNA gene\\t99.0\\t80\\t1e-90\\t7955\\t1500\\t1600\\t1\\t1600\\t1\\t1600\\t6997\\tq1\\n'",
    "    exit 0",
    "    ;;",
    "esac"
  ), stub)
  Sys.chmod(stub, "0755")

  panel <- file.path(panel_dir, "ribopilot_rrna")
  file.create(paste0(panel, ".nin"))

  db_path <- file.path(root, "test.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER,
    entrez_query TEXT, extra_opts TEXT, max_target_seqs INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',1,'','',50)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT, ignore INT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL, blast_source TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',0,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,2,
    'TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA',0,NULL,NULL,NULL,NULL,NULL,NULL)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, taxid INT, lineage TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',7955,NULL,NULL)")
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, bitscore REAL, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    evalue REAL, bitscore REAL, is_nuclear_rrna INT, time_stamp INT)")
  DBI::dbDisconnect(con)

  writeLines(c(
    sprintf("include { blast_genbank } from '%s'", blast_nf),
    "workflow {",
    "    blast_genbank(Channel.of('S1'))",
    "}"
  ), file.path(wf_dir, "main.nf"))
  writeLines(c(
    sprintf("params.db = '%s'", db_path),
    sprintf("params.blastdb = '%s'", panel),
    "params.taxdump_dir = ''",
    sprintf("params.ribopilot_src = '%s'", src),
    "process {",
    "    resourceLimits = [cpus: 1]",
    "    withName: 'blast_genbank' {",
    "        errorStrategy = 'terminate'",
    "        maxRetries = 0",
    "    }",
    "}"
  ), file.path(wf_dir, "nextflow.config"))

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(fake_bin, bin_dir, old_path, sep = .Platform$path.sep))
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)

  withr::local_dir(wf_dir)
  status <- system2("nextflow", c("run", "main.nf", "-w", shQuote(work_dir)),
    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0)

  con2 <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con2), add = TRUE)
  searched <- DBI::dbGetQuery(con2,
    "SELECT scaffold, source FROM blast_searched WHERE ID='S1' ORDER BY scaffold")
  expect_equal(nrow(searched), 2L)
  expect_equal(searched$source[searched$scaffold == 1], "remote",
    info = "unit 1's local search failed, so it must fall back to remote")
  expect_equal(searched$source[searched$scaffold == 2], "local",
    info = "unit 2's local search would have succeeded; unit 1's failure must not leak into it")
})

# User-supplied assemblies enter WF1 where an assembly would, rather than as a
# second workflow. See dev/design/20_USER_ASSEMBLIES.md.

test_that("query_samples carries the columns the assembly branch reads", {
  start <- grep("^process query_samples \\{", pipeline_r_nf)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]
  hdr <- grep("printf 'ID", body, value = TRUE)

  for (col in c("Assembly", "find_rrna", "rrna_tool", "rrna_kingdom")) {
    expect_match(hdr, col)
    expect_true(any(grepl(paste0("AS ", col), body)))
  }
  expect_true(any(grepl("JOIN samples s", body)))
})

samples_sql <- function() {
  start <- grep("^process query_samples \\{", pipeline_r_nf)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]
  sql <- body[grep("SELECT p.ID", body):grep(";\" >> samplesheet.tsv", body, fixed = TRUE)]
  sql <- sub(";\" >> samplesheet.tsv", "", paste(sql, collapse = "\n"), fixed = TRUE)
  gsub("\\$\\{params\\.db\\}", "", sql)
}

test_that("query_samples skips samples whose taxonomy is not locked", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Assembly TEXT, taxonomy_lock INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, assemble_lock INT,
    assemble_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE preprocess(ID TEXT, R1 TEXT, R2 TEXT, pre_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE pre_opts(pre_opts TEXT, cpus INT, memory INT, fastp TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, cpus INT, memory INT,
    getOrganelle TEXT, max_paths INT, seeds_db TEXT, labels_db TEXT, assembler TEXT, find_rrna INT,
    bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate_opts(annotate_opts TEXT, rrna_tool TEXT, rrna_kingdom TEXT,
    barrnap_opts TEXT, pybarrnap_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO pre_opts VALUES('default',4,8,'')")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default',4,16,'--flag',10,NULL,NULL,'getorganelle',1,NULL,NULL)")
  DBI::dbExecute(con, "INSERT INTO annotate_opts VALUES('default','pybarrnap','euk','','')")
  for (id in c("LOCKED", "UNLOCKED")) {
    DBI::dbExecute(con, "INSERT INTO samples VALUES(?,NULL,?)",
      params = list(id, if (id == "LOCKED") 1L else 0L))
    DBI::dbExecute(con, "INSERT INTO preprocess VALUES(?,'a.fq','b.fq','default')", params = list(id))
    DBI::dbExecute(con, "INSERT INTO assemble VALUES(?,1,0,'default')", params = list(id))
  }
  got <- DBI::dbGetQuery(con, samples_sql())$ID
  expect_equal(got, "LOCKED")
})

test_that("main.nf routes supplied assemblies away from preprocess and assembly", {
  br <- grep("sample_rows = rows\\.branch", main_nf)
  expect_length(br, 1)
  block <- paste(main_nf[br:(br + 4)], collapse = "\n")
  expect_match(block, "user: r\\.Assembly")
  expect_match(block, "asis: r\\.Assembly")
  # `reads: true`, not an equality test, for the same reason the assembler
  # branch uses it: branch silently drops an item matching no condition.
  expect_match(block, "reads: true")

  # preprocess only for rows that actually name reads
  expect_true(any(grepl("filter \\{ r -> r\\.R1\\?\\.trim\\(\\) && r\\.R2\\?\\.trim\\(\\) \\}", main_nf)))
  # ...and assembly only for rows assembling from reads, so a supplied assembly
  # that also names reads (for coverage) is not assembled as well
  expect_true(any(grepl("read_ids = sample_rows\\.reads", main_nf)))
  expect_true(any(grepl("\\.join\\(read_ids\\)", main_nf)))
})

test_that("find_rrna emits the assembler's fasta contract", {
  nf <- readLines(system.file("nextflow", "modules", "find_rrna.nf", package = "RiboPilot"))
  # published to out/<id>/assemble, named what score() expects
  expect_true(any(grepl("\\$\\{params\\.publishDir\\}/\\$\\{id\\}/assemble", nf)))
  expect_true(any(grepl("\\$\\{id\\}_assembly_\\*\\.fasta", nf)))
  # the sentinel name load_assemblies() skips, for a sample with nothing found
  expect_true(any(grepl("\\$\\{id\\}_assembly_0\\.fasta", nf)))
  # chunked so a draft genome is never read whole
  expect_true(any(grepl("splitFasta\\(by: params\\.find_rrna_chunk_size", nf)))
  # the default must live in THIS module: includes resolve before main.nf's body
  # runs, so a default set there reaches the splitFasta above as null
  expect_true(any(grepl("^params\\.find_rrna_chunk_size = [0-9]+", nf)))
  expect_false(any(grepl("params\\.find_rrna_chunk_size", main_nf)))
})

test_that("find_rrna_pick stages every scan chunk under its own name", {
  nf <- readLines(system.file("nextflow", "modules", "find_rrna.nf", package = "RiboPilot"))
  expect_true(any(grepl("tuple val(id), path(assembly), path(gffs, stageAs: 'chunk_*.gff')", nf, fixed = TRUE)))
})

test_that("query_samples asks GetOrganelle for one path over max_paths, and sets no cap for 0", {
  start <- grep("^process query_samples \\{", pipeline_r_nf)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  line <- grep("--max-paths-num", pipeline_r_nf[start:end], value = TRUE)
  expect_length(line, 1)
  expr <- sub("\\s+AS getOrganelle,\\s*$", "", trimws(line))
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE assemble_opts (getOrganelle TEXT, max_paths INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES ('-F anonym', 10), ('-F anonym', 0), ('-F anonym', NULL)")
  got <- DBI::dbGetQuery(con, paste("SELECT", expr, "AS go FROM assemble_opts ao"))$go
  expect_equal(got, c("-F anonym --max-paths-num 11", "-F anonym", "-F anonym"))
})

# ---- rawDir is only required by samples that actually have reads ----

mk_rawdir_project <- function(dir, raw = "NA", r1 = NA_character_) {
  writeLines(c("params {", paste0("  rawDir = '", raw, "'"), "}"), file.path(dir, ".config"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE preprocess (ID TEXT, R1 TEXT, R2 TEXT)")
  if (!is.na(r1)) {
    DBI::dbExecute(con, "INSERT INTO preprocess VALUES ('S1', ?, 'b.fq')", params = list(r1))
  }
  DBI::dbDisconnect(con)
  dir
}

test_that("an assembly-only project runs without a rawDir", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_rawdir_project(proj)
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_false("--rawDir" %in% a)
})

test_that("a project with reads still refuses to run without a rawDir", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_rawdir_project(proj, raw = "NA", r1 = "a.fq")
  expect_error(RiboPilot:::nf_command(proj, resume = FALSE), "No rawDir set")
})

test_that("a project with reads and a rawDir passes it through", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  mk_rawdir_project(proj, raw = "/data/raw", r1 = "a.fq")
  a <- RiboPilot:::nf_command(proj, resume = FALSE)
  expect_equal(a[which(a == "--rawDir") + 1], "/data/raw")
})

test_that("a reference folder with a space stays one mount argument", {
  skip_on_os("windows")
  d <- "/x/Library/Application Support/ref_dbs"
  for (engine in c("docker", "singularity")) {
    v <- RiboPilot:::container_bind_opts(RiboPilot:::container_engine_block(engine), d)[2]
    words <- system2("bash", c("-c", shQuote(paste("printf '%s\\n'", v))), stdout = TRUE)
    expect_equal(words[2], if (engine == "docker") paste0(d, ":", d) else d, info = engine)
  }
})

test_that("samples below minDepth are marked failed with their read count", {
  pre_nf <- readLines(system.file("nextflow", "modules", "preprocess.nf", package = "RiboPilot"))
  expect_true(any(grepl("include \\{ preprocess; below_min_depth \\}", main_nf)))
  start <- grep("^\\s*below_min_depth\\(", main_nf)
  expect_length(start, 1)
  block <- paste(main_nf[start:(start + 4)], collapse = "\n")
  expect_match(block, "(after as Integer) < (params.minDepth as Integer)", fixed = TRUE)
  # only samples assembling from reads; a supplied assembly with shallow reads is fine
  expect_match(block, ".join(read_ids)", fixed = TRUE)
  p <- grep("^process below_min_depth \\{", pre_nf)
  expect_length(p, 1)
  body <- paste(pre_nf[p:length(pre_nf)], collapse = "\n")
  expect_match(body, "assemble_switch=3", fixed = TRUE)
  expect_match(body, "below min_depth: ${after} reads after trimming", fixed = TRUE)
  expect_match(body, "cache false", fixed = TRUE)
  expect_false(grepl("'terminate'", body, fixed = TRUE))
})

test_that("score reads the assembly files staged into its task", {
  start <- grep("^process score \\{", pipeline_r_nf)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- paste(pipeline_r_nf[start:end], collapse = "\n")
  expect_match(body, "rp_score.R '${params.db}' '${id}' . '${params.taxdump_dir}'", fixed = TRUE)
  expect_false(grepl("publishDir", body, fixed = TRUE))
  asm <- readLines(system.file("nextflow", "modules", "assemble.nf", package = "RiboPilot"))
  expect_true(any(grepl('tuple val(id), path("{${id}_assembly_*.fasta,*.gfa}"), emit: paths',
    asm, fixed = TRUE)))
})

test_that("a fallback assembly's databases and tool options are task inputs", {
  fb <- readLines(system.file("nextflow", "modules", "fallback.nf", package = "RiboPilot"))
  expect_true(any(grepl(paste("tuple val(id), path(reads), val(nthreads), val(mem), val(method),",
    "path(seeds, stageAs: 'seeds/*'), path(labels, stageAs: 'labels/*'), val(seeds_url), val(labels_url),",
    "val(bowtie2), val(spades)"), fb, fixed = TRUE)))
  expect_false(any(grepl("params.db", fb, fixed = TRUE)))
  txt <- paste(main_nf, collapse = "\n")
  expect_match(txt, "def s = method == 'bait' ? seeds : ''", fixed = TRUE)
  src <- readLines(system.file("nextflow", "bin", "rp_fallback.R", package = "RiboPilot"))
  expect_false(any(grepl("open_project_db|dbGetQuery", src)))
})

test_that("a database URL reaches the fallback as text, so a dead one fails only its sample", {
  txt <- paste(main_nf, collapse = "\n")
  expect_match(txt, "def stage_db(p) { p && !is_url(p) ? file(p) : [] }", fixed = TRUE)
  expect_match(txt, "stage_db(s), stage_db(labels), url_db(s), url_db(labels), bt2, spades", fixed = TRUE)
  fb <- paste(readLines(system.file("nextflow", "modules", "fallback.nf", package = "RiboPilot")),
    collapse = "\n")
  expect_match(fb, "'${seeds}' '${labels}' '${bowtie2}' '${spades}' '${seeds_url}' '${labels_url}'",
    fixed = TRUE)
  src <- paste(readLines(system.file("nextflow", "bin", "rp_fallback.R", package = "RiboPilot")),
    collapse = "\n")
  expect_match(src, "reference = db(7), labels_db = db(8)", fixed = TRUE)
})

test_that("seeds and labels with the same file name do not collide", {
  for (f in c("fallback.nf", "assemble.nf")) {
    nf <- paste(readLines(system.file("nextflow", "modules", f, package = "RiboPilot")), collapse = "\n")
    expect_match(nf, "path(seeds, stageAs: 'seeds/*'), path(labels, stageAs: 'labels/*')",
      fixed = TRUE, info = f)
  }
})

test_that("a scaffold never searched is sent to BLAST and the screen when it is annotated", {
  con <- mk_units_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, "PUTBACK")
  add_unit(con, "SEARCHED")
  DBI::dbExecute(con, "INSERT INTO blast_searched VALUES('SEARCHED',1,1)")
  got <- DBI::dbGetQuery(con, paste(annotate_units_sql(), "ORDER BY an.ID"))
  expect_equal(got$searched, c(0L, 1L))
  code <- paste(main_nf[!grepl("^\\s*//", main_nf)], collapse = "\n")
  expect_match(code, "filter { r -> r.searched == '0' }", fixed = TRUE)
  expect_match(code, "blast_genbank(score.out.mix(unsearched))", fixed = TRUE)
})

test_that("find_rrna_pick and stage_assembly always run, so their database writes do", {
  nf <- readLines(system.file("nextflow", "modules", "find_rrna.nf", package = "RiboPilot"))
  for (p in c("find_rrna_pick", "stage_assembly")) {
    s <- grep(paste0("^process ", p, " \\{"), nf)
    body <- nf[s:(s - 1 + which(nf[s:length(nf)] == "}")[1])]
    expect_true(any(grepl("^\\s*cache false$", body)), info = p)
  }
})

test_that("the pick step runs even when every scan part failed, and scan parts retry", {
  txt <- paste(readLines(system.file("nextflow", "modules", "find_rrna.nf", package = "RiboPilot")),
    collapse = "\n")
  expect_match(txt, ".join(gffs_by_id, remainder: true)", fixed = TRUE)
  expect_match(txt, "tuple(id, assembly, gffs ?: [])", fixed = TRUE)
  expect_match(txt, "${params.find_rrna_chunk_size} ${gffs}", fixed = TRUE)
  cfg <- paste(readLines(system.file("nextflow", "nextflow.config", package = "RiboPilot")),
    collapse = " ")
  expect_match(cfg, "withName: 'assemble|fallback_assemble|preprocess|annotate_tools|find_rrna_scan|recheck_tools'",
    fixed = TRUE)
})

test_that("user rRNA caller flags reach both the per-unit and the scan commands", {
  nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))
  expect_true(any(grepl("barrnap_opts ELSE ao.pybarrnap_opts", nf, fixed = TRUE)))
  expect_true(any(grepl('"\\$RRNAO" ${task.cpus}', nf, fixed = TRUE)))
  expect_true(any(grepl("AS rrna_opts", nf, fixed = TRUE)))
  scan <- readLines(system.file("nextflow", "modules", "find_rrna.nf", package = "RiboPilot"))
  expect_true(any(grepl("--quiet ${opts} ${chunk}", scan, fixed = TRUE)))
  expect_true(any(grepl("-q ${opts} ${chunk}", scan, fixed = TRUE)))
  main <- readLines(system.file("nextflow", "main.nf", package = "RiboPilot"))
  expect_true(any(grepl("r.rrna_kingdom, r.rrna_opts)", main, fixed = TRUE)))
})

# ---- the recheck workflow: main_recheck.nf and modules/recheck.nf ----

recheck_nf <- readLines(system.file("nextflow", "modules", "recheck.nf", package = "RiboPilot"))
main_recheck_nf <- readLines(system.file("nextflow", "main_recheck.nf", package = "RiboPilot"))

recheck_body <- function(lines, proc) {
  start <- grep(paste0("^process ", proc, " \\{"), lines)
  expect_length(start, 1)
  lines[start:(start - 1 + which(lines[start:length(lines)] == "}")[1])]
}

test_that("main_recheck.nf feeds query_recheck through recheck_tools into recheck", {
  expect_true(any(main_recheck_nf ==
    "include { query_recheck; recheck_tools; recheck } from './modules/recheck.nf'"))
  code <- paste(main_recheck_nf[!grepl("^\\s*//", main_recheck_nf)], collapse = "\n")
  expect_match(code, "query_recheck().splitCsv(header: true, sep: '\\t')", fixed = TRUE)
  expect_match(code, "tuple(r.ID, r.path, r.scaffold, r.rev, r.cpus, r.mem, r.mode)", fixed = TRUE)
  expect_match(code, "recheck_tools(units)", fixed = TRUE)
  expect_match(code, "recheck(recheck_tools.out)", fixed = TRUE)
  # every field the workflow reads is a column query_recheck writes
  hdr <- grep("printf 'ID", recheck_body(recheck_nf, "query_recheck"), value = TRUE)
  expect_length(hdr, 1)
  cols <- strsplit(sub(".*printf '", "", sub("\\\\n' > .*", "", hdr)), "\\\\t")[[1]]
  expect_equal(cols, c("ID", "path", "scaffold", "rev", "cpus", "mem", "mode"))
  used <- sub("^r\\.", "", regmatches(code, gregexpr("\\br\\.[A-Za-z_]+", code))[[1]])
  expect_true(all(used %in% cols))
  # uncached, and the per-contig processes are tagged by contig
  for (p in c("query_recheck", "recheck_tools", "recheck")) {
    expect_true(any(grepl("^\\s*cache false$", recheck_body(recheck_nf, p))), info = p)
  }
  for (p in c("recheck_tools", "recheck")) {
    body <- recheck_body(recheck_nf, p)
    expect_true(any(grepl('tag "${id}.${upath}.${uscaf}"', body, fixed = TRUE)), info = p)
    expect_true(any(grepl("errorStrategy 'ignore'", body, fixed = TRUE)), info = p)
  }
  expect_true(any(grepl("executor 'local'", recheck_body(recheck_nf, "query_recheck"), fixed = TRUE)))
  expect_true(all(c("params.blastdb = ''", "params.taxdump_dir = ''") %in% recheck_nf))
})

test_that("query_recheck selects RECHECK_QUEUE_SQL word for word and marks it running", {
  body <- paste(recheck_body(recheck_nf, "query_recheck"), collapse = "\n")
  sql <- regmatches(body, regexpr("CREATE TEMP TABLE q AS\\s*\\K[^;]*", body, perl = TRUE))
  norm <- function(x) gsub("\\s+", " ", trimws(x))
  expect_equal(norm(sql), norm(RiboPilot:::RECHECK_QUEUE_SQL))
  # it travels inside a Nextflow string and a double-quoted shell argument
  expect_false(any(strsplit(RiboPilot:::RECHECK_QUEUE_SQL, "")[[1]] %in% c("$", "\"", "\\", ";")))
  expect_match(body, "BEGIN IMMEDIATE;", fixed = TRUE)
  expect_match(body, "UPDATE review_units SET check_switch = 4", fixed = TRUE)
  expect_match(body, "COMMIT;\" >> recheck_units.tsv", fixed = TRUE)
})

test_that("recheck_tools reads the tool options exactly as annotate_tools does", {
  opts <- function(lines, proc) {
    body <- recheck_body(lines, proc)
    from <- grep("read TOOL KINGDOM ITSXT ITSXO <<<", body, fixed = TRUE)
    to <- grep("RRNAO=\"", body, fixed = TRUE)
    expect_length(from, 1)
    expect_length(to, 1)
    body[from:(to - 1 + grep("\")\"$", body[to:length(body)])[1])]
  }
  got <- opts(recheck_nf, "recheck_tools")
  expect_gt(length(got), 5)
  expect_identical(got, opts(pipeline_r_nf, "annotate_tools"))
})

test_that("recheck_tools runs the tools only on the queued rev and fails only on out-of-memory exits", {
  body <- recheck_body(recheck_nf, "recheck_tools")
  txt <- paste(body, collapse = "\n")
  expect_match(txt, "'>'||ID||'.'||path||'.'||scaffold||char(10)||sequence FROM review_units", fixed = TRUE)
  expect_match(txt, "AND rev=${rev}", fixed = TRUE)
  expect_match(txt, "check_switch IN (1,4) AND sequence IS NOT NULL", fixed = TRUE)
  expect_match(txt, "if [ -s unit.fasta ]; then", fixed = TRUE)
  call <- grep("rp_rrna_tools.sh ", body, value = TRUE, fixed = TRUE)
  expect_length(call, 1)
  expect_equal(trimws(call), paste('rp_rrna_tools.sh unit.fasta "\\$TOOL" "\\$KINGDOM" "\\$ITSXT"',
    '"\\$ITSXO" "\\$RRNAO" ${task.cpus} || rc=\\$?'))
  codes <- grep('case "\\$rc" in', body, value = TRUE, fixed = TRUE)
  expect_length(codes, 1)
  codes <- as.integer(strsplit(sub(".* in ([0-9|]+)\\).*", "\\1", codes), "|", fixed = TRUE)[[1]])
  expect_setequal(codes, RiboPilot:::NF_RETRY_CODES)
  expect_match(txt, "path('unit.itsx_check*'), optional: true", fixed = TRUE)
})

test_that("the recheck process resolves the panel and labels as the chimera process does", {
  block <- function(lines, proc) {
    body <- recheck_body(lines, proc)
    from <- grep('DBARG=""', body, fixed = TRUE)
    to <- grep("2>/dev/null || echo \"\")", body, fixed = TRUE)
    expect_length(from, 1)
    body[from:to[length(to)]]
  }
  expect_identical(block(recheck_nf, "recheck"), block(pipeline_r_nf, "chimera"))
  call <- gsub("\\s+", " ", paste(recheck_body(recheck_nf, "recheck"), collapse = " "))
  expect_match(call, paste("RP_THREADS=${task.cpus} RIBOPILOT_SRC='${params.ribopilot_src}' rp_recheck.R \\",
    "'${params.db}' '${id}' ${upath} ${uscaf} ${rev} ${mode} . \"\\$DBARG\" '${params.taxdump_dir}' \"\\$LABELS\""),
    fixed = TRUE)
})

test_that("nf_recheck_command runs main_recheck.nf in work_recheck, without -resume", {
  skip_if_not_installed("withr")
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir()
  refs <- withr::local_tempdir()
  mk_ref_config_project(proj, refs, blastdb = "/opt/rp/ribopilot_rrna")
  writeLines(c("params {", "  rawDir = '/data/raw'", "  minDepth = 5",
    "  ncbi_api_key = 'abc123'", "}"), file.path(proj, ".config"))
  a <- RiboPilot:::nf_recheck_command(proj)
  val <- function(flag) a[which(a == flag) + 1]
  expect_equal(a[1:2], c("run", system.file("nextflow", "main_recheck.nf", package = "RiboPilot")))
  expect_true(file.exists(a[2]))
  expect_equal(val("-c"), file.path(proj, ".config"))
  expect_equal(val("-work-dir"), file.path(proj, "work_recheck"))
  expect_equal(val("--db"), file.path(proj, ".sqlite"))
  expect_equal(val("--projectDir"), proj)
  expect_equal(val("--rawDir"), "/data/raw")
  expect_equal(val("--taxdump_dir"), file.path(proj, ".taxonomy"))
  expect_equal(val("--blastdb"), "/opt/rp/ribopilot_rrna")
  expect_equal(val("--ribopilot_version"), as.character(utils::packageVersion("RiboPilot")))
  for (f in c("-resume", "--minDepth", "--ref_dir", "--ncbi_api_key")) expect_false(f %in% a, info = f)
  expect_true(all(nzchar(a)))
  # the main pipeline's command for the same project does carry them
  expect_true(all(c("--minDepth", "--ref_dir", "--ncbi_api_key") %in%
    RiboPilot:::nf_command(proj, resume = FALSE)))
})

test_that("nf_recheck_command leaves out a blank rawDir and blastdb instead of refusing", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  # reads but no rawDir: the main pipeline refuses, the recheck reads no reads
  mk_rawdir_project(proj, raw = "NA", r1 = "a.fq")
  expect_error(RiboPilot:::nf_command(proj, resume = FALSE), "No rawDir set")
  a <- RiboPilot:::nf_recheck_command(proj)
  expect_false("--rawDir" %in% a)
  expect_false("--blastdb" %in% a)
  expect_true(all(nzchar(a)))
})

test_that("nf_recheck_command mounts the reference folders as nf_command does", {
  skip_if_not_installed("withr")
  withr::local_envvar(RIBOPILOT_REF_DIR = tempfile("no-ref-dir-"))
  proj <- withr::local_tempdir(); refs <- withr::local_tempdir()
  dir.create(file.path(refs, "seeds")); dir.create(file.path(refs, "labels"))
  s <- file.path(refs, "seeds", "x_seeds.fasta"); l <- file.path(refs, "labels", "x_labels.fasta")
  writeLines(">18S", s); writeLines(">18S", l)
  mk_mount_project(proj, s, l)
  a <- RiboPilot:::nf_recheck_command(proj)
  b <- RiboPilot:::nf_command(proj, resume = FALSE)
  opt <- "-process.containerOptions"
  expect_length(which(a == opt), 1L)
  expect_equal(a[which(a == opt) + 1], b[which(b == opt) + 1])
  expect_true(grepl(file.path(refs, "seeds"), a[which(a == opt) + 1], fixed = TRUE))
})

test_that("building the recheck command records no panel provenance", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  writeLines(c("params {", "  rawDir = '/data/raw'", "}"), file.path(d, ".config"))
  panel_dir <- withr::local_tempdir()
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}',
    file.path(panel_dir, "manifest.json"))
  prefix <- file.path(panel_dir, "ribopilot_rrna")
  file.create(paste0(prefix, ".nin"))
  mk_ref_config_project(d, NA_character_, blastdb = prefix)
  expect_true("--blastdb" %in% RiboPilot:::nf_recheck_command(d))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::blastdb_provenance(con), "")
})
