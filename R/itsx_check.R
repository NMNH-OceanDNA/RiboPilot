# Final pre-export ITSx check, run the way NCBI's Submission Portal runs it:
# stock ITSx, default parameters, no organism-group restriction. The pipeline's
# own ITSx pass (annotation step) uses the project's itsx_taxon / itsx_opts and
# happens before curation, so it cannot stand in for what the Portal will see.
#
# Portal review flags chimeric sequences, gaps between ITS regions, and
# potential misassemblies. Each is mapped to what stock ITSx reports:
#   chimeric      "Chimeric!" (5.8S overlaps SSU/LSU) or domains out of order
#   gaps          an N-run inside the ITSx-detected ITS span
#   misassembly   missing/partial 5.8S, "ITS region too long!", missing ITS1/2,
#                 or ITS boundaries far from RiboPilot's own annotation

# max distance (bp) between an ITSx and a RiboPilot ITS boundary before it is flagged
ITSX_CHECK_SHIFT_BP <- 25L

# ITSx origin letter -> group name
ITSX_ORIGINS <- c(
  A = "Alveolata", B = "Bryophyta", C = "Bacillariophyta", D = "Amoebozoa",
  E = "Euglenozoa", F = "Fungi", G = "Chlorophyta", H = "Rhodophyta",
  I = "Phaeophyceae", L = "Marchantiophyta", M = "Metazoa", O = "Oomycota",
  P = "Haptophyceae", Q = "Raphidophyceae", R = "Rhizaria", S = "Synurophyceae",
  T = "Tracheophyta", U = "Eustigmatophyceae", X = "Apusozoa", Y = "Parabasalia")

# lineage name -> ITSx origin letters that are consistent with it
ITSX_LINEAGE_ORIGINS <- list(
  Metazoa = "M", Fungi = "F", Streptophyta = c("T", "B", "L"),
  Chlorophyta = "G", Rhodophyta = "H", Alveolata = "A", Amoebozoa = "D",
  Euglenozoa = "E", Oomycota = "O", Rhizaria = "R", Bacillariophyta = "C",
  Phaeophyceae = "I", Haptophyta = "P", Apusozoa = "X", Parabasalia = "Y")

#' How to invoke ITSx for a project: through its container, or from PATH
#'
#' The pipeline's tools live in the project's container (docker / singularity
#' per `.config`), so the check runs ITSx there. Only a project on the conda
#' executor (no container enabled) falls back to an ITSx on PATH.
#' @param project Project directory (holds `.config`).
#' @param workdir Directory the container must see: the check's input FASTA
#'   and output land here. Bind-mounted at the same path.
#' @param which Lookup for executables (injectable for tests).
#' @return list(bin, prefix, message): `bin` is the executable to call and
#'   `prefix` the arguments that go before ITSx's own; `bin` is "" and
#'   `message` says why when the check cannot run.
#' @noRd
itsx_check_command <- function(project, workdir, which = Sys.which) {
  cfg <- file.path(project, ".config")
  lines <- if (file.exists(cfg)) readLines(cfg, warn = FALSE) else character(0)
  engine <- config_container_engine(lines)
  none <- function(msg) list(bin = "", prefix = character(0), message = msg)
  if (is.na(engine)) {
    bin <- unname(which("ITSx"))
    if (!nzchar(bin)) return(none("ITSx not on PATH (no container enabled in .config)"))
    return(list(bin = bin, prefix = character(0), message = ""))
  }
  img <- config_container_image(lines)
  if (!nzchar(img)) return(none("no process.container image in .config"))
  wd <- normalizePath(workdir, mustWork = FALSE)
  if (identical(engine, "docker")) {
    bin <- unname(which("docker"))
    if (!nzchar(bin)) return(none("docker not on PATH"))
    ids <- vapply(c("-u", "-g"), function(f) {
      tryCatch(trimws(system2("id", f, stdout = TRUE, stderr = FALSE))[1], error = function(e) "")
    }, character(1))
    user <- if (all(nzchar(ids))) c("-u", paste(ids, collapse = ":")) else character(0)
    return(list(bin = bin,
      prefix = c("run", "--rm", user, "-v", paste0(wd, ":", wd), "-w", wd, img, "ITSx"),
      message = ""))
  }
  w <- which(c("singularity", "apptainer"))
  w <- w[nzchar(w)]
  if (length(w) < 1) return(none("singularity/apptainer not on PATH"))
  list(bin = unname(w[1]),
    prefix = c("exec", "--bind", wd, "--pwd", wd, img, "ITSx"),
    message = "")
}

#' Run stock ITSx on a FASTA
#'
#' Only output-shaping options are passed (`--preserve`, `--detailed_results`,
#' `--save_regions`, `--graphical`, `--cpu`). Detection parameters are left at
#' their defaults on purpose: that is what the Portal runs.
#' @param fasta Input FASTA; must sit inside the directory `cmd` was built for.
#' @param out_prefix Output prefix (ITSx appends its own suffixes); same rule.
#' @param cpu Threads.
#' @param cmd list(bin, prefix, message) from [itsx_check_command()].
#' @return list(ok, positions, results, problematic, no_detections, chimeric,
#'   message); `ok = FALSE` and a message when ITSx could not run.
#' @noRd
run_itsx_check <- function(fasta, out_prefix, cpu = 2L, cmd) {
  empty <- list(ok = FALSE, positions = character(0), results = character(0),
    problematic = character(0), no_detections = character(0), chimeric = character(0),
    message = "")
  if (!nzchar(cmd$bin %||% "")) {
    empty$message <- cmd$message %||% "ITSx unavailable"
    return(empty)
  }
  dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
  err <- tempfile("rp_itsx_check_err_")
  on.exit(unlink(err), add = TRUE)
  fasta <- normalizePath(fasta, mustWork = FALSE)
  out_prefix <- file.path(normalizePath(dirname(out_prefix), mustWork = FALSE), basename(out_prefix))
  st <- suppressWarnings(system2(cmd$bin, c(cmd$prefix,
    "-i", shQuote(fasta), "-o", shQuote(out_prefix),
    "--preserve", "T", "--detailed_results", "T", "--save_regions", "none",
    "--graphical", "F", "--cpu", as.integer(cpu)),
    stdout = FALSE, stderr = err))
  if (!identical(as.integer(st), 0L)) {
    empty$message <- paste("ITSx failed:",
      paste(utils::tail(tryCatch(readLines(err, warn = FALSE), error = function(e) ""), 3L),
        collapse = " | "))
    return(empty)
  }
  rd <- function(suffix) {
    f <- paste0(out_prefix, suffix)
    if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
  }
  chim <- rd(".chimeric.fasta")
  list(ok = TRUE,
    positions = rd(".positions.txt"),
    results = rd(".extraction.results"),
    problematic = rd(".problematic.txt"),
    no_detections = trimws(rd("_no_detections.txt")),
    chimeric = sub("^>", "", chim[startsWith(chim, ">")]),
    message = "")
}

#' Per-record notes from an ITSx positions file
#'
#' ITSx appends free-text warnings after the region columns, each ending in
#' "!" (e.g. "Broken or partial sequence, no 5.8S! ", "Chimeric! ").
#' @return data.frame(seqid, notes) with one row per line; `notes` joins the
#'   warnings with "; " (empty when none).
#' @noRd
parse_itsx_notes <- function(pos_lines) {
  pos_lines <- pos_lines[nzchar(trimws(pos_lines))]
  if (length(pos_lines) == 0) {
    return(data.frame(seqid = character(0), notes = character(0), stringsAsFactors = FALSE))
  }
  fields <- strsplit(pos_lines, "\t", fixed = TRUE)
  seqid <- vapply(fields, `[`, character(1), 1)
  notes <- vapply(fields, function(f) {
    rest <- f[-1]
    rest <- rest[!grepl("^\\s*([0-9]+ bp\\.|(SSU|ITS1|5\\.8S|ITS2|LSU):)", rest)]
    tail_txt <- paste(rest, collapse = " ")
    m <- regmatches(tail_txt, gregexpr("[A-Za-z0-9.,' ]+!", tail_txt))[[1]]
    paste(trimws(m), collapse = "; ")
  }, character(1))
  data.frame(seqid = seqid, notes = notes, stringsAsFactors = FALSE)
}

#' Parse an ITSx `.extraction.results` table
#' @return data.frame(seqid, origin, strand, domains, chimeric, domain_order)
#' @noRd
parse_itsx_results <- function(res_lines) {
  res_lines <- res_lines[nzchar(trimws(res_lines))]
  out <- data.frame(seqid = character(0), origin = character(0), strand = integer(0),
    domains = integer(0), chimeric = logical(0), domain_order = character(0),
    stringsAsFactors = FALSE)
  if (length(res_lines) == 0) return(out)
  fields <- strsplit(res_lines, "\t", fixed = TRUE)
  fields <- Filter(function(f) length(f) >= 13, fields)
  if (length(fields) == 0) return(out)
  data.frame(
    seqid = vapply(fields, `[`, character(1), 1),
    origin = vapply(fields, `[`, character(1), 3),
    strand = suppressWarnings(as.integer(vapply(fields, `[`, character(1), 4))),
    domains = suppressWarnings(as.integer(vapply(fields, `[`, character(1), 5))),
    chimeric = vapply(fields, function(f) identical(trimws(f[13]), "Chimeric"), logical(1)),
    domain_order = vapply(fields, function(f) if (length(f) >= 14) trimws(f[14]) else "", character(1)),
    stringsAsFactors = FALSE)
}

#' Parse an ITSx `.problematic.txt` file
#' @return data.frame(seqid, problem)
#' @noRd
parse_itsx_problematic <- function(prob_lines) {
  prob_lines <- prob_lines[nzchar(trimws(prob_lines))]
  if (length(prob_lines) == 0) {
    return(data.frame(seqid = character(0), problem = character(0), stringsAsFactors = FALSE))
  }
  fields <- strsplit(prob_lines, "\t", fixed = TRUE)
  data.frame(seqid = vapply(fields, `[`, character(1), 1),
    problem = vapply(fields, function(f) trimws(paste(f[-1], collapse = " ")), character(1)),
    stringsAsFactors = FALSE)
}

#' ITSx origin letters consistent with a sample's lineage string
#' @return character vector of letters, or NULL when the lineage names no known group.
#' @noRd
itsx_expected_origins <- function(lineage) {
  if (is.null(lineage) || is.na(lineage) || !nzchar(lineage)) return(NULL)
  parts <- trimws(strsplit(lineage, ";", fixed = TRUE)[[1]])
  hit <- names(ITSX_LINEAGE_ORIGINS)[names(ITSX_LINEAGE_ORIGINS) %in% parts]
  if (length(hit) == 0) return(NULL)
  unique(unlist(ITSX_LINEAGE_ORIGINS[hit], use.names = FALSE))
}

#' Portal-parity flags for one record from stock-ITSx output
#'
#' Pure: takes already-parsed ITSx output plus the record's sequence, RiboPilot
#' annotation, and sample lineage.
#' @param seqid Record ID as written in the checked FASTA.
#' @param seq Record sequence (or NA).
#' @param ann data.frame(element, start, stop) of RiboPilot's annotation (or NULL).
#' @param its data.frame from [parse_itsx()] over the positions file.
#' @param notes data.frame from [parse_itsx_notes()].
#' @param res data.frame from [parse_itsx_results()].
#' @param prob data.frame from [parse_itsx_problematic()].
#' @param no_detections IDs ITSx listed as not containing an ITS.
#' @param chimeric_ids IDs ITSx wrote to its chimeric FASTA (domains out of order).
#' @param lineage Sample lineage string ("Metazoa;Chordata;...") or NA.
#' @return list(flags = character vector, note = flags joined with ";").
#' @noRd
itsx_check_flags <- function(seqid, seq, ann, its, notes, res, prob,
                             no_detections = character(0), chimeric_ids = character(0),
                             lineage = NA_character_) {
  flags <- character(0)
  my_its <- if (nrow(its)) its[its$seqid == seqid, , drop = FALSE] else its
  my_note <- if (nrow(notes)) notes$notes[notes$seqid == seqid] else character(0)
  my_res <- if (nrow(res)) res[res$seqid == seqid, , drop = FALSE] else res
  my_prob <- if (nrow(prob)) prob$problem[prob$seqid == seqid] else character(0)
  note_txt <- paste(c(my_note, my_prob), collapse = "; ")

  if (seqid %in% no_detections || (nrow(my_its) == 0 && nrow(my_res) == 0)) {
    return(list(flags = "itsx_not_detected", note = "itsx_not_detected"))
  }

  # chimeric: overlap-based "Chimeric!" note, domain-order chimera, or the
  # problematic-file wording
  if (grepl("Chimeric!", note_txt, fixed = TRUE) || seqid %in% chimeric_ids ||
      any(my_res$chimeric) || grepl("may be chimeric", note_txt, fixed = TRUE)) {
    flags <- c(flags, "itsx_chimeric")
  }
  # misassembly signals reported by ITSx itself
  if (grepl("no 5.8S!", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_no_5.8S")
  if (grepl("only partial 5.8S!", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_partial_5.8S")
  if (grepl("overlap with SSU!", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_5.8S_overlaps_SSU")
  if (grepl("overlap with LSU!", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_5.8S_overlaps_LSU")
  if (grepl("too long!", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_ITS_too_long")
  if (grepl("No ITS1 sequence", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_no_ITS1")
  if (grepl("No ITS2 sequence", note_txt, fixed = TRUE)) flags <- c(flags, "itsx_no_ITS2")

  # gap between ITS regions: an N-run inside the span ITSx would annotate
  span <- my_its[my_its$element %in% c("ITS1", "5.8S", "ITS2"), , drop = FALSE]
  if (nrow(span) > 0 && length(seq) == 1 && !is.na(seq)) {
    gaps <- find_n_runs(seq)
    if (nrow(gaps) > 0 &&
        any(gaps$stop >= min(span$start) & gaps$start <= max(span$stop))) {
      flags <- c(flags, "itsx_gap_in_ITS")
    }
  }

  # boundary disagreement with RiboPilot's own annotation (ITS1/5.8S/ITS2 only;
  # ITSx's SSU/LSU spans run to the sequence ends, so they carry no boundary)
  if (nrow(span) > 0 && !is.null(ann) && nrow(ann) > 0) {
    shifts <- character(0)
    for (el in c("ITS1", "5.8S", "ITS2")) {
      a <- ann[ann$element == el, , drop = FALSE]
      x <- span[span$element == el, , drop = FALSE]
      if (nrow(a) == 0 || nrow(x) == 0) next
      d <- max(abs(x$start[1] - a$start[1]), abs(x$stop[1] - a$stop[1]))
      if (d > ITSX_CHECK_SHIFT_BP) shifts <- c(shifts, sprintf("%s:%dbp", el, d))
    }
    if (length(shifts)) {
      flags <- c(flags, sprintf("itsx_boundary_shift(%s)", paste(shifts, collapse = ",")))
    }
  }

  # origin group vs. sample lineage
  expected <- itsx_expected_origins(lineage)
  if (nrow(my_res) > 0 && !is.null(expected)) {
    o <- my_res$origin[1]
    if (!is.na(o) && nzchar(o) && !(o %in% expected)) {
      flags <- c(flags, sprintf("itsx_origin(%s)", unname(ITSX_ORIGINS[o]) %|NA|% o))
    }
  }
  list(flags = flags, note = paste(flags, collapse = ";"))
}

#' Run the final ITSx check over an export's records
#'
#' Writes the ITSx output under `<gdir>/itsx_check/` and a per-record CSV
#' `<gdir>/<label>_itsx_check.csv`, and returns one row per unit.
#' @param con Project DB connection.
#' @param units data.frame(ID, path, scaffold, seqid) as built by [export_bundle()].
#' @param gdir Bundle directory.
#' @param label Bundle label.
#' @param cpu Threads for ITSx.
#' @param project Project directory; its `.config` decides the container.
#' @param cmd Override for [itsx_check_command()] (tests).
#' @return data.frame(seqid, itsx_check, itsx_origin, itsx_ITS1, itsx_5.8S,
#'   itsx_ITS2, itsx_notes). `itsx_check` is "" when clean, the flag list when
#'   not, or "not_run(<why>)" when ITSx could not be run.
#' @noRd
itsx_export_check <- function(con, units, gdir, label, cpu = 2L, project = ".", cmd = NULL) {
  n <- nrow(units)
  blank <- data.frame(seqid = units$seqid, itsx_check = rep("", n),
    itsx_origin = rep(NA_character_, n), itsx_ITS1 = rep(NA_character_, n),
    `itsx_5.8S` = rep(NA_character_, n), itsx_ITS2 = rep(NA_character_, n),
    itsx_notes = rep("", n), stringsAsFactors = FALSE, check.names = FALSE)
  if (n == 0) return(blank)

  seqs <- vapply(seq_len(n), function(i) {
    r <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
    if (nrow(r)) r$sequence[1] else NA_character_
  }, character(1))
  cdir <- file.path(gdir, "itsx_check")
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  fa <- file.path(cdir, "records.fasta")
  ok <- !is.na(seqs) & nzchar(seqs)
  writeLines(as.vector(rbind(paste0(">", units$seqid[ok]), seqs[ok])), fa)

  if (is.null(cmd)) cmd <- itsx_check_command(project, cdir)
  out <- run_itsx_check(fa, file.path(cdir, "itsx"), cpu = cpu, cmd = cmd)
  if (!isTRUE(out$ok)) {
    blank$itsx_check <- sprintf("not_run(%s)", out$message)
    return(blank)
  }
  its <- parse_itsx(out$positions)
  if (nrow(its) == 0) {
    its <- data.frame(seqid = character(0), element = character(0),
      start = integer(0), stop = integer(0), stringsAsFactors = FALSE)
  }
  notes <- parse_itsx_notes(out$positions)
  res <- parse_itsx_results(out$results)
  prob <- parse_itsx_problematic(out$problematic)

  ids <- unique(as.character(units$ID))
  lin <- tryCatch(DBI::dbGetQuery(con,
    sprintf("SELECT ID, lineage FROM samples WHERE ID IN (%s)",
      paste(rep("?", length(ids)), collapse = ",")), params = as.list(ids)),
    error = function(e) data.frame(ID = character(0), lineage = character(0)))

  rng <- function(sid, el) {
    x <- its[its$seqid == sid & its$element == el, , drop = FALSE]
    if (nrow(x)) sprintf("%d-%d", x$start[1], x$stop[1]) else NA_character_
  }
  for (i in seq_len(n)) {
    sid <- units$seqid[i]
    ann <- DBI::dbGetQuery(con,
      "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(units$ID[i], units$path[i], units$scaffold[i]))
    lineage <- lin$lineage[match(units$ID[i], lin$ID)]
    fl <- itsx_check_flags(sid, seqs[i], ann, its, notes, res, prob,
      no_detections = out$no_detections, chimeric_ids = out$chimeric,
      lineage = if (length(lineage)) lineage else NA_character_)
    blank$itsx_check[i] <- fl$note
    o <- res$origin[match(sid, res$seqid)]
    blank$itsx_origin[i] <- if (length(o) && !is.na(o)) unname(ITSX_ORIGINS[o]) %|NA|% o else NA_character_
    blank$itsx_ITS1[i] <- rng(sid, "ITS1")
    blank$`itsx_5.8S`[i] <- rng(sid, "5.8S")
    blank$itsx_ITS2[i] <- rng(sid, "ITS2")
    blank$itsx_notes[i] <- paste(c(notes$notes[notes$seqid == sid], prob$problem[prob$seqid == sid]),
      collapse = "; ")
  }
  utils::write.csv(blank, file.path(gdir, paste0(label, "_itsx_check.csv")), row.names = FALSE)
  blank
}
