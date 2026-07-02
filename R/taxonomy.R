# NCBI taxonomy resolution for sample Taxon strings.
#
# Free-text `Taxon` values are resolved to an NCBI taxid + full lineage so the
# scaffold/path scorer can measure how closely a candidate's top-BLAST-hit
# lineage matches the sample.
#
# Reproducibility: resolution uses a per-project snapshot of the NCBI taxonomy
# (taxdump) downloaded at project init and pinned in the project directory, so a
# project always resolves against the same taxonomy version even as NCBI changes.
# The confirmed taxid + lineage are also frozen in the project database, so
# downstream steps never re-query taxonomy at all.
#
# Resolution is typo-tolerant: exact names resolve directly; a mistyped
# "Genus species" falls back to genus-anchored fuzzy matching (edit distance
# against that genus's species); genuine non-matches are flagged. EVERY proposal
# is surfaced for user confirmation, not only the failures.

RRNA_TAX_RANKS <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

#' Normalize a free-text taxon string to a query-able name
#'
#' Strips open-nomenclature qualifiers (sp., cf., aff., nr.), parenthetical
#' content, and non-letter junk, keeping "Genus species".
#' @noRd
clean_taxon <- function(x) {
  x <- trimws(x)
  x <- gsub("[_]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- sub("\\b(sp|spp|cf|aff|nr|indet)\\.?\\b.*$", "", x, ignore.case = TRUE)
  x <- gsub("\\([^)]*\\)", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  toks <- strsplit(x, " ")[[1]]
  toks <- toks[grepl("^[A-Za-z.-]+$", toks)]
  paste(utils::head(toks, 2), collapse = " ")
}

# ---- Per-project local taxonomy (taxdump + taxonkit) ------------------------

#' Directory holding a project's pinned NCBI taxdump
#' @param project Project directory.
#' @noRd
project_taxdump_dir <- function(project) file.path(project, ".taxonomy")

#' Path to a project's pinned taxdump tarball
#' @noRd
taxdump_tarball <- function(dir) file.path(dir, "taxdump.tar.gz")

#' Is a usable taxdump present in `dir`? (pinned tarball, or pre-extracted dmps)
#' @noRd
have_taxdump <- function(dir) {
  if (is.null(dir)) {
    return(FALSE)
  }
  file.exists(taxdump_tarball(dir)) ||
    (file.exists(file.path(dir, "names.dmp")) && file.exists(file.path(dir, "nodes.dmp")))
}

#' Extract the dmp files taxonkit needs from a pinned tarball into a temp dir
#'
#' The ~500 MB of extracted dmp files are not persisted in the project; they are
#' unpacked to a temporary location for the duration of one resolution call and
#' then deleted by the caller.
#' @return temp directory containing names.dmp / nodes.dmp / merged.dmp / delnodes.dmp
#' @noRd
extract_taxdump <- function(dir) {
  tmp <- tempfile("ribopilot_taxdump_")
  dir.create(tmp)
  utils::untar(taxdump_tarball(dir),
    files = c("names.dmp", "nodes.dmp", "merged.dmp", "delnodes.dmp"),
    exdir = tmp)
  tmp
}

#' Locate the taxonkit binary
#' @noRd
taxonkit_bin <- function() {
  getOption("ribopilot.taxonkit_bin", unname(Sys.which("taxonkit")))
}

#' Download a pinned NCBI taxdump tarball into `dir` (once per project)
#'
#' Stores only the ~72 MB compressed tarball plus a provenance file. The dmp
#' files are extracted to a temp location at resolution time (see
#' [extract_taxdump()]), not persisted in the project.
#'
#' @param dir Destination directory (e.g. `<project>/.taxonomy`).
#' @param force Re-download even if a tarball already exists.
#' @return `dir`.
#' @export
ensure_taxdump <- function(dir, force = FALSE) {
  if (file.exists(taxdump_tarball(dir)) && !force) {
    return(dir)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tgz <- taxdump_tarball(dir)
  url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  message("Downloading NCBI taxdump snapshot (~72 MB) to ", dir, " ...")
  if (utils::download.file(url, tgz, mode = "wb", quiet = TRUE) != 0) {
    stop("Failed to download taxdump from ", url, call. = FALSE)
  }
  writeLines(
    c(
      paste0("source: ", url),
      paste0("downloaded: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste0("tarball_mtime: ", format(file.info(tgz)$mtime, "%Y-%m-%d %H:%M:%S %Z"))
    ),
    file.path(dir, "TAXDUMP_VERSION.txt")
  )
  dir
}

#' Run taxonkit with a given data-dir, feeding `input` on stdin
#' @noRd
tk_run <- function(args, input, dir) {
  bin <- taxonkit_bin()
  if (!nzchar(bin)) {
    stop("'taxonkit' not found. Install taxonkit (see environment.yml).", call. = FALSE)
  }
  # Always route through /bin/sh (system2 uses a shell iff `input` is non-NULL)
  # so shell-quoting of args is consistent. The reformat format string contains
  # ';' and '{}', which must be quoted.
  suppressWarnings(system2(
    bin, c(shQuote(args), "--data-dir", shQuote(dir)),
    input = if (is.null(input)) "" else input, stdout = TRUE, stderr = FALSE
  ))
}

#' Map names -> taxids locally (exact; includes synonyms). First taxid per name.
#' @return named integer vector (input name -> taxid, NA if unmatched)
#' @noRd
tk_name2taxid <- function(names, dir) {
  names <- unique(names)
  if (length(names) == 0) {
    return(integer(0))
  }
  out <- tk_run("name2taxid", paste(names, collapse = "\n"), dir)
  res <- setNames(rep(NA_integer_, length(names)), names)
  for (p in strsplit(out, "\t")) {
    nm <- p[1]
    tx <- suppressWarnings(as.integer(p[2] %||% NA))
    if (!is.na(nm) && nm %in% names && !is.na(tx) && is.na(res[[nm]])) res[[nm]] <- tx
  }
  res
}

#' Local taxid info: sci_name, rank, ranked lineage (batched: two taxonkit calls)
#' @return data.frame(taxid, sci_name, rank, lineage)
#' @noRd
tk_taxid_info <- function(taxids, dir) {
  taxids <- unique(stats::na.omit(taxids))
  if (length(taxids) == 0) {
    return(data.frame(taxid = integer(0), sci_name = character(0), rank = character(0), lineage = character(0)))
  }
  ids <- paste(taxids, collapse = "\n")
  nr <- tk_run(c("lineage", "-n", "-r"), ids, dir) # taxid \t lineage \t name \t rank
  nr <- do.call(rbind, lapply(strsplit(nr, "\t"), function(p) {
    if (length(p) < 4) {
      return(NULL)
    }
    data.frame(taxid = as.integer(p[1]), sci_name = p[length(p) - 1], rank = p[length(p)], stringsAsFactors = FALSE)
  }))
  lin <- tk_run("lineage", ids, dir) # taxid \t lineage
  ref <- tk_run(c("reformat", "-f", "{K};{p};{c};{o};{f};{g};{s}"), paste(lin, collapse = "\n"), dir)
  ref <- do.call(rbind, lapply(strsplit(ref, "\t"), function(p) {
    if (length(p) < 3) {
      return(NULL)
    }
    data.frame(taxid = as.integer(p[1]),
      lineage = trimws(gsub("; *;", ";", gsub("(^; *|; *$)", "", gsub(";", "; ", p[length(p)])))),
      stringsAsFactors = FALSE)
  }))
  merge(nr, ref, by = "taxid", all.x = TRUE)
}

#' Local species list under a genus taxid (for fuzzy epithet matching)
#' @return data.frame(taxid, name)
#' @noRd
tk_genus_species <- function(genus_taxid, dir) {
  out <- tk_run(c("list", "--ids", as.character(genus_taxid), "--indent", "", "--show-name", "--show-rank"), NULL, dir)
  out <- out[grepl("\\[species\\]", out)]
  if (length(out) == 0) {
    return(data.frame(taxid = integer(0), name = character(0)))
  }
  m <- regmatches(out, regexec("^(\\d+)\\s+\\[species\\]\\s+(.*)$", out))
  m <- Filter(function(x) length(x) == 3, m)
  data.frame(
    taxid = as.integer(vapply(m, `[`, character(1), 2)),
    name = vapply(m, `[`, character(1), 3),
    stringsAsFactors = FALSE
  )
}

#' Batched local resolution of many Taxon strings
#' @noRd
resolve_taxonomy_local <- function(taxa, dir) {
  clean <- vapply(taxa, clean_taxon, character(1), USE.NAMES = FALSE)
  res <- data.frame(
    Taxon = taxa, clean = clean, match_type = "unresolved",
    taxid = NA_integer_, sci_name = NA_character_, rank = NA_character_,
    lineage = NA_character_, alternatives = NA_character_, stringsAsFactors = FALSE
  )
  ok <- nzchar(clean)
  if (!any(ok)) {
    return(res)
  }

  # --- exact pass (one name2taxid + one info call for ALL exact names) ---
  n2t <- tk_name2taxid(clean[ok], dir)
  exact_tx <- unname(n2t[clean])
  hit <- ok & !is.na(exact_tx)
  if (any(hit)) {
    info <- tk_taxid_info(exact_tx[hit], dir)
    m <- info[match(exact_tx[hit], info$taxid), ]
    res$match_type[hit] <- "exact"
    res$taxid[hit] <- exact_tx[hit]
    res$sci_name[hit] <- m$sci_name
    res$rank[hit] <- m$rank
    res$lineage[hit] <- m$lineage
  }

  # --- fuzzy pass (genus-anchored) for two-token misses ---
  miss <- which(ok & !hit)
  if (length(miss) > 0) {
    genus <- vapply(clean[miss], function(x) {
      t <- strsplit(x, " ")[[1]]
      if (length(t) == 2) t[1] else NA_character_
    }, character(1), USE.NAMES = FALSE)
    gtx_map <- tk_name2taxid(stats::na.omit(unique(genus)), dir)
    chosen <- integer(length(miss)) # taxid to look up lineage for, at the end
    chosen[] <- NA_integer_
    for (k in seq_along(miss)) {
      g <- genus[k]
      if (is.na(g) || is.na(gtx_map[[g]] %||% NA)) next
      gt <- gtx_map[[g]]
      epithet <- strsplit(clean[miss[k]], " ")[[1]][2]
      sp <- tk_genus_species(gt, dir)
      if (nrow(sp) > 0) {
        sp$dist <- as.integer(adist(tolower(clean[miss[k]]), tolower(sp$name)))
        sp <- sp[order(sp$dist), ]
        if (sp$dist[1] <= max(2, ceiling(nchar(epithet) / 3))) {
          res$match_type[miss[k]] <- "fuzzy_species"
          res$sci_name[miss[k]] <- sp$name[1]
          res$alternatives[miss[k]] <- paste0(utils::head(sp$name, 5), " [", utils::head(sp$taxid, 5), "]", collapse = " | ")
          chosen[k] <- sp$taxid[1]
          next
        }
      }
      res$match_type[miss[k]] <- "fuzzy_genus"
      res$sci_name[miss[k]] <- g
      chosen[k] <- gt
    }
    # one info call for all fuzzy picks
    pick <- !is.na(chosen)
    if (any(pick)) {
      info <- tk_taxid_info(chosen[pick], dir)
      m <- info[match(chosen[pick], info$taxid), ]
      idx <- miss[pick]
      res$taxid[idx] <- chosen[pick]
      res$rank[idx] <- m$rank
      res$lineage[idx] <- m$lineage
    }
  }
  res
}

# ---- Network fallback (NCBI `datasets`, no local taxdump) --------------------

#' Locate the NCBI `datasets` binary
#' @noRd
datasets_bin <- function() {
  getOption("ribopilot.datasets_bin", unname(Sys.which("datasets")))
}

#' Run `datasets summary taxonomy taxon <query>`, stdout + stderr separately
#' @noRd
run_datasets_taxon <- function(query) {
  bin <- datasets_bin()
  if (!nzchar(bin)) {
    stop("NCBI 'datasets' CLI not found. Install ncbi-datasets-cli (see environment.yml).", call. = FALSE)
  }
  errf <- tempfile()
  on.exit(unlink(errf), add = TRUE)
  out <- suppressWarnings(system2(bin, c("summary", "taxonomy", "taxon", shQuote(query)),
    stdout = TRUE, stderr = errf))
  list(stdout = paste(out, collapse = "\n"), stderr = paste(readLines(errf, warn = FALSE), collapse = "\n"))
}

#' Lineage data.frame (rank, id, name) from a datasets classification
#' @noRd
classification_to_lineage <- function(cl) {
  ranks <- intersect(RRNA_TAX_RANKS, names(cl))
  if (length(ranks) == 0) {
    return(data.frame(rank = character(0), id = integer(0), name = character(0)))
  }
  data.frame(
    rank = ranks,
    id = vapply(ranks, function(r) as.integer(cl[[r]]$id %||% NA), integer(1)),
    name = vapply(ranks, function(r) as.character(cl[[r]]$name %||% NA), character(1)),
    stringsAsFactors = FALSE
  )
}

#' Fetch classification for an exact name or a taxid via datasets
#' @noRd
datasets_exact <- function(query) {
  raw <- run_datasets_taxon(query)$stdout
  raw <- sub("^[^{]*", "", raw)
  parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$reports)) {
    return(NULL)
  }
  tx <- parsed$reports[[1]]$taxonomy
  lin <- classification_to_lineage(tx$classification)
  finest <- if (nrow(lin) > 0) lin[nrow(lin), ] else NULL
  list(
    taxid = if (!is.null(finest)) finest$id else NA_integer_,
    sci_name = tx$current_scientific_name$name %||% (if (!is.null(finest)) finest$name else NA_character_),
    rank = if (!is.null(finest)) finest$rank else NA_character_,
    lineage = lin
  )
}

#' Resolve one Taxon via network datasets (exact -> suggestions -> unresolved)
#' @noRd
resolve_taxon_one <- function(taxon) {
  clean <- clean_taxon(taxon)
  empty <- data.frame(
    Taxon = taxon, clean = clean, match_type = "unresolved",
    taxid = NA_integer_, sci_name = NA_character_, rank = NA_character_,
    lineage = NA_character_, alternatives = NA_character_, stringsAsFactors = FALSE
  )
  if (!nzchar(clean)) {
    return(empty)
  }
  raw <- run_datasets_taxon(clean)
  lineage_str <- function(lin) if (is.data.frame(lin) && nrow(lin) > 0) paste(lin$name, collapse = "; ") else NA_character_
  if (grepl("\"reports\"", raw$stdout, fixed = TRUE)) {
    ex <- datasets_exact(clean)
    if (!is.null(ex) && !is.na(ex$taxid)) {
      out <- empty
      out$match_type <- "exact"
      out$taxid <- ex$taxid
      out$sci_name <- ex$sci_name
      out$rank <- ex$rank
      out$lineage <- lineage_str(ex$lineage)
      return(out)
    }
  }
  empty # network fuzzy fallback intentionally minimal; prefer the local path
}

# ---- Public entry point -----------------------------------------------------

#' Resolve sample Taxon strings to NCBI taxids + lineages (for confirmation)
#'
#' Prefers a pinned local taxdump (fast, offline, reproducible) when
#' `taxdump_dir` has one; otherwise falls back to live NCBI `datasets`. EVERY row
#' is returned for user confirmation (see `match_type`), not only failures.
#'
#' @param taxa Character vector of Taxon strings (duplicates de-duplicated).
#' @param taxdump_dir Directory of a pinned taxdump (see [ensure_taxdump()]). If
#'   it contains a taxdump, resolution is local; otherwise network `datasets` is
#'   used.
#' @return data.frame, one row per unique Taxon: `match_type` in
#'   "exact"/"fuzzy_species"/"fuzzy_genus"/"unresolved", plus proposed `taxid`,
#'   `sci_name`, `rank`, `lineage`, `alternatives`.
#' @export
resolve_taxonomy <- function(taxa, taxdump_dir = NULL) {
  uniq <- unique(stats::na.omit(taxa))
  if (length(uniq) == 0) {
    return(data.frame())
  }
  if (have_taxdump(taxdump_dir)) {
    # already-extracted dmps -> use directly; otherwise unpack the pinned tarball
    # to a temp dir for this call only, then clean up.
    if (file.exists(file.path(taxdump_dir, "names.dmp"))) {
      data_dir <- taxdump_dir
    } else {
      data_dir <- extract_taxdump(taxdump_dir)
      on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)
    }
    return(resolve_taxonomy_local(uniq, data_dir))
  }
  do.call(rbind, lapply(uniq, resolve_taxon_one))
}
