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

# ---- Per-project local taxonomy (pinned taxdump) ----------------------------

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

#' Extract the dmp files from a pinned tarball into a temp dir
#'
#' The ~500 MB of extracted dmp files are not persisted in the project; they are
#' unpacked to a temporary location only while the cached index is built (see
#' [build_tax_index()]) and then deleted by the caller.
#' @return temp directory containing names.dmp / nodes.dmp / merged.dmp / delnodes.dmp
#' @noRd
extract_taxdump <- function(dir) {
  tmp <- tempfile("ribopilot_taxdump_")
  dir.create(tmp)
  st <- suppressWarnings(utils::untar(taxdump_tarball(dir),
    files = c("names.dmp", "nodes.dmp", "merged.dmp", "delnodes.dmp"), exdir = tmp))
  if (!identical(as.integer(st), 0L)) {
    unlink(tmp, recursive = TRUE)
    stop("Could not unpack ", taxdump_tarball(dir), " (tar status ", st,
      "); the file is incomplete or corrupt. Delete it and try again.", call. = FALSE)
  }
  tmp
}

#' Download a pinned NCBI taxdump tarball into `dir` (once per project)
#'
#' Stores only the about 80 MB compressed tarball plus a provenance file. The
#' dmp files are extracted to a temp location once, when the project's
#' taxonomy index is built (see [build_tax_index()]), and are not persisted.
#'
#' @param dir Destination directory (e.g. `<project>/.taxonomy`).
#' @param force Re-download even if a tarball already exists.
#' @return `dir`.
#' @noRd
ensure_taxdump <- function(dir, force = FALSE) {
  tgz <- taxdump_tarball(dir)
  if (file.exists(tgz) && !force) {
    return(dir)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  message("Downloading the NCBI taxonomy (about 80 MB) to ", dir, " ...")
  part <- paste0(tgz, ".part")
  on.exit(unlink(c(part, paste0(part, ".md5"))), add = TRUE)
  download_or_stop(url, part, "the NCBI taxonomy")
  download_or_stop(paste0(url, ".md5"), paste0(part, ".md5"), "the NCBI taxonomy checksum")
  want <- tolower(sub("\\s.*$", "", readLines(paste0(part, ".md5"), warn = FALSE)[1]))
  got <- unname(tools::md5sum(part))
  if (!identical(want, got)) {
    stop("The NCBI taxonomy download is incomplete (checksum mismatch). Try again.", call. = FALSE)
  }
  file.rename(part, tgz)
  unlink(tax_index_path(dir))
  writeLines(
    c(
      paste0("source: ", url),
      paste0("md5: ", got),
      paste0("downloaded: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste0("tarball_mtime: ", format(file.info(tgz)$mtime, "%Y-%m-%d %H:%M:%S %Z"))
    ),
    file.path(dir, "TAXDUMP_VERSION.txt")
  )
  dir
}

#' Kingdom name of each rank-aligned lineage string (NA where there is none)
#' @noRd
rank_kingdom <- function(lineage) {
  vapply(lineage, function(l) {
    if (is.na(l)) NA_character_ else unname(parse_rank_names(l)[["kingdom"]])
  }, character(1), USE.NAMES = FALSE)
}

#' Describe homonym candidates for review: "name [taxid] (kingdom)" per taxid
#' @noRd
describe_candidates <- function(taxids, con) {
  info <- td_taxid_info(con, taxids)
  m <- info[match(taxids, info$taxid), ]
  k <- rank_kingdom(m$lineage)
  lab <- paste0(ifelse(is.na(m$sci_name), "?", m$sci_name), " [", taxids, "]")
  paste(ifelse(is.na(k), lab, paste0(lab, " (", k, ")")), collapse = " | ")
}

#' Local taxid info: sci_name, rank, ranked lineage (one batched index query)
#'
#' The lineage is the same seven rank-aligned fields [tk_ranked_lineage()]
#' produces, names and taxids alike. Alignment is the whole point: consumers read
#' these strings positionally ([parse_rank_names()], [parse_rank_taxids()]), so a
#' taxon missing a rank must carry a filled placeholder there rather than have
#' the slot collapsed away and every deeper rank shifted up one.
#' @return data.frame(taxid, sci_name, rank, lineage, lineage_taxids)
#' @noRd
td_taxid_info <- function(con, taxids) {
  taxids <- unique(stats::na.omit(suppressWarnings(as.integer(taxids))))
  empty <- data.frame(taxid = integer(0), sci_name = character(0), rank = character(0),
    lineage = character(0), lineage_taxids = character(0), stringsAsFactors = FALSE)
  if (length(taxids) == 0) {
    return(empty)
  }
  paths <- td_paths(con, taxids)
  gone <- setdiff(as.character(taxids), names(paths))
  if (length(gone) > 0) {
    warning("No lineage returned for taxid(s): ", paste(gone, collapse = ", "), call. = FALSE)
  }
  if (length(paths) == 0) {
    return(empty)
  }
  rows <- lapply(names(paths), function(k) {
    p <- paths[[k]]
    tip <- p[nrow(p), ]
    al <- td_rank_align(p)
    data.frame(taxid = as.integer(k), sci_name = tip$name, rank = tip$rank,
      lineage = al$names, lineage_taxids = al$taxids, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Local taxid info by taxdump directory
#' @noRd
tk_taxid_info <- function(taxids, taxdump_dir) {
  empty <- data.frame(taxid = integer(0), sci_name = character(0), rank = character(0),
    lineage = character(0), lineage_taxids = character(0), stringsAsFactors = FALSE)
  with_tax_index(taxdump_dir, function(con) td_taxid_info(con, taxids), default = empty)
}

#' Local species list under a genus taxid (for fuzzy epithet matching)
#' @return data.frame(taxid, name)
#' @noRd
tk_genus_species <- function(genus_taxid, con) {
  td_descendant_species(con, genus_taxid)
}

#' Batched local resolution of many Taxon strings
#' @param con Open connection to the project's taxonomy index.
#' @noRd
resolve_taxonomy_local <- function(taxa, con) {
  clean <- vapply(taxa, clean_taxon, character(1), USE.NAMES = FALSE)
  res <- data.frame(
    Taxon = taxa, clean = clean, match_type = "unresolved",
    taxid = NA_integer_, sci_name = NA_character_, rank = NA_character_,
    lineage = NA_character_, lineage_taxids = NA_character_,
    alternatives = NA_character_, stringsAsFactors = FALSE
  )
  ok <- nzchar(clean)
  if (!any(ok)) {
    return(res)
  }

  # --- exact pass ---
  cand <- td_name2taxid(con, clean[ok])
  n_cand <- vapply(clean, function(x) length(cand[[x]]), integer(1), USE.NAMES = FALSE)
  n_cand[!ok] <- 0L
  amb <- n_cand > 1
  if (any(amb)) {
    res$match_type[amb] <- "ambiguous"
    res$alternatives[amb] <- vapply(clean[amb],
      function(x) describe_candidates(cand[[x]], con), character(1), USE.NAMES = FALSE)
  }
  hit <- n_cand == 1
  if (any(hit)) {
    tx <- vapply(clean[hit], function(x) cand[[x]], integer(1), USE.NAMES = FALSE)
    info <- td_taxid_info(con, tx)
    m <- info[match(tx, info$taxid), ]
    res$match_type[hit] <- "exact"
    res$taxid[hit] <- tx
    res$sci_name[hit] <- m$sci_name
    res$rank[hit] <- m$rank
    res$lineage[hit] <- m$lineage
    res$lineage_taxids[hit] <- m$lineage_taxids
  }

  # --- fuzzy pass (genus-anchored) for two-token misses ---
  miss <- which(ok & n_cand == 0)
  if (length(miss) > 0) {
    genus <- vapply(clean[miss], function(x) {
      t <- strsplit(x, " ")[[1]]
      if (length(t) == 2) t[1] else NA_character_
    }, character(1), USE.NAMES = FALSE)
    gcand <- td_name2taxid(con, stats::na.omit(unique(genus)))
    chosen <- rep(NA_integer_, length(miss))
    for (k in seq_along(miss)) {
      g <- genus[k]
      gts <- if (is.na(g)) integer(0) else gcand[[g]] %||% integer(0)
      if (length(gts) == 0) next
      epithet <- strsplit(clean[miss[k]], " ")[[1]][2]
      sp <- do.call(rbind, lapply(gts, tk_genus_species, con = con))
      if (nrow(sp) > 0) {
        sp$dist <- as.integer(utils::adist(tolower(clean[miss[k]]), tolower(sp$name)))
        sp <- sp[order(sp$dist), ]
        if (sp$dist[1] <= max(2, ceiling(nchar(epithet) / 3))) {
          res$match_type[miss[k]] <- "fuzzy_species"
          res$sci_name[miss[k]] <- sp$name[1]
          res$alternatives[miss[k]] <- paste0(utils::head(sp$name, 5), " [", utils::head(sp$taxid, 5), "]", collapse = " | ")
          chosen[k] <- sp$taxid[1]
          next
        }
      }
      if (length(gts) > 1) {
        res$match_type[miss[k]] <- "ambiguous"
        res$alternatives[miss[k]] <- describe_candidates(gts, con)
        next
      }
      res$match_type[miss[k]] <- "fuzzy_genus"
      res$sci_name[miss[k]] <- g
      chosen[k] <- gts
    }
    pick <- !is.na(chosen)
    if (any(pick)) {
      info <- td_taxid_info(con, chosen[pick])
      m <- info[match(chosen[pick], info$taxid), ]
      idx <- miss[pick]
      res$taxid[idx] <- chosen[pick]
      res$rank[idx] <- m$rank
      res$lineage[idx] <- m$lineage
      res$lineage_taxids[idx] <- m$lineage_taxids
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

#' Lay a (rank, id, name) classification onto the seven rank-aligned slots
#'
#' The network path used to join whatever ranks NCBI happened to return, so a
#' taxon with no class produced a shorter string and every rank below it shifted
#' up one slot. Consumers read these strings positionally, so absent ranks must
#' be left empty rather than closed up. Domain is dropped: [RRNA_LINEAGE_RANKS]
#' starts at kingdom.
#' @param lin data.frame(rank, id, name) from [classification_to_lineage()].
#' @return list(names, taxids), each a ";"-joined 7-field string.
#' @noRd
rank_aligned_lineage <- function(lin) {
  if (!is.data.frame(lin) || nrow(lin) == 0) {
    return(list(names = NA_character_, taxids = NA_character_))
  }
  i <- match(RRNA_LINEAGE_RANKS, lin$rank)
  fill <- function(x) paste(ifelse(is.na(i) | is.na(x[i]), "", as.character(x[i])), collapse = ";")
  list(names = fill(lin$name), taxids = fill(lin$id))
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
    lineage = NA_character_, lineage_taxids = NA_character_,
    alternatives = NA_character_, stringsAsFactors = FALSE
  )
  if (!nzchar(clean)) {
    return(empty)
  }
  raw <- run_datasets_taxon(clean)
  if (grepl("\"reports\"", raw$stdout, fixed = TRUE)) {
    ex <- datasets_exact(clean)
    if (!is.null(ex) && !is.na(ex$taxid)) {
      out <- empty
      out$match_type <- "exact"
      out$taxid <- ex$taxid
      out$sci_name <- ex$sci_name
      out$rank <- ex$rank
      aligned <- rank_aligned_lineage(ex$lineage)
      out$lineage <- aligned$names
      out$lineage_taxids <- aligned$taxids
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
#'   it contains a taxdump, resolution is local (an unreadable taxdump is an
#'   error); otherwise network `datasets` is used.
#' @return data.frame, one row per unique Taxon: `match_type` in
#'   "exact"/"ambiguous"/"fuzzy_species"/"fuzzy_genus"/"unresolved", plus `taxid`,
#'   `sci_name`, `rank`, `lineage`, `alternatives`. An ambiguous name (several
#'   taxa share it) gets no taxid; \code{alternatives} lists each candidate
#'   with its kingdom.
#' @export
resolve_taxonomy <- function(taxa, taxdump_dir = NULL) {
  uniq <- unique(stats::na.omit(taxa))
  if (length(uniq) == 0) {
    return(data.frame())
  }
  if (have_taxdump(taxdump_dir)) {
    local <- with_tax_index(taxdump_dir, function(con) resolve_taxonomy_local(uniq, con))
    if (is.null(local)) {
      stop("Could not read the NCBI taxonomy in ", taxdump_dir, ". Delete that folder, then ",
        "press Update in the Taxonomy step to download a fresh copy.", call. = FALSE)
    }
    return(local)
  }
  do.call(rbind, lapply(uniq, resolve_taxon_one))
}
