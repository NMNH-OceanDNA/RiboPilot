# GenBank submission bundle export (mirrors MitoPilot's group-based export).
#
# Units are organized into export groups (annotate.export_group). Each group /
# selection is written to <out_dir>/export/<label>/:
#
#   <label>.fasta          combined records, deflines from the header template
#   <label>.tbl            concatenated 5-column feature tables
#   <label>_summary.csv    one row per exported unit
#   <label>_report.html    self-contained HTML report
#   README.txt             how the bundle maps onto NCBI's submission routes
#   samples/               <seqid>.fasta and <seqid>.tbl, one pair per unit
#   markers/<ELEMENT>/     <seqid>_<ELEMENT>.fasta and all_<ELEMENT>.fasta
#
# Deflines come from the named header configurations in export_opts.R.
# table2asn is run per sample file if available.

#' Assign selected units to an export group (export_group is per unit now)
#'
#' @param path Project dir.
#' @param units data.frame(ID, path, scaffold) OR a character vector of sample
#'   IDs (assigns all of each sample's units).
#' @param group Group name.
#' @noRd
assign_export_group <- function(path, units, group) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  set_export_group(con, units, group)
  invisible()
}

#' Clear the export group of selected units
#' @param path Project dir.
#' @param units data.frame(ID, path, scaffold) or a character vector of IDs.
#' @noRd
clear_export_group <- function(path, units) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  set_export_group(con, units, NA)
  invisible()
}

#' Write annotate.export_group for the given units (NA clears). Accepts a
#' data.frame(ID,path,scaffold) or a character vector of sample IDs. `group` is
#' one name for all units, or one name per unit.
#' @noRd
set_export_group <- function(con, units, group) {
  n <- if (is.data.frame(units)) nrow(units) else length(units)
  if (n == 0) return(invisible())
  g <- if (length(group) == 0) NA_character_ else as.character(group)
  g <- if (length(g) == 1) rep(g, n) else g
  stopifnot(length(g) == n)
  # all-or-nothing: a partial group assignment is worse than none
  DBI::dbBegin(con)
  done <- FALSE
  on.exit(if (!done) DBI::dbRollback(con))
  if (is.data.frame(units)) {
    DBI::dbExecute(con, "UPDATE annotate SET export_group = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(g, as.character(units$ID), as.integer(units$path), as.integer(units$scaffold)))
  } else {
    DBI::dbExecute(con, "UPDATE annotate SET export_group = ? WHERE ID = ?",
      params = list(g, as.character(units)))
  }
  DBI::dbCommit(con)
  done <- TRUE
  invisible()
}

#' Assembler that produced each sample (assemble.assemble_opts -> assemble_opts.assembler)
#'
#' A GenBank submission carries one assembly method, so a group should not mix
#' assemblers. Missing rows fall back to "getorganelle", the workflow default.
#' @param path Project dir. @param ids Sample IDs.
#' @return character vector aligned to `ids`.
#' @noRd
sample_assembler <- function(path, ids) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  a <- tryCatch(DBI::dbGetQuery(con,
    "SELECT a.ID, COALESCE(ao.assembler, 'getorganelle') AS assembler
     FROM assemble a LEFT JOIN assemble_opts ao ON ao.assemble_opts = a.assemble_opts"),
    error = function(e) data.frame(ID = character(0), assembler = character(0)))
  out <- a$assembler[match(ids, a$ID)]
  out[is.na(out)] <- "getorganelle"
  out
}

# Minimum run of Ns annotated as an assembly_gap; shorter runs are treated as
# ordinary ambiguous bases, not gaps (per INSDC convention).
MIN_ASSEMBLY_GAP <- 10L

#' 1-based start/stop of every N-run of at least `min_len` in a sequence
#' @noRd
find_n_runs <- function(seq, min_len = MIN_ASSEMBLY_GAP) {
  empty <- data.frame(start = integer(0), stop = integer(0))
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) return(empty)
  m <- gregexpr(sprintf("[Nn]{%d,}", min_len), seq, perl = TRUE)[[1]]
  if (length(m) == 1 && m[1] == -1) return(empty)
  len <- attr(m, "match.length")
  data.frame(start = as.integer(m), stop = as.integer(m + len - 1))
}

#' Feature table (.tbl) text for one assembly unit from an open connection
#'
#' Emits the rRNA/spacer features plus an `assembly_gap` feature for each N-run
#' (>= [MIN_ASSEMBLY_GAP]) so table2asn does not reject fabricated gap runs, and
#' marks the first/last real feature partial (`<`/`>`) when it reaches a
#' sequence end or when `force_partial` is set (a curator-flagged partial repeat,
#' so it is not submitted as a complete cistron). `seqid` is the FASTA/feature id.
#' @noRd
feature_table_text_con <- function(con, id, path = 1L, scaffold = 1L, seqid = id, force_partial = FALSE) {
  d <- DBI::dbGetQuery(con,
    "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  hdr <- paste0(">Feature ", seqid)
  if (nrow(d) == 0) {
    return(paste0(hdr, "\n"))
  }
  st <- DBI::dbGetQuery(con, "SELECT sequence, length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  seq <- if (nrow(st)) st$sequence[1] else NA_character_
  seqlen <- if (nrow(st) && !is.na(st$length[1])) st$length[1] else max(d$stop)

  gaps <- find_n_runs(seq)
  feats <- rbind(
    data.frame(start = d$start, stop = d$stop, type = d$type, element = d$element,
      kind = "rna", trunc = is_truncated(d), stringsAsFactors = FALSE),
    if (nrow(gaps) > 0) data.frame(start = gaps$start, stop = gaps$stop,
      type = "assembly_gap", element = NA_character_, kind = "gap", trunc = FALSE,
      stringsAsFactors = FALSE))
  feats <- feats[order(feats$start), , drop = FALSE]
  rna <- which(feats$kind == "rna")
  first_rna <- if (length(rna)) min(rna) else NA_integer_
  last_rna <- if (length(rna)) max(rna) else NA_integer_

  lines <- hdr
  for (i in seq_len(nrow(feats))) {
    if (identical(feats$kind[i], "gap")) {
      lines <- c(lines,
        paste0(feats$start[i], "\t", feats$stop[i], "\tassembly_gap"),
        paste0("\t\t\testimated_length\t", feats$stop[i] - feats$start[i] + 1),
        "\t\t\tgap_type\twithin scaffold",
        "\t\t\tlinkage_evidence\tunspecified")
      next
    }
    # a truncated gene is partial whichever end is missing, so it gets both marks
    s <- if (feats$trunc[i] || (i == first_rna && (feats$start[i] <= 1 || isTRUE(force_partial)))) {
      paste0("<", feats$start[i])
    } else {
      as.character(feats$start[i])
    }
    e <- if (feats$trunc[i] || (i == last_rna && (feats$stop[i] >= seqlen || isTRUE(force_partial)))) {
      paste0(">", feats$stop[i])
    } else {
      as.character(feats$stop[i])
    }
    lines <- c(lines, paste0(s, "\t", e, "\t", feats$type[i]),
      paste0("\t\t\tproduct\t", rrna_element_product(feats$element[i])))
    if (identical(feats$type[i], "misc_RNA")) lines <- c(lines, paste0("\t\t\tnote\t", feats$element[i]))
  }
  paste(lines, collapse = "\n")
}

#' Portal-parity submission-risk flags for one export sequence
#'
#' NCBI's Submission Portal accepts FASTA only for rRNA-ITS and runs ITSx (stock
#' parameters) on it, silently removing sequences it judges chimeric, gapped
#' between ITS regions, or misassembled
#' (https://submit.ncbi.nlm.nih.gov/help/genbank/#annotation_rRNAITS_faq). This
#' re-checks each sequence the same way so a curator sees the risk at export
#' rather than as a rejection. Warn-only: nothing is dropped from the bundle.
#' @return list(chimeric, gap_in_rdna, misassembly, note); `note` is a
#'   `;`-joined string of the flags that fired ("" when clean).
#' @noRd
submission_risks <- function(con, id, path = 1L, scaffold = 1L) {
  chim <- DBI::dbGetQuery(con,
    "SELECT chimeric FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$chimeric[1]
  chm <- tryCatch(DBI::dbGetQuery(con,
    "SELECT label, reason FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  ann <- DBI::dbGetQuery(con,
    "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  seq <- NA_character_
  if (nrow(ann) > 0) {
    st <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))
    if (nrow(st)) seq <- st$sequence[1]
  }
  risk_flags(chim, ann, seq,
    chimera_label = if (nrow(chm)) chm$label[1] else NA_character_,
    chimera_reason = if (nrow(chm)) chm$reason[1] else NA_character_)
}

#' The risk computation of [submission_risks()] over already-fetched facts
#'
#' `ann` is that unit's annotations ordered by start (or NULL/0 rows), `seq` its
#' assembled sequence (length 1, may be NA). `chimera_label`/`chimera_reason`
#' come from the `chimera` table ([chimera_unit()]); `suspect` and `likely` both
#' raise the flag, warn-only, never blocking export.
#' @noRd
risk_flags <- function(chimeric, ann, seq, chimera_label = NA_character_, chimera_reason = NA_character_) {
  chim <- isTRUE(chimeric[1] == 1L) || isTRUE(chimera_label %in% c("suspect", "likely"))
  gap_in_rdna <- FALSE
  misassembly <- FALSE
  if (!is.null(ann) && nrow(ann) > 0) {
    gaps <- find_n_runs(if (length(seq)) seq[1] else NA_character_)
    if (nrow(gaps) > 0) {
      # an assembly gap inside the transcribed span is NCBI's "gap between ITS regions"
      gap_in_rdna <- any(gaps$stop >= min(ann$start) & gaps$start <= max(ann$stop))
    }
    # features must run strictly 5'->3' in canonical order and not overlap
    idx <- match(ann$element, RRNA_ELEMENT_ORDER)
    misassembly <- any(diff(idx) <= 0, na.rm = TRUE) ||
      (nrow(ann) > 1 && any(ann$start[-1] <= ann$stop[-nrow(ann)]))
  }
  chim_tag <- if (chim) {
    if (!is.na(chimera_reason) && nzchar(chimera_reason)) {
      sprintf("chimeric(%s)", gsub(";", ",", chimera_reason, fixed = TRUE))
    } else "chimeric"
  }
  flags <- c(chim_tag, if (gap_in_rdna) "gap_in_rDNA", if (misassembly) "misassembly")
  list(chimeric = chim, gap_in_rdna = gap_in_rdna, misassembly = misassembly,
    note = paste(flags, collapse = ";"))
}

#' Which annotation rows the curation step flagged truncated (too little of the
#' rRNA model recovered). `d` is annotation rows; FALSE when it has no flags column.
#' @noRd
is_truncated <- function(d) {
  if (is.null(d$flags)) return(rep(FALSE, nrow(d)))
  grepl("(^|,)truncated(,|$)", ifelse(is.na(d$flags), "", d$flags))
}

#' TRUE where a submission_risk string carries a contamination tag
#' (off_target or contaminant), the ones the export window asks about.
#' @noRd
has_screen_risk <- function(risk) {
  grepl("(^|;)(off_target\\(|contaminant(;|$))", ifelse(is.na(risk), "", risk))
}

#' README text mapping the bundle to NCBI's two submission routes
#' @noRd
submission_readme <- function(label, itsx_check = TRUE) {
  paste(c(
    paste0("RiboPilot GenBank submission bundle: ", label),
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    "HOW NCBI HANDLES rRNA-ITS SUBMISSIONS",
    "-------------------------------------",
    "The NCBI Submission Portal (https://submit.ncbi.nlm.nih.gov) accepts FASTA",
    "ONLY for rRNA-ITS sequences. It runs ITSx (stock parameters, no organism",
    "tuning) and writes the annotation itself as a SINGLE misc_RNA feature over",
    "the whole region, e.g. (accession PZ692210):",
    "",
    "  misc_RNA  <1..>N",
    "     /note=\"contains small subunit ribosomal RNA, internal transcribed",
    "            spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2,",
    "            and large subunit ribosomal RNA\"",
    "",
    "Per-element annotation is NOT submitted through the Portal, and the .tbl",
    "feature tables in this bundle are NOT used there.",
    "",
    paste0("  -> Portal route: upload ", label, ".fasta only; NCBI annotates for you."),
    "",
    "The .tbl files are only for NCBI's command-line tool table2asn",
    "(https://www.ncbi.nlm.nih.gov/genbank/table2asn/). Run it with -V v, read",
    "the .val file it writes, and send the .sqn file the way NCBI's GenBank",
    "submission page says.",
    "",
    "WHAT IS IN THIS FOLDER",
    "----------------------",
    paste0("  ", label, ".fasta          every exported record, one FASTA file"),
    paste0("  ", label, ".tbl            every feature table, concatenated"),
    paste0("  ", label, "_summary.csv    one row per exported record"),
    paste0("  ", label, "_report.html    HTML report (feature maps and alignments)"),
    "  samples/             <seqid>.fasta and <seqid>.tbl, one pair per record",
    "  markers/<ELEMENT>/   one FASTA per record per element, plus all_<ELEMENT>.fasta",
    "",
    "SUBMISSION-RISK CHECKS",
    "----------------------",
    "NCBI's ITSx evaluation may silently remove sequences it judges chimeric,",
    "gapped between ITS regions, or misassembled. RiboPilot re-checks each",
    paste0("sequence the same way; see the submission_risk column in ", label, "_summary.csv."),
    "Flagged sequences are still exported (warn-only); review before submitting:",
    "  chimeric     the chimera detector or ITSx flagged mixed origins",
    "  gap_in_rDNA  an assembly-gap (N-run) falls inside the 18S..28S span",
    "  misassembly  features overlap or are out of 5'->3' order",
    "  off_target(<taxon>)  the contamination screen's GenBank hits agree on",
    "               <taxon>, which is not on the sample's lineage",
    "  contaminant  the best GenBank hit shares no phylum with the sample",
    "  taxon_differs(<name>)  your Taxon differs from <name>, the name set in",
    "               the Taxonomy step. The organism column shows which one was",
    "               exported: <name> when your Taxon matched it exactly or you",
    "               picked a species there, otherwise your Taxon, so labels",
    "               like 'Genus sp. A' are kept.",
    "Check off_target and contaminant records before submitting: the sequence",
    "may be another organism's rDNA, exported under this sample's name.",
    "",
    if (itsx_check) c(
    "FINAL ITSx CHECK",
    "----------------",
    "Stock ITSx (default parameters, no organism group, the Portal's setup) was",
    paste0("run over ", label, ".fasta; see the itsx_check column in ", label, "_summary.csv"),
    paste0("and the per-record detail in ", label, "_itsx_check.csv (raw ITSx output in"),
    "itsx_check/). It ran inside the project's container, and the same flags",
    "were written back to the app's ITSx Warnings column. Empty means ITSx saw",
    "nothing to flag; not_run(...) means ITSx was unavailable. Flags, and the",
    "Portal problem each maps to:",
    "  itsx_chimeric             chimera: 5.8S overlaps SSU/LSU, or domains out of order",
    "  itsx_gap_in_ITS           gap: an N-run inside the ITSx-detected ITS span",
    "  itsx_no_5.8S / itsx_partial_5.8S   misassembly: 5.8S missing or truncated",
    "  itsx_5.8S_overlaps_SSU / _LSU      misassembly: conserved domains collide",
    "  itsx_ITS_too_long         misassembly: ITS region longer than ITSx allows",
    "  itsx_no_ITS1 / itsx_no_ITS2        misassembly: a spacer could not be placed",
    "  itsx_boundary_shift(...)  ITSx places an ITS boundary >25 bp from RiboPilot's",
    "  itsx_origin(<group>)      ITSx assigns the record to a group other than the sample's lineage",
    "  itsx_not_detected         ITSx found no ITS; the Portal cannot annotate it",
    "")), collapse = "\n")
}

#' Per-unit export summary data.frame
#'
#' `units` is a data.frame(ID, path, scaffold, seqid). Includes a
#' `submission_risk` column ([submission_risks]) so the CSV and HTML report both
#' surface portal-parity risk flags.
#' @noRd
export_summary_df <- function(con, units) {
  if (nrow(units) == 0) return(NULL)
  n <- nrow(units)
  ids <- unique(as.character(units$ID))
  qs <- paste(rep("?", length(ids)), collapse = ",")
  # one query per fact table for the whole unit set, joined in R
  fetch <- function(sql) tryCatch(DBI::dbGetQuery(con, sprintf(sql, qs), params = as.list(ids)),
    error = function(e) data.frame())
  m <- fetch(
    "SELECT s.ID, an.path, an.scaffold, asm.length, an.n_features,
       an.missing, an.reviewed, an.partial, an.problematic
     FROM samples s
     JOIN annotate an ON an.ID = s.ID
     JOIN assemblies asm ON asm.ID = s.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE s.ID IN (%s)")
  chi <- fetch("SELECT ID, path, scaffold, chimeric FROM annotate WHERE ID IN (%s)")
  chm <- fetch("SELECT ID, path, scaffold, label, reason FROM chimera WHERE ID IN (%s)")
  ann <- fetch("SELECT ID, path, scaffold, element, start, stop FROM annotations
                WHERE ID IN (%s) ORDER BY ID, path, scaffold, start")
  asm <- fetch("SELECT ID, path, scaffold, sequence FROM assemblies WHERE ID IN (%s)")
  scr <- fetch("SELECT ID, path, scaffold, variant_class, variant_group, variant_identity
                FROM screen WHERE ID IN (%s)")
  sv <- fetch("SELECT ID, path, scaffold, screen_verdict, lca_taxon FROM screen WHERE ID IN (%s)")
  nts <- fetch("SELECT ID, path, scaffold, annotate_notes FROM annotate WHERE ID IN (%s)")
  smp <- fetch("SELECT * FROM samples WHERE ID IN (%s)")

  ukey <- unit_key(units$ID, units$path, units$scaffold)
  mi <- if (nrow(m)) match(ukey, unit_key(m$ID, m$path, m$scaffold)) else rep(NA_integer_, n)
  pick <- function(d, nm, na) {
    if (is.null(d[[nm]]) || length(d[[nm]]) == 0) rep(na, n) else d[[nm]][mi]
  }
  anni <- if (nrow(ann)) split(seq_len(nrow(ann)), unit_key(ann$ID, ann$path, ann$scaffold)) else list()
  seqv <- if (nrow(asm)) asm$sequence[match(ukey, unit_key(asm$ID, asm$path, asm$scaffold))] else rep(NA_character_, n)
  chiv <- if (nrow(chi)) chi$chimeric[match(ukey, unit_key(chi$ID, chi$path, chi$scaffold))] else rep(NA_integer_, n)
  cmi <- if (nrow(chm)) match(ukey, unit_key(chm$ID, chm$path, chm$scaffold)) else rep(NA_integer_, n)
  chm_label <- if (nrow(chm)) chm$label[cmi] else rep(NA_character_, n)
  chm_reason <- if (nrow(chm)) chm$reason[cmi] else rep(NA_character_, n)
  v <- unit_variant_facts_df(scr, units)
  svi <- if (nrow(sv)) match(ukey, unit_key(sv$ID, sv$path, sv$scaffold)) else rep(NA_integer_, n)
  nti <- if (nrow(nts)) match(ukey, unit_key(nts$ID, nts$path, nts$scaffold)) else rep(NA_integer_, n)

  si <- match(units$ID, smp$ID)
  taxon <- (smp$Taxon %||% NA)[si]
  sci <- (smp$sci_name %||% NA)[si]
  organism <- export_organism(taxon, sci, (smp$match_type %||% NA)[si], (smp$rank %||% NA)[si])
  differs <- !is.na(sci) & nzchar(sci) & !is.na(taxon) & trimws(taxon) != sci
  risk <- vapply(seq_len(n), function(i) {
    note <- risk_flags(chiv[i], ann[anni[ukey[i]][[1]], , drop = FALSE], seqv[i],
      chimera_label = chm_label[i], chimera_reason = chm_reason[i])$note
    # the sequence may be another organism's rDNA, exported under this sample's name
    off <- isTRUE(sv$screen_verdict[svi[i]] == "off_target")
    notes <- if (is.na(nti[i])) character(0) else strsplit(nts$annotate_notes[nti[i]] %|NA|% "", "; ", fixed = TRUE)[[1]]
    paste(c(if (nzchar(note)) note,
      if (off) sprintf("off_target(%s)", sv$lca_taxon[svi[i]] %|NA|% "other taxon"),
      if ("contaminant" %in% notes) "contaminant",
      if (differs[i]) sprintf("taxon_differs(%s)", sci[i]),
      # NCBI rejects a record with no [organism=]; flag it rather than drop the record
      if (is.na(organism[i]) || !nzchar(trimws(organism[i]))) "missing_organism"),
      collapse = ";")
  }, character(1))

  data.frame(seqid = units$seqid, ID = units$ID, path = units$path, scaffold = units$scaffold,
    organism = organism,
    length = pick(m, "length", NA_integer_),
    n_features = pick(m, "n_features", NA_integer_),
    missing = pick(m, "missing", NA_character_),
    reviewed = pick(m, "reviewed", NA_character_),
    partial = pick(m, "partial", NA_character_),
    problematic = pick(m, "problematic", NA_character_),
    variant_class = v$variant_class, variant_group = v$variant_group,
    variant_identity = v$variant_identity, n_variants = v$n_variants,
    submission_risk = risk,
    stringsAsFactors = FALSE)
}

#' Join key for an assembly unit
#' @noRd
unit_key <- function(id, path, scaffold) {
  paste(as.character(id), as.integer(path), as.integer(scaffold), sep = "\r")
}

#' The intragenomic variant facts for a set of exported units
#'
#' Records that a unit is one of N intragenomic variants detected in its sample
#' and the anchor identity between them. It answers doc 08's open question 7
#' without imposing a submission policy: the fact is recorded, the choice stays
#' the curator's. Deliberately NOT routed through [submission_risks()] -- a
#' paralog is real biology, and warning on export is the behaviour this
#' sub-project exists to remove.
#'
#' `n_variants` counts the sample's COMPARABLE on-target units, so it agrees with
#' what the grouping actually measured rather than with how many units exist.
#'
#' `scr` is the `screen` rows for the units' samples (0 rows when the project has
#' no screen table).
#' @noRd
unit_variant_facts_df <- function(scr, units) {
  n <- nrow(units)
  out <- list(variant_class = rep(NA_character_, n), variant_group = rep(NA_character_, n),
    variant_identity = rep(NA_real_, n), n_variants = rep(NA_integer_, n))
  if (is.null(scr) || nrow(scr) == 0) return(out)
  i <- match(unit_key(units$ID, units$path, units$scaffold),
    unit_key(scr$ID, scr$path, scr$scaffold))
  hit <- !is.na(i)
  if (!any(hit)) return(out)
  out$variant_class[hit] <- as.character(scr$variant_class[i[hit]])
  out$variant_group[hit] <- as.character(scr$variant_group[i[hit]])
  out$variant_identity[hit] <- as.numeric(scr$variant_identity[i[hit]])
  cnt <- tapply(!is.na(scr$variant_class), as.character(scr$ID), sum)
  out$n_variants[hit] <- as.integer(cnt[as.character(units$ID)[hit]])
  out
}

#' Inline SVG feature-map for one unit (for the HTML report)
#' @noRd
feature_map_svg <- function(con, id, path = 1L, scaffold = 1L, seqid = id) {
  d <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  if (nrow(d) == 0) {
    return("<em>no annotations</em>")
  }
  feature_map_svg_text(d)
}

#' @noRd
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Hand-written SVG for the five-box rRNA feature diagram
#'
#' Right-pointing gene arrows on one row, filled from [RRNA_ELEMENT_COLORS],
#' with centred element labels and a bp axis. Written by hand rather than
#' rendered: the HTML report embeds one of these for every exported unit.
#'
#' `d` is that unit's annotations (element, start, stop) ordered by start.
#' @noRd
feature_map_svg_text <- function(d, width = 900, height = 110) {
  x0 <- min(d$start, na.rm = TRUE)
  x1 <- max(d$stop, na.rm = TRUE)
  span <- max(x1 - x0, 1)
  pad <- 26
  inner <- width - 2 * pad
  sx <- function(v) pad + (v - x0) / span * inner
  ytop <- 18
  hbody <- 34
  ybot <- ytop + hbody
  ymid <- ytop + hbody / 2
  axis_y <- ybot + 22

  ticks <- pretty(c(x0, x1), n = 5)
  ticks <- ticks[ticks >= x0 & ticks <= x1]
  grid <- vapply(ticks, function(t) sprintf(
    "<line x1='%.1f' y1='%d' x2='%.1f' y2='%d' stroke='#EBEBEB' stroke-width='1'/>",
    sx(t), ytop - 8, sx(t), ybot + 8), character(1))
  axis <- vapply(ticks, function(t) sprintf(
    "<text x='%.1f' y='%d' text-anchor='middle' font-size='10' fill='#4D4D4D'>%s</text>",
    sx(t), axis_y, format(t, scientific = FALSE, trim = TRUE)), character(1))

  feats <- vapply(seq_len(nrow(d)), function(i) {
    xs <- sx(d$start[i])
    xe <- sx(d$stop[i])
    head <- min(10, max(0, (xe - xs) * 0.35))
    fill <- unname(RRNA_ELEMENT_COLORS[d$element[i]])
    if (is.na(fill)) fill <- "#BBBBBB"
    # gene arrow: body, then a head standing slightly proud of it
    pts <- sprintf("%.1f,%d %.1f,%d %.1f,%d %.1f,%.1f %.1f,%d %.1f,%d %.1f,%d",
      xs, ytop, xe - head, ytop, xe - head, ytop - 2, xe, ymid,
      xe - head, ybot + 2, xe - head, ybot, xs, ybot)
    lab <- xml_escape(d$element[i])
    paste0(
      sprintf("<polygon points='%s' fill='%s' stroke='#333333' stroke-width='1'/>", pts, fill),
      sprintf("<text x='%.1f' y='%.1f' text-anchor='middle' dominant-baseline='central' font-size='11' fill='#000000'>%s</text>",
        (xs + xe - head) / 2, ymid, lab))
  }, character(1))

  paste0(sprintf("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 %d %d' width='%d' height='%d'>",
    width, height, width, height),
    paste(grid, collapse = ""), paste(feats, collapse = ""), paste(axis, collapse = ""),
    "</svg>")
}

#' Write a self-contained HTML report for an export group / selection
#'
#' `units` is a data.frame(ID, path, scaffold, seqid) - one row per exported unit.
#' `smry` is the already-computed [export_summary_df()] for those units.
#' @noRd
export_report_html <- function(con, units, out_html, label, smry = export_summary_df(con, units)) {
  hdr <- c("seqid", "ID", "organism", "length", "n_features", "missing", "reviewed", "submission_risk", "itsx_check")
  rows <- apply(smry[, intersect(hdr, names(smry)), drop = FALSE], 1, function(r) {
    paste0("<tr>", paste0("<td>", ifelse(is.na(r), "", r), "</td>", collapse = ""), "</tr>")
  })
  maps <- vapply(seq_len(nrow(units)), function(i) {
    u <- units[i, ]
    paste0("<div class='sample'><h3>", u$seqid, " <small>",
      (smry$organism[i] %||% NA) %|NA|% "", "</small></h3>",
      feature_map_svg(con, u$ID, u$path, u$scaffold, u$seqid), "</div>")
  }, character(1))
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>RiboPilot export: ", label,
    "</title><style>body{font-family:system-ui,Arial,sans-serif;margin:2em;color:#222;}",
    "table{border-collapse:collapse;margin:1em 0;}th,td{border:1px solid #ddd;padding:4px 8px;font-size:0.9em;}",
    "th{background:#f4f6f8;text-align:left;}.sample{margin:1.2em 0;border-top:1px solid #eee;padding-top:0.5em;}",
    "small{color:#888;font-weight:normal;}svg{max-width:100%;height:auto;}</style></head><body>",
    "<h1>RiboPilot export report</h1><p><strong>Group:</strong> ", label,
    " &middot; <strong>Units:</strong> ", nrow(units),
    " &middot; ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>",
    "<h2>Summary</h2><table><thead><tr>",
    paste0("<th>", intersect(hdr, names(smry)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table>",
    "<h2>Feature models</h2>", paste(maps, collapse = ""),
    "</body></html>")
  writeLines(html, out_html)
}

#' Build a GenBank submission bundle for a set of samples (one export group)
#'
#' @param path Project directory.
#' @param ids Sample IDs to export (default: every exportable unit). If `group`
#'   is given, restricted to that group. Sample-grain: every exportable unit of
#'   a listed sample is included. Use `units` to pick individual scaffolds.
#' @param units Optional data.frame(ID, path, scaffold) restricting the export to
#'   those exact assembly units. Applied on top of `ids` / `group`.
#' @param group Export-group name to filter by and label the bundle.
#' @param out_dir Project "out" directory (default `<path>/out`); files land in
#'   `<out_dir>/export/<label>`, mirroring MitoPilot.
#' @param label Bundle label (default: `group` or "export").
#' @param markers,report Also write per-marker fastas / an HTML report.
#' @param template Optional path to a GenBank submission template (`.sbt`);
#'   copied into the bundle as `<label>.sbt` so table2asn picks it up.
#' @param run_table2asn Run table2asn if on PATH.
#' @param zip Also write `<out_dir>/<label>.zip`.
#' @param export_opts Name of the saved FASTA-header configuration to use.
#' @param fasta_header,fasta_header_marker Optional header templates
#'   (\code{\{token\}} glue syntax) overriding the saved configuration.
#' @param itsx_check Run stock ITSx (default parameters, no organism group) over
#'   the exported records, the way NCBI's Submission Portal does, and report
#'   what it would flag in the `itsx_check` column of the summary CSV and in
#'   `<label>_itsx_check.csv`. Warn-only. ITSx runs inside the project's
#'   container (docker / singularity per `.config`); a conda-executor project
#'   uses the ITSx on PATH. Otherwise the column reads `not_run(...)`.
#' @param itsx_cpu Threads for the ITSx check.
#' @return (invisibly) the group output directory, with the summary data.frame
#'   attached as attribute `summary`.
#' @export
export_bundle <- function(path = ".", ids = NULL, units = NULL, group = NULL,
                          out_dir = file.path(path, "out"),
                          label = NULL, markers = TRUE, report = TRUE, template = NULL,
                          run_table2asn = FALSE, zip = TRUE,
                          export_opts = "default", fasta_header = NULL,
                          fasta_header_marker = NULL, itsx_check = TRUE,
                          itsx_cpu = 2L) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Export annotation-locked, non-ignored, non-problematic UNITS (one GenBank
  # record per (ID, path, scaffold)). Filtered by group / ids / units when given.
  all_units <- DBI::dbGetQuery(con,
    "SELECT an.ID, an.path, an.scaffold, an.export_group
     FROM annotate an
     JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
     WHERE an.annotate_lock = 1 AND COALESCE(asm.ignore,0) = 0
       AND COALESCE(an.problematic,'') != 'yes'
     ORDER BY an.ID, an.path, an.scaffold")
  sel <- all_units
  if (!is.null(group)) sel <- sel[!is.na(sel$export_group) & sel$export_group == group, , drop = FALSE]
  if (!is.null(ids)) sel <- sel[sel$ID %in% ids, , drop = FALSE]
  if (!is.null(units) && nrow(units) > 0) {
    key <- function(d) paste(d$ID, as.integer(d$path), as.integer(d$scaffold))
    sel <- sel[key(sel) %in% key(units), , drop = FALSE]
  }
  if (nrow(sel) == 0) {
    stop("No annotation-locked units to export.", call. = FALSE)
  }
  # SeqID: the sample's ID when it exports a single unit, else ID_p<path>_s<scaffold>.
  # Counted over ALL exportable units of the sample, not over this bundle's
  # subset, so a record keeps the same SeqID across a group export and a
  # whole-project export. Renaming a record after an accession is bound to it
  # cannot be undone.
  n_per_id <- table(all_units$ID)
  units <- sel
  n_units <- as.integer(n_per_id[units$ID])
  units$seqid <- ifelse(!is.na(n_units) & n_units > 1,
    sprintf("%s_p%s_s%s", units$ID, units$path, units$scaffold), units$ID)
  label <- label %||% group %||% "export"

  # Headers are resolved and checked before anything is written: a bad template
  # or a header that repeats across records is a rejected submission, not a
  # cosmetic problem, and a half-written bundle would hide it.
  opts <- get_export_opts(con, export_opts)
  tpl_record <- fasta_header %||% opts$fasta_header
  tpl_marker <- fasta_header_marker %||% opts$fasta_header_marker
  hdata <- export_header_data(con, units)
  marker_data <- function(i, element) {
    r <- hdata[i, , drop = FALSE]
    r$element <- element
    r$product <- if (identical(element, "ITS")) {
      paste("internal transcribed spacer 1, 5.8S ribosomal RNA gene,",
        "and internal transcribed spacer 2")
    } else {
      rrna_element_product(element)
    }
    r
  }
  first_el <- vapply(strsplit(hdata$elements %||% "", "-"),
    function(x) if (length(x)) x[1] else "18S", character(1))
  for (chk in list(
    list(what = "FASTA header", res = validate_fasta_header(tpl_record, hdata)),
    list(what = "marker FASTA header",
      res = validate_fasta_header(tpl_marker,
        do.call(rbind, lapply(seq_len(nrow(hdata)), function(i) marker_data(i, first_el[i]))))))) {
    if (identical(chk$res$level, "error")) {
      stop("The ", chk$what, " template cannot be used: ", chk$res$message, call. = FALSE)
    }
    if (identical(chk$res$level, "warn")) {
      warning("The ", chk$what, " template: ", chk$res$message, call. = FALSE)
    }
  }
  headers <- render_fasta_headers(tpl_record, hdata)
  dup <- unique(headers[duplicated(headers)])
  if (length(dup) > 0) {
    stop("The FASTA header template gives ", length(dup),
      " duplicated header(s); add {seqid} so every record is named once: ",
      paste(dup, collapse = "; "), call. = FALSE)
  }

  gdir <- file.path(out_dir, "export", label)
  unlink(gdir, recursive = TRUE)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  sdir <- file.path(gdir, "samples")
  dir.create(sdir, showWarnings = FALSE)
  if (markers) dir.create(file.path(gdir, "markers"), showWarnings = FALSE)

  group_fsa <- character(0)
  group_tbl <- character(0)
  marker_pool <- list()

  for (ui in seq_len(nrow(units))) {
    id <- units$ID[ui]; upath <- units$path[ui]; uscaf <- units$scaffold[ui]; seqid <- units$seqid[ui]
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, upath, uscaf))$sequence[1]
    defline <- paste0(">", headers[ui])
    writeLines(c(defline, seq), file.path(sdir, paste0(seqid, ".fasta")))
    is_partial <- identical(DBI::dbGetQuery(con,
      "SELECT partial FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, upath, uscaf))$partial[1], "yes")
    tbl <- feature_table_text_con(con, id, upath, uscaf, seqid = seqid, force_partial = is_partial)
    writeLines(tbl, file.path(sdir, paste0(seqid, ".tbl")))
    group_fsa <- c(group_fsa, defline, seq)
    group_tbl <- c(group_tbl, tbl)

    if (markers) {
      ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
        params = list(id, upath, uscaf))
      wm <- function(el, s, e) {
        if (is.na(s) || is.na(e) || e < s) {
          return(invisible())
        }
        el_clean <- gsub("[^A-Za-z0-9]", "", el)
        el_dir <- file.path(gdir, "markers", el_clean)
        dir.create(el_dir, recursive = TRUE, showWarnings = FALSE)
        rec <- c(paste0(">", render_fasta_headers(tpl_marker, marker_data(ui, el))),
          substr(seq, s, e))
        writeLines(rec, file.path(el_dir, paste0(seqid, "_", el_clean, ".fasta")))
        marker_pool[[el]] <<- c(marker_pool[[el]], rec)
      }
      for (i in seq_len(nrow(ann))) wm(ann$element[i], ann$start[i], ann$stop[i])
      i1 <- ann[ann$element == "ITS1", ]
      i2 <- ann[ann$element == "ITS2", ]
      if (nrow(i1) && nrow(i2)) wm("ITS", i1$start[1], i2$stop[1])
    }
  }

  # combined group files
  writeLines(group_fsa, file.path(gdir, paste0(label, ".fasta")))
  writeLines(paste(group_tbl, collapse = "\n"), file.path(gdir, paste0(label, ".tbl")))
  if (markers) {
    for (el in names(marker_pool)) {
      el_clean <- gsub("[^A-Za-z0-9]", "", el)
      writeLines(marker_pool[[el]],
        file.path(gdir, "markers", el_clean, paste0("all_", el_clean, ".fasta")))
    }
  }
  # per-unit summary CSV (+ warn-only portal-parity submission-risk gate)
  smry <- export_summary_df(con, units)
  risky <- which(nzchar(smry$submission_risk))
  if (length(risky) > 0) {
    warning("Submission-risk flags on ", length(risky), " of ", nrow(smry),
      " unit(s) (still exported; see submission_risk column): ",
      paste(sprintf("%s[%s]", smry$seqid[risky], smry$submission_risk[risky]), collapse = ", "),
      call. = FALSE)
  }
  # final ITSx check: stock run, as the Portal does it (warn-only)
  if (isTRUE(itsx_check)) {
    chk <- itsx_export_check(con, units, gdir, label, cpu = itsx_cpu, project = path)
    smry$itsx_check <- chk$itsx_check[match(smry$seqid, chk$seqid)]
    not_run <- startsWith(smry$itsx_check, "not_run(")
    if (any(not_run)) {
      message("ITSx check skipped: ", sub("^not_run\\((.*)\\)$", "\\1", smry$itsx_check[not_run][1]))
    }
    flagged <- which(nzchar(smry$itsx_check) & !not_run)
    if (length(flagged) > 0) {
      warning("ITSx check flags on ", length(flagged), " of ", nrow(smry),
        " unit(s) (still exported; see itsx_check column): ",
        paste(sprintf("%s[%s]", smry$seqid[flagged], smry$itsx_check[flagged]), collapse = ", "),
        call. = FALSE)
    }
  }
  utils::write.csv(smry, file.path(gdir, paste0(label, "_summary.csv")), row.names = FALSE)
  writeLines(submission_readme(label, itsx_check = isTRUE(itsx_check)), file.path(gdir, "README.txt"))
  # HTML report
  if (isTRUE(report)) {
    tryCatch(export_report(con, units,
      out_file = file.path(gdir, paste0(label, "_report.html")),
      label = label, alignments = TRUE),
      error = function(e) message("Report generation skipped: ", conditionMessage(e)))
  }
  # submission template (.sbt): copy in so table2asn (-t) picks it up
  if (!is.null(template)) {
    if (!file.exists(template)) {
      warning("Submission template not found: ", template)
    } else {
      file.copy(template, file.path(gdir, paste0(label, ".sbt")), overwrite = TRUE)
    }
  }

  # optional table2asn, per file: its -indir mode only picks up .fsa
  if (isTRUE(run_table2asn) && nzchar(unname(Sys.which("table2asn")))) {
    sbt <- Sys.glob(file.path(gdir, "*.sbt"))
    for (f in Sys.glob(file.path(sdir, "*.fasta"))) {
      args <- c("-i", shQuote(f), "-a", "s")
      if (length(sbt) > 0) args <- c(args, "-t", shQuote(sbt[1]))
      tryCatch(system2("table2asn", args, stdout = TRUE, stderr = TRUE),
        error = function(e) message("table2asn skipped for ", basename(f), ": ",
          conditionMessage(e)))
    }
  }

  message("Exported ", nrow(units), " unit(s) to ", gdir,
    "\n  GenBank Portal takes FASTA only: upload ", label,
    ".fasta; NCBI auto-annotates the rRNA-ITS region (see README.txt).")
  if (isTRUE(zip)) {
    zf <- paste0(gdir, ".zip")
    old <- setwd(dirname(gdir))
    on.exit(setwd(old), add = TRUE)
    utils::zip(zf, label, flags = "-rq")
  }
  attr(gdir, "summary") <- smry
  invisible(gdir)
}

#' Export every assigned export group in a project
#' @param path Project directory.
#' @param out_dir Project "out" directory; bundles land in `<out_dir>/export`.
#' @return (invisibly) the export groups exported.
#' @export
export_project <- function(path = ".", out_dir = file.path(path, "out")) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  groups <- DBI::dbGetQuery(con, "SELECT DISTINCT export_group FROM annotate WHERE export_group IS NOT NULL")$export_group
  DBI::dbDisconnect(con)
  if (length(groups) == 0) {
    stop("No export groups assigned. Assign samples to a group first.", call. = FALSE)
  }
  for (g in groups) export_bundle(path, group = g, out_dir = out_dir, zip = TRUE)
  invisible(groups)
}

#' Plain FASTA for a set of assembly units, for the app's copy-to-clipboard
#'
#' Deliberately not the submission defline built by [export_bundle()]: this is
#' for pasting into BLAST or an alignment, so the header is the three things
#' that identify the sequence to a person -- sample, the taxon the user gave,
#' and the unit. `samples.Taxon` is the user's own label; `sci_name` is what
#' taxonomy resolution assigned and is not used here.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param units data.frame with `path` and `scaffold` columns, in output order.
#' @param width Sequence line width.
#' @return Character vector of FASTA lines (empty when nothing to write).
#' @noRd
unit_fasta_text <- function(con, id, units, width = 60) {
  if (is.null(units) || nrow(units) == 0) return(character(0))
  taxon <- tryCatch(DBI::dbGetQuery(con, "SELECT Taxon FROM samples WHERE ID = ?",
    params = list(id))$Taxon[1], error = function(e) NA_character_)
  out <- character(0)
  for (i in seq_len(nrow(units))) {
    p <- as.integer(units$path[i]); s <- as.integer(units$scaffold[i])
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))$sequence[1]
    if (length(seq) == 0 || is.na(seq) || !nzchar(seq)) next
    hdr <- paste(c(id, if (!is.na(taxon) && nzchar(taxon)) taxon else NULL,
      sprintf("p%d.s%d", p, s)), collapse = " ")
    starts <- seq.int(1L, nchar(seq), by = width)
    out <- c(out, paste0(">", hdr), substring(seq, starts, starts + width - 1L))
  }
  out
}
