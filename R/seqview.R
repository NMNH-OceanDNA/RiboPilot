# Sequence viewer (inst/app/www/seqviewer.js), ported from MitoPilot.

#' Sequence viewer payload
#'
#' Pure: no Shiny, no database. Features are keyed by element name, which is
#' unique within a unit.
#' @param seq Unit sequence. @param feats Annotations (element, start, stop,
#'   strand, partial_start, partial_stop). @param depths,clips Per-position
#'   read depth and soft-clip counts, or NULL. @param chim [chimera_origin()]
#'   keys to add, or NULL.
#' @noRd
seqview_payload <- function(seq, feats, depths = NULL, clips = NULL, selected = NULL, chim = NULL) {
  seq <- toupper(as.character(seq)[1] %|NA|% "")
  n <- nchar(seq)
  ok <- which(!is.na(feats$start) & !is.na(feats$stop))
  out <- list(len = n, seq = seq,
    cutoff = CHIMERA_CLIP_FRAC,
    features = lapply(ok, function(i) list(
      row = feats$element[i], gene = feats$element[i],
      pos1 = as.integer(feats$start[i]), pos2 = as.integer(feats$stop[i]),
      dir = if (identical(feats$strand[i], "-")) "-" else "+",
      partial5 = isTRUE(feats$partial_start[i] == 1), partial3 = isTRUE(feats$partial_stop[i] == 1),
      color = unname(RRNA_ELEMENT_COLORS[feats$element[i]] %|NA|% "#888888"))))
  if (length(selected) == 1L) out$selected <- selected
  if (length(depths) && n) {
    d <- as.numeric(depths)[seq_len(n)]
    out$depth <- as.integer(d)
    if (length(clips)) {
      # low-depth fractions are noise; the chimera check ignores them too
      f <- round(as.numeric(clips)[seq_len(n)] / d, 3)
      f[is.na(d) | d < 10] <- NA
      if (any(!is.na(f))) out$clip <- f
    }
  }
  out[names(chim)] <- chim
  out
}

#' Sequence viewer section
#' @param prefix Id prefix: the canvas is `<prefix>-canvas`.
#' @param summary Section heading. @param open Start expanded?
#' @param pick Does clicking a gene select it?
#' @noRd
seqview_ui <- function(prefix, summary, open = TRUE, pick = TRUE) {
  id <- function(x) paste0(prefix, "-", x)
  btn <- function(action, label, icon = NULL, title = NULL) {
    tags$button(type = "button", class = "btn btn-default btn-sm", `data-rpseq` = action,
      title = title, `aria-label` = title %||% label, icon, label)
  }
  goto <- tags$input(id = id("goto"), type = "number", min = 1, step = 1, class = "form-control input-sm",
    style = "width:130px;", placeholder = "Go to position", `aria-label` = "Go to position")
  box <- function(x, label) tags$label(class = "rp-seqview-check",
    tags$input(type = "checkbox", id = id(x), checked = NA), label)
  rp_details(summary, open = open,
    div(class = "rp-seqview-controls",
      goto,
      btn("fit", "Fit gene", title = "Zoom to the selected gene"),
      btn("whole", "Whole sequence"),
      btn("zoom_in", NULL, icon("magnifying-glass-plus"), "Zoom in"),
      btn("zoom_out", NULL, icon("magnifying-glass-minus"), "Zoom out"),
      btn("copy", "Copy sequence", title = "Copy the bases in view as FASTA"),
      box("show_cov", "Coverage"), box("show_clip", "Clipped reads"),
      box("show_ident", "Window identity"), box("show_origin", "Origin"), box("show_nt", "Nucleotides")),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      paste0("Drag or scroll sideways to pan, scroll to zoom", if (pick) "; click a gene to select it",
        ". Letters appear when zoomed in. Copy sequence copies the bases in view. "),
      sprintf("Clipped reads is the share of reads soft-clipped at each base (shown where depth is at least 10); bars turn orange at %d%%, the chimera check's cutoff.",
        round(100 * CHIMERA_CLIP_FRAC)),
      " Origin colours each stretch by the class its BLAST windows match: grey is your lineage, white is no hit, ",
      "other colours are foreign; a dashed line marks a possible junction. Window identity is the best hit's identity. ",
      "Hover for the match."),
    div(class = "rp-seqview",
      tags$canvas(id = id("canvas"), class = "rp-seqview-canvas"),
      tags$div(id = id("tip"), class = "rp-seqview-tip", hidden = NA)))
}

#' Send one unit to a viewer
#' @param canvas Canvas id. @param input Input id for gene clicks, or NULL.
#' @param review TRUE (Annotation Review): the review layer, with the copy's
#'   coverage window and trimmed ends. FALSE (Assemble step): the originals.
#' @noRd
seqview_send <- function(session, canvas, input, id, path, scaffold, selected = NULL, review = FALSE) {
  con <- app_con(); on.exit(DBI::dbDisconnect(con))
  q <- function(sql, tail = "") DBI::dbGetQuery(con,
    paste(sql, "WHERE ID = ? AND path = ? AND scaffold = ?", tail), params = list(id, path, scaffold))
  s <- q(paste("SELECT sequence FROM", if (review) "review_units_v" else "assemblies"))
  if (!nrow(s) || is.na(s$sequence[1])) return(invisible())
  f <- q(paste("SELECT element, start, stop, strand,", if (review) REVIEW_EFF_PARTIALS else "partial_start, partial_stop",
    "FROM", if (review) "review_annotations_v" else "annotations"), "ORDER BY start")
  if (nrow(f)) {
    e <- partial_ends(f, nchar(s$sequence[1]), if (review) unit_trimmed(con, id, path, scaffold) else c(FALSE, FALSE))
    f$partial_start <- as.integer(e$p5); f$partial_stop <- as.integer(e$p3)
  }
  tr <- if (review) review_coverage(con, id, path, scaffold) else {
    cov <- q("SELECT * FROM coverage")
    list(depths = if (nrow(cov)) parse_depths(cov$depths[1]),
      clips = if (nrow(cov) && !is.null(cov$clips)) parse_depths(cov$clips[1]))
  }
  chim <- seqview_chimera(con, id, path, scaffold, review, nchar(s$sequence[1]))
  pl <- seqview_payload(s$sequence[1], f, tr$depths, tr$clips, selected, chim)
  pl$id <- canvas; pl$input <- input; pl$unit <- sprintf("%s_p%s_s%s", id, path, scaffold)
  session$sendCustomMessage("rpseq", pl)
}

#' Chimera evidence for the viewer: de-overlapped origin segments, junctions,
#' and per-base window identity. Pure.
#' @noRd
chimera_origin <- function(seg, win, junctions, label, channels, sample_lt, len) {
  out <- list(origin = list(), junctions = list(), ident = NULL)
  if (NROW(seg)) {
    seg <- seg[order(seg$start), , drop = FALSE]
    # stored JSON drops fields that were NA in every row
    for (col in c("species", "pident")) if (is.null(seg[[col]])) seg[[col]] <- NA
    seg$col <- chimera_seg_colors(seg, sample_lt)
    stop <- pmin(seg$stop, c(seg$start[-1] - 1L, Inf))
    out$origin <- lapply(seq_len(nrow(seg)), function(i) list(
      pos1 = as.integer(seg$start[i]), pos2 = as.integer(stop[i]), color = unname(seg$col[i]),
      taxon = seg$taxon[i], species = seg$species[i] %|NA|% "", pident = seg$pident[i],
      host = identical(unname(seg$col[i]), "#B3B3B3")))
  }
  if (NROW(win) && any(!is.na(win$pident))) {
    id <- rep(NA_real_, len)
    for (i in which(!is.na(win$pident))) {
      r <- max(1L, win$start[i]):min(len, win$stop[i])
      id[r] <- pmax(id[r], win$pident[i], na.rm = TRUE)
    }
    out$ident <- id
  }
  j <- chimera_junction_table(junctions, channels)
  if (nrow(j)) out$junctions <- lapply(seq_len(nrow(j)), function(i) list(
    pos = as.integer(j$position[i]), severe = identical(label, "severe"), checks = j$checks[i]))
  out
}

#' Chimera evidence of one unit's layer; empty before the Annotation Review
#' checks are done or without a chimera row
#' @noRd
seqview_chimera <- function(con, id, path, scaffold, review, len) {
  if (review && !identical(review_unit_state(con,
    data.frame(ID = id, path = as.integer(path), scaffold = as.integer(scaffold)))$state, "done")) return(list())
  tbl <- function(x) if (review) paste0("review_", x, "_v") else x
  q <- function(sql) DBI::dbGetQuery(con, paste(sql, "WHERE ID = ? AND path = ? AND scaffold = ?"),
    params = list(id, path, scaffold))
  ch <- q(paste("SELECT label, junctions, segments, channels FROM", tbl("chimera")))
  if (!nrow(ch)) return(list())
  win <- q(paste("SELECT start, stop, pident FROM", tbl("chimera_windows")))
  seg <- if (!is.na(ch$segments[1])) jsonlite::fromJSON(ch$segments[1]) else data.frame()
  samp <- DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))
  o <- chimera_origin(seg, win, ch$junctions[1], ch$label[1], ch$channels[1],
    if (nrow(samp)) samp$lineage_taxids[1] else NA_character_, len)
  Filter(length, o)
}

#' Zoom one viewer to a range
#' @noRd
seqview_goto <- function(session, canvas, pos1, pos2)
  session$sendCustomMessage("rpseq_goto", list(id = canvas, pos1 = as.integer(pos1), pos2 = as.integer(pos2)))
