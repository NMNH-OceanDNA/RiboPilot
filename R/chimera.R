# Chimera detection: per-gene taxonomy, windowed local BLAST, ribotyper,
# read evidence. Design: dev/design/21_CHIMERA_DETECTION.md.

CHIMERA_WIN <- 150L
CHIMERA_WIN_STEP <- 150L
# a class-level switch needs this many windows on each side; phylum or deeper needs 2
CHIMERA_CLASS_MIN_WINDOWS <- 3L
# read-evidence step (depth medians, clip pileup margin)
CHIMERA_STEP <- 250L
CHIMERA_PID_DROP <- 90
CHIMERA_DEPTH_STEP <- 3
CHIMERA_CLIP_FRAC <- 0.4
CHIMERA_JUNCTION_TOL <- 300L
CHIMERA_MIN_DEPTH <- 20
CHIMERA_REPEAT_K <- 50L
CHIMERA_TIE_FRAC <- 0.02
# window searches keep many hits: genome annotations list hundreds of identical
# rRNA copies per species, and a 50-hit list fills with one of them
CHIMERA_WIN_MAX_TARGETS <- 5000L

#' @noRd
chimera_tile <- function(seq_len, win = CHIMERA_WIN, step = CHIMERA_WIN_STEP) {
  seq_len <- as.integer(seq_len)
  if (seq_len <= win) return(data.frame(start = 1L, stop = seq_len))
  starts <- seq.int(1L, seq_len - win + 1L, by = step)
  stops <- starts + win - 1L
  tail <- seq_len - stops[length(stops)]
  # leftover tail: overlapping tiles get a right-anchored window; end-to-end
  # tiles get the tail as its own window, or join it to the last when under 50 bp
  if (tail > 0L && step < win) { starts <- c(starts, seq_len - win + 1L); stops <- c(stops, seq_len) }
  else if (tail >= 50L) { starts <- c(starts, stops[length(stops)] + 1L); stops <- c(stops, seq_len) }
  else if (tail > 0L) stops[length(stops)] <- seq_len
  data.frame(start = as.integer(starts), stop = as.integer(stops))
}

#' @noRd
chimera_windows_fasta <- function(seq, tiles, file) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  subs <- substring(seq, tiles$start, tiles$stop)
  writeLines(as.vector(rbind(paste0(">", ids), subs)), file)
  invisible(ids)
}

#' Bitscore-weighted LCA per annotated element from stored HSPs
#' @noRd
chimera_gene_lca <- function(hits, feats, min_cov = 0.5) {
  out <- lapply(seq_len(nrow(feats)), function(i) {
    s <- feats$start[i]; e <- feats$stop[i]; len <- e - s + 1
    ov <- pmin(hits$qend, e) - pmax(hits$qstart, s) + 1
    h <- hits[!is.na(ov) & ov >= min_cov * len & !is.na(hits$lineage_taxids), , drop = FALSE]
    h <- h[order(-h$bitscore), , drop = FALSE]
    h <- h[!duplicated(h$accession), , drop = FALSE]
    base <- data.frame(element = feats$element[i], taxon = NA_character_, rank = NA_character_,
      taxid = NA_integer_, species = NA_character_, pident = NA_real_, n_hits = nrow(h),
      lineage_taxids = NA_character_, stringsAsFactors = FALSE)
    if (nrow(h) < 3) return(base)
    lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
    if (is.null(lca) || is.na(lca$taxid %||% NA)) return(base)
    base$taxon <- lca$taxon; base$rank <- lca$rank; base$taxid <- as.integer(lca$taxid)
    base$species <- h$species[1]; base$pident <- mean(h$pident, na.rm = TRUE)
    base$lineage_taxids <- chimera_cut_lineage(h$lineage_taxids[1], lca$rank)
    base
  })
  do.call(rbind, out)
}

#' Blank the fields of a ';'-joined 7-rank lineage below `rank` ("none" blanks all)
#' @noRd
chimera_cut_lineage <- function(x, rank) {
  if (is.na(x)) return(x)
  p <- strsplit(x, ";", fixed = TRUE)[[1]]
  k <- match(rank, RRNA_LINEAGE_RANKS, nomatch = 0L)
  if (k < length(p)) p[(k + 1L):length(p)] <- ""
  paste(p, collapse = ";")
}

#' Shallowest rank (1 = kingdom .. 7 = species) where both lineages resolve and
#' differ; NA when one is an ancestor of the other
#' @noRd
chimera_split_rank <- function(a, b) {
  x <- parse_rank_taxids(a); y <- parse_rank_taxids(b)
  k <- which(!is.na(x) & !is.na(y) & x != y)
  if (length(k)) min(k) else NA_integer_
}

#' @noRd
chimera_unavailable <- function(availability = "unavailable") {
  list(score = 0, junctions = integer(0), reason = NA_character_, availability = availability)
}

#' Channel 1: two resolved elements whose LCAs sit in different classes
#' @noRd
chimera_channel_genes <- function(gene_taxa, feats = NULL) {
  g <- gene_taxa[!is.na(gene_taxa$taxon), , drop = FALSE]
  if (nrow(g) < 2) return(chimera_unavailable("uninformative"))
  worst <- NULL
  shallow <- NULL
  for (i in seq_len(nrow(g) - 1)) for (j in (i + 1):nrow(g)) {
    k <- chimera_split_rank(g$lineage_taxids[i], g$lineage_taxids[j])
    sv <- isTRUE(k <= 2L && min(g$pident[i], g$pident[j]) >= CHIMERA_PID_DROP)
    if (!is.na(k) && k <= 3L && (is.null(worst) || sv > worst$sv || (sv == worst$sv && k < worst$k)))
      worst <- list(i = i, j = j, k = k, sv = sv)
    if (is.na(k) && is.null(shallow)) {
      x <- parse_rank_taxids(g$lineage_taxids[i]); y <- parse_rank_taxids(g$lineage_taxids[j])
      shared <- deepest_shared_rank(x, y)
      if (shared <= 2L && !identical(x, y)) shallow <- list(i = i, j = j, shared = shared)
    }
  }
  if (is.null(worst) && !is.null(shallow)) {
    a <- g[shallow$i, ]; b <- g[shallow$j, ]
    return(list(score = 0.3, junctions = integer(0),
      reason = sprintf("%s resolves to %s (%s, %.1f%%) and %s to %s (%s, %.1f%%): agree only to %s",
        a$element, a$taxon, a$species, a$pident, b$element, b$taxon, b$species, b$pident,
        RRNA_LINEAGE_RANKS[max(shallow$shared, 1L)]),
      availability = "available", severe = FALSE))
  }
  if (is.null(worst)) return(list(score = 0, junctions = integer(0),
    reason = sprintf("elements agree: %s", paste(sprintf("%s %s", g$element, g$taxon), collapse = ", ")),
    availability = "available"))
  a <- g[worst$i, ]; b <- g[worst$j, ]
  junction <- NA_integer_
  if (!is.null(feats)) {
    fa <- feats[feats$element == a$element, ]; fb <- feats[feats$element == b$element, ]
    if (nrow(fa) && nrow(fb)) junction <- as.integer(round((min(fa$stop, fb$stop) + max(fa$start, fb$start)) / 2))
  }
  list(score = 0.6, junctions = junction[!is.na(junction)],
    reason = sprintf("%s resolves to %s (%s, %.1f%%) but %s to %s (%s, %.1f%%): different %s",
      a$element, a$taxon, a$species, a$pident, b$element, b$taxon, b$species, b$pident,
      RRNA_LINEAGE_RANKS[worst$k]),
    availability = "available", severe = worst$sv)
}

#' Best hit per tile from a parse_blast_hits_file() frame
#' @noRd
chimera_window_hits <- function(hits, tiles, sample_taxids = NA_character_,
                                taxdump_dir = NULL, db_source = NA_character_, cap = CHIMERA_WIN_MAX_TARGETS) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  if (!nrow(hits)) {
    return(data.frame(start = tiles$start, stop = tiles$stop,
      accession = NA_character_, species = NA_character_,
      staxid = NA_integer_, lineage = NA_character_,
      lineage_taxids = NA_character_, rank_distance = NA_integer_,
      hit_class = NA_character_, pident = NA_real_, bitscore = NA_real_,
      db_source = db_source, stringsAsFactors = FALSE))
  }
  # one record per taxon, so a tie counts taxa rather than copies
  n_raw <- table(hits$qseqid)
  hits <- hits[order(hits$qseqid, -hits$bitscore), , drop = FALSE]
  hits <- hits[!duplicated(hits[c("qseqid", "staxid")]), , drop = FALSE]
  n_taxa <- table(hits$qseqid)
  lin <- tryCatch(tk_ranked_lineage(hits$staxid, taxdump_dir), error = function(e) data.frame())
  li <- if (nrow(lin)) match(hits$staxid, lin$taxid) else rep(NA_integer_, nrow(hits))
  hits$lineage <- if (nrow(lin)) lin$lineage[li] else NA_character_
  hits$lineage_taxids <- if (nrow(lin)) lin$lineage_taxids[li] else NA_character_
  best <- hits[!duplicated(hits$qseqid), , drop = FALSE]
  tied <- integer(nrow(best))
  for (k in seq_len(nrow(best))) {
    h <- hits[hits$qseqid == best$qseqid[k] & !is.na(hits$lineage_taxids), , drop = FALSE]
    if (nrow(h) < 2) next
    lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, frac = CHIMERA_TIE_FRAC, min_hits = 1L)
    if (lca$n_hits < 2 || identical(lca$rank, "species")) next
    best$lineage[k] <- chimera_cut_lineage(best$lineage[k], lca$rank)
    best$lineage_taxids[k] <- chimera_cut_lineage(best$lineage_taxids[k], lca$rank)
    tied[k] <- lca$n_hits - 1L
  }
  i <- match(ids, best$qseqid)
  lt <- best$lineage_taxids[i]
  rd <- if (!is.na(sample_taxids) && nzchar(sample_taxids)) {
    vapply(lt, function(x) if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, sample_taxids)),
      integer(1), USE.NAMES = FALSE)
  } else rep(NA_integer_, length(i))
  sp <- extract_species(best$title[i])
  # a tie names the clade the tied hits share, not one arbitrary hit
  clade <- vapply(strsplit(ifelse(is.na(best$lineage[i]), "", best$lineage[i]), ";", fixed = TRUE),
    function(p) utils::tail(c(NA_character_, p[nzchar(p)]), 1), character(1))
  sp <- ifelse(!is.na(tied[i]) & tied[i] > 0L & !is.na(clade),
    sprintf("%s (%d tied taxa)", clade, tied[i] + 1L), sp)
  # a full hit list of one taxon says nothing about who else matches
  sat <- !is.na(i) & as.integer(n_raw[best$qseqid[i]]) >= cap & as.integer(n_taxa[best$qseqid[i]]) == 1L
  if (any(sat)) {
    sp[sat] <- sprintf("%s fills the hit list", extract_species(best$title[i[sat]]))
    cut1 <- function(x) vapply(x, chimera_cut_lineage, character(1), rank = RRNA_LINEAGE_RANKS[1], USE.NAMES = FALSE)
    best$lineage[i[sat]] <- cut1(best$lineage[i[sat]]); lt[sat] <- cut1(lt[sat])
  }
  data.frame(start = tiles$start, stop = tiles$stop,
    accession = best$accession[i], species = sp,
    staxid = best$staxid[i], lineage = best$lineage[i],
    lineage_taxids = lt, rank_distance = rd,
    hit_class = classify_hit_title(best$title[i]),
    pident = best$pident[i], bitscore = best$bitscore[i], db_source = db_source,
    stringsAsFactors = FALSE)
}

#' Name at a rank from a ';'-joined 7-field lineage
#' @noRd
chimera_rank_name <- function(lineage, rank = "class") {
  k <- match(rank, RRNA_LINEAGE_RANKS)
  vapply(strsplit(lineage %||% NA_character_, ";", fixed = TRUE), function(p) {
    if (length(p) >= k && nzchar(p[k])) p[k] else NA_character_
  }, character(1))
}

#' Collapse windows into runs sharing a name at the deciding rank
#' @noRd
chimera_segments <- function(windows, deciding = "class") {
  if (!nrow(windows)) return(data.frame())
  key <- chimera_rank_name(windows$lineage, deciding)
  key[is.na(key) & !is.na(windows$lineage)] <- "unresolved"
  key[is.na(key)] <- ifelse(is.na(windows$pident[is.na(key)]), "nohit",
    ifelse(windows$pident[is.na(key)] >= CHIMERA_PID_DROP, "high_id", "low_id"))
  run <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
  do.call(rbind, lapply(split(seq_along(run), run), function(ix) {
    w <- windows[ix, , drop = FALSE]
    top <- ix[which.max(w$bitscore %||% rep(0, nrow(w)))]
    data.frame(start = w$start[1], stop = w$stop[nrow(w)],
      taxon = key[ix[1]], rank = deciding,
      species = windows$species[top] %||% NA_character_,
      pident = if (all(is.na(w$pident))) NA_real_ else mean(w$pident, na.rm = TRUE), n_windows = nrow(w),
      db_source = w$db_source[1], lineage_taxids = windows$lineage_taxids[top] %||% NA_character_,
      stringsAsFactors = FALSE)
  }))
}

#' Channel 2: adjacent segments of comparable length in different classes
#' @noRd
chimera_channel_windows <- function(segments, identity_only = FALSE) {
  if (!nrow(segments)) return(chimera_unavailable())
  s <- segments[segments$n_windows >= 2 & !segments$taxon %in% c("nohit", "unresolved"), , drop = FALSE]
  if (nrow(s) < 2) return(list(score = 0, junctions = integer(0),
    reason = sprintf("one lineage along the contig: %s", s$taxon[1] %|NA|% segments$taxon[1]),
    availability = "available"))
  best <- NULL
  for (i in seq_len(nrow(s) - 1)) {
    a <- s[i, ]; b <- s[i + 1, ]
    if (identical(a$taxon, b$taxon)) next
    la <- a$stop - a$start + 1; lb <- b$stop - b$start + 1
    if (min(la, lb) < 0.25 * max(la, lb)) next
    k <- if (identity_only) NA_integer_ else chimera_split_rank(a$lineage_taxids, b$lineage_taxids)
    if (!identity_only && (is.na(k) || k > 3L)) next
    if (!identity_only && k == 3L && min(a$n_windows, b$n_windows) < CHIMERA_CLASS_MIN_WINDOWS) next
    sv <- !identity_only && isTRUE(k <= 2L && min(a$pident, b$pident) >= CHIMERA_PID_DROP)
    if (is.null(best) || sv > best$sv || (sv == best$sv && isTRUE(k < best$k)))
      best <- list(a = a, b = b, k = k, sv = sv, junction = (a$stop + b$start) %/% 2L)
  }
  if (is.null(best)) return(list(score = 0, junctions = integer(0),
    reason = "no class-level switch between comparable segments", availability = "available"))
  desc <- function(x) sprintf("%d-%d %s (%s, %.1f%%)", x$start, x$stop, x$taxon, x$species %||% "", x$pident)
  list(score = if (identity_only) 0 else 0.5, junctions = best$junction,
    reason = paste(c(desc(best$a), desc(best$b),
      if (!identity_only) sprintf("different %s", RRNA_LINEAGE_RANKS[best$k])), collapse = "; "),
    availability = if (identity_only) "identity_only" else "available",
    severe = best$sv)
}

CHIMERA_RIBOTYPER_SCORES <- c(InconsistentHits = 0.8, DuplicateRegion = 0.8, BothStrands = 0.8,
  MultipleFamilies = 0.6, LowScoreDifference = 0.2)

#' @noRd
chimera_parse_ribotyper <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  ln <- readLines(file, warn = FALSE)
  ln <- ln[!grepl("^#", ln) & nzchar(trimws(ln))]
  if (!length(ln)) return(data.frame(target = character(0), pass = logical(0), codes = I(list())))
  f <- strsplit(trimws(ln), "\\s+")
  target <- vapply(f, `[`, character(1), 2)
  pass <- vapply(f, `[`, character(1), 5) == "PASS"
  feat <- vapply(f, function(x) if (length(x) >= 6) paste(x[6:length(x)], collapse = "") else "-", character(1))
  codes <- lapply(feat, function(x) {
    if (x == "-") return(character(0))
    p <- strsplit(x, ";", fixed = TRUE)[[1]]
    p <- p[grepl("^\\*", p)]
    sub("^\\*([A-Za-z]+).*$", "\\1", p)
  })
  data.frame(target = target, pass = pass, codes = I(codes), stringsAsFactors = FALSE)
}

#' Channel 3: ribotyper structural failures, ITSx chimeric as a weak add-on
#' spans is accepted for a future model-coordinate junction mapping; unused in v1.
#' @noRd
chimera_channel_ribotyper <- function(rt, spans, itsx_chimeric = FALSE) {
  score <- if (isTRUE(itsx_chimeric)) 0.2 else 0
  parts <- if (isTRUE(itsx_chimeric)) "ITSx chimeric" else character(0)
  if (is.null(rt) || !nrow(rt)) {
    if (score == 0) return(chimera_unavailable())
    return(list(score = score, junctions = integer(0), reason = parts, availability = "available"))
  }
  for (i in seq_len(nrow(rt))) {
    cd <- rt$codes[[i]]; cd <- cd[cd %in% names(CHIMERA_RIBOTYPER_SCORES)]
    if (!length(cd)) next
    score <- max(score, CHIMERA_RIBOTYPER_SCORES[cd])
    parts <- c(parts, sprintf("%s %s", sub("^u", "", rt$target[i]), paste(cd, collapse = ";")))
  }
  list(score = unname(score), junctions = integer(0),
    reason = if (length(parts)) paste(parts, collapse = "; ") else "ribotyper passed",
    availability = "available")
}

#' Positions inside a sequence repeated elsewhere on the scaffold (exact k-mers)
#' @noRd
chimera_repeat_mask <- function(seq, k = CHIMERA_REPEAT_K) {
  n <- nchar(seq); m <- logical(n)
  if (n < k) return(m)
  km <- substring(seq, seq_len(n - k + 1L), k:n)
  for (i in which(duplicated(km) | duplicated(km, fromLast = TRUE))) m[i:(i + k - 1L)] <- TRUE
  m
}

#' @noRd
chimera_depth_steps <- function(depths, feats, seq = NULL, step = CHIMERA_STEP, ratio = CHIMERA_DEPTH_STEP,
                                edge = 100L, min_depth = CHIMERA_MIN_DEPTH) {
  n <- length(depths)
  if (n < 4L * step) return(integer(0))
  if (!is.null(seq)) depths[which(chimera_repeat_mask(seq)[seq_len(n)])] <- NA
  med <- function(a, z) stats::median(depths[a:z], na.rm = TRUE)
  edges <- c(feats$start, feats$stop)
  out <- integer(0)
  for (p in seq.int(2L * step, n - 2L * step, by = step)) {
    l <- med(p - step + 1L, p); r <- med(p + 1L, p + step)
    if (is.na(l) || is.na(r) || min(l, r) <= 0 || max(l, r) < min_depth) next
    if (max(l, r) / min(l, r) < ratio) next
    if (length(edges) && any(abs(edges - p) <= edge)) next
    out <- c(out, as.integer(p))
  }
  out
}

#' @noRd
chimera_clip_pileups <- function(depths, clips, min_depth = 10, frac = CHIMERA_CLIP_FRAC, end = CHIMERA_STEP) {
  n <- min(length(depths), length(clips))
  if (n <= 2 * end) return(integer(0))
  i <- (end + 1):(n - end)
  ok <- depths[i] >= min_depth & clips[i] / pmax(depths[i], 1) >= frac
  as.integer(i[which(ok)])
}

#' Channel 4: depth steps and clip pileups
#' @noRd
chimera_channel_reads <- function(depths, clips, feats, seq = NULL) {
  if (is.null(depths) || !length(depths)) return(chimera_unavailable())
  steps <- chimera_depth_steps(depths, feats, seq)
  pile <- if (!is.null(clips) && length(clips)) chimera_clip_pileups(depths, clips) else integer(0)
  pile <- chimera_merge_positions(pile, 50L)
  score <- max(if (length(steps)) 0.4 else 0, if (length(pile)) 0.7 else 0)
  parts <- c(if (length(steps)) sprintf("depth step at %s", paste(steps, collapse = ",")),
             if (length(pile)) sprintf("clip pileup at %s", paste(pile, collapse = ",")))
  list(score = score, junctions = sort(unique(c(steps, pile))),
    reason = if (length(parts)) paste(parts, collapse = "; ") else "no read-evidence breakpoints",
    availability = if (is.null(clips) || !length(clips)) "depth_only" else "available")
}

#' Collapse positions closer than `gap` into their first member
#' @noRd
chimera_merge_positions <- function(x, gap) {
  if (!length(x)) return(integer(0))
  x <- sort(unique(x)); keep <- c(TRUE, diff(x) > gap); x[keep]
}

#' Chimera label: "severe" when a taxonomy channel saw a phylum- or
#' kingdom-level switch, "minor" for any other evidence, else "clean"
#' @noRd
chimera_label <- function(score, severe = FALSE) {
  if (is.na(score)) return(NA_character_)
  if (isTRUE(severe)) "severe" else if (score < 0.3) "clean" else "minor"
}

#' Combine channels: strongest score, severe flag, corroborated junctions
#' @noRd
chimera_score <- function(channels) {
  sc <- vapply(channels, function(x) x$score %||% 0, numeric(1))
  js <- lapply(channels, function(x) as.integer(x$junctions %||% integer(0)))
  corroborated <- integer(0)
  for (i in seq_along(js)) for (j in seq_along(js)) if (i < j) for (a in js[[i]]) {
    near <- js[[j]][abs(js[[j]] - a) <= CHIMERA_JUNCTION_TOL]
    if (length(near)) corroborated <- c(corroborated, min(a, near))
  }
  score <- max(sc, 0)
  severe <- any(vapply(channels, function(x) isTRUE(x$severe), logical(1)))
  list(score = score, label = chimera_label(score, severe),
    junctions = if (length(corroborated)) chimera_merge_positions(corroborated, CHIMERA_JUNCTION_TOL)
      else sort(unique(unlist(js))))
}

#' One row per possible junction: the checks that point to it and why, in words
#' @param junctions Comma-separated positions. @param channels Channels JSON.
#' @noRd
chimera_junction_table <- function(junctions, channels) {
  j <- suppressWarnings(as.integer(strsplit(junctions %|NA|% "", ",", fixed = TRUE)[[1]]))
  j <- j[!is.na(j)]
  out <- data.frame(position = j, checks = character(length(j)), why = character(length(j)))
  if (!length(j) || is.na(channels)) return(out)
  ch <- jsonlite::fromJSON(channels, simplifyVector = FALSE)
  near <- function(x, p) any(abs(as.integer(unlist(x)) - p) <= CHIMERA_JUNCTION_TOL)
  at <- function(reason, what) {
    m <- regmatches(reason, regexpr(paste0(what, " at [0-9,]+"), reason))
    if (length(m)) as.integer(strsplit(sub(".* at ", "", m), ",", fixed = TRUE)[[1]]) else integer(0)
  }
  for (i in seq_along(j)) {
    p <- j[i]
    who <- names(ch)[vapply(ch, function(x) near(x$junctions, p), logical(1))]
    why <- character(0)
    if ("genes" %in% who) why <- c(why, paste0("The genes' own BLAST hits disagree: ", ch$genes$reason, "."))
    if ("windows" %in% who) why <- c(why, paste0("BLAST windows switch lineage here: ", ch$windows$reason, "."))
    if ("reads" %in% who) {
      r <- ch$reads$reason %||% ""
      for (d in at(r, "depth step")[abs(at(r, "depth step") - p) <= CHIMERA_JUNCTION_TOL])
        why <- c(why, sprintf("Read depth changes at least %gx at %d bp (median of %d bp on each side).",
          CHIMERA_DEPTH_STEP, d, CHIMERA_STEP))
      for (d in at(r, "clip pileup")[abs(at(r, "clip pileup") - p) <= CHIMERA_JUNCTION_TOL])
        why <- c(why, sprintf("At least %d%% of reads are clipped at %d bp, so reads stop matching there.",
          round(100 * CHIMERA_CLIP_FRAC), d))
    }
    if (identical(who, "reads")) why <- c(why,
      "Only read evidence points here; read evidence alone gives a minor call at most.")
    out$checks[i] <- paste(who, collapse = ", ")
    out$why[i] <- paste(why, collapse = " ")
  }
  out
}

#' The call in words: one "check: verdict" clause per check, the check that
#' raised the call first, the junction last. The segments table carries the
#' per-segment detail, so it is not repeated here.
#' @noRd
chimera_reason <- function(segments, channels, junctions) {
  checks <- c("genes", "windows", "ribotyper", "reads")
  txt <- vapply(checks, function(k) {
    r <- channels[[k]]$reason
    if (is.null(r) || is.na(r)) {
      if (k == "windows" && !is.null(segments) && nrow(segments)) r <- paste(unique(segments$taxon), collapse = ", ")
      else return(NA_character_)
    }
    if (startsWith(r, k)) r else paste0(k, ": ", r)
  }, character(1))
  score <- vapply(checks, function(k) channels[[k]]$score %||% 0, numeric(1))
  txt <- txt[order(-score)]
  txt <- txt[!is.na(txt)]
  jtxt <- if (length(junctions)) {
    who <- vapply(junctions, function(j) paste(names(channels)[vapply(channels, function(x)
      any(abs((x$junctions %||% integer(0)) - j) <= CHIMERA_JUNCTION_TOL), logical(1))], collapse = ","), character(1))
    paste(sprintf("junction near %d (%s)", junctions, who), collapse = "; ")
  } else "no junction"
  paste(c(txt, jtxt), collapse = "; ")
}

#' Split a stored reason back into its clauses
#' @noRd
chimera_reason_parts <- function(reason) {
  trimws(strsplit(reason, "; (?=(genes|windows|ribotyper|reads|junction near|junctions trimmed|no junction))", perl = TRUE)[[1]])
}

#' @noRd
chimera_run_blastn <- function(query, dbarg, task, out, threads = 1L, blastn = "blastn", max_targets = 50L) {
  system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query), "-task", task,
    "-outfmt", shQuote("6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"),
    "-max_target_seqs", as.character(max_targets), "-max_hsps", "1", "-num_threads", threads), stdout = out, stderr = FALSE)
}

#' @noRd
chimera_run_ribotyper <- function(fasta, outdir, ribotyper = "ribotyper") {
  system2(ribotyper, c("-f", shQuote(fasta), shQuote(outdir)), stdout = FALSE, stderr = FALSE)
}

#' Turn a labels FASTA into a usable blastn -db argument, building the BLAST
#' db when one is not already sitting beside the FASTA.
#' @noRd
chimera_label_prefix <- function(label_db, workdir) {
  if (is.null(label_db) || !nzchar(label_db)) return(label_db)
  if (!grepl("\\.(fasta|fa)$", label_db, ignore.case = TRUE)) return(label_db)
  if (file.exists(paste0(label_db, ".nin"))) return(label_db)
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  prefix <- file.path(workdir, "labels")
  system2("makeblastdb", c("-in", shQuote(label_db), "-dbtype", "nucl", "-out", shQuote(prefix)),
    stdout = FALSE, stderr = FALSE)
  prefix
}

#' Run every channel for one unit and write chimera + chimera_windows
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param dbarg blastn -db argument for the windowed local BLAST channel, or
#'   NULL/"" to fall back to `label_db` in identity-only mode.
#' @param taxdump_dir Local NCBI taxdump directory (or NULL for the bundled one).
#' @param label_db Path to a labels FASTA/BLAST db used when `dbarg` is absent.
#' @param threads blastn thread count.
#' @param workdir Scratch directory for intermediate files.
#' @param review TRUE checks the contig's review copy: its sequence, features,
#'   hits, and ITSx chimera flag come from the review tables, coverage is cut
#'   to the copy's window, and the result goes to `review_chimera` and
#'   `review_chimera_windows`.
#' @return (invisibly) the written `chimera` row, as a list.
#' @export
chimera_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                         label_db = NULL, threads = 1L, workdir = tempfile("chimera_"), review = FALSE) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  label_db <- chimera_label_prefix(label_db, workdir)
  seq <- DBI::dbGetQuery(con, sprintf("SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("assemblies", review)), params = list(id, path, scaffold))$sequence[1]
  if (is.na(seq %||% NA) || !nzchar(seq)) return(invisible(NULL))
  feats <- DBI::dbGetQuery(con, sprintf("SELECT element, start, stop FROM %s WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    review_tbl("annotations", review)), params = list(id, path, scaffold))
  hits <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT accession, bitscore, qstart, qend, lineage_taxids, lineage, species, pident FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("blast_hits", review)), params = list(id, path, scaffold)), error = function(e) data.frame())
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  itsx <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT chimeric FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    review_tbl("annotate", review)), params = list(id, path, scaffold))$chimeric[1], error = function(e) NA)
  track <- function(x) if (length(x) && !is.na(x) && nzchar(x)) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]) else NULL
  cov <- if (review) review_coverage(con, id, path, scaffold) else {
    cv <- tryCatch(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold)), error = function(e) data.frame())
    list(depths = if (nrow(cv)) track(cv$depths[1]), clips = if (nrow(cv) && "clips" %in% names(cv)) track(cv$clips[1]))
  }

  # channel 1
  gene_taxa <- if (nrow(hits) && nrow(feats)) chimera_gene_lca(hits, feats) else NULL
  ch_genes <- if (is.null(gene_taxa)) chimera_unavailable() else chimera_channel_genes(gene_taxa, feats)

  # channel 2
  tiles <- chimera_tile(nchar(seq))
  qf <- file.path(workdir, "windows.fa"); chimera_windows_fasta(seq, tiles, qf)
  identity_only <- is.null(dbarg) || !nzchar(dbarg)
  db_use <- if (identity_only) label_db else dbarg
  windows <- NULL; ch_windows <- chimera_unavailable(); segments <- NULL
  if (!is.null(db_use) && nzchar(db_use)) {
    out1 <- file.path(workdir, "mega.txt")
    chimera_run_blastn(qf, db_use, "megablast", out1, threads, max_targets = CHIMERA_WIN_MAX_TARGETS)
    h <- parse_blast_hits_file(out1)
    missing <- setdiff(sprintf("w%d_%d", tiles$start, tiles$stop), unique(h$qseqid))
    if (length(missing)) {
      keep <- tiles[sprintf("w%d_%d", tiles$start, tiles$stop) %in% missing, , drop = FALSE]
      qf2 <- file.path(workdir, "windows2.fa"); chimera_windows_fasta(seq, keep, qf2)
      out2 <- file.path(workdir, "dc.txt")
      chimera_run_blastn(qf2, db_use, "dc-megablast", out2, threads, max_targets = CHIMERA_WIN_MAX_TARGETS)
      h <- rbind(h, parse_blast_hits_file(out2))
    }
    windows <- chimera_window_hits(h, tiles, samp, taxdump_dir, db_source = db_use)
    no_taxonomy <- !identity_only && all(is.na(windows$lineage))
    if (no_taxonomy) identity_only <- TRUE
    segments <- chimera_segments(windows)
    ch_windows <- chimera_channel_windows(segments, identity_only = identity_only)
    if (no_taxonomy) ch_windows$availability <- "no_taxonomy"
  }

  # channel 3
  ch_rt <- chimera_unavailable(); codes <- NA_character_
  rt_feats <- feats[feats$element %in% c("18S", "28S"), , drop = FALSE]
  if (nrow(rt_feats) && nzchar(Sys.which("ribotyper"))) {
    rf <- file.path(workdir, "rt.fa")
    st <- pmax(1L, rt_feats$start - 50L); en <- pmin(nchar(seq), rt_feats$stop + 50L)
    writeLines(as.vector(rbind(paste0(">u", rt_feats$element), substring(seq, st, en))), rf)
    rdir <- file.path(workdir, "rt"); dir.create(rdir, showWarnings = FALSE)
    chimera_run_ribotyper(rf, rdir)
    rt <- chimera_parse_ribotyper(list.files(rdir, "ribotyper\\.short\\.out$", full.names = TRUE)[1])
    ch_rt <- chimera_channel_ribotyper(rt, data.frame(target = paste0("u", rt_feats$element), start = st),
      itsx_chimeric = isTRUE(itsx == 1L))
    if (!is.null(rt)) codes <- paste(unlist(rt$codes), collapse = ";")
  } else if (isTRUE(itsx == 1L)) {
    ch_rt <- chimera_channel_ribotyper(NULL, NULL, itsx_chimeric = TRUE)
  }

  # channel 4
  depths <- cov$depths
  clips <- cov$clips
  ch_reads <- chimera_channel_reads(depths, clips, feats, seq)

  channels <- list(genes = ch_genes, windows = ch_windows, ribotyper = ch_rt, reads = ch_reads)
  sc <- chimera_score(channels)
  reason <- chimera_reason(segments, channels, sc$junctions)

  if (!is.null(windows) && nrow(windows) && !is.null(depths)) {
    windows$depth_median <- vapply(seq_len(nrow(windows)), function(i)
      stats::median(depths[windows$start[i]:min(windows$stop[i], length(depths))], na.rm = TRUE), numeric(1))
    windows$clip_frac <- if (!is.null(clips)) vapply(seq_len(nrow(windows)), function(i) {
      ix <- windows$start[i]:min(windows$stop[i], length(clips))
      max(clips[ix] / pmax(depths[ix], 1), na.rm = TRUE) }, numeric(1)) else NA_real_
  } else if (!is.null(windows)) { windows$depth_median <- NA_real_; windows$clip_frac <- NA_real_ }

  ts <- as.integer(Sys.time())
  row <- data.frame(ID = id, path = path, scaffold = scaffold, score = sc$score, label = sc$label,
    junctions = paste(sc$junctions, collapse = ","), reason = reason,
    segments = if (!is.null(segments) && nrow(segments)) jsonlite::toJSON(segments, dataframe = "rows") else NA_character_,
    gene_taxa = if (!is.null(gene_taxa)) jsonlite::toJSON(gene_taxa, dataframe = "rows") else NA_character_,
    channels = jsonlite::toJSON(channels, auto_unbox = TRUE),
    ribotyper_codes = codes, db_source = db_use %||% NA_character_, time_stamp = ts,
    stringsAsFactors = FALSE)
  ch_tbl <- review_tbl("chimera", review)
  win_tbl <- review_tbl("chimera_windows", review)
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", ch_tbl),
      params = list(id, path, scaffold))
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", win_tbl),
      params = list(id, path, scaffold))
    DBI::dbAppendTable(con, ch_tbl, row)
    if (!is.null(windows) && nrow(windows)) {
      windows$ID <- id; windows$path <- path; windows$scaffold <- scaffold
      DBI::dbAppendTable(con, win_tbl, windows[, c("ID", "path", "scaffold", "start", "stop", "accession",
        "species", "staxid", "lineage", "lineage_taxids", "rank_distance", "hit_class", "pident", "bitscore",
        "depth_median", "clip_frac", "db_source")])
    }
  })
  as.list(row)
}

#' Move a unit's stored chimera coordinates after its sequence was trimmed or restored
#'
#' New position = old position - `offset`. Windows and junctions that fall
#' outside 1..`len` are dropped, and the segments are rebuilt from the windows kept.
#' @noRd
chimera_shift <- function(con, id, path, scaffold, offset, len) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  ch <- DBI::dbGetQuery(con, "SELECT junctions, label, score FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?",
    params = p)
  if (!nrow(ch)) return(invisible())
  w <- DBI::dbGetQuery(con, "SELECT * FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = p)
  w$start <- w$start - offset; w$stop <- w$stop - offset
  w <- w[w$start >= 1 & w$stop <= len, , drop = FALSE]
  j0 <- suppressWarnings(as.integer(strsplit(ch$junctions[1] %|NA|% "", ",", fixed = TRUE)[[1]]))
  j <- j0 - offset
  j <- j[!is.na(j) & j >= 1 & j <= len]
  seg <- chimera_segments(w)
  gone <- any(!is.na(j0)) && !length(j)
  reason <- gsub(" ()", "", chimera_reason(seg, list(), j), fixed = TRUE)
  if (gone) reason <- sub("no junction", "junctions trimmed away", reason, fixed = TRUE)
  with_savepoint(con, {
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    if (nrow(w)) DBI::dbAppendTable(con, "chimera_windows", w)
    DBI::dbExecute(con,
      "UPDATE chimera SET junctions = ?, segments = ?, reason = ?, label = ?, score = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(paste(j, collapse = ","),
        if (nrow(seg)) as.character(jsonlite::toJSON(seg, dataframe = "rows")) else NA_character_,
        reason, if (gone) "clean" else ch$label[1], if (gone) 0 else ch$score[1]), p))
  })
  invisible()
}

#' Put back a chimera result and windows saved as JSON by an older trim
#' @noRd
chimera_restore <- function(con, id, path, scaffold, snap) {
  s <- jsonlite::fromJSON(snap)
  p <- list(id, as.integer(path), as.integer(scaffold))
  with_savepoint(con, {
    DBI::dbExecute(con, "DELETE FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    DBI::dbAppendTable(con, "chimera", as.data.frame(s$chimera))
    if (is.data.frame(s$windows) && nrow(s$windows)) DBI::dbAppendTable(con, "chimera_windows", s$windows)
  })
  invisible()
}
