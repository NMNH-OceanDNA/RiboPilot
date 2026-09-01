# Chimera detection: per-gene taxonomy, windowed local BLAST, ribotyper,
# read evidence. Design: dev/design/21_CHIMERA_DETECTION.md.

CHIMERA_WIN <- 500L
CHIMERA_STEP <- 250L
CHIMERA_PID_DROP <- 90
CHIMERA_DEPTH_STEP <- 3
CHIMERA_CLIP_FRAC <- 0.3
CHIMERA_JUNCTION_TOL <- 300L
CHIMERA_TAXONOMY_CAP <- 0.7
CHIMERA_LABELS <- c("clean", "suspect", "likely")

#' @noRd
chimera_tile <- function(seq_len, win = CHIMERA_WIN, step = CHIMERA_STEP) {
  seq_len <- as.integer(seq_len)
  if (seq_len <= win) return(data.frame(start = 1L, stop = seq_len))
  starts <- seq.int(1L, seq_len - win + 1L, by = step)
  if (utils::tail(starts, 1) + win - 1L < seq_len) starts <- c(starts, seq_len - win + 1L)
  data.frame(start = as.integer(starts), stop = as.integer(starts + win - 1L))
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
    parts <- strsplit(h$lineage_taxids[1], ";", fixed = TRUE)[[1]]
    k <- match(lca$rank, RRNA_LINEAGE_RANKS)
    if (!is.na(k) && k < length(parts)) parts[(k + 1L):length(parts)] <- ""
    base$lineage_taxids <- paste(parts, collapse = ";")
    base
  })
  do.call(rbind, out)
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
  for (i in seq_len(nrow(g) - 1)) for (j in (i + 1):nrow(g)) {
    d <- rank_distance(g$lineage_taxids[i], g$lineage_taxids[j])
    if (!is.na(d) && d >= 5 && (is.null(worst) || d > worst$d)) worst <- list(i = i, j = j, d = d)
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
    reason = sprintf("%s resolves to %s (%s, %.1f%%) but %s to %s (%s, %.1f%%)",
      a$element, a$taxon, a$species, a$pident, b$element, b$taxon, b$species, b$pident),
    availability = "available")
}

#' Best hit per tile from a parse_blast_hits_file() frame
#' @noRd
chimera_window_hits <- function(hits, tiles, sample_taxids = NA_character_,
                                taxdump_dir = NULL, db_source = NA_character_) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  if (!nrow(hits)) {
    return(data.frame(start = tiles$start, stop = tiles$stop,
      accession = NA_character_, species = NA_character_,
      staxid = NA_integer_, lineage = NA_character_,
      lineage_taxids = NA_character_, rank_distance = NA_integer_,
      hit_class = NA_character_, pident = NA_real_, bitscore = NA_real_,
      db_source = db_source, stringsAsFactors = FALSE))
  }
  best <- hits[order(hits$qseqid, -hits$bitscore), , drop = FALSE]
  best <- best[!duplicated(best$qseqid), , drop = FALSE]
  i <- match(ids, best$qseqid)
  lin <- tryCatch(tk_ranked_lineage(best$staxid, taxdump_dir), error = function(e) data.frame())
  li <- if (nrow(lin)) match(best$staxid[i], lin$taxid) else rep(NA_integer_, length(i))
  lt <- if (nrow(lin)) lin$lineage_taxids[li] else rep(NA_character_, length(i))
  rd <- if (!is.na(sample_taxids) && nzchar(sample_taxids)) {
    vapply(lt, function(x) if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, sample_taxids)),
      integer(1), USE.NAMES = FALSE)
  } else rep(NA_integer_, length(i))
  data.frame(start = tiles$start, stop = tiles$stop,
    accession = best$accession[i], species = extract_species(best$title[i]),
    staxid = best$staxid[i], lineage = if (nrow(lin)) lin$lineage[li] else NA_character_,
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
  key[is.na(key)] <- ifelse(is.na(windows$pident[is.na(key)]), "nohit",
    ifelse(windows$pident[is.na(key)] >= CHIMERA_PID_DROP, "high_id", "low_id"))
  run <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
  do.call(rbind, lapply(split(seq_along(run), run), function(ix) {
    w <- windows[ix, , drop = FALSE]
    top <- ix[which.max(w$bitscore %||% rep(0, nrow(w)))]
    data.frame(start = w$start[1], stop = w$stop[nrow(w)],
      taxon = key[ix[1]], rank = deciding,
      species = windows$species[top] %||% NA_character_,
      pident = mean(w$pident, na.rm = TRUE), n_windows = nrow(w),
      db_source = w$db_source[1], lineage_taxids = windows$lineage_taxids[top] %||% NA_character_,
      stringsAsFactors = FALSE)
  }))
}

#' Channel 2: adjacent segments of comparable length in different classes
#' @noRd
chimera_channel_windows <- function(segments, identity_only = FALSE) {
  if (!nrow(segments)) return(chimera_unavailable())
  s <- segments[segments$n_windows >= 2, , drop = FALSE]
  if (nrow(s) < 2) return(list(score = 0, junctions = integer(0),
    reason = sprintf("one lineage along the scaffold: %s", segments$taxon[1]), availability = "available"))
  best <- NULL
  for (i in seq_len(nrow(s) - 1)) {
    a <- s[i, ]; b <- s[i + 1, ]
    if (identical(a$taxon, b$taxon)) next
    if (a$taxon %in% c("nohit") || b$taxon %in% c("nohit")) next
    la <- a$stop - a$start + 1; lb <- b$stop - b$start + 1
    if (min(la, lb) < 0.25 * max(la, lb)) next
    if (!identity_only) {
      if (is.na(a$lineage_taxids) || is.na(b$lineage_taxids)) next
      if (rank_distance(a$lineage_taxids, b$lineage_taxids) < 5) next
    }
    best <- list(a = a, b = b, junction = (a$stop + b$start) %/% 2L)
  }
  if (is.null(best)) return(list(score = 0, junctions = integer(0),
    reason = "no class-level switch between comparable segments", availability = "available"))
  desc <- function(x) sprintf("%d-%d %s (%s, %.1f%%)", x$start, x$stop, x$taxon, x$species %||% "", x$pident)
  list(score = if (identity_only) 0 else 0.5, junctions = best$junction,
    reason = paste(desc(best$a), desc(best$b), sep = "; "),
    availability = if (identity_only) "identity_only" else "available")
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

#' @noRd
chimera_depth_steps <- function(depths, feats, step = CHIMERA_STEP, ratio = CHIMERA_DEPTH_STEP, edge = 100L) {
  n <- length(depths)
  if (n < 2 * step) return(integer(0))
  b <- seq.int(step, n - step, by = step)
  med <- function(a, z) stats::median(depths[a:z], na.rm = TRUE)
  out <- integer(0)
  for (p in b) {
    l <- med(max(1, p - step + 1), p); r <- med(p + 1, min(n, p + step))
    if (is.na(l) || is.na(r) || min(l, r) <= 0) next
    if (max(l, r) / min(l, r) < ratio) next
    edges <- c(feats$start, feats$stop)
    if (length(edges) && any(abs(edges - p) <= edge)) next
    out <- c(out, as.integer(p))
  }
  out
}

#' @noRd
chimera_clip_pileups <- function(depths, clips, min_depth = 10, frac = CHIMERA_CLIP_FRAC, end = 50L) {
  n <- min(length(depths), length(clips))
  if (n <= 2 * end) return(integer(0))
  i <- (end + 1):(n - end)
  ok <- depths[i] >= min_depth & clips[i] / pmax(depths[i], 1) >= frac
  as.integer(i[which(ok)])
}

#' Channel 4: depth steps and clip pileups
#' @noRd
chimera_channel_reads <- function(depths, clips, feats) {
  if (is.null(depths) || !length(depths)) return(chimera_unavailable())
  steps <- chimera_depth_steps(depths, feats)
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

#' @noRd
chimera_label <- function(score) {
  if (is.na(score)) return(NA_character_)
  if (score < 0.3) "clean" else if (score <= 0.7) "suspect" else "likely"
}

#' Combine channels: max, corroboration bonus, taxonomy-only cap
#' @noRd
chimera_score <- function(channels) {
  sc <- vapply(channels, function(x) x$score %||% 0, numeric(1))
  score <- max(sc, 0)
  js <- lapply(channels, function(x) as.integer(x$junctions %||% integer(0)))
  corroborated <- integer(0)
  nm <- names(js)
  for (i in seq_along(js)) for (j in seq_along(js)) if (i < j) for (a in js[[i]]) {
    near <- js[[j]][abs(js[[j]] - a) <= CHIMERA_JUNCTION_TOL]
    if (length(near)) corroborated <- c(corroborated, min(a, near))
  }
  if (length(corroborated)) score <- min(1, score + 0.25)
  strong <- sc[setdiff(nm, c("genes", "windows"))]
  if (!length(strong) || all(strong == 0)) score <- min(score, CHIMERA_TAXONOMY_CAP)
  all_j <- sort(unique(unlist(js)))
  list(score = score, label = chimera_label(score),
    junctions = if (length(corroborated)) chimera_merge_positions(corroborated, CHIMERA_JUNCTION_TOL) else all_j)
}

#' @noRd
chimera_reason <- function(segments, channels, junctions) {
  segtxt <- if (!is.null(segments) && nrow(segments)) paste(sprintf("%d-%d %s (%s, %.1f%%)",
    segments$start, segments$stop, segments$taxon, segments$species %||% "", segments$pident), collapse = "; ")
    else "no window segments"
  jtxt <- if (length(junctions)) {
    who <- vapply(junctions, function(j) paste(names(channels)[vapply(channels, function(x)
      any(abs((x$junctions %||% integer(0)) - j) <= CHIMERA_JUNCTION_TOL), logical(1))], collapse = ","), character(1))
    paste(sprintf("junction near %d (%s)", junctions, who), collapse = "; ")
  } else "no junction"
  extra <- vapply(channels[c("ribotyper", "reads")], function(x) x$reason %||% NA_character_, character(1))
  extra <- extra[!is.na(extra)]
  paste(c(segtxt, jtxt, extra), collapse = "; ")
}

#' @noRd
chimera_run_blastn <- function(query, dbarg, task, out, threads = 1L, blastn = "blastn") {
  system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query), "-task", task,
    "-outfmt", shQuote("6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"),
    "-max_target_seqs", "5", "-max_hsps", "1", "-num_threads", threads), stdout = out, stderr = FALSE)
}

#' @noRd
chimera_run_ribotyper <- function(fasta, outdir, ribotyper = "ribotyper") {
  system2(ribotyper, c("-f", shQuote(fasta), shQuote(outdir)), stdout = FALSE, stderr = FALSE)
}

#' Run every channel for one unit and write chimera + chimera_windows
#' @noRd
chimera_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                         blastn = "blastn", ribotyper = "ribotyper", workdir = tempfile()) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  seq <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$sequence[1]
  if (is.na(seq %||% NA) || !nzchar(seq)) return(invisible(NULL))
  feats <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  hits <- tryCatch(DBI::dbGetQuery(con, "SELECT accession, bitscore, qstart, qend, lineage_taxids, lineage, species, pident FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  cov <- tryCatch(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  itsx <- tryCatch(DBI::dbGetQuery(con, "SELECT chimeric FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$chimeric[1], error = function(e) NA)
  track <- function(x) if (length(x) && !is.na(x) && nzchar(x)) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]) else NULL

  # channel 1
  gene_taxa <- if (nrow(hits) && nrow(feats)) chimera_gene_lca(hits, feats) else NULL
  ch_genes <- if (is.null(gene_taxa)) chimera_unavailable() else chimera_channel_genes(gene_taxa, feats)

  # channel 2
  tiles <- chimera_tile(nchar(seq))
  qf <- file.path(workdir, "windows.fa"); chimera_windows_fasta(seq, tiles, qf)
  identity_only <- is.null(dbarg) || !nzchar(dbarg)
  db_use <- if (identity_only) NULL else dbarg
  windows <- NULL; ch_windows <- chimera_unavailable(); segments <- NULL
  if (!is.null(db_use) && nzchar(db_use)) {
    out1 <- file.path(workdir, "mega.txt")
    chimera_run_blastn(qf, db_use, "megablast", out1, threads = 1L, blastn = blastn)
    h <- parse_blast_hits_file(out1)
    missing <- setdiff(sprintf("w%d_%d", tiles$start, tiles$stop), unique(h$qseqid))
    if (length(missing)) {
      keep <- tiles[sprintf("w%d_%d", tiles$start, tiles$stop) %in% missing, , drop = FALSE]
      qf2 <- file.path(workdir, "windows2.fa"); chimera_windows_fasta(seq, keep, qf2)
      out2 <- file.path(workdir, "dc.txt")
      chimera_run_blastn(qf2, db_use, "dc-megablast", out2, threads = 1L, blastn = blastn)
      h <- rbind(h, parse_blast_hits_file(out2))
    }
    windows <- chimera_window_hits(h, tiles, samp, taxdump_dir, db_source = db_use)
    segments <- chimera_segments(windows)
    ch_windows <- chimera_channel_windows(segments, identity_only = identity_only)
  }

  # channel 3
  ch_rt <- chimera_unavailable(); codes <- NA_character_
  rt_feats <- feats[feats$element %in% c("18S", "28S"), , drop = FALSE]
  if (nrow(rt_feats) && nzchar(Sys.which(ribotyper))) {
    rf <- file.path(workdir, "rt.fa")
    st <- pmax(1L, rt_feats$start - 50L); en <- pmin(nchar(seq), rt_feats$stop + 50L)
    writeLines(as.vector(rbind(paste0(">u", rt_feats$element), substring(seq, st, en))), rf)
    rdir <- file.path(workdir, "rt"); dir.create(rdir, showWarnings = FALSE)
    chimera_run_ribotyper(rf, rdir, ribotyper = ribotyper)
    rt <- chimera_parse_ribotyper(list.files(rdir, "ribotyper\\.short\\.out$", full.names = TRUE)[1])
    ch_rt <- chimera_channel_ribotyper(rt, data.frame(target = paste0("u", rt_feats$element), start = st),
      itsx_chimeric = isTRUE(itsx == 1L))
    if (!is.null(rt)) codes <- paste(unlist(rt$codes), collapse = ";")
  } else if (isTRUE(itsx == 1L)) {
    ch_rt <- chimera_channel_ribotyper(NULL, NULL, itsx_chimeric = TRUE)
  }

  # channel 4
  depths <- if (nrow(cov)) track(cov$depths[1]) else NULL
  clips <- if (nrow(cov) && "clips" %in% names(cov)) track(cov$clips[1]) else NULL
  ch_reads <- chimera_channel_reads(depths, clips, feats)

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
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbAppendTable(con, "chimera", row)
    if (!is.null(windows) && nrow(windows)) {
      windows$ID <- id; windows$path <- path; windows$scaffold <- scaffold
      DBI::dbAppendTable(con, "chimera_windows", windows[, c("ID", "path", "scaffold", "start", "stop", "accession",
        "species", "staxid", "lineage", "lineage_taxids", "rank_distance", "hit_class", "pident", "bitscore",
        "depth_median", "clip_frac", "db_source")])
    }
  })
  as.list(row)
}
