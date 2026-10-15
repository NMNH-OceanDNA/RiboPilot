# Paths alignment panel: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** In the Assemble details window, align the homologous contigs of a multi-path sample and show the alignment in the sequence viewer with features per row and differing columns highlighted.

**Architecture:** R builds the homologous sets from annotations, aligns each set with `DECIPHER::AlignSeqs` on demand, and sends a `rows` payload to the existing canvas viewer. `seqviewer.js` gains an alignment mode (one row per contig, a Differences track, prev/next difference) behind `if (this.rows)`. Nothing is stored in the database.

**Tech Stack:** R, Shiny, DECIPHER, Biostrings, testthat; hand-rolled canvas JS in `inst/app/www/seqviewer.js` (no JS test harness; the app is verified by driving it with chromote).

**Spec:** `dev/specs/2026-10-15-paths-alignment-design.md`

## Global Constraints

- Work on branch `feat/paths-alignment` (already created off `devel`). Never push.
- No new dependencies. No non-ASCII in R or JS code. Minimal comments.
- Alignment call is exactly `DECIPHER::AlignSeqs(dna, verbose = FALSE)` inside `tryCatch`, as in `compute_feature_msa()` (`R/blast_ref_align.R`).
- Anchors are `SCREEN_ANCHOR_ELEMENTS` (`R/screen_variants.R`): `18S`, `5.8S`, `28S`.
- Feature colours come from `RRNA_ELEMENT_COLORS` (`R/app_modals.R:890`); unknown element is `#888888`.
- The existing per-unit viewer mode must not change behaviour. All JS changes are behind `if (this.rows)` or in new functions.
- Run tests with `Rscript -e 'testthat::test_local(filter = "<name>")'` from the repo root (the package loads with `devtools::load_all()`).
- Commit messages: brief, no attribution lines.

## Review Focus

1. Sample with three single-contig paths: one set of three rows, no dropdown. Test in Task 1 and Task 4.
2. A contig whose sequence is NA in `assemblies` (failed path): skipped before aligning; if fewer than two remain the summary says "Could not align these contigs." Test in Task 3.
3. A feature that spans an alignment gap in its own row: `pos2` lands on the column of its last base, not on the gap. Test in Task 2.
4. Alignment where one row is entirely gaps at the ends (fragment vs full length): the Differences track counts the terminal gap columns as differences (spec: gaps count). Test in Task 2 (`paths_diff_count`).
5. Reopening the modal for the same sample: the cached payload is replayed and no realignment happens. The cache is inline in the server, so this is checked by hand in Task 7 (the summary must appear at once on reopen).

---

### Task 1: Homologous sets (`paths_homology_sets`)

**Files:**
- Create: `R/paths_align.R`
- Test: `tests/testthat/test-paths-align.R`

**Interfaces:**
- Produces: `paths_homology_sets(feats)` where `feats` is a data frame with columns `path`, `scaffold`, `element`. Returns `list()` or a list of `list(label = <chr>, units = data.frame(path = <int>, scaffold = <int>))`.

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-paths-align.R
feats_df <- function(...) {
  rows <- list(...)
  data.frame(path = vapply(rows, `[[`, 1, 1), scaffold = vapply(rows, `[[`, 1, 2),
    element = vapply(rows, `[[`, "", 3), stringsAsFactors = FALSE)
}

test_that("contigs sharing an anchor form one set; stray contigs are dropped", {
  f <- feats_df(list(1, 1, "18S"), list(1, 1, "28S"), list(2, 1, "18S"), list(2, 2, "28S"),
    list(3, 1, "ITS1"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(s[[1]]$units, data.frame(path = c(1L, 2L, 2L), scaffold = c(1L, 1L, 2L)), ignore_attr = TRUE)
  expect_equal(s[[1]]$label, "18S-28S (3 contigs)")
})

test_that("three single-contig paths are one set of three", {
  f <- feats_df(list(1, 1, "18S"), list(2, 1, "18S"), list(3, 1, "18S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(nrow(s[[1]]$units), 3L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
})

test_that("disjoint sets come out largest first and singletons vanish", {
  f <- feats_df(list(1, 1, "28S"), list(2, 1, "28S"), list(1, 2, "18S"), list(2, 2, "18S"),
    list(3, 2, "18S"), list(4, 1, "5.8S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 2L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
  expect_equal(s[[2]]$label, "28S (2 contigs)")
  expect_equal(s[[2]]$units$path, c(1L, 2L))
})

test_that("no anchors or empty input gives an empty list", {
  expect_equal(RiboPilot:::paths_homology_sets(feats_df(list(1, 1, "ITS1"), list(2, 1, "ITS2"))), list())
  expect_equal(RiboPilot:::paths_homology_sets(data.frame(path = integer(), scaffold = integer(),
    element = character())), list())
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: FAIL, `paths_homology_sets` not found.

- [ ] **Step 3: Implement**

```r
# R/paths_align.R
# Paths alignment panel: homologous contigs across a sample's assembly paths,
# aligned with DECIPHER and shown in the sequence viewer's alignment mode.

#' Homologous sets of contigs: connected components of "shares an anchor"
#' @param feats data.frame(path, scaffold, element), one row per annotation.
#' @return list of list(label, units = data.frame(path, scaffold)); only sets
#'   of two or more contigs, largest first.
#' @noRd
paths_homology_sets <- function(feats) {
  f <- feats[feats$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (!nrow(f)) return(list())
  f$unit <- paste(f$path, f$scaffold)
  units <- unique(f$unit)
  comp <- stats::setNames(seq_along(units), units)
  for (el in unique(f$element)) {
    members <- unique(f$unit[f$element == el])
    comp[comp %in% comp[members]] <- min(comp[members])
  }
  sets <- split(names(comp), comp)
  sets <- Filter(function(u) length(u) >= 2, sets)
  if (!length(sets)) return(list())
  out <- lapply(sets, function(u) {
    pu <- do.call(rbind, strsplit(u, " ", fixed = TRUE))
    units <- data.frame(path = as.integer(pu[, 1]), scaffold = as.integer(pu[, 2]))
    units <- units[order(units$path, units$scaffold), , drop = FALSE]
    rownames(units) <- NULL
    anchors <- SCREEN_ANCHOR_ELEMENTS[SCREEN_ANCHOR_ELEMENTS %in% f$element[f$unit %in% u]]
    list(label = sprintf("%s (%d contigs)", paste(anchors, collapse = "-"), nrow(units)), units = units)
  })
  n <- vapply(out, function(s) nrow(s$units), 1L)
  first <- vapply(out, function(s) s$units$path[1] * 1e6 + s$units$scaffold[1], 1)
  unname(out[order(-n, first)])
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add R/paths_align.R tests/testthat/test-paths-align.R
git commit -m "feat: homologous contig sets across assembly paths"
```

---

### Task 2: Alignment, payload, and difference count

**Files:**
- Modify: `R/paths_align.R`
- Test: `tests/testthat/test-paths-align.R`

**Interfaces:**
- Produces:
  - `paths_align(seqs)`: named character vector in; `Biostrings::DNAStringSet` (aligned, names kept) or `NULL`.
  - `paths_payload(aln, feats, labels)`: `aln` is the aligned `DNAStringSet`; `feats` is a data frame with `unit`, `element`, `start`, `stop`, `strand`, `partial_start`, `partial_stop`; `labels` is a named character vector keyed by unit. Returns `list(len, rows)`; each row is `list(unit, label, seq, features)`; each feature is `list(row, gene, pos1, pos2, dir, partial5, partial3, color)`.
  - `paths_diff_count(rows)`: integer count of columns where the rows' characters (gap included) are not all equal.

- [ ] **Step 1: Write the failing tests**

```r
# append to tests/testthat/test-paths-align.R
aln_fixture <- function() {
  Biostrings::DNAStringSet(c(`p1.s1` = "ACGTACGTAC", `p2.s1` = "ACG--CGTAC"))
}
no_feats <- function() data.frame(unit = character(), element = character(), start = numeric(),
  stop = numeric(), strand = character(), partial_start = numeric(), partial_stop = numeric())

test_that("the payload maps features from contig to alignment columns across gaps", {
  f <- data.frame(unit = c("p1.s1", "p2.s1", "p2.s1"), element = c("18S", "18S", "ITS1"),
    start = c(2, 2, NA), stop = c(5, 5, 7), strand = "+", partial_start = c(0, 1, 0),
    partial_stop = c(0, 0, 0), stringsAsFactors = FALSE)
  p <- RiboPilot:::paths_payload(aln_fixture(), f, c(p1.s1 = "p1.s1 keep", p2.s1 = "p2.s1 drop"))
  expect_equal(p$len, 10L)
  expect_length(p$rows, 2L)
  expect_equal(p$rows[[2]]$label, "p2.s1 drop")
  expect_equal(p$rows[[2]]$seq, "ACG--CGTAC")
  r1 <- p$rows[[1]]$features[[1]]; r2 <- p$rows[[2]]$features[[1]]
  expect_equal(c(r1$pos1, r1$pos2), c(2L, 5L))
  expect_equal(c(r2$pos1, r2$pos2), c(2L, 7L))   # bases 2..5 of p2.s1 sit at columns 2,3,6,7
  expect_true(r2$partial5); expect_false(r1$partial5)
  expect_equal(r1$color, unname(RiboPilot:::RRNA_ELEMENT_COLORS["18S"]))
  expect_length(p$rows[[2]]$features, 1L)        # NA start is skipped
})

test_that("the difference count includes gap columns and is zero for identical rows", {
  lab <- c(a = "a", b = "b", p1.s1 = "a", p2.s1 = "b")
  rows_of <- function(x) RiboPilot:::paths_payload(x, no_feats(), lab)$rows
  expect_equal(RiboPilot:::paths_diff_count(rows_of(aln_fixture())), 2L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGT", b = "ACGT")))), 0L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGTAC", b = "--GTAC")))), 2L)
})

test_that("paths_align aligns two near-identical sequences to equal width", {
  skip_if_not_installed("DECIPHER")
  a <- strrep("ACGTTGCAAC", 8)
  b <- paste0(substr(a, 1, 40), substr(a, 44, 80))
  m <- RiboPilot:::paths_align(c(`p1.s1` = a, `p2.s1` = b))
  expect_s4_class(m, "DNAStringSet")
  expect_equal(names(m), c("p1.s1", "p2.s1"))
  expect_equal(length(unique(Biostrings::width(m))), 1L)
  expect_null(RiboPilot:::paths_align(c(`p1.s1` = a)))
})
```

`paths_payload` must accept the zero-row data frame from `no_feats()`.

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: FAIL, `paths_payload` not found.

- [ ] **Step 3: Implement**

```r
# append to R/paths_align.R

#' Align the contigs of one set
#' @param seqs Named character vector (unit keys as names).
#' @return Aligned DNAStringSet or NULL.
#' @noRd
paths_align <- function(seqs) {
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  if (length(seqs) < 2) return(NULL)
  # ponytail: no cap on contig count; DECIPHER on a few 10 kb contigs takes seconds
  dna <- Biostrings::DNAStringSet(toupper(seqs))
  tryCatch(DECIPHER::AlignSeqs(dna, verbose = FALSE), error = function(e) NULL)
}

#' Viewer payload for an alignment: one row per contig, features in columns
#' @noRd
paths_payload <- function(aln, feats, labels) {
  seqs <- as.character(aln)
  rows <- lapply(names(seqs), function(u) {
    s <- seqs[[u]]
    cols <- which(strsplit(s, "", fixed = TRUE)[[1]] != "-")
    f <- feats[feats$unit == u & !is.na(feats$start) & !is.na(feats$stop), , drop = FALSE]
    list(unit = u, label = unname(labels[u]) %|NA|% u, seq = s,
      features = lapply(seq_len(nrow(f)), function(i) list(
        row = f$element[i], gene = f$element[i],
        pos1 = as.integer(cols[f$start[i]]), pos2 = as.integer(cols[f$stop[i]]),
        dir = if (identical(f$strand[i], "-")) "-" else "+",
        partial5 = isTRUE(f$partial_start[i] == 1), partial3 = isTRUE(f$partial_stop[i] == 1),
        color = unname(RRNA_ELEMENT_COLORS[f$element[i]] %|NA|% "#888888"))))
  })
  list(len = as.integer(nchar(seqs[1])), rows = rows)
}

#' Columns where the rows disagree (gaps count)
#' @noRd
paths_diff_count <- function(rows) {
  if (length(rows) < 2) return(0L)
  m <- do.call(rbind, lapply(rows, function(r) strsplit(r$seq, "", fixed = TRUE)[[1]]))
  sum(apply(m, 2, function(col) length(unique(col)) > 1))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: PASS (7 tests).

- [ ] **Step 5: Commit**

```bash
git add R/paths_align.R tests/testthat/test-paths-align.R
git commit -m "feat: align a contig set and build the viewer payload"
```

---

### Task 3: Sender and viewer UI mode

**Files:**
- Modify: `R/paths_align.R`
- Modify: `R/seqview.R:44-76` (`seqview_ui`)
- Test: `tests/testthat/test-paths-align.R`

**Interfaces:**
- Consumes: `paths_align`, `paths_payload` from Task 2.
- Produces:
  - `paths_send(session, id, set)`: `set` is one element of `paths_homology_sets()`. Queries the DB through `app_con()`, sends `rpseq` with `id = "seqv_paths-canvas"`, returns the payload or `NULL`.
  - `seqview_ui(prefix, summary, open = TRUE, pick = TRUE, mode = c("unit", "align"))`: in `align` mode emits buttons `data-rpseq` = `whole`, `zoom_in`, `zoom_out`, `prev_diff`, `next_diff`, `copy`, the goto input, and the `show_nt` box only; with `summary = NULL` returns a plain `div` (no `<details>`).

- [ ] **Step 1: Write the failing tests**

```r
# append to tests/testthat/test-paths-align.R
test_that("paths_send aligns the set's contigs and sends rows to the paths canvas", {
  skip_if_not_installed("withr"); skip_if_not_installed("DECIPHER")
  dir <- mk_app_project(); withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  a <- strrep("ACGTTGCAAC", 8); b <- paste0(substr(a, 1, 40), substr(a, 44, 80))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, triage = 'keep' WHERE path = 1 AND scaffold = 1", params = list(a))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, triage = 'drop' WHERE path = 2 AND scaffold = 1", params = list(b))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,30,'+'), ('S1',2,1,'18S',1,30,'+')")
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  set <- list(label = "18S (2 contigs)", units = data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  p <- RiboPilot:::paths_send(session, "S1", set)
  expect_equal(sent$id, "seqv_paths-canvas")
  expect_equal(sent$unit, "S1")
  expect_equal(vapply(sent$rows, `[[`, "", "label"), c("p1.s1 keep", "p2.s1 drop"))
  expect_equal(sent$rows[[1]]$features[[1]]$gene, "18S")
  expect_identical(p$rows, sent$rows)
})

test_that("paths_send returns NULL when fewer than two contigs have a sequence", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = NULL WHERE path = 2")
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  set <- list(label = "x", units = data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  expect_null(RiboPilot:::paths_send(session, "S1", set))
  expect_null(sent)
})

test_that("the viewer's align mode has difference buttons and no coverage boxes", {
  html <- paste(as.character(RiboPilot:::seqview_ui("seqv_paths", NULL, mode = "align")), collapse = " ")
  expect_match(html, 'data-rpseq="prev_diff"', fixed = TRUE)
  expect_match(html, 'data-rpseq="next_diff"', fixed = TRUE)
  expect_match(html, "seqv_paths-canvas", fixed = TRUE)
  expect_match(html, "seqv_paths-show_nt", fixed = TRUE)
  expect_false(grepl("show_cov", html, fixed = TRUE))
  expect_false(grepl('data-rpseq="fit"', html, fixed = TRUE))
  expect_false(grepl("<details", html, fixed = TRUE))
  unit <- paste(as.character(RiboPilot:::seqview_ui("seqv_asm", "Cov")), collapse = " ")
  expect_match(unit, "show_cov", fixed = TRUE)
  expect_false(grepl("prev_diff", unit, fixed = TRUE))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: FAIL, `paths_send` not found; align-mode assertions fail.

- [ ] **Step 3: Implement `paths_send`**

```r
# append to R/paths_align.R

#' Align one homologous set and send it to the paths viewer
#' @return The payload, or NULL when the set could not be aligned.
#' @noRd
paths_send <- function(session, id, set) {
  con <- app_con(); on.exit(DBI::dbDisconnect(con))
  u <- set$units
  keys <- sprintf("p%d.s%d", u$path, u$scaffold)
  seqs <- labels <- stats::setNames(rep(NA_character_, nrow(u)), keys)
  feats <- list()
  for (i in seq_len(nrow(u))) {
    a <- DBI::dbGetQuery(con, "SELECT sequence, triage FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, u$path[i], u$scaffold[i]))
    if (!nrow(a)) next
    seqs[i] <- a$sequence[1]
    labels[i] <- trimws(paste(keys[i], a$triage[1] %|NA|% ""))
    f <- DBI::dbGetQuery(con, "SELECT element, start, stop, strand, partial_start, partial_stop
      FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(id, u$path[i], u$scaffold[i]))
    if (nrow(f)) feats[[keys[i]]] <- cbind(unit = keys[i], f, stringsAsFactors = FALSE)
  }
  aln <- paths_align(seqs)
  if (is.null(aln)) return(NULL)
  feats <- if (length(feats)) do.call(rbind, feats) else
    data.frame(unit = character(), element = character(), start = numeric(), stop = numeric(),
      strand = character(), partial_start = numeric(), partial_stop = numeric())
  pl <- paths_payload(aln, feats, labels)
  pl$id <- "seqv_paths-canvas"; pl$unit <- id
  session$sendCustomMessage("rpseq", pl)
  pl
}
```

- [ ] **Step 4: Add `mode` to `seqview_ui`**

Replace the whole `seqview_ui` function in `R/seqview.R` with:

```r
seqview_ui <- function(prefix, summary, open = TRUE, pick = TRUE, mode = c("unit", "align")) {
  mode <- match.arg(mode)
  id <- function(x) paste0(prefix, "-", x)
  btn <- function(action, label, icon = NULL, title = NULL) {
    tags$button(type = "button", class = "btn btn-default btn-sm", `data-rpseq` = action,
      title = title, `aria-label` = title %||% label, icon, label)
  }
  goto <- tags$input(id = id("goto"), type = "number", min = 1, step = 1, class = "form-control input-sm",
    style = "width:130px;", placeholder = if (mode == "align") "Go to column" else "Go to position",
    `aria-label` = "Go to position")
  box <- function(x, label) tags$label(class = "rp-seqview-check",
    tags$input(type = "checkbox", id = id(x), checked = NA), label)
  controls <- if (mode == "align") div(class = "rp-seqview-controls",
      goto,
      btn("whole", "Whole alignment"),
      btn("zoom_in", NULL, icon("magnifying-glass-plus"), "Zoom in"),
      btn("zoom_out", NULL, icon("magnifying-glass-minus"), "Zoom out"),
      btn("prev_diff", NULL, icon("backward-step"), "Previous difference"),
      btn("next_diff", NULL, icon("forward-step"), "Next difference"),
      btn("copy", "Copy sequence", title = "Copy the rows in view as gapped FASTA"),
      box("show_nt", "Nucleotides"))
    else div(class = "rp-seqview-controls",
      goto,
      btn("fit", "Fit gene", title = "Zoom to the selected gene"),
      btn("whole", "Whole sequence"),
      btn("zoom_in", NULL, icon("magnifying-glass-plus"), "Zoom in"),
      btn("zoom_out", NULL, icon("magnifying-glass-minus"), "Zoom out"),
      btn("copy", "Copy sequence", title = "Copy the bases in view as FASTA"),
      box("show_cov", "Coverage"), box("show_clip", "Clipped reads"),
      box("show_ident", "BLAST identity"), box("show_origin", "BLAST origin"), box("show_nt", "Nucleotides"))
  help <- if (mode == "unit") tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      paste0("Drag or scroll sideways to pan, scroll to zoom", if (pick) "; click a gene to select it",
        ". Letters appear when zoomed in. Copy sequence copies the bases in view. "),
      sprintf("Clipped reads is the share of reads soft-clipped at each base (shown where depth is at least 10); bars turn orange at %d%%, the chimera check's cutoff.",
        round(100 * CHIMERA_CLIP_FRAC)),
      " BLAST origin colours each stretch by the class its BLAST windows match: grey is your lineage, white is no hit, ",
      "other colours are foreign; a dashed line marks a possible junction. BLAST identity is the best hit's identity in each window, coloured like the origin band. ",
      "Hover for the match.")
  body <- tagList(controls, help,
    div(class = "rp-seqview",
      tags$canvas(id = id("canvas"), class = "rp-seqview-canvas"),
      tags$div(id = id("tip"), class = "rp-seqview-tip", hidden = NA)))
  if (is.null(summary)) div(body) else rp_details(summary, open = open, body)
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align|seqview")'`
Expected: PASS. The existing `test-seqview.R` tests still pass (the unit mode markup is unchanged).

- [ ] **Step 6: Commit**

```bash
git add R/paths_align.R R/seqview.R tests/testthat/test-paths-align.R
git commit -m "feat: paths viewer sender and align-mode controls"
```

---

### Task 4: Modal section

**Files:**
- Modify: `R/app_modals.R:171-215` (`assemble_modal`), add `paths_section()` next to `rejected_contigs_section()`
- Test: `tests/testthat/test-paths-align.R`

**Interfaces:**
- Consumes: `paths_homology_sets` (Task 1), `seqview_ui(mode = "align")` (Task 3).
- Produces: `paths_section(id)` returning `NULL` or a `<details>` with `selectInput("paths_set")` (only for two or more sets), `uiOutput("paths_summary")`, and the `seqv_paths` viewer. Also `paths_sets_for(id)` returning `paths_homology_sets()` of the sample's annotations, used by the server in Task 5.

- [ ] **Step 1: Write the failing tests**

```r
# append to tests/testthat/test-paths-align.R
test_that("the Assemble window shows the paths section only for a multi-path sample with a set", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  html <- function() paste(as.character(RiboPilot:::assemble_modal("S1")), collapse = " ")
  expect_false(grepl("seqv_paths-canvas", html(), fixed = TRUE))   # no annotations yet
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,30,'+'), ('S1',2,1,'18S',1,30,'+')")
  h <- html()
  expect_match(h, "seqv_paths-canvas", fixed = TRUE)
  expect_match(h, "Compare assembly paths", fixed = TRUE)
  expect_false(grepl("paths_set", h, fixed = TRUE))                 # one set: no dropdown
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,2,'28S',1,30,'+'), ('S1',2,2,'28S',1,30,'+')")
  expect_match(html(), 'id="paths_set"', fixed = TRUE)
  DBI::dbExecute(con, "DELETE FROM assemblies WHERE path = 2")
  DBI::dbExecute(con, "DELETE FROM annotations WHERE path = 2")
  DBI::dbDisconnect(con)
  expect_false(grepl("seqv_paths-canvas", html(), fixed = TRUE))   # one path left
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: FAIL on `expect_match(h, "seqv_paths-canvas")`.

- [ ] **Step 3: Implement**

Add after `rejected_contigs_section()` in `R/app_modals.R`:

```r
#' The sample's homologous contig sets, from its annotations
#' @noRd
paths_sets_for <- function(id) {
  f <- app_query("SELECT path, scaffold, element FROM annotations WHERE ID = ?", params = list(id))
  paths_homology_sets(f)
}

#' Paths alignment section of the Assemble window, or nothing for a
#' single-path sample
#' @noRd
paths_section <- function(id) {
  n_paths <- app_query("SELECT COUNT(DISTINCT path) AS n FROM assemblies WHERE ID = ?", params = list(id))$n[1]
  if (is.na(n_paths) || n_paths < 2) return(NULL)
  sets <- paths_sets_for(id)
  if (!length(sets)) return(NULL)
  labels <- vapply(sets, `[[`, "", "label")
  picker <- if (length(sets) > 1) div(style = "max-width:320px; margin-bottom:4px;",
    selectInput("paths_set", NULL, choices = stats::setNames(seq_along(sets), labels), width = "100%"))
  rp_details("Compare assembly paths", open = TRUE,
    tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
      "Rows are contigs from each path aligned together. Red ticks on the Differences track mark ",
      "columns where the paths disagree; use the arrows to step through them. Letters appear when ",
      "zoomed in; bases that match in every row are greyed."),
    picker,
    uiOutput("paths_summary"),
    seqview_ui("seqv_paths", NULL, pick = FALSE, mode = "align"))
}
```

In `assemble_modal()`, insert `paths_section(id),` right after the Copy FASTA `div(...)` and before `tags$hr(style = "margin:10px 0;"),`:

```r
    div(style = "margin-top:6px;",
      actionButton("copy_fasta", "Copy FASTA", icon = icon("clipboard"), class = "btn-sm"),
      tags$span(class = "text-muted", style = "margin-left:8px; font-size:0.82em;",
        "Copies the ticked contigs as FASTA, headed by sample, taxon and contig.")),
    paths_section(id),
    tags$hr(style = "margin:10px 0;"),
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align|seqview|app-table")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app_modals.R tests/testthat/test-paths-align.R
git commit -m "feat: Compare assembly paths section in the Assemble window"
```

---

### Task 5: Server wiring with per-sample cache

**Files:**
- Modify: `R/app_server.R` (Assemble modal block, after the `observe` at line ~2368 that calls `seqview_send(session, "seqv_asm-canvas", ...)`)
- Test: `tests/testthat/test-paths-align.R`

**Interfaces:**
- Consumes: `paths_sets_for(id)` (Task 4), `paths_send(session, id, set)` (Task 3), `paths_diff_count(rows)` (Task 2).
- Produces: `paths_summary_text(pl)`: pure, `"Aligning contigs..."` is never its output; returns `"Could not align these contigs."` for `NULL`, else `"<n> contigs, <len> columns, <d> differ"` with thousands separators.

- [ ] **Step 1: Write the failing test**

```r
# append to tests/testthat/test-paths-align.R
test_that("the summary line counts rows, columns, and differences", {
  rows <- RiboPilot:::paths_payload(aln_fixture(), no_feats(), c(p1.s1 = "a", p2.s1 = "b"))
  expect_equal(RiboPilot:::paths_summary_text(rows), "2 contigs, 10 columns, 2 differ")
  expect_equal(RiboPilot:::paths_summary_text(NULL), "Could not align these contigs.")
  big <- list(len = 7812L, rows = rows$rows)
  expect_match(RiboPilot:::paths_summary_text(big), "7,812 columns", fixed = TRUE)
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align")'`
Expected: FAIL, `paths_summary_text` not found.

- [ ] **Step 3: Implement the pure helper**

Append to `R/paths_align.R`:

```r
#' Summary line above the paths viewer
#' @noRd
paths_summary_text <- function(pl) {
  if (is.null(pl)) return("Could not align these contigs.")
  sprintf("%d contigs, %s columns, %s differ", length(pl$rows),
    format(pl$len, big.mark = ","), format(paths_diff_count(pl$rows), big.mark = ","))
}
```

- [ ] **Step 4: Wire the server**

In `R/app_server.R`, directly after the block

```r
  observe({
    edit_tick(); refresh()
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    try(seqview_send(session, "seqv_asm-canvas", NULL, id, h$path, h$scaffold), silent = TRUE)
  })
```

add:

```r
  # ---- paths alignment: homologous contigs across paths, aligned once per set ----
  paths_cache <- new.env(parent = emptyenv())
  paths_result <- reactiveVal(NULL)
  observe({
    refresh()
    id <- sel_id(); req(id)
    # one sample's alignments at a time; the dropdown value can be stale from
    # the previous sample, so it only counts when this sample has a dropdown
    if (!identical(paths_cache$id, id)) { rm(list = ls(paths_cache), envir = paths_cache); paths_cache$id <- id; paths_result(NULL) }
    sets <- paths_sets_for(id)
    k <- if (length(sets) > 1) suppressWarnings(as.integer(input$paths_set %||% 1L)) else 1L
    if (is.na(k) || k > length(sets)) k <- 1L
    if (!length(sets)) { paths_result(NULL); return() }
    key <- paste("set", k)
    if (!exists(key, envir = paths_cache, inherits = FALSE)) {
      assign(key, try(paths_send(session, id, sets[[k]]), silent = TRUE), envir = paths_cache)
    }
    r <- get(key, envir = paths_cache)
    if (inherits(r, "try-error")) r <- NULL
    if (!is.null(r)) session$sendCustomMessage("rpseq", c(r, list(id = "seqv_paths-canvas", unit = id)))
    paths_result(list(id = id, pl = r))
  })
  output$paths_summary <- renderUI({
    r <- paths_result()
    txt <- if (is.null(r)) "Aligning contigs..." else paths_summary_text(r$pl)
    span(txt, style = "color:#555; display:block; margin-bottom:4px;")
  })
```

Notes for the implementer:
- `paths_send()` already sends the message the first time; the explicit `sendCustomMessage` on cache hits replays the payload when the modal is reopened for the same sample (the `observe` re-runs because `refresh()` ticks on open). Sending twice on the first run is harmless: the viewer reloads the same rows.
- `%||%` already exists in the package (used in `seqview_ui`).

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "paths-align|seqview|app-")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add R/app_server.R R/paths_align.R tests/testthat/test-paths-align.R
git commit -m "feat: align paths on open and cache per sample and set"
```

---

### Task 6: Viewer alignment mode (JS)

**Files:**
- Modify: `inst/app/www/seqviewer.js`
- Modify: `inst/app/www/custom.css` (only if the controls wrap badly; otherwise untouched)

**Interfaces:**
- Consumes: the `rpseq` message with `rows` (Task 3 payload shape).
- Produces: `window.rpseq.state(id)` gains `rows` and `diff` counts in align mode; `window.rpseq.diffAt(id, col)` (boolean) for driving checks.

No R test harness covers JS. Verification is by driving the app (Task 7). Work in these steps so each is small and reviewable.

- [ ] **Step 1: Pure helpers in the geometry block**

In the first IIFE (after `segOwns`), add and export:

```js
  // Columns where the rows disagree; a gap is a symbol like any other.
  function diffColumns(seqs) {
    var out = [], n = seqs.length ? seqs[0].length : 0;
    for (var c = 0; c < n; c++) {
      var ch = seqs[0].charAt(c);
      for (var r = 1; r < seqs.length; r++) if (seqs[r].charAt(c) !== ch) { out.push(c + 1); break; }
    }
    return out;
  }
  // Native position of column `col` in a gapped string: count of non-gap
  // characters up to and including it, or null when the column is a gap.
  function nativePos(seq, col) {
    if (seq.charAt(col - 1) === '-') return null;
    var n = 0;
    for (var i = 0; i < col; i++) if (seq.charAt(i) !== '-') n++;
    return n;
  }
```

and extend the export line:

```js
  window.rpseq.geom = { span: span, lanes: lanes, partialEdge: partialEdge, segOwns: segOwns,
                        diffColumns: diffColumns, nativePos: nativePos };
```

- [ ] **Step 2: Load rows**

In `Viewer.prototype.load`, after `var sameSeq = ...` insert alignment handling so the whole function reads:

```js
  Viewer.prototype.load = function (p) {
    var seq = p.rows ? (p.rows[0] ? p.rows[0].seq : '') : (p.seq || '');
    var sameSeq = this.seq === seq && !!this.rows === !!p.rows;
    this.len = p.len; this.seq = seq; this.input = p.input; this.unit = p.unit || 'unit';
    this.rows = p.rows || null;
    if (this.rows) {
      this.rows.forEach(function (r) {
        r.feats = (r.features || []).map(function (f) { return Object.assign({}, f); });
        r.nLanes = Math.max(1, G.lanes(r.feats));
      });
      this.feats = [];
      this.diff = G.diffColumns(this.rows.map(function (r) { return r.seq; }));
      this.diffSet = {}; this.diff.forEach(function (c) { this.diffSet[c] = true; }, this);
      this.nLanes = 0; this.depth = null; this.clip = null; this.origin = null; this.junctions = null; this.ident = null;
    } else {
      this.feats = (p.features || []).map(function (f) { return Object.assign({}, f); });
      this.nLanes = G.lanes(this.feats);
      this.depth = p.depth || null; this.clip = p.clip || null;
      this.origin = p.origin || null; this.junctions = p.junctions || null; this.ident = p.ident || null;
      this.diff = null; this.diffSet = null;
    }
    if (p.cutoff) this.clipFlag = p.cutoff;
    this.depthMax = (this.depth || []).reduce(function (m, d) { return d !== null && d > m ? d : m; }, 1);
    this.selected = (p.selected === undefined || p.selected === null) ? null : p.selected;
    if (!sameSeq) this.whole();
    this.draw();
  };
```

- [ ] **Step 3: Layout and draw in align mode**

Add constants next to the existing ones:

```js
  var DIFF_H = 14, ROW_GAP = 10, DIFF_RED = '#c0392b', MATCH_FILL = '#d9d9d9', MATCH_INK = '#888888';
```

Add row geometry helpers after `Viewer.prototype.height`:

```js
  // align mode: y of each row's first lane; the nt strip sits under its lanes
  Viewer.prototype.rowTop = function (i) {
    var y = RULER_H + PAD + DIFF_H + GAP;
    for (var k = 0; k < i; k++) y += this.rows[k].nLanes * LANE_H + (this.ntOn() ? NT_H : 0) + ROW_GAP;
    return y;
  };
  Viewer.prototype.rowNtY = function (i) { return this.rowTop(i) + this.rows[i].nLanes * LANE_H; };
  Viewer.prototype.alignHeight = function () {
    return this.rowTop(this.rows.length) - ROW_GAP + PAD;
  };
```

Change `height()` to:

```js
  Viewer.prototype.height = function () {
    if (this.rows) return this.alignHeight();
    return this.chimY() + (this.identOn() ? IDENT_H + GAP : 0) + (this.originOn() ? ORIGIN_H + GAP : 0);
  };
```

In `Viewer.prototype.draw`, right after `this.drawRuler(c);` add:

```js
    if (this.rows) { this.drawAlign(c, W); return; }
```

Add `drawAlign`, `drawDiffTrack`, and `drawRowNt`:

```js
  Viewer.prototype.drawAlign = function (c, W) {
    var y = RULER_H + PAD;
    this.drawDiffTrack(c, W, y);
    var self = this, H = this.alignHeight();
    // red wash over every differing column in view, behind all rows
    if (this.ppb >= NT_LETTER) {
      c.fillStyle = 'rgba(192, 57, 43, 0.12)';
      var vs = Math.floor(this.viewStart), ve = Math.ceil(this.viewStart + this.viewLen());
      for (var col = vs; col <= ve; col++) {
        if (!this.diffSet[col]) continue;
        var x = this.x(col); if (x < GUTTER) continue;
        c.fillRect(x, y + DIFF_H + GAP, Math.max(1, this.ppb), H - (y + DIFF_H + GAP));
      }
    }
    this.rows.forEach(function (r, i) {
      var top = self.rowTop(i);
      self.drawLanesAt(c, r.feats, top);
      if (self.ntOn()) self.drawRowNt(c, r, self.rowNtY(i));
      c.fillStyle = cssVar('--mp-text', '#333'); c.textAlign = 'right'; c.textBaseline = 'middle';
      c.fillText(r.label, GUTTER - 6, (self.ntOn() ? self.rowNtY(i) + NT_H / 2 : top + LANE_H / 2), GUTTER - 10);
    });
  };
  Viewer.prototype.drawDiffTrack = function (c, W, y) {
    c.fillStyle = cssVar('--mp-surface-alt', '#f5f5f5'); c.fillRect(GUTTER, y, W - GUTTER, DIFF_H);
    c.fillStyle = DIFF_RED;
    var w = Math.max(1, this.ppb), vs = this.viewStart, ve = vs + this.viewLen();
    for (var i = 0; i < this.diff.length; i++) {
      var col = this.diff[i]; if (col < vs - 1 || col > ve + 1) continue;
      var x = this.x(col); if (x < GUTTER || x > W) continue;
      c.fillRect(x, y, w, DIFF_H);
    }
    c.strokeStyle = cssVar('--mp-border', '#ccc'); c.beginPath(); c.moveTo(GUTTER, y + DIFF_H + 0.5); c.lineTo(W, y + DIFF_H + 0.5); c.stroke();
    c.fillStyle = cssVar('--mp-text-muted', '#6a6a6a'); c.textAlign = 'right'; c.textBaseline = 'middle';
    c.fillText('Differences', GUTTER - 6, y + DIFF_H / 2, GUTTER - 10);
    this.diffTrack = { y0: y, y1: y + DIFF_H };
  };
  Viewer.prototype.drawRowNt = function (c, r, y) {
    var vs = Math.floor(this.viewStart), ve = Math.ceil(this.viewStart + this.viewLen());
    var letter = this.ppb >= NT_LETTER, w = Math.max(1, this.ppb - (letter ? 1 : 0));
    c.textAlign = 'center'; c.textBaseline = 'middle';
    for (var col = vs; col <= ve; col++) {
      if (col < 1 || col > this.len) continue;
      var x = this.x(col); if (x < GUTTER) continue;
      var b = r.seq.charAt(col - 1) || 'N', same = letter && !this.diffSet[col];
      if (b === '-') {
        c.fillStyle = MATCH_INK; c.fillRect(x + w * 0.2, y + NT_H / 2 - 1, w * 0.6, 2);
        continue;
      }
      c.fillStyle = same ? MATCH_FILL : (BASE[b] || BASE.N) + (letter ? '99' : 'cc');
      c.fillRect(x, y + 2, w, NT_H - 4);
      if (letter) { c.fillStyle = same ? MATCH_INK : '#000000'; c.fillText(b, x + this.ppb / 2, y + NT_H / 2); }
    }
  };
```

Refactor `drawLanes` so the per-unit mode and the align mode share one body. Rename the existing function to `drawLanesAt(c, feats, top)` with the lane y computed as `top + f.lane * LANE_H` instead of `this.laneY(f.lane)`, and add:

```js
  Viewer.prototype.drawLanes = function (c) { this.drawLanesAt(c, this.feats, this.topH()); };
```

Inside `drawLanesAt`, the `this.feats.forEach(function (f) {` becomes `feats.forEach(function (f) {` and `var ... y = this.laneY(f.lane)` becomes `y = top + f.lane * LANE_H`. Everything else in the body stays.

- [ ] **Step 4: Navigation, click, tooltip, copy**

Add after `Viewer.prototype.fitRange`:

```js
  // nearest differing column outside the view on the given side; wraps
  Viewer.prototype.stepDiff = function (dir) {
    if (!this.diff || !this.diff.length) return;
    var vs = Math.ceil(this.viewStart), ve = Math.floor(this.viewStart + this.viewLen()), col = null;
    if (dir > 0) { for (var i = 0; i < this.diff.length; i++) if (this.diff[i] > ve) { col = this.diff[i]; break; } if (col === null) col = this.diff[0]; }
    else { for (var j = this.diff.length - 1; j >= 0; j--) if (this.diff[j] < vs) { col = this.diff[j]; break; } if (col === null) col = this.diff[this.diff.length - 1]; }
    this.goto(col);
  };
```

In `bindControls`, extend the button handler:

```js
        else if (a === 'prev_diff') self.stepDiff(-1); else if (a === 'next_diff') self.stepDiff(1);
```

In `Viewer.prototype.click`, before `var h = this.hitTest(px, py);` add:

```js
    if (this.rows && this.diffTrack && py >= this.diffTrack.y0 && py <= this.diffTrack.y1 && px >= GUTTER) {
      this.goto(Math.floor(this.viewStart + (px - GUTTER) / this.ppb)); return;
    }
```

In the canvas `mousemove` tooltip handler, right after `if (!b || px < GUTTER) { self.tip.hidden = true; return; }` add:

```js
      if (self.rows) {
        var t2 = 'Column ' + b.pos.toLocaleString();
        self.rows.forEach(function (r) {
          var np = G.nativePos(r.seq, b.pos), ch = r.seq.charAt(b.pos - 1);
          t2 += ' | ' + r.label.split(' ')[0] + ' ' + (np === null ? 'gap' : np.toLocaleString() + ' ' + ch);
        });
        if (h) t2 += ' | ' + h.f.gene + ' ' + h.f.pos1.toLocaleString() + '-' + h.f.pos2.toLocaleString() + ' (' + h.f.dir + ')';
        self.tip.textContent = t2; self.tip.hidden = false;
        var tw2 = self.tip.offsetWidth, th2 = self.tip.offsetHeight;
        self.tip.style.left = Math.max(0, Math.min(px + 12, cv.clientWidth - tw2)) + 'px';
        self.tip.style.top = (py + 12 + th2 > cv.clientHeight ? py - th2 - 8 : py + 12) + 'px';
        return;
      }
```

Change `visibleFasta` to handle rows:

```js
  Viewer.prototype.visibleFasta = function () {
    var a = Math.max(1, Math.ceil(this.viewStart)), b = Math.min(this.len, Math.floor(this.viewStart + this.viewLen()));
    var fa = function (name, s) { return '>' + name + ':' + a + '-' + b + '\n' + (s.match(/.{1,60}/g) || []).join('\n') + '\n'; };
    if (this.rows) {
      var n = 0, text = this.rows.map(function (r) { var s = r.seq.substring(a - 1, b); n += s.length; return fa(r.label.replace(/\s+/g, '_'), s); }).join('');
      return { n: n, text: text };
    }
    var s = this.seq.substring(a - 1, b);
    return { n: s.length, text: fa(this.unit, s) };
  };
```

Extend `rpseq.state` and add `diffAt`:

```js
  window.rpseq.state = function (id) {
    var v = viewers[id]; if (!v) return null;
    return { len: v.len, viewStart: v.viewStart, ppb: v.ppb, selected: v.selected, height: v.height(), nLanes: v.nLanes,
             origin: (v.origin || []).length, junctions: (v.junctions || []).length,
             rows: (v.rows || []).length, diff: (v.diff || []).length,
             features: v.feats.map(function (f) { return { row: f.row, gene: f.gene, pos1: f.pos1, pos2: f.pos2, lane: f.lane }; }) };
  };
  window.rpseq.diffAt = function (id, col) { var v = viewers[id]; return !!(v && v.diffSet && v.diffSet[col]); };
```

- [ ] **Step 5: Syntax check**

Run: `node --check inst/app/www/seqviewer.js` if `node` is on PATH; otherwise open the app (Task 7) and check the browser console for errors.

- [ ] **Step 6: Run the R suite for regressions**

Run: `Rscript -e 'testthat::test_local(filter = "seqview|paths-align")'`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add inst/app/www/seqviewer.js
git commit -m "feat: alignment mode in the sequence viewer"
```

---

### Task 7: Drive the app and check both viewers

**Files:**
- None modified unless a defect is found.

Follow the project's headless drive approach (memory: chromote with flatpak Chrome, `headless = FALSE`, short `TMPDIR`). Use the test project from `init_test_project()` or any project with a multi-path sample.

- [ ] **Step 1: Open an Assemble details window for a sample with two or more paths**

Check: the "Compare assembly paths" section is present, open, and the summary reads "N contigs, M columns, D differ" after a moment.

- [ ] **Step 2: Inspect the viewer**

Run in the page: `rpseq.state('seqv_paths-canvas')`. Expect `rows >= 2`, `diff >= 0`, `height > 0`. Click "Next difference" and confirm the view jumps and letters appear with matching columns greyed. Hover and confirm the tooltip lists each row's native position.

- [ ] **Step 3: Check the per-contig viewer still works**

Click a contig row, expand "Coverage and windowed BLAST for p1.s1", run `rpseq.state('seqv_asm-canvas')` and confirm `rows === 0` and the coverage track draws as before.

- [ ] **Step 4: Open a single-path sample**

Confirm no "Compare assembly paths" section.

- [ ] **Step 5: Fix anything found, re-run the full suite, commit**

Run: `Rscript -e 'testthat::test_local()'`
Expected: PASS (the known-flaky remote BLAST tests aside).

```bash
git add -A
git commit -m "fix: paths alignment drive fixes"
```

(Skip the commit when nothing changed.)

---

### Task 8: Docs and NEWS

**Files:**
- Modify: `vignettes/Test-Project-Assemble.Rmd` (section describing the details window)
- Modify: `NEWS.md` (top unreleased section)

- [ ] **Step 1: Add a short paragraph to the Assemble vignette**

After the paragraph that explains clicking a contig row for its evidence, add:

```markdown
When a sample has more than one assembly path, the details window also shows
**Compare assembly paths**: the contigs that cover the same part of the rDNA
repeat, aligned together. Red ticks on the Differences track mark the columns
where the paths disagree; use the arrow buttons to step through them. Zoom in
to see the bases, with the matching ones greyed out. Each row shows its own
annotated features, so you can see whether a difference falls inside a gene.
```

- [ ] **Step 2: Add a NEWS entry**

Under the current unreleased heading in `NEWS.md`:

```markdown
- The Assemble details window aligns a sample's assembly paths and shows the
  alignment in the sequence viewer, with differing columns marked.
```

- [ ] **Step 3: Commit**

```bash
git add vignettes/Test-Project-Assemble.Rmd NEWS.md
git commit -m "docs: Compare assembly paths section"
```
