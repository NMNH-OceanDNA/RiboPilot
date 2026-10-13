# Chimera check in the sequence viewer: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the chimera check's lineage segments, junctions, and window identity out of the static ggplot into the canvas sequence viewer, make the chimera tables click-to-zoom, and show the verdict inline.

**Architecture:** R builds one pure payload (`chimera_origin()`) from the stored segments, windows, and junctions; `seqview_send()` adds it to the message the viewer already receives. `seqviewer.js` draws two more tracks (origin band, % identity) and junction lines, with hover and a `rpseq_goto` message for table clicks. `plot_chimera()` and `output$chimera_plot` are deleted. The chimera section keeps the pill, an inline reason, the segments table (windows as expandable rows), and the junctions table.

**Tech Stack:** R, Shiny, reactable, testthat; hand-rolled canvas JS in `inst/app/www/seqviewer.js` (no JS test harness; the app is verified by driving it with chromote).

**Spec:** `dev/review/2026-10-13/CHIMERA_UI_REVIEW.md` (sections "Addendum" and "Labels and tables after the merge").

## Global Constraints

- Branch off `devel`: `git switch devel && git switch -c feat/chimera-seqview`. Never push.
- No new dependencies. No non-ASCII in R code. Minimal comments.
- Keep `chimera_seg_colors()` as the single place that picks segment colours.
- Coordinates sent to the viewer are the layer's own (review copy on Annotation Review, original on Assemble), exactly as `seqview_send()` already does for genes and depth.
- Trim modal junction lines / presets and the action-button verdict card are NOT in this plan (follow-up).
- Commit messages: brief, no attribution lines beyond the session link the harness adds.

## Review Focus

1. Contig with a chimera row but `segments` NULL (clean call, no windows stored): viewer must draw no origin band and no error. Test in Task 1.
2. Overlapping segments (200 bp windows at step 100 always overlap by one window): the band must show a clean edge, tables keep true coordinates. Test in Task 1.
3. Annotation Review right after a trim (checks not done): the viewer gets no origin data, the pill says "not run", nothing errors. Test in Task 2 and Task 5.
4. Junction exactly at a segment edge or at position 1 / `len`: dashed line still drawn inside the canvas; JS clamps x to the gutter. Covered by Task 3 manual check.
5. Clicking a table row when the viewer section is collapsed: the section opens and the viewer zooms. Covered by Task 4 JS (`section.open = true`) and the Task 7 drive.

---

### Task 1: `chimera_origin()` payload builder and the 3-colour palette

**Files:**
- Modify: `R/app_modals.R:1003-1023` (palette, `chimera_seg_colors`)
- Create: nothing new; add `chimera_origin()` to `R/seqview.R` after `seqview_payload()`
- Test: `tests/testthat/test-seqview.R`

**Interfaces:**
- Produces: `chimera_origin(seg, win, junctions, label, channels, sample_lt, len)` returning `list(origin = <list of lists>, junctions = <list of lists>, ident = <numeric vector length len or NULL>)`. Each origin item: `pos1, pos2, color, taxon, species, pident, host` (host = logical, own lineage). Each junction item: `pos, severe, checks` (checks = comma string from `chimera_junction_table()`).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-seqview.R`:

```r
test_that("chimera_origin de-overlaps segments for drawing and keeps the host grey", {
  seg <- data.frame(start = c(1, 501, 2251), stop = c(750, 2500, 4500),
    taxon = c("nohit", "Actinopterygii", "Eurotiomycetes"), species = c(NA, "Salmo salar", "Penicillium sp."),
    pident = c(NA, 99, 97), n_windows = c(3, 10, 9),
    lineage_taxids = c(NA, "2759;33208;7711;7898;186623", "2759;4751;4890;147545;451871"),
    stringsAsFactors = FALSE)
  win <- data.frame(start = c(1, 2251), stop = c(500, 2750), pident = c(NA, 97))
  o <- RiboPilot:::chimera_origin(seg, win, "2375", "severe", NA, "2759;33208;7711;7898;186623", 4500L)
  expect_equal(vapply(o$origin, `[[`, 0, "pos2"), c(500, 2250, 4500))
  expect_true(o$origin[[2]]$host); expect_false(o$origin[[3]]$host)
  expect_equal(o$origin[[2]]$color, "#B3B3B3")
  expect_equal(o$origin[[1]]$color, "#FFFFFF")
  expect_equal(o$origin[[3]]$color, RiboPilot:::CHIMERA_LINEAGE_COLORS[1])
  expect_equal(o$junctions[[1]]$pos, 2375L); expect_true(o$junctions[[1]]$severe)
  expect_length(o$ident, 4500L)
  expect_true(is.na(o$ident[1])); expect_equal(o$ident[2300], 97)
})

test_that("chimera_origin with no segments or junctions sends nothing", {
  o <- RiboPilot:::chimera_origin(data.frame(), data.frame(), NA, "clean", NA, NA_character_, 100L)
  expect_length(o$origin, 0L); expect_length(o$junctions, 0L); expect_null(o$ident)
})

test_that("the lineage palette has three colours", {
  expect_equal(RiboPilot:::CHIMERA_LINEAGE_COLORS, c("#EEDD88", "#77AADD", "#AAAA00"))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: 3 failures, "could not find function chimera_origin" and palette mismatch.

- [ ] **Step 3: Implement**

In `R/app_modals.R` replace the palette line:

```r
CHIMERA_LINEAGE_COLORS <- c("#EEDD88", "#77AADD", "#AAAA00")
```

Append to `R/seqview.R`:

```r
#' Chimera evidence for the viewer: de-overlapped origin segments, junctions,
#' and per-base window identity. Pure.
#' @noRd
chimera_origin <- function(seg, win, junctions, label, channels, sample_lt, len) {
  out <- list(origin = list(), junctions = list(), ident = NULL)
  if (NROW(seg)) {
    seg <- seg[order(seg$start), , drop = FALSE]
    seg$col <- chimera_seg_colors(seg, sample_lt)
    n <- nrow(seg)
    stop <- pmin(seg$stop, c(seg$start[-1] - 1L, Inf))
    out$origin <- lapply(seq_len(n), function(i) list(
      pos1 = as.integer(seg$start[i]), pos2 = as.integer(stop[i]), color = seg$col[i],
      taxon = seg$taxon[i], species = seg$species[i] %|NA|% "", pident = seg$pident[i] %|NA|% NA,
      host = identical(seg$col[i], "#B3B3B3")))
  }
  if (NROW(win) && any(!is.na(win$pident))) {
    id <- rep(NA_real_, len)
    for (i in seq_len(nrow(win))) if (!is.na(win$pident[i])) {
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
```

Note `chimera_junction_table(NA, NA)` returns a zero-row frame, so a NULL/NA junctions string is safe.

- [ ] **Step 4: Run tests**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: all pass. Also run `test-chimera.R`; the `chimera_seg_colors` behaviour is unchanged apart from the palette.

- [ ] **Step 5: Commit**

```bash
git add R/seqview.R R/app_modals.R tests/testthat/test-seqview.R
git commit -m "Chimera origin payload and 3-colour lineage palette"
```

---

### Task 2: `seqview_send()` carries the chimera evidence

**Files:**
- Modify: `R/seqview.R` (`seqview_payload`, `seqview_send`)
- Test: `tests/testthat/test-seqview.R`

**Interfaces:**
- Consumes: `chimera_origin()` from Task 1; `review_unit_state()` (R/review.R:428).
- Produces: message fields `origin`, `junctions`, `ident` on the `rpseq` payload. Absent (not NULL-valued keys; simply missing) when there is no chimera row or the review checks are not done.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-seqview.R`:

```r
test_that("the viewer gets the chimera evidence of its own layer, and none before a recheck", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 100 WHERE path = 1 AND scaffold = 1",
    params = list(strrep("ACGTACGTAC", 10)))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions, segments) VALUES ('S1',1,1,'severe','50',
    '[{\"start\":1,\"stop\":60,\"taxon\":\"Actinopterygii\",\"species\":\"Salmo salar\",\"pident\":99,\"n_windows\":3},
      {\"start\":41,\"stop\":100,\"taxon\":\"Eurotiomycetes\",\"species\":\"Penicillium sp.\",\"pident\":97,\"n_windows\":3}]')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident) VALUES ('S1',1,1,1,20,99), ('S1',1,1,41,60,97)")
  fixture_mark_done(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)

  RiboPilot:::seqview_send(session, "seqv_asm-canvas", NULL, "S1", 1L, 1L)
  expect_equal(vapply(sent$origin, `[[`, 0, "pos2"), c(40, 100))
  expect_equal(sent$junctions[[1]]$pos, 50L)
  expect_equal(sent$ident[50], 97)

  RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 90)
  RiboPilot:::seqview_send(session, "seqv_ann-canvas", NULL, "S1", 1L, 1L, review = TRUE)
  expect_null(sent$origin); expect_null(sent$junctions); expect_null(sent$ident)

  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label, junctions) VALUES ('S1',1,1,'minor','30')")
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2")
  RiboPilot:::seqview_send(session, "seqv_ann-canvas", NULL, "S1", 1L, 1L, review = TRUE)
  expect_equal(sent$junctions[[1]]$pos, 30L); expect_false(sent$junctions[[1]]$severe)
  expect_null(sent$origin)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-seqview.R")'`
Expected: FAIL, `sent$origin` is NULL on the first block.

- [ ] **Step 3: Implement**

In `seqview_payload()` add `chim = NULL` as the last parameter and, before `out`'s return, `out[names(chim)] <- chim` (only keys present in `chim` land in the payload). In `seqview_send()`, after `pl <- seqview_payload(...)` becomes:

```r
  chim <- seqview_chimera(con, id, path, scaffold, review, nchar(s$sequence[1]))
  pl <- seqview_payload(s$sequence[1], f, tr$depths, tr$clips, selected, chim)
```

Add to `R/seqview.R`:

```r
#' Chimera evidence of one unit's layer, or an empty list when the Annotation
#' Review checks are not done or there is no chimera row
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
```

`Filter(length, o)` drops empty `origin`/`junctions` lists and a NULL `ident`, so the JS sees missing keys, not empty arrays.

- [ ] **Step 4: Run tests**

Run the seqview and chimera test files. Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add R/seqview.R tests/testthat/test-seqview.R
git commit -m "Send chimera origin, junctions, and identity to the viewer"
```

---

### Task 3: Viewer draws origin band, identity track, junctions, hover

**Files:**
- Modify: `inst/app/www/seqviewer.js` (constants line 44; `load` 75-87; `topH` 134-136; `draw` 141-157; tooltip 321-334; checkbox wiring 348-351)
- Modify: `R/seqview.R` `seqview_ui()` (two more checkboxes, help sentence)

**Interfaces:**
- Consumes: payload keys `origin` (array of `{pos1,pos2,color,taxon,species,pident,host}`), `junctions` (array of `{pos,severe,checks}`), `ident` (array of numbers or null per base).
- Produces: `Viewer.prototype.fitRange(pos1, pos2)`; checkboxes `<prefix>-show_origin`, `<prefix>-show_ident`.

- [ ] **Step 1: UI checkboxes**

In `seqview_ui()` change the controls line to:

```r
      box("show_cov", "Coverage"), box("show_clip", "Clipped reads"),
      box("show_ident", "Window identity"), box("show_origin", "Origin"), box("show_nt", "Nucleotides")),
```

and extend the help paragraph with one sentence after the clipped-reads one:

```r
      " Origin colours each stretch by the class its BLAST windows match: grey is your lineage, white is no hit, ",
      "other colours are foreign; a dashed line marks a possible junction. Hover for the match."
```

- [ ] **Step 2: JS constants and state**

Line 44, add `ORIGIN_H = 18, IDENT_H = 36` to the var list. In `load()` after the `clip` line add:

```js
    this.origin = p.origin || null; this.junctions = p.junctions || null; this.ident = p.ident || null;
```

Add after `clipOn`:

```js
  Viewer.prototype.identOn = function () { return this.showIdent && !!this.ident; };
  Viewer.prototype.originOn = function () { return this.showOrigin && !!this.origin; };
```

`topH()` becomes:

```js
  Viewer.prototype.topH = function () {
    return RULER_H + PAD + (this.covOn() ? COV_H + GAP : 0) + (this.clipOn() ? CLIP_H + GAP : 0)
      + (this.identOn() ? IDENT_H + GAP : 0) + (this.originOn() ? ORIGIN_H + GAP : 0);
  };
```

In `draw()`, after the clip line and before `drawLanes`:

```js
    if (this.identOn()) { this.drawTrack(c, W, y, IDENT_H, this.ident, 100, '% id', false); y += IDENT_H + GAP; }
    if (this.originOn()) { this.drawOrigin(c, W, y); y += ORIGIN_H + GAP; }
    this.drawLanes(c);
    if (this.ntOn()) this.drawNt(c, this.topH() + this.nLanes * LANE_H + PAD);
    this.drawJunctions(c, H);
```

Checkbox wiring (line 348) gains two pairs: `['show_ident', 'showIdent'], ['show_origin', 'showOrigin']`.

- [ ] **Step 3: drawOrigin and drawJunctions**

Add after `drawTrack`:

```js
  // One rect per segment in the segment's colour, labelled when the label
  // fits the visible part (same rule as the gene lanes).
  Viewer.prototype.drawOrigin = function (c, W, y) {
    var muted = cssVar('--mp-text-muted', '#6a6a6a');
    c.textAlign = 'center'; c.textBaseline = 'middle';
    this.origin.forEach(function (s) {
      if (!this.inView(s)) return;
      var x0 = Math.max(GUTTER, this.x(s.pos1)), x1 = Math.min(W, this.x(s.pos2 + 1));
      if (x1 - x0 < 1) return;
      c.fillStyle = s.color; c.fillRect(x0, y, x1 - x0, ORIGIN_H);
      c.strokeStyle = cssVar('--mp-border', '#ccc'); c.strokeRect(x0 + 0.5, y + 0.5, x1 - x0 - 1, ORIGIN_H - 1);
      var full = s.host ? s.taxon + ' (your lineage)' : s.pident ? s.taxon + ', ' + s.species + ' ' + s.pident.toFixed(1) + '%' : s.taxon;
      var lab = x1 - x0 > c.measureText(full).width + 8 ? full : x1 - x0 > c.measureText(s.taxon).width + 8 ? s.taxon : null;
      if (lab) { c.fillStyle = cssVar('--mp-text', '#333'); c.fillText(lab, (x0 + x1) / 2, y + ORIGIN_H / 2); }
    }, this);
    c.fillStyle = muted; c.textAlign = 'right'; c.fillText('origin', GUTTER - 6, y + ORIGIN_H / 2, GUTTER - 10);
  };
  Viewer.prototype.drawJunctions = function (c, H) {
    if (!this.junctions) return;
    c.save(); c.setLineDash([5, 3]); c.lineWidth = 1.5;
    this.junctions.forEach(function (j) {
      var x = this.x(j.pos); if (x < GUTTER || x > this.wrap.clientWidth) return;
      c.strokeStyle = j.severe ? '#c0392b' : '#b9770e';
      c.beginPath(); c.moveTo(x, RULER_H); c.lineTo(x, H); c.stroke();
    }, this);
    c.restore();
  };
  Viewer.prototype.originAt = function (pos) {
    return (this.origin || []).find(function (s) { return pos >= s.pos1 && pos <= s.pos2; }) || null;
  };
```

`inView(f)` works on any object with `pos1`/`pos2`.

- [ ] **Step 4: Tooltip**

In the canvas `mousemove` handler, after the `h` line:

```js
      var o = self.originAt(b.pos);
      if (o) t += ' | ' + (o.host ? o.taxon + ' (your lineage)' : o.taxon === 'nohit' ? 'no hit' : o.taxon + (o.species ? ', ' + o.species : '') + (o.pident ? ' ' + o.pident.toFixed(1) + '%' : ''));
      var j = (self.junctions || []).find(function (q) { return Math.abs(self.x(q.pos) - px) <= 6; });
      if (j) t += ' | junction ' + j.pos.toLocaleString() + (j.checks ? ' (' + j.checks + ')' : '');
```

- [ ] **Step 5: fitRange and the `rpseq_goto` handler**

Replace `fit()` body so it delegates:

```js
  Viewer.prototype.fit = function (row) {
    var f = this.feats.find(function (x) { return x.row === row; });
    if (!f) return;
    this.selected = row; this.fitRange(f.pos1, f.pos2);
  };
  Viewer.prototype.fitRange = function (pos1, pos2) {
    this.isWhole = false;
    var s = pos2 - pos1 + 1, margin = Math.max(3, s * 0.05);
    this.ppb = Math.min(MAX_PPB, this.width() / (s + 2 * margin));
    this.viewStart = pos1 - margin; this.clamp();
    if (this.section) this.section.open = true;
    this.draw();
  };
```

Next to the `rpseq_select` handler add:

```js
    window.Shiny.addCustomMessageHandler('rpseq_goto', function (p) {
      withViewer(p.id, function (v) { v.fitRange(p.pos1, p.pos2); });
    });
```

- [ ] **Step 6: Syntax check**

Run: `node --check inst/app/www/seqviewer.js` (if node is absent, skip; the Task 7 drive catches errors in the browser console).

- [ ] **Step 7: Commit**

```bash
git add inst/app/www/seqviewer.js R/seqview.R
git commit -m "Viewer: origin band, window identity track, junction lines"
```

---

### Task 4: Tables click to zoom; windows fold into the segments table

**Files:**
- Modify: `R/utils_reactable.R` (add `rt_zoom_row()` after `rt_examine_row`, line ~762)
- Modify: `R/app_server.R:2714-2733` (`chimera_segments`), `2785-2795` (`chimera_junctions`), `2796-2833` (`chimera_windows` becomes a helper used as `details`)
- Modify: `R/app_modals.R:149-170` (`chimera_panel`: drop the "All windows" subsection)
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Consumes: `rpseq_goto` message from Task 3.
- Produces: Shiny input `chim_zoom` = `list(pos1, pos2, nonce)`; `rt_zoom_row(from, to)` returning a reactable `onClick` JS.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-chimera.R`:

```r
test_that("clicking a chimera table row zooms the open modal's viewer", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  msgs <- list()
  shiny::testServer(RiboPilot:::app_server, {
    session$userData$msgs <- NULL
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    local_mocked_bindings(seqview_goto = function(session, canvas, pos1, pos2) msgs[[length(msgs) + 1]] <<- list(canvas, pos1, pos2))
    session$setInputs(chim_zoom = list(pos1 = 100, pos2 = 300, nonce = 1))
    session$setInputs(mode = "Assemble")
    session$setInputs(chim_zoom = list(pos1 = 2625, pos2 = 2625, nonce = 2))
  })
  expect_equal(msgs[[1]], list("seqv_ann-canvas", 100L, 300L))
  expect_equal(msgs[[2]], list("seqv_asm-canvas", 2625L, 2625L))
  js <- as.character(RiboPilot:::rt_zoom_row("start", "stop"))
  expect_match(js, "chim_zoom", fixed = TRUE)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-chimera.R")'`
Expected: FAIL, `seqview_goto` / `rt_zoom_row` not found.

- [ ] **Step 3: Implement**

`R/utils_reactable.R`, after `rt_examine_row`:

```r
#' Row click that zooms the open modal's sequence viewer to from..to
#' @noRd
rt_zoom_row <- function(from, to) htmlwidgets::JS(sprintf(
  "function(rowInfo, column) {
    if (column.id === '.details') return;
    Shiny.setInputValue('chim_zoom', {pos1: rowInfo.row['%s'], pos2: rowInfo.row['%s'], nonce: Date.now()}, {priority: 'event'});
  }", from, to))
```

`R/seqview.R`:

```r
#' Zoom one viewer to a range
#' @noRd
seqview_goto <- function(session, canvas, pos1, pos2)
  session$sendCustomMessage("rpseq_goto", list(id = canvas, pos1 = as.integer(pos1), pos2 = as.integer(pos2)))
```

`R/app_server.R`, next to the `chim_unit` reactive:

```r
  observeEvent(input$chim_zoom, {
    z <- input$chim_zoom
    seqview_goto(session, if (identical(input$mode, "Annotate")) "seqv_ann-canvas" else "seqv_asm-canvas",
      z$pos1, z$pos2)
  })
```

Segments table (`output$chimera_segments`): keep the current column list, add `onClick = rt_zoom_row("start", "stop")`, and add a `details` that renders the segment's windows. Extract the existing windows-table body into a function in app_server.R (same file, above the output):

```r
  chimera_windows_tbl <- function(w) {
    w$lineage_short <- vapply(w$lineage, function(l) {
      if (is.na(l)) return(NA_character_)
      p <- strsplit(l, ";", fixed = TRUE)[[1]]; p <- p[nzchar(p)]
      paste(utils::tail(p, 3), collapse = " > ")
    }, character(1), USE.NAMES = FALSE)
    w <- w[, c("start", "stop", "species", "lineage_short", "pident", "bitscore", "rank_distance", "depth_median", "clip_frac", "lineage")]
    reactable::reactable(w, compact = TRUE, pagination = FALSE, onClick = rt_zoom_row("start", "stop"),
      defaultColDef = reactable::colDef(na = "", minWidth = 70, headerStyle = list(whiteSpace = "nowrap")),
      columns = <the existing `columns = list(...)` of output$chimera_windows, unchanged>)
  }
```

Then in `output$chimera_segments`, query the windows once (`chim_query("start, stop, species, lineage, pident, bitscore, rank_distance, depth_median, clip_frac", "chimera_windows", " ORDER BY start")`) and pass:

```r
      onClick = rt_zoom_row("start", "stop"),
      details = function(index) {
        ww <- w[findInterval(w$start, s$start) == index, , drop = FALSE]
        if (!nrow(ww)) return(NULL)
        chimera_windows_tbl(ww)
      },
```

Delete `output$chimera_windows` entirely. Junctions table: add `onClick = rt_zoom_row("position", "position")`.

`chimera_panel()`: delete the `rp_details("All windows", ...)` subsection and the surrounding `div(class = "rp-sub-list", ...)` wrapper if only one subsection remains (keep the "Possible junctions" details as is). Update the segments help: add `tags$p(class = "mp-help-text", "Click a row to zoom the viewer to it; expand a row for its windows.")` above the segments table.

- [ ] **Step 4: Run tests**

Run `test-chimera.R`, `test-seqview.R`, `test-batch-trim.R`. Expected: the new test passes; any test that referenced `output$chimera_windows` is updated to read `output$chimera_segments` instead (grep `chimera_windows` in tests first).

- [ ] **Step 5: Commit**

```bash
git add R/utils_reactable.R R/seqview.R R/app_server.R R/app_modals.R tests/testthat/test-chimera.R
git commit -m "Chimera tables: click to zoom, windows nested under segments"
```

---

### Task 5: Verdict inline, "not run" state, retire the static plot

**Files:**
- Modify: `R/app_modals.R:149-170` (`chimera_panel`), `1025-1118` (delete `chimera_result_text`'s hover form, delete `plot_chimera`), `917` (delete `fixed_gutter` if no other caller: `grep -n fixed_gutter R/`), `971` (comment mentions plot_chimera; reword to "same x extent as the viewer")
- Modify: `R/app_server.R:2700-2713` (delete `output$chimera_plot`, `output$chimera_db_source`; rewrite `output$chimera_result`)
- Modify: `tests/testthat/test-chimera.R:435-460, 492-501, 592-700`; `tests/testthat/test-batch-trim.R:250-262`

**Interfaces:**
- Produces: `chimera_result_text(label, reason, state = "done")`: pill plus the reason as visible help text; when `state != "done"`, a grey `mp-pill-default` "not run" pill with the state note.

- [ ] **Step 1: Rewrite the tests**

Replace the `chimera_result_text` test (test-chimera.R:492-501) with:

```r
test_that("the chimera panel shows the level and the reason in the open, and not run before a recheck", {
  expect_match(as.character(RiboPilot:::chimera_panel()), "chimera_result", fixed = TRUE)
  expect_false(grepl("chimera_plot", as.character(RiboPilot:::chimera_panel()), fixed = TRUE))
  x <- as.character(RiboPilot:::chimera_result_text("minor", "18S resolves to Chordata: agree only to phylum"))
  expect_match(x, '<span class="mp-pill mp-pill-warning">minor</span>', fixed = TRUE)
  expect_match(x, "agree only to phylum", fixed = TRUE)
  expect_false(grepl("title=", x, fixed = TRUE))
  x <- as.character(RiboPilot:::chimera_result_text(NA, NA, state = "not run"))
  expect_match(x, "not run", fixed = TRUE)
})
```

Delete the two `plot_chimera` tests at 435-460 and the state test at 592-650 (the viewer test in Task 2 covers the layer/state logic). In the panel test at 652-700: remove the `local_mocked_bindings(plot_chimera = ...)`, the `output$chimera_plot` lines, the `expect_true(drew)` / `expect_false(drew)` lines, and the two `chimera_db_source` expectations; the "not run" block changes `expect_error(output$chimera_result)` to `expect_match(output$chimera_result$html, "not run", fixed = TRUE)`.

In test-batch-trim.R:250-262, the test "the chimera plot redraws after a trim" becomes "the viewer is re-sent after a trim": mock `seqview_send = function(...) { drawn <<- drawn + 1L; invisible() }` and keep the rest of the test's assertions on `drawn`.

- [ ] **Step 2: Run to verify failures**

Run `test-chimera.R` and `test-batch-trim.R`. Expected: failures on the `title=` absence, the `state` argument, and the panel's `chimera_plot` presence.

- [ ] **Step 3: Implement**

`chimera_result_text`:

```r
chimera_result_text <- function(label, reason, state = "done") {
  if (!identical(state, "done")) return(tags$span("Chimera call: ",
    tags$span(class = "mp-pill mp-pill-default", "not run"),
    tags$span(class = "mp-help-text", style = "font-weight:normal;",
      if (identical(state, "waiting")) " (waiting for the Assemble step's Update)" else " (press Update to run the checks)")))
  tone <- switch(label, clean = "success", severe = "danger", "warning")
  tagList(tags$span("Chimera call: ", tags$span(class = paste0("mp-pill mp-pill-", tone), label)),
    if (!is.na(reason) && nzchar(reason)) tags$div(class = "mp-help-text", style = "font-weight:normal; margin-top:2px;", reason))
}
```

Check `mp-pill-default` exists in the CSS (`grep -n "mp-pill-" inst/app/www/*.css`); if not, use the class the "not run" badge in `rt_gated()` uses (`R/utils_reactable.R:274-290`).

`chimera_panel()` becomes:

```r
chimera_panel <- function() {
  tagList(
    tags$div(style = "margin:2px 0 6px; font-weight:600;", uiOutput("chimera_result")),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      "Looks for pieces of different origin joined in one contig, using the genes' BLAST hits, ",
      "200 bp BLAST windows against the bundled rDNA database, ribotyper (checks that each gene is one whole gene), ",
      "and read depth and clipping. Severe (different phyla or kingdoms) is flagged at export; minor is only reported. ",
      "Unresolved windows match several classes equally and do not count. ",
      "The Origin track and the dashed junction lines in the sequence viewer above show where."),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      "Click a row to zoom the viewer to it; expand a row for its windows."),
    reactable::reactableOutput("chimera_segments"),
    htmltools::tagAppendAttributes(class = "rp-subsection", rp_details("Possible junctions",
      tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
        "Where the checks think two pieces may be joined, which checks point there, and why."),
      uiOutput("chimera_junctions"))))
}
```

Server: delete `output$chimera_plot` and `output$chimera_db_source`. `output$chimera_result`:

```r
  output$chimera_result <- renderUI({
    edit_tick(); refresh(); u <- hit_unit()
    if (identical(input$mode, "Annotate") && !identical(sel_check_state(), "done"))
      return(chimera_result_text(NA, NA, state = sel_check_state()))
    ch <- app_query(sprintf("SELECT label, reason FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
      if (identical(input$mode, "Annotate")) "review_chimera_v" else "chimera"), list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$label[1]))
    chimera_result_text(ch$label[1], ch$reason[1])
  })
```

`output$chimera_junctions`: when `chim_unit()` is not done it `req()`s out, which is fine now that the pill says so.

Delete `plot_chimera()`, and `fixed_gutter()` if `grep` finds no other caller. Fix the "(, NA%)" fragment in `chimera_reason()` (R/chimera.R:402-403): `sprintf("%d-%d %s (%s, %.1f%%)", ...)` becomes a per-row paste that drops the parenthetical when `species` is NA/empty and `pident` is NA:

```r
  segtxt <- if (!is.null(segments) && nrow(segments)) paste(mapply(function(a, b, t, sp, p) {
    det <- c(if (!is.na(sp) && nzchar(sp)) sp, if (!is.na(p)) sprintf("%.1f%%", p))
    paste0(a, "-", b, " ", t, if (length(det)) paste0(" (", paste(det, collapse = ", "), ")"))
  }, segments$start, segments$stop, segments$taxon, segments$species %||% NA_character_, segments$pident), collapse = "; ")
    else "no window segments"
```

Grep tests for `"(, NA%)"` or exact reason strings (`grep -n "nohit (" tests/testthat/*.R`) and update expectations.

- [ ] **Step 4: Run the whole suite**

Run: `Rscript -e 'devtools::test()'` (long; or at least test-chimera, test-seqview, test-batch-trim, test-triage-modal, test-review*). Expected: all pass. `R CMD check`-level: run `Rscript -e 'devtools::document(roclets = "rd")'` only if a roxygen block changed its signature (none here are exported).

- [ ] **Step 5: Commit**

```bash
git add R/app_modals.R R/app_server.R R/chimera.R tests/testthat/test-chimera.R tests/testthat/test-batch-trim.R
git commit -m "Chimera panel: verdict inline, not-run state, static plot removed"
```

---

### Task 6: Open the viewer on a flagged contig in the Assemble modal

**Files:**
- Modify: `R/app_server.R:2258` (the examine observer that collapses `#asm_unit_panels details`)
- Test: `tests/testthat/test-chimera.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("examining a flagged contig opens its viewer section", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1',1,1,'severe')")
  DBI::dbDisconnect(con)
  js <- character(0)
  local_mocked_bindings(runjs = function(code) js <<- c(js, code), .package = "shinyjs")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()
  })
  expect_true(any(grepl("seqv_asm-canvas", js, fixed = TRUE) & grepl("open = true", js, fixed = TRUE)))
})
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL, no runjs call mentions `seqv_asm-canvas`.

- [ ] **Step 3: Implement**

After the existing collapse `runjs` at R/app_server.R:2258 add:

```r
    lab <- app_query("SELECT label FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?",
      list(id, h$path, h$scaffold))$label
    if (length(lab) && !lab[1] %in% c(NA, "clean")) shinyjs::runjs(
      "var d = document.getElementById('seqv_asm-canvas').closest('details'); if (d) d.open = true;")
```

Use whatever names the observer already has for the id and the hit (read the surrounding lines; the pattern is `h <- sel_hit_unit(); id <- sel_id()`).

- [ ] **Step 4: Run tests, commit**

```bash
git add R/app_server.R tests/testthat/test-chimera.R
git commit -m "Assemble modal opens the viewer for a flagged contig"
```

---

### Task 7: Drive the app and look

**Files:** none changed unless the drive finds a bug.

- [ ] **Step 1: Launch on the test project**

Follow the recipe in the `drive-app-headless` memory: Chrome at the flatpak path, `RiboPilot(proj, port = P, launch.browser = FALSE, headless = FALSE)`, a keeper chromote session, short `TMPDIR` under `/dev/shm`. Project: `/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/` (CHIM01 is the severe chimera).

- [ ] **Step 2: Check, screenshot to the scratchpad**

1. Assemble modal, CHIM01: viewer section is open, origin band shows grey then yellow with a red dashed line near 2625, hover on the band names Eurotiomycetes / Penicillium, the "% id" track has a gap over no-hit windows.
2. Zoom in until the 1-window segments get their taxon text; zoom out until long labels drop to the short form.
3. Click the second segments row: viewer zooms to it. Expand the row: windows listed; click a window: viewer zooms. Click a junction row: zooms around 2625.
4. Annotation Review, CHIM01: same, pill with the reason visible under it. Trim the contig: pill says "not run", origin band gone, no console errors. Press Update later is out of scope; Restore and confirm the band returns.
5. Browser console: no errors (`chromote` `Runtime.consoleAPICalled` or read via `b$Log`).

- [ ] **Step 3: Fix anything found, commit as separate small commits.**

---

### Task 8: Docs and NEWS

**Files:**
- Modify: `vignettes/Curation-Flags.Rmd:198-222`, `vignettes/Test-Project-Assemble.Rmd:93-104`, `vignettes/Test-Project-Annotate.Rmd:47`
- Modify: `NEWS.md` (development version section)

- [ ] **Step 1: Text**

Curation-Flags.Rmd: replace "open the **Chimera check** panel, and hover over the contig's badge in the **Chimera** column to read the evidence behind the call" with "open **Sequence and read coverage** and switch on **Origin**: each stretch is coloured by the class its BLAST windows match (grey is your lineage), and a dashed line marks a possible junction. The **Chimera check** panel below gives the call, the reason, and the segments; click a row to zoom the viewer to it." Keep the two screenshot references but note in the commit that they need re-shooting (`assemble-chimera-chim01.png`, `chimera-mild.png`); re-shoot during Task 7 if the drive is already up, saving to the same paths.

Test-Project-Assemble.Rmd:101-104: reword "the **Chimera check** panel shows more detail. You can see where ..." to point at the Origin track in the viewer and the segments table.

- [ ] **Step 2: NEWS**

Under "# RiboPilot (development version)" add a "## Chimera check" section:

```
## Chimera check

* The chimera evidence now lives in the sequence viewer: an Origin track colours
  each stretch by the class its BLAST windows match, a Window identity track
  shows the hit identity, and dashed lines mark possible junctions. The static
  chimera plot is gone.
* The Chimera check panel shows the reason under the call instead of on hover,
  says "not run" after a trim, and lists windows under their segment. Clicking a
  segment, window, or junction row zooms the viewer to it.
* The lineage palette is down to three colours that stay apart under red-green
  colour blindness.
```

- [ ] **Step 3: Commit**

```bash
git add vignettes/Curation-Flags.Rmd vignettes/Test-Project-Assemble.Rmd vignettes/Test-Project-Annotate.Rmd NEWS.md
git commit -m "Docs: chimera evidence in the sequence viewer"
```

---

## Self-review notes

- Spec coverage: origin track, identity track, junction lines, hover, tiered labels, 3-colour palette, de-overlap, tables kept with windows nested, click-to-zoom on all three tables, verdict inline, not-run state, static plot removed, Assemble auto-open, docs. Trim modal and action buttons deferred by design.
- Names used across tasks: `chimera_origin`, `seqview_chimera`, `seqview_goto`, `rt_zoom_row`, `fitRange`, `rpseq_goto`, input `chim_zoom`, checkboxes `show_origin` / `show_ident`. Consistent.
- Review Focus items 1-3 have tests (Tasks 1, 2, 5); 4 and 5 are visual (Task 7).
