# Paths alignment panel in the Assemble details window

Date: 2026-10-15. Branch: `feat/paths-alignment` (off `devel`).

## Purpose

A sample often yields several assembly paths: alternative reconstructions of
the same rDNA repeat. They are usually near-identical, and the curator has to
pick one with no view of where they differ. This panel aligns the homologous
contigs across paths and shows the alignment in the sequence viewer, with the
annotated features on every row and the disagreeing columns made obvious.

## Decisions (user-approved)

- Shown only in the Assemble details window, and only for a sample that has
  two or more paths and at least one homologous set of two or more contigs.
- Homology rule: two contigs are homologous when they share at least one
  anchor element (18S, 5.8S, 28S) in `annotations` and the two element copies
  overlap: an ends-free pairwise alignment covers at least `REF_ALN_MIN_COV` of
  the shorter copy and at least `REF_ALN_MIN_BP` bases (constants from
  `R/blast_ref_align.R`). Sharing the element name alone is not enough: a
  28S-only fragment that holds a different part of 28S than the full contig
  has no shared sequence and gets smeared across the alignment. A homologous
  set is a
  connected component under that relation. Sets of one contig are dropped.
  Contigs with no anchor element are never in a set.
- When a sample has two or more sets, a dropdown picks the set. Otherwise no
  dropdown.
- Dropped contigs stay in the alignment; the row label carries the triage
  verdict.
- Alignment engine: `DECIPHER::AlignSeqs(dna, verbose = FALSE)`, the same call
  `compute_feature_msa()` makes. Pure R; no external tool.
- Viewer: the existing `inst/app/www/seqviewer.js` gains an alignment mode
  (one row per contig) instead of a second widget. msaR is not used here.
- Alignment is computed when the window opens, cached per (sample, set) for
  the session, never stored in the database. No schema change, no pipeline
  change.

## R side

All in a new file `R/paths_align.R`. The first three are pure (no Shiny, no
database) and unit-tested.

### `paths_homology_sets(feats)`

- Input: data frame with `path`, `scaffold`, `element`, one row per annotation
  of the sample (all contigs, all paths).
- Only anchors count: `SCREEN_ANCHOR_ELEMENTS` (`18S`, `5.8S`, `28S`).
- Union-find or repeated merging over contigs: two contigs join when they share
  an anchor whose copies overlap. Sets are the components with two or more contigs.
- Output: list of sets. Each set is `list(label, units)` where `units` is a
  data frame `(path, scaffold)` ordered by path then scaffold, and `label` is
  the anchors present in the set in transcription order joined by `-` plus the
  count, e.g. `"18S-5.8S-28S (4 contigs)"`. Sets are ordered by size
  descending, then by first unit.
- Empty input or no qualifying set: `list()`.

### `paths_align(seqs)`

- Input: named character vector of contig sequences (names are unit keys
  `p<path>.s<scaffold>`).
- Returns the `DNAStringSet` from `DECIPHER::AlignSeqs`, or `NULL` on error or
  with fewer than two sequences. `tryCatch` as in `compute_feature_msa()`.

### `paths_payload(aln, feats, labels)`

- Input: the aligned `DNAStringSet` (names are unit keys); `feats` data frame
  with `unit`, `element`, `start`, `stop`, `strand`, `partial_start`,
  `partial_stop` in native contig coordinates; `labels` named character vector
  of row labels keyed by unit (e.g. `"p1.s1 keep"`).
- Column mapping: for each row, the k-th non-gap character sits at column
  `which(chars != "-")[k]`. A feature with native `start`/`stop` gets
  `pos1`/`pos2` set to the columns of its start and stop bases.
- Output: `list(len = <columns>, rows = list(...))`. Each row is
  `list(unit, label, seq = <gapped, upper case>, features = list(...))` where
  each feature has the shape the viewer already takes: `row`, `gene`, `pos1`,
  `pos2`, `dir`, `partial5`, `partial3`, `color` (from `RRNA_ELEMENT_COLORS`).
  Features with NA start or stop are skipped, as in `seqview_payload()`.
- Rows follow the order of `aln`.

### `paths_send(session, id, set)`

- Shiny-side. Queries `assemblies.sequence` and `assemblies.triage` for the
  set's units and their `annotations`, builds labels (`p1.s1 keep`,
  `p2.s1 drop`, `p2.s2` when triage is NA), calls `paths_align()` then
  `paths_payload()`, and sends `rpseq` with `id = "seqv_paths-canvas"`,
  `unit = <sample ID>`, and the payload. Returns the payload (or `NULL`) so
  the server can cache it and render the summary line.
- The cache is a `reactiveValues`-free plain environment keyed by
  `paste(id, set_index)` held in the server function; cleared when the sample
  changes. Switching sets back and forth never realigns.

## Modal UI (`R/app_modals.R`)

- `assemble_modal(id)` calls `paths_section(id)`, which returns `NULL` unless
  the sample has two or more distinct `path` values in `assemblies` and
  `paths_homology_sets()` on its annotations is non-empty. The section sits
  directly after the Copy FASTA row and before the `<hr>`, open by default.
- Section title: "Compare assembly paths".
- Contents, top to bottom:
  - Help text: "Rows are contigs from each path aligned together. Red ticks on
    the Differences track mark columns where the paths disagree; use the
    arrows to step through them. Letters appear when zoomed in; bases that
    match in every row are greyed."
  - `selectInput("paths_set", NULL, choices = <set labels>)` only when there
    are two or more sets (inline, narrow).
  - `uiOutput("paths_summary")`: "4 contigs, 7,812 columns, 23 differ", or
    "Aligning contigs..." before the payload lands, or "Could not align these
    contigs." when `paths_align()` returns `NULL`.
  - `seqview_ui("seqv_paths", NULL, open = TRUE, pick = FALSE, mode = "align")`
    (see below).
- `seqview_ui()` gains `mode = c("unit", "align")`. In `align` mode the
  controls are: Go to column, Whole alignment, Zoom in, Zoom out, Previous
  difference, Next difference, Copy alignment, and the Nucleotides box. The
  coverage, clipped, BLAST identity, BLAST origin boxes and the Fit gene
  button are not rendered. The help paragraph is omitted (the section supplies
  its own). When `summary` is `NULL`, no `<details>` wrapper is emitted: the
  controls and canvas are returned in a plain `div`.

## Server (`R/app_server.R`, Assemble modal block)

- `paths_sets <- reactive()`: `scaffold_rows()` dependency plus `refresh()`;
  returns `paths_homology_sets()` of the sample's annotations.
- `observe()` on `sel_id()`, `input$paths_set`, and `paths_sets()`: picks the
  set index (1 when there is no dropdown), looks up the cache, otherwise calls
  `paths_send()` and stores the result; then sets a `reactiveVal` that the
  summary output reads. Wrapped in `try(..., silent = TRUE)` like the existing
  `seqview_send` call.
- `output$paths_summary` renders the summary string from that `reactiveVal`.
  The count of differing columns comes from R (`paths_diff_count(rows)`: a
  row's leading and trailing gaps mean "not covered"; a column differs only
  when, among the rows that cover it, the characters (internal gaps included)
  are not all equal; fewer than two covering rows is never a difference) so
  the summary and the viewer agree.

## JS: alignment mode (`inst/app/www/seqviewer.js`)

- `Viewer.prototype.load(p)` sets `this.rows = p.rows || null`. When `rows`
  is set, `this.len = p.len`, `this.seq` is the first row's gapped string (so
  `baseAt` and the ruler keep working on columns), and `this.feats` is the
  concatenation of every row's features with a `rowIdx` added. Lanes are
  packed per row, not globally.
- `this.diff`: an array of differing column numbers computed on load (a
  row's leading and trailing gaps mean "not covered"; a column differs only
  when, among the rows that cover it, the characters, internal gaps included,
  are not all equal; fewer than two covering rows is never a difference).
  Exposed in `rpseq.state(id)` as `diff` (count).
- Layout in align mode, top to bottom: ruler (`RULER_H`), Differences track
  (`DIFF_H = 14`), gap, then for each row: its feature lanes (`LANE_H` each,
  at least one lane so the label has a home) and a nucleotide strip (`NT_H`),
  then a row gap. `height()` is the sum. The gutter shows the row label,
  right-aligned, next to the nucleotide strip.
- Differences track: a pale background across the width, one red
  (`#c0392b`) 1 px tick per differing column in view (`Math.max(1, ppb)`
  wide). Label "Differences" in the gutter.
- Nucleotide strip per row: as `drawNt` today, with two changes. Terminal gaps draw
  nothing; an internal gap is a short grey dash centred in the cell. When letters are visible
  (`ppb >= NT_LETTER`), a column that is not in `diff` is drawn in pale grey
  (`#d9d9d9` fill, `#888` letter); a differing column keeps the base colours
  and gets a light red wash (`rgba(192, 57, 43, 0.12)`) spanning all rows.
  Below the letter threshold every column keeps the base colours so the
  zoomed-out strip still reads as sequence.
- Feature lanes per row: `drawLanes` unchanged in shape, offset by the row's
  y. Partial-end dashing and labels as today.
- Prev/Next difference (`data-rpseq="prev_diff"` / `"next_diff"`): the nearest
  differing column strictly left of the current view's first column (or right
  of its last). Centres it with `goto()` (which raises `ppb` to letter zoom).
  Wraps around at either end. Does nothing when `diff` is empty.
- Clicking in the Differences track zooms to the column under the pointer via
  `goto()`.
- Tooltip in align mode: `Column 1,234 | p1.s1 1,204 A | p2.s1 1,201 G`, each
  row's native position being the count of non-gap characters up to that
  column (a `-` prints as `gap`), then ` | <gene> pos1-pos2 (dir)` when a
  feature is under the pointer.
- Copy alignment: every row in view as gapped FASTA, header `><label>:<a>-<b>`
  with the column range.
- The existing per-unit mode is unchanged: every new branch is behind
  `if (this.rows)`.
- Message handling, `withViewer`, modal-reopen replay, and `rpseq.state` work
  unchanged; `state` adds `rows` (count) and `diff` (count) when in align mode.

## Edge cases

- One path, or no two contigs (counted across all paths) sharing an anchor
  with overlapping copies:
  no section. Three paths of one contig each form one set of three.
- `paths_align()` returns `NULL`: summary says "Could not align these
  contigs."; the canvas stays empty.
- A contig whose sequence is NA is skipped from the set before aligning; if
  fewer than two remain, treated as alignment failure.
- Sequences are used as stored. The annotate step has already oriented contigs,
  so no reverse complement here.
- Size: a handful of 5 to 10 kb contigs aligns in a few seconds and blocks the
  session meanwhile. The summary shows "Aligning contigs..." first. No cap on
  contig count; a `ponytail:` comment names the ceiling.

## Tests (`tests/testthat/test-paths-align.R`)

- `paths_homology_sets()`: one full-length contig plus a split path yields one
  set of three; a stray 28S-only contig with no overlap is dropped; two
  disjoint sets come out ordered by size; no anchors gives `list()`.
- `paths_payload()`: a hand-built two-row alignment with gaps maps a feature
  that spans a gap to the right columns; a feature with NA start is skipped;
  row labels and `len` come through.
- `paths_diff_count()`: counts mismatch and internal gap columns, ignores
  terminal gaps, zero for identical rows.
- `assemble_modal()` on the fixture project renders no `seqv_paths-canvas`
  with one path, and renders it (and no dropdown) once a second path with a
  shared 18S annotation is inserted; a second disjoint set renders the
  dropdown.
- `seqview_ui(mode = "align")` emits the prev/next difference buttons and not
  the coverage boxes.
- `paths_align()` on two short near-identical sequences returns an aligned set
  with equal widths (skip if DECIPHER is not installed).

## Out of scope

- Storing alignments in the database.
- Showing the panel in Annotation Review.
- Reverse-complementing contigs that annotate on the minus strand.
- Per-row coverage or BLAST tracks in the alignment.
