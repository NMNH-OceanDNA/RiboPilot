# RiboPilot export revamp

Date: 2026-07-31
Branch: `export-revamp` (from `wf2-decommission`)

## Goal

Give the curator control over exported FASTA headers, save those choices as named
configurations, check them before anything is written, and reorganise the export
folder. Extend the HTML report with per-feature alignments. Patterns are ported
from MitoPilot where they already work.

## 1. Named header configurations

New table, created by `init_db()` and back-filled by `migrate_db()`:

```sql
CREATE TABLE export_opts (
  export_opts TEXT NOT NULL,     -- configuration name
  fasta_header TEXT,             -- whole-record template
  fasta_header_marker TEXT,      -- per-marker template
  PRIMARY KEY (export_opts)
);
```

Seeded with an `export_opts = 'default'` row holding the package defaults. The
selected configuration is session state in the app (default `'default'`), as in
MitoPilot; editing a header box upserts into the selected name.

Accessors mirror `triage_opts.R`, defensive against a missing table or row:

- `list_export_templates(con)` -> `"default"` first, then the rest sorted
- `get_export_opts(con, name)` -> both templates, package default per empty field
- `set_export_opts(con, fasta_header, fasta_header_marker, name)` -> upsert

## 2. Header templates

`{token}` syntax resolved with `stringr::str_glue_data()` against one row per
exported unit, built by `export_header_data(con, units)`:

- every `samples` column, so every extra column of the user's mapping file is a
  valid token with no extra plumbing
- unit facts: `path`, `scaffold`, `length`, `topology`, `n_features`, `missing`,
  `partial`, `reviewed`, `export_group`
- computed: `seqid`, `completeness`, `elements`, `rrna_definition`
- marker rows additionally carry `element` and `product`

`seqid` keeps its current rule (`ID`, or `ID_p<path>_s<scaffold>` when the sample
exports more than one unit) computed over all of the sample's exportable units, so
IDs stay stable between a partial and a whole export.

`rrna_definition` builds the NCBI definition phrase from the features actually
annotated and their partial marks, e.g. "18S ribosomal RNA gene, partial sequence;
ITS1, 5.8S ribosomal RNA gene, and ITS2, complete sequence; 28S ribosomal RNA gene,
partial sequence".

Defaults:

```
record: {seqid} [organism={Taxon}] [moltype=genomic DNA] {Taxon} {rrna_definition}
marker: {seqid}_{element} [organism={Taxon}] {Taxon} {product}
```

Allowlisted INSDC source modifiers (`SOURCE_MOD_ALLOW`) are still appended to the
record header after the template, unchanged from today.

## 3. Verification

`validate_fasta_header(template, data, require_seqid = TRUE)` returns
`list(ok, level, message)` with `level` in `ok` / `warn` / `error`:

1. empty template -> error
2. unbalanced brace -> error, naming the position and surrounding text
3. dry-run glue against a real row -> error, rewriting glue's "object not found"
   into `column "X" not found in the database`
4. missing `{seqid}` -> error when any sample exports more than one unit (the
   FASTA header and the `>Feature` line would disagree and NCBI would reject the
   submission), warning otherwise

In the app: debounced status line under each box, green/amber/red. On Export, an
invalid template blocks with a centred alert naming the offending box.

At write time `export_bundle()` stops when two units resolve to the same header,
which is what a template without `{seqid}` produces.

## 4. Output layout

```
out/export/<GROUP>/
    <GROUP>.fasta
    <GROUP>.tbl
    <GROUP>_summary.csv
    <GROUP>_report.html
    README.txt
    samples/    <seqid>.fasta, <seqid>.tbl
    markers/    <ELEMENT>/<seqid>_<ELEMENT>.fasta, all_<ELEMENT>.fasta
```

Changes from today: per-unit files move out of the group root into `samples/`,
`.fsa` becomes `.fasta`, `README_SUBMISSION.txt` becomes `README.txt`. The optional
table2asn step is invoked per file rather than over the directory, because its
directory mode only picks up `.fsa`.

## 5. HTML report

`inst/export_report.Rmd`, rendered with `rmarkdown::render()`, self-contained,
floating table of contents. When rmarkdown is unavailable the report step is
skipped with a warning rather than failing the export.

Sections:

1. Summary table (as today)
2. Per-unit feature maps (existing `feature_map_svg()`)
3. Per-feature alignments: for each of 18S, 5.8S, 28S, ITS1, ITS2 present in the
   group, a nucleotide `DECIPHER::AlignSeqs()` alignment across exported units,
   rows ordered by complete-linkage clustering on the DECIPHER distance matrix,
   embedded as an interactive `msaR::msaR` widget with a stat line (units, length
   range, mean identity). Features with fewer than two units say so and are
   skipped.

Sequence gathering lives in `R/export_report.R` (`gather_feature_seqs()`,
`feature_alignment()`) so the Rmd stays thin and the helpers are testable.

## 6. Export dialog

Group picker (unchanged), plus a template-name selectize (typing a new name
creates it), the two header boxes with live status and the list of available
columns, and checkboxes for per-marker files and for the report. Header edits
auto-save into the selected template, debounced, only while both templates
validate.

## 7. Tests

- header data assembly, including a custom mapping-file column
- `rrna_definition` phrasing across feature combinations and partial marks
- every validator branch
- header resolution end to end, plus duplicate-header blocking
- folder layout: `samples/` contents, `.fasta` suffix, group files, `README.txt`
- alignment gathering (per feature, skips singletons)
- report render smoke test, skipped when rmarkdown is unavailable

## Known breakage

Scripts referencing `<seqid>.fsa` or `README_SUBMISSION.txt`, and any expectation
that per-unit files sit in the group root. Old export folders are not migrated;
the next export rebuilds the group directory as it already does.
