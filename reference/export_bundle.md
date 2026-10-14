# Build a GenBank submission bundle for a set of samples (one export group)

Build a GenBank submission bundle for a set of samples (one export
group)

## Usage

``` r
export_bundle(
  path = ".",
  ids = NULL,
  units = NULL,
  group = NULL,
  out_dir = file.path(path, "out"),
  label = NULL,
  markers = TRUE,
  report = TRUE,
  template = NULL,
  run_table2asn = FALSE,
  zip = TRUE,
  export_opts = "default",
  fasta_header = NULL,
  fasta_header_marker = NULL
)
```

## Arguments

- path:

  Project directory.

- ids:

  Sample IDs to export (default: every exportable unit). If `group` is
  given, restricted to that group. Sample-grain: every exportable unit
  of a listed sample is included. Use `units` to pick individual
  contigs.

- units:

  Optional data.frame(ID, path, scaffold) restricting the export to
  those exact assembly units. Applied on top of `ids` / `group`.

- group:

  Export-group name to filter by and label the bundle.

- out_dir:

  Project "out" directory (default `<path>/out`); files land in
  `<out_dir>/export/<label>`, mirroring MitoPilot.

- label:

  Bundle label (default: `group` or "export").

- markers, report:

  Also write per-marker fastas / an HTML report.

- template:

  Optional path to a GenBank submission template (`.sbt`); copied into
  the bundle as `<label>.sbt` so table2asn picks it up.

- run_table2asn:

  Run table2asn if on PATH.

- zip:

  Also write `<out_dir>/export/<label>.zip`.

- export_opts:

  Name of the saved FASTA-header configuration to use.

- fasta_header, fasta_header_marker:

  Optional header templates (`{column}` tokens) overriding the saved
  configuration.

## Value

(invisibly) the group output directory, with the summary data.frame
attached as attribute `summary` and the contigs left out because their
checks are not finished as attribute `not_done` (data.frame(ID, path,
scaffold)).

## Details

A unit is exported when its sample's assembly and the unit's annotation
are both locked, its contig is kept, it is not excluded from export, and
its checks are done (press Update on Annotation Review after an edit).
Every value written comes from the reviewed copy of the contig when it
has one. The `itsx_check` column of the summary CSV reports the stock
ITSx check stored for each contig (`not_run(not checked)` when none is
stored); export does not run ITSx.
