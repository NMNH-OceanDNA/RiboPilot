# Build a GenBank submission bundle for a set of samples (one export group)

Build a GenBank submission bundle for a set of samples (one export
group)

## Usage

``` r
export_bundle(
  path = ".",
  ids = NULL,
  group = NULL,
  out_dir = file.path(path, "export"),
  label = NULL,
  markers = TRUE,
  report = TRUE,
  template = NULL,
  run_table2asn = FALSE,
  zip = TRUE
)
```

## Arguments

- path:

  Project directory.

- ids:

  Sample IDs to export (default: all with a stitched sequence). If
  \`group\` is given, restricted to that group.

- group:

  Export-group name to filter by and label the bundle.

- out_dir:

  Output directory (default \`\<path\>/export\`).

- label:

  Bundle label (default: \`group\` or "export").

- markers, report:

  Also write per-marker fastas / an HTML report.

- template:

  Optional path to a GenBank submission template (\`.sbt\`); copied into
  the bundle as \`\<label\>.sbt\` so table2asn picks it up.

- run_table2asn:

  Run table2asn if on PATH.

- zip:

  Also write \`\<out_dir\>/\<label\>.zip\`.

## Value

(invisibly) the group output directory.
