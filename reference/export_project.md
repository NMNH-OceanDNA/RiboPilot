# Export a project's locked units as GenBank bundles

Writes one bundle per export group made in the app's Export step, using
the default header template. When no group was made, writes every unit
with both its assembly and its annotation locked into one bundle named
"export". Only contigs whose checks are done are written; a group with
none is skipped with a message.

## Usage

``` r
export_project(path = ".", out_dir = file.path(path, "out"))
```

## Arguments

- path:

  Project directory.

- out_dir:

  Project "out" directory; bundles land in `<out_dir>/export`.

## Value

(invisibly) the bundle names written.
