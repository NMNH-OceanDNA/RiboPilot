# Download a pinned NCBI taxdump tarball into \`dir\` (once per project)

Stores only the ~72 MB compressed tarball plus a provenance file. The
dmp files are extracted to a temp location at resolution time (see
\[extract_taxdump()\]), not persisted in the project.

## Usage

``` r
ensure_taxdump(dir, force = FALSE)
```

## Arguments

- dir:

  Destination directory (e.g. \`\<project\>/.taxonomy\`).

- force:

  Re-download even if a tarball already exists.

## Value

\`dir\`.
