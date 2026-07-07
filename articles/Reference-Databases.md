# Using RiboPilot with your group of organisms

To assemble a sample, RiboPilot compares its reads to known ribosomal
sequences from related organisms. It comes ready to use for **ray-finned
fishes**. For any other group, you build that reference data once, and
it is reused for every future project on that group.

You do not need to find or format any sequences yourself. RiboPilot
downloads them from GenBank for you.

## Build reference data for a group

Use
[`build_reference_db()`](https://dmacguigan.github.io/RiboPilot/reference/build_reference_db.md).
It takes a standard NCBI search (the same kind of query you would type
on the [NCBI website](https://www.ncbi.nlm.nih.gov/)) and a short label
for the group. It then downloads the relevant 18S, 5.8S, 28S, and ITS
sequences and prepares the reference data.

``` r

# 1. First, just check how many sequences are available. Nothing is downloaded.
RiboPilot::build_reference_db(
  '"Gastropoda"[Organism]', "gastropoda",
  count_only = TRUE, conda_env = "ribopilot"
)

# 2. Then build the reference data and save it somewhere you can reuse.
db <- RiboPilot::build_reference_db(
  '"Gastropoda"[Organism]', "gastropoda",
  out_dir = "~/ribopilot_refs", conda_env = "ribopilot"
)
```

Replace `Gastropoda` with your group. The `"[Organism]"` part tells NCBI
to search by organism name, so `'"Gastropoda"[Organism]'` means “all
gastropod sequences”. The second value, `"gastropoda"`, is just a short
label used in the file names.

Running this needs the tools you installed with conda, so start R from a
terminal where you have run `conda activate ribopilot`, or keep the
`conda_env = “ribopilot”` setting shown above.

## Use it in a project

When you set up a project, point it at the reference data you built by
adding two lines to
[`new_project()`](https://dmacguigan.github.io/RiboPilot/reference/new_project.md):

``` r

RiboPilot::new_project(
  path       = "~/my_project",
  mapping_fn = "~/my_samples.csv",
  data_path  = "~/my_reads",
  executor   = "conda",
  custom_seeds_db  = db$seeds,
  custom_labels_db = db$labels
)
```

## Good to know

- **Build once, reuse forever.** If you run
  [`build_reference_db()`](https://dmacguigan.github.io/RiboPilot/reference/build_reference_db.md)
  again for the same group and location, it reuses what you already
  built instead of re-downloading. To force a fresh rebuild, add
  `overwrite = TRUE`.
- **Where it is saved.** If you leave out `out_dir`, RiboPilot saves the
  reference data in a default folder on your computer and finds it
  automatically next time. Setting `out_dir` to a shared location lets
  your whole lab reuse the same data.
- **Groups with few sequences.** Some groups have very little published
  ribosomal data. If your search returns only a handful of records
  (check with `count_only = TRUE` first), broaden it to a larger group.
  For example, search a whole family or order rather than a single
  genus, so RiboPilot has enough reference sequence to work with.
- **The 5.8S gene.** GenBank rarely stores 5.8S on its own; it is
  usually bundled inside ITS records. That is fine: RiboPilot still
  recovers the 5.8S region from the full repeat, and its quality checks
  anchor on the 18S and 28S genes.
