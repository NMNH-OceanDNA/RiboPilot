# Write fetched records and their annotations into the shared store

Best effort. A reference folder that cannot be written leaves the run
unaffected: the project keeps its own copy either way.

## Usage

``` r
push_to_ref_store(ref_dir, fasta_file, ann, tool)
```

## Arguments

- ref_dir:

  Reference folder.

- fasta_file:

  The multi-record fasta that was fetched.

- ann:

  annotate_records() output, or NULL.

- tool:

  Pipeline-level provenance label for this write, e.g.
  `"ribopilot-annotate"`. Not the name of one external binary: anchors
  come from whichever rRNA caller is configured (pybarrnap or barrnap)
  and spacers from ITSx, so a single per-accession string cannot name
  both.

## Value

(invisibly) the number of accessions written.
