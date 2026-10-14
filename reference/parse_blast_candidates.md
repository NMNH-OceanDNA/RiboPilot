# Parse blastn hits into ranked reference candidates + a contamination screen

Reads a blastn `-outfmt "6 saccver stitle pident qcovs evalue staxids"`
file for one sample and splits it two ways:

## Usage

``` r
parse_blast_candidates(
  con,
  id,
  path = 1L,
  scaffold = 1L,
  hits_file,
  taxdump_dir = NULL,
  review = FALSE
)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- hits_file:

  blastn outfmt-6 output (or a "NO_SIGNIFICANT_HITS" sentinel).

- taxdump_dir:

  Optional pinned taxdump for resolving the screen hit's lineage.
  Without it the lineage is left NA, which curation flags rather than
  treats as clean.

- review:

  TRUE stores the hits of the contig's review copy in
  `review_blast_hits` and calls derive_blast_products() in review mode.
  The sample's lineage taxids are never filled in.

## Value

(invisibly) the ranked candidate data.frame.

## Details

- **Reference panel** - hits whose title is nuclear rDNA. Deduped by
  accession (best pident\*qcovs), ranked into `blast_ref_candidates`,
  rank-1 written to the unit's `assemblies.blast_*` columns. Same
  records the restricted search returned before.

- **Contamination screen** - the top hit of the whole search (the
  bundled rDNA database with its decoys, or NCBI core_nt), unfiltered,
  into `blast_screen`. This is the only signal that can distinguish "the
  host" from "a symbiont" or "not nuclear rDNA at all", so it is
  deliberately taken before the nuclear-rDNA filter.
