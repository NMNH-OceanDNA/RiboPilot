# RiboPilot (development version)

## Fixes

* Tables in the Run reports, Assemble, and other detail windows showed only
  a sliver of their rows in Safari: the horizontal scrollbar covered the
  last row and a vertical scrollbar appeared. Tables that size to their rows
  now include the bar.
## GenBank BLAST check

* A new check predicts the two BLAST flags GenBank's rRNA/ITS curation raises
  before it assigns accessions: top-hit query coverage under 75% ("Query
  coverage is less than 75%") and a top hit aligned in more than one piece
  ("HSP_problem"). RiboPilot runs NCBI's own searches (default megablast
  against `ITS_eukaryote_sequences` and the SSU/LSU eukaryote rRNA sets) on
  each contig in the annotation step and again on the recheck, and stores the
  result as `genbank_check`. It reproduces NCBI's reply exactly on the ten
  records of SUB16547976.
* The result shows in the Checks column (amber, never a gate), in a **GenBank
  BLAST** column on the scaffold and Annotation Review tables, and at export:
  a `genbank_check` column in `<group>_summary.csv` and the HTML report, a
  count of flagged records in the Export done window, and a README section
  with NCBI's wording and the reply to send gb-admin. See "The GenBank BLAST
  check" in the GenBank Submission guide.
* Most fish records trip both flags because the ITS database holds few fishes;
  the flags predict GenBank's reply, not a bad assembly.
* The check needs `ITS_eukaryote_sequences` beside the panel. It ships with the
  next panel release; until then every contig stores
  `not_run(no ITS_eukaryote_sequences database)`.
* Project database version 7: `genbank_check` columns on `annotate` and
  `review_units`. Contigs already locked keep `not_run(not checked)` until
  they are unlocked and rechecked.

# RiboPilot 0.1.1

## Assemble

* The Assemble details window aligns a sample's assembly paths and shows the
  alignment in the sequence viewer, with differing columns marked.
* New projects default to a maximum of 20 assembly paths and 20 contigs per
  path (was 10). Existing projects keep their settings.

# RiboPilot 0.1.0

## Before you upgrade

* Make your existing projects again with `new_project()`. Projects made before
  this version do not record their database version and cannot be updated.
* From this version on, a project records its database version. After a later
  upgrade, run `backwards_compatibility()` on each project instead of making it
  again. See
  ["Upgrading RiboPilot"](https://nmnh-oceandna.github.io/RiboPilot/articles/Installation.html#upgrading-ribopilot)
  in the Installation guide.

## Reference database

* The rDNA BLAST panel is 2026.08.18.1. It leaves out six Canis aureus
  PREDICTED rRNA records (XR_013375787, 788, 790, 791, 794, and 795) whose
  sequences are flatfish contamination; they produced a false "PREDICTED: Canis"
  chimera flag on a flounder sample. Known-bad GenBank records are listed in
  `dev/refdb/exclusions.csv` and dropped at build time.
