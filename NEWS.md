# RiboPilot 0.1.2

## Browser support

* Opening the app in Safari now shows a warning on every launch that
  Safari has known display and interaction issues and suggests Chrome or
  Firefox.

## Help

* Longer help text across the app moved behind "?" icons and was shortened:
  options windows, export window, NMNH voucher report, details windows, triage,
  trim, reference settings, and the cluster job script.
* Every "?" icon now looks and behaves the same: click to open, one open at a
  time, closes on click away or Escape. Tool documentation keeps its wide
  scrolling panel.
* The Chimera check panel says "features" where it said "genes".

## Annotate

* The Partial column is now Completeness and is worked out automatically: a
  contig is "complete" when all five features are present and no gene end is
  partial, otherwise "partial". It can be sorted and filtered. The Mark
  Partial buttons are gone; set a gene end with the 5' and 3' buttons in the
  contig window instead. The `{partial}` export token follows the same rule.
* The two sliding-window BLAST tracks in the sequence viewer are now labeled
  "window" (Window BLAST identity, Window BLAST taxon), grouped by a bracket, to tell them apart
  from Feature BLAST.
* The sequence viewer has a Feature BLAST track above the window BLAST tracks:
  one block per feature, green for OK, orange for a weak or low-confidence
  match, red for off-target or non-rDNA, grey for no hit or not searched. In
  the contig window, clicking a block opens the Feature BLAST table.
* Hovering the sequence viewer now shows only the value of the track under
  the pointer, and its help moved behind a "?" icon beside the viewer title.
* Contigs that are both locked and excluded have their own red lock icon, sort
  position, and Lock filter option.

## Assemble

* The annotate and curate option links grey out when the sample's annotation
  is locked, like the other option links do for a locked assembly.

## Export

* Assign Group now offers the groups that already exist, and warns when a
  selected contig is already in another group.
* Clear Group asks for confirmation first.
* The Export table has an Exported filter (Yes / No).

## Fixes

* Tables in the Run reports, Assemble, and other detail windows showed only
  a sliver of their rows in Safari: the horizontal scrollbar covered the
  last row and a vertical scrollbar appeared. Tables that size to their rows
  now include the bar.
* A new project's `.config` file is readable by your group (mode 0640, was
  0600), so others in the group can work with a shared project.

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
  with NCBI's wording. See "The GenBank BLAST
  check" in the GenBank Submission guide.
* Most fish records trip both flags because the ITS database holds few fishes;
  the flags predict GenBank's reply, not a bad assembly.
* The check needs `ITS_eukaryote_sequences` beside the panel. rDNA panel
  2026.08.18.2 ships it and is included in the 0.1.2 container image. With an
  older panel every contig stores
  `not_run(no ITS_eukaryote_sequences database)`.
* Project database version 7: `genbank_check` columns on `annotate` and
  `review_units`. Contigs already locked keep `not_run(not checked)` until
  they are unlocked and rechecked.

## Reference database

* The rDNA BLAST panel is 2026.08.18.2. The rDNA sequences are the same as in
  2026.08.18.1; this release adds NCBI's `ITS_eukaryote_sequences` database
  for the GenBank BLAST check and refreshes the NCBI 16S and taxonomy
  databases shipped beside the panel.

## Reference folder

* `clear_ref_cache()` deletes the cached reference data in the shared
  reference folder (clade databases, record counts, and stored reference
  sequences). It affects every project that uses the folder, so it lists what
  it found, asks you to type `clear`, and refuses to run from a script. The
  BLAST panel and container images are left alone.

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
