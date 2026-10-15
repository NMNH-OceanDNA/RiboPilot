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
