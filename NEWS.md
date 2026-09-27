# RiboPilot (development version)

RiboPilot is in development. Its version number stays at 0.0.0.9018 for now, so
read this page, not the number, to see what changed.

## Before you upgrade

* Recreate your projects. The project database changed, and RiboPilot does not
  yet convert an existing project. See "Upgrading RiboPilot" in the Installation
  guide.
* Use a container image built from this code. The published 0.0.0.9018 image
  holds older code, and a run with it stops before doing any work.

## GenBank records

* Before it writes anything, the Export window lists every unit the BLAST screen
  or the chimera check flagged that nobody marked Reviewed, and leaves them out
  unless you choose Export anyway. Lock asks about the same units, on both tables
  and in the unit window, and the Export complete window names each flagged
  record.
* A gene the curation step finds truncated is written as partial.
* `[organism=]` uses the name confirmed in the Taxonomy step. A record whose Taxon
  differs gets `taxon_differs(<name>)` in `submission_risk`.
* A gap of 5 bp or less between a gene and its spacer is closed; a larger one
  gets a `gap` warning. Assembly gaps use "unspecified" linkage evidence. An
  empty source modifier in a header template gives a warning.
* Header templates are filled by column name and never run as code, and export
  group names must work as folder names.
* `export_project()` exports every exportable unit as one "export" bundle when
  there are no export groups, by the same rule as the app.

## Quality checks

* A Checks column on the Assemble and Annotation Review tables sums up each
  sample or unit: red for off target, contaminant, severe chimera, or ITSx
  chimeric; amber for things worth a look; green "ok".
* Chimera calls are clean, minor, or severe. Only severe calls are a submission
  risk; minor calls are listed in the report.
* Changing a sample's taxon updates its lineage and screens it again. A name
  that matches more than one organism is never picked automatically.
* Isomer and paralog classes are computed after annotation finishes.

## Re-running samples

* Assembling a sample again clears the old genes, screen, chimera call, BLAST
  results, lock, and marks of every scaffold whose sequence changed. Your own
  keep and drop calls follow their sequence.
* Changing annotation options re-runs only annotation, and keeps trims.
* Changing a BLAST setting searches again.
* A rescue with the baited assembler is annotated in the same run.
* An rRNA scan of a supplied assembly that does not finish fails the sample
  instead of keeping part of the result.

## Runs

* The end of a run says "Run finished with no failures", or names each failed
  sample and step with the exit code and the likely cause. **Copy report**
  copies the report, which is also saved in `nextflow_logs/`.
* A failed sample's window says why it failed and what to try next.
* Closing the browser tab, or pressing Stop / Interrupt, stops a run cleanly.
  Only one run goes at a time per project.
* Each task is capped at this computer's CPUs and memory.
* The rescue for a failed sample is the baited SPAdes choice under Assembler in
  Assembly Options; `run_fallback()` is removed.
* Cluster job scripts for SLURM, SGE, PBS, and LSF take the queue, account, and
  walltime from your profile. The awsbatch profile is removed.

## Projects

* `new_project()` and `new_test_project()` work outside RStudio and use Docker
  by default. Conda is unsupported and warns.
* The mapping file, IDs, and paths are checked before anything is written.
* New `add_samples()` and `update_sample_metadata()`.
* A locked unit is read only. Problematic is now Exclude from export. Lock &
  Close only locks; Mark Reviewed is separate. Only one assembly path per sample
  can be locked.

## Working at scale, reproducibility, and clusters

* The Lock and State filters now drop rows from the underlying table instead of
  hiding them with CSS, so the table's paging and counts always match what the
  filters allow; refreshing the table keeps its page and selection. The header
  tick box, bulk actions, and Triage all act only on rows the search and column
  filters still show, and Triage skips assembly-locked samples with a toast
  rather than silently including them.
* Every tool chain package in the container image is pinned to the version
  actually installed, and the pins match `environment.yml`.
* A run stops before doing any work when the RiboPilot that launched it and the
  RiboPilot installed in the container differ in version or in the pipeline
  scripts themselves, instead of only warning.
* The database that answered a unit's BLAST search (the bundled rDNA database,
  with its version, or NCBI core_nt) is recorded for every unit, and the BLAST
  evidence panel names it.
* On the NMNH Hydra cluster, the rRNA search of a supplied assembly, and the
  other per-sample memory steps, ask Hydra for memory the way it expects.

## Documentation

* A new Get started guide that follows the app, and new pages on installation,
  your own project, GenBank submission, scripting, clusters, caches and offline
  use, troubleshooting, and citing RiboPilot.
* This changelog, `citation("RiboPilot")`, and the tool versions in the
  container image, for methods sections.
* The docs and the app use the same names: "step", "BLAST screen", and "bundled
  rDNA database". RiboPilot does not ship a reference database; build one for
  every group in the Reference step.
* Help pages show their formatting and links.
