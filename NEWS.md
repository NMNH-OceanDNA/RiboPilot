# RiboPilot (development version)

## Before you upgrade

* Make your existing projects again with `new_project()`. Projects made before
  this version do not record their database version and cannot be updated.
* From this version on, a project records its database version. After a later
  upgrade, run `backwards_compatibility()` on each project instead of making it
  again. See
  ["Upgrading RiboPilot"](https://dmacguigan.github.io/RiboPilot/articles/Installation.html#upgrading-ribopilot)
  in the Installation guide.

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
  results, lock, and marks of every contig whose sequence changed. Your own
  keep and drop calls follow their sequence.
* Changing annotation options re-runs only annotation, and keeps trims.
* Annotation options take extra command-line flags for the rRNA caller, one
  field for pybarrnap and one for barrnap, since the two take different flags.
  The flags go to both per-contig annotation and the rRNA scan of a supplied
  assembly. Existing projects need `backwards_compatibility()` once.
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
* New `backwards_compatibility()` updates a project's database and points the
  project at the container image for the installed RiboPilot, after saving a
  copy of both. The next run downloads that image by itself. The app does not
  open a project that needs this update, and says how to run it.
* A locked unit is read only. Problematic is now Exclude from export. Lock &
  Close only locks; Mark Reviewed is separate. Only one assembly path per sample
  can be locked.

## Working at scale, reproducibility, and clusters

* The Assemble details window has a Notes box for your own notes on a sample.
  They show in the Assemble table's User notes column, and a unit's notes show
  in the Annotation Review table's User notes column.
* The Partial column reads "no" for every unit not marked Partial, so a filter
  on "no" lists only complete units.
* The Checks column can be dragged wider.
* The chimera check plot, segments, and windows redraw after a trim, a Restore,
  or a feature boundary edit, instead of only when another unit is picked.
* Trim and Restore buttons on the Assemble and Annotation Review toolbars trim
  many units at once, to their gene ends or gene ends plus a flank, and undo it.
  Only units with all five features are trimmed; the window lists the rest and
  why before you confirm.
* The Lock and State filters now drop rows from the underlying table instead of
  hiding them with CSS, so the table's paging and counts always match what the
  filters allow; refreshing the table keeps its page and selection. The header
  tick box, bulk actions, and Triage all act only on rows the search and column
  filters still show, and Triage skips assembly-locked samples with a toast
  rather than silently including them.
* Clicking a row's lock or state icon on the Taxonomy, Assemble, and Annotation
  Review tables locks, unlocks, or sets the state of that row alone, with the
  same questions as the toolbar buttons. The Unlock and "Edit beyond the
  selection" questions have a "Don't ask again this session" box.
* Every tool chain package in the container image is pinned to the version
  actually installed, and the pins match `environment.yml`.
* Opening a project whose container image is for another RiboPilot version shows
  a message that says so and how to switch; runs carry on with that image.
* The database that answered a unit's BLAST search (the bundled rDNA database,
  with its version, or NCBI core_nt) is recorded for every unit, and the BLAST
  evidence panel names it.
* On the NMNH Hydra cluster, the rRNA search of a supplied assembly, and the
  other per-sample memory steps, ask Hydra for memory the way it expects.

## Sample metadata

* Link samples to GEOME, GBIF, and NCBI BioSample records through mapping-file
  columns `GEOME_BCID`, `GBIF_ID`, and `BioSample`, or other names via
  `mapping_geome`, `mapping_gbif`, and `mapping_biosample`.
* `new_project()`, `new_test_project()`, `add_samples()`, and
  `update_sample_metadata()` fetch the records, unless `fetch_geome`,
  `fetch_gbif`, or `fetch_biosample` is FALSE. `link_sources = TRUE` also
  follows links between databases.
* New `fetch_geome()`, `fetch_gbif()`, `fetch_biosample()`, `remove_metadata()`,
  and `set_metadata_columns()`, plus `geome_normalize_bcid()`,
  `gbif_normalize_id()`, and `ncbi_normalize_id()`.
* The NCBI API key is read from the project `.config`, then `NCBI_API_KEY`, then
  `ENTREZ_KEY`.
* A Metadata column in the Taxonomy, Assemble, Annotation Review, and Export
  tables opens a viewer with per-source tabs and a Compare tab. A Metadata
  button picks which fields show as table columns.
* Set Export Metadata picks fields for header templates and
  `<group>_summary.csv`. Export asks before writing when sources disagree on
  something a template uses. A mapping-file column with the same name as a
  ticked field wins.

## Documentation

* A new Get started guide that follows the app: a short page on setup, then a
  test project walkthrough with one page per step (Taxonomy, Reference,
  Assemble, Annotation Review, Export). Boxes tagged HYDRA mark where the NMNH
  Hydra cluster differs.
* New pages on installation, your own project, GenBank submission, scripting,
  clusters, caches and offline use, troubleshooting, and citing RiboPilot.
* This changelog, `citation("RiboPilot")`, and the tool versions in the
  container image, for methods sections.
* The docs and the app say "contig" where they said "scaffold": every assembled
  sequence is gapless. Names you type or read in files are unchanged: the
  `scaffold` column and `{scaffold}` header token, `max_scaffolds`, and record
  names such as `SRR33849683_p1_s2`.
* The bundled example is called the test project, as in MitoPilot.
* The docs and the app use the same names: "step", "BLAST screen", and "bundled
  rDNA database". RiboPilot does not ship a reference database; build one for
  every group in the Reference step.
* Help pages show their formatting and links.
