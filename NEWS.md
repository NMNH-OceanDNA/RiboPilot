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
* NMNH user switch in Export data. Adds the `[specimen_voucher=]` and
  `[voucherURI=]` modifiers NMNH requires on GenBank submissions, fills them per
  sample from your mapping file, GBIF, NCBI, or GEOME, checks the voucher codes
  and EZIDs (a tissue EZID is replaced with its parent specimen), and reports
  samples that are missing or wrong. Fix them in the report or with a CSV round
  trip. See `vignette("Sample-Metadata")`.
* Click a sample ID in the NMNH voucher report to open its metadata, then
  return to the report. The metadata window now has a copy button on every
  value and Export CSV to save all of a sample's metadata.
* Tick Not NMNH in the voucher report for samples that are not NMNH
  specimens: their voucher and URI are used as found, without NMNH checks. For
  any sample, pick the mapping-file column or metadata field its
  `specimen_voucher` or `voucherURI` comes from.
* The NMNH user switch is remembered per project. The voucher report's field
  lists are sorted and offer None to leave a modifier out.
* A `[modifier={field}]` with no value for a sample is left out of that
  record's FASTA header instead of being written empty.

## Chimera check in the sequence viewer

* The chimera evidence now lives in the sequence viewer of the Assemble and
  Annotation Review windows: an **Origin** track colours each stretch by the
  class its BLAST windows match (grey is the sample's own lineage), a **Window
  identity** track shows the hit identity, and dashed lines mark possible
  junctions. Hover for the match. The static chimera plot is gone.
* The Chimera check panel shows the reason under the call instead of on hover,
  says "not run" after a trim, and lists windows under their segment. Clicking
  a segment, window, or junction row zooms the viewer to it. On the Assemble
  window, a flagged contig opens the viewer by itself.
* The lineage palette is down to three colours that stay apart under red-green
  colour blindness.

## Quality checks

* A Checks column on the Assemble and Annotation Review tables sums up each
  sample or unit: red for off target, contaminant, severe chimera, or ITSx
  chimeric; amber for things worth a look; green "ok".
* Chimera calls are clean, minor, or severe. Only severe calls are a submission
  risk; minor calls are listed in the report.
* Each core feature (18S, ITS1, 5.8S, ITS2, 28S) is now searched on its own
  against the bundled rDNA database. A new Feature BLAST warnings column in
  Annotation Review lists features with no hit, a non-rDNA hit, a hit outside
  the sample's phylum, or too short to search (under 50 bp), and a Feature
  BLAST panel in the details window shows each feature's matches. Its Search
  column names the algorithm that found the hit; a hit found only by the more
  sensitive dc-megablast retry is shown in amber and listed as a weak hit.
* The chimera check now searches 200 bp windows every 100 bp (was 500 bp every
  250 bp), so it catches shorter foreign stretches. A switch between classes
  still needs about 500 bp on each side, which keeps conserved stretches from
  raising false alarms.
* Changing a sample's taxon updates its lineage and screens it again. A name
  that matches more than one organism is never picked automatically.
* Isomer and paralog classes are computed after annotation finishes.

## Annotation Review edits and Update

* Trims, gene boundary edits, and 5'/3' choices on Annotation Review change a
  working copy of the contig. The original contig and its check results stay as
  the pipeline made them, and the Assemble step shows only the original.
* A trim or a boundary edit marks the contig's checks "not run". The new
  **Update** button on Annotation Review reruns every offline check on those
  contigs, one contig at a time: feature counts, gene completeness, ITSx, the
  BLAST hits, the BLAST screen, the reference alignment, the quality flags, and
  the chimera check. It never re-annotates and never contacts NCBI. Quality
  flags are no longer recalculated in the app right after an edit.
* State on Annotation Review offers only Ready to run and On hold, and queues
  contigs for that Update. A contig whose original still waits for the
  Assemble step's Update says so.
* Lock on Annotation Review refuses the whole click when a selected contig's
  checks are not finished, and lists those contigs with what to do next.
  Exclude & Lock still works. Export leaves such contigs out and counts them.
* Restore deletes the working copy, so the contig shows its original sequence,
  features, and checks again. It undoes trims, boundary edits, and 5'/3'
  choices; marks, notes, and export group stay. Restore on the toolbar works on
  any contig with a working copy, not only trimmed ones.
* Export no longer runs ITSx. The `itsx_check` column of `<group>_summary.csv`
  and the Export done window report the ITSx check stored for each contig, and
  `<group>_itsx_check.csv` and the `itsx_check/` folder are gone. An ITSx run
  that crashed reads `not_run(ITSx did not finish)`, never clean, and the
  Checks column shows it as "ITSx not checked".
* Existing projects need `backwards_compatibility()` once. It turns earlier
  trims and edits into working copies and puts the original contigs back.

## Re-running samples

* Assembling a sample again clears the old genes, screen, chimera call, BLAST
  results, lock, and marks of every contig whose sequence changed. Your own
  keep and drop calls follow their sequence.
* Changing annotation options re-runs only annotation. A contig you trimmed or
  edited on Annotation Review keeps its edits, and its checks are queued for
  Update on Annotation Review; Restore takes the new annotation.
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
  copies the report.
* **Run Reports** in the toolbar lists every run, reference build, and cluster
  job, with its result, and shows each one's report. The app says when a run it
  did not watch, such as a cluster job, has finished. While open, the window
  checks every 15 seconds for new runs and for runs that have finished.
  The run list fits its runs and scrolls past 10.
* Run files go in the project's `.runs/` folder: job scripts, scheduler logs,
  Nextflow logs, and reports. Runs from `run_assembly()` get their own log
  there too.
* A failed sample's window says why it failed and what to try next.
* Closing the browser tab, or pressing Stop or Stop and Close in the run
  window, stops a run cleanly. Only one run goes at a time per project.
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
* The Reference groups table's record count is now "# GenBank Hits", and a new
  "# Database Sequences" column shows how many sequences a built database holds.
  The build drops duplicate and short records, so it is usually much smaller.
* After you change the reference folder, Propose groups points each group at
  that folder's database, or marks it as needing one. Databases you picked
  yourself from elsewhere are kept.

## Working at scale, reproducibility, and clusters

* The Assemble details window has a Notes box for your own notes on a sample.
  They show in the Assemble table's User notes column, and a unit's notes show
  in the Annotation Review table's User notes column.
* The Partial column reads "no" for every unit not marked Partial, so a filter
  on "no" lists only complete units.
* The Checks column can be dragged wider.
* The chimera check panel follows the contig's edits: after a trim or a feature
  boundary edit it says "Not run yet. Press Update." until the recheck has run,
  and after a Restore it shows the original check again. A contig waiting for the
  Assemble step's Update shows "Waiting for the Assemble step's Update." there.
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

* `new_test_project()` now fetches each sample's NCBI metadata from its SRA
  run (CHIM00 and CHIM01 use the Upeneus parvus run they were built from)
  and follows links to GBIF and GEOME, with linking saved as the project
  setting. Pass `fetch_metadata = FALSE` when offline.
* Link samples to GEOME, GBIF, and NCBI BioSample records through mapping-file
  columns `GEOME_BCID`, `GBIF_ID`, and `BioSample`, or other names via
  `mapping_geome`, `mapping_gbif`, and `mapping_ncbi`.
* `new_project()`, `new_test_project()`, `add_samples()`, and
  `update_sample_metadata()` fetch the records, unless `fetch_geome`,
  `fetch_gbif`, or `fetch_ncbi` is FALSE. `link_sources = TRUE` also
  follows links between databases.
* New `fetch_geome()`, `fetch_gbif()`, `fetch_ncbi()`, `remove_metadata()`,
  and `set_metadata_columns()`, plus `geome_normalize_bcid()`,
  `gbif_normalize_id()`, and `ncbi_normalize_id()`.
* The NCBI API key is read from the project `.config`, then `NCBI_API_KEY`, then
  `ENTREZ_KEY`.
* A Metadata column in the Taxonomy, Assemble, Annotation Review, and Export
  tables opens a viewer with per-source tabs and a Compare tab.
* One Metadata fields window, opened by the Metadata button, picks both which
  metadata fields show as table columns (Show) and which can be used at export
  (Export). Tick both together sets both with one click, header boxes tick every
  row in view, GenBank-ready and Ticked only narrow the list, sample counts show
  gaps (for example `12/20`), and a click copies a field's template text. Choose
  fields in Export Data opens the same window filtered to that source. Export
  fields go into header templates and `<group>_summary.csv`. Export asks before
  writing when sources disagree on
  something a template uses. A mapping-file column with the same name as a
  ticked field wins.
* Linking to GBIF by voucher follows Darwin Core. When a GEOME or NCBI record
  has no Smithsonian EZID, RiboPilot searches GBIF for the specimen voucher in
  every field that can hold one (`specimen_voucher`, `genbankSpecimenVoucher`,
  `materialSampleID`, `voucherCatalogNumber`, `otherCatalogNumbers`,
  `bio_material`, or a bare catalog number with the institution code). It reads
  triplets, doublets, `urn:catalog:` prefixes, and `USNM 419933`-style spacing,
  keeps the only match even when its collection code differs, and skips GBIF
  occurrences copied from GenBank or BOLD.
* `geome_specimen_voucher` uses `genbankSpecimenVoucher` when GEOME has it, and
  `gbif_specimen_voucher` drops a catalog number's leading repeat of the
  institution code.

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
