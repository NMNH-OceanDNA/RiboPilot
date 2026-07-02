# 01 - Architecture

Mirror MitoPilot's structure closely. Where a component is identical, reuse it;
where the biology differs, swap the internals but keep the interface.

## Package layout (proposed)

```
RiboPilot/
  DESCRIPTION            # same Imports set as MitoPilot (Biostrings, DECIPHER, DBI,
                         #   dbplyr, golem, shiny, processx, ...)
  R/
    init_project.R       # new_project(): same signature as MitoPilot
    init_db.R            # rDNA-specific schema (see below)
    add_samples.R        # reuse ~verbatim
    generate_config.R    # reuse verbatim (executor/config system)
    constants.R          # rDNA element names, expected length ranges per clade
    assemble_*.R         # assembly driver + path scoring for the repeat
    annotate_barrnap.R   # rRNA gene boundaries (18S/5.8S/28S)
    annotate_itsx.R      # ITS1/5.8S/ITS2 boundaries
    annotate_utils.R     # merge annotators into one feature model
    curate_rrna_core.R   # shared curation core
    curate_<clade>_rrna.R# thin per-clade wrappers
    params_<clade>_rrna.R# per-clade expected features + length tolerances
    export.R             # GenBank submission prep
    app_*.R              # Shiny modules (assemble/annotate/export/config)
    ...
  inst/
    nextflow/            # main.nf + modules/ (rDNA workflow)
    config.*             # copy MitoPilot's executor templates unchanged
    app/                 # golem app assets
    ref_dbs/             # seed + label + profile DBs (see 04)
  ref_dbs/               # packaged/downloadable reference sets
  tests/
  docs/
```

## Input map file (unchanged)

Same CSV contract as MitoPilot: required `ID`, `Taxon`, `R1`, `R2`; any number of
extra metadata columns preserved. `new_project(mapping_fn=, data_path=, ...)`
keeps the same arguments (`min_depth`, `executor`, `container`, `ncbi_api_key`,
etc). `genetic_code` is irrelevant for rDNA and can be dropped or ignored.

Practical win: a user points RiboPilot and MitoPilot at the **same map file and
same `data_path`**, producing mito + rDNA markers from one skim.

References are built per-clade up front (see `02_ASSEMBLY`), not per-sample, so no
extra map-file columns are required and the shared MitoPilot map file works
unchanged. Optional per-sample `seed_ref` / `gene_ref` override columns can be
supported for odd samples, matching MitoPilot's `custom_seeds_db` /
`custom_labels_db` mechanism.

## Database schema

Keep MitoPilot's spine (`samples`, `preprocess`, `pre_opts`, plus opts tables
keyed by an opts hash) and replace the mito-specific tables. Proposed tables:

- `samples` - unchanged (`ID`, `Taxon`, metadata JSON, hidden flag, ...).
- `preprocess` / `pre_opts` - unchanged (fastp trimming, read counts, depth).
- `assemble` - per-sample assembly summary. Columns like MitoPilot's but
  rDNA-flavored: `length`, `topology` (linear monomer vs circularized repeat
  unit), `paths`, `scaffolds`, `assemble_notes`, `assemble_switch/lock`,
  `assemble_opts`, `blast_opts`, best-hit BLAST columns
  (`blast_accession/species/pident/qcovs/evalue/lineage`), `time_stamp`.
- `assemble_opts` - assembler settings. Instead of MitoFinder fields:
  `getOrganelle`, `seeds_db`, `labels_db`, `assembler` (getorganelle | rnabloom |
  spades-region | reference-guided), `k_list`, `target` (`nr` seed set),
  `max_paths`, `max_scaffolds`, `min_length`, `expected_length`.
- `assemblies` - per path/scaffold sequence rows (like MitoPilot: `ID`, `path`,
  `scaffold`, `sequence`, `length`, `depth`, `gc`, `ignore`, `edited`,
  `edit_positions`, BLAST hit cols). One repeat unit per path.
- `blast_opts` / `assembly_blast` - reuse MitoPilot's candidate-BLAST machinery
  (recent branch work) but against rDNA reference DBs (SILVA/EUKARYOME/PR2).
- `annotate` - per-assembly feature table (see 03). Rows are features
  (18S, ITS1, 5.8S, ITS2, 28S, optional ETS/28S D-domains) with `start`, `stop`,
  `strand`, `type` (rRNA | misc_RNA/spacer), `annotator` (barrnap | itsx |
  reference), `edited`, `notes`, curation flags.
- `annotate_opts` - annotator settings (barrnap kingdom/evalue, ITSx taxon set,
  reference profile choice).
- `curate_opts` - per-clade ruleset selection + tolerances.

Same DB-write discipline as MitoPilot: **writes happen driver-side over JDBC/DBI,
not from a head-node local executor**, and curate/validate commits assemblies +
annotations atomically in one transaction (carry over MitoPilot's
`curate-validate atomic writer` and `DB writes driver-side` lessons).

## Nextflow workflow

`inst/nextflow/main.nf` with modules mirroring MitoPilot's stage split:

```
preprocess  -> fastp trim, read/depth QC              (reuse MitoPilot module)
assemble    -> targeted nr assembly (see 02)
             -> orient + rotate to canonical start (18S 5' end)
             -> BLAST vs rDNA ref DB, pick candidate refs
annotate    -> barrnap (18S/5.8S/28S) + ITSx (ITS1/ITS2)
             -> reconcile boundaries into one feature model
curate      -> per-clade ruleset: anchor genes, derive spacers, flag oddities
validate    -> length/structure checks; write atomic to DB
export       -> GenBank tbl/fasta/source-modifiers
```

Executor/config system is **reused unchanged** from MitoPilot
(`generate_config()`, `list_configs()`, saved profiles, the `clusterOptions`
gotcha noted in MitoPilot memory). Container is a RiboPilot image bundling
GetOrganelle, barrnap/HMMER, ITSx, BLAST+, seqkit, etc.

## Shiny app

Same three-tab curation flow. Differences:

- Assemble tab: show the linear repeat-unit map with element blocks (gggenes),
  coverage track, path selector, BLAST candidate refs modal (reuse MitoPilot's
  "All BLAST Hits" modal work).
- Annotate tab: draggable boundary editing for the 4 internal breakpoints
  (18S|ITS1, ITS1|5.8S, 5.8S|ITS2, ITS2|28S). Because spacers are defined by the
  anchors, editing an anchor edge auto-updates the adjacent spacer.
- Export tab: same as MitoPilot, feature-table generation tuned for rRNA +
  misc_RNA/spacer features and rDNA source modifiers.

## Look-and-feel and naming parity (hard requirement)

RiboPilot should feel like MitoPilot to an existing user. Mirror MitoPilot
wherever the biology allows:

- **Visual styling** - same theme/CSS, colors, fonts, `bslib`/shiny theme,
  logo treatment, button styles, icon set (same `shiny::icon`/font-awesome
  symbols), table styling (same DT/reactable look), pickers and inputs
  (same `selectizeInput`/`pickerInput` choices), modal layout, spacing.
  Copy MitoPilot's `inst/app/www` CSS/assets and golden theme as the starting
  point and only diverge where a control has no rDNA analogue.
- **Component choices** - same widget for the same job (e.g. the path picker,
  the coverage plot, the "All BLAST Hits" modal, lock/switch/hide toggles,
  the config editor) so muscle memory transfers.
- **Function names** - keep MitoPilot's exported R function names where the role
  matches: `new_project()`, `new_test_project()`, `add_samples()`,
  `generate_config()`, `list_configs()`, `export()`/export helpers, the
  `curate_*` / `params_*` family. Same argument names too (`mapping_fn`,
  `data_path`, `executor`, `min_depth`, `container`, `full_size`, ...).
- **GUI labels / button text / column headers** - reuse MitoPilot's wording
  where it still applies (Assemble/Annotate/Export tabs; buttons like
  "Run", "Lock", "Hide", "Ignore", "Edit"; DB/table column names shared with
  MitoPilot such as `ID`, `Taxon`, `length`, `topology`, `paths`, `scaffolds`,
  `blast_*`, `assemble_switch/lock`, `edited`, `notes`, `time_stamp`).
- Only rename when the rDNA feature genuinely differs (e.g. rRNA/spacer feature
  labels, ITS boundary controls), and even then follow MitoPilot's naming style.

## What we can copy nearly verbatim from MitoPilot

- `generate_config.R`, all `inst/config.*` templates, executor resolution.
- `add_samples.R`, project scaffolding in `init_project.R`.
- Preprocess module (fastp) and coverage/depth plumbing.
- BLAST candidate-ref machinery (branch `wf1-multi-blast-ref`), retargeted to
  rDNA DBs.
- Shiny shell, golem config, DB-access helpers, atomic-writer pattern.
- Backwards-compatibility / DB-migration helper pattern.
