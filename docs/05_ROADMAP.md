# 05 - Roadmap and open questions

## Phased build plan

### Phase 0 - Scaffolding (reuse MitoPilot)
- Fork MitoPilot's non-biological spine: `new_project()`, `add_samples()`,
  `generate_config()` + all `inst/config.*`, golem Shiny shell, DB helpers,
  atomic writer, preprocess (fastp) Nextflow module.
- Stand up the rDNA DB schema (`01_ARCHITECTURE.md`).
- Build the RiboPilot Docker image: GetOrganelle, HMMER, barrnap/pyBarrnap, ITSx,
  BLAST+, seqkit, (optional Infernal), R + package deps.
- Deliverable: `new_project()` produces a DB + config; empty pipeline runs
  through preprocess.

### Phase 1 - Assembly (Option B1)
- Wrap **`go_fetch.py`** as a reference-prep step: per-sample seed + gene DBs from
  NCBI taxid, with a packaged SILVA/EUKARYOME fallback for missing taxids/offline
  runs and an Entrez-response cache. (Borrowed from skim2rrna; preferred over
  maintaining one giant pan-eukaryote DB.)
- GetOrganelle nr assembly module; orient + rotate to 18S 5' start; multi-path
  output + path scoring (port MitoPilot).
- BLAST candidate-ref module vs rDNA DB (port MitoPilot's multi-blast-ref branch);
  lineage/contaminant flags.
- Assemble tab in Shiny: element map, coverage, path picker, BLAST modal.
- Deliverable: reads -> curated repeat-unit candidate(s) in DB, reviewable.

### Phase 2 - Annotation
- barrnap module (18S/5.8S/28S) + ITSx module (ITS1/ITS2); reconcile into one
  gap-free feature model (`03`).
- Optional cmsearch/Rfam high-accuracy anchor mode.
- Annotate tab: draggable 4 breakpoints, reference-alignment/identity overlay.
- Deliverable: full feature model per assembly, editable.

### Phase 3 - Reference-guided curation
- `curate_rrna_core` + first clade wrappers/params (start with a well-sampled
  metazoan clade, e.g. fish or a mollusc group already used in MitoPilot tests).
- Anchor/length/divergence/contaminant flags; boundary snapping.
- Deliverable: automated pass with human-in-the-loop flag resolution.

### Phase 4 - Export + fallbacks
- GenBank tbl/fasta/source-modifier export; per-marker fasta bundles.
- Option B2 reference-guided baiting fallback; Option A whole-assembly cross-check.
- Backwards-compat/DB-migration helper.
- Deliverable: end-to-end run to submission-ready output.

### Phase 5 - Broaden taxa
- Add clade rulesets + reference sets outward (broad metazoa, then optionally
  plants/fungi/protists). Same "params + wrapper" extension pattern as MitoPilot.
- Optional: 5S array, ETS/IGS recovery, 28S D-domain sub-annotation.

## Open questions (decide before/early in build)

1. **Primary taxonomic scope for v1?** Broad metazoa is the natural first target
   (matches MitoPilot's user base). Plants/fungi already served partly by
   GetOrganelle presets + UNITE; include or defer?
2. **Circular vs linear output?** Lean linear 18S-start monomer, IGS excluded.
   Confirm this is what downstream barcoding submissions want.
3. **How much spacer curation is worth it?** Proposal: bound-and-trim only, no
   base-level spacer curation. Confirm acceptable scientifically.
4. **28S extent:** single feature vs D-domain sub-annotation by default?
   Proposal: single feature, D-domains optional.
5. **Reference DB packaging:** *leaning resolved* - use `go_fetch` per-sample
   fetch by taxid as primary (no giant DB to ship), plus a small packaged
   EUKARYOME/SILVA fallback in the image for missing-taxid/offline runs. Confirm
   this over shipping full trimmed subsets.
6. **Intragenomic variants:** report top-N paths (like MitoPilot) or force a
   single consensus? Proposal: top-N with dominant flagged primary.
7. **Shared repo vs sibling repo:** keep RiboPilot fully separate (this repo) or
   eventually factor MitoPilot's shared spine into a common base package both
   depend on? Start separate; revisit if drift is painful.
8. **5S rRNA:** in or out? Proposal: out for v1 (separate array).

## Guardrails carried from MitoPilot (do not relearn the hard way)

- DB writes happen **driver-side over DBI/JDBC**, never from a head-node local
  executor.
- curate/validate commits **assemblies + annotations atomically** in one txn.
- Watch the `INSERT OR REPLACE` column-wipe trap (list every column or unlisted
  ones go NULL).
- No reused gargoyle flag name as a `session$userData` payload key.
- Scope fixes narrowly (per-focal-feature recompute, not blanket rebuilds).
- `clusterOptions` config gotcha in the executor system.
