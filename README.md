# RiboPilot

Nuclear ribosomal repeat (nrDNA) assembly, annotation, and curation from genome
skimming data. Sister package to
[MitoPilot](https://github.com/Smithsonian/MitoPilot).

RiboPilot targets the eukaryotic nuclear ribosomal cistron and its spacers:

```
   ETS   18S    ITS1   5.8S   ITS2    28S     IGS
 --|===|=======|~~~~~|=====|~~~~~|===========|------ (tandem repeat, ~100s of copies)
       SSU              LSU
```

Goal: recover the `18S - ITS1 - 5.8S - ITS2 - 28S` region (optionally with ETS/IGS
flanks), annotate gene and spacer boundaries, curate against references, and
export GenBank-ready submissions. Designed to run on the same input map file,
project layout, SQLite tracking database, Nextflow execution, and HPC executor
configs as MitoPilot, so the two can be run side by side on the same genome
skims.

**Status:** planning only. No code yet. See [`docs/`](docs/).

## Planning documents

- [`docs/00_OVERVIEW.md`](docs/00_OVERVIEW.md) - vision, scope, biology primer
- [`docs/01_ARCHITECTURE.md`](docs/01_ARCHITECTURE.md) - package layout, DB schema, Nextflow, reuse of MitoPilot
- [`docs/02_ASSEMBLY.md`](docs/02_ASSEMBLY.md) - assembly strategies and trade-offs
- [`docs/03_ANNOTATION_CURATION.md`](docs/03_ANNOTATION_CURATION.md) - annotation + reference-guided curation
- [`docs/04_EXISTING_TOOLS.md`](docs/04_EXISTING_TOOLS.md) - prior art, tools, reference databases
- [`docs/05_ROADMAP.md`](docs/05_ROADMAP.md) - phased build plan and open questions
