
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RiboPilot <a href="https://dmacguigan.github.io/RiboPilot/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

RiboPilot recovers the **nuclear ribosomal repeat (nrDNA)** from genome
skimming data, then helps you check and clean up the result and package
it for submission to GenBank. It is a companion to
[MitoPilot](https://github.com/Smithsonian/MitoPilot), which does the
same thing for the mitochondrial genome. You can run both on the same
samples and get both marker sets from one sequencing run.

You do not need to be a programmer to use RiboPilot. Most of your work
happens in a **point-and-click app** that runs in your web browser. The
step-by-step [Get started
guide](https://dmacguigan.github.io/RiboPilot/articles/RiboPilot.html)
walks through a complete project from start to finish.

# What RiboPilot recovers

Every cell carries the ribosomal RNA genes in a repeated unit that is
copied hundreds to thousands of times in a row. Because there are so
many copies, even a shallow, inexpensive “genome skim” sequences this
region very deeply, essentially for free. RiboPilot reconstructs that
repeat unit:

``` 
   ETS   18S    ITS1   5.8S   ITS2    28S     IGS
 --|===|=======|~~~~~|=====|~~~~~|===========|------ (repeated hundreds of times)
       SSU              LSU
```

From it you get the most widely used nuclear markers for phylogenetics
and DNA barcoding: **18S**, **28S**, **ITS1**, **ITS2**, and the full
**ITS** region.

The three ribosomal genes (18S, 5.8S, 28S) are highly conserved and act
as anchors. The two spacers (ITS1, ITS2) sit between them, evolve
quickly, and are the regions that distinguish closely related species.

# What RiboPilot does for you

Starting from raw sequencing reads, RiboPilot:

1.  **Confirms each sample’s organism.** RiboPilot resolves the `Taxon`
    label you provide to an NCBI taxon (taxid + full lineage) using a
    local copy of the NCBI taxonomy, then lets you **verify or correct
    it in the app** before anything else runs. This is a dedicated first
    stage (the **Taxonomy** tab): samples must be reviewed and locked
    here before they move on to assembly, so the taxonomy that guides
    scoring and curation is one you have signed off on. The taxonomy
    snapshot (\~72 MB) downloads once per project and is indexed the
    first time it is used, which takes about a minute and needs no
    software beyond the R package.
2.  **Assembles** the ribosomal repeat for each sample, cleaning up the
    reads and piecing the region back together, then picking the best
    reconstruction automatically.
3.  **Annotates** it by locating the boundaries of all five features
    (18S, ITS1, 5.8S, ITS2, 28S).
4.  **Checks (curates)** each result against a reference and flags
    anything that looks off, such as a missing gene, an unusual length,
    or a possible contaminant, so you can focus your attention where it
    is needed.
5.  **Exports** a tidy, GenBank-ready submission package.

Behind the scenes RiboPilot uses established bioinformatics tools
([fastp](https://github.com/OpenGene/fastp),
[GetOrganelle](https://github.com/Kinggerm/GetOrganelle),
[barrnap](https://github.com/tseemann/barrnap) /
[pybarrnap](https://github.com/moshi4/pybarrnap),
[ITSx](https://microbiology.se/software/itsx/),
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), and others) and
stitches them together for you. You do not need to learn or run any of
them yourself.

# Which organisms it works with

RiboPilot was first built for ray-finned fishes and comes ready to use
for them out of the box. For any other group of organisms, it can build
the reference data it needs automatically by downloading published
sequences from GenBank. This is a one-time setup step per group; see the
[reference database
guide](https://dmacguigan.github.io/RiboPilot/articles/Reference-Databases.html).

It also runs automatic quality checks tuned for a taxonomic group.
Checks currently ship for a generic metazoan (animal) setting and a
fish-specific setting; see the [curation
guide](https://dmacguigan.github.io/RiboPilot/articles/Curation-Flags.html).

# Installation

**Supported execution: Docker or Singularity.** The analysis steps are
supported only inside the RiboPilot container image, which bundles every
bioinformatics tool the pipeline calls at the versions it was tested
with. The image is published on Docker Hub as
[`macguigand/ribopilot`](https://hub.docker.com/r/macguigand/ribopilot),
tagged with the RiboPilot version, and `new_project()` picks the
matching tag for you. On a workstation, create projects with `executor =
"docker"`. On a cluster, use Singularity: the built-in `NMNH_Hydra` and
`NOAA_SEDNA` profiles already do, and `generate_config()` writes a
Singularity config for any other scheduler.

Installing the tools yourself, whether one at a time or through the
`conda` environment described below, is **not supported going forward**.
Those instructions are left in place because they still work, but they
are not tested and problems that come from them cannot be supported.

RiboPilot needs three things installed on your computer (or your lab’s
server):

1.  **[R](https://www.r-project.org/)** (version 4.4.0 or newer), the
    language RiboPilot is written in. You will type a few short commands
    into R, but you do not need to know R.
2.  **[Nextflow](https://www.nextflow.io/docs/latest/install.html)**
    (version 24.10.x through 25.10.x), which runs the heavy analysis
    steps. RiboPilot refuses to load on an older Nextflow and prompts
    you to update; on a newer, unverified Nextflow it automatically pins
    runs to a supported version via `NXF_VER` (the `conda` environment
    below already installs a compatible one).
3.  **The bioinformatics tools** the pipeline uses, and only if you are
    not running in a container. The container already carries them, so
    with Docker or Singularity you can skip this step. **This conda
    route is no longer supported** (see above); it is kept for setups
    that already rely on it. The easiest way to get all of the tools at
    once is
    [conda](https://docs.conda.io/projects/miniconda/en/latest/), a tool
    that installs a bundle of programs into a self-contained
    “environment” without interfering with anything else on your
    machine. Once conda is installed, one command sets up everything
    RiboPilot needs:

<!-- end list -->

``` bash
conda env create -f environment.yml
conda activate ribopilot
```

Then, in R, install RiboPilot itself:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("dmacguigan/RiboPilot")
```

If any of this is unfamiliar, the [Get started
guide](https://dmacguigan.github.io/RiboPilot/articles/RiboPilot.html)
explains each piece in more detail, and your local bioinformatics
support staff can usually help with the one-time setup.

# Getting started

The best way to learn RiboPilot is to run the built-in practice project.
RiboPilot ships a small set of fish samples and a single command that
downloads them and sets everything up:

``` r
RiboPilot::new_test_project(path = "~/ribopilot_test", executor = "docker")
```

`executor` is passed explicitly here because `new_test_project()` still
defaults to `executor = "conda"`, which is no longer a supported way to
run.

From there, the [Get started
guide](https://dmacguigan.github.io/RiboPilot/articles/RiboPilot.html)
walks you through assembling the samples, reviewing and cleaning them up
in the app, and exporting the results.

## On the NMNH Hydra cluster

RStudio Server sessions on Hydra start without the grid engine or Java
on `PATH`, so Nextflow cannot find `qsub` and every job fails at
submission. Run this once per session, before launching the app:

``` r
RiboPilot::hydra_setup()
```

# Getting help

  - Full documentation: <https://dmacguigan.github.io/RiboPilot/>
  - Questions or problems: please open an
    [issue](https://github.com/dmacguigan/RiboPilot/issues).
