# RiboPilot Docker image

This directory builds the container image that runs the RiboPilot Nextflow
pipeline. The image bundles R + the RiboPilot package + the full nrDNA tool chain
(fastp, GetOrganelle, SPAdes, BLAST, bowtie2, minimap2, samtools, barrnap, ITSx,
HMMER, seqkit/seqtk, Entrez Direct, sqlite3, ...) on a single `PATH`. Unlike the
MitoPilot image, everything lives in one conda environment, because the RiboPilot
pipeline calls each tool directly rather than through `conda run -n <env>`.

## Files

- `Dockerfile` - the image definition.
- `Dockerfile.dockerignore` - trims the build context to just the package tarball.
- `deploy-local.sh` - build the image locally.
- `deploy-dockerhub.sh` - build and push to Docker Hub.
- `RiboPilot_*.tar.gz` - the package tarball the deploy scripts build for the
  Dockerfile to install (git-ignored; regenerated on each build).

## Build locally

From the repository root:

```bash
docker/deploy-local.sh            # tags ribopilot:latest
docker/deploy-local.sh v0.1.0     # tags ribopilot:v0.1.0
```

Each script first builds the R package tarball into this directory
(`devtools::build(path = "docker")`), then runs `docker build`. Note it does
**not** run `devtools::document()`: RiboPilot's `NAMESPACE` is hand-maintained.

## Push to Docker Hub

```bash
docker/deploy-dockerhub.sh v0.1.0
```

It prompts for your Docker Hub credentials, then builds, tags
`<username>/ribopilot:v0.1.0`, and pushes.

## Use the image in a project

Create the project with the `docker` executor and point it at the image tag you
built or pushed:

```r
RiboPilot::new_project(
  path       = "~/my_project",
  mapping_fn = "~/my_samples.csv",
  data_path  = "~/my_reads",
  executor   = "docker",
  container  = "macguigand/ribopilot:v0.1.0"   # or a local tag, e.g. "ribopilot:latest"
)
```

This writes a `.config` from `inst/config.docker` that enables Docker and runs
every process inside the image. Because the pipeline opens the project SQLite
database directly by path, the config bind-mounts the project directory and the
raw-read directory into the container (`docker.runOptions`). If your reference
database lives outside the project directory, add another `-v` mount to
`docker.runOptions` in the project's `.config`.

## Notes

- Reference databases are not baked into the image. The default fish databases
  are fetched from GitHub at run time; user-built databases are read from disk
  (mount them, see above).
- Nextflow itself runs on the host (it launches the containers), so it is not
  installed in the image.
