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
docker/deploy-local.sh            # tags ribopilot:<DESCRIPTION Version>
docker/deploy-local.sh v0.1.0     # tags ribopilot:v0.1.0
docker/deploy-local.sh --no-db    # build without the rDNA panel
```

Each script first builds the R package tarball into this directory
(`devtools::build(path = "docker")`), then runs `docker build`. Note it does
**not** run `devtools::document()`: RiboPilot's `NAMESPACE` is hand-maintained.

## The baked-in rDNA panel

By default both scripts stage the rDNA panel and its decoys from the
`db-<RP_DB_VERSION>` GitHub release (with `gh`, since the repo is private) and
bake them into the image at `RIBOPILOT_DB_DIR`, so a Docker run searches locally
with no configuration. If that release does not exist yet, the scripts stop and
name the commands that fix it.

`--no-db` builds without the panel. That image leaves `RIBOPILOT_DB_DIR` unset,
so nothing is discovered and every BLAST search goes to NCBI.

Setting `RIBOPILOT_NO_LOCAL_BLAST=1` in the shell that launches a run forces the
remote search even on an image that carries the panel.

## Push to Docker Hub

Publishes to `macguigand/ribopilot`.

```bash
docker/deploy-dockerhub.sh                  # tag = the DESCRIPTION Version
docker/deploy-dockerhub.sh v0.1.0           # explicit tag
docker/deploy-dockerhub.sh --dry-run        # build and test, push nothing
docker/deploy-dockerhub.sh --no-db          # build without the rDNA panel
```

The tag defaults to the DESCRIPTION `Version` on purpose: `new_project()`
derives its default `container` from the same place, so a tag chosen by hand can
be a tag no project ever asks for.

There is no `:latest`, and the scripts cannot publish one. A floating tag lets
two machines run different code while naming the same image, which is exactly
the ambiguity a pinned version removes.

Before pushing, the script:

- clears stale tarballs from `docker/` (the Dockerfile COPYs
  `RiboPilot_*.tar.gz` to one destination, so two of them break the build);
- builds for `linux/amd64` explicitly, because not every tool in this image has
  an arm64 conda build and an image that quietly matches the host architecture
  is an image that fails on someone else's cluster;
- labels the image with the commit it came from, and warns if the tree is dirty;
- runs a smoke test in the built image: every tool the pipeline calls must be on
  `PATH` and the package must load. Several pipeline steps degrade to a SILENT
  no-op when a tool is missing, so this is worth catching before publishing
  rather than hours into somebody's run;
- asks for confirmation before it pushes anything.

Set `DOCKERHUB_TOKEN` to skip the login prompt, or `DOCKERHUB_USER` to publish
under a different account. Use an access token rather than your password.

## Use the image in a project

Create the project with the `docker` executor and point it at the image tag you
built or pushed:

```r
RiboPilot::new_project(
  path       = "~/my_project",
  mapping_fn = "~/my_samples.csv",
  data_path  = "~/my_reads",
  executor   = "docker",
  container  = "macguigand/ribopilot:0.1.0"   # or a local tag, e.g. "ribopilot:0.0.0.9004"
)
```

This writes a `.config` from `inst/config.docker` that enables Docker and runs
every process inside the image. Because the pipeline opens the project SQLite
database directly by path, the config bind-mounts the project directory and the
raw-read directory into the container (`docker.runOptions`). If your reference
database lives outside the project directory, add another `-v` mount to
`docker.runOptions` in the project's `.config`.

## Notes

- Clade reference databases are not baked into the image (the rDNA BLAST panel
  is, see above). The default fish databases are fetched from GitHub at run
  time; user-built databases are read from disk (mount them, see above).
- Nextflow itself runs on the host (it launches the containers), so it is not
  installed in the image.
