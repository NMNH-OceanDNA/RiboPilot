# Running on a cluster

On a computing cluster, RiboPilot sends each step of the pipeline to the
cluster’s scheduler as a job, and runs the analysis tools in Singularity
(or Apptainer). You still use the app to review results and export.

## Choose a profile

A profile tells RiboPilot/Nextflow how to submit jobs on your cluster.
Two are built in and ready to use: `executor = "NMNH_Hydra"` and
`executor = "NOAA_SEDNA"`. Each has its own setup page: [Smithsonian
Hydra
setup](https://nmnh-oceandna.github.io/RiboPilot/articles/NMNH-Hydra.md)
and [NOAA SEDNA
setup](https://nmnh-oceandna.github.io/RiboPilot/articles/NOAA-SEDNA.md).

For any other Slurm, SGE, PBS, or LSF cluster, make a Nextflow profile
once with
[`RiboPilot::generate_config()`](https://nmnh-oceandna.github.io/RiboPilot/reference/generate_config.md),
then pass its name to
[`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md):

``` r

RiboPilot::generate_config(
  "mylab_slurm",
  scheduler       = "slurm",
  queue           = "general",
  account         = "mylab",
  walltime        = "48h",
  container_cache = "/scratch/mylab/singularity"
)
RiboPilot::new_project(path = "/shared/mylab/my_project",
  mapping_fn = "/shared/mylab/my_samples.csv",
  data_path = "/shared/mylab/reads", 
  executor = "mylab_slurm")
```

`container_cache` is the folder where Singularity keeps the RiboPilot
image, shared by every project that uses the profile. If you leave it
out, RiboPilot uses `NXF_SINGULARITY_CACHEDIR` when that is set, and
otherwise a `containers` folder inside your reference folder. Many
clusters need a `walltime`: some reject a job that asks for none.

Profiles are saved in your own settings folder (on Linux,
`~/.config/R/RiboPilot`). To share one with your lab, save it in a
shared folder with `profile_dir`, and pass the same `profile_dir` to
[`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md):

``` r

RiboPilot::generate_config(
  "mylab_slurm",
  scheduler   = "slurm",
  queue       = "general",
  walltime    = "48h",
  profile_dir = "/shared/mylab/ribopilot_profiles"
)
RiboPilot::new_project(
  path        = "/shared/mylab/my_project",
  mapping_fn  = "/shared/mylab/my_samples.csv",
  data_path   = "/shared/mylab/reads",
  executor    = "mylab_slurm",
  profile_dir = "/shared/mylab/ribopilot_profiles"
)
```

**Manual profile edits.** You can manually adjust profiles in the
`profile_dir`. This may be necessary for clusters that use a customized
job scheduler. If you are unsure how to do this, please consult your
cluster’s job submission documentation or contact a cluster admin. For
additional guidance, check [NF-core](https://nf-co.re/configs/) to see
if your cluster already has a custom config or instructions.

## Put everything on shared storage

The project folder, the reads, the reference folder, and the R library
that holds RiboPilot must be visible from the login node and from every
compute node.

## Open the app

On **RStudio Server**, just run
`RiboPilot::RiboPilot("/shared/mylab/my_project")`.

Without RStudio Server, run R and start the app on the login node or an
interactive session to reach it from your own computer through an SSH
tunnel:

``` r

RiboPilot::RiboPilot("/shared/mylab/my_project", port = 8123,
                     launch.browser = FALSE)
```

The app prints the address to open and the tunnel command. On your own
computer, in a terminal:

``` bash
ssh -N -L 8123:localhost:8123 you@login.cluster.edu
```

Then open <http://localhost:8123> in your browser.

When a batch scheduler (sbatch, qsub, or bsub) is on the PATH, or the
app runs on RStudio Server, the run window offers a job script to submit
instead of **Run from App**. The `headless` argument of
[`RiboPilot()`](https://nmnh-oceandna.github.io/RiboPilot/reference/run_app.md)
overrides this:

``` r

# Always offer a job script, even when no scheduler is detected
RiboPilot::RiboPilot("/shared/mylab/my_project", headless = TRUE)

# Always offer Run from App: Nextflow runs inside the R session that
# started the app, not through the scheduler
RiboPilot::RiboPilot("/shared/mylab/my_project", headless = FALSE)
```

`headless = FALSE` only makes sense when the machine running R is
allowed to do the work itself, for example an interactive compute
session. Most clusters forbid heavy jobs on login nodes.

## Submit a run as a job

Press **Update**. The run window shows the job script for the Nextflow
head job. The head job is small but long-lived: 2 CPUs, 8 GB of memory,
and 48 hours, in the queue and with the account from your profile. Check
the partition, account, walltime, and memory, and edit the script in the
window if you need to. Then press **Submit as job**, or **Save script**
to submit it yourself.

The job runs on its own, so you can close the browser and watch the job
with your scheduler’s tools (e.g. `qstat` on Hydra). Once job is done,
use the reload button above the table to see new results. The app does
not track a submitted job, so do not start another run for the same
project until the job has finished.

Each submission leaves three kinds of files in the project’s `.runs`
folder:

- `.runs/jobs/`: the job scripts;
- `.runs/logs/`: the scheduler’s log for each job;
- `.runs/nextflow/`: one Nextflow log per run.

Once a job has finished, use **Run Reports** in the toolbar to see how
it went.

In the Assemble step, **Work Dirs** finds the pipeline’s working folders
for a sample, from both app runs and cluster jobs.
