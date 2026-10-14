# Smithsonian Hydra setup

This page covers what is specific to the Smithsonian Hydra cluster. Once
you are set up here, continue with the [Get
started](https://nmnh-oceandna.github.io/RiboPilot/articles/RiboPilot.md)
walkthrough, but using `executor = "NMNH_Hydra"` for the
[`new_test_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_test_project.md)
function. [Running on a
cluster](https://nmnh-oceandna.github.io/RiboPilot/articles/Clusters.md)
has general instructions for using RiboPilot on a cluster.

If you need a Hydra account, see the [Hydra policies
page](https://confluence.si.edu/display/HPC/Hydra+Policies).

## Get an R session on Hydra

The simplest route is the dedicated RStudio Server at
<https://galaxy.si.edu/R4>. Log in from a browser with your Hydra
credentials. This is a shared resource, so use it for interactive work
and submit larger runs as jobs. Details are on the [Hydra RStudio
page](https://confluence.si.edu/display/HPC/Dedicated+RStudio+Server).

## Install Nextflow

Hydra has no Nextflow module, so install your own copy once. Run these
commands in a **terminal** connected to Hydra:

``` bash
cd ~
module load tools/java/21.0.2
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p ~/bin
mv nextflow ~/bin/nextflow
echo 'export PATH="${HOME}/bin:${PATH}"' >> ~/.bashrc
source ~/.bashrc
```

Job scripts written by the app load Java and read `~/.bashrc` for you.
If you want to use Nextflow from a shell yourself, run
`module load tools/java/21.0.2` first.

## Install RiboPilot

In the R console of your RStudio Server session, follow step 5 of
[Install the
tools](https://nmnh-oceandna.github.io/RiboPilot/articles/Installation.html#install-the-tools).
Singularity is already on Hydra, so there is nothing else to install.

## Create a project

``` r

RiboPilot::new_test_project(
  path = "/pool/public/genomics/<USER>/RiboPilot/test_project",
  executor = "NMNH_Hydra"
)
```

Use the same `executor = "NMNH_Hydra"` with
[`RiboPilot::new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md)
for your own data. Work under `/pool` or `/scratch`. Home folders are
too small, and projects cannot run from `/store`.

By default Singularity’s download cache and the reference folder also
sit in your home folder. [Disk space, caches, and offline
use](https://nmnh-oceandna.github.io/RiboPilot/articles/Caches-and-Offline.md)
shows how to move them.

## Open the app and run

``` r

RiboPilot::RiboPilot("/pool/public/genomics/<USER>/RiboPilot/test_project")
```

The **Update** run window shows a job script that is sent to Hydra’s
queue for workflow managers. Press **Submit as job**, or **Save script**
to submit it yourself with `qsub`. The job keeps going after you close
the browser. Watch it with `qstat` in a terminal, then press the reload
button above the table to see new results. [Submit a run as a
job](https://nmnh-oceandna.github.io/RiboPilot/articles/Clusters.html#submit-a-run-as-a-job)
has the details.

The [Get
started](https://nmnh-oceandna.github.io/RiboPilot/articles/RiboPilot.md)
walkthrough marks each place where Hydra differs, in boxes tagged HYDRA.
