# Package index

## Start a project

Create a project, open the app, and add or update samples.

- [`new_test_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_test_project.md)
  : Initialize a RiboPilot test project
- [`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md)
  : Initialize a new RiboPilot project
- [`run_app()`](https://nmnh-oceandna.github.io/RiboPilot/reference/run_app.md)
  [`RiboPilot()`](https://nmnh-oceandna.github.io/RiboPilot/reference/run_app.md)
  : Open the RiboPilot app
- [`add_samples()`](https://nmnh-oceandna.github.io/RiboPilot/reference/add_samples.md)
  : Add samples to an existing project
- [`update_sample_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/update_sample_metadata.md)
  : Update metadata for existing samples
- [`backwards_compatibility()`](https://nmnh-oceandna.github.io/RiboPilot/reference/backwards_compatibility.md)
  : Update an old project to the installed RiboPilot

## Sample metadata

Fetch specimen records from GEOME, GBIF, and NCBI, and choose what to
compare.

- [`fetch_geome()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_geome.md)
  : Fetch GEOME metadata for project samples
- [`fetch_gbif()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_gbif.md)
  : Fetch GBIF occurrence metadata for project samples
- [`fetch_ncbi()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_ncbi.md)
  : Fetch NCBI BioSample and BioProject metadata for project samples
- [`remove_metadata()`](https://nmnh-oceandna.github.io/RiboPilot/reference/remove_metadata.md)
  : Remove fetched GEOME, GBIF, and NCBI metadata
- [`set_metadata_columns()`](https://nmnh-oceandna.github.io/RiboPilot/reference/set_metadata_columns.md)
  : Choose which mapping-file columns are compared with GEOME, GBIF, and
  NCBI

## Reference databases

Build a reference database from R, and find the shared reference folder.

- [`build_reference_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/build_reference_db.md)
  : Build a clade nuclear-rDNA reference database (GetOrganelle seeds +
  labels)
- [`ribopilot_ref_dir()`](https://nmnh-oceandna.github.io/RiboPilot/reference/ribopilot_ref_dir.md)
  : Default directory for user-built reference databases

## Clusters

Set up a profile for your cluster’s scheduler.

- [`generate_config()`](https://nmnh-oceandna.github.io/RiboPilot/reference/generate_config.md)
  : Generate and save a reusable cluster config profile
- [`hydra_setup()`](https://nmnh-oceandna.github.io/RiboPilot/reference/hydra_setup.md)
  : Configure the R session environment for the NMNH Hydra cluster

## Export

Write GenBank submission files without the app.

- [`export_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/export_project.md)
  : Export a project's locked units as GenBank bundles
