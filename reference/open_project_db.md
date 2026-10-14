# Open a project DB connection with a busy timeout and a filesystem-appropriate journal mode

The busy timeout is set first: `journal_mode` takes an exclusive lock,
so setting it before the timeout would run that one statement with no
timeout.

## Usage

``` r
open_project_db(db)
```

## Arguments

- db:

  Path to the project SQLite file.
