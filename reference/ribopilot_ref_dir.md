# Default directory for user-built reference databases

A per-user, cross-project location (override with the
\`RIBOPILOT_REF_DIR\` environment variable). Kept outside project
directories so a clade DB built once is reused by every project. Uses
\[tools::R_user_dir()\].

## Usage

``` r
ribopilot_ref_dir()
```
