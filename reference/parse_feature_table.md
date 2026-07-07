# Parse an NCBI feature table (efetch -format ft) into rRNA/ITS features

Parse an NCBI feature table (efetch -format ft) into rRNA/ITS features

## Usage

``` r
parse_feature_table(ft_lines)
```

## Arguments

- ft_lines:

  Character vector (lines of the feature-table file).

## Value

data.frame(element, type, start, stop, strand); rRNA/ITS rows only.
