# Batch Circadian Analysis of Multiple AGD Files

Reads a set of ActiGraph `.agd` files, runs the full circadian analysis
on each, and returns one combined summary row per file, the package
equivalent of uploading a batch of files. Files that fail to read or
analyse are reported in an `error` column rather than aborting the run.
Optionally writes a single Excel workbook whose Summary sheet has one
row per file.

## Usage

``` r
circadian.batch(
  files,
  file = NULL,
  metric = c("axis1", "vm"),
  include_period_ci = FALSE,
  n_boot = 200,
  epoch_length = NULL,
  verbose = TRUE
)
```

## Arguments

- files:

  Character vector of `.agd` paths, or a single directory (all its
  `.agd` files are used).

- file:

  Optional output `.xlsx` path for a combined-summary workbook.

- metric:

  Activity metric: `"axis1"` (default) or `"vm"` (vector magnitude).

- include_period_ci:

  Bootstrap the period CI per file (default `FALSE` for batch speed).

- n_boot:

  Bootstrap replicates when `include_period_ci` is `TRUE`.

- epoch_length:

  Epoch length in seconds; `NULL` (default) infers it per file from the
  timestamps.

- verbose:

  Print per-file progress (default `TRUE`).

## Value

A data frame with one row per file (a `file` column, an `error` column,
and every summary metric). Returned invisibly when a workbook `file` is
written.

## See also

[`circadian.workbook`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)

## Examples

``` r
# \donttest{
# Pass a folder to analyse a whole batch; a single file is used here.
batch <- circadian.batch(example_agd(), verbose = FALSE)
batch[, c("file", "IS", "IV", "RA", "rhythm_p_value")]
#>                      file     IS     IV   RA rhythm_p_value
#> 1 MOS2E39230594_60sec.agd 0.2279 1.0008 0.98    0.002396228
# }
```
