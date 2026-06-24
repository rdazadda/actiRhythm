# Activity Counts from a Raw .gt3x File

Compute activity counts from a raw ActiGraph `.gt3x` accelerometer file
using the agcounts implementation of the ActiGraph count algorithm
(Neishabouri 2022). Requires the agcounts and read.gt3x packages.

## Usage

``` r
gt3x.counts(path, epoch = 60, lfe = FALSE, tz = "UTC")
```

## Arguments

- path:

  Path to a `.gt3x` file.

- epoch:

  Epoch length in seconds (default 60).

- lfe:

  Use the low-frequency extension filter (default `FALSE`).

- tz:

  Time zone for the timestamps (default `"UTC"`).

## Value

Data frame with `time`, `axis1`, `axis2`, `axis3` and `vm`, one row per
epoch.

## References

Neishabouri A, et al. (2022). Quantification of acceleration as activity
counts in ActiGraph wearable. *Scientific Reports*, 12:11958.

## See also

[`axivity.counts`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md),
[`geneactiv.counts`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md),
[`read.raw`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
