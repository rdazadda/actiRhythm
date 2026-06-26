# Activity Counts from a Raw Axivity .cwa File

Reads a raw Axivity (`.cwa`) accelerometer file and converts it to
ActiGraph-equivalent activity counts via the agcounts band-pass filter
(the Brond 2017 method). Requires the GGIRread and agcounts packages.

## Usage

``` r
axivity.counts(path, epoch = 60, lfe = FALSE, tz = "UTC")
```

## Arguments

- path:

  Path to a `.cwa` file.

- epoch:

  Epoch length in seconds (default 60).

- lfe:

  Use the low-frequency extension filter (default `FALSE`).

- tz:

  Time zone for the timestamps (default `"UTC"`).

## Value

Data frame with `time`, `axis1`, `axis2`, `axis3` and `vm`, one row per
epoch (the same shape as
[`gt3x.counts`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md)).

## Cross-brand counts

These are an *approximation* of ActiGraph counts, not native ActiGraph
output. Axivity-to-count conversion has been directly validated (Brond
et al. 2017). The result is appropriate for the relative and normalized
circadian metrics in this package (IS, IV, RA, L5, M10, SRI); it should
*not* be used to apply ActiGraph intensity cut-points or to compare
absolute counts across device brands.

## References

Brond JC, Andersen LB, Arvidsson D (2017). “Generating ActiGraph counts
from raw acceleration recorded by an alternative monitor.” *Medicine &
Science in Sports & Exercise*, **49**(11), 2351–2360.
[doi:10.1249/MSS.0000000000001344](https://doi.org/10.1249/MSS.0000000000001344)
.

## See also

[`geneactiv.counts`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md),
[`gt3x.counts`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md),
[`read.raw`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
