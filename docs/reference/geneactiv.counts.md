# Activity Counts from a Raw GENEActiv .bin File

Reads a raw GENEActiv (`.bin`) accelerometer file and converts it to
ActiGraph-equivalent activity counts via the agcounts band-pass filter
(the Brond 2017 method). Requires the GGIRread and agcounts packages.

## Usage

``` r
geneactiv.counts(path, epoch = 60, lfe = FALSE, tz = "UTC")
```

## Arguments

- path:

  Path to a `.bin` file.

- epoch:

  Epoch length in seconds (default 60).

- lfe:

  Use the low-frequency extension filter (default `FALSE`).

- tz:

  Time zone for the timestamps (default `"UTC"`).

## Value

Data frame with `time`, `axis1`, `axis2`, `axis3` and `vm`, one row per
epoch.

## Cross-brand counts

These are an *approximation* of ActiGraph counts, not native ActiGraph
output, and GENEActiv-to-count conversion is *not* empirically validated
(only theoretically motivated by the shared filter). GENEActiv records
at about 85.7 Hz, which is resampled to 30 Hz before the count filter.
Appropriate for the relative and normalized circadian metrics here (IS,
IV, RA, L5, M10, SRI); do *not* apply ActiGraph cut-points or compare
absolute counts across brands.

## References

Brond JC, Andersen LB, Arvidsson D (2017). Generating ActiGraph counts
from raw acceleration recorded by an alternative monitor. *Medicine &
Science in Sports & Exercise*, 49(11):2351-2360.
[doi:10.1249/MSS.0000000000001344](https://doi.org/10.1249/MSS.0000000000001344)

## See also

[`axivity.counts`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md),
[`gt3x.counts`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md),
[`read.raw`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
