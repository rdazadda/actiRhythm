# Activity Counts from a Raw Accelerometer File (Any Supported Brand)

A single entry point that dispatches to the brand-specific reader by
file extension (or an explicit `device`) and returns
ActiGraph-equivalent counts: ActiGraph `.gt3x`, Axivity `.cwa`, or
GENEActiv `.bin`.

## Usage

``` r
read.raw(
  path,
  device = c("auto", "gt3x", "axivity", "geneactiv"),
  epoch = 60,
  lfe = FALSE,
  tz = "UTC"
)
```

## Arguments

- path:

  Path to a raw file.

- device:

  One of `"auto"` (infer from the extension), `"gt3x"`, `"axivity"`,
  `"geneactiv"`.

- epoch:

  Epoch length in seconds (default 60).

- lfe:

  Use the low-frequency extension filter (default `FALSE`).

- tz:

  Time zone for the timestamps (default `"UTC"`).

## Value

Data frame with `time`, `axis1`, `axis2`, `axis3` and `vm`, one row per
epoch. Counts from non-ActiGraph devices are an approximation; see
[`axivity.counts`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md)
/
[`geneactiv.counts`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md).

## See also

[`gt3x.counts`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md),
[`axivity.counts`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md),
[`geneactiv.counts`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md)
