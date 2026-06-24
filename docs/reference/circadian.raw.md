# Circadian Analysis Directly from Raw Acceleration

From a raw accelerometer file (or raw data frame), computes a per-epoch
raw metric (ENMO or MAD) with auto-calibration, then runs the
count-style circadian analysis
([`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md))
on it. The metrics (IS, IV, RA, L5/M10, cosinor, DFA, periodograms, SRI,
...) run on the raw metric the same way they run on counts, so raw data
needs no new per-method code.

## Usage

``` r
circadian.raw(
  x,
  metric = c("ENMO", "MAD"),
  device = "auto",
  epoch = 60,
  calibrate = TRUE,
  tz = "UTC",
  ...
)
```

## Arguments

- x:

  A path to a raw file (`.gt3x`, `.cwa`, `.bin`) or a raw data frame
  (see
  [`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
  /
  [`example_raw`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)).

- metric:

  Raw metric to analyse, `"ENMO"` (default) or `"MAD"`.

- device:

  Device brand or `"auto"` (default; used for file input).

- epoch:

  Epoch length in seconds (default 60).

- calibrate:

  Apply van Hees auto-calibration first (default `TRUE`).

- tz:

  Time zone (default `"UTC"`).

- ...:

  Passed to
  [`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md).

## Value

The
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
result computed on the raw metric.

## See also

[`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md),
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
[`example_raw`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)

## Examples

``` r
# The full count-style battery on synthetic raw ENMO
# \donttest{
cr <- circadian.raw(example_raw(days = 2), metric = "ENMO")
c(IS = cr$IS, IV = cr$IV, RA = cr$RA)
#>     IS     IV     RA 
#> 1.0000 0.3823 0.9873 
# }
```
