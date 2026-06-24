# Cosinor Rhythmicity Test

Tests whether a rhythm of the given period is present, using the Halberg
zero-amplitude F-test (H0: amplitude = 0) on the single-component
cosinor fit, and reports the percent rhythm (the proportion of variance
the cosinor explains, R^2). It reuses
[`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
as the single fitting engine, so the MESOR/amplitude/acrophase and
valid-day gating are identical.

## Usage

``` r
rhythmicity.test(
  counts,
  timestamps,
  period = 24,
  alpha = 0.05,
  wear_time = NULL,
  min_valid_hours = 10,
  cosinor_result = NULL
)
```

## Arguments

- counts:

  Numeric vector of activity counts (or any activity measure).

- timestamps:

  POSIXct timestamps, one per count.

- period:

  Rhythm period in hours (default 24).

- alpha:

  Significance level for the `rhythmic` flag (default 0.05).

- wear_time:

  Optional logical wear-time mask passed to the cosinor fit.

- min_valid_hours:

  Minimum valid hours per day for the cosinor fit.

- cosinor_result:

  Optional existing
  [`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  result to test without refitting; when supplied, `counts`/`timestamps`
  are ignored.

## Value

An object of class `actiRhythm_rhythmicity`: a list with the F
statistic, numerator/denominator degrees of freedom (`df1`, `df2`),
`p_value`, `percent_rhythm` (and `r_squared`), a logical `rhythmic`
flag, and the cosinor parameters.

## References

Nelson W, Tong YL, Lee JK, Halberg F (1979). Methods for
cosinor-rhythmometry. Chronobiologia, 6(4), 305-323.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440)
h <- as.numeric(format(ts, "%H"))
counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
rhythmicity.test(counts, ts)
#> Cosinor Rhythmicity Test (Halberg zero-amplitude F-test)
#> 
#> H0: amplitude = 0 (no rhythm)
#>   Period:         24 h
#>   F(2, 21):        Inf
#>   P-value:        <1e-300
#>   Percent rhythm: 100.0% (R-squared = 1.0000)
#>   Rhythmic:       YES (alpha = 0.05)
```
