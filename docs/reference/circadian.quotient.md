# Circadian Quotient and Relative Amplitude from a Cosinor Fit

Computes the circadian quotient (amplitude divided by MESOR) and a
cosinor relative amplitude (amplitude divided by the overall mean) from
a cosinor result. The circadian quotient is a dimensionless measure of
rhythm strength against the rhythm-adjusted mean.

## Usage

``` r
circadian.quotient(cosinor_result)
```

## Arguments

- cosinor_result:

  A list returned by
  [`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  (class `"actiRhythm_cosinor"`) or any list providing `amplitude` and
  `mesor`. An optional `overall_mean` (raw arithmetic mean of the
  counts) is used for the relative-amplitude denominator when present;
  otherwise `mesor` is used.

## Value

A list with class `"actiRhythm_circadian_quotient"` containing:

- circadian_quotient:

  `amplitude / mesor`.

- relative_amplitude:

  `amplitude / overall_mean` (falls back to `amplitude / mesor` when
  `overall_mean` is absent).

Returns `NA` entries when the inputs are missing or the denominator is
non-positive.

## References

Nelson W, Tong YL, Lee JK, Halberg F (1979). Methods for
cosinor-rhythmometry. *Chronobiologia*, 6(4):305-323.

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440 * 3)
counts <- 100 + 80 * cos(2 * pi * (as.numeric(format(ts, "%H")) - 14) / 24)
cos <- cosinor.analysis(counts, ts)
circadian.quotient(cos)
#> $circadian_quotient
#> [1] 0.8
#> 
#> $relative_amplitude
#> [1] 0.8
#> 
#> attr(,"class")
#> [1] "actiRhythm_circadian_quotient"
# }
```
