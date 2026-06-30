# Compare Cosinor Rhythms Between Two Groups

Tests whether the rest-activity rhythm differs between two groups. Fits
each subject's cosinor with the same engine as
[`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md).
An Bingham et al. (1982) comparison is a 2-variate Hotelling's \\T^2\\
on the cosine and sine coefficients (testing equal amplitude and
acrophase, accounting for their correlation), with the MESOR compared by
a separate pooled-variance F-test. Auxiliary per-parameter two-sample
(Welch) t-tests break the result down for the MESOR, amplitude, and
acrophase; the acrophase test is circular-aware, unwrapping the
per-subject acrophases about their common circular mean.

## Usage

``` r
cosinor.compare(
  activity,
  timestamps,
  subject,
  group,
  period = 24,
  level = 0.95,
  min_valid_hours = 12
)
```

## Arguments

- activity:

  Numeric activity vector with all subjects stacked.

- timestamps:

  POSIXct timestamps, one per value.

- subject:

  Subject identifier, one per value.

- group:

  Group identifier (exactly two levels), one per value.

- period:

  Rhythm period in hours (default 24).

- level:

  Confidence level for the difference CIs (default 0.95).

- min_valid_hours:

  Minimum profile hours for a subject to be included.

## Value

An object of class `actiRhythm_cosinor_compare`: a `joint` list (the
amplitude/acrophase Hotelling \\T^2\\ with its F, degrees of freedom and
p-value, plus a separate MESOR F-test), a `tests` data frame of
auxiliary per-parameter Welch comparisons (each group's estimate, the
difference, the t statistic, p-value and CI), and the per-subject
coefficients.

## References

Bingham C, Arbogast B, Cornelissen Guillaume G, Lee JK, Halberg F
(1982). “Inferential statistical methods for estimating and comparing
cosinor parameters.” *Chronobiologia*, **9**(4), 397–439.

## See also

[`population.cosinor`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)

## Examples

``` r
set.seed(1); hrs <- 0:23
act <- ts <- subj <- grp <- NULL
for (g in c("A", "B")) for (i in 1:5) {
  acro <- if (g == "A") 8 else 11
  y <- 100 + 40 * cos(2 * pi * (hrs - acro) / 24) + rnorm(24, 0, 4)
  act <- c(act, y); subj <- c(subj, rep(paste0(g, i), 24)); grp <- c(grp, rep(g, 24))
  ts <- c(ts, as.POSIXct("2024-01-01", tz = "UTC") + hrs * 3600)
}
cosinor.compare(act, as.POSIXct(ts, tz = "UTC", origin = "1970-01-01"), subj, grp)
#> Cosinor Comparison Between Two Groups
#> 
#>   Groups:  A (n=5)  vs  B (n=5)
#>   Period:  24 h
#> 
#>   Amplitude/acrophase (Bingham T2):  F(2,7) = 994.58   p = 2.55e-09
#>   MESOR (Bingham F):                 F(1,8) = 1.70   p = 0.228
#> 
#>   mesor     diff = +0.70   t(5.6) = 1.31   p = 0.243
#>   amplitude diff = -0.22   t(7.8) = -0.34   p = 0.743
#>   acrophase diff = -2.92   t(6.8) = -39.92   p = 2.77e-09
```
