# Confidence Interval for the Endogenous Circadian Period

Attaches a bootstrap confidence interval (and standard error) to the
Lomb-Scargle period estimate from
[`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md).
Because activity is strongly autocorrelated, an i.i.d. resample would
destroy the rhythm and give an invalid (too-narrow) interval, so a
*circular moving-block* bootstrap of the cosinor residuals is used: the
rhythm fitted at the point estimate is held fixed while blocks of
residuals are resampled, the periodogram peak is re-located on each
replicate (with sub-grid parabolic refinement), and the interval is
taken from the quantiles of the replicate periods.

## Usage

``` r
period.ci(
  counts,
  timestamps,
  from = 18,
  to = 30,
  ofac = 4,
  n_boot = 200,
  block_hours = 24,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- counts:

  Numeric activity vector; `NA` are dropped with their times.

- timestamps:

  POSIXct (or numeric-coercible) timestamps, same length.

- from, to:

  Period search window in hours (default 18, 30).

- ofac:

  Integer oversampling factor for the period grid (default 4).

- n_boot:

  Number of bootstrap replicates (default 200).

- block_hours:

  Moving-block length in hours (default 24).

- level:

  Confidence level (default 0.95).

- seed:

  Optional integer seed passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) for reproducible
  bootstrap draws.

## Value

An object of class `actiRhythm_period_ci`: a list with `tau` (the point
estimate, hours), `ci_lower`/`ci_upper`, `se`, the `level`, the number
of valid replicates, and `tau_boot` (the vector of bootstrap replicate
periods).

## References

Kunsch HR (1989). “The jackknife and the bootstrap for general
stationary observations.” *The Annals of Statistics*, **17**(3),
1217–1241.
[doi:10.1214/aos/1176347265](https://doi.org/10.1214/aos/1176347265) .

Politis DN, Romano JP (1992). “A circular block-resampling procedure for
stationary data.” In LePage R, Billard L (eds.), *Exploring the Limits
of Bootstrap*, 263–270. Wiley, New York.

## See also

[`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)

## Examples

``` r
# \donttest{
t_hours <- seq(0, 5 * 24 - 1/60, by = 1/60)
ts <- as.POSIXct("2024-01-01", tz = "UTC") + t_hours * 3600
counts <- 100 + 80 * cos(2 * pi * (t_hours - 8) / 24) + rnorm(length(t_hours), 0, 20)
period.ci(counts, ts, n_boot = 50, seed = 1)
#> Circadian Period with Bootstrap Confidence Interval
#> 
#>   tau:      24.088 h
#>   95% CI:   [23.975, 24.051] h
#>   SE:       0.020 h
#>   Method:   circular block residual bootstrap (50/50 valid reps)
# }
```
