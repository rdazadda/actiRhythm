# Crespo Rest/Activity Period Detection

Detects the main rest and activity periods across a recording with the
Crespo algorithm (Crespo et al. 2012): a rank-order (median) filter at
the `alpha` (~rest-length) scale, an `alpha`/24h percentile threshold,
and a separate minutes-scale morphology pass that turn activity counts
into a binary rest/activity series read bout by bout from its
transitions. The wide smoothing window suppresses short within-period
transitions, so it reports the main daily rest rather than every nap.

## Usage

``` r
rest.crespo(
  counts,
  timestamps,
  epoch_length = 60,
  zeta = 15,
  zeta_r = 30,
  zeta_a = 2,
  t = 0.33,
  alpha = 8 * 3600,
  beta = 1 * 3600
)
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- epoch_length:

  Epoch length in seconds (default 60).

- zeta, zeta_r, zeta_a:

  Maximum valid consecutive-zero run, in epochs, for the global
  pre-conditioning, within rest segments, and within active segments
  (defaults 15, 30, 2). Longer zero runs are treated as non-wear.

- t:

  Quantile of the counts used to replace non-wear epochs (default 0.33).

- alpha:

  Expected daily rest length, in seconds (default 28800, 8 hours); sets
  the rank-order smoothing window (Eq 4), the threshold percentile
  (alpha / 24 h), and the minimum data required.

- beta:

  Morphology scale, in seconds (default 3600, 1 hour); sets the
  structuring-element size that consolidates the thresholded series.

## Value

An object of class `actiRhythm_crespo`: a `rest_periods` data frame (one
row per bout, with onset, offset, and duration), the per-epoch
`rest_state` ("R"/"A"), and per-bout counts. The function never errors;
with no resolvable bout it returns an empty `rest_periods`.

## References

Crespo C, Aboy M, Fernandez JR, Mojon A (2012). “Automatic
identification of activity-rest periods based on actigraphy.” *Medical &
Biological Engineering & Computing*, **50**(4), 329–340.
[doi:10.1007/s11517-012-0875-y](https://doi.org/10.1007/s11517-012-0875-y)
.

Hammad G, Reyt M, Beliy N, Baillet M, Deantoni M, Lesoinne A, Muto V,
Schmidt C (2021). “pyActigraphy: open-source python package for
actigraphy data visualization and analysis.” *PLOS Computational
Biology*, **17**(10), e1009514.
[doi:10.1371/journal.pcbi.1009514](https://doi.org/10.1371/journal.pcbi.1009514)
.

## See also

[`rest.periods`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md),
[`sleep.changepoints`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)

## Examples

``` r
# Two nights of rest with a daytime nap between them
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
day <- as.integer(format(ts, "%j"))
counts <- ifelse(h >= 23 | h < 7, 5, 300)
counts[h >= 14 & h < 15 & day == min(day) + 1L] <- 5
rest.crespo(counts, ts)
#> Crespo Rest/Activity Periods
#> 
#>   Rest bouts:   3  (1.5 per day over 2 days)
#>   Total rest:   23.0 h
#> 
#>   First bouts:
#>     01-01 00:01 -> 01-01 06:59  (419 min)
#>     01-01 23:00 -> 01-02 06:59  (480 min)
#>     01-02 23:00 -> 01-03 06:59  (480 min)
#> 
#>   Reference: Crespo et al. (2012)
#> 
```
