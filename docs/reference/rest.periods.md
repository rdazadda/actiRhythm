# Consolidated Rest-Period Detection (Roenneberg / MASDA)

Detects all consolidated rest bouts across a recording with the
Roenneberg consolidation algorithm (the Munich Actimetry Sleep Detection
Algorithm). Each epoch is compared to a fraction of its own 24-hour
activity trend, runs below that threshold seed candidate rest bouts, and
a correlation procedure grows each seed into a consolidated bout. Unlike
[`sleep.changepoints`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
(one nightly bout per cycle), this returns any number of bouts,
including daytime naps and fragmented or polyphasic rest.

## Usage

``` r
rest.periods(
  counts,
  timestamps,
  epoch_length = 60,
  trend_period = 86400,
  min_trend_period = 43200,
  threshold = 0.15,
  min_seed_period = 1800,
  max_test_period = 43200,
  r_consec_below = 1800,
  nap_max_minutes = 180
)
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- epoch_length:

  Epoch length in seconds (default 60).

- trend_period:

  Moving-average window for the activity trend, in seconds (default
  86400, 24 hours).

- min_trend_period:

  Minimum non-missing seconds required in the trend window (default
  43200, 12 hours).

- threshold:

  Fraction of the local trend below which an epoch is candidate rest
  (default 0.15).

- min_seed_period:

  Minimum below-threshold run to seed a bout, in seconds (default 1800,
  30 minutes).

- max_test_period:

  Maximum bout length the consolidation tests, in seconds (default
  43200, 12 hours).

- r_consec_below:

  Above-threshold tolerance during consolidation, in seconds (default
  1800, 30 minutes).

- nap_max_minutes:

  A non-main bout shorter than this is labelled a nap (default 180).

## Value

An object of class `actiRhythm_roenneberg`: a `rest_periods` data frame
(one row per consolidated bout, with onset, offset, duration, and a
main/nap label), the per-epoch `rest` vector, the activity `trend`, and
per-bout counts. The function never errors; with no resolvable bout it
returns an empty `rest_periods`.

## References

Roenneberg T, Keller LK, Fischer D, Matera JL, Vetter C, Winnebeck EC
(2015). Human activity and rest in situ. *Methods in Enzymology*,
552:257-283.
[doi:10.1016/bs.mie.2014.11.028](https://doi.org/10.1016/bs.mie.2014.11.028)

Loock A-S, Khan Sullivan A, Reis C, et al. (2021). Validation of the
Munich Actimetry Sleep Detection Algorithm for estimating sleep-wake
patterns from activity recordings. *Journal of Sleep Research*,
30(6):e13371. [doi:10.1111/jsr.13371](https://doi.org/10.1111/jsr.13371)

Hammad G, Reyt M, Beliy N, et al. (2021). pyActigraphy: open-source
python package for actigraphy data visualization and analysis. *PLOS
Computational Biology*, 17(10):e1009514.
[doi:10.1371/journal.pcbi.1009514](https://doi.org/10.1371/journal.pcbi.1009514)

## See also

[`sleep.changepoints`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md),
[`sleep.cole.kripke`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)

## Examples

``` r
# Three nights of rest plus one daytime nap
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
day <- as.integer(format(ts, "%j"))
counts <- ifelse(h >= 23 | h < 7, 5, 300)
counts[h >= 14 & h < 15.5 & day == min(day) + 1L] <- 5   # a nap on day 2
rest.periods(counts, ts)
#> Consolidated Rest Periods (Roenneberg/MASDA)
#> 
#>   Bouts detected:   4  (2.0 per day over 2 days)
#>   Total rest:       24.5 h
#>   Main bouts / naps: 3 / 1
#> 
#>   First bouts:
#>     main 01-01 00:00 -> 01-01 06:59  (420 min)
#>     main 01-01 23:00 -> 01-02 06:59  (480 min)
#>     nap  01-02 14:00 -> 01-02 15:29  (90 min)
#> 
#>   Reference: Roenneberg et al. (2015); Loock et al. (2021)
#> 
```
