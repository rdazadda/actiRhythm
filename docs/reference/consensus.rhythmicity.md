# Consensus Rhythmicity Across Methods

Combines the rhythmicity verdicts of three tests into one call: the
cosinor zero-amplitude F-test
([`rhythmicity.test`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)),
the Lomb-Scargle Baluev false-alarm probability
([`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)),
and the chi-square (Sokolove-Bushell) periodogram
([`chi.sq.periodogram`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md)).
The three share one series, so their p-values are pooled by the Cauchy
combination, which is valid under arbitrary dependence, and a majority
vote is also reported.

## Usage

``` r
consensus.rhythmicity(
  counts,
  timestamps,
  period = 24,
  alpha = 0.05,
  wear_time = NULL
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- period:

  Rhythm period in hours for the cosinor tests (default 24).

- alpha:

  Significance level (default 0.05).

- wear_time:

  Optional logical wear-time mask for the cosinor fit.

## Value

An object of class `actiRhythm_consensus`: the combined p-value and
consensus call, the vote count, a per-method `tests` data frame, and the
agreement fraction.

## References

Liu Y, Xie J (2020). “Cauchy combination test: a powerful test with
analytic p-value calculation under arbitrary dependency structures.”
*Journal of the American Statistical Association*, **115**(529),
393–402.
[doi:10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485)
.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h <- as.numeric(format(ts, "%H"))
counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
consensus.rhythmicity(counts, ts)
#> Consensus Rhythmicity (multi-method)
#> 
#>   cosinor F-test         p=<2e-16  rhythmic
#>   Lomb-Scargle FAP       p=<2e-16  rhythmic
#>   chi-square periodogram p=<2e-16  rhythmic
#> 
#>   Votes:        3 / 3 methods
#>   Combined p:   <2e-16
#>   Consensus:    RHYTHMIC (alpha = 0.05)
```
