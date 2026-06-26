# Consensus Rhythmicity Across Methods

Combines the rhythmicity verdicts of several independent tests into one
call: the cosinor zero-amplitude F-test
([`rhythmicity.test`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)),
the Bingham confidence ellipse
([`cosinor.confidence.ellipse`](https://rdazadda.github.io/actiRhythm/reference/cosinor.confidence.ellipse.md)),
the Lomb-Scargle Baluev false-alarm probability
([`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)),
and the chi-square (Sokolove-Bushell) periodogram
([`chi.sq.periodogram`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md)).
The available p-values are pooled by Fisher's method, and a majority
vote across all methods is also reported.

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

An object of class `actiRhythm_consensus`: the Fisher-combined p-value
and consensus call, the vote count, a per-method `tests` data frame, and
the agreement fraction.

## References

Fisher RA (1925). *Statistical Methods for Research Workers*. Oliver and
Boyd, Edinburgh.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h <- as.numeric(format(ts, "%H"))
counts <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24))
consensus.rhythmicity(counts, ts)
#> Consensus Rhythmicity (multi-method)
#> 
#>   cosinor F-test         p=<2e-16  rhythmic
#>   Bingham ellipse                  no
#>   Lomb-Scargle FAP       p=<2e-16  rhythmic
#>   chi-square periodogram p=<2e-16  rhythmic
#> 
#>   Votes:        3 / 4 methods
#>   Fisher p:     <2e-16
#>   Consensus:    RHYTHMIC (alpha = 0.05)
```
