# Multi-resolution Interdaily Stability (IS)

Computes interdaily stability at several within-day bin resolutions. IS
is classically reported at 60-min bins, but finer bins (30/15 min)
capture higher-frequency day-to-day regularity. For each bin width,
`IS = (N * sum_h (xbar_h - xbar)^2) / (p * sum_i (x_i - xbar)^2)` where
`p` is bins-per-day, `xbar_h` the mean for bin-of-day `h`, and `N` the
total number of bins.

## Usage

``` r
circadian.is.multiscale(
  counts,
  timestamps,
  bin_minutes = (1:60)[1440L%%(1:60) == 0L]
)
```

## Arguments

- counts:

  Numeric vector of epoch-level activity.

- timestamps:

  POSIXct timestamps (one per epoch).

- bin_minutes:

  Integer bin widths in minutes that divide 1440 (default the divisors
  of 1440 from 1 to 60 min, per Goncalves et al. 2014).

## Value

An object of class `actiRhythm_ism`: a per-bin `IS` table and the
averaged `ISm`.

## References

Witting W, Kwa IH, Eikelenboom P, Mirmiran M, Swaab DF (1990).
“Alterations in the circadian rest-activity rhythm in aging and
Alzheimer's disease.” *Biological Psychiatry*, **27**(6), 563–572.
[doi:10.1016/0006-3223(90)90523-5](https://doi.org/10.1016/0006-3223%2890%2990523-5)
.

Goncalves BSB, Cavalcanti PRA, Tavares GR, Campos TF, Araujo JF (2014).
“Nonparametric methods in actigraphy: an update.” *Sleep Science*,
**7**(3), 158–164.
[doi:10.1016/j.slsci.2014.09.013](https://doi.org/10.1016/j.slsci.2014.09.013)
.
