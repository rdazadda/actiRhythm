# Multi-resolution Interdaily Stability (IS)

Computes interdaily stability at several within-day bin resolutions. IS
is classically reported at 60-min bins, but finer bins (30/15 min)
capture higher-frequency day-to-day regularity. For each bin width,
`IS = (N * sum_h (xbar_h - xbar)^2) / (p * sum_i (x_i - xbar)^2)` where
`p` is bins-per-day, `xbar_h` the mean for bin-of-day `h`, and `N` the
total number of bins.

## Usage

``` r
circadian.is.multiscale(counts, timestamps, bin_minutes = c(60, 30, 15))
```

## Arguments

- counts:

  Numeric vector of epoch-level activity.

- timestamps:

  POSIXct timestamps (one per epoch).

- bin_minutes:

  Integer vector of bin widths in minutes (default `c(60, 30, 15)`).

## Value

Data frame with columns `bin_minutes` and `IS`.

## References

Witting W, et al. (1990). Biol Psychiatry 27(6):563-572.
