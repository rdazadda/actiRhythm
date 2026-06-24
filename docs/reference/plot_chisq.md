# Plot the Chi-Square (Sokolove-Bushell) Periodogram

Draws the chi-square periodogram from
[`chi.sq.periodogram`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md):
the \\Q_P\\ statistic against trial period, with the per-period
chi-square significance threshold overlaid as a dashed line. The
estimated period (the \\Q_P\\ peak) is marked and a 24-hour reference
line drawn. Unlike the Lomb-Scargle plot, the explicit threshold shows
directly whether a rhythm is significant at a given period.

## Usage

``` r
plot_chisq(
  counts,
  timestamps,
  from = 18,
  to = 30,
  alpha = 0.05,
  epoch_length = NULL
)
```

## Arguments

- counts:

  Numeric vector of activity counts on a regular epoch grid.

- timestamps:

  A `POSIXct` vector the same length as `counts`.

- from, to:

  Period search window in hours (default 18 to 30).

- alpha:

  Significance level for the chi-square threshold (default 0.05).

- epoch_length:

  Epoch length in seconds; inferred from `timestamps` when `NULL`.

## Value

A `ggplot` object. Never errors; returns an annotated empty plot on
insufficient data.

## See also

[`plot_periodogram`](https://rdazadda.github.io/actiRhythm/reference/plot_periodogram.md)
for the Lomb-Scargle spectrum.
