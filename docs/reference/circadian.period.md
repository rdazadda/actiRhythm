# Estimate Endogenous Circadian Period via the Lomb-Scargle Periodogram

Computes the Lomb-Scargle periodogram of an activity `counts` series
sampled at the supplied `timestamps` and returns the period (tau, in
hours) of the strongest spectral peak within the search window
`[from, to]`. The Lomb-Scargle method (Lomb 1976; Scargle 1982) is the
least-squares-equivalent spectral estimator for unevenly sampled time
series and is therefore appropriate for actigraphy data containing gaps,
which the FFT cannot accommodate.

## Usage

``` r
circadian.period(counts, timestamps, from = 18, to = 30, ofac = 4)
```

## Arguments

- counts:

  Numeric vector of activity counts (minute-level recommended). `NA`
  values (e.g. non-wear epochs) are dropped together with their
  timestamps before estimation.

- timestamps:

  A `POSIXct` vector (or anything coercible by `as.numeric`) of epoch
  timestamps, the same length as `counts`. Internally converted to hours
  elapsed since the first timestamp.

- from:

  Numeric. Lower bound of the period search window, in hours (default
  `18`).

- to:

  Numeric. Upper bound of the period search window, in hours (default
  `30`).

- ofac:

  Integer oversampling factor controlling the period-grid resolution.
  Higher values give a finer period grid and a more precise peak
  location at the cost of computation (default `4`).

## Value

A named `list` with elements:

- tau:

  Numeric. Period (hours) of the strongest Lomb-Scargle peak in
  `[from, to]`, i.e. the estimated endogenous circadian period.
  `NA_real_` when the data are insufficient.

- peak_power:

  Numeric. Normalized power of that peak (the Lomb-Scargle peak
  statistic). `NA_real_` when insufficient.

- p_value:

  Numeric. P-value of the peak under the null hypothesis of Gaussian
  noise (Baluev 2008 analytic false-alarm probability). `NA_real_` when
  insufficient.

- oversampling:

  The `ofac` oversampling factor used.

- n_used:

  Integer. Number of non-`NA` observations actually passed to the
  periodogram (`NA_integer_` when not run).

- span_days:

  Numeric. Total recording span in days (max minus min timestamp), used
  for the \>= 2-day guard.

- scanned:

  Numeric vector of trial periods (hours) of the full Lomb-Scargle
  spectrum (`numeric(0)` when not run).

- power:

  Numeric vector of standard-normalized Lomb-Scargle power, aligned to
  `scanned` (`numeric(0)` when not run).

On any edge case (too few points, too short a span, degenerate input, or
an internal numerical failure) the function returns this same structure
with `tau`, `peak_power` and `p_value` set to `NA`; it never throws.

## Details

Processing steps:

1.  Timestamps are converted to hours since the first sample
    (`t_hours = (as.numeric(timestamps) - min) / 3600`).

2.  Pairs with a missing `count` or a missing/non-finite time are
    dropped.

3.  Two guards are applied so the estimate is never based on
    insufficient data: the recording must span at least **2 days**
    (otherwise a 18-30 h period cannot be resolved) and at least **10
    non-`NA`** observations must remain.

4.  The standard-normalized Lomb-Scargle periodogram is evaluated over
    the period window. Its strongest peak gives the period (hours), and
    the Baluev (2008) analytic false-alarm probability gives the
    p-value.

The Lomb-Scargle periodogram is chosen specifically because actigraphy
series are rarely gap-free: the FFT assumes uniform sampling, whereas
Lomb-Scargle fits sinusoids by least squares at each trial frequency and
is unbiased for irregular sampling.

## References

Lomb NR (1976). Least-squares frequency analysis of unequally spaced
data. *Astrophysics and Space Science*, 39(2):447-462.

Scargle JD (1982). Studies in astronomical time series analysis. II.
Statistical aspects of spectral analysis of unevenly spaced data. *The
Astrophysical Journal*, 263:835-853.

Ruf T (1999). The Lomb-Scargle periodogram in biological rhythm
research: analysis of incomplete and unequally spaced time-series.
*Biological Rhythm Research*, 30(2):178-201.

Refinetti R, Cornelissen G, Halberg F (2007). Procedures for numerical
analysis of circadian rhythms. *Biological Rhythm Research*,
38(4):275-325.

## See also

[`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
for parametric (fixed-period) rhythm estimation,
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
for non-parametric L5/M10/IS/IV metrics.

## Examples

``` r
# \donttest{
# Seven days of minute-level data with a 24 h rhythm
t_hours <- seq(0, 7 * 24 - 1/60, by = 1/60)
ts <- as.POSIXct("2024-01-01 00:00:00") + t_hours * 3600
counts <- 100 + 80 * cos(2 * pi * (t_hours - 8) / 24) + rnorm(length(t_hours), 0, 5)
circadian.period(counts, ts)$tau   # about 24.0
#> [1] 23.99762
# }
```
