# Chi-square (Sokolove-Bushell) Periodogram

Estimates the dominant period of an activity time series with the
Sokolove-Bushell (1978) chi-square periodogram, the periodogram most
widely used in chronobiology and actigraphy. It pairs with the
Lomb-Scargle estimator in
[`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md):
Lomb-Scargle is the least-squares spectral method for unevenly sampled
data; the chi-square periodogram is the analysis-of-variance method for
regularly sampled data and reports a significance threshold.

## Usage

``` r
chi.sq.periodogram(
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

  Numeric vector of activity counts on a regular epoch grid. `NA` values
  (e.g. non-wear) are handled by removal from the phase and grand means.

- timestamps:

  A `POSIXct` vector (or numeric seconds) the same length as `counts`;
  used to infer the epoch length.

- from, to:

  Period search window in hours (default `18` to `30`).

- alpha:

  Significance level for the chi-square threshold (default `0.05`).

- epoch_length:

  Epoch length in seconds. If `NULL` (default) it is inferred from the
  median spacing of `timestamps`.

## Value

A named list with `period` (hours, the \\Q_P\\ peak), `Qp_peak`,
`p_value` (family-wise Sidak p-value of the peak across the scanned
periods), `significant` (logical, `p_value < alpha`), `scanned` (trial
periods in hours), `Qp` (the periodogram aligned to `scanned`),
`critical` (the per-period chi-square threshold), `epoch_length` and
`alpha`. On insufficient/degenerate data the same shape is returned with
`period`/`Qp_peak` `NA` and empty vectors; the function never throws.

## Details

For each trial period \\P\\ (an integer number of epochs) the series is
folded into \\P\\ phase bins over \\K = \lfloor N/P \rfloor\\ complete
cycles and the statistic \$\$Q_P = \frac{K\\N
\sum\_{h=1}^{P}(\bar{A}\_h - \bar{A})^2}{\sum\_{i=1}^{N}(A_i -
\bar{A})^2}\$\$ is computed, where \\\bar{A}\_h\\ is the mean at phase
\\h\\ and \\\bar{A}\\ the grand mean over the \\N = KP\\ retained
points. Under the null hypothesis of no rhythm at \\P\\, \\Q_P\\ follows
a chi-square distribution with \\P-1\\ degrees of freedom, giving the
significance line \\\chi^2\_{P-1,\\1-\alpha}\\. The estimated period is
the \\P\\ maximising \\Q_P\\ within the search window.

## References

Sokolove PG, Bushell WN (1978). “The chi square periodogram: its utility
for analysis of circadian rhythms.” *Journal of Theoretical Biology*,
**72**(1), 131–160.
[doi:10.1016/0022-5193(78)90022-X](https://doi.org/10.1016/0022-5193%2878%2990022-X)
.

Refinetti R, Cornelissen G, Halberg F (2007). “Procedures for numerical
analysis of circadian rhythms.” *Biological Rhythm Research*, **38**(4),
275–325.
[doi:10.1080/09291010600903692](https://doi.org/10.1080/09291010600903692)
.

Sidak Z (1967). “Rectangular confidence regions for the means of
multivariate normal distributions.” *Journal of the American Statistical
Association*, **62**(318), 626–633.
[doi:10.1080/01621459.1967.10482935](https://doi.org/10.1080/01621459.1967.10482935)
.

## See also

[`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
for the Lomb-Scargle estimator.
