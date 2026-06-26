# Singular Spectrum Analysis of an Activity Series

Decomposes an activity-count series into additive components (trend, a
circadian component, ultradian components, and noise) with Basic
Singular Spectrum Analysis: embed the series into a Hankel trajectory
matrix, take its singular value decomposition, group the resulting
elementary series, and reconstruct each group by diagonal averaging.
pyActigraphy implements this model-free decomposition for actigraphy.

## Usage

``` r
circadian.ssa(
  counts,
  timestamps,
  window_hours = 24,
  n_components = 10,
  groups = NULL,
  w_components = NULL,
  period_range = c(20, 28),
  detrend = FALSE
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- window_hours:

  Embedding window length in hours (default 24). The window in epochs is
  `L = round(window_hours * 3600 / epoch_seconds)`.

- n_components:

  Number of leading elementary components to reconstruct and keep
  (default 10).

- groups:

  Optional named list of 1-based component indices, e.g.
  `list(trend = 1, circadian = 2:3)`. When `NULL` the grouping is chosen
  automatically (component 1 is the trend; the circadian pair is the
  largest component whose period falls in `period_range`).

- w_components:

  Number of leading components used for the w-correlation matrix
  (default `min(n_components, 10)`).

- period_range:

  Two-element period window in hours used to identify the circadian
  component (default `c(20, 28)`).

- detrend:

  If `TRUE`, the trend component is removed from `reconstructed`.

## Value

An object of class `actiRhythm_ssa`: the singular values and partial
variances, the reconstructed component series (`trend`, `circadian`,
`ultradian`), the w-correlation matrix, the circadian
`fundamental_period`, and the share of variance the circadian component
carries. The function never errors; on insufficient data it returns the
same structure with `insufficient = TRUE`.

## Details

Singular Spectrum Analysis on a long minute-level series builds a large
trajectory matrix; for multi-day recordings, resampling to a coarser
epoch (for example 10 to 30 minutes) before calling keeps the
decomposition fast.

## References

Golyandina N, Zhigljavsky A (2013). “Singular Spectrum Analysis for Time
Series.” *SpringerBriefs in Statistics*.
[doi:10.1007/978-3-642-34913-3](https://doi.org/10.1007/978-3-642-34913-3)
.

Vautard R, Yiou P, Ghil M (1992). “Singular-spectrum analysis: a toolkit
for short, noisy chaotic signals.” *Physica D: Nonlinear Phenomena*,
**58**(1-4), 95–126.
[doi:10.1016/0167-2789(92)90103-T](https://doi.org/10.1016/0167-2789%2892%2990103-T)
.

Hammad G, Reyt M, Beliy N, Baillet M, Deantoni M, Lesoinne A, Muto V,
Schmidt C (2021). “pyActigraphy: open-source python package for
actigraphy data visualization and analysis.” *PLOS Computational
Biology*, **17**(10), e1009514.
[doi:10.1371/journal.pcbi.1009514](https://doi.org/10.1371/journal.pcbi.1009514)
.

## See also

[`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md),
[`circadian.flm`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 1800, length.out = 7 * 48)
h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
counts <- 50 + 0.5 * seq_along(ts) + 60 * cos(2 * pi * (h - 14) / 24)
ssa <- circadian.ssa(counts, ts)
ssa
#> Singular Spectrum Analysis (Basic SSA)
#> 
#>   Window length L:    48 epochs (24.0 h)
#>   Series length n:    336 (K = 289, span 7.0 days)
#>   Components kept:     4 of rank 4
#> 
#>   Variance explained (leading components):
#>     ET1   lambda = 0.9171  (cumulative 0.9171)
#>     ET2   lambda = 0.0418  (cumulative 0.9589)
#>     ET3   lambda = 0.0410  (cumulative 0.9999)
#>     ET4   lambda = 0.0001  (cumulative 1.0000)
#>     ET5   lambda = 0.0000  (cumulative 1.0000)
#> 
#>   Grouping:
#>     Trend:      components 1
#>     Circadian:  components 2, 3 (8.3% of variance)
#>     Fundamental period: 23.93 h
#> 
#>   Reference: Golyandina and Zhigljavsky (2013)
#> 
```
