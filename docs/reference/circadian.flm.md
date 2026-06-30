# Functional Linear Model of the 24-Hour Activity Profile

Fits the averaged 24-hour activity profile with a periodic basis
expansion (Fourier by default, B-spline alternative) by weighted least
squares, giving a smooth functional form of the daily activity pattern.
The single-component cosinor is the one-harmonic special case; adding
harmonics fits the non-sinusoidal shape of a real rest-activity profile.
This follows the functional form of Wang et al. (2011), as implemented
for actigraphy by pyActigraphy.

## Usage

``` r
circadian.flm(
  counts,
  timestamps,
  basis = c("fourier", "bspline"),
  n_harmonics = 4,
  nbasis = 9,
  spline_order = 4,
  period = 24,
  wear_time = NULL,
  min_valid_hours = 10,
  weights = c("n", "none"),
  n_eval = 1440
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- basis:

  Basis type: `"fourier"` (default) or `"bspline"`.

- n_harmonics:

  Number of Fourier harmonics (default 4, giving a 9-term expansion as
  in Wang et al. (2011)).

- nbasis:

  Number of B-spline basis functions (default 9).

- spline_order:

  B-spline order (default 4, cubic).

- period:

  Profile period in hours (default 24).

- wear_time:

  Optional logical wear-time mask.

- min_valid_hours:

  Minimum valid hours per day (default 10).

- weights:

  Profile weighting: `"n"` (square root of the per-hour observation
  count, default, matching
  [`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md))
  or `"none"` (plain least squares, matching pyActigraphy).

- n_eval:

  Length of the dense within-day evaluation grid (default 1440).

## Value

An object of class `actiRhythm_flm`: the fitted `coefficients`, the
smooth daily curve (`smooth_curve`), the fitted profile, per-harmonic
amplitudes and acrophases (`harmonics`, Fourier only), the peak and
trough, and the fit statistics (`r_squared`, `aic`, `f_statistic`,
`p_value`). The function never errors; on insufficient data it returns
the same structure with `r_squared` NA.

## References

Wang J, Xian H, Licis A, Deych E, Ding J, McLeland J, Toedebusch C, Li
T, Duntley S, Shannon W (2011). “Measuring the impact of apnea and
obesity on circadian activity patterns using functional linear modeling
of actigraphy data.” *Journal of Circadian Rhythms*, **9**, 11.
[doi:10.1186/1740-3391-9-11](https://doi.org/10.1186/1740-3391-9-11) .

Ramsay JO, Silverman BW (2005). *Functional Data Analysis*, 2nd edition.
Springer. [doi:10.1007/b98888](https://doi.org/10.1007/b98888) .

Hammad G, Reyt M, Beliy N, Baillet M, Deantoni M, Lesoinne A, Muto V,
Schmidt C (2021). “pyActigraphy: open-source python package for
actigraphy data visualization and analysis.” *PLOS Computational
Biology*, **17**(10), e1009514.
[doi:10.1371/journal.pcbi.1009514](https://doi.org/10.1371/journal.pcbi.1009514)
.

## See also

[`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
[`circadian.ssa`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
counts <- 100 + 60 * cos(2 * pi * (h - 14) / 24) + 25 * cos(2 * pi * (h - 6) / 12)
circadian.flm(counts, ts)
#> Functional Linear Model (24-hour activity profile)
#> 
#>   Basis:        fourier (order 4, 9 terms)
#>   Period:       24 hours
#>   Days / profile points: 4 / 24
#> 
#>   Model fit:
#>     R-squared:   1.0000
#>     AIC:         -1371.37
#>     F-statistic: 13059180580519173999266600266066.00 (p = <2e-16)
#> 
#>   Dominant harmonic: H1, amplitude 59.83, acrophase 14.01 h
#> 
#>   Peak 165.0 at 16.5 h; trough 20.5 at 0.8 h
#> 
#>   Basis fit follows Wang et al. (2011); F/p are ordinary regression
#>   statistics (not Wang's between-subject permutation test).
#> 
```
