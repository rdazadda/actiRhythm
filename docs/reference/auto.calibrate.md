# Auto-Calibrate Raw Acceleration to the Unit Gravity Sphere

Estimates per-axis gain and offset corrections by the van Hees et al.
(2014) method: during non-movement windows the acceleration vector
should lie on the 1 g sphere, so gain/offset are fit by iteratively
projecting non-movement window means onto the closest point of the unit
sphere. Apply as `(raw - offset) * scale`.

## Usage

``` r
auto.calibrate(
  xyz,
  fs,
  sphere_crit = 0.3,
  sd_crit = 0.013,
  max_iter = 1000,
  tol = 1e-09
)
```

## Arguments

- xyz:

  A data frame or matrix of raw acceleration in g (columns x, y, z).

- fs:

  Sample rate in Hz.

- sphere_crit:

  Minimum coverage (g) each axis must span on both sides of zero for a
  stable fit (default 0.3).

- sd_crit:

  Per-axis rolling SD (g) below which a 10-second window counts as
  non-movement (default 0.013).

- max_iter:

  Maximum refinement iterations (default 1000).

- tol:

  Convergence tolerance on the calibration error (default 1e-9).

## Value

A list with `scale` and `offset` (length-3), the calibration error
before and after (`cal_error_start`, `cal_error_end`, mean absolute
deviation from 1 g), `npoints`, and a `calibrated` flag. When there is
too little non-movement data the identity correction is returned with
`calibrated = FALSE`.

## References

van Hees VT, et al. (2014). Autocalibration of accelerometer data for
free-living physical activity assessment using local gravity and
temperature. *Journal of Applied Physiology*, 117(7):738-744.
[doi:10.1152/japplphysiol.00421.2014](https://doi.org/10.1152/japplphysiol.00421.2014)

## Examples

``` r
# Recover a known per-axis gain and offset from non-movement windows
set.seed(1)
u <- matrix(rnorm(40 * 3), 40, 3); u <- u / sqrt(rowSums(u^2))   # sphere directions
raw <- do.call(rbind, lapply(seq_len(40), function(i)
  matrix(rep(u[i, ] / c(1.03, 0.97, 1.01) + c(0.04, -0.03, 0.02), each = 300),
         300, 3) + rnorm(900, 0, 0.004)))
auto.calibrate(data.frame(x = raw[, 1], y = raw[, 2], z = raw[, 3]), fs = 30)$scale
#> [1] 1.0301368 0.9701036 1.0099054
```
