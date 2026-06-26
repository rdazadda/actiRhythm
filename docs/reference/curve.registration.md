# Curve Registration of Daily Activity Profiles

Aligns each day's 24-hour activity profile on its active-phase landmark
(the M10 centre), separating the horizontal phase variation (how the
timing shifts day to day) from the vertical amplitude variation (the
registered mean profile, sharper than the plain average because phase
jitter no longer blurs it). The per-day landmark times are a
scale-invariant chronotype phase, unchanged by any rescaling of the
counts (Ramsay and Silverman 2005).

## Usage

``` r
curve.registration(counts, timestamps, n_grid = 144L, period = 24)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- n_grid:

  Bins per day for the profile (default 144 = 10-minute bins).

- period:

  Period in hours (default 24).

## Value

An object of class `actiRhythm_registration`: a per-day landmark table
(L5 and M10 centre hours), the circular-mean landmarks, the phase
variability (circular SD of the M10 landmark), and the registered mean
profile. Never errors.

## References

Ramsay JO, Silverman BW (2005). *Functional Data Analysis*, 2nd edition.
Springer. [doi:10.1007/b98888](https://doi.org/10.1007/b98888) .

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
h  <- as.numeric(format(ts, "%H"))
curve.registration(ifelse(h >= 23 | h < 7, 5, 300), ts)
#> Curve Registration (landmark)
#> 
#>   Days registered:   5
#>   Mean L5 centre:    02.50 h
#>   Mean M10 centre:   12.00 h
#>   Phase variability: -0.00 h (circular SD of M10)
#> 
```
