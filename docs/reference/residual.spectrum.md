# Residual Circadian Spectrum

Removes the fitted cosinor mean and estimates the spectrum of what is
left, the residual circadian spectrum, integrating it into frequency
bands. It measures the ultradian and noise structure that the cosinor
does not fit (Krafty et al. 2019), so two recordings with the same
24-hour rhythm but different residual fragmentation are told apart.

## Usage

``` r
residual.spectrum(
  counts,
  timestamps,
  period = 24,
  bands = list(ultradian = c(2, 8), high_freq = c(0.5, 2))
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- period:

  Cosinor period in hours (default 24).

- bands:

  Named list of `c(low_hours, high_hours)` period bands for the residual
  power; defaults to ultradian (2-8 h) and high-frequency (\<2 h).

## Value

An object of class `actiRhythm_rcs`: the residual variance, a per-band
power table, and the spectrum. Never errors.

## References

Krafty RT, Fu H, Graves JL, et al. (2019). Measuring variability in
rest-activity rhythms from actigraphy with application to characterizing
symptoms of depression. *Statistics in Biosciences*, 11:314-333.
[doi:10.1007/s12561-018-09230-2](https://doi.org/10.1007/s12561-018-09230-2)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
residual.spectrum(100 + 80 * cos(2 * pi * th / 24) + 20 * sin(2 * pi * th / 4), ts,
                  period = 24)
#> Residual Circadian Spectrum
#> 
#>   Residual variance: 200.2
#> 
#>       band low_h high_h power fraction
#>  ultradian   2.0      8 86400        1
#>  high_freq   0.5      2     0        0
#> 
```
