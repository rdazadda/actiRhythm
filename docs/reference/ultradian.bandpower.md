# Ultradian Wavelet Band Power

Partitions the activity variance into dyadic period bands with an
undecimated (MODWT) Haar wavelet transform, isolating ultradian bands
such as the about-90-minute, about-4-hour, and about-8-hour rhythms
(Percival and Walden 2000; Leise 2013). Reports the energy, power, and
fraction of total variance in each band.

## Usage

``` r
ultradian.bandpower(
  counts,
  timestamps,
  bands = list(`90min` = c(1, 2), `4h` = c(2, 6), `8h` = c(6, 12)),
  epoch_length = 60
)
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- bands:

  Named list of `c(low_hours, high_hours)` period bands (default about
  90 min, about 4 h, about 8 h).

- epoch_length:

  Epoch length in seconds (default 60).

## Value

An object of class `actiRhythm_bandpower`: a per-band table (energy,
power, fraction) and the per-level wavelet variance. Never errors.

## References

Percival DB, Walden AT (2000). *Wavelet Methods for Time Series
Analysis*. Cambridge University Press.
[doi:10.1017/CBO9780511841040](https://doi.org/10.1017/CBO9780511841040)
.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
ultradian.bandpower(100 + 50 * sin(2 * pi * th / 4), ts)   # a 4-hour rhythm
#> Ultradian Wavelet Band Power
#> 
#>   band low_h high_h  power fraction
#>  90min     1      2  248.2    0.199
#>     4h     2      6 1173.8    0.939
#>     8h     6     12  435.8    0.349
#> 
```
