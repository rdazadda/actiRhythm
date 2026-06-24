# Generalized Least- and Most-Active Periods (MX / LX)

The least-active `LX` and most-active `MX` periods for arbitrary window
lengths, generalizing L5 and M10 (Van Someren et al. 1999). For each
window the function returns the mean activity over the window and its
onset and midpoint clock times. A fixed L5/M10 pair does not give these
phase markers.

## Usage

``` r
activity.extrema(counts, timestamps, windows = c(5, 10))
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- windows:

  Window lengths in hours (default `c(5, 10)` = L5, M10).

## Value

An object of class `actiRhythm_mxlx`: a data frame with, for each
window, the least- and most-active mean level, onset hour, and midpoint
hour. Never errors.

## References

Van Someren EJW, Swaab DF, Colenda CC, Cohen W, McCall WV, Rosenquist PB
(1999). Bright light therapy: improved sensitivity to its effects on
rest-activity rhythms in Alzheimer patients by application of
nonparametric methods. *Chronobiology International*, 16(4):505-518.
[doi:10.3109/07420529908998724](https://doi.org/10.3109/07420529908998724)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H"))
activity.extrema(ifelse(h >= 23 | h < 7, 5, 300), ts, windows = c(5, 8, 10))
#> Least / Most Active Periods (MX / LX)
#> 
#>   L5: 5 at 00.0h   M5: 300 at 07.0h
#>   L8: 5 at 23.0h   M8: 300 at 07.0h
#>   L10: 64 at 21.0h   M10: 300 at 07.0h
#> 
```
