# Per-Day Nonparametric Metrics

The intraday nonparametric metrics computed for each day separately,
showing within-recording drift that the pooled values hide. Interdaily
stability is a between-day measure and so is not a per-day quantity; the
per-day table reports L5, M10, their onset times, relative amplitude,
and intradaily variability (Goncalves et al. 2014).

## Usage

``` r
circadian.daily(counts, timestamps, min_hours = 12)
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- min_hours:

  Minimum hours of data a day needs to be reported (default 12).

## Value

An object of class `actiRhythm_daily`: a `daily` data frame, one row per
day. Never errors.

## References

Goncalves BSB, Cavalcanti PRA, Tavares GR, Campos TF, Araujo JF (2014).
Nonparametric methods in actigraphy: an update. *Sleep Science*,
7(3):158-164.
[doi:10.1016/j.slsci.2014.09.013](https://doi.org/10.1016/j.slsci.2014.09.013)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H"))
circadian.daily(ifelse(h >= 23 | h < 7, 5, 300), ts)
#> Per-Day Nonparametric Metrics (4 days)
#> 
#>        date L5 L5_onset_h M10 M10_onset_h    RA    IV  total
#>  2024-01-01  5          0 300           7 0.967 0.391 290400
#>  2024-01-02  5          0 300           7 0.967 0.391 290400
#>  2024-01-03  5          0 300           7 0.967 0.391 290400
#>  2024-01-04  5          0 300           7 0.967 0.391 290400
#> 
```
