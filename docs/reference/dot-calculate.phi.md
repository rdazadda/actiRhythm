# Calculate Phi (First-Order Autocorrelation)

Phi is the correlation of the activity time series with itself at a
1-hour lag. GGIR uses it as a measure of activity predictability.

## Usage

``` r
.calculate.phi(counts, timestamps)
```

## Arguments

- counts:

  Numeric vector of activity counts

- timestamps:

  POSIXct timestamps

## Value

Numeric phi value (-1 to 1, higher = more autocorrelated)

## Details

Higher phi means more consistent, predictable activity. Low or negative
phi means fragmented, unpredictable activity.

This is the lag-1 autocorrelation of the hourly activity series. The
series is placed on a regular hourly grid (NA for missing/non-wear
hours) so the "1-hour lag" only correlates genuinely adjacent clock
hours rather than bridging multi-hour gaps. It equals the AR(1)
coefficient only for a true AR(1) process. It is computed directly with
[`stats::acf`](https://rdrr.io/r/stats/acf.html) and has no separate
method reference.
