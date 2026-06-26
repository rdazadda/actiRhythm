# Non-Wear Detection from Raw Acceleration

Flags non-wear time from raw acceleration by the van Hees et al. (2011)
standard-deviation-and-range rule: a block is non-wear when, over a
window centred on it (default 60 minutes), at least two of the three
axes have both a standard deviation below `sd_crit` and a value range
below `range_crit`. A stationary, taken-off device reads as non-wear.
This lets the z-angle sleep detector tell device-off periods from real
sleep (a still arm that keeps micro-movement). Pass the result as the
`wear` argument of
[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md).

## Usage

``` r
detect.nonwear.raw(
  x,
  device = "auto",
  epoch = 5,
  block = 300,
  window = 3600,
  sd_crit = 0.013,
  range_crit = 0.05,
  tz = "UTC"
)
```

## Arguments

- x:

  A path to a raw file or a raw data frame (see
  [`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)).

- device:

  Device brand or `"auto"` (file input only).

- epoch:

  Output epoch length in seconds for the returned mask (default 5,
  matching
  [`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)).

- block:

  Internal classification block length in seconds (default 300).

- window:

  Window in seconds over which the SD and range are taken (default
  3600).

- sd_crit:

  Per-axis SD threshold in g (default 0.013).

- range_crit:

  Per-axis range threshold in g (default 0.050).

- tz:

  Time zone (default `"UTC"`).

## Value

A logical vector, one per epoch: `TRUE` = wear, `FALSE` = non-wear.
Never errors.

## References

van Hees VT, Renstrom F, Wright A, Gradmark A, Catt M, Chen KY, Lof M,
Bluck L, Pomeroy J, Wareham NJ, Ekelund U, Brage S, Franks PW (2011).
“Estimation of daily energy expenditure in pregnant and non-pregnant
women using a wrist-worn tri-axial accelerometer.” *PLoS ONE*, **6**(7),
e22922.
[doi:10.1371/journal.pone.0022922](https://doi.org/10.1371/journal.pone.0022922)
.

## See also

[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md),
[`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)

## Examples

``` r
# \donttest{
raw  <- example_raw(days = 2, device_off = 1)   # two worn days + one device-off day
mean(detect.nonwear.raw(raw, epoch = 60))       # fraction of epochs worn
#> [1] 0.6736111
# }
```
