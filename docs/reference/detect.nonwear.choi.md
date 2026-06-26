# Choi (2011) Non-Wear Detection

Classifies each epoch as wear or non-wear with the Choi et al. (2011)
algorithm: a run of consecutive zero-count epochs of at least `frame`
minutes is non-wear, tolerating a short nonzero spike only when it is
flanked by a fully zero window of `stream` minutes both before and
after. The returned mask can be passed as the `wear_time` argument to
the rest-activity and sleep functions.

## Usage

``` r
detect.nonwear.choi(
  counts,
  epoch_length = 60,
  frame = 90,
  spike_tolerance = 2,
  stream = 30
)
```

## Arguments

- counts:

  Numeric activity vector (vertical axis), minute epochs assumed.

- epoch_length:

  Epoch length in seconds (default 60). Window lengths are given in
  minutes and scaled to epochs by this value.

- frame:

  Minimum non-wear window in minutes (default 90).

- spike_tolerance:

  Maximum tolerated nonzero spike in minutes (default 2).

- stream:

  Flanking all-zero window required around a tolerated spike, in minutes
  (default 30).

## Value

A logical vector, one per epoch: `TRUE` = wear, `FALSE` = non-wear.
Never errors.

## References

Choi L, Liu Z, Matthews CE, Buchowski MS (2011). “Validation of
accelerometer wear and nonwear time classification algorithm.” *Medicine
& Science in Sports & Exercise*, **43**(2), 357–364.
[doi:10.1249/MSS.0b013e3181ed61a3](https://doi.org/10.1249/MSS.0b013e3181ed61a3)
.

## See also

[`detect.nonwear.troiano`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.troiano.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 600)
counts <- c(rep(200, 200), rep(0, 200), rep(200, 200))   # 200-min non-wear gap
table(detect.nonwear.choi(counts))
#> 
#> FALSE  TRUE 
#>   200   400 
```
