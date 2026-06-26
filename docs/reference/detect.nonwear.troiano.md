# Troiano (2008) Non-Wear Detection

Classifies each epoch as wear or non-wear with the Troiano et al. (2008,
NHANES) algorithm: a run of at least `frame` minutes of zero counts is
non-wear, tolerating up to `spike_tolerance` consecutive nonzero minutes
only if every one of them is at or below `stoplevel` counts. Unlike Choi
it has no flanking-window requirement but applies a count ceiling on
spikes.

## Usage

``` r
detect.nonwear.troiano(
  counts,
  epoch_length = 60,
  frame = 60,
  spike_tolerance = 2,
  stoplevel = 100
)
```

## Arguments

- counts:

  Numeric activity vector (vertical axis), minute epochs assumed.

- epoch_length:

  Epoch length in seconds (default 60).

- frame:

  Minimum non-wear window in minutes (default 60).

- spike_tolerance:

  Maximum tolerated nonzero spike in minutes (default 2).

- stoplevel:

  Count above which a spike ends the non-wear bout (default 100).

## Value

A logical vector, one per epoch: `TRUE` = wear, `FALSE` = non-wear.
Never errors.

## References

Troiano RP, Berrigan D, Dodd KW, Masse LC, Tilert T, McDowell M (2008).
“Physical activity in the United States measured by accelerometer.”
*Medicine & Science in Sports & Exercise*, **40**(1), 181–188.
[doi:10.1249/mss.0b013e31815a51b3](https://doi.org/10.1249/mss.0b013e31815a51b3)
.

## See also

[`detect.nonwear.choi`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.choi.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 400)
counts <- c(rep(150, 100), rep(0, 200), rep(150, 100))
table(detect.nonwear.troiano(counts))
#> 
#> FALSE  TRUE 
#>   200   200 
```
