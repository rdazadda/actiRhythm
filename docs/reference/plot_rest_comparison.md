# Rest-Detector Comparison Strip

Runs the four rest/sleep detectors on one recording and stacks their
detected rest bands over the activity series on a shared time axis, so
the user can see that the differing bout counts reflect different
questions (main night vs every bout vs latent state), not contradictory
accuracy. Returns a `ggplot` object and never errors.

## Usage

``` r
plot_rest_comparison(counts, timestamps)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

## Value

A `ggplot` object.

## See also

[`sleep.changepoints`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md),
[`rest.periods`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md),
[`rest.crespo`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md),
[`rest.hmm`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H"))
plot_rest_comparison(ifelse(h >= 23 | h < 7, 5, 300), ts)

# }
```
