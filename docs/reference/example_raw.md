# Synthetic Raw Acceleration Recording

Generates a deterministic synthetic triaxial raw acceleration recording
(in g) with a day/night posture cycle, daytime movement and posture
changes, still (sleeping) nights, and a slight built-in miscalibration.
It is the file-free stand-in for a real raw file in the examples: feed
it to
[`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md),
[`circadian.raw`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md),
or the z-angle sleep pipeline. Raw acceleration is far too large to ship
as data, so it is generated on demand rather than bundled.

## Usage

``` r
example_raw(days = 2, fs = 30, device_off = 0)
```

## Arguments

- days:

  Recording length in days (default 2, giving two nights).

- fs:

  Sample rate in Hz (default 30).

- device_off:

  Days of a still, taken-off device to append at the end (default 0).
  Use it to exercise
  [`detect.nonwear.raw`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
  and the non-wear gate of
  [`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md).

## Value

A data frame with a POSIXct `time` column and `x`/`y`/ `z` acceleration
in g, with the sample rate in the `"fs"` attribute.

## See also

[`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md),
[`circadian.raw`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md),
[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)

## Examples

``` r
raw <- example_raw(days = 1)
str(raw)
#> 'data.frame':    2592000 obs. of  4 variables:
#>  $ time: POSIXct, format: "2024-01-01 12:00:00" "2024-01-01 12:00:00" ...
#>  $ x   : num  0.1094 0.149 0.0827 0.108 0.1341 ...
#>  $ y   : num  -1.045 -1.01 -0.971 -0.972 -1.016 ...
#>  $ z   : num  0.261 0.236 0.247 0.299 0.231 ...
#>  - attr(*, "fs")= num 30
```
