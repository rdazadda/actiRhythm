# Activity Onset and Offset (Relative-Difference Phase Markers)

Finds the daily activity onset and offset by a relative-difference
contrast on the averaged 24-hour profile: the onset is where mean
activity rises most sharply (the relative difference of the window after
versus before is largest) and the offset is where it falls most sharply.
These are non-cosinor, non-changepoint phase markers, a
normalized-contrast edge detector on the daily profile rather than a
published actigraphy algorithm.

## Usage

``` r
activity.onset.offset(counts, timestamps, window_hours = 6)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- window_hours:

  Half-window, in hours, compared before versus after each minute
  (default 6).

## Value

An object of class `actiRhythm_aont`: `onset_h` and `offset_h` (clock
hours) and the relative-difference profile. Never errors.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H"))
activity.onset.offset(ifelse(h >= 23 | h < 7, 5, 300), ts)
#> Activity Onset / Offset
#> 
#>   Activity onset:  07.00 h
#>   Activity offset: 23.00 h
#> 
```
