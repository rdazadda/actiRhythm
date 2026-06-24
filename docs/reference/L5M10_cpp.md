# Calculate L5/M10 Circadian Metrics

Finds least active 5 hours (L5) and most active 10 hours (M10) using van
Someren (1999) average-profile method.

## Usage

``` r
L5M10_cpp(minute_data, window_L5 = 300L, window_M10 = 600L, start_minute = 0L)
```

## Arguments

- minute_data:

  Numeric vector of minute-level activity

- window_L5:

  L5 window in minutes (default: 300)

- window_M10:

  M10 window in minutes (default: 600)

- start_minute:

  Start minute of day (0-1439)

## Value

List with L5_value, L5_onset, M10_value, M10_onset, RA
