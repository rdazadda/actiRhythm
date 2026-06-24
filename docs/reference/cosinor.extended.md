# Multi-Component Cosinor Analysis with Harmonics

Fits a multi-component cosinor model with the fundamental (24h) and
harmonics (12h, 8h, etc.) to capture patterns the single cosinor misses,
such as bimodal rhythms.

## Usage

``` r
cosinor.extended(
  counts,
  timestamps,
  harmonics = c(24, 12),
  wear_time = NULL,
  min_valid_hours = 10
)
```

## Arguments

- counts:

  Numeric vector of activity counts

- timestamps:

  POSIXct vector of timestamps

- harmonics:

  Vector of periods to include (default: c(24, 12) for 24h + 12h)

- wear_time:

  Optional logical vector indicating valid wear time

- min_valid_hours:

  Numeric. Valid-day criterion (GGIR includedaycrit): minimum wear hours
  for a day to count; default 10, 0/NULL to disable.

## Value

A list with class "actiRhythm_cosinor_extended" containing:

- mesor:

  Rhythm-adjusted mean

- components:

  Data frame with amplitude, acrophase, acrophase_time for each harmonic

- dominant_period:

  Period with largest amplitude

- dominant_acrophase:

  Acrophase of dominant component

- r_squared:

  Model R-squared (should be higher than single-component)

- r_squared_improvement:

  Improvement over single 24h cosinor

## Details

The multi-component model is:

`Y(t) = M + sum_k[A_k * cos(2*pi*t/T_k - phi_k)]`

Where T_k are the periods (24h, 12h, 8h, etc.)

Use it for:

- Bimodal patterns (morning + evening peaks), captured by the 12h
  harmonic

- Complex daily routines with multiple activity bouts

- Shift workers or irregular schedules

## References

Cornelissen G. (2014). Cosinor-based rhythmometry. Theor Biol Med Model,
11:16. Refinetti R, et al. (2007). Procedures for numerical analysis of
circadian rhythms.

## Examples

``` r
# \donttest{
counts <- agd.counts(read.agd(example_agd()))
#> AGD file tables: awakenings, capsense, crouterEpoch, crouterMinute, data, filterCategory, filters, logDiaryTimes, logEventHistory, logEventType, proximity, settings, sleep, sqlite_sequence, wtvBouts 
#> Epochs loaded: 9919 
#> Sleep periods found: 4 
#> Awakenings found: 41 
#> Wear time bouts: 5 
#> Capsense samples: 9912 
result <- cosinor.extended(counts$axis1, counts$timestamp, harmonics = c(24, 12, 8))
print(result)
#> Multi-Component Cosinor Analysis
#> 
#> MESOR (mean):     330.12
#> Pattern type:     Strong 24h
#> R-squared:        0.5938 (59.4% of variance explained)
#> Improvement:      +15.7% over single 24h cosinor
#> 
#> Harmonic Components:
#>   24h: Amplitude=300.2, Peak=16:58 (17.0h), Power=74.0% *DOMINANT*
#>   12h: Amplitude=171.2, Peak=07:24 (7.4h), Power=24.1%
#>   8h: Amplitude=48.6, Peak=04:31 (4.5h), Power=1.9%
#> 
#> [!] Bimodal pattern: two activity peaks per day
#>     Estimated peaks at: around 07:00 and around 19:00
#> 
#> F-statistic: 4.14, p-value: 9.58e-03
# }
```
