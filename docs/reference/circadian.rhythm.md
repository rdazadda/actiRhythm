# Circadian Rhythm Analysis

Computes non-parametric circadian metrics as in GGIR, ActCR, and
nparACT. Uses a minute-level sliding window for L5/M10, plus the Sleep
Regularity Index and several phase metrics.

## Usage

``` r
circadian.rhythm(
  counts,
  timestamps,
  sleep_state = NULL,
  wear_time = NULL,
  min_valid_hours = 10,
  epoch_length = 60,
  calculate_sri = TRUE,
  use_cpp = TRUE
)
```

## Arguments

- counts:

  Numeric vector of activity counts (minute-level recommended for
  accuracy)

- timestamps:

  POSIXct vector of epoch timestamps

- sleep_state:

  Optional character vector of sleep states ("S" or "W") for SRI
  calculation

- wear_time:

  Optional logical vector indicating wear time (TRUE = worn)

- min_valid_hours:

  Numeric. Valid-day criterion (GGIR includedaycrit): minimum wear hours
  for a day to count. Applied only when `wear_time` is given; default
  `10`. Set `0`/`NULL` to disable.

- epoch_length:

  Numeric. Epoch length in seconds (default: 60)

- calculate_sri:

  Logical. Calculate Sleep Regularity Index? (default: TRUE if
  sleep_state provided)

- use_cpp:

  Logical. Use C++ backend for faster computation? (default: TRUE)

## Value

List with class 'actiRhythm_circadian' containing:

- L5, L5_start, L5_start_hour:

  Least active 5-hour average and timing (sliding window)

- M10, M10_start, M10_start_hour:

  Most active 10-hour average and timing (sliding window)

- L1, M1:

  Least/most active 1-hour for additional granularity

- RA:

  Relative amplitude: (M10-L5)/(M10+L5), range 0-1

- IS:

  Interdaily stability: day-to-day consistency (0-1, higher=more stable)

- IV:

  Intradaily variability: within-day fragmentation (0 for a sine wave,
  near 2 for noise)

- phi:

  First-order autocorrelation at 1-hour lag

- SRI:

  Sleep Regularity Index (-100 to 100, higher=more regular)

- onset_timing_variability:

  Mean circular SD of daily L5/M10 onset times (day-to-day phase
  variability). NOT the published Fischer/Roenneberg CPD.

- hourly_profile:

  Mean activity by hour of day

- daily_metrics:

  Per-day L5, M10, RA values

## Details

**Non-Parametric Metrics (van Someren et al., 1999):**

L5/M10 use a **minute-level sliding window**, not hourly aggregation, so
timing is resolved to the minute rather than the hour.

- **L5** = Average activity during the least active 5 consecutive hours.
  Uses sliding window across minute-level data for precise timing.

- **M10** = Average activity during the most active 10 consecutive
  hours.

- **RA** = Relative amplitude = (M10 - L5) / (M10 + L5). Ranges 0 to 1;
  higher values mean stronger day-night amplitude.

- **IS** = Interdaily stability. Measures coupling of the rhythm to
  stable zeitgebers. Range 0 (Gaussian noise) to 1 (perfect stability).

- **IV** = Intradaily variability. Measures rhythm fragmentation. About
  0 for a perfect sine wave, about 2 for Gaussian noise, and can exceed
  2 with ultradian rhythms.

**Sleep Regularity Index (Phillips et al., 2017):**

SRI = probability of being in same sleep/wake state at any two time
points 24 hours apart. Range -100 to +100, with 100 indicating perfect
regularity.

## See also

[`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
for standalone SRI calculation,
[`social.jet.lag`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md)
for social jet lag calculation

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
result <- circadian.rhythm(counts$axis1, counts$timestamp)
print(result)
#> 
#> Circadian Rhythm Analysis
#> 
#> Data Summary
#>   Days analyzed:            8
#>   Valid circadian days:     6
#>   Epoch length:             60 seconds
#> 
#> Non-Parametric Metrics (IS/IV: Witting et al. 1990; RA/L5/M10: van Someren et al. 1999)
#>   L5 (least active 5h):     6.10 counts/min, onset 01:28
#>   M10 (most active 10h):    602.93 counts/min, onset 12:02
#>   L1 (least active 1h):     3.10 counts/min, onset 02:35
#>   M1 (most active 1h):      1175.51 counts/min, onset 17:43
#>   Relative Amplitude (RA):  0.9800 (range 0-1, higher=stronger rhythm)
#>   Interdaily Stability (IS): 0.2279 (range 0-1, higher=more consistent)
#>   Intradaily Variability (IV): 1.0008 (near 0 = sine, near 2 = noise)
#>   Phi (autocorrelation):    0.4151 (higher=more predictable)
#> 
#> Sleep-Based & Variability Metrics
#>   Sleep Regularity Index:   Not calculated (requires sleep_state input)
#>   Onset timing variability: 2.23 hours
#>   L5 timing variability:    0.77 hours (circular SD)
#>   M10 timing variability:   3.70 hours (circular SD)
#> 
#> References: Witting (1990), van Someren (1999)
#> 
# }
```
