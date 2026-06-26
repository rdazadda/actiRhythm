# Cosinor Analysis for Circadian Rhythm

Fits a single-component cosinor model to activity data and returns the
parametric circadian parameters.

## Usage

``` r
cosinor.analysis(
  counts,
  timestamps,
  period = 24,
  wear_time = NULL,
  min_valid_hours = 10,
  transform = c("none", "log1p")
)
```

## Arguments

- counts:

  Numeric vector of activity counts

- timestamps:

  POSIXct vector of timestamps

- period:

  Period in hours (default: 24 for circadian rhythm)

- wear_time:

  Optional logical vector indicating wear time

- min_valid_hours:

  Numeric. Valid-day criterion (GGIR includedaycrit): minimum wear hours
  for a day to count; default 10, 0/NULL to disable.

- transform:

  Input transform: `"none"` (default) or `"log1p"`, GGIR's
  `log(ENMO_mg + 1)` pre-transform for raw acceleration metrics (g-unit
  input is scaled to mg first).

## Value

List with class 'actiRhythm_cosinor' containing:

- mesor:

  Rhythm-adjusted mean (MESOR)

- amplitude:

  Half the peak-to-trough difference

- acrophase:

  Time of peak activity (hours, 0-24)

- acrophase_time:

  Acrophase as HH:MM string

- r_squared:

  Goodness of fit (coefficient of determination)

- f_statistic:

  F-statistic for model significance

- p_value:

  P-value for model significance

## Details

The single-component cosinor model fits: \$\$Y(t) = M + A \cdot cos(2\pi
t / T + \phi)\$\$

Where:

- M = mesor (rhythm-adjusted mean)

- A = amplitude (half peak-to-trough difference)

- T = period (typically 24 hours)

- φ = acrophase (phase angle, converted to time)

The model is fit using linear least squares on the linearized form: Y(t)
= M + β₁·cos(2πt/T) + β₂·sin(2πt/T)

## References

Nelson W, Tong YL, Lee JK, Halberg F (1979). “Methods for
cosinor-rhythmometry.” *Chronobiologia*, **6**(4), 305–323.

Cornelissen G (2014). “Cosinor-based rhythmometry.” *Theoretical Biology
and Medical Modelling*, **11**, 16.
[doi:10.1186/1742-4682-11-16](https://doi.org/10.1186/1742-4682-11-16) .

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
result <- cosinor.analysis(counts$axis1, counts$timestamp)
print(result)
#> Cosinor Analysis Results
#> 
#> Period:     24 hours
#> N obs:      9919 (8.0 days)
#> 
#> Parameters:
#>   MESOR:      327.25 (rhythm-adjusted mean)
#>   Amplitude:  295.63 (half peak-to-trough)
#>   Acrophase:  16:55 (16.92 h, time of peak)
#> 
#> Model Fit:
#>   R-squared:  0.4371
#>   F-statistic: 8.15
#>   P-value:    2.4e-03
# }
```
