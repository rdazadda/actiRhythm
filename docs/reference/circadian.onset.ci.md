# Confidence Intervals for L5/M10 Onset Timing

Percentile bootstrap confidence interval for the mean daily onset time
of a circadian phase marker (e.g. L5 or M10), using circular resampling
of the per-day onsets.

## Usage

``` r
circadian.onset.ci(onset_hours, level = 0.95, n_boot = 2000)
```

## Arguments

- onset_hours:

  Numeric vector of daily onset times in decimal hours.

- level:

  Confidence level (default 0.95).

- n_boot:

  Bootstrap replicates (default 2000).

## Value

List with `mean_onset`, `ci_lower`, `ci_upper` (decimal hours) and
`n_days`.
