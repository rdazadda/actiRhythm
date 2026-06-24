# Calculate Onset Timing Variability

Measures day-to-day variability in the timing of the L5 and M10 activity
windows using circular statistics.

## Usage

``` r
.onset.timing.variability(daily_metrics)
```

## Arguments

- daily_metrics:

  Data frame with daily L5_start, M10_start

## Value

List with onset_timing_variability and component (L5/M10) variabilities

## Details

Uses circular standard deviation to properly handle midnight wraparound.
For example, if L5 starts at 23:00 one day and 01:00 the next, the
actual variability is 2 hours, not 22 hours.

NOTE: This is the mean of the circular SD of daily L5 and M10 onset
times. It is NOT the published Composite Phase Deviation (CPD) of
Fischer & Roenneberg (2016), which combines each day's precision
(deviation from the individual's own mean phase) and accuracy (deviation
from a reference phase) as mean(sqrt(precision^2 + accuracy^2)). It is
named here to reflect exactly what it computes so it is not mistaken for
the established CPD metric.
