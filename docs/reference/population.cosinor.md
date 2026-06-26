# Population-Mean Cosinor (Bingham)

Pools single-subject cosinor fits into a group-mean rhythm with
confidence intervals, following Bingham et al. (1982). Fits each
subject's averaged 24-hour profile with the same weighted-least-squares
engine as
[`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
averages the linearized cos/sin coefficients across subjects, and
returns the group MESOR, amplitude, and acrophase with Bingham
confidence intervals.

## Usage

``` r
population.cosinor(
  activity,
  timestamps,
  subject,
  group = NULL,
  period = 24,
  level = 0.95,
  min_valid_hours = 12
)
```

## Arguments

- activity:

  Numeric activity vector with all subjects stacked together.

- timestamps:

  POSIXct timestamps, one per value.

- subject:

  Subject identifier, one per value.

- group:

  Optional group identifier, one per value; when supplied a population
  cosinor is returned for each group.

- period:

  Rhythm period in hours (default 24).

- level:

  Confidence level (default 0.95).

- min_valid_hours:

  Minimum profile hours for a subject to be included (default 12).

## Value

An object of class `actiRhythm_population_cosinor` (group MESOR,
amplitude, acrophase with Bingham CIs and a `conf_interval_valid` flag),
or a named list of them (class `actiRhythm_population_cosinor_list`)
when `group` is supplied.

## References

Bingham C, Arbogast B, Cornelissen Guillaume G, Lee JK, Halberg F
(1982). “Inferential statistical methods for estimating and comparing
cosinor parameters.” *Chronobiologia*, **9**(4), 397–439.

## Examples

``` r
set.seed(1)
hrs <- 0:23
act <- ts <- subj <- NULL
for (i in 1:6) {
  y <- 100 + 40 * cos(2 * pi * (hrs - (8 + i / 3)) / 24) + rnorm(24, 0, 4)
  act <- c(act, y); subj <- c(subj, rep(paste0("S", i), 24))
  ts <- c(ts, as.POSIXct("2024-01-01", tz = "UTC") + hrs * 3600)
}
population.cosinor(act, as.POSIXct(ts, tz = "UTC", origin = "1970-01-01"), subj)
#> Population-Mean Cosinor (Bingham)
#> 
#>   Subjects:   6
#>   Period:     24 h
#>   MESOR:      100.17  [99.32, 101.01]
#>   Amplitude:  39.68  [38.54, 40.81]
#>   Acrophase:  9.71 h  [9.07, 10.36]
```
