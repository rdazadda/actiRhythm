# Export a Full Circadian Analysis to an Excel Workbook

Runs the complete actiRhythm circadian analysis on one activity time
series and writes a multi-sheet `.xlsx` workbook: a one-row **Summary**
of every metric (nonparametric IS/IV/RA/L5/M10, cosinor, the rhythmicity
F-test, the Lomb-Scargle period with its bootstrap confidence interval,
the chi-square periodogram, DFA, multifractal DFA, multiscale entropy,
rest-activity state transitions, and, when sleep is supplied, the Sleep
Regularity Index, social jet lag and LIDS), plus detail sheets for the
hourly profile, both periodograms, the fluctuation and multifractal
spectra, the transition curves, the LIDS fits, and a data dictionary.

## Usage

``` r
circadian.workbook(
  activity,
  timestamps,
  file = NULL,
  sleep_state = NULL,
  sleep_periods = NULL,
  wear_time = NULL,
  epoch_length = 60,
  include_period_ci = TRUE,
  n_boot = 200
)
```

## Arguments

- activity:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- file:

  Output `.xlsx` path; if `NULL` the workbook object is returned without
  writing.

- sleep_state:

  Optional per-epoch sleep/wake vector (turns on the SRI).

- sleep_periods:

  Optional sleep-period data frame with `in_bed_time`/`out_bed_time`
  (turns on social jet lag and LIDS).

- wear_time:

  Optional logical wear-time mask.

- epoch_length:

  Epoch length in seconds (default 60).

- include_period_ci:

  Whether to bootstrap the period CI (default TRUE).

- n_boot:

  Bootstrap replicates for the period CI (default 200).

## Value

The `openxlsx` workbook object, invisibly.

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
h  <- as.numeric(format(ts, "%H"))
act <- pmax(0, 100 + 80 * cos(2 * pi * (h - 14) / 24) + rnorm(length(ts), 0, 20))
circadian.workbook(act, ts, file = tempfile(fileext = ".xlsx"), include_period_ci = FALSE)
# }
```
