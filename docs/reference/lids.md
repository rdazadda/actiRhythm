# Locomotor Inactivity During Sleep (LIDS)

Quantifies the ultradian (sleep-cycle) oscillation of inactivity during
sleep, following Winnebeck et al. (2018). Within each sleep period the
activity is transformed to LIDS (`100 / (activity + 1)`), smoothed with
a centered 30-minute moving average, and fit with an ultradian cosine
over a scan of candidate periods; the best period is the one maximising
the Munich Rhythmicity Index (MRI).

## Usage

``` r
lids(
  counts,
  timestamps,
  sleep_periods,
  epoch_length = 60,
  smooth_minutes = 30,
  period_min = 30,
  period_max = 180,
  period_step = 5
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct (or parseable) timestamps, one per count.

- sleep_periods:

  Data frame with `in_bed_time` and `out_bed_time` (same schema as
  [`social.jet.lag`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md));
  one LIDS fit per period.

- epoch_length:

  Epoch length in seconds (default 60).

- smooth_minutes:

  Centered moving-average window in minutes (default 30).

- period_min, period_max, period_step:

  Period scan grid in minutes (defaults 30, 180, 5).

## Value

An object of class `actiRhythm_lids`: a per-period data frame (period in
minutes, MRI, amplitude, Pearson r, ...) and the mean period and MRI
across periods.

## References

Winnebeck EC, Fischer D, Leise T, Roenneberg T (2018). “Dynamics and
ultradian structure of human sleep in real life.” *Current Biology*,
**28**(1), 49–59.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)
.
