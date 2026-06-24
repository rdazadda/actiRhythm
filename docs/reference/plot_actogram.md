# Double-Plotted Actogram

Draws a classic double-plotted actogram: one row per calendar day, each
row showing 48 hours (the day itself on the left, the following day on
the right) so circadian phase can be traced down the diagonal. Activity
is rendered as a per-minute raster (darker = more active). Non-wear (via
`wear_time`) and missing time are left blank, and any skipped calendar
days appear as empty rows.

## Usage

``` r
plot_actogram(
  counts,
  timestamps,
  epoch_length = NULL,
  wear_time = NULL,
  double_plot = TRUE,
  L5_onset = NULL,
  M10_onset = NULL,
  sleep_mask = NULL,
  scale = c("linear", "sqrt")
)
```

## Arguments

- counts:

  Numeric vector of activity counts on a regular epoch grid.

- timestamps:

  A `POSIXct` vector the same length as `counts`.

- epoch_length:

  Epoch length in seconds. If `NULL` (default) it is inferred from the
  median spacing of `timestamps`.

- wear_time:

  Optional logical vector (`TRUE` = worn) the same length as `counts`;
  non-wear epochs are blanked.

- double_plot:

  Logical; draw the 48-hour double plot (default `TRUE`) or a single
  24-hour plot.

- L5_onset, M10_onset:

  Optional L5 / M10 onset to overlay as a dashed vertical phase line.
  Accepts a decimal hour, an `"HH:MM"` string, or a `POSIXct`.

- sleep_mask:

  Optional logical/character vector (`TRUE`/"S" = asleep) the same
  length as `counts`; sleep is shaded over the raster.

- scale:

  Fill scaling: `"linear"` (default) or `"sqrt"` to compress a
  heavy-tailed activity range.

## Value

A `ggplot` object. Never errors; returns an annotated empty plot on
insufficient data.
