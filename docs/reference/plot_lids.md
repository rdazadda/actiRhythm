# LIDS Ultradian Cycle Curve

Plots the smoothed Locomotor Inactivity During Sleep (LIDS) signal
across one sleep period with the best-fit ultradian cosine overlaid, the
companion figure to
[`lids`](https://rdazadda.github.io/actiRhythm/reference/lids.md). The
cosine's period and Munich Rhythmicity Index summarise the sleep-cycle
oscillation. Returns a `ggplot` object and never errors.

## Usage

``` r
plot_lids(counts, timestamps, sleep_periods, period = 1L, ...)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per count.

- sleep_periods:

  Data frame with `in_bed_time` and `out_bed_time` (as in
  [`lids`](https://rdazadda.github.io/actiRhythm/reference/lids.md)).

- period:

  Which sleep period to plot (default 1).

- ...:

  Passed to
  [`lids`](https://rdazadda.github.io/actiRhythm/reference/lids.md).

## Value

A `ggplot` object.

## References

Winnebeck EC, Fischer D, Leise T, Roenneberg T (2018). “Dynamics and
ultradian structure of human sleep in real life.” *Current Biology*,
**28**(1), 49–59.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)
.

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01 22:00", tz = "UTC"), by = 60, length.out = 480)
sp <- data.frame(in_bed_time = ts[1], out_bed_time = ts[480])
plot_lids(50 + 45 * cos(2 * pi * seq_along(ts) / 110), ts, sp)

# }
```
