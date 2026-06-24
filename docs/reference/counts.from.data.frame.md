# Activity Counts from a Data Frame

A device-neutral entry point: pull a count column (and optionally a time
column) out of any data frame and return a tidy `timestamp`/`counts`
frame ready for the analysis functions. Use it for non-ActiGraph counts
or any pre-extracted series.

## Usage

``` r
counts.from.data.frame(
  df,
  count_col = "axis1",
  time_col = NULL,
  epoch_length = 60,
  tz = "UTC",
  start = "1970-01-01"
)
```

## Arguments

- df:

  A data frame.

- count_col:

  Name or index of the count column (default `"axis1"`).

- time_col:

  Name or index of the timestamp column; if `NULL`, timestamps are
  synthesized from `start` by `epoch_length` (with a warning).

- epoch_length:

  Epoch length in seconds for synthesized timestamps (default 60).

- tz:

  Time zone (default `"UTC"`).

- start:

  Start time for synthesized timestamps (default `"1970-01-01"`).

## Value

A data frame with `timestamp` (POSIXct) and `counts`.

## See also

[`read.actigraph.csv`](https://rdazadda.github.io/actiRhythm/reference/read.actigraph.csv.md)

## Examples

``` r
df <- data.frame(activity = c(0, 50, 300), clock = c("2024-01-01 00:00:00",
  "2024-01-01 00:01:00", "2024-01-01 00:02:00"))
counts.from.data.frame(df, count_col = "activity", time_col = "clock")
#>             timestamp counts
#> 1 2024-01-01 00:00:00      0
#> 2 2024-01-01 00:01:00     50
#> 3 2024-01-01 00:02:00    300
```
