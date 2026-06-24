# Write Activity Counts to an .agd File

Writes a per-epoch counts data frame to a SQLite `.agd` file (the
standard `settings` + `data` schema), readable by
[`read.agd`](https://rdazadda.github.io/actiRhythm/reference/read.agd.md).

## Usage

``` r
write.agd(
  data,
  path,
  epoch_length = 60,
  start_time = NULL,
  device_serial = "",
  device_name = "wGT3XBT",
  subject_name = "",
  sample_rate = 30,
  mode = 61
)
```

## Arguments

- data:

  Data frame with `timestamp` (POSIXct) and `axis1`, `axis2`, `axis3`;
  optional `steps`, `lux` and `inclineOff`/`Standing`/`Sitting`/`Lying`.

- path:

  Output `.agd` path.

- epoch_length:

  Epoch length in seconds (default 60).

- start_time:

  Recording start (POSIXct); defaults to the first timestamp.

- device_serial, device_name, subject_name, sample_rate, mode:

  Metadata for the settings table.

## Value

The output `path`, invisibly.
