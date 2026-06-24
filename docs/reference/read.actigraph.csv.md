# Read an ActiGraph (ActiLife) Epoch CSV

Reads an ActiLife epoch (count) CSV export. The metadata header is
parsed for the start time and epoch length, the count columns are
recognized by name or position, and a regular timestamp grid is built.
The result is a data frame with a POSIXct `timestamp` and the available
count columns, the same shape as
[`agd.counts`](https://rdazadda.github.io/actiRhythm/reference/agd.counts.md),
so it feeds the analysis functions directly.

## Usage

``` r
read.actigraph.csv(
  filepath,
  tz = "UTC",
  date_format = "%m/%d/%Y",
  epoch_length = 60
)
```

## Arguments

- filepath:

  Path to the ActiLife CSV.

- tz:

  Time zone for the timestamps (default `"UTC"`; ActiLife stores local
  clock time without a zone).

- date_format:

  Date format of the header Start Date (default `"%m/%d/%Y"`, the
  ActiLife M/d/yyyy default).

- epoch_length:

  Fallback epoch length in seconds if the header lacks an Epoch Period
  line (default 60).

## Value

A data frame with `timestamp` (POSIXct) and the count columns present
(`axis1`, `axis2`, `axis3`, `vm`, `steps`, `lux`).

## See also

[`counts.from.data.frame`](https://rdazadda.github.io/actiRhythm/reference/counts.from.data.frame.md),
[`read.agd`](https://rdazadda.github.io/actiRhythm/reference/read.agd.md)
