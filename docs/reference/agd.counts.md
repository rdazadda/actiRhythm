# Extract Counts from .agd Data

Extracts activity counts, timestamps, steps, and inclinometer data (if
available) from AGD file data. Inclinometer data can be used to enhance
sleep detection for AGD files that lack raw acceleration.

## Usage

``` r
agd.counts(agd_data, convert.timestamps = TRUE, include_inclinometer = TRUE)
```

## Arguments

- agd_data:

  List returned from read.agd()

- convert.timestamps:

  Logical. Convert ActiGraph timestamps to POSIXct (default: TRUE)

- include_inclinometer:

  Logical. Include inclinometer data if available? (default: TRUE)

## Value

Data frame with counts per minute, timestamps, and optional inclinometer
data:

- `timestamp` - POSIXct timestamp (if convert.timestamps = TRUE)

- `axis1` - Vertical axis counts (Y-axis)

- `axis2` - Horizontal axis counts (X-axis)

- `axis3` - Anterior-posterior axis counts (Z-axis)

- `steps` - Step counts

- `lux` - Light intensity (if available)

- `inclineOff` - Off-body indicator (if available)

- `inclineStanding` - Standing indicator (if available)

- `inclineSitting` - Sitting indicator (if available)

- `inclineLying` - Lying indicator (if available)

## Details

Inclinometer data is particularly useful for sleep detection in AGD
files:

- `inclineLying = 1` indicates the device detected lying posture

- Combined with low activity (Cole-Kripke sleep), this sharpens
  specificity

- Not all AGD files contain inclinometer data (depends on ActiLife
  processing settings)
