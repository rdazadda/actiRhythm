# Read ActiGraph .agd File

Reads epoch-level AGD files pre-processed in ActiLife. AGD files hold
activity counts already computed from raw acceleration data.

## Usage

``` r
read.agd(
  filepath,
  include_sleep = TRUE,
  include_wear_time = TRUE,
  verbose = TRUE
)
```

## Arguments

- filepath:

  Path to .agd file

- include_sleep:

  Logical. Include sleep analysis data if available? (default: TRUE)

- include_wear_time:

  Logical. Include wear time validation data if available? (default:
  TRUE)

- verbose:

  Logical. Print progress messages? (default: TRUE)

## Value

List containing:

- `data` - Data frame with epoch-level activity counts (axis1, axis2,
  axis3, steps, etc.)

- `settings` - Data frame with device and subject settings

- `sleep` - Sleep periods detected by ActiLife (if available and
  requested)

- `awakenings` - Awakening events during sleep (if available and
  requested)

- `wear_time` - Wear time validation bouts (if available and requested)

- `capsense` - Capacitive sensor data for on-body detection (if
  available)

- `tables` - Character vector of all tables in the AGD file

## Details

AGD files are SQLite databases created by ActiLife software. They
contain pre-processed activity counts at user-specified epoch lengths
(e.g., 60 seconds).
