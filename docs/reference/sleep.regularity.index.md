# Calculate Sleep Regularity Index (SRI) - Exported Version

Calculate Sleep Regularity Index (SRI) - Exported Version

## Usage

``` r
sleep.regularity.index(sleep_state, timestamps, epoch_length = 60)
```

## Arguments

- sleep_state:

  Character vector of sleep states ("S" or "W")

- timestamps:

  POSIXct timestamps (must be regular epochs)

- epoch_length:

  Epoch length in seconds (default 60)

## Value

Numeric SRI value (-100 to 100)

## References

Phillips AJK, et al. (2017). Scientific Reports, 7(1):3216.
