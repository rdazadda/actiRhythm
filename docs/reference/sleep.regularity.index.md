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

Phillips AJK, Clerx WM, O'Brien CS, Sano A, Barger LK, Picard RW,
Lockley SW, Klerman EB, Czeisler CA (2017). “Irregular sleep/wake
patterns are associated with poorer academic performance and delayed
circadian and sleep/wake timing.” *Scientific Reports*, **7**(1), 3216.
[doi:10.1038/s41598-017-03171-4](https://doi.org/10.1038/s41598-017-03171-4)
.
