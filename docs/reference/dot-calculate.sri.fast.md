# Calculate Sleep Regularity Index (SRI) - Fast Vectorized Version

The probability of being in the same sleep/wake state at any two time
points 24 hours apart, averaged across the recording period.

## Usage

``` r
.calculate.sri.fast(sleep_state, timestamps, epoch_length = 60)
```

## Arguments

- sleep_state:

  Character vector of sleep states ("S" or "W")

- timestamps:

  POSIXct timestamps

- epoch_length:

  Epoch length in seconds

## Value

Numeric SRI value (-100 to 100, higher = more regular)

## Details

SRI = 200 \* P(same state at t and t-24h) - 100

A perfectly regular sleeper (same schedule every day) scores 100. Random
sleep/wake patterns score about 0. Perfectly anti-regular patterns
(opposite states) score -100.

## References

Phillips AJK, Clerx WM, O'Brien CS, Sano A, Barger LK, Picard RW,
Lockley SW, Klerman EB, Czeisler CA (2017). “Irregular sleep/wake
patterns are associated with poorer academic performance and delayed
circadian and sleep/wake timing.” *Scientific Reports*, **7**(1), 3216.
[doi:10.1038/s41598-017-03171-4](https://doi.org/10.1038/s41598-017-03171-4)
.
