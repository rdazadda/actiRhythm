# Sleep Regularity Index (SRI) - Phillips (2017) Epoch-of-Day x Day Matrix

Computes the Sleep Regularity Index exactly as defined by Phillips et
al. (2017) using the full epoch-of-day by day concordance matrix. This
is the published form of the SRI. Unlike a single fixed 24-hour-lag
comparison, every epoch is binned by its real clock time (epoch-of-day)
and by calendar date, a binary sleep/wake matrix is built, and the index
is the average agreement between the SAME clock time on CONSECUTIVE
calendar days, rescaled to the interval \[-100, 100\].

## Usage

``` r
sri.matrix(sleep_state, timestamps, epoch_length = 60)
```

## Arguments

- sleep_state:

  Sleep/wake state per epoch. Either a character vector of `"S"` (sleep)
  / `"W"` (wake), or a numeric/integer/logical vector where 1/TRUE =
  sleep and 0/FALSE = wake. `NA` marks unscored / non-wear epochs and is
  excluded from the concordance count. Must be the same length as
  `timestamps`.

- timestamps:

  POSIXct vector of epoch start times, one per element of `sleep_state`.
  Need not start at midnight and need not be perfectly contiguous;
  alignment is by real clock time.

- epoch_length:

  Epoch length in seconds (default 60). Must divide 86400 evenly;
  `M = 86400 / epoch_length` is the number of epochs per day.

## Value

A list with components:

- SRI:

  Numeric Sleep Regularity Index in \[-100, 100\] (rounded to 2 dp), or
  `NA_real_` if it cannot be computed.

- n_days:

  Integer number of distinct calendar days (N) spanned.

- n_valid_pairs:

  Integer number of consecutive-day epoch pairs that were actually
  compared (both cells non-NA).

- method:

  Character constant `"phillips_matrix"`.

## Details

Let `M = 86400 / epoch_length` be the number of epochs in a day and let
`N` be the number of calendar days spanned by the recording. Each epoch
is assigned an epoch-of-day index `i` (1..M) from its clock time
(hour/minute/second) and a day index `j` (1..N) from its calendar date.
These two indices populate a binary matrix `s[i, j]` where 1 = sleep and
0 = wake. The SRI is then

\$\$SRI = -100 + \frac{200}{M (N - 1)} \sum\_{i=1}^{M} \sum\_{j=1}^{N-1}
\mathbf{1}\\ s\_{i,j} = s\_{i,j+1} \\\$\$

i.e. the proportion of clock-time epochs that hold the same state on two
consecutive days, mapped linearly so that perfect regularity scores
+100, chance-level (independent) scoring about 0, and perfect
anti-regularity -100.

Robustness to gaps and partial wear: any consecutive-day pair
`(i, j) -> (i, j + 1)` in which either cell is `NA` (missing epoch,
non-wear, or a calendar day not represented in the data) is skipped and
NOT counted. The sum is divided by the actual number of valid pairs
rather than the theoretical `M * (N - 1)`, so missing data reduce
statistical power but do not bias the estimate.

Days are aligned by true clock time using
[`as.POSIXlt()`](https://rdrr.io/r/base/as.POSIXlt.html) (hour, minute,
second give the epoch-of-day; the calendar date gives the day index), so
the result is correct even when a recording does not start at midnight.

A `SRI` of `NA_real_` is returned (never an error) when there are fewer
than two days, no valid consecutive-day pairs, or no usable input. The
return structure is always the same list.

## References

Phillips AJK, Clerx WM, O'Brien CS, Sano A, Barger LK, Picard RW,
Lockley SW, Klerman EB, Czeisler CA (2017). Irregular sleep/wake
patterns are associated with poorer academic performance and delayed
circadian and sleep/wake timing. *Scientific Reports*, 7(1):3216.
[doi:10.1038/s41598-017-03171-4](https://doi.org/10.1038/s41598-017-03171-4)

## Examples

``` r
# A perfectly regular sleeper: same 8h sleep block every day -> SRI near 100
ts <- seq(as.POSIXct("2024-01-01 00:00:00", tz = "UTC"),
          by = 60, length.out = 1440 * 4)
tod <- as.POSIXlt(ts)$hour
state <- ifelse(tod < 8, "S", "W")   # asleep 00:00-08:00 every day
sri.matrix(state, ts, epoch_length = 60)$SRI
#> [1] 100
```
