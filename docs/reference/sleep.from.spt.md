# Sleep Parameters from SPT and Sustained-Inactivity Bouts

Intersects the per-epoch sustained-inactivity sleep score
([`sib.vanhees`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md))
with the sleep-period-time window
([`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md))
to derive per-night sleep parameters: total sleep time, onset and wake,
wake-after-sleep-onset, efficiency, and awakenings.

## Usage

``` r
sleep.from.spt(spt, sib, timestamps, epoch_length = 5)
```

## Arguments

- spt:

  An `actiRhythm_spt` object from
  [`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md).

- sib:

  A character `"S"`/`"W"` vector from
  [`sib.vanhees`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md).

- timestamps:

  POSIXct timestamps, one per epoch (matching `sib`).

- epoch_length:

  Epoch length in seconds (default 5).

## Value

An object of class `actiRhythm_sleep`: a data frame with one row per
night (`date`, sleep `onset`, `offset`, `tst` hours, `waso` minutes,
`efficiency`, `n_awakenings`, `mid_sleep`).

## References

van Hees VT, Sabia S, Jones SE, Wood AR, Anderson KN, Kivimaki M,
Frayling TM, Pack AI, Bucan M, Trenell MI, Mazzotti DR, Gehrman PR,
Singh-Manoux BA, Weedon MN (2018). “Estimating sleep parameters using an
accelerometer without sleep diary.” *Scientific Reports*, **8**, 12975.
[doi:10.1038/s41598-018-31266-z](https://doi.org/10.1038/s41598-018-31266-z)
.

## See also

[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md),
[`sib.vanhees`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01 12:00", tz = "UTC"), by = 5, length.out = 17280)
h <- as.numeric(format(ts, "%H")); night <- h >= 23 | h < 7
set.seed(1)
anglez <- ifelse(night, -60, -30) + rnorm(17280, 0, ifelse(night, 0.02, 20))
spt <- rest.spt(anglez, ts, epoch_length = 5)
sib <- sib.vanhees(anglez, epoch_length = 5)
sleep.from.spt(spt, sib, ts, epoch_length = 5)
#>         date               onset              offset      tst waso efficiency
#> 1 2024-01-01 2024-01-01 23:00:05 2024-01-02 06:59:55 7.998611    0          1
#>   n_awakenings           mid_sleep
#> 1            0 2024-01-02 03:00:00
```
