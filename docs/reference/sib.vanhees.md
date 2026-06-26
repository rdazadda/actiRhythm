# Sustained-Inactivity-Bout Sleep Scoring from the z-Angle

Scores each epoch sleep (`"S"`) or wake (`"W"`) from the z-angle by the
van Hees et al. (2015) sustained-inactivity-bout rule: an interval with
no arm-posture change exceeding `angle_thresh` degrees for at least
`time_thresh` minutes is sustained inactivity (sleep). The output has
the same shape as
[`sleep.cole.kripke`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md),
so it feeds
[`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md),
[`lids`](https://rdazadda.github.io/actiRhythm/reference/lids.md) and
the rest consumers directly. Intersect it with
[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
via
[`sleep.from.spt`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md).

## Usage

``` r
sib.vanhees(anglez, angle_thresh = 5, time_thresh = 5, epoch_length = 5)
```

## Arguments

- anglez:

  Numeric z-angle (degrees) per epoch.

- angle_thresh:

  Posture-change threshold in degrees (default 5).

- time_thresh:

  Minimum sustained-inactivity duration in minutes (default 5).

- epoch_length:

  Epoch length in seconds (default 5).

## Value

Character vector of `"S"`/`"W"`, one per epoch. A day with fewer than 10
posture changes is scored all `"S"` (sustained inactivity, to be gated
by the SPT window and wear time).

## References

van Hees VT, Sabia S, Anderson KN, Denton SJ, Oliver J, Catt M, Abell
JG, Kivimaki M, Trenell MI, Singh-Manoux A (2015). “A novel, open access
method to assess sleep duration using a wrist-worn accelerometer.” *PLoS
ONE*, **10**(11), e0142533.
[doi:10.1371/journal.pone.0142533](https://doi.org/10.1371/journal.pone.0142533)
.

## See also

[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md),
[`sleep.from.spt`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)

## Examples

``` r
# A long still stretch scores as sustained inactivity (sleep)
set.seed(1)
anglez <- c(rnorm(3000, -60, 0.02), rnorm(3000, 0, 20))   # still, then active
table(sib.vanhees(anglez, epoch_length = 5))
#> 
#>    W 
#> 6000 
```
