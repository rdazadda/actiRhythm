# Sadeh Sleep/Wake Scoring

Classifies each epoch as sleep or wake from activity counts with the
Sadeh algorithm (Sadeh et al. 1994), validated on adults and adolescents
on one-minute epochs.

## Usage

``` r
sleep.sadeh(
  counts,
  epoch_length = 60,
  wake_threshold = 0,
  clip = NULL,
  na_action = c("na", "wake", "zero")
)
```

## Arguments

- counts:

  Numeric vector of activity counts (vertical axis).

- epoch_length:

  Epoch length in seconds (default 60). The algorithm was validated on
  60-second epochs; other lengths raise a warning.

- wake_threshold:

  Sleep-index cut: an epoch is sleep when \\SI \ge\\ `wake_threshold`
  (default 0, Sadeh 1994; ActiLife uses -4).

- clip:

  Optional upper cap on counts before scoring (default `NULL`, no cap;
  ActiLife and pyActigraphy use 300). Affects AVG, SD, and LG.

- na_action:

  How NA-count epochs appear in the output: `"na"` (default) emits `NA`,
  so non-wear gaps are not read as sleep; `"wake"` scores them wake;
  `"zero"` scores them from a zero count.

## Value

Character vector of states, `"S"` (sleep) or `"W"` (wake), the same
length as `counts`.

## Details

The sleep index uses an eleven-epoch window (five before, the current
epoch, and five after): \$\$SI = 7.601 - 0.065 \cdot AVG - 1.08 \cdot
NATS - 0.056 \cdot SD - 0.703 \cdot LG\$\$ where \\AVG\\ is the window
mean, \\NATS\\ the number of epochs with counts in \[50, 100), \\SD\\
the standard deviation over the current and five preceding epochs, and
\\LG = \log(\mathrm{count} + 1)\\. An epoch is scored sleep when \\SI
\ge 0\\ (Sadeh et al. 1994). For ActiLife/ActiGraph parity, set
`wake_threshold = -4` and `clip = 300`.

## References

Sadeh A, Sharkey KM, Carskadon MA (1994). “Activity-based sleep-wake
identification: an empirical test of methodological issues.” *Sleep*,
**17**(3), 201–207.
[doi:10.1093/sleep/17.3.201](https://doi.org/10.1093/sleep/17.3.201) .

## See also

[`sleep.cole.kripke`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md),
[`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)

## Examples

``` r
# \donttest{
agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
state <- sleep.sadeh(agd$axis1)
table(state)
#> state
#>    S    W 
#> 6757 3162 
# }
```
