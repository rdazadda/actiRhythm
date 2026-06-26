# Cole-Kripke Sleep/Wake Scoring

Classifies each epoch as sleep or wake from activity counts with the
Cole-Kripke algorithm (Cole et al. 1992), validated for adults on
one-minute epochs, with the optional Webster rescoring rules (Webster et
al. 1982). Returns the per-epoch `sleep_state` that the regularity and
locomotor-sleep metrics read
([`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md),
[`lids`](https://rdazadda.github.io/actiRhythm/reference/lids.md)).

## Usage

``` r
sleep.cole.kripke(
  counts,
  apply_rescoring = TRUE,
  epoch_length = 60,
  na_action = c("na", "wake", "zero")
)
```

## Arguments

- counts:

  Numeric vector of activity counts (vertical axis).

- apply_rescoring:

  Apply Webster's rescoring rules (default `TRUE`).

- epoch_length:

  Epoch length in seconds (default 60). The algorithm was validated on
  60-second epochs; other lengths raise a warning.

- na_action:

  How NA-count epochs appear in the output: `"na"` (default) emits `NA`,
  so non-wear gaps are not read as sleep; `"wake"` scores them wake;
  `"zero"` scores them from a zero count.

## Value

Character vector of states, `"S"` (sleep) or `"W"` (wake), the same
length as `counts`.

## Details

The sleep index uses a seven-epoch window (four before, the current
epoch, and two after), with counts divided by 100 and capped at 300:
\$\$D = 0.001 (106 P_4 + 54 P_3 + 58 P_2 + 76 P_1 + 230 C + 74 N_1 + 67
N_2)\$\$ An epoch is scored sleep when \\D \< 1\\. Webster's rescoring
then re-labels short sleep bouts that follow or are surrounded by
sustained wake as wake.

## References

Cole RJ, Kripke DF, Gruen W, Mullaney DJ, Gillin JC (1992). “Automatic
sleep/wake identification from wrist activity.” *Sleep*, **15**(5),
461–469.
[doi:10.1093/sleep/15.5.461](https://doi.org/10.1093/sleep/15.5.461) .

Webster JB, Kripke DF, Messin S, Mullaney DJ, Wyborney G (1982). “An
activity based sleep monitor system for ambulatory use.” *Sleep*,
**5**(4), 389–399.
[doi:10.1093/sleep/5.4.389](https://doi.org/10.1093/sleep/5.4.389) .

## See also

[`sleep.sadeh`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md),
[`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)

## Examples

``` r
# \donttest{
agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
state <- sleep.cole.kripke(agd$axis1)
table(state)
#> state
#>    S    W 
#> 7646 2273 
# }
```
