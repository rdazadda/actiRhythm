# Multifractal Detrended Fluctuation Analysis (MF-DFA)

Generalizes detrended fluctuation analysis
([`fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md))
to a spectrum of moment orders q (Kantelhardt et al. 2002). A flat
generalized Hurst exponent h(q) indicates a monofractal signal; a
decreasing h(q) and a wide multifractal spectrum indicate
multifractality. h(2) equals the standard DFA scaling exponent.

## Usage

``` r
mfdfa(
  x,
  scale_min = 8L,
  scale_max = NULL,
  q_values = seq(-5, 5, by = 0.5),
  both_ends = TRUE
)
```

## Arguments

- x:

  Numeric series (e.g. activity counts); the longest gap-free run is
  analysed.

- scale_min:

  Smallest window size in samples (default 8).

- scale_max:

  Largest window size (default `floor(N/4)`).

- q_values:

  Moment orders to evaluate (default `seq(-5, 5, 0.5)`).

- both_ends:

  If `TRUE` (default) windows are taken from both ends so the tail of
  the profile is not discarded; `FALSE` reproduces the start-only
  convention of
  [`fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md).

## Value

An object of class `actiRhythm_mfdfa`: a list with `q_values`, `h_q`
(generalized Hurst exponent), `tau_q` (mass exponent), `alpha`/`f_alpha`
(the multifractal spectrum), `alpha_dfa` (`= h(2)`) and the spectrum
`width`. Returns an NA structure on insufficient data; never errors.

## References

Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, Havlin S, Bunde A,
Stanley HE (2002). Multifractal detrended fluctuation analysis of
nonstationary time series. Physica A, 316(1-4):87-114.

## See also

[`fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)

## Examples

``` r
set.seed(1)
mfdfa(stats::rnorm(4096))$alpha_dfa   # near 0.5 for white noise
#> [1] 0.4772896
```
