# Multiscale Sample Entropy (MSE)

Computes Multiscale Sample Entropy (Costa et al., 2002, 2005), which
applies Sample Entropy (Richman & Moorman, 2000) to coarse-grained
versions of an activity time series across a range of temporal scales.
MSE distinguishes genuinely complex signals (whose entropy is sustained
or increases across scales, e.g. 1/f-like physiological signals) from
uncorrelated random signals (whose entropy falls monotonically with
scale).

## Usage

``` r
multiscale.entropy(x, scales = 1:20, m = 2, r = 0.15)
```

## Arguments

- x:

  Numeric vector of activity counts (minute-level recommended).
  Internally analyzed on the longest continuous non-NA segment.

- scales:

  Integer vector of coarse-graining scale factors tau. Default `1:20`.

- m:

  Integer embedding dimension (template length) for Sample Entropy.
  Default 2.

- r:

  Numeric tolerance as a fraction of the standard deviation of the
  ORIGINAL (scale-1) series. The Chebyshev matching radius used at every
  scale is `r * sd(x_original)`, following the standard MSE convention
  of holding r fixed in absolute units across scales. Default 0.15.

## Value

A list with class `"actiRhythm_mse"` containing:

- mse:

  Numeric vector of SampEn values, one per requested scale (NA where the
  coarse-grained series is too short or yields no matches).

- scales:

  Integer vector of the scale factors used.

- area:

  Sum of `mse` over scales (complexity index), `sum(mse, na.rm = TRUE)`.

- slope:

  Slope of `lm(mse ~ scales)` over the non-NA points (negative =\>
  entropy declines with scale, typical of noise).

- r_absolute:

  The absolute tolerance `r * sd(x_original)` used.

- n_used:

  Length of the analyzed (longest non-NA) segment.

On an unusable series the entropy vector is all NA with the same
structure (never an error).

## Details

For each scale tau the series is coarse-grained into non-overlapping
means of length tau, then Sample Entropy is computed on the
coarse-grained series with embedding dimension `m` and the FIXED
absolute tolerance `r * sd(x_original)`. Sample Entropy is the negative
natural log of the conditional probability that two sub-sequences
matching for `m` points (within the tolerance, Chebyshev distance,
self-matches excluded) also match for `m + 1` points. The implementation
is fully self-contained (base R / stats only).

## See also

[`fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md),
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)

## Examples

``` r
# \donttest{
set.seed(1)
# White noise: entropy decreases with scale
multiscale.entropy(rnorm(2000))$mse
#>  [1] 2.446324 2.105090 1.935192 1.777846 1.621644 1.552057 1.569350 1.427566
#>  [9] 1.357506 1.303940 1.305857 1.384235 1.256591 1.199179 1.138336 1.120196
#> [17] 1.254367 1.146115 1.051494 0.892656
# 1/f-like (random walk): flatter / sustained entropy
multiscale.entropy(cumsum(rnorm(2000)))$mse
#>  [1] 0.1147112 0.1594668 0.2091053 0.2474783 0.2828905 0.3114901 0.3333711
#>  [8] 0.3486852 0.3554676 0.3622740 0.4004146 0.3876029 0.4276109 0.4148211
#> [15] 0.4246965 0.4166384 0.4553923 0.4447822 0.5024574 0.4326161
# }
```
