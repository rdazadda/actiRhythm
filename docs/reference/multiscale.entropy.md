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
#>  [1] 2.4475282 2.1095652 1.9424928 1.7839571 1.6271744 1.5529452 1.5719688
#>  [8] 1.4297197 1.3623616 1.3134330 1.3228068 1.4076692 1.2764495 1.2212146
#> [15] 1.1643755 1.1213814 1.2828048 1.1881730 1.0746168 0.9322216
# 1/f-like (random walk): flatter / sustained entropy
multiscale.entropy(cumsum(rnorm(2000)))$mse
#>  [1] 0.1155946 0.1613639 0.2117052 0.2512892 0.2873474 0.3162274 0.3391400
#>  [8] 0.3552696 0.3638887 0.3696590 0.4063044 0.3947861 0.4363337 0.4262172
#> [15] 0.4326016 0.4276882 0.4768064 0.4615090 0.5163658 0.4473874
# }
```
