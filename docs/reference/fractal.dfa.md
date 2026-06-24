# Detrended Fluctuation Analysis (DFA)

Computes the scaling exponent alpha of an activity time series using
Detrended Fluctuation Analysis (Peng et al., 1994). DFA quantifies
long-range temporal correlations: alpha approximately 0.5 indicates
uncorrelated (white) noise, alpha approximately 1.0 indicates 1/f (pink)
noise, and alpha approximately 1.5 indicates Brownian (random-walk /
brown) noise. Healthy human activity fluctuations typically show alpha
in the 0.9 to 1.0 range, with reductions reported in aging and
Alzheimer's disease (Hu et al., 2009).

## Usage

``` r
fractal.dfa(x, scale_min = 4, scale_max = NULL, breakpoint_min = 90)
```

## Arguments

- x:

  Numeric vector of activity counts (minute-level recommended).
  Internally analyzed on the longest continuous non-NA segment.

- scale_min:

  Integer. Smallest window size (box length) in samples. Must be \>= 4
  so that a line can be detrended with residual degrees of freedom.
  Default 4.

- scale_max:

  Integer or NULL. Largest window size in samples. If NULL (default) it
  is set to floor(N / 4) where N is the length of the analyzed segment,
  ensuring at least four windows at the largest scale.

- breakpoint_min:

  Numeric. Window-size boundary (in samples / minutes for minute-level
  data) separating the short-timescale exponent `alpha1` (scales \<
  breakpoint_min) from the long-timescale exponent `alpha2` (scales \>=
  breakpoint_min). Default 90.

## Value

A list with class `"actiRhythm_dfa"` containing:

- alpha:

  Overall scaling exponent: slope of `lm(log10(F) ~ log10(n))` across
  all scales.

- alpha1:

  Short-timescale exponent (scales \< `breakpoint_min`). NA if fewer
  than two qualifying scales.

- alpha2:

  Long-timescale exponent (scales \>= `breakpoint_min`). NA if fewer
  than two qualifying scales.

- scales:

  Integer vector of window sizes n that were used.

- fluctuations:

  Numeric vector of fluctuation magnitudes F(n) corresponding to
  `scales`.

- n_used:

  Length of the analyzed (longest non-NA) segment.

- breakpoint_min:

  The breakpoint value used.

On an unusable series (too short, all NA, or constant) the numeric
scaling outputs are returned as NA with the same structure (never an
error).

## Details

Algorithm (integrated, non-overlapping, linear DFA):

1.  Extract the longest continuous non-NA segment of `x`.

2.  Integrate the mean-centred signal: `y = cumsum(x - mean(x))`.

3.  For each window size n (log-spaced from `scale_min` to `scale_max`),
    split y into floor(N / n) non-overlapping windows of length n, fit
    and remove a least-squares line within each window, and pool the
    residuals. The fluctuation is `F(n) = sqrt(mean(residuals^2))` over
    all pooled residuals.

4.  The scaling exponent is the slope of `log10(F(n))` regressed on
    `log10(n)`.

Window sizes are unique integers chosen on a base-10 log grid, which
gives approximately even spacing in the log-log fit and matches the
convention used by reference implementations such as nonlinearTseries.

## See also

[`multiscale.entropy`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md),
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)

## Examples

``` r
# \donttest{
set.seed(1)
# White noise -> alpha near 0.5
fractal.dfa(rnorm(10000))$alpha
#> [1] 0.4953119
# Brown noise -> alpha near 1.5
fractal.dfa(cumsum(rnorm(10000)))$alpha
#> [1] 1.494906
# }
```
