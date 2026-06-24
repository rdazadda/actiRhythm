# Activity Balance Index

The Activity Balance Index (Danilevicz et al. 2024), a 0 to 1 transform
of a detrended fluctuation analysis scaling exponent that peaks at the
healthy \\\alpha = 1\\ (1/f) balance: \\ABI(\alpha) = \exp(-\|\alpha -
1\|\\e^{-2})\\.

## Usage

``` r
activity.balance.index(x)
```

## Arguments

- x:

  Either a numeric scaling exponent, or a fractal object with `alpha`
  (and optionally `alpha1`, `alpha2`) such as the result of the
  package's detrended fluctuation analysis.

## Value

If `x` is numeric, the scalar ABI. If `x` is a fractal object, a list
with `ABI_overall`, `ABI_short`, `ABI_long`.

## References

Danilevicz IM, et al. (2024). Measures of fragmentation of rest activity
patterns: mathematical properties and interpretability. *BMC Medical
Research Methodology*, 24:132.
[doi:10.1186/s12874-024-02255-w](https://doi.org/10.1186/s12874-024-02255-w)

## Examples

``` r
activity.balance.index(1.0)   # perfect 1/f balance -> 1
#> [1] 1
activity.balance.index(0.7)
#> [1] 0.9602126
```
