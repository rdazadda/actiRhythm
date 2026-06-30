# Multiscale Entropy Curve

Sample entropy from
[`multiscale.entropy`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
against the coarse-graining scale, with the complexity index (area)
shaded and the across-scale slope drawn. A curve that holds its entropy
across scales marks a complex signal; one that falls is closer to noise.
Returns a `ggplot` object and never errors.

## Usage

``` r
plot_mse(counts, scales = 1:20, m = 2L, r = 0.15)
```

## Arguments

- counts:

  Numeric activity vector.

- scales:

  Integer coarse-graining scales (default `1:20`).

- m:

  Embedding dimension for sample entropy (default 2).

- r:

  Tolerance as a fraction of the series SD (default 0.15).

## Value

A `ggplot` object.

## References

Costa M, Goldberger AL, Peng C (2002). “Multiscale entropy analysis of
complex physiologic time series.” *Physical Review Letters*, **89**(6),
068102.
[doi:10.1103/PhysRevLett.89.068102](https://doi.org/10.1103/PhysRevLett.89.068102)
.

## Examples

``` r
# \donttest{
set.seed(1)
plot_mse(rnorm(3000))

# }
```
