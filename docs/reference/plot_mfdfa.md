# Multifractal DFA Spectrum

Two panels from
[`mfdfa`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md): the
generalized Hurst exponent h(q) against q (flat = monofractal,
decreasing = multifractal, with h(2) marked) and the singularity
spectrum f(alpha) (the arch whose width measures multifractality).
Returns a `ggplot` object and never errors.

## Usage

``` r
plot_mfdfa(counts, q_values = seq(-5, 5, by = 0.5), detrend_order = 1L)
```

## Arguments

- counts:

  Numeric activity vector.

- q_values:

  Moment orders to evaluate (default `seq(-5, 5, by = 0.5)`).

- detrend_order:

  Within-window polynomial detrend order (default 1).

## Value

A `ggplot` object.

## References

Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, Havlin S, Bunde A,
Stanley HE (2002). “Multifractal detrended fluctuation analysis of
nonstationary time series.” *Physica A: Statistical Mechanics and its
Applications*, **316**(1-4), 87–114.
[doi:10.1016/S0378-4371(02)01383-3](https://doi.org/10.1016/S0378-4371%2802%2901383-3)
.

## Examples

``` r
# \donttest{
set.seed(1)
plot_mfdfa(cumsum(rnorm(8000)))

# }
```
