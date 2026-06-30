# Rest-Active Transition Probabilities

Maximum-likelihood and Bayesian estimates of the rest-to-active and
active-to-rest transition probabilities (Danilevicz et al. 2024): the
transitions out of a state divided by its epochs at risk.

## Usage

``` r
transition.probability(counts, threshold = 1, eps = 0.5)
```

## Arguments

- counts:

  Numeric activity vector; `NA` are dropped.

- threshold:

  Counts at or above which an epoch is active (default 1).

- eps:

  Bayesian Beta pseudo-count added to the transition count and the
  epochs at risk (default 0.5).

## Value

A list with the maximum-likelihood and Bayesian `tp_ra` (rest to active)
and `tp_ar` (active to rest), the active bout count, and the mean active
bout length.

## References

Danilevicz IM, van Hees VT, van der Heide F, Jacob L, Landre B,
Benadjaoud MA, Sabia S (2024). “Measures of fragmentation of rest
activity patterns: mathematical properties and interpretability.” *BMC
Medical Research Methodology*, **24**, 132.
[doi:10.1186/s12874-024-02255-w](https://doi.org/10.1186/s12874-024-02255-w)
.

## See also

[`state.transitions`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md),
[`activity.balance.index`](https://rdazadda.github.io/actiRhythm/reference/activity.balance.index.md)

## Examples

``` r
counts <- c(rep(0, 50), rep(100, 20), rep(0, 40), rep(80, 30), rep(0, 60))
transition.probability(counts)
#> $tp_ar_mle
#> [1] 0.04
#> 
#> $tp_ra_mle
#> [1] 0.01342282
#> 
#> $tp_ar_bayes
#> [1] 0.04950495
#> 
#> $tp_ra_bayes
#> [1] 0.01672241
#> 
#> $n_active_bouts
#> [1] 2
#> 
#> $mean_active_bout
#> [1] 25
#> 
```
