# Rest-Activity State Transition Rates (kRA, kAR)

Computes the rest-to-activity and activity-to-rest transition rates from
a binarized activity series, following the survival-curve method used by
pyActigraphy (Lim et al. 2011). Thresholds the series into rest/active
epochs, builds the per-lag transition probability (hazard) of each bout
type from the bout-length survival curve, and takes a single rate over
the LOWESS "sustained" plateau of that curve.

## Usage

``` r
state.transitions(counts, threshold = 1, frac = 0.3, iter = 0)
```

## Arguments

- counts:

  Numeric activity vector; `NA` are dropped.

- threshold:

  Activity level at or above which an epoch is "active" (default 1, i.e.
  any non-zero count).

- frac:

  LOWESS smoother span for the sustained-region search (default 0.3).

- iter:

  LOWESS robustifying iterations (default 0).

## Value

An object of class `actiRhythm_transitions`: a list with `kRA`/`kAR`
(sustained rest-to-active / active-to-rest rates), `pRA`/`pAR` (overall
per-epoch transition probabilities, 1 / mean bout length), bout counts,
and the two transition curves.

## References

Lim ASP, Yu L, Costa MD, Leurgans SE, Buchman AS, Bennett DA, Saper CB
(2011). “Quantification of the fragmentation of rest-activity patterns
in elderly individuals using a state transition analysis.” *Sleep*,
**34**(11), 1569–1581.
[doi:10.5665/sleep.1400](https://doi.org/10.5665/sleep.1400) .

## Examples

``` r
set.seed(1)
counts <- as.integer(stats::runif(5000) < 0.1) * 100
state.transitions(counts)
#> Rest-Activity State Transitions
#> 
#>   Threshold:  >= 1 counts = active
#>   kRA (rest->active): 0.1051   (471 rest bouts)
#>   kAR (active->rest): 0.8747   (470 active bouts)
#>   pRA / pAR:          0.1055 / 0.8752
```
