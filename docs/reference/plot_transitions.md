# Rest-Activity Transition Curves

The rest-to-active and active-to-rest transition hazards against bout
length, with a LOWESS fit and the sustained rates kRA and kAR marked.
Returns a `ggplot` object and never errors.

## Usage

``` r
plot_transitions(counts, threshold = 1, frac = 0.3)
```

## Arguments

- counts:

  Numeric activity vector.

- threshold:

  Counts at or above which an epoch is active (default 1).

- frac:

  LOWESS span (default 0.3).

## Value

A `ggplot` object.

## References

Lim ASP, Yu L, Costa MD, Leurgans SE, Buchman AS, Bennett DA, Saper CB
(2011). “Quantification of the fragmentation of rest-activity patterns
in elderly individuals using a state transition analysis.” *Sleep*,
**34**(11), 1569–1581.
[doi:10.5665/sleep.1400](https://doi.org/10.5665/sleep.1400) .

## Examples

``` r
set.seed(1)
plot_transitions(as.integer(stats::runif(8000) < 0.1) * 100)

```
