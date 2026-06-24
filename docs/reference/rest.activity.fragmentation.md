# Rest-Activity Bout Fragmentation

Summarizes how broken up the rest-activity rhythm is, from a per-epoch
rest/active state: the mean and median rest and active bout durations,
the number of state transitions, and transitions per day. These add a
bout-length view of fragmentation to the transition probabilities of
[`state.transitions`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
(kRA/kAR) (Lim et al. 2011). This covers rest-activity-rhythm
fragmentation only; it omits the sedentary-behaviour bout distribution
(Gini, power law, hazard), which is a physical-activity-epidemiology
concern, not a circadian one.

## Usage

``` r
rest.activity.fragmentation(state, timestamps, epoch_length = 60)
```

## Arguments

- state:

  Per-epoch state: a logical vector (TRUE = active) or a character
  vector where `"R"`/`"S"`/`"sleep"`/`"rest"` mark rest.

- timestamps:

  POSIXct timestamps, one per value.

- epoch_length:

  Epoch length in seconds (default 60).

## Value

An object of class `actiRhythm_rafrag`: mean/median rest and active bout
durations (minutes), bout counts, transition count, and transitions per
day. Never errors.

## References

Lim ASP, Yu L, Costa MD, et al. (2011). Quantification of the
fragmentation of rest-activity patterns in elderly individuals using a
state transition analysis. *Sleep*, 34(11):1569-1581.
[doi:10.5665/sleep.1400](https://doi.org/10.5665/sleep.1400)

## See also

[`state.transitions`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h  <- as.numeric(format(ts, "%H"))
rest.activity.fragmentation(h >= 7 & h < 23, ts)
#> Rest-Activity Bout Fragmentation
#> 
#>   Active bouts: 3, mean 960 min (median 960)
#>   Rest bouts:   4, mean 360 min (median 450)
#>   Transitions:  6 (2.0 per day)
#> 
```
