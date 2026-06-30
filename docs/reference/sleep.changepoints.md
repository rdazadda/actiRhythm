# Change-Point Detection of Sleep and Wake Onsets

Locates the sleep-onset and wake-onset time of each circadian cycle with
a cosinor-anchored mean-shift change point. A fixed 24-hour cosinor
bounds each rest and active span roughly (the cosinor anchoring follows
CircaCP, Chen and Sun 2024), and the precise transition inside each
bound is then placed with a single least-squares mean-shift change point
on the raw counts. The result is a per-night sleep-onset / wake-onset
table. A single rest-activity transition rate (such as
[`state.transitions`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md))
cannot localise this timing.

## Usage

``` r
sleep.changepoints(
  counts,
  timestamps,
  period = 1440,
  thr = 0.2,
  window_minutes = 240
)
```

## Arguments

- counts:

  Numeric activity vector (minute epochs recommended).

- timestamps:

  POSIXct timestamps, one per value.

- period:

  Cosinor period in minutes (default 1440, one day).

- thr:

  Dichotomisation threshold on the range-scaled cosine, in \\\[0, 1\]\\
  (default 0.2, approximating CircaCP's lower-20\\ fitted curve above
  `thr` is the rough active span, below it the rough rest span.

- window_minutes:

  Half-width of the search window, in minutes, in which each rough
  boundary is refined to a change point (default 240).

## Value

An object of class `actiRhythm_changepoints`: the cosinor summary, a
`changepoints` data frame (time and type, "sleep onset" or "wake
onset"), a `sleep_episodes` data frame (sleep onset, wake onset, and
duration in hours), and the mean sleep duration. The function never
errors; on insufficient data it returns the structure with
`insufficient = TRUE`.

## References

Chen S, Sun X (2024). “Validating CircaCP: a generic sleep-wake cycle
detection algorithm for unlabelled actigraphy data.” *Royal Society Open
Science*, **11**(5), 231468.
[doi:10.1098/rsos.231468](https://doi.org/10.1098/rsos.231468) .

## See also

[`state.transitions`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md),
[`sleep.regularity.index`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)

## Examples

``` r
# Five days with a clear active day / restful night
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 5 * 1440)
h  <- as.numeric(format(ts, "%H"))
counts <- ifelse(h >= 8 & h < 23, 300, 5) + pmax(0, stats::rnorm(length(ts), 0, 5))
sleep.changepoints(counts, ts)
#> Change-Point Sleep/Wake Detection
#> 
#>   Span:           5.0 days (7200 epochs)
#>   Cosinor acrophase: 15.5 h
#>   Change points:  10 (4 sleep episodes)
#>   Mean sleep duration: 9.0 h
#> 
#>   First sleep episodes:
#>     sleep 01-01 22:59  ->  wake 01-02 07:59  (9.0 h)
#>     sleep 01-02 22:59  ->  wake 01-03 07:59  (9.0 h)
#>     sleep 01-03 22:59  ->  wake 01-04 07:59  (9.0 h)
#> 
#>   Reference: Chen and Sun (2024)
#> 
```
