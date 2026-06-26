# State-Space (Hidden Markov) Rest-Activity Model

Fits an unsupervised Gaussian hidden Markov model to the activity
counts, inferring latent rest and active states and the rhythm with
which the subject moves between them (Huang et al. 2018). A threshold
scorer labels each epoch independently; the HMM uses the persistence of
states. Its decoded path gives a 24-hour state-occupancy profile, the
probability of being at rest across the day.

## Usage

``` r
rest.hmm(
  counts,
  timestamps,
  states = 2L,
  transform = c("sqrt", "log", "none"),
  max_iter = 200L,
  tol = 1e-06,
  n_starts = 5L,
  seed = NULL
)
```

## Arguments

- counts:

  Numeric activity vector (a coarse epoch is recommended for speed).

- timestamps:

  POSIXct timestamps, one per value.

- states:

  Number of hidden states (default 2 = rest/active; 3 adds a moderate
  state).

- transform:

  Variance-stabilizing transform of the counts: `"sqrt"` (default),
  `"log"`, or `"none"`.

- max_iter, tol:

  EM iteration cap and log-likelihood tolerance.

- n_starts:

  Random restarts; the highest-likelihood fit is kept (default 5).

- seed:

  Optional seed for the restarts.

## Value

An object of class `actiRhythm_hmm`: the per-state emission means and
SDs, the transition matrix, the decoded `state_path` and a `sleep_state`
("S"/"W") vector, a 24-hour state-occupancy profile, and the
log-likelihood with AIC/BIC. Never errors.

## References

Huang Q, Cohen D, Komarzynski S, Li X, Innominato P, Levi F, Finkenstadt
B (2018). “Hidden Markov models for monitoring circadian rhythmicity in
telemetric activity data.” *Journal of the Royal Society Interface*,
**15**(139), 20170885.
[doi:10.1098/rsif.2017.0885](https://doi.org/10.1098/rsif.2017.0885) .

## See also

[`sleep.changepoints`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md),
[`sleep.cole.kripke`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
h  <- as.numeric(format(ts, "%H"))
counts <- ifelse(h >= 23 | h < 7, 2, 250) + pmax(0, stats::rnorm(length(ts), 0, 10))
rest.hmm(counts, ts)
#> State-Space Rest-Activity Model (Gaussian HMM)
#> 
#>   States: 2   log-likelihood: -349.0   AIC: 710.1
#> 
#>  state mean_transformed sd_transformed  label
#>      1             2.26           1.07   rest
#>      2            15.94           0.19 active
#> 
#>   Time at rest: 33%   persistence (rest self-transition): 0.98
#> 
```
