# Composite Phase Deviation (Fischer & Roenneberg, 2016)

Combines, for each day, the deviation of the phase marker from the
individual's own mean phase (precision) and from a reference phase
(accuracy) as `CPD = mean(sqrt(precision^2 + accuracy^2))`. When no
external reference phase is supplied, accuracy is taken relative to the
individual's own mean phase (accuracy = 0), so CPD reduces to the mean
absolute phase deviation - still a valid, published measure of phase
instability, and distinct from the circular SD reported by
`onset_timing_variability`.

## Usage

``` r
composite.phase.deviation(onset_hours, reference_phase = NULL)
```

## Arguments

- onset_hours:

  Numeric vector of daily phase-marker onset times in decimal hours
  (e.g. daily L5 onsets).

- reference_phase:

  Optional reference phase in decimal hours (e.g. a scheduled/expected
  time, or a group mean). Default `NULL` uses the individual's own mean
  phase (accuracy term = 0).

## Value

List with `CPD`, `precision` (mean absolute deviation from own mean
phase, hours), `accuracy` (mean absolute deviation from the reference,
hours), `reference_phase`, and `n_days`.

## References

Fischer D, Vetter C, Roenneberg T (2016). “A novel method to visualise
and quantify circadian misalignment.” *Scientific Reports*, **6**,
38601. [doi:10.1038/srep38601](https://doi.org/10.1038/srep38601) .
