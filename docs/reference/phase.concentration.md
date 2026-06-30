# Phase Concentration Tests

Tests whether a set of daily phase markers (acrophases, onsets, L5/M10
times) are concentrated rather than scattered around the clock. Reports
the mean resultant vector, the Rayleigh test of uniformity (Fisher
1993), and the Hermans-Rasson test (Landler et al. 2019), which catches
multimodal clustering that Rayleigh misses.

## Usage

``` r
phase.concentration(times_h, period = 24, n_perm = 2000)
```

## Arguments

- times_h:

  Numeric vector of clock times (hours, 0-24), one per day.

- period:

  Period the times wrap on, in hours (default 24).

- n_perm:

  Permutations for the Hermans-Rasson p-value (default 2000, a speed
  tradeoff; Landler et al. 2019 use 9999 for a finer minimum p-value).

## Value

An object of class `actiRhythm_phasetest`: mean direction, mean
resultant length R, and the Rayleigh and Hermans-Rasson statistics and
p-values. Never errors.

## References

Fisher NI (1993). *Statistical Analysis of Circular Data*. Cambridge
University Press.
[doi:10.1017/CBO9780511564345](https://doi.org/10.1017/CBO9780511564345)
.

Landler L, Ruxton GD, Malkemper EP (2019). “The Hermans-Rasson test as a
powerful alternative to the Rayleigh test for circular statistics in
biology.” *BMC Ecology*, **19**, 30.
[doi:10.1186/s12898-019-0246-8](https://doi.org/10.1186/s12898-019-0246-8)
.

## Examples

``` r
set.seed(1)
onsets <- 23 + stats::rnorm(10, 0, 0.5)      # tightly clustered near 23:00
phase.concentration(onsets %% 24)
#> Phase Concentration Tests
#> 
#>   n days:          10
#>   Mean direction:  23.07 h    R: 0.995
#>   Rayleigh:        Z = 9.91, p = 0.0000
#>   Hermans-Rasson:  V = 27.22, p = 0.0005
#> 
```
