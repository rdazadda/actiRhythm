# Hilbert-Huang Instantaneous Phase and Frequency

Computes the instantaneous amplitude, phase, and period of a circadian
intrinsic mode function via its analytic signal (Huang et al. 1998).
Where the cosinor gives one acrophase for the whole recording, the
instantaneous phase tracks the rhythm cycle by cycle, and the spread of
the instantaneous period measures how stationary the circadian band is.

## Usage

``` r
hilbert.huang(x, timestamps = NULL, imf = NULL, epoch_length = 60)
```

## Arguments

- x:

  A `actiRhythm_emd` object, or a numeric activity vector (then
  `timestamps` is required).

- timestamps:

  POSIXct timestamps (required when `x` is numeric).

- imf:

  Which IMF to analyse; default is the circadian IMF of the EMD object.

- epoch_length:

  Epoch length in seconds (default 60; used when `x` is numeric).

## Value

An object of class `actiRhythm_hht`: the instantaneous amplitude, phase,
and period series, plus the mean instantaneous period and its SD, mean
amplitude and its CV, and the fraction of time the period stays in 20-28
h. Never errors.

## References

Huang NE, Shen Z, Long SR, et al. (1998). The empirical mode
decomposition and the Hilbert spectrum. *Proceedings of the Royal
Society A*, 454:903-995.
[doi:10.1098/rspa.1998.0193](https://doi.org/10.1098/rspa.1998.0193)

## See also

[`circadian.emd`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
e <- circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
hilbert.huang(e)
#> Hilbert-Huang Instantaneous Dynamics
#> 
#>   Instantaneous period: 24.85 h (SD 3.60)
#>   Instantaneous amplitude: 58.8 (CV 0.25)
#>   Time in 20-28 h band: 76%
#> 
```
