# Empirical Mode Decomposition of the Activity Rhythm

Decomposes the activity series into intrinsic mode functions by
empirical mode decomposition, a data-adaptive, nonlinear alternative to
the linear SSA decomposition (Huang et al. 1998). The intrinsic mode
function whose period is nearest 24 hours is taken as the circadian
component. Optional ensemble EMD adds noise to reduce mode mixing (Wu
and Huang 2009).

## Usage

``` r
circadian.emd(
  counts,
  timestamps,
  max_imf = 10L,
  ensemble = 1L,
  noise_sd = 0.2,
  period_range = c(20, 28),
  epoch_length = 60,
  seed = NULL
)
```

## Arguments

- counts:

  Numeric activity vector (a coarse epoch is recommended for speed).

- timestamps:

  POSIXct timestamps, one per value.

- max_imf:

  Maximum number of modes to extract (default 10).

- ensemble:

  Ensemble size for EEMD (default 1 = plain EMD).

- noise_sd:

  Added-noise SD as a fraction of the series SD (default 0.2).

- period_range:

  Period window (hours) for the circadian mode (default 20 to 28).

- epoch_length:

  Epoch length in seconds (default 60).

- seed:

  Optional seed for the EEMD noise.

## Value

An object of class `actiRhythm_emd`: the IMF matrix, the residual trend,
per-IMF period and variance share, the circadian IMF index, and the
reconstruction error. Never errors.

## References

Huang NE, Shen Z, Long SR, Wu MC, Shih HH, Zheng Q, Yen NC, Tung CC, Liu
HH (1998). “The empirical mode decomposition and the Hilbert spectrum
for nonlinear and non-stationary time series analysis.” *Proceedings of
the Royal Society A*, **454**(1971), 903–995.
[doi:10.1098/rspa.1998.0193](https://doi.org/10.1098/rspa.1998.0193) .

Wu Z, Huang NE (2009). “Ensemble empirical mode decomposition: a
noise-assisted data analysis method.” *Advances in Adaptive Data
Analysis*, **1**(1), 1–41.
[doi:10.1142/S1793536909000047](https://doi.org/10.1142/S1793536909000047)
.

## See also

[`hilbert.huang`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md),
[`circadian.ssa`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
circadian.emd(100 + 60 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#> Empirical Mode Decomposition
#> 
#>   IMFs: 2   reconstruction error: 7.11e-15
#>  IMF period_h var_share
#>    1    22.15     0.932
#>    2    48.00     0.068
#> 
#>   Circadian IMF: 1 (period 22.15 h)
#> 
```
