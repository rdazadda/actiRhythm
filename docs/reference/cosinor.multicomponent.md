# Multicomponent Cosinor with Model Selection

Fits the multi-component (harmonic) cosinor of Cornelissen (2014) with
one to several harmonics of the fundamental period and picks the number
of harmonics by an information criterion (AIC or BIC, a package choice),
so it captures a bimodal or asymmetric daily shape without your choosing
the order by hand. The single cosinor is the one-harmonic special case.

## Usage

``` r
cosinor.multicomponent(
  counts,
  timestamps,
  period = 24,
  max_harmonics = 3,
  criterion = c("AIC", "BIC")
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- period:

  Fundamental period in hours (default 24).

- max_harmonics:

  Largest number of harmonics to consider (default 3).

- criterion:

  `"AIC"` (default) or `"BIC"` for model selection.

## Value

An object of class `actiRhythm_multicosinor`: the selected number of
harmonics, the per-harmonic amplitude and acrophase, the MESOR,
R-squared, and the full model-comparison table. Never errors.

## References

Cornelissen G (2014). “Cosinor-based rhythmometry.” *Theoretical Biology
and Medical Modelling*, **11**, 16.
[doi:10.1186/1742-4682-11-16](https://doi.org/10.1186/1742-4682-11-16) .

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
y  <- 100 + 40 * cos(2 * pi * th / 24) + 20 * cos(2 * pi * 2 * th / 24)
cosinor.multicomponent(y, ts)
#> Multicomponent Cosinor
#> 
#>   Selected harmonics: 2 (by AIC)   MESOR: 100.0   R-squared: 1.000
#> 
#>  harmonic amplitude acrophase_h
#>         1     39.89        0.01
#>         2     19.77        0.01
#> 
```
