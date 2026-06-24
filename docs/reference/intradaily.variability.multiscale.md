# Multiscale Intradaily Variability (IVm)

Intradaily variability computed across a set of bin sizes and averaged,
the counterpart to multiscale interdaily stability
([`circadian.is.multiscale`](https://rdazadda.github.io/actiRhythm/reference/circadian.is.multiscale.md)).
Fragmentation that is invisible at the hourly scale often shows at finer
bins, so the averaged IVm varies less across recordings than the single
hourly IV (Goncalves et al. 2014).

## Usage

``` r
intradaily.variability.multiscale(
  counts,
  timestamps,
  bin_minutes = c(5, 10, 15, 30, 60)
)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- bin_minutes:

  Bin sizes in minutes (default 5 to 60).

## Value

An object of class `actiRhythm_ivm`: a per-bin `IV` table and the
averaged `IVm`. Never errors; returns `NA` on insufficient data.

## References

Goncalves BSB, Cavalcanti PRA, Tavares GR, Campos TF, Araujo JF (2014).
Nonparametric methods in actigraphy: an update. *Sleep Science*,
7(3):158-164.
[doi:10.1016/j.slsci.2014.09.013](https://doi.org/10.1016/j.slsci.2014.09.013)

## See also

[`circadian.is.multiscale`](https://rdazadda.github.io/actiRhythm/reference/circadian.is.multiscale.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h  <- as.numeric(format(ts, "%H"))
intradaily.variability.multiscale(ifelse(h >= 23 | h < 7, 5, 300), ts)
#> Multiscale Intradaily Variability
#> 
#>   IVm (averaged): 0.151
#> 
#>  bin_minutes         IV
#>            5 0.03128621
#>           10 0.06264501
#>           15 0.09407666
#>           30 0.18881119
#>           60 0.38028169
#> 
```
