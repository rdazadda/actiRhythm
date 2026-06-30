# Multiscale Intradaily Variability (IVm)

Intradaily variability computed across a set of bin sizes and averaged,
the counterpart to multiscale interdaily stability
([`circadian.is.multiscale`](https://rdazadda.github.io/actiRhythm/reference/circadian.is.multiscale.md)).
Fragmentation that is invisible at the hourly scale often shows at finer
bins, so the averaged IVm varies less across recordings than the single
hourly IV (Goncalves et al. 2014).

## Usage

``` r
intradaily.variability.multiscale(counts, timestamps, bin_minutes = 1:60)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- bin_minutes:

  Bin sizes in minutes (default 1 to 60, per Goncalves et al. 2014).

## Value

An object of class `actiRhythm_ivm`: a per-bin `IV` table and the
averaged `IVm`. Never errors; returns `NA` on insufficient data.

## References

Goncalves BSB, Cavalcanti PRA, Tavares GR, Campos TF, Araujo JF (2014).
“Nonparametric methods in actigraphy: an update.” *Sleep Science*,
**7**(3), 158–164.
[doi:10.1016/j.slsci.2014.09.013](https://doi.org/10.1016/j.slsci.2014.09.013)
.

## See also

[`circadian.is.multiscale`](https://rdazadda.github.io/actiRhythm/reference/circadian.is.multiscale.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 3 * 1440)
h  <- as.numeric(format(ts, "%H"))
intradaily.variability.multiscale(ifelse(h >= 23 | h < 7, 5, 300), ts)
#> Multiscale Intradaily Variability
#> 
#>   IVm (averaged): 0.137
#> 
#>  bin_minutes          IV
#>            1 0.006251447
#>            2 0.012505790
#>            3 0.018763030
#>            4 0.025023170
#>            5 0.031286211
#>            6 0.037552156
#>            7 0.030231324
#>            8 0.025363425
#>            9 0.031711636
#>           10 0.062645012
#>           11 0.048476586
#>           12 0.075208914
#>           13 0.059410690
#>           14 0.058635458
#>           15 0.094076655
#>           16 0.063931049
#>           17 0.070327503
#>           18 0.064370776
#>           19 0.078840244
#>           20 0.125581395
#>           21 0.093998547
#>           22 0.097726887
#>           23 0.100205157
#>           24 0.078357397
#>           25 0.109348619
#>           26 0.116119329
#>           27 0.114449213
#>           28 0.125546284
#>           29 0.122368742
#>           30 0.188811189
#>           31 0.131938583
#>           32 0.160936531
#>           33 0.149314613
#>           34 0.146035800
#>           35 0.149708487
#>           36 0.132684653
#>           37 0.155550423
#>           38 0.187181810
#>           39 0.150143505
#>           40 0.134579439
#>           41 0.165639191
#>           42 0.189378057
#>           43 0.191448217
#>           44 0.198615186
#>           45 0.168421053
#>           46 0.202910241
#>           47 0.213954877
#>           48 0.200907806
#>           49 0.217924586
#>           50 0.215230312
#>           51 0.210418259
#>           52 0.227165187
#>           53 0.255578212
#>           54 0.237481797
#>           55 0.252917248
#>           56 0.244180296
#>           57 0.235365600
#>           58 0.258376506
#>           59 0.262059810
#>           60 0.380281690
#> 
```
