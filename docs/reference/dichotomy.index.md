# Dichotomy Index (I \< O)

The fraction of in-rest activity counts that fall below the median of
the out-of-rest (active) counts, a rest/active separation index used in
cancer chronobiology (Mormont et al. 2000). A high I\<O means rest is
quiet relative to the active day, marking a well-separated rhythm.

## Usage

``` r
dichotomy.index(counts, timestamps, rest)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- rest:

  A logical vector (TRUE = rest / in-bed), or a character vector of
  states where `"R"`/`"S"`/`"sleep"` mark rest. Same length as `counts`.

## Value

An object of class `actiRhythm_dichotomy`: the index `IO` (percent), the
active-period median, and epoch counts. Never errors.

## References

Mormont MC, Waterhouse J, Bleuzen P, et al. (2000). Marked 24-h
rest/activity rhythms are associated with better quality of life, better
response, and longer survival in patients with metastatic colorectal
cancer. *Clinical Cancer Research*, 6(8):3038-3045.

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 2 * 1440)
h  <- as.numeric(format(ts, "%H"))
counts <- ifelse(h >= 23 | h < 7, 5, 300)
dichotomy.index(counts, ts, rest = h >= 23 | h < 7)
#> Dichotomy Index (I<O)
#> 
#>   I<O:             100.0%
#>   Active median:   300.0 counts
#>   Rest / active epochs: 960 / 1920
#> 
```
