# Wavelet Cone of Influence

The Morlet cone of influence (Torrence and Compo 1998): the largest
period, in hours, at each time point below which the wavelet power is
free of edge effects. Power at periods above this curve, near the start
and end of the recording, is unreliable.
[`circadian.wavelet`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)
returns this curve in its result; this helper computes it directly for
plotting or masking.

## Usage

``` r
wavelet.coi(timestamps, omega0 = 6, epoch_length = 60)
```

## Arguments

- timestamps:

  POSIXct timestamps, one per epoch.

- omega0:

  Morlet central frequency (default 6).

- epoch_length:

  Epoch length in seconds (default 60).

## Value

Numeric vector, one per epoch, of the maximum reliable period (hours).

## References

Torrence C, Compo GP (1998). A practical guide to wavelet analysis.
*Bulletin of the American Meteorological Society*, 79(1):61-78.
[doi:10.1175/1520-0477(1998)079\<0061:APGTWA\>2.0.CO;2](https://doi.org/10.1175/1520-0477%281998%29079%3C0061%3AAPGTWA%3E2.0.CO%3B2)

## See also

[`circadian.wavelet`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
range(wavelet.coi(ts, epoch_length = 600))
#> [1] 1.217454e-06 5.247225e+01
```
