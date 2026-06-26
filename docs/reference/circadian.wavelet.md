# Continuous Wavelet Transform of the Activity Rhythm

Runs a Morlet continuous wavelet transform on the activity series and
returns the time-frequency power surface across circadian and ultradian
periods and a dominant-period-over-time track (Torrence and Compo 1998;
Leise 2013). Unlike the sliding chi-square spectrogram, it localizes
period drift at every time point, so a lengthening or fragmenting rhythm
shows up directly.

## Usage

``` r
circadian.wavelet(counts, timestamps, dj = 1/12, omega0 = 6, epoch_length = 60)
```

## Arguments

- counts:

  Numeric activity vector (a coarse epoch, e.g. 10-minute bins, is
  recommended for speed; see the example).

- timestamps:

  POSIXct timestamps, one per value.

- dj:

  Scale resolution in voices per octave step (default 1/12).

- omega0:

  Morlet central frequency (default 6).

- epoch_length:

  Epoch length in seconds (default 60).

## Value

An object of class `actiRhythm_wavelet`: the period grid (hours), the
time-averaged global power spectrum, the per-time dominant period, the
overall peak period, the power matrix, and the cone of influence
(`coi_period_h`, the largest reliable period at each time, with the
`in_coi` mask of edge-affected cells). The global spectrum and dominant
period are computed outside the cone so edge effects do not bias them. A
`significant` logical matrix flags cells whose power exceeds the 95\\
confidence level against an AR(1) red-noise background (with the
per-scale threshold `sig_power` and the lag-1 autocorrelation `phi`),
which separates a real rhythm from background. Never errors.

## References

Torrence C, Compo GP (1998). “A practical guide to wavelet analysis.”
*Bulletin of the American Meteorological Society*, **79**(1), 61–78.
[doi:10.1175/1520-0477(1998)079\<0061:APGTWA\>2.0.CO;2](https://doi.org/10.1175/1520-0477%281998%29079%3C0061%3AAPGTWA%3E2.0.CO%3B2)
.

Leise TL (2013). “Wavelet analysis of circadian and ultradian behavioral
rhythms.” *Journal of Circadian Rhythms*, **11**, 5.
[doi:10.1186/1740-3391-11-5](https://doi.org/10.1186/1740-3391-11-5) .

## Examples

``` r
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 600, length.out = 6 * 144)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
circadian.wavelet(100 + 80 * cos(2 * pi * th / 24), ts, epoch_length = 600)
#> Continuous Wavelet Transform (Morlet)
#> 
#>   Periods scanned: 0.3 to 47.9 h (87 scales)
#>   Peak global period: 23.95 h
#>   Dominant period: median 23.95 h (IQR 12.68-23.95)
#>   Cone of influence: 81% of cells edge-reliable
#>   Significant power: 7.2% of reliable cells (95% vs AR(1) red noise, phi=0.99)
#> 
```
