# Beyond the basics: profile, decomposition, and finer metrics

``` r

library(actiRhythm)
library(ggplot2)
```

This article continues from the get-started vignette, picking up the
same bundled recording and going past the core workflow. It covers a
fuller functional profile of the day and a singular-spectrum
decomposition of the whole series, then a wider set of methods: finer
nonparametric metrics, the shape and phase of the day, time-frequency
and adaptive decomposition, and registration, residual structure, and
states. Every object below is rebuilt from the installed package so the
article stands on its own.

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))

bin      <- as.integer(as.numeric(agd$timestamp) %/% 600)   # 10-minute bins
counts10 <- as.numeric(tapply(agd$axis1, bin, sum))
t10      <- as.POSIXct(as.numeric(names(tapply(agd$axis1, bin, sum))) * 600,
                       origin = "1970-01-01", tz = "UTC")

state <- sleep.cole.kripke(agd$axis1)
```

## A fuller profile, and a decomposition

The cosinor fits one cosine to the day.
[`circadian.flm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md)
fits a functional model of the same 24-hour profile with several Fourier
harmonics, so it follows the real, asymmetric shape of the day instead
of forcing a symmetric curve ([Wang et al., 2011](#ref-wang2011flm));
the single cosine is its one-harmonic special case.

``` r

flm <- circadian.flm(agd$axis1, agd$timestamp)
flm
#> Functional Linear Model (24-hour activity profile)
#> 
#>   Basis:        fourier (order 4, 9 terms)
#>   Period:       24 hours
#>   Days / profile points: 8 / 24
#> 
#>   Model fit:
#>     R-squared:   0.7352
#>     AIC:         406.84
#>     F-statistic: 5.21 (p = 0.003)
#> 
#>   Dominant harmonic: H1, amplitude 303.40, acrophase 16.99 h
#> 
#>   Peak 956.3 at 19.2 h; trough -117.7 at 3.8 h
#> 
#>   Basis fit follows Wang et al. (2011); F/p are ordinary regression
#>   statistics (not Wang's between-subject permutation test).
```

``` r

ggplot(flm$smooth_curve, aes(t, activity)) +
  geom_point(data = flm$fitted_profile, aes(t, observed), colour = "grey55") +
  geom_line(colour = "#236192", linewidth = 1) +
  labs(x = "Hour of day", y = "Activity") +
  theme_actiRhythm()
```

![The functional fit (line) over the averaged hourly profile (points).
The extra Fourier harmonics let the curve follow the asymmetric morning
rise and evening decline that a single cosine smooths
over.](beyond-the-basics_files/figure-html/flm-plot-1.png)

The functional fit (line) over the averaged hourly profile (points). The
extra Fourier harmonics let the curve follow the asymmetric morning rise
and evening decline that a single cosine smooths over.

The extra harmonics buy back the variance the single cosine left on the
table, and the curve now bends through the morning rise the one cosine
ran straight past.

[`circadian.ssa()`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)
takes the opposite view. Instead of an average day it decomposes the
whole recording with singular spectrum analysis into additive
components: a slow trend, a circadian pair, and noise ([Golyandina &
Zhigljavsky, 2013](#ref-golyandina2013)). The trajectory matrix grows
with the series, so bin a long minute recording to a coarser epoch
first.

``` r

bin      <- as.integer(as.numeric(agd$timestamp) %/% 600)   # 10-minute bins
counts10 <- tapply(agd$axis1, bin, sum)
t10      <- as.POSIXct(as.numeric(names(counts10)) * 600, origin = "1970-01-01", tz = "UTC")
ssa <- circadian.ssa(as.numeric(counts10), t10)
ssa
#> Singular Spectrum Analysis (Basic SSA)
#> 
#>   Window length L:    144 epochs (24.0 h)
#>   Series length n:    993 (K = 850, span 6.9 days)
#>   Components kept:     10 of rank 144
#> 
#>   Variance explained (leading components):
#>     ET1   lambda = 0.2929  (cumulative 0.2929)
#>     ET2   lambda = 0.0814  (cumulative 0.3743)
#>     ET3   lambda = 0.0578  (cumulative 0.4321)
#>     ET4   lambda = 0.0295  (cumulative 0.4616)
#>     ET5   lambda = 0.0293  (cumulative 0.4909)
#> 
#>   Grouping:
#>     Trend:      components 1
#>     Circadian:  components 2, 3 (13.9% of variance)
#>     Fundamental period: 24.49 h
#> 
#>   Reference: Golyandina and Zhigljavsky (2013)
```

``` r

ggplot(data.frame(time = t10, circadian = ssa$circadian), aes(time, circadian)) +
  geom_line(colour = "#236192") +
  labs(x = "Time", y = "Circadian component") +
  theme_actiRhythm()
```

![The circadian component singular spectrum analysis pulls out of the
whole recording: the daily rhythm separated from the slow trend and the
noise.](beyond-the-basics_files/figure-html/ssa-plot-1.png)

The circadian component singular spectrum analysis pulls out of the
whole recording: the daily rhythm separated from the slow trend and the
noise.

For this recording the circadian pair carries about 14 percent of the
variance with a fundamental period near 24.5 hours, the same
slightly-long period the periodogram found earlier, now isolated as a
clean component of the raw series.

## A wider set of methods

The core above answers most questions. actiRhythm also carries a wider
set of circadian methods for finer or more specialized work, shown here
on the same recording.

### Finer nonparametric metrics

[`intradaily.variability.multiscale()`](https://rdazadda.github.io/actiRhythm/reference/intradaily.variability.multiscale.md)
averages intradaily variability across bin sizes,
[`activity.extrema()`](https://rdazadda.github.io/actiRhythm/reference/activity.extrema.md)
generalizes L5 and M10 to any window with onset and midpoint times,
[`dichotomy.index()`](https://rdazadda.github.io/actiRhythm/reference/dichotomy.index.md)
measures how cleanly rest separates from the active day ([Mormont et
al., 2000](#ref-mormont2000)), and
[`circadian.daily()`](https://rdazadda.github.io/actiRhythm/reference/circadian.daily.md)
reports each day on its own so within-recording drift shows.
[`rest.activity.fragmentation()`](https://rdazadda.github.io/actiRhythm/reference/rest.activity.fragmentation.md)
adds the rest and active bout-length view that complements the kRA/kAR
rates above.

``` r

daily <- circadian.daily(agd$axis1, agd$timestamp)
dichotomy.index(agd$axis1, rest = state == "S")
#> Dichotomy Index (I<O)
#> 
#>   I<O:             100.0%
#>   Active median:   651.5 counts
#>   Rest / active epochs: 7641 / 2278
```

``` r

ggplot(daily$daily, aes(date, M10_onset_h)) +
  geom_line(colour = "grey60") + geom_point(colour = "#236192", size = 2) +
  labs(x = NULL, y = "M10 onset (h)") +
  theme_actiRhythm()
```

![The M10 onset (the start of the most active ten hours) for each day.
Reporting it per day, rather than pooled, shows how steady the active
phase is across the
recording.](beyond-the-basics_files/figure-html/daily-plot-1.png)

The M10 onset (the start of the most active ten hours) for each day.
Reporting it per day, rather than pooled, shows how steady the active
phase is across the recording.

### The shape and phase of the day

[`cosinor.multicomponent()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md)
adds harmonics to the single cosine and picks how many by information
criterion,
[`activity.onset.offset()`](https://rdazadda.github.io/actiRhythm/reference/activity.onset.offset.md)
marks the daily activity onset and offset by a relative-difference
contrast, and
[`phase.concentration()`](https://rdazadda.github.io/actiRhythm/reference/phase.concentration.md)
tests whether the daily phase markers cluster (Rayleigh and
Hermans-Rasson).

``` r

cosinor.multicomponent(agd$axis1, agd$timestamp)
#> Multicomponent Cosinor
#> 
#>   Selected harmonics: 3 (by AIC)   MESOR: 330.1   R-squared: 0.640
#> 
#>  harmonic amplitude acrophase_h
#>         1    294.33       20.56
#>         2    173.67       11.03
#>         3     42.81        0.14
phase.concentration(daily$daily$M10_onset_h)
#> Phase Concentration Tests
#> 
#>   n days:          7
#>   Mean direction:  13.86 h    R: 0.679
#>   Rayleigh:        Z = 3.23, p = 0.0332
#>   Hermans-Rasson:  V = 6.21, p = 0.0175
```

### Time-frequency and adaptive decomposition

On the ten-minute series from earlier,
[`circadian.wavelet()`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)
gives the time-frequency power spectrum,
[`ultradian.bandpower()`](https://rdazadda.github.io/actiRhythm/reference/ultradian.bandpower.md)
splits the variance into ultradian bands, and
[`circadian.emd()`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)
with
[`hilbert.huang()`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md)
extracts the circadian component data-adaptively and tracks its
instantaneous period.

``` r

wav <- circadian.wavelet(as.numeric(counts10), t10, epoch_length = 600)
ggplot(data.frame(period = wav$period_hours, power = wav$global_power),
       aes(period, power)) +
  geom_line(colour = "#236192") +
  geom_vline(xintercept = 24, linetype = 2, colour = "grey50") +
  scale_x_continuous(trans = "log2") +
  labs(x = "Period (h, log scale)", y = "Wavelet power (scale-rectified)") +
  theme_actiRhythm()
```

![The scale-rectified global wavelet power spectrum (\|W\|^2/scale, Liu
et al. 2007). The peak sits at the circadian period (dashed line at 24
hours), recovered without assuming a fixed cosine
shape.](beyond-the-basics_files/figure-html/wavelet-plot-1.png)

The scale-rectified global wavelet power spectrum (\|W\|^2/scale, Liu et
al. 2007). The peak sits at the circadian period (dashed line at 24
hours), recovered without assuming a fixed cosine shape.

``` r

plot_wavelet(as.numeric(counts10), t10, epoch_length = 600)
```

![The full wavelet power scalogram: scale-rectified power across period
(log axis) and time, with the cone of influence faded at the ends and
the dominant-period ridge traced in white. The circadian band runs the
length of the
recording.](beyond-the-basics_files/figure-html/wavelet-scalogram-1.png)

The full wavelet power scalogram: scale-rectified power across period
(log axis) and time, with the cone of influence faded at the ends and
the dominant-period ridge traced in white. The circadian band runs the
length of the recording.

``` r

emd <- circadian.emd(as.numeric(counts10), t10, epoch_length = 600)
hilbert.huang(emd)
#> Hilbert-Huang Instantaneous Dynamics
#> 
#>   No circadian IMF to analyse
```

``` r

plot_emd(as.numeric(counts10), t10, epoch_length = 600)
```

![The empirical mode decomposition: intrinsic mode functions from finest
at the top to the slow trend at the bottom, with the circadian mode
(orange) the one nearest 24
hours.](beyond-the-basics_files/figure-html/emd-stack-1.png)

The empirical mode decomposition: intrinsic mode functions from finest
at the top to the slow trend at the bottom, with the circadian mode
(orange) the one nearest 24 hours.

### Registration, residual structure, and states

[`curve.registration()`](https://rdazadda.github.io/actiRhythm/reference/curve.registration.md)
aligns the days on their active-phase landmark and reports a
scale-invariant phase marker (the M10-window centre),
[`residual.spectrum()`](https://rdazadda.github.io/actiRhythm/reference/residual.spectrum.md)
examines what the cosinor leaves behind, and
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
fits a state-space model whose decoded path gives the probability of
rest across the day.

``` r

hmm <- rest.hmm(as.numeric(counts10), t10)
ggplot(hmm$tod_profile, aes(hour, p_rest)) +
  geom_col(fill = "#236192") +
  labs(x = "Hour of day", y = "P(rest)") +
  theme_actiRhythm()
```

![The 24-hour rest-probability profile from the state-space model: the
share of each clock hour the decoded path spends in the rest
state.](beyond-the-basics_files/figure-html/hmm-plot-1.png)

The 24-hour rest-probability profile from the state-space model: the
share of each clock hour the decoded path spends in the rest state.

## References

Golyandina, N., & Zhigljavsky, A. (2013). Singular spectrum analysis for
time series. *SpringerBriefs in Statistics*.
<https://doi.org/10.1007/978-3-642-34913-3>

Mormont, M. C., Waterhouse, J., Bleuzen, P., Giacchetti, S., Jami, A.,
Bogdan, A., Lellouch, J., Misset, J. L., Touitou, Y., & Levi, F. (2000).
Marked 24-h rest/activity rhythms are associated with better quality of
life, better response, and longer survival in patients with metastatic
colorectal cancer and good performance status. *Clinical Cancer
Research*, *6*(8), 3038–3045.

Wang, J., Xian, H., Licis, A., Deych, E., Ding, J., McLeland, J.,
Toedebusch, C., Li, T., Duntley, S., & Shannon, W. (2011). Measuring the
impact of apnea and obesity on circadian activity patterns using
functional linear modeling of actigraphy data. *Journal of Circadian
Rhythms*, *9*, 11. <https://doi.org/10.1186/1740-3391-9-11>
