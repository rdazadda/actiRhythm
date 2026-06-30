# Period and spectral analysis

``` r

library(actiRhythm)
library(ggplot2)
```

## The idea, and a rule to remember

The cosinor assumes the period; period estimation **finds** it. Given a
recording, these tools ask “what is the dominant cycle length?”, not
“how well does a 24-hour cosine fit?” That matters whenever the rhythm
may be free-running (a tau that is not exactly 24), as it is in jet lag,
shift work, dim-light protocols, and many clinical and animal
recordings. The two workhorses are the **Lomb-Scargle** periodogram, the
least-squares spectral estimator for unevenly sampled series ([Lomb,
1976](#ref-lomb1976); [Scargle, 1982](#ref-scargle1982)), and the
**chi-square (Sokolove-Bushell)** periodogram, the analysis-of-variance
method long standard in chronobiology ([Sokolove & Bushell,
1978](#ref-sokolove1978)).

The rule to carry through this article: **the chi-square periodogram
needs near-regular sampling, while Lomb-Scargle tolerates gaps.**
Actigraphy is rarely gap-free (non-wear, dropped epochs, mixed epoch
lengths), and a single device-off stretch can pull the chi-square peak
off the true period while Lomb-Scargle, which works from the actual
timestamps, stays on it. When in doubt with gappy data, trust
Lomb-Scargle.

## The math

Write the activity as values $`y_i`$ (mean-subtracted) sampled at times
$`t_i`$, $`i = 1, \dots, n`$, and let $`\omega = 2\pi / \tau`$ be the
angular frequency of a trial period $`\tau`$. The **Lomb-Scargle** power
at $`\tau`$ is

``` math
P(\omega) = \frac{1}{\sum_i y_i^2}\left[
  \frac{\bigl(\sum_i y_i \cos\omega(t_i - \delta)\bigr)^2}{\sum_i \cos^2\omega(t_i - \delta)}
  + \frac{\bigl(\sum_i y_i \sin\omega(t_i - \delta)\bigr)^2}{\sum_i \sin^2\omega(t_i - \delta)}
\right],
```

where the time offset $`\delta`$ is fixed by
$`\tan(2\omega\delta) = \sum_i \sin 2\omega t_i \big/ \sum_i \cos 2\omega t_i`$.
This offset is exactly what makes the estimate invariant to the
irregular spacing of the $`t_i`$: Lomb-Scargle is the least-squares fit
of a sinusoid at each trial period. Standard-normalised this way, the
power lies in $`[0, 1]`$, and its significance is read from the Baluev
([2008](#ref-baluev2008)) analytic false-alarm probability.

The **chi-square** periodogram instead folds the series at each integer
trial period $`P`$ (in epochs) into $`P`$ phase bins over
$`K = \lfloor N/P \rfloor`$ whole cycles and measures how much the phase
means $`\bar A_h`$ separate from the grand mean $`\bar A`$([Sokolove &
Bushell, 1978](#ref-sokolove1978)):

``` math
Q_P = \frac{K\,N \sum_{h=1}^{P}(\bar A_h - \bar A)^2}{\sum_{i=1}^{N}(A_i - \bar A)^2}.
```

Under the null of no rhythm at $`P`$, $`Q_P \sim \chi^2_{P-1}`$, which
supplies the significance line. The folding is where the regular-grid
assumption enters: it treats the $`i`$-th retained sample as if it sat
at phase $`i \bmod P`$, so a gap that shifts later samples off their
true phase distorts every $`Q_P`$.

## Assumptions, and when they break

- **A long enough span.** A period in the 18-30 h band cannot be
  resolved from less than two days;
  [`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
  enforces a \>= 2-day span and \>= 10 valid points and returns `NA`
  rather than guessing.
- **Stationarity.** Both periodograms assume one period holds across the
  recording. When tau drifts (re-entrainment, fragmentation), a single
  global number hides it. Use
  [`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
  to see the period over time.
- **Regular sampling, for the chi-square periodogram only.** Its
  phase-folding needs a near-uniform grid; gaps mis-align the fold.
  Lomb-Scargle carries the real times and is unbiased under irregular
  sampling, the centrepiece of this article.
- **A search window.** The peak is taken within `[from, to]` (default
  18-30 h). A rhythm outside the window cannot be found; widen it for
  ultradian or infradian work.

## Recovering known truth

Before trusting these tools on real data, watch them recover a period we
plant. We build a ten-day minute-level recording from a cosine whose
period is **25 hours, not 24** (a deliberately free-running tau), then
add noise.

``` r

ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 10 * 1440)
th <- as.numeric(difftime(ts, ts[1], units = "hours"))
set.seed(1)
tau_true <- 25
counts <- pmax(0, 100 + 80 * cos(2 * pi * (th - 8) / tau_true) + rnorm(length(th), 0, 15))

ls0  <- circadian.period(counts, ts)
chi0 <- chi.sq.periodogram(counts, ts)
knitr::kable(
  data.frame(method   = c("Lomb-Scargle", "chi-square"),
             planted  = c(tau_true, tau_true),
             tau_h    = c(ls0$tau, chi0$period)),
  digits = 2, caption = "Both periodograms recover the planted 25-hour period."
)
```

| method       | planted | tau_h |
|:-------------|--------:|------:|
| Lomb-Scargle |      25 | 25.26 |
| chi-square   |      25 | 25.05 |

Both periodograms recover the planted 25-hour period. {.table}

Both land on 25, not the conventional 24: the method reports the period
we put in, not the one we expected. The confidence interval should then
**bracket** the planted value.
[`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md)
attaches a circular moving-block bootstrap interval that respects the
autocorrelation of activity ([Kunsch, 1989](#ref-kunsch1989); [Politis &
Romano, 1992](#ref-politis1992)).

``` r

ci <- period.ci(counts, ts, n_boot = 200, seed = 1)
ci
#> Circadian Period with Bootstrap Confidence Interval
#> 
#>   tau:      25.027 h
#>   95% CI:   [24.987, 25.012] h
#>   SE:       0.006 h
#>   Method:   circular block residual bootstrap (200/200 valid reps)
tau_true >= ci$ci_lower && tau_true <= ci$ci_upper   # the interval brackets the truth
#> [1] TRUE
```

## On a real recording

The bundled recording runs the same way.
[`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
returns a list whose `tau` is the estimated endogenous period and whose
`p_value` is the Baluev false-alarm probability.

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
ls  <- circadian.period(agd$axis1, agd$timestamp)
c(tau = ls$tau, peak_power = ls$peak_power, p_value = ls$p_value,
  n_used = ls$n_used, span_days = round(ls$span_days, 2))
#>          tau   peak_power      p_value       n_used    span_days 
#> 2.448889e+01 4.418582e-02 4.916458e-98 9.919000e+03 6.890000e+00
```

``` r

plot_periodogram(agd$axis1, agd$timestamp)
```

![The Lomb-Scargle periodogram of the bundled recording over the 18-30 h
window. The orange line marks the dominant period, the grey dashed line
is 24 hours, and the red dashed line is the 0.05 Baluev false-alarm
threshold, which is flat across periods (unlike the chi-square critical
curve below). The peak sits to the right of 24, a longer-than-day
dominant cycle.](period-estimation_files/figure-html/ls-plot-1.png)

The Lomb-Scargle periodogram of the bundled recording over the 18-30 h
window. The orange line marks the dominant period, the grey dashed line
is 24 hours, and the red dashed line is the 0.05 Baluev false-alarm
threshold, which is flat across periods (unlike the chi-square critical
curve below). The peak sits to the right of 24, a longer-than-day
dominant cycle.

The bootstrap interval shows how firmly that period is pinned down. On
this recording the replicate periods spread widely across the search
window, so the interval spans several hours: a single tau on a noisy
free-running record carries real uncertainty.

``` r

plot_period_ci(agd$axis1, agd$timestamp, n_boot = 200, seed = 1)
```

![Bootstrap replicate periods (histogram) with the point estimate
(orange) and the 95% confidence-interval band. A broad spread means the
dominant period is only weakly
determined.](period-estimation_files/figure-html/ci-plot-1.png)

Bootstrap replicate periods (histogram) with the point estimate (orange)
and the 95% confidence-interval band. A broad spread means the dominant
period is only weakly determined.

## Reading the numbers

State each output in human terms:

- **tau** is the dominant period in hours. Near 24 is entrained to the
  day; clearly above or below 24 suggests a free-running or
  mis-entrained rhythm. Read it against the 24-hour line, not in
  isolation.
- **peak_power** (Lomb-Scargle) runs 0 to 1; higher is a sharper, more
  dominant cycle. **Qp_peak** (chi-square) is on a chi-square scale and
  is judged against its per-period critical line.
- **p_value** is the chance a peak this tall arose from noise.
  Lomb-Scargle uses the Baluev analytic false-alarm probability; the
  chi-square periodogram reports a family-wise (Sidak) p-value across
  the scanned periods, so a single tall bar is not called real purely
  from multiple testing.
- **The confidence interval** is the honest precision of tau. A wide
  interval that straddles 24 is the data telling you it cannot separate
  an entrained from a free-running rhythm.

## When the sampling is gappy

This is the rule in action. We take the same planted-25 h series and
**knock out a contiguous block of epochs** (a roughly two-day device-off
gap), then re-estimate. Lomb-Scargle reads the real timestamps and
should hold at 25; the chi-square fold, which assumes a regular grid,
should slip.

``` r

gap   <- th >= 5 * 24 & th < 7.3 * 24          # remove ~2.3 days in the middle
ts_g  <- ts[!gap]
cnt_g <- counts[!gap]

ls_g  <- circadian.period(cnt_g, ts_g)
chi_g <- chi.sq.periodogram(cnt_g, ts_g)
knitr::kable(
  data.frame(method = c("Lomb-Scargle", "chi-square"),
             truth  = c(tau_true, tau_true),
             full   = c(ls0$tau,  chi0$period),
             gappy  = c(ls_g$tau, chi_g$period)),
  digits = 2,
  caption = "After a multi-day gap, Lomb-Scargle stays on 25 h; the chi-square peak slips off."
)
```

| method       | truth |  full | gappy |
|:-------------|------:|------:|------:|
| Lomb-Scargle |    25 | 25.26 | 25.26 |
| chi-square   |    25 | 25.05 | 24.13 |

After a multi-day gap, Lomb-Scargle stays on 25 h; the chi-square peak
slips off. {.table}

The chi-square period moves off 25 (pulled toward a wrong value because
the samples after the gap no longer fold onto their true phase), while
Lomb-Scargle, anchored to the actual times, returns essentially the same
tau as before. The lesson is the rule: with gappy actigraphy, reach for
Lomb-Scargle and treat a chi-square peak with suspicion until you have
confirmed the sampling is near-regular.

## The wider period family

The same spectral idea generates a family of tools that look past a
single global number.

**Period over time.**
[`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
slides a window across the recording and runs the chi-square periodogram
in each, producing a period-by-time heat map that shows drift,
fragmentation, and re-entrainment a single fit averages away.

``` r

sg <- circadian.spectrogram(agd$axis1, agd$timestamp, step_hours = 12)
sg$plot
```

![Sliding-window chi-square spectrogram of the bundled recording,
coloured by Qp relative to its per-period significance threshold so
equal colour means equal significance across periods. The white line
traces the per-window peak period and the dashed line is 24
h.](period-estimation_files/figure-html/spectrogram-1.png)

Sliding-window chi-square spectrogram of the bundled recording, coloured
by Qp relative to its per-period significance threshold so equal colour
means equal significance across periods. The white line traces the
per-window peak period and the dashed line is 24 h.

**Ultradian band power.** Rhythms faster than a day (the roughly
90-minute, 4-hour, and 8-hour bands) carry real physiology.
[`ultradian.bandpower()`](https://rdazadda.github.io/actiRhythm/reference/ultradian.bandpower.md)
partitions the activity variance into dyadic period bands with an
undecimated Haar wavelet and reports the fraction of variance in each
([Percival & Walden, 2000](#ref-percival2000)).

``` r

ultradian.bandpower(agd$axis1, agd$timestamp)
#> Ultradian Wavelet Band Power
#> 
#>   band low_h high_h    power fraction
#>  90min     1      2 186918.1    0.178
#>     4h     2      6 269890.2    0.257
#>     8h     6     12 122031.3    0.116
```

**What the cosinor leaves behind.**
[`residual.spectrum()`](https://rdazadda.github.io/actiRhythm/reference/residual.spectrum.md)
removes the fitted 24-hour cosinor mean and estimates the spectrum of
the residual, splitting it into ultradian and high-frequency bands. Two
recordings with the same 24-hour rhythm but different residual
fragmentation are told apart here.

``` r

residual.spectrum(agd$axis1, agd$timestamp, period = 24)
#> Residual Circadian Spectrum
#> 
#>   Residual variance: 1009357.4
#> 
#>       band low_h high_h     power fraction
#>  ultradian   2.0      8 995304256    0.197
#>  high_freq   0.5      2 925775830    0.183
```

## Limitations

- **Resolution scales with span.** The period grid is set by the
  recording length; short records give coarse, wide intervals.
  Distinguishing 24.0 from 24.5 h needs many days.
- **The chi-square periodogram is fragile to gaps.** As shown above,
  treat its peak with care on irregular sampling; prefer Lomb-Scargle,
  or gate and report wear time first.
- **Stationarity is assumed.** A single tau hides a drifting period; use
  the spectrogram when re-entrainment or fragmentation is plausible.
- **The window bounds the answer.** A rhythm outside `[from, to]` is
  invisible; set the window to the science, not the default.

## Reference and validation

The Lomb-Scargle estimator follows Lomb ([1976](#ref-lomb1976)) and
Scargle ([1982](#ref-scargle1982)), with the analytic false-alarm
probability of Baluev ([2008](#ref-baluev2008)) and the actigraphy
application of Ruf ([1999](#ref-ruf1999)); the chi-square periodogram is
that of Sokolove & Bushell ([1978](#ref-sokolove1978)). The period
confidence interval uses the moving-block bootstrap of Kunsch
([1989](#ref-kunsch1989)) in its circular form ([Politis & Romano,
1992](#ref-politis1992)). actiRhythm’s Lomb-Scargle power is
cross-checked against a direct least-squares periodogram (the exact
estimator: it uses the double-angle time offset of Scargle
([1982](#ref-scargle1982)), where the `lomb` package uses a single-angle
offset), and the chi-square $`Q_P`$ statistic is validated against its
chi-square null, in the package’s test suite.

## References

Baluev, R. V. (2008). Assessing the statistical significance of
periodogram peaks. *Monthly Notices of the Royal Astronomical Society*,
*385*(3), 1279–1285. <https://doi.org/10.1111/j.1365-2966.2008.12689.x>

Kunsch, H. R. (1989). The jackknife and the bootstrap for general
stationary observations. *The Annals of Statistics*, *17*(3), 1217–1241.
<https://doi.org/10.1214/aos/1176347265>

Lomb, N. R. (1976). Least-squares frequency analysis of unequally spaced
data. *Astrophysics and Space Science*, *39*(2), 447–462.
<https://doi.org/10.1007/BF00648343>

Percival, D. B., & Walden, A. T. (2000). *Wavelet methods for time
series analysis*. Cambridge University Press.
<https://doi.org/10.1017/CBO9780511841040>

Politis, D. N., & Romano, J. P. (1992). A circular block-resampling
procedure for stationary data. In R. LePage & L. Billard (Eds.),
*Exploring the limits of bootstrap* (pp. 263–270). Wiley.

Ruf, T. (1999). The Lomb-Scargle periodogram in biological rhythm
research: Analysis of incomplete and unequally spaced time-series.
*Biological Rhythm Research*, *30*(2), 178–201.
<https://doi.org/10.1076/brhm.30.2.178.1422>

Scargle, J. D. (1982). Studies in astronomical time series analysis. II.
Statistical aspects of spectral analysis of unevenly spaced data. *The
Astrophysical Journal*, *263*, 835–853. <https://doi.org/10.1086/160554>

Sokolove, P. G., & Bushell, W. N. (1978). The chi square periodogram:
Its utility for analysis of circadian rhythms. *Journal of Theoretical
Biology*, *72*(1), 131–160.
<https://doi.org/10.1016/0022-5193(78)90022-X>
