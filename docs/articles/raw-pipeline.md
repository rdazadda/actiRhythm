# From raw acceleration

``` r

library(actiRhythm)
library(ggplot2)
```

This article continues from the get-started vignette and works from raw
accelerometer data rather than activity counts. It shows how actiRhythm
reads raw ActiGraph, Axivity, and GENEActiv files, auto-calibrates them,
derives the ENMO, MAD, and z-angle metrics, runs the full circadian
analysis and a posture-based sleep detector directly from a file, and
validates the whole chain against GGIR. The objects below are rebuilt
here so the article stands on its own.

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
raw <- example_raw(days = 2)        # synthetic 2-day recording (or a real file path)
```

## From raw acceleration

Everything so far runs on activity counts. Raw accelerometer files also
record the gravity component, and with it body posture, which counts
cannot carry. actiRhythm reads raw ActiGraph `.gt3x`, Axivity `.cwa`,
and GENEActiv `.bin` files, auto-calibrates them, and derives the three
raw metrics used across the field: ENMO, MAD, and the z-angle.

[`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
takes a file path. Raw acceleration is far too large to bundle, so for a
reproducible illustration
[`example_raw()`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)
synthesises a recording; a real file is read the same way, with
`raw.metrics("recording.cwa")`.

``` r

raw <- example_raw(days = 2)        # synthetic 2-day recording (or a real file path)
m <- raw.metrics(raw, epoch = 60)   # per-epoch ENMO (mg), MAD, and the z-angle
head(m)
#>                  time     ENMO      MAD    anglez
#> 1 2024-01-01 12:00:00 43.95883 24.78791  8.280445
#> 2 2024-01-01 12:01:00 43.13860 24.89003 14.666095
#> 3 2024-01-01 12:02:00 44.56626 25.55921  8.287697
#> 4 2024-01-01 12:03:00 44.25995 24.90778 14.705513
#> 5 2024-01-01 12:04:00 44.79032 25.36485  8.292428
#> 6 2024-01-01 12:05:00 44.17071 24.72920 14.755796
```

ENMO is the activity signal, and every method above takes it directly:
pass it to
[`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
or let
[`circadian.raw()`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md)
run the whole analysis from the file.

``` r

cr <- circadian.rhythm(m$ENMO, m$time)
c(IS = cr$IS, IV = cr$IV, RA = cr$RA)
#>     IS     IV     RA 
#> 1.0000 0.3823 0.9873
```

Calibration is applied first inside
[`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md).
[`auto.calibrate()`](https://rdazadda.github.io/actiRhythm/reference/auto.calibrate.md)
finds the per-axis gain and offset that return still periods to the 1 g
sphere ([Hees et al., 2014](#ref-vanhees2014)); on data with a known
distortion it returns the original gain and offset to within a fraction
of a percent.

``` r

set.seed(1)
u <- matrix(rnorm(40 * 3), 40, 3); u <- u / sqrt(rowSums(u^2))
cal <- do.call(rbind, lapply(seq_len(40), function(i)
  matrix(rep(u[i, ] / c(1.03, 0.97, 1.01) + c(0.04, -0.03, 0.02), each = 300),
         300, 3) + rnorm(900, 0, 0.004)))
auto.calibrate(data.frame(x = cal[, 1], y = cal[, 2], z = cal[, 3]),
               fs = 30)[c("scale", "offset")]
#> $scale
#> [1] 1.0301368 0.9701036 1.0099054
#> 
#> $offset
#> [1]  0.04000814 -0.02992012  0.01990532
```

The z-angle supports a sleep detector that needs no diary and works from
posture, which the counts cannot do.
[`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
finds the nightly sleep-period-time window from the distribution of
angle change \[HDCZA; Hees et al. ([2018](#ref-vanhees2018))\],
[`sib.vanhees()`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md)
scores sustained-inactivity bouts ([Hees et al.,
2015](#ref-vanhees2015)), and
[`sleep.from.spt()`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)
combines them into onset, wake and efficiency.
[`detect.nonwear.raw()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
flags a stationary, taken-off device \[low standard deviation and range
over the hour; Hees et al. ([2011](#ref-vanhees2011))\], and passing its
mask as `wear` keeps a device-off stretch from being read as one long
night. Computed here on the synthetic recording at a 5-second epoch,
with the detected windows shaded.

``` r

m5    <- raw.metrics(raw, epoch = 5, metrics = "anglez")
wear  <- detect.nonwear.raw(raw, epoch = 5)
spt   <- rest.spt(m5$anglez, m5$time, epoch_length = 5, wear = wear)
sib   <- sib.vanhees(m5$anglez, epoch_length = 5)
sleep <- sleep.from.spt(spt, sib, m5$time, epoch_length = 5)
thin <- seq(1, nrow(m5), by = 12)
ggplot() +
  geom_rect(data = spt, aes(xmin = onset, xmax = offset, ymin = -90, ymax = 90),
            fill = "#236192", alpha = 0.18) +
  geom_line(data = m5[thin, ], aes(time, anglez), linewidth = 0.2, colour = "grey30") +
  labs(x = NULL, y = "z-angle (deg)") +
  theme_actiRhythm()
```

![The z-angle of the synthetic recording at 5-second epochs, with the
two nightly sleep-period-time windows detected from the angle alone
(HDCZA), gated by raw non-wear,
shaded.](raw-pipeline_files/figure-html/spt-1.png)

The z-angle of the synthetic recording at 5-second epochs, with the two
nightly sleep-period-time windows detected from the angle alone (HDCZA),
gated by raw non-wear, shaded.

``` r

sleep[, c("date", "onset", "offset", "tst", "efficiency")]
#>         date               onset              offset      tst efficiency
#> 1 2024-01-01 2024-01-01 23:00:00 2024-01-02 07:00:00 8.001389  0.9968853
#> 2 2024-01-02 2024-01-02 23:00:00 2024-01-03 07:00:00 8.001389  0.9968853
```

Two more estimators measure fragmentation, and both run on the counts
already loaded:
[`activity.balance.index()`](https://rdazadda.github.io/actiRhythm/reference/activity.balance.index.md)
maps a detrended-fluctuation exponent to a 0 to 1 score that peaks at
the healthy 1/f balance ([Danilevicz et al.,
2024](#ref-danilevicz2024)), and
[`transition.probability()`](https://rdazadda.github.io/actiRhythm/reference/transition.probability.md)
gives the closed-form rest-to-active and active-to-rest transition
probabilities.

``` r

activity.balance.index(fractal.dfa(agd$axis1))
#> $ABI_overall
#> [1] 0.9900827
#> 
#> $ABI_short
#> [1] 0.9943423
#> 
#> $ABI_long
#> [1] 0.9820051
transition.probability(agd$axis1)[c("tp_ra_mle", "tp_ar_mle")]
#> $tp_ra_mle
#> [1] 0.04853705
#> 
#> $tp_ar_mle
#> [1] 0.1366603
```

### Agreement with GGIR

The raw metrics and the z-angle sleep detector reimplement the van Hees
algorithms that GGIR established ([Migueles et al.,
2019](#ref-migueles2019)), so the whole chain should reproduce GGIR’s
output, not just be internally consistent: calibration, the raw metrics,
and the z-angle sleep detection. Run side by side on the same real 7-day
wrist recording (118,800 five-second epochs), they agree at every stage.

Calibration matches to about three decimals, with scale 1.001 / 0.995 /
0.999 and the same 0.006 g calibration error, and the offsets agree in
magnitude once the two sign conventions are lined up: GGIR applies
`scale * x + offset`, actiRhythm `(x - offset) * scale`. Per-epoch ENMO
and z-angle both correlate at r = 0.99, with mean absolute differences
of 3.6 mg and 2.2 degrees. The HDCZA sleep-period window lands within
about five to ten seconds of GGIR each night, and total sleep time and
WASO agree to about a minute on the valid nights. Non-wear matches to
within 0.05 per worn day, and the wider gap between 0.41 and 0.52 is
entirely the parked-device tail that both tools already exclude.

The comparison is reproducible on your own file with GGIR installed:

``` r

library(GGIR)
P <- load_params(); P$params_general$windowsizes <- c(5, 900, 3600)
I <- g.inspectfile(file, params_rawdata = P$params_rawdata, params_general = P$params_general)
ggir <- g.getmeta(file, params_rawdata = P$params_rawdata, params_general = P$params_general,
                  params_cleaning = P$params_cleaning, inspectfileobject = I)$metashort
acti <- raw.metrics(file, epoch = 5)

# Align on wall-clock time: the .gt3x stores local time and the two tools anchor it
# slightly differently, so match the clock label rather than absolute seconds.
gkey <- gsub("T", " ", substr(as.character(ggir$timestamp), 1, 19))
akey <- format(acti$time, "%Y-%m-%d %H:%M:%S")
m <- merge(data.frame(k = gkey, ge = ggir$ENMO * 1000, ga = ggir$anglez),
           data.frame(k = akey, ae = acti$ENMO,        aa = acti$anglez), by = "k")
c(ENMO = cor(m$ge, m$ae), anglez = cor(m$ga, m$aa))
```

The sleep figures come from
[`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
and
[`sleep.from.spt()`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)
run against GGIR’s full part 1 to 4 pipeline (`HASPT.algo = "HDCZA"`).
On every clean night the window and the parameters match to seconds and
to about a minute. They diverge only on nights GGIR flags as heavily
invalid, 40 to 53% non-wear, where the two handle wear differently. That
is an edge case on badly degraded data, not a difference in the
algorithm.

ActiGraph idle-sleep gaps are imputed as zeros by `read.gt3x`, which
would otherwise collapse the z-angle to a constant during quiescent
periods. Like GGIR, actiRhythm carries the last gravity vector forward
through those gaps, so the angle stays correct through the still periods
the sleep detector relies on.

## References

Danilevicz, I. M., Hees, V. T. van, Heide, F. van der, Jacob, L.,
Landre, B., Benadjaoud, M. A., & Sabia, S. (2024). Measures of
fragmentation of rest activity patterns: Mathematical properties and
interpretability. *BMC Medical Research Methodology*, *24*, 132.
<https://doi.org/10.1186/s12874-024-02255-w>

Hees, V. T. van, Fang, Z., Langford, J., Assah, F., Mohammad, A., Silva,
I. C. M. da, Trenell, M. I., White, T., Wareham, N. J., & Brage, S.
(2014). Autocalibration of accelerometer data for free-living physical
activity assessment using local gravity and temperature: An evaluation
on four continents. *Journal of Applied Physiology*, *117*(7), 738–744.
<https://doi.org/10.1152/japplphysiol.00421.2014>

Hees, V. T. van, Renstrom, F., Wright, A., Gradmark, A., Catt, M., Chen,
K. Y., Lof, M., Bluck, L., Pomeroy, J., Wareham, N. J., Ekelund, U.,
Brage, S., & Franks, P. W. (2011). Estimation of daily energy
expenditure in pregnant and non-pregnant women using a wrist-worn
tri-axial accelerometer. *PLoS ONE*, *6*(7), e22922.
<https://doi.org/10.1371/journal.pone.0022922>

Hees, V. T. van, Sabia, S., Anderson, K. N., Denton, S. J., Oliver, J.,
Catt, M., Abell, J. G., Kivimaki, M., Trenell, M. I., & Singh-Manoux, A.
(2015). A novel, open access method to assess sleep duration using a
wrist-worn accelerometer. *PLoS ONE*, *10*(11), e0142533.
<https://doi.org/10.1371/journal.pone.0142533>

Hees, V. T. van, Sabia, S., Jones, S. E., Wood, A. R., Anderson, K. N.,
Kivimaki, M., Frayling, T. M., Pack, A. I., Bucan, M., Trenell, M. I.,
Mazzotti, D. R., Gehrman, P. R., Singh-Manoux, B. A., & Weedon, M. N.
(2018). Estimating sleep parameters using an accelerometer without sleep
diary. *Scientific Reports*, *8*, 12975.
<https://doi.org/10.1038/s41598-018-31266-z>

Migueles, J. H., Rowlands, A. V., Huber, F., Sabia, S., & Hees, V. T.
van. (2019). GGIR: A research community-driven open source R package for
generating physical activity and sleep outcomes from multi-day raw
accelerometer data. *Journal for the Measurement of Physical Behaviour*,
*2*(3), 188–196. <https://doi.org/10.1123/jmpb.2018-0063>
