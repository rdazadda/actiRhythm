# Getting started with actiRhythm

You have an ActiGraph recording for one person and a practical question.
Does this person have a stable daily rhythm, how strong is it, and is it
real or could it be noise? This article answers that question end to end
on one bundled recording. By the end you will be able to read the file,
describe the rhythm without assuming a shape, test whether it is real,
locate and bound its period, probe its finer structure, and write out a
one-row summary you can stack across people.

The method assumes a multi-day activity-count series with timestamps,
which is what an ActiGraph `.agd` file holds and what counts exported
from any source can supply. Every number and figure below is
reproducible from the installed package.

``` r

library(actiRhythm)
library(ggplot2)
```

## The recording

[`read.agd()`](https://rdazadda.github.io/actiRhythm/reference/read.agd.md)
opens the `.agd` SQLite file and auto-detects the epoch length and axes;
[`agd.counts()`](https://rdazadda.github.io/actiRhythm/reference/agd.counts.md)
returns a tidy data frame with `$axis1` (the vertical-axis counts) and
`$timestamp`. One recording ships with the package, reached through
[`example_agd()`](https://rdazadda.github.io/actiRhythm/reference/example_agd.md).

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
head(agd[, c("timestamp", "axis1")])
#>             timestamp axis1
#> 1 2025-10-07 20:28:00  2676
#> 2 2025-10-07 20:29:00  1380
#> 3 2025-10-07 20:30:00  4126
#> 4 2025-10-07 20:31:00  5763
#> 5 2025-10-07 20:32:00  4076
#> 6 2025-10-07 20:33:00  4043
```

Read the actogram top to bottom before computing anything. It is the
fastest way to see whether there is a rhythm to measure at all.

``` r

plot_actogram(agd$axis1, agd$timestamp)
```

![Double-plotted actogram: each row is one day shown twice across 48
hours, so a stable rhythm forms a vertical band of activity at the same
clock time each day; scattered fill marks a fragmented or shifting
rhythm.](actiRhythm_files/figure-html/actogram-1.png)

Double-plotted actogram: each row is one day shown twice across 48
hours, so a stable rhythm forms a vertical band of activity at the same
clock time each day; scattered fill marks a fragmented or shifting
rhythm.

In this recording the active rows do not stack into one clean vertical
band; they slide a little later or earlier from day to day, which is the
low interdaily stability the next section puts a number on.

## The nonparametric rhythm

[`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
describes the rest-activity pattern without assuming a waveform. It
returns interdaily stability (IS) and intradaily variability (IV),
introduced by Witting et al. ([1990](#ref-witting1990)), together with
the relative amplitude (RA) and the least- and most-active windows L5
and M10 of Van Someren et al. ([1999](#ref-vansomeren1999)), each with
its onset time.

``` r

rhythm <- circadian.rhythm(agd$axis1, agd$timestamp)
rhythm
#> 
#> Circadian Rhythm Analysis
#> 
#> Data Summary
#>   Days analyzed:            8
#>   Valid circadian days:     6
#>   Epoch length:             60 seconds
#> 
#> Non-Parametric Metrics (IS/IV: Witting et al. 1990; RA/L5/M10: van Someren et al. 1999)
#>   L5 (least active 5h):     6.10 counts/min, onset 01:28
#>   M10 (most active 10h):    602.93 counts/min, onset 12:02
#>   L1 (least active 1h):     3.10 counts/min, onset 02:35
#>   M1 (most active 1h):      1175.51 counts/min, onset 17:43
#>   Relative Amplitude (RA):  0.9800 (range 0-1, higher=stronger rhythm)
#>   Interdaily Stability (IS): 0.2279 (range 0-1, higher=more consistent)
#>   Intradaily Variability (IV): 1.0008 (near 0 = sine, near 2 = noise)
#>   Phi (autocorrelation):    0.4151 (higher=more predictable)
#> 
#> Sleep-Based & Variability Metrics
#>   Sleep Regularity Index:   Not calculated (requires sleep_state input)
#>   Onset timing variability: 2.23 hours
#>   L5 timing variability:    0.77 hours (circular SD)
#>   M10 timing variability:   3.70 hours (circular SD)
#> 
#> References: Witting (1990), van Someren (1999)
```

Read the three headline numbers together. IS ranges 0 to 1 and rises as
the daily pattern repeats more exactly from day to day; IV runs from
about 0 (a smooth sine) to about 2 (noisy or split into ultradian
bouts); RA ranges 0 to 1 and grows with the contrast between an active
day and a restful night. Here a high RA (about 0.98) sits with a low IS
(about 0.23): a strong day-night contrast carried on irregular timing.
The header also notes that six of the eight calendar days qualified,
because
[`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
sets aside days with too little valid recording before computing the
metrics.

## Cosinor and the rhythmicity test

[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
fits a 24-hour cosine and returns the MESOR (rhythm-adjusted mean),
amplitude, and acrophase, the clock time of the peak ([Cornelissen,
2014](#ref-cornelissen2014)).

``` r

cos <- cosinor.analysis(agd$axis1, agd$timestamp, period = 24)
c(mesor = unname(cos$mesor), amplitude = unname(cos$amplitude),
  acrophase = unname(cos$acrophase))
#>     mesor amplitude acrophase 
#>    327.25    295.63     16.92
```

The acrophase is **circular**: it is a clock time, so 23.5 and 0.5 are
an hour apart, not 23 hours apart, and it is only interpretable once the
amplitude is distinguishable from zero.
[`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
is that check, the zero-amplitude F-test of Nelson et al.
([1979](#ref-nelson1979)), reported with the percent of variance the
rhythm explains.

``` r

rhythmicity.test(agd$axis1, agd$timestamp, cosinor_result = cos)
#> Cosinor Rhythmicity Test (Halberg zero-amplitude F-test)
#> 
#> H0: amplitude = 0 (no rhythm)
#>   Period:         24 h
#>   F(2, 21):        8.15
#>   P-value:        0.0024
#>   Percent rhythm: 43.7% (R-squared = 0.4371)
#>   Rhythmic:       YES (alpha = 0.05)
```

A small p-value rejects the flat-line null. Here the test is clearly
significant (p well under 0.05) and the single cosine explains a little
over forty percent of the variance, so the acrophase above is worth
reading. The percent-rhythm also reminds you that the rest of the
variance, more than half, is a mix of noise and structure a single
cosine cannot capture. Overlay the fit to see what it is and is not
catching.

``` r

plot_extended_cosinor(agd$axis1, agd$timestamp)
```

![Hourly activity with the fitted 24-hour cosine overlaid; the peak of
the curve is the acrophase and its height above the MESOR line is the
amplitude. Systematic departures from the curve are where a single
cosine is too smooth for real rest-activity
data.](actiRhythm_files/figure-html/cosinor-plot-1.png)

Hourly activity with the fitted 24-hour cosine overlaid; the peak of the
curve is the acrophase and its height above the MESOR line is the
amplitude. Systematic departures from the curve are where a single
cosine is too smooth for real rest-activity data.

Look where the fitted curve and the hourly points part company: the
single cosine runs too smooth through the early hours, and that gap is
the structure the percent-rhythm says one cosine leaves on the table.
When you need the joint uncertainty in amplitude and acrophase together,
`cosinor.confidence.ellipse(cos)` returns the Bingham et al.
([1982](#ref-bingham1982)) ellipse; an ellipse that encloses the origin
is the geometric form of “no detectable rhythm”.

## The period

A cosinor assumes exactly 24 hours.
[`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
does not: it runs a Lomb-Scargle periodogram over the 18 to 30 hour band
and returns the dominant period `tau`, with a false-alarm probability in
`$p_value` ([Baluev, 2008](#ref-baluev2008); [Lomb,
1976](#ref-lomb1976); [Scargle, 1982](#ref-scargle1982)).

``` r

per <- circadian.period(agd$axis1, agd$timestamp)
c(tau = per$tau, p_value = per$p_value)
#>          tau      p_value 
#> 2.448889e+01 4.916458e-98
```

The strongest cycle here sits near 25.4 hours, meaningfully longer than
a solar day, and the false-alarm probability is essentially zero, so the
cycle itself is not in doubt. Its exact length is less certain.
[`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md)
refines that peak and bootstraps it with a circular block bootstrap
([Kunsch, 1989](#ref-kunsch1989); [Politis & Romano,
1992](#ref-politis1992)), which respects the autocorrelation in activity
data. The interval it returns is wide enough to span 24 hours. On a
single recording the period is clearly present but not pinned to the
hour. Passing `seed = 1` makes the bootstrap reproducible; the 50
replicates here keep the vignette fast, and a real interval would use
more.

``` r

period.ci(agd$axis1, agd$timestamp, n_boot = 50, seed = 1)
#> Circadian Period with Bootstrap Confidence Interval
#> 
#>   tau:      24.895 h
#>   95% CI:   [23.373, 26.193] h
#>   SE:       1.195 h
#>   Method:   circular block residual bootstrap (50/50 valid reps)
```

``` r

plot_periodogram(agd$axis1, agd$timestamp)
```

![Lomb-Scargle spectral power across candidate periods; the tallest peak
is the dominant rhythm, and a peak rising above the significance line
rejects the no-rhythm
null.](actiRhythm_files/figure-html/periodogram-1.png)

Lomb-Scargle spectral power across candidate periods; the tallest peak
is the dominant rhythm, and a peak rising above the significance line
rejects the no-rhythm null.

The chi-square periodogram of Sokolove & Bushell
([1978](#ref-sokolove1978)) asks the same question a different way,
scoring each trial period by how much of the variance folds onto it. A
second estimator that clears its own threshold at the same period is
reassurance that the peak is real rather than an artifact of one method.

``` r

plot_chisq(agd$axis1, agd$timestamp)
```

![Chi-square (Sokolove-Bushell) periodogram: the Qp statistic against
trial period, with the significance threshold. The peak should fall at
the period the Lomb-Scargle periodogram
found.](actiRhythm_files/figure-html/chisq-1.png)

Chi-square (Sokolove-Bushell) periodogram: the Qp statistic against
trial period, with the significance threshold. The peak should fall at
the period the Lomb-Scargle periodogram found.

## Watching the period drift

A single periodogram averages over the whole recording.
[`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
slides a window across it and recomputes the periodogram in each, so you
can see the dominant period move, the chronobiology equivalent of a
music spectrogram.

``` r

circadian.spectrogram(agd$axis1, agd$timestamp, step_hours = 24)$plot
```

![Period (vertical) against time (horizontal), coloured by spectral
power. A horizontal band of colour is a stable period; a band that bends
upward or downward is a rhythm lengthening or shortening across the
recording.](actiRhythm_files/figure-html/spectrogram-1.png)

Period (vertical) against time (horizontal), coloured by spectral power.
A horizontal band of colour is a stable period; a band that bends upward
or downward is a rhythm lengthening or shortening across the recording.

Follow the band of colour left to right: here it drifts upward, the
dominant period lengthening across the recording rather than holding at
24 hours. Read the later windows with some caution, since they ride the
upper edge of the search band, where a weakening rhythm and a genuinely
longer period can look alike.

## Nonlinear structure

Two recordings can share the same IS and amplitude yet differ in their
moment-to-moment dynamics.
[`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
estimates the detrended-fluctuation exponent alpha, the long-range
temporal correlation in the series ([Hu et al., 2001](#ref-hu2001);
[Peng et al., 1994](#ref-peng1994)).

``` r

fractal.dfa(agd$axis1)$alpha
#> [1] 0.9263546
```

An alpha near 1 is the “1/f” scaling typical of healthy physiological
output; near 0.5 is uncorrelated noise; well above 1 is a smoother, more
random-walk-like signal. Here alpha comes out near 0.93, close to the
1/f end, so this series carries genuine long-range correlation rather
than looking like noise. The log-log fit shows whether one exponent
holds across all time scales.

``` r

plot_dfa(agd$axis1)
```

![Detrended fluctuation on log-log axes; the slope is the DFA exponent
alpha. One straight line is monofractal scaling, a visible kink is a
crossover between short- and long-range
dynamics.](actiRhythm_files/figure-html/dfa-plot-1.png)

Detrended fluctuation on log-log axes; the slope is the DFA exponent
alpha. One straight line is monofractal scaling, a visible kink is a
crossover between short- and long-range dynamics.

Two further questions use the same series and a single call each.
[`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
asks how the complexity of the signal changes as you coarse-grain it
across time scales ([Costa et al., 2002](#ref-costa2002)), and
[`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md)
asks whether one scaling exponent is enough or the series is genuinely
multifractal ([Kantelhardt et al., 2002](#ref-kantelhardt2002)). Their
help pages carry the math, and the [Beyond the
basics](https://rdazadda.github.io/actiRhythm/articles/beyond-the-basics.html)
article runs them on this recording; here it is enough to know they
refine, not replace, the single alpha above.

## Rest-activity transitions

[`state.transitions()`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
summarises how readily the subject switches between rest and activity,
the kRA (rest-to-active) and kAR (active-to-rest) rates that capture
fragmentation a single amplitude cannot ([Lim et al.,
2011](#ref-lim2011)).

``` r

st <- state.transitions(agd$axis1)
c(kRA = st$kRA, kAR = st$kAR)
#>        kRA        kAR 
#> 0.03920806 0.10681627
```

For this subject kAR runs well above kRA: once active they settle back
to rest quickly, but once at rest they stay there. That asymmetry is the
fragmentation a single amplitude cannot see.

## Pinpointing sleep and wake

[`state.transitions()`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
gives one pair of rates for the whole recording.
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
instead locates the sleep-onset and wake-onset of each night directly
from the counts, with no scored sleep required: a 24-hour cosinor bounds
each rest and active span roughly, and a change point inside each bound
places the precise transition ([Chen & Sun, 2024](#ref-chensun2024)).

``` r

cp <- sleep.changepoints(agd$axis1, agd$timestamp)
cp
#> Change-Point Sleep/Wake Detection
#> 
#>   Span:           6.9 days (9919 epochs)
#>   Cosinor acrophase: 20.5 h
#>   Change points:  14 (7 sleep episodes)
#>   Mean sleep duration: 11.1 h
#> 
#>   First sleep episodes:
#>     sleep 10-07 22:36  ->  wake 10-08 07:09  (8.6 h)
#>     sleep 10-08 23:30  ->  wake 10-09 07:50  (8.3 h)
#>     sleep 10-09 21:30  ->  wake 10-10 12:11  (14.7 h)
#> 
#>   Reference: Chen and Sun (2024)
```

``` r

cps <- cp$changepoints
cps$clock <- as.numeric(format(cps$time, "%H")) + as.numeric(format(cps$time, "%M")) / 60
ggplot(cps, aes(time, clock, colour = type)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("sleep onset" = "#236192", "wake onset" = "#E69F00")) +
  labs(x = "Date", y = "Clock hour", colour = NULL) +
  theme_actiRhythm()
```

![Each detected sleep onset (blue) and wake onset (orange) by date and
clock time. The detector reads the per-night sleep and wake timing
directly from the counts.](actiRhythm_files/figure-html/cp-plot-1.png)

Each detected sleep onset (blue) and wake onset (orange) by date and
clock time. The detector reads the per-night sleep and wake timing
directly from the counts.

For this recording the detector finds seven nights, sleep onset late in
the evening and wake onset around seven in the morning, with most nights
near eight to nine hours of rest. That is the per-night sleep and wake
timing the single transition rate above cannot localise.

## Every rest bout, naps included

[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
finds the one main rest bout of each cycle.
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
takes the complementary view, consolidating every spell of low activity
into a rest bout, however many a day holds, daytime naps and fragmented
rest included ([Roenneberg et al., 2015](#ref-roenneberg2015)). Each
epoch is compared to a fraction of its own 24-hour activity level, and
runs that stay below it grow into consolidated bouts.

``` r

rp <- rest.periods(agd$axis1, agd$timestamp)
rp
#> Consolidated Rest Periods (Roenneberg/MASDA)
#> 
#>   Bouts detected:   9  (1.5 per day over 6 days)
#>   Total rest:       61.3 h
#>   Main bouts / naps: 9 / 0
#> 
#>   First bouts:
#>     main 10-08 00:49 -> 10-08 07:09  (381 min)
#>     main 10-08 23:31 -> 10-09 07:50  (500 min)
#>     main 10-10 00:12 -> 10-10 06:57  (406 min)
#> 
#>   Reference: Roenneberg et al. (2015); Loock et al. (2021)
```

``` r

ggplot(rp$rest_periods, aes(onset, bout, xend = offset, yend = bout, colour = type)) +
  geom_segment(linewidth = 3) +
  scale_colour_manual(values = c(main = "#236192", nap = "#E69F00")) +
  labs(x = "Time", y = "Rest bout", colour = NULL) +
  theme_actiRhythm()
```

![Each consolidated rest bout as a bar from onset to offset. Detecting
more than one bout a day, including naps, is what separates this from
the one-per-cycle change-point
detector.](actiRhythm_files/figure-html/restperiods-plot-1.png)

Each consolidated rest bout as a bar from onset to offset. Detecting
more than one bout a day, including naps, is what separates this from
the one-per-cycle change-point detector.

For this recording it consolidates nine rest bouts across six days, more
than one a day, so a few days carry a second rest beyond the main night.

[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
reaches a comparable result by a different route, rank-order filtering
and binary morphology rather than seed-and-grow consolidation ([Crespo
et al., 2012](#ref-crespo2012)); running both gives an independent
cross-check of where a recording’s rest bouts fall. Neither applies a
wear-time filter, so gate the counts on valid wear first if a recording
has device-off stretches.

## Scoring sleep and its regularity

The metrics so far describe the rhythm; the next need a per-epoch sleep
or wake label.
[`sleep.cole.kripke()`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)
produces that label directly from the counts, the standard count-based
classifier for adults ([Cole et al., 1992](#ref-cole1992)), with
[`sleep.sadeh()`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md)
as the children-and-adolescents alternative.

``` r

state <- sleep.cole.kripke(agd$axis1)
table(state)
#> state
#>    S    W 
#> 7641 2278
```

Count-based scoring labels every low-activity epoch as sleep, so a
mostly sedentary recording scores a high sleep fraction; what matters
downstream is how regular that label is from day to day.
[`sleep.regularity.index()`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
measures exactly that, the probability that two epochs 24 hours apart
share the same state ([Phillips et al., 2017](#ref-phillips2017)).

``` r

sleep.regularity.index(state, agd$timestamp)
#> [1] 58.27
```

When you also have explicit sleep periods,
[`social.jet.lag()`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md)
returns the work-day to free-day mid-sleep difference ([Wittmann et al.,
2006](#ref-wittmann2006)) and
[`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md) the
locomotor inactivity during sleep that tracks ultradian sleep structure
([Winnebeck et al., 2018](#ref-winnebeck2018)).

## Pooling the evidence

Each test above answers “is there a rhythm?” in its own language.
[`consensus.rhythmicity()`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
runs three of them, the cosinor F-test, the Lomb-Scargle false-alarm
probability, and the chi-square periodogram, and reports both a majority
vote and a Cauchy-combined p-value ([Liu & Xie, 2020](#ref-liu2020)).

``` r

consensus.rhythmicity(agd$axis1, agd$timestamp)
#> Consensus Rhythmicity (multi-method)
#> 
#>   cosinor F-test         p=0.0024  rhythmic
#>   Lomb-Scargle FAP       p=<2e-16  rhythmic
#>   chi-square periodogram p=1       no
#> 
#>   Votes:        2 / 3 methods
#>   Combined p:   <2e-16
#>   Consensus:    RHYTHMIC (alpha = 0.05)
```

For this recording the tests agree that a rhythm is present and the
combined p-value is effectively zero, so the evidence leans firmly
toward a real rhythm even where one method on its own might hesitate.
The three tests run on the same series, so they are not independent; the
Cauchy combination stays valid under that dependence, where Fisher’s
method would not.

## A one-row summary

A one-row summary is then easy to assemble and
[`rbind()`](https://rdrr.io/r/base/cbind.html) across subjects.

``` r

data.frame(
  IS        = rhythm$IS,
  IV        = rhythm$IV,
  RA        = rhythm$RA,
  tau       = per$tau,
  amplitude = cos$amplitude,
  acrophase = cos$acrophase,
  dfa_alpha = fractal.dfa(agd$axis1)$alpha
)
#>              IS     IV   RA      tau amplitude acrophase dfa_alpha
#> cos_term 0.2279 1.0008 0.98 24.48889    295.63     16.92 0.9263546
```

Read across, this row is one person: a strongly amplitude-modulated
rhythm whose timing is irregular from day to day, running on a period a
little over 24 hours.

## A whole study at once

You have now run the whole arc on one recording, ending with a single
summary row. Most studies, though, are a folder of recordings rather
than one, and that is what
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
is for. It runs this entire arc over every `.agd` file in a directory
and returns one row per recording, so a cohort becomes a single data
frame you can model. A file that fails to read is reported in an `error`
column instead of stopping the run. The package ships two recordings, so
this runs as a real, if small, batch.

``` r

batch <- circadian.batch(
  system.file("extdata", package = "actiRhythm"),
  verbose = FALSE
)
batch[, c("file", "IS", "IV", "RA", "period_tau", "rhythm_p_value")]
#>                      file     IS     IV     RA period_tau rhythm_p_value
#> 1 MOS2E39230594_60sec.agd 0.2279 1.0008 0.9800   24.48889    0.002396228
#> 2  MOS2E3923063660sec.agd 0.4230 1.2687 0.7706   24.09425    0.050151893
```

Each row is one subject, described by the same metrics as the
walkthrough above, and already the two differ: the second is steadier
from day to day (higher IS) but less amplitude-modulated (lower RA), and
its rhythm clears significance only at the margin. Point the call at
your own folder and nothing else changes.

When you want the full analysis for a single recording rather than a
summary row,
[`circadian.workbook()`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)
writes every table to a multi-sheet Excel file: a summary row plus the
hourly profile, both periodograms, the fluctuation and multifractal
spectra, and a data dictionary.

``` r

circadian.workbook(agd$axis1, agd$timestamp, file = "subject01.xlsx")
```

## Where next

- Decomposition, the multi-harmonic profile, and the finer nonparametric
  and phase metrics are demonstrated in [Beyond the
  basics](https://rdazadda.github.io/actiRhythm/articles/beyond-the-basics.html).
- Reading raw `.gt3x`/`.cwa`/`.bin` files, calibration, the
  ENMO/MAD/z-angle metrics, and diary-free posture-based sleep are in
  [From raw
  acceleration](https://rdazadda.github.io/actiRhythm/articles/raw-pipeline.html).
- When a rhythm drifts, shifts phase, or fragments across the recording,
  see [Nonstationary and complex
  rhythms](https://rdazadda.github.io/actiRhythm/articles/nonstationary-rhythms.html).
- [`?circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
  [`?cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
  [`?circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md),
  [`?fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md),
  [`?consensus.rhythmicity`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
  give the per-function help with full references.
- [`?circadian.batch`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
  and
  [`?circadian.workbook`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)
  are the study-scale entry points.
- For a multi-subject design,
  [`population.cosinor()`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)
  and
  [`cosinor.compare()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.compare.md)
  give the mean rhythm and a two-group test (Bingham/Hotelling); for
  covariates, repeated measures, or nesting, feed the per-subject
  metrics this package emits (MESOR, amplitude, acrophase, IS, IV, RA,
  period) into a mixed model with `lme4`/`nlme` alongside your design.
  The package gives you the metrics; you supply the multilevel model.
- Cite the package with `citation("actiRhythm")`, and report
  `packageVersion("actiRhythm")` with your results so the analysis stays
  reproducible across releases.
- Questions and bugs: <https://github.com/rdazadda/actiRhythm/issues>.

## References

Baluev, R. V. (2008). Assessing the statistical significance of
periodogram peaks. *Monthly Notices of the Royal Astronomical Society*,
*385*(3), 1279–1285. <https://doi.org/10.1111/j.1365-2966.2008.12689.x>

Bingham, C., Arbogast, B., Cornelissen Guillaume, G., Lee, J. K., &
Halberg, F. (1982). Inferential statistical methods for estimating and
comparing cosinor parameters. *Chronobiologia*, *9*(4), 397–439.

Chen, S., & Sun, X. (2024). Validating CircaCP: A generic sleep-wake
cycle detection algorithm for unlabelled actigraphy data. *Royal Society
Open Science*, *11*(5), 231468. <https://doi.org/10.1098/rsos.231468>

Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
(1992). Automatic sleep/wake identification from wrist activity.
*Sleep*, *15*(5), 461–469. <https://doi.org/10.1093/sleep/15.5.461>

Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical Biology
and Medical Modelling*, *11*, 16.
<https://doi.org/10.1186/1742-4682-11-16>

Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Multiscale entropy
analysis of complex physiologic time series. *Physical Review Letters*,
*89*(6), 068102. <https://doi.org/10.1103/PhysRevLett.89.068102>

Crespo, C., Aboy, M., Fernandez, J. R., & Mojon, A. (2012). Automatic
identification of activity-rest periods based on actigraphy. *Medical &
Biological Engineering & Computing*, *50*(4), 329–340.
<https://doi.org/10.1007/s11517-012-0875-y>

Hu, K., Ivanov, P. C., Chen, Z., Carpena, P., & Stanley, H. E. (2001).
Effect of trends on detrended fluctuation analysis. *Physical Review E*,
*64*(1), 011114. <https://doi.org/10.1103/PhysRevE.64.011114>

Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation
analysis of nonstationary time series. *Physica A: Statistical Mechanics
and Its Applications*, *316*(1-4), 87–114.
<https://doi.org/10.1016/S0378-4371(02)01383-3>

Kunsch, H. R. (1989). The jackknife and the bootstrap for general
stationary observations. *The Annals of Statistics*, *17*(3), 1217–1241.
<https://doi.org/10.1214/aos/1176347265>

Lim, A. S. P., Yu, L., Costa, M. D., Leurgans, S. E., Buchman, A. S.,
Bennett, D. A., & Saper, C. B. (2011). Quantification of the
fragmentation of rest-activity patterns in elderly individuals using a
state transition analysis. *Sleep*, *34*(11), 1569–1581.
<https://doi.org/10.5665/sleep.1400>

Liu, Y., & Xie, J. (2020). Cauchy combination test: A powerful test with
analytic p-value calculation under arbitrary dependency structures.
*Journal of the American Statistical Association*, *115*(529), 393–402.
<https://doi.org/10.1080/01621459.2018.1554485>

Lomb, N. R. (1976). Least-squares frequency analysis of unequally spaced
data. *Astrophysics and Space Science*, *39*(2), 447–462.
<https://doi.org/10.1007/BF00648343>

Nelson, W., Tong, Y. L., Lee, J. K., & Halberg, F. (1979). Methods for
cosinor-rhythmometry. *Chronobiologia*, *6*(4), 305–323.

Peng, C. K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E., &
Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
*Physical Review E*, *49*(2), 1685–1689.
<https://doi.org/10.1103/PhysRevE.49.1685>

Phillips, A. J. K., Clerx, W. M., O’Brien, C. S., Sano, A., Barger, L.
K., Picard, R. W., Lockley, S. W., Klerman, E. B., & Czeisler, C. A.
(2017). Irregular sleep/wake patterns are associated with poorer
academic performance and delayed circadian and sleep/wake timing.
*Scientific Reports*, *7*(1), 3216.
<https://doi.org/10.1038/s41598-017-03171-4>

Politis, D. N., & Romano, J. P. (1992). A circular block-resampling
procedure for stationary data. In R. LePage & L. Billard (Eds.),
*Exploring the limits of bootstrap* (pp. 263–270). Wiley.

Roenneberg, T., Keller, L. K., Fischer, D., Matera, J. L., Vetter, C., &
Winnebeck, E. C. (2015). Human activity and rest in situ. *Methods in
Enzymology*, *552*, 257–283.
<https://doi.org/10.1016/bs.mie.2014.11.028>

Scargle, J. D. (1982). Studies in astronomical time series analysis. II.
Statistical aspects of spectral analysis of unevenly spaced data. *The
Astrophysical Journal*, *263*, 835–853. <https://doi.org/10.1086/160554>

Sokolove, P. G., & Bushell, W. N. (1978). The chi square periodogram:
Its utility for analysis of circadian rhythms. *Journal of Theoretical
Biology*, *72*(1), 131–160.
<https://doi.org/10.1016/0022-5193(78)90022-X>

Van Someren, E. J. W., Swaab, D. F., Colenda, C. C., Cohen, W., McCall,
W. V., & Rosenquist, P. B. (1999). Bright light therapy: Improved
sensitivity to its effects on rest-activity rhythms in Alzheimer
patients by application of nonparametric methods. *Chronobiology
International*, *16*(4), 505–518.
<https://doi.org/10.3109/07420529908998724>

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and ultradian structure of human sleep in real life. *Current
Biology*, *28*(1), 49–59. <https://doi.org/10.1016/j.cub.2017.11.063>

Witting, W., Kwa, I. H., Eikelenboom, P., Mirmiran, M., & Swaab, D. F.
(1990). Alterations in the circadian rest-activity rhythm in aging and
Alzheimer’s disease. *Biological Psychiatry*, *27*(6), 563–572.
<https://doi.org/10.1016/0006-3223(90)90523-5>

Wittmann, M., Dinich, J., Merrow, M., & Roenneberg, T. (2006). Social
jetlag: Misalignment of biological and social time. *Chronobiology
International*, *23*(1-2), 497–509.
<https://doi.org/10.1080/07420520500545979>
