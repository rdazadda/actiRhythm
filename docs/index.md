# actiRhythm

**actiRhythm measures circadian rest-activity rhythms from activity
counts.** It needs only a count vector and its timestamps, and it reads
those counts straight out of an ActiGraph `.agd` file. From one
recording it produces the metrics a chronobiology analysis reports, from
the nonparametric measures through cosinor, periodograms, and fractal
structure. Every analysis returns a typed object that prints its own
metrics and plots with one call.

You can usually see the answer before computing it. Plot the recording
as an actogram and a regular sleeper forms a single vertical band of
activity at the same clock time each day, while a drifting rhythm
scatters across the panel. The methods that follow put numbers to what
that picture shows, and each carries the reference that defined it, from
Witting (1990) and van Someren (1999) through the cosinor and the
periodogram.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("rdazadda/actiRhythm")
```

actiRhythm needs R 4.1 or newer and a C++17 compiler (Rtools on Windows,
the Xcode command-line tools on macOS) for the small backend behind its
nonparametric metrics.

## A first analysis

One ActiGraph recording ships with the package, so this runs without any
data of your own. The bundled file is a single de-identified
60-second-epoch recording, included only so the examples run.

``` r

library(actiRhythm)

# Read the bundled .agd and pull the count series
agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))

# Describe the rest-activity rhythm without assuming a waveform
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
#>   Intradaily Variability (IV): 1.0008 (~0=sine, ~2=noise)
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

The actogram is usually the first thing worth plotting.

``` r

plot_actogram(agd$axis1, agd$timestamp)
```

![Double-plotted actogram: each row is one day shown twice across 48
hours, so a stable rhythm forms a vertical band of activity at the same
clock time each day, while scattered fill marks a fragmented or shifting
rhythm.](reference/figures/README-actogram-1.png)

Double-plotted actogram: each row is one day shown twice across 48
hours, so a stable rhythm forms a vertical band of activity at the same
clock time each day, while scattered fill marks a fragmented or shifting
rhythm.

For this recording a high relative amplitude near 0.98 sits with a low
interdaily stability near 0.23: the days are strongly active, but the
pattern does not land at the same clock time from one day to the next.

Every number and figure here regenerates from the installed package and
its bundled recording, and the calculations behind them are exercised by
the package’s test suite.

## Going further

The same recording runs through the rest of the package.
[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
fits a 24-hour cosine (Cornelissen 2014) and
[`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
checks whether the rhythm is statistically real.
[`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
estimates the free-running period from a Lomb-Scargle periodogram (Lomb
1976),
[`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md)
puts a bootstrap confidence interval on it, and
[`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
shows how that period drifts across the recording.
[`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
measures the long-range correlation in the series (Peng 1994), and
[`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
and
[`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md)
describe its finer nonlinear structure.
[`circadian.flm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md)
fits a functional (multi-harmonic) model of the daily profile and
[`circadian.ssa()`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)
decomposes the recording into trend, circadian, and noise components.
[`sleep.regularity.index()`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
(Phillips 2017),
[`social.jet.lag()`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md),
and
[`state.transitions()`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
summarise sleep timing and rest-activity fragmentation, and
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
locates each night’s sleep and wake onset directly from the counts (Chen
and Sun 2024).
[`sleep.cole.kripke()`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)
and
[`sleep.sadeh()`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md)
score each epoch as sleep or wake (Cole et al. 1992; Sadeh et al. 1994),
and
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
(Roenneberg et al. 2015) and
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
(Crespo et al. 2012) consolidate every rest bout across the recording,
naps included, by two independent algorithms.

A further set of methods covers the harder cases.
[`circadian.wavelet()`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)
gives the time-frequency power surface,
[`ultradian.bandpower()`](https://rdazadda.github.io/actiRhythm/reference/ultradian.bandpower.md)
the dyadic ultradian bands, and
[`circadian.emd()`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)
with
[`hilbert.huang()`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md)
a data-adaptive decomposition and instantaneous phase.
[`cosinor.multicomponent()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md)
fits a multi-harmonic shape,
[`activity.onset.offset()`](https://rdazadda.github.io/actiRhythm/reference/activity.onset.offset.md)
marks daily onset and offset, and
[`phase.concentration()`](https://rdazadda.github.io/actiRhythm/reference/phase.concentration.md)
tests day-to-day phase clustering.
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
is a state-space rest-activity model,
[`curve.registration()`](https://rdazadda.github.io/actiRhythm/reference/curve.registration.md)
aligns daily profiles for a scale-invariant chronotype phase, and
[`residual.spectrum()`](https://rdazadda.github.io/actiRhythm/reference/residual.spectrum.md)
spectra what the cosinor leaves behind.
[`circadian.daily()`](https://rdazadda.github.io/actiRhythm/reference/circadian.daily.md),
[`intradaily.variability.multiscale()`](https://rdazadda.github.io/actiRhythm/reference/intradaily.variability.multiscale.md),
[`activity.extrema()`](https://rdazadda.github.io/actiRhythm/reference/activity.extrema.md),
and
[`dichotomy.index()`](https://rdazadda.github.io/actiRhythm/reference/dichotomy.index.md)
complete the nonparametric set.

Counts can come from more than `.agd` files:
[`read.actigraph.csv()`](https://rdazadda.github.io/actiRhythm/reference/read.actigraph.csv.md)
reads ActiLife epoch CSVs and
[`counts.from.data.frame()`](https://rdazadda.github.io/actiRhythm/reference/counts.from.data.frame.md)
takes a count series from any data frame.
[`read.raw()`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
turns raw ActiGraph `.gt3x`, Axivity `.cwa`, and GENEActiv `.bin` files
into counts with the agcounts band-pass filter, and
[`gt3x.counts()`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md),
[`axivity.counts()`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md),
and
[`geneactiv.counts()`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md)
do the same per brand. Cross-brand counts are an approximation rather
than native ActiGraph output.
[`detect.nonwear.choi()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.choi.md)
and
[`detect.nonwear.troiano()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.troiano.md)
flag non-wear time directly from the counts to use as the wear filter
(Choi et al. 2011; Troiano et al. 2008).

Raw acceleration also gives the gravity-based metrics that counts cannot
carry.
[`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
returns per-epoch ENMO, MAD and the z-angle with van Hees (2014)
auto-calibration, and
[`circadian.raw()`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md)
runs every method above on ENMO from a single call. The z-angle also
supports a sleep pipeline that needs no diary, which counts cannot do.
[`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
finds the nightly sleep-period window (HDCZA, van Hees 2018) and
[`sib.vanhees()`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md)
scores sustained-inactivity bouts (van Hees 2015);
[`sleep.from.spt()`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)
turns them into onset, wake, WASO and efficiency.
[`detect.nonwear.raw()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
gates out a device taken off and left still, so it is not read as sleep
(van Hees 2011).
[`activity.balance.index()`](https://rdazadda.github.io/actiRhythm/reference/activity.balance.index.md)
and
[`transition.probability()`](https://rdazadda.github.io/actiRhythm/reference/transition.probability.md)
add two fragmentation estimators.

When a study has many files,
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
runs the whole pipeline over a folder and returns one row per recording,
and
[`circadian.workbook()`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)
writes a full analysis to a multi-sheet Excel file. Every plotting
function returns a `ggplot` object you can theme and save with
[`save.circadian.plot()`](https://rdazadda.github.io/actiRhythm/reference/save.circadian.plot.md).

## Documentation

[`vignette("actiRhythm")`](https://rdazadda.github.io/actiRhythm/articles/actiRhythm.md)
walks through a full analysis, from reading a file to exporting the
results. Each function’s help page carries its method reference, so
[`?circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
[`?cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
and
[`?circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
are good places to read more. Questions and bug reports go to the issue
tracker at <https://github.com/rdazadda/actiRhythm/issues>, and
contributions are welcome (see `CONTRIBUTING.md`).

## Citation

Run `citation("actiRhythm")` for the entry to use in published work.
When you publish results, report the version you ran, from
`packageVersion("actiRhythm")`, so the analysis stays reproducible
across releases.

## License

MIT, see [LICENSE](https://rdazadda.github.io/actiRhythm/LICENSE).

actiRhythm is developed and maintained by the Center for Alaska Native
Health Research at the University of Alaska Fairbanks.
