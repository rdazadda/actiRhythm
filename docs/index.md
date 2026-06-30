# actiRhythm

**actiRhythm turns activity counts or raw acceleration into a complete
picture of the circadian rest-activity rhythm.** It reads ActiGraph
`.agd` counts and raw `.gt3x`, `.cwa`, and `.bin` recordings,
auto-calibrates the raw signal, and derives the ENMO, MAD, and z-angle
metrics that agree with GGIR. From a single recording it computes what a
chronobiology analysis reports, from the nonparametric measures and
cosinor through periodograms and fractal structure to wavelet and
empirical-mode decomposition, sleep scoring, and a state-space
rest-activity model. Methods that usually live in separate packages, or
in no R package at all, sit behind one consistent interface, and every
analysis returns a typed object that prints its own summary and draws
its own plots.

Each method carries the reference that defined it, and where a trusted
implementation already exists, actiRhythm is checked against it: the raw
metrics against GGIR, the single-component cosinor against the `cosinor`
package, and the nonparametric measures against `nparACT` and `ActCR`.
Those cross-checks run live in the package’s [validation
article](https://rdazadda.github.io/actiRhythm/articles/validation.html).

You can usually see the answer before you compute it. Plot a recording
as an actogram and a regular sleeper forms a single vertical band of
activity at the same clock time each day, while a drifting rhythm
scatters across the panel. The numbers that follow put a figure to what
the picture shows, from Witting (1990) and van Someren (1999) through
Cornelissen (2014) and the periodogram.

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

## What the package covers

**Cosinor and the spectrum.**
[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
fits a 24-hour cosine (Cornelissen 2014) and
[`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
checks whether the rhythm is statistically real;
[`cosinor.multicomponent()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md)
adds harmonics and
[`population.cosinor()`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)
compares groups.
[`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
estimates the free-running period from a Lomb-Scargle periodogram (Lomb
1976) and
[`chi.sq.periodogram()`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md)
from the chi-square periodogram (Sokolove 1978);
[`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md)
puts a bootstrap interval on it, and
[`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
shows how the period drifts across the recording.
[`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
measures the long-range correlation in the series (Peng 1994), and
[`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md)
and
[`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
describe its finer nonlinear structure.

**Rhythms that move rather than repeat.** A set of nonstationary methods
reads the change directly.
[`circadian.wavelet()`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)
returns the time-frequency power surface with its cone of influence,
[`circadian.emd()`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)
with
[`hilbert.huang()`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md)
give a data-adaptive decomposition and an instantaneous period, and
[`circadian.ssa()`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)
separates the series into trend, circadian, and residual components.
[`circadian.flm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md)
fits a functional model of the daily profile,
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
infers rest and active states with a hidden Markov model,
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
locates each night’s onset and wake straight from the counts (Chen and
Sun 2024), and
[`curve.registration()`](https://rdazadda.github.io/actiRhythm/reference/curve.registration.md)
aligns the days on their active-phase landmark for a scale-invariant
phase marker.

**Sleep, regularity, and fragmentation.**
[`sleep.cole.kripke()`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)
and
[`sleep.sadeh()`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md)
score each epoch as sleep or wake (Cole et al. 1992; Sadeh et al. 1994).
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
(Roenneberg et al. 2015) returns every consolidated rest bout, naps
included, while
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
(Crespo et al. 2012) detects the main daily rest periods.
[`sleep.regularity.index()`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
(Phillips 2017),
[`social.jet.lag()`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md),
and [`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md)
(Winnebeck 2018) summarise sleep regularity, chronotype misalignment,
and the ultradian sleep cycle, and
[`phase.concentration()`](https://rdazadda.github.io/actiRhythm/reference/phase.concentration.md)
and
[`state.transitions()`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
test day-to-day phase clustering and rest-activity fragmentation.

**From files to counts to raw acceleration.** Counts can come from more
than `.agd` files:
[`read.actigraph.csv()`](https://rdazadda.github.io/actiRhythm/reference/read.actigraph.csv.md)
reads ActiLife epoch CSVs, and
[`read.raw()`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
turns raw `.gt3x`, `.cwa`, and `.bin` files into counts with the
agcounts implementation of the ActiGraph count algorithm (Neishabouri
2022);
[`gt3x.counts()`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md),
[`axivity.counts()`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md),
and
[`geneactiv.counts()`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md)
do the same per brand, with cross-brand counts treated as an
approximation rather than native ActiGraph output.
[`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
returns per-epoch ENMO, MAD, and the z-angle with van Hees (2014)
auto-calibration, and
[`circadian.raw()`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md)
runs every method above on ENMO from one call. The z-angle drives a
diary-free sleep pipeline that counts cannot:
[`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
finds the nightly sleep-period window (HDCZA, van Hees 2018),
[`sib.vanhees()`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md)
scores sustained-inactivity bouts (van Hees 2015), and
[`sleep.from.spt()`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)
reports onset, wake, WASO, and efficiency.
[`detect.nonwear.choi()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.choi.md),
[`detect.nonwear.troiano()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.troiano.md),
and
[`detect.nonwear.raw()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
flag non-wear time (Choi et al. 2011; Troiano et al. 2008; van Hees et
al. 2013).

**Plots and batch reports.** Every analysis prints a readable summary
and draws a themed `ggplot`: the actogram and periodogram, the
spectrogram, the cosinor fit and confidence ellipse, the wavelet
scalogram, the singular-spectrum components, a circular phase rose, and
a rest-detector comparison strip, each one you can recolour with
[`theme_actiRhythm()`](https://rdazadda.github.io/actiRhythm/reference/theme_actiRhythm.md)
and save with
[`save.circadian.plot()`](https://rdazadda.github.io/actiRhythm/reference/save.circadian.plot.md).
When a study has many files,
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
runs the whole pipeline over a folder and returns one row per recording,
and
[`circadian.workbook()`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)
writes the full analysis to a multi-sheet Excel file.

## Documentation

The package website, <https://rdazadda.github.io/actiRhythm/>, collects
the full documentation.
[`vignette("actiRhythm")`](https://rdazadda.github.io/actiRhythm/articles/actiRhythm.md)
is the get-started walkthrough and
[`vignette("output-codebook")`](https://rdazadda.github.io/actiRhythm/articles/output-codebook.md)
is a data dictionary for every metric the package reports. A set of
articles goes deeper: [Choosing a
method](https://rdazadda.github.io/actiRhythm/articles/choosing-a-method.html)
maps each question to the right function, [Validation against reference
packages](https://rdazadda.github.io/actiRhythm/articles/validation.html)
cross-checks the metrics live against GGIR, `nparACT`, `ActCR`, and the
`cosinor` package, and the method-family articles (nonparametric,
cosinor, period, fractal, nonstationary, phase and regularity, sleep,
raw acceleration, batch) each work one family from the math to a worked
example. Every function’s help page also carries its method reference,
so
[`?circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
[`?cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
and
[`?circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
are good places to read more. Questions and bug reports go to the issue
tracker at <https://github.com/rdazadda/actiRhythm/issues>, and
contributions are welcome.

## Citation

Run `citation("actiRhythm")` for the entry to use in published work.
When you publish results, report the version you ran, from
`packageVersion("actiRhythm")`, so the analysis stays reproducible
across releases.

## License

MIT, see [LICENSE](https://rdazadda.github.io/actiRhythm/LICENSE).

actiRhythm is developed and maintained by the Center for Alaska Native
Health Research at the University of Alaska Fairbanks.
