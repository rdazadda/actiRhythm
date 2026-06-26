# Package index

## Reading & preparing data

Read ActiGraph counts and raw multi-brand accelerometer files, compute
and calibrate raw metrics, detect non-wear, and access the bundled
example data.

- [`read.agd()`](https://rdazadda.github.io/actiRhythm/reference/read.agd.md)
  : Read ActiGraph .agd File
- [`agd.counts()`](https://rdazadda.github.io/actiRhythm/reference/agd.counts.md)
  : Extract Counts from .agd Data
- [`write.agd()`](https://rdazadda.github.io/actiRhythm/reference/write.agd.md)
  : Write Activity Counts to an .agd File
- [`gt3x.to.agd()`](https://rdazadda.github.io/actiRhythm/reference/gt3x.to.agd.md)
  : Convert a Raw .gt3x File to an .agd File
- [`has.inclinometer()`](https://rdazadda.github.io/actiRhythm/reference/has.inclinometer.md)
  : Check if AGD Data Has Inclinometer Information
- [`read.actigraph.csv()`](https://rdazadda.github.io/actiRhythm/reference/read.actigraph.csv.md)
  : Read an ActiGraph (ActiLife) Epoch CSV
- [`counts.from.data.frame()`](https://rdazadda.github.io/actiRhythm/reference/counts.from.data.frame.md)
  : Activity Counts from a Data Frame
- [`read.raw()`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
  : Activity Counts from a Raw Accelerometer File (Any Supported Brand)
- [`gt3x.counts()`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md)
  : Activity Counts from a Raw .gt3x File
- [`axivity.counts()`](https://rdazadda.github.io/actiRhythm/reference/axivity.counts.md)
  : Activity Counts from a Raw Axivity .cwa File
- [`geneactiv.counts()`](https://rdazadda.github.io/actiRhythm/reference/geneactiv.counts.md)
  : Activity Counts from a Raw GENEActiv .bin File
- [`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
  : Raw-Acceleration Epoch Metrics (ENMO, MAD, z-Angle)
- [`auto.calibrate()`](https://rdazadda.github.io/actiRhythm/reference/auto.calibrate.md)
  : Auto-Calibrate Raw Acceleration to the Unit Gravity Sphere
- [`detect.nonwear.choi()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.choi.md)
  : Choi (2011) Non-Wear Detection
- [`detect.nonwear.troiano()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.troiano.md)
  : Troiano (2008) Non-Wear Detection
- [`detect.nonwear.raw()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
  : Non-Wear Detection from Raw Acceleration
- [`example_agd()`](https://rdazadda.github.io/actiRhythm/reference/example_agd.md)
  : Get Path to Example AGD Files
- [`example_raw()`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)
  : Synthetic Raw Acceleration Recording

## Nonparametric rhythm metrics

Shape-free descriptors of the rest-activity rhythm: IS, IV, RA, L5/M10,
and fragmentation, with no assumed waveform.

- [`circadian-analysis`](https://rdazadda.github.io/actiRhythm/reference/circadian-analysis.md)
  : Circadian Rhythm Analysis
- [`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
  : Circadian Rhythm Analysis
- [`intradaily.variability.multiscale()`](https://rdazadda.github.io/actiRhythm/reference/intradaily.variability.multiscale.md)
  : Multiscale Intradaily Variability (IVm)
- [`circadian.is.multiscale()`](https://rdazadda.github.io/actiRhythm/reference/circadian.is.multiscale.md)
  : Multi-resolution Interdaily Stability (IS)
- [`dichotomy.index()`](https://rdazadda.github.io/actiRhythm/reference/dichotomy.index.md)
  : Dichotomy Index (I \< O)
- [`circadian.quotient()`](https://rdazadda.github.io/actiRhythm/reference/circadian.quotient.md)
  : Circadian Quotient and Relative Amplitude from a Cosinor Fit
- [`activity.balance.index()`](https://rdazadda.github.io/actiRhythm/reference/activity.balance.index.md)
  : Activity Balance Index
- [`rest.activity.fragmentation()`](https://rdazadda.github.io/actiRhythm/reference/rest.activity.fragmentation.md)
  : Rest-Activity Bout Fragmentation
- [`transition.probability()`](https://rdazadda.github.io/actiRhythm/reference/transition.probability.md)
  : Rest-Active Transition Probabilities (ASTP / SATP)
- [`state.transitions()`](https://rdazadda.github.io/actiRhythm/reference/state.transitions.md)
  : Rest-Activity State Transition Rates (kRA, kAR)
- [`activity.extrema()`](https://rdazadda.github.io/actiRhythm/reference/activity.extrema.md)
  : Generalized Least- and Most-Active Periods (MX / LX)

## Cosinor models

Sinusoidal rhythm models from the single cosinor and its rhythmicity
test through extended, multicomponent, population-mean, and
group-comparison fits.

- [`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  : Cosinor Analysis for Circadian Rhythm
- [`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
  : Cosinor Rhythmicity Test
- [`cosinor.confidence.ellipse()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.confidence.ellipse.md)
  : Joint Amplitude-Acrophase Confidence Ellipse (Bingham et al. 1982)
- [`cosinor.extended()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.extended.md)
  : Multi-Component Cosinor Analysis with Harmonics
- [`cosinor.antilogistic()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.antilogistic.md)
  : Anti-Logistic (Extended) Cosinor Analysis (Marler et al. 2006)
- [`cosinor.multicomponent()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md)
  : Multicomponent Cosinor with Model Selection
- [`population.cosinor()`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)
  : Population-Mean Cosinor (Bingham)
- [`cosinor.compare()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.compare.md)
  : Compare Cosinor Rhythms Between Two Groups
- [`consensus.rhythmicity()`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
  : Consensus Rhythmicity Across Methods

## Period & spectral analysis

Estimate the free-running period and track how it drifts: Lomb-Scargle
and chi-square periodograms, a bootstrap CI, a spectrogram, and
ultradian power.

- [`circadian-period`](https://rdazadda.github.io/actiRhythm/reference/circadian-period.md)
  : Endogenous Circadian Period Estimation (Lomb-Scargle)
- [`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md)
  : Estimate Endogenous Circadian Period via the Lomb-Scargle
  Periodogram
- [`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md)
  : Confidence Interval for the Endogenous Circadian Period
- [`chi.sq.periodogram()`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md)
  : Chi-square (Sokolove-Bushell) Periodogram
- [`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md)
  : Circadian Spectrogram (Period over Time)
- [`ultradian.bandpower()`](https://rdazadda.github.io/actiRhythm/reference/ultradian.bandpower.md)
  : Ultradian Wavelet Band Power
- [`residual.spectrum()`](https://rdazadda.github.io/actiRhythm/reference/residual.spectrum.md)
  : Residual Circadian Spectrum

## Fractal & nonlinear dynamics

Long-range correlation and complexity of the activity series.

- [`circadian-fractal`](https://rdazadda.github.io/actiRhythm/reference/circadian-fractal.md)
  : Fractal and Complexity Metrics for Activity Time Series
- [`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
  : Detrended Fluctuation Analysis (DFA)
- [`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md)
  : Multifractal Detrended Fluctuation Analysis (MF-DFA)
- [`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
  : Multiscale Sample Entropy (MSE)

## Decomposition & time-frequency

Data-adaptive decompositions of the recording into trend, rhythm, and
noise.

- [`circadian.wavelet()`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md)
  : Continuous Wavelet Transform of the Activity Rhythm
- [`wavelet.coi()`](https://rdazadda.github.io/actiRhythm/reference/wavelet.coi.md)
  : Wavelet Cone of Influence
- [`circadian.emd()`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)
  : Empirical Mode Decomposition of the Activity Rhythm
- [`hilbert.huang()`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md)
  : Hilbert-Huang Instantaneous Phase and Frequency
- [`circadian.ssa()`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md)
  : Singular Spectrum Analysis of an Activity Series
- [`circadian.flm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md)
  : Functional Linear Model of the 24-Hour Activity Profile
- [`curve.registration()`](https://rdazadda.github.io/actiRhythm/reference/curve.registration.md)
  : Curve Registration of Daily Activity Profiles

## Sleep & rest detection

Epoch sleep/wake scoring, rest-period and sleep-period-time detection,
and locomotor inactivity during sleep.

- [`sleep.cole.kripke()`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)
  : Cole-Kripke Sleep/Wake Scoring
- [`sleep.sadeh()`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md)
  : Sadeh Sleep/Wake Scoring
- [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
  : Change-Point Detection of Sleep and Wake Onsets
- [`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
  : Consolidated Rest-Period Detection (Roenneberg / MASDA)
- [`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
  : Crespo Rest/Activity Period Detection
- [`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
  : State-Space (Hidden Markov) Rest-Activity Model
- [`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md)
  : Sleep-Period-Time Window from the z-Angle (HDCZA)
- [`sib.vanhees()`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md)
  : Sustained-Inactivity-Bout Sleep Scoring from the z-Angle
- [`sleep.from.spt()`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md)
  : Sleep Parameters from SPT and Sustained-Inactivity Bouts
- [`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md) :
  Locomotor Inactivity During Sleep (LIDS)

## Phase, onset & regularity

Sleep timing, regularity, and day-to-day phase stability.

- [`sleep.regularity.index()`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
  : Calculate Sleep Regularity Index (SRI) - Exported Version
- [`sri.matrix()`](https://rdazadda.github.io/actiRhythm/reference/sri.matrix.md)
  : Sleep Regularity Index (SRI) - Phillips (2017) Epoch-of-Day x Day
  Matrix
- [`social.jet.lag()`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md)
  : Calculate Social Jet Lag
- [`composite.phase.deviation()`](https://rdazadda.github.io/actiRhythm/reference/composite.phase.deviation.md)
  : Composite Phase Deviation (Fischer & Roenneberg, 2016)
- [`phase.concentration()`](https://rdazadda.github.io/actiRhythm/reference/phase.concentration.md)
  : Phase Concentration Tests
- [`activity.onset.offset()`](https://rdazadda.github.io/actiRhythm/reference/activity.onset.offset.md)
  : Activity Onset and Offset (Relative-Difference Phase Markers)
- [`circadian.onset.ci()`](https://rdazadda.github.io/actiRhythm/reference/circadian.onset.ci.md)
  : Confidence Intervals for L5/M10 Onset Timing

## Batch & reporting workflow

Run the pipeline over many recordings and export full analyses.

- [`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
  : Batch Circadian Analysis of Multiple AGD Files
- [`circadian.daily()`](https://rdazadda.github.io/actiRhythm/reference/circadian.daily.md)
  : Per-Day Nonparametric Metrics
- [`circadian.raw()`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md)
  : Circadian Analysis Directly from Raw Acceleration
- [`circadian.workbook()`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md)
  : Export a Full Circadian Analysis to an Excel Workbook
- [`save.circadian.plot()`](https://rdazadda.github.io/actiRhythm/reference/save.circadian.plot.md)
  : Save a actiRhythm Plot to a File

## Plotting & theme

ggplot2 plots and the package theme/palette.

- [`circadian-plots`](https://rdazadda.github.io/actiRhythm/reference/circadian-plots.md)
  : Circadian Rhythm Visualizations
- [`plot_actogram()`](https://rdazadda.github.io/actiRhythm/reference/plot_actogram.md)
  : Double-Plotted Actogram
- [`plot_periodogram()`](https://rdazadda.github.io/actiRhythm/reference/plot_periodogram.md)
  : Plot the Lomb-Scargle Periodogram with the Endogenous Period
- [`plot_chisq()`](https://rdazadda.github.io/actiRhythm/reference/plot_chisq.md)
  : Plot the Chi-Square (Sokolove-Bushell) Periodogram
- [`plot_dfa()`](https://rdazadda.github.io/actiRhythm/reference/plot_dfa.md)
  : Plot the Detrended Fluctuation Analysis (DFA) Scaling Relationship
- [`plot_extended_cosinor()`](https://rdazadda.github.io/actiRhythm/reference/plot_extended_cosinor.md)
  : Plot the Extended (Marler) Cosinor Fit on the 24-Hour Activity
  Profile
- [`theme_actiRhythm()`](https://rdazadda.github.io/actiRhythm/reference/theme_actiRhythm.md)
  : actiRhythm ggplot2 Theme
- [`set_actiRhythm_theme()`](https://rdazadda.github.io/actiRhythm/reference/set_actiRhythm_theme.md)
  : Set actiRhythm Theme as Default
- [`scale_color_actiRhythm()`](https://rdazadda.github.io/actiRhythm/reference/scale_color_actiRhythm.md)
  : ggplot2 Scale for actiRhythm Colors
- [`scale_fill_actiRhythm()`](https://rdazadda.github.io/actiRhythm/reference/scale_fill_actiRhythm.md)
  : ggplot2 Fill Scale for actiRhythm Colors
- [`actiRhythm_colors()`](https://rdazadda.github.io/actiRhythm/reference/actiRhythm_colors.md)
  : actiRhythm Color Generator

## C++ backend

The compiled fast path behind the nonparametric metrics, and its
controls.

- [`cpp_available()`](https://rdazadda.github.io/actiRhythm/reference/cpp_available.md)
  : Check C++ Availability
- [`backend_info()`](https://rdazadda.github.io/actiRhythm/reference/backend_info.md)
  : Backend Information
- [`circadian_cpp()`](https://rdazadda.github.io/actiRhythm/reference/circadian_cpp.md)
  : Complete Circadian Analysis
- [`IS_cpp()`](https://rdazadda.github.io/actiRhythm/reference/IS_cpp.md)
  : Interdaily Stability (IS)
- [`IV_cpp()`](https://rdazadda.github.io/actiRhythm/reference/IV_cpp.md)
  : Intradaily Variability (IV)
- [`L5M10_cpp()`](https://rdazadda.github.io/actiRhythm/reference/L5M10_cpp.md)
  : Calculate L5/M10 Circadian Metrics
- [`rolling_mean()`](https://rdazadda.github.io/actiRhythm/reference/rolling_mean.md)
  : Rolling Mean
- [`rolling_sd()`](https://rdazadda.github.io/actiRhythm/reference/rolling_sd.md)
  : Rolling Standard Deviation
- [`rolling_sum()`](https://rdazadda.github.io/actiRhythm/reference/rolling_sum.md)
  : Rolling Sum

## Print & plot methods

S3 methods dispatched automatically when you print or plot a result
object.

- [`print(`*`<actiRhythm_circadian>`*`)`](https://rdazadda.github.io/actiRhythm/reference/print.actiRhythm_circadian.md)
  : Print Method for Circadian Rhythm Results
- [`print(`*`<actiRhythm_cosinor>`*`)`](https://rdazadda.github.io/actiRhythm/reference/print.actiRhythm_cosinor.md)
  : Print method for cosinor analysis
- [`print(`*`<actiRhythm_cosinor_extended>`*`)`](https://rdazadda.github.io/actiRhythm/reference/print.actiRhythm_cosinor_extended.md)
  : Print method for extended cosinor analysis
- [`print(`*`<actiRhythm_dfa>`*`)`](https://rdazadda.github.io/actiRhythm/reference/print.actiRhythm_dfa.md)
  : Print Method for DFA Results
- [`print(`*`<actiRhythm_mse>`*`)`](https://rdazadda.github.io/actiRhythm/reference/print.actiRhythm_mse.md)
  : Print Method for Multiscale Entropy Results
- [`plot(`*`<actiRhythm_circadian>`*`)`](https://rdazadda.github.io/actiRhythm/reference/plot.actiRhythm_circadian.md)
  : Plot Circadian Rhythm Profile
- [`plot_actogram()`](https://rdazadda.github.io/actiRhythm/reference/plot_actogram.md)
  : Double-Plotted Actogram
- [`plot_chisq()`](https://rdazadda.github.io/actiRhythm/reference/plot_chisq.md)
  : Plot the Chi-Square (Sokolove-Bushell) Periodogram
- [`plot_dfa()`](https://rdazadda.github.io/actiRhythm/reference/plot_dfa.md)
  : Plot the Detrended Fluctuation Analysis (DFA) Scaling Relationship
- [`plot_extended_cosinor()`](https://rdazadda.github.io/actiRhythm/reference/plot_extended_cosinor.md)
  : Plot the Extended (Marler) Cosinor Fit on the 24-Hour Activity
  Profile
- [`plot_periodogram()`](https://rdazadda.github.io/actiRhythm/reference/plot_periodogram.md)
  : Plot the Lomb-Scargle Periodogram with the Endogenous Period
