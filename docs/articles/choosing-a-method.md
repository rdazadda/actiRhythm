# Choosing a method

actiRhythm gives you many functions, and several can look like they
answer the same question. This guide maps the question you have to the
function that answers it, and pulls apart the methods that overlap. The
rule throughout: start simple, and escalate only when the data demand
it. The actogram and the nonparametric summary answer most questions on
their own; the rest of the package is there for the recordings that are
not so tidy.

## Start here, every time

Plot the recording first with
[`plot_actogram()`](https://rdazadda.github.io/actiRhythm/reference/plot_actogram.md),
then run
[`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md),
[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md),
and
[`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
(or
[`consensus.rhythmicity()`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
for all of them at once). If the active band stacks into a straight
vertical column and the numbers sit comfortably, you are done. Reach
further only when the band bends, scatters, or fragments, or when a
specific question below is yours.

## Is there a rhythm, and how strong?

Two complementary answers:

- **Without assuming a shape**:
  [`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
  returns interdaily stability (IS, how tightly the pattern repeats),
  intradaily variability (IV, how fragmented it is), and relative
  amplitude (RA, the day-night contrast) ([Van Someren et al.,
  1999](#ref-vansomeren1999); [Witting et al., 1990](#ref-witting1990)).
- **With a fitted cosine**:
  [`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  gives an interpretable MESOR, amplitude, and acrophase, and
  [`rhythmicity.test()`](https://rdazadda.github.io/actiRhythm/reference/rhythmicity.test.md)
  is the zero-amplitude F-test of whether that rhythm is real
  ([Cornelissen, 2014](#ref-cornelissen2014); [Nelson et al.,
  1979](#ref-nelson1979)).

Use the nonparametric route when you do not want to assume a sinusoid,
and the cosinor route when you want interpretable parameters. They
complement rather than compete: a high RA with a low IS, for instance,
is a strong rhythm carried on irregular timing, something neither number
says alone. To pool every test into one verdict,
[`consensus.rhythmicity()`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
runs the cosinor F-test, the Bingham ellipse, the Lomb-Scargle
false-alarm probability, and the chi-square periodogram, and reports a
majority vote and a Fisher-combined p-value ([Fisher,
1925](#ref-fisher1925)).

## What is the period, and is it really 24 hours?

A cosinor assumes a fixed period (24 h by default). If the rhythm might
run long or short (a free-running or drifting period), estimate it
instead of assuming it:

- [`circadian.period()`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md):
  a Lomb-Scargle periodogram; tolerates missing epochs and uneven
  sampling ([Lomb, 1976](#ref-lomb1976); [Ruf, 1999](#ref-ruf1999)).
- [`chi.sq.periodogram()`](https://rdazadda.github.io/actiRhythm/reference/chi.sq.periodogram.md):
  the classic Sokolove-Bushell periodogram; expects near-regular
  sampling ([Sokolove & Bushell, 1978](#ref-sokolove1978)).
- [`period.ci()`](https://rdazadda.github.io/actiRhythm/reference/period.ci.md):
  a circular block-bootstrap confidence interval on the period ([Politis
  & Romano, 1992](#ref-politis1992)).

**Lomb-Scargle vs chi-square:** prefer Lomb-Scargle when the recording
has gaps or irregular epochs; the chi-square periodogram is the
traditional choice for regularly sampled actigraphy. Running both and
seeing the peak land in the same place is reassurance that it is real
and not an artifact of one method.
[`consensus.rhythmicity()`](https://rdazadda.github.io/actiRhythm/reference/consensus.rhythmicity.md)
uses both for exactly this reason.

## How is the rhythm shaped?

In order of increasing waveform flexibility, choose the simplest the
percent-rhythm says is enough:

- [`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md):
  one symmetric cosine; simplest and most interpretable.
- [`cosinor.extended()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.extended.md)
  /
  [`cosinor.antilogistic()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.antilogistic.md):
  adds nonlinear shape for the asymmetric, squared-off real day ([Marler
  et al., 2006](#ref-marler2006)).
- [`cosinor.multicomponent()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md):
  several harmonics, with the number chosen by AIC or BIC.
- [`circadian.flm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.flm.md):
  a functional model of the whole daily profile that follows the true
  morning rise and evening decline ([Wang et al.,
  2011](#ref-wang2011flm)).

The last two are demonstrated in the *Beyond the basics* article.

## Does the rhythm drift, shift, or fragment?

When the actogram band bends or breaks up, a single summary averages the
change away. Switch to the time-resolved tools:

- [`circadian.spectrogram()`](https://rdazadda.github.io/actiRhythm/reference/circadian.spectrogram.md):
  period and its strength as a surface over time.
- [`circadian.ssa()`](https://rdazadda.github.io/actiRhythm/reference/circadian.ssa.md):
  separates the series into trend, circadian, and noise ([Golyandina &
  Zhigljavsky, 2013](#ref-golyandina2013)).
- [`circadian.emd()`](https://rdazadda.github.io/actiRhythm/reference/circadian.emd.md)
  with
  [`hilbert.huang()`](https://rdazadda.github.io/actiRhythm/reference/hilbert.huang.md):
  an instantaneous, cycle-by-cycle period ([N. E. Huang et al.,
  1998](#ref-huang1998)).
- [`circadian.wavelet()`](https://rdazadda.github.io/actiRhythm/reference/circadian.wavelet.md):
  the full time-frequency power surface.
- [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md):
  dates a regime change directly from the counts.
- [`curve.registration()`](https://rdazadda.github.io/actiRhythm/reference/curve.registration.md):
  a chronotype phase that does not blur when the days disagree ([Krafty
  et al., 2019](#ref-krafty2019)).

These have their own walkthrough in the *Nonstationary and complex
rhythms* article.

## How complex or fractal is the signal?

Two recordings can share the same IS and amplitude yet differ in their
moment-to-moment dynamics:

- [`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md):
  the detrended-fluctuation exponent (long-range correlation) ([Peng et
  al., 1994](#ref-peng1994)).
- [`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md):
  multifractality, summarised by the spectrum width ([Kantelhardt et
  al., 2002](#ref-kantelhardt2002)).
- [`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md):
  how complexity changes across coarse-graining scales ([Costa et al.,
  2002](#ref-costa2002)).

Demonstrated in *Beyond the basics*.

## When does the person sleep, and how regular?

- **A per-epoch sleep/wake label from counts**:
  [`sleep.cole.kripke()`](https://rdazadda.github.io/actiRhythm/reference/sleep.cole.kripke.md)
  for adults ([Cole et al., 1992](#ref-cole1992)),
  [`sleep.sadeh()`](https://rdazadda.github.io/actiRhythm/reference/sleep.sadeh.md)
  for children and adolescents ([Sadeh et al., 1994](#ref-sadeh1994)).
- **Regularity of that label day to day**:
  [`sleep.regularity.index()`](https://rdazadda.github.io/actiRhythm/reference/sleep.regularity.index.md)
  ([Phillips et al., 2017](#ref-phillips2017)), with
  [`sri.matrix()`](https://rdazadda.github.io/actiRhythm/reference/sri.matrix.md)
  for the day-by-day heatmap.
- **Each night’s onset and wake, straight from the counts**:
  [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
  ([Chen & Sun, 2024](#ref-chensun2024)).
- **Every rest bout, naps included**:
  [`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
  (Roenneberg) ([Roenneberg et al., 2015](#ref-roenneberg2015)) and
  [`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
  (Crespo) ([Crespo et al., 2012](#ref-crespo2012)), two independent
  algorithms;
  [`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
  is a probabilistic state-space alternative ([Q. Huang et al.,
  2018](#ref-huang2018hmm)).
- **Timing against the social clock**:
  [`social.jet.lag()`](https://rdazadda.github.io/actiRhythm/reference/social.jet.lag.md)
  ([Wittmann et al., 2006](#ref-wittmann2006)); the ultradian structure
  within sleep,
  [`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md)
  ([Winnebeck et al., 2018](#ref-winnebeck2018)).

**Which rest detector?**
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
finds the *one main* rest bout of each night (sleep and wake timing).
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
and
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
consolidate *every* rest bout across the recording, daytime naps
included. Run both for an independent cross-check.
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
models rest probabilistically rather than as hard bouts.

## Counts or raw acceleration?

- You have count files (`.agd`) or pre-extracted counts: analyse them
  directly; everything above takes a count series.
- You have raw files (`.gt3x`, `.cwa`, `.bin`):
  [`read.raw()`](https://rdazadda.github.io/actiRhythm/reference/read.raw.md)
  and
  [`raw.metrics()`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)
  give ENMO, MAD, and the **z-angle**, the gravity/posture signal counts
  cannot carry. The z-angle drives diary-free, posture-based sleep
  detection
  ([`rest.spt()`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md),
  [`sib.vanhees()`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md))
  that needs no scored sleep at all. See the *From raw acceleration*
  article.

Cross-brand counts (converting Axivity or GENEActiv to
ActiGraph-equivalent counts) are an approximation ([Brond et al.,
2017](#ref-brond2017)): fine for the relative and normalised metrics
here (IS, IV, RA, the periodograms, SRI), not for absolute-count
comparisons or ActiGraph cut-points.

## One subject or a group?

Every function above describes one recording. For a study:

- **Group mean rhythm**:
  [`population.cosinor()`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)
  pools per-subject cosinor fits into a Bingham population-mean rhythm
  with confidence intervals ([Bingham et al., 1982](#ref-bingham1982)).
- **Two-group comparison**:
  [`cosinor.compare()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.compare.md)
  tests whether the rhythm differs between groups (an omnibus Hotelling
  T^2 plus per-parameter tests) ([Bingham et al.,
  1982](#ref-bingham1982)).
- **Covariates, repeated measures, or nesting**: run
  [`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
  over the folder to get one row of metrics per recording. Then feed
  those per-subject metrics (IS, IV, RA, MESOR, amplitude, acrophase,
  period) into a mixed model with `lme4` or `nlme` alongside your
  design. The package supplies the metrics; you supply the multilevel
  model.

## A quick map

| Your question | Reach for | Walkthrough |
|----|----|----|
| Is there a rhythm, and how strong? | `circadian.rhythm`, `cosinor.analysis` + `rhythmicity.test`, `consensus.rhythmicity` | Get started |
| What is the period? | `circadian.period`, `chi.sq.periodogram`, `period.ci` | Get started |
| What shape is the day? | `cosinor.extended`, `cosinor.multicomponent`, `circadian.flm` | Beyond the basics |
| Does it drift or fragment? | `circadian.spectrogram`, `circadian.ssa`, `circadian.emd`, `sleep.changepoints` | Nonstationary and complex rhythms |
| How fractal/complex? | `fractal.dfa`, `mfdfa`, `multiscale.entropy` | Beyond the basics |
| When is sleep, how regular? | `sleep.cole.kripke`/`sleep.sadeh`, `sleep.regularity.index`, `rest.periods`/`rest.crespo` | Get started |
| Raw files, posture, diary-free sleep | `read.raw`, `raw.metrics`, `rest.spt` + `sib.vanhees` | From raw acceleration |
| A group, not one subject | `population.cosinor`, `cosinor.compare`, then `lme4`/`nlme` | Get started |
| What does column *X* mean? | the metric tables | Output codebook |

## Assumptions to respect

These hold across the package and decide whether any of the above is
trustworthy:

- **Evenly spaced epochs and enough valid days.**
  [`circadian.rhythm()`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md)
  sets aside days with too little valid recording before computing; very
  short records give unstable IS/IV and unreliable periods.
- **Wear time.** Gate the counts on valid wear, or supply a mask; if you
  have none,
  [`detect.nonwear.choi()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.choi.md)
  /
  [`detect.nonwear.troiano()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.troiano.md)
  (counts) or
  [`detect.nonwear.raw()`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md)
  (raw) build one.
  [`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
  and
  [`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
  apply no wear filter themselves.
- **Epoch length.** It affects IV and the epoch sleep scorers, which
  were validated at 60-second epochs.

## References

Bingham, C., Arbogast, B., Cornelissen Guillaume, G., Lee, J. K., &
Halberg, F. (1982). Inferential statistical methods for estimating and
comparing cosinor parameters. *Chronobiologia*, *9*(4), 397–439.

Brond, J. C., Andersen, L. B., & Arvidsson, D. (2017). Generating
ActiGraph counts from raw acceleration recorded by an alternative
monitor. *Medicine & Science in Sports & Exercise*, *49*(11), 2351–2360.
<https://doi.org/10.1249/MSS.0000000000001344>

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

Fisher, R. A. (1925). *Statistical methods for research workers*.
Oliver; Boyd.

Golyandina, N., & Zhigljavsky, A. (2013). Singular spectrum analysis for
time series. *SpringerBriefs in Statistics*.
<https://doi.org/10.1007/978-3-642-34913-3>

Huang, N. E., Shen, Z., Long, S. R., Wu, M. C., Shih, H. H., Zheng, Q.,
Yen, N. C., Tung, C. C., & Liu, H. H. (1998). The empirical mode
decomposition and the Hilbert spectrum for nonlinear and non-stationary
time series analysis. *Proceedings of the Royal Society A*, *454*(1971),
903–995. <https://doi.org/10.1098/rspa.1998.0193>

Huang, Q., Cohen, D., Komarzynski, S., Li, X.-M., Innominato, P., Levi,
F., & Finkenstadt, B. (2018). Hidden Markov models for monitoring
circadian rhythmicity in telemetric activity data. *Journal of the Royal
Society Interface*, *15*(139), 20170885.
<https://doi.org/10.1098/rsif.2017.0885>

Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation
analysis of nonstationary time series. *Physica A: Statistical Mechanics
and Its Applications*, *316*(1-4), 87–114.
<https://doi.org/10.1016/S0378-4371(02)01383-3>

Krafty, R. T., Fu, H., Graves, J. L., Bruce, S. A., Hall, M. H., &
Smagula, S. F. (2019). Measuring variability in rest-activity rhythms
from actigraphy with application to characterizing symptoms of
depression. *Statistics in Biosciences*, *11*, 314–333.
<https://doi.org/10.1007/s12561-018-09230-2>

Lomb, N. R. (1976). Least-squares frequency analysis of unequally spaced
data. *Astrophysics and Space Science*, *39*(2), 447–462.
<https://doi.org/10.1007/BF00648343>

Marler, M. R., Gehrman, P., Martin, J. L., & Ancoli-Israel, S. (2006).
The sigmoidally transformed cosine curve: A mathematical model for
circadian rhythms with symmetric non-sinusoidal shapes. *Statistics in
Medicine*, *25*(22), 3893–3904. <https://doi.org/10.1002/sim.2466>

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

Ruf, T. (1999). The Lomb-Scargle periodogram in biological rhythm
research: Analysis of incomplete and unequally spaced time-series.
*Biological Rhythm Research*, *30*(2), 178–201.
<https://doi.org/10.1076/brhm.30.2.178.1422>

Sadeh, A., Sharkey, K. M., & Carskadon, M. A. (1994). Activity-based
sleep-wake identification: An empirical test of methodological issues.
*Sleep*, *17*(3), 201–207. <https://doi.org/10.1093/sleep/17.3.201>

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

Wang, J., Xian, H., Licis, A., Deych, E., Ding, J., McLeland, J.,
Toedebusch, C., Li, T., Duntley, S., & Shannon, W. (2011). Measuring the
impact of apnea and obesity on circadian activity patterns using
functional linear modeling of actigraphy data. *Journal of Circadian
Rhythms*, *9*, 11. <https://doi.org/10.1186/1740-3391-9-11>

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
