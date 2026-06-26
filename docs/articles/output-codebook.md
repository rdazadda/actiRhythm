# Output codebook (metrics data dictionary)

This is the lookup for every metric `actiRhythm` emits in its one-row
summary – the columns of
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md),
the **Summary** sheet of
[`circadian.workbook()`](https://rdazadda.github.io/actiRhythm/reference/circadian.workbook.md),
and the slots of the individual result objects. Find any output here:
what it is, its formula or units, the range you should expect, the
defining reference, and where it lives in the result object.

The table below, the workbook’s *Data Dictionary* sheet, and the
package’s tests all read one file –
`inst/extdata/metrics_dictionary.csv` – so the codebook can never drift
from the code (a test checks that the documented metrics are exactly
those
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
returns). Individual result objects carry further fields beyond this
summary; see each function’s help page (for example
[`?circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md))
for the complete object.

## Reading this codebook

Each metric is named exactly as it appears as a column of the
[`circadian.batch()`](https://rdazadda.github.io/actiRhythm/reference/circadian.batch.md)
data frame and the workbook **Summary** sheet. The **Output object**
column shows where the same value lives when you call the producing
function directly – `circadian.rhythm()$IS`, `cosinor.analysis()$mesor`,
and so on. The **Reference** column gives the citation key for the
method; full references are listed at the end. A metric whose **Range**
notes “requires `sleep_state`” or “requires `sleep_periods`” is `NA`
unless that input is supplied.

## Recording descriptors

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| n_epochs | Number of epochs analysed in the recording | count | non-negative integer |  | summary$`n_epochs       |
|epoch_length_s |Epoch length of the recording              |seconds         |typically 15-60 s      |          |summary`$epoch_length_s |

## Nonparametric metrics

Shape-free descriptors of the rest-activity rhythm, computed with no
assumed waveform ([Van Someren et al., 1999](#ref-vansomeren1999);
[Witting et al., 1990](#ref-witting1990)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| L5 | Mean activity over the least-active 5 hours of the average day | activity counts per min | non-negative; lower = deeper rest | vansomeren1999 | circadian.rhythm()$`L5  |
|M10    |Mean activity over the most-active 10 hours of the average day           |activity counts per min                                          |non-negative; higher = stronger active phase     |vansomeren1999 |circadian.rhythm()`$M10 |
| RA | Relative amplitude of the rest-activity rhythm | (M10 - L5) / (M10 + L5); unitless | 0 to 1; higher = stronger rhythm | vansomeren1999 | circadian.rhythm()$`RA  |
|IS     |Interdaily stability: strength of coupling of the rhythm to the 24 h day |n*sum((xbar_h - xbar)^2) / (p*sum((x_i - xbar)^2)); unitless     |0 to 1; near 1 = highly stable, < 0.3 = weak     |witting1990    |circadian.rhythm()`$IS |
| IV | Intradaily variability: fragmentation of the rest-activity rhythm | n*sum((x_i - x\_{i-1})^2) / ((n-1)*sum((x_i - xbar)^2)); unitless | 0 to ~2; ~0 = smooth sine, ~2 = noisy/fragmented | witting1990 | circadian.rhythm()$`IV  |
|phi    |Autocorrelation-based predictability of the rhythm at the 24 h lag       |unitless                                                         |higher = more predictable day-to-day             |               |circadian.rhythm()`$phi |

## Cosinor and rhythmicity

The single-component cosinor and its zero-amplitude rhythmicity test
([Cornelissen, 2014](#ref-cornelissen2014); [Nelson et al.,
1979](#ref-nelson1979)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| cosinor_mesor | MESOR: the rhythm-adjusted mean of the fitted 24 h cosine | activity counts | non-negative | cornelissen2014 | cosinor.analysis()$`mesor          |
|cosinor_amplitude |Amplitude: half the peak-to-trough extent of the fitted cosine   |activity counts    |non-negative; larger = stronger sinusoidal component |cornelissen2014 |cosinor.analysis()`$amplitude |
| cosinor_acrophase | Acrophase: clock time of the cosine peak (circular) | clock hours (0-24) | 0 to 24 | cornelissen2014 | cosinor.analysis()$`acrophase      |
|cosinor_r_squared |Proportion of variance explained by the single cosine            |unitless           |0 to 1                                               |cornelissen2014 |cosinor.analysis()`$r_squared |
| rhythm_F | F statistic of the zero-amplitude (no-rhythm) test | unitless | non-negative; larger = stronger evidence of a rhythm | nelson1979 | rhythmicity.test()$`F              |
|rhythm_df2        |Denominator degrees of freedom of the rhythmicity F-test         |count              |non-negative                                         |nelson1979      |rhythmicity.test()`$df2 |
| rhythm_p_value | p-value of the zero-amplitude F-test | unitless | 0 to 1; \< 0.05 = detectable rhythm | nelson1979 | rhythmicity.test()$`p_value        |
|percent_rhythm    |Percent of variance the single cosine explains (100 * R-squared) |percent            |0 to 100                                             |nelson1979      |rhythmicity.test()`$percent_rhythm |
| rhythmic | Whether a rhythm is detected at the 0.05 level | logical | TRUE / FALSE | nelson1979 | rhythmicity.test()\$rhythmic |

## Period and spectral

The endogenous period from the Lomb-Scargle and chi-square periodograms,
with a bootstrap confidence interval ([Lomb, 1976](#ref-lomb1976);
[Sokolove & Bushell, 1978](#ref-sokolove1978)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| period_tau | Dominant free-running period from the Lomb-Scargle periodogram | hours | within the search band (default 18-30 h) | lomb1976 | circadian.period()$`tau        |
|period_peak_power |Lomb-Scargle spectral power at the dominant period                 |unitless        |non-negative; higher = stronger periodicity |lomb1976     |circadian.period()`$peak_power |
| period_p_value | False-alarm probability of the periodogram peak (Baluev) | unitless | 0 to 1; \< 0.05 = significant period | baluev2008 | circadian.period()$`p_value    |
|period_ci_lower   |Lower bound of the bootstrap period confidence interval            |hours           |<= period_tau                               |politis1992  |period.ci()`$ci_lower |
| period_ci_upper | Upper bound of the bootstrap period confidence interval | hours | \>= period_tau | politis1992 | period.ci()$`ci_upper          |
|period_ci_se      |Bootstrap standard error of the period estimate                    |hours           |non-negative                                |politis1992  |period.ci()`$se |
| chisq_period | Dominant period from the chi-square (Sokolove-Bushell) periodogram | hours | within the search band | sokolove1978 | chi.sq.periodogram()$`period   |
|chisq_Qp_peak     |Peak Qp statistic of the chi-square periodogram                    |unitless        |non-negative; above threshold = significant |sokolove1978 |chi.sq.periodogram()`$Qp_peak |
| chisq_p_value | Significance of the chi-square periodogram peak (Sidak-corrected) | unitless | 0 to 1 | sokolove1978 | chi.sq.periodogram()\$p_value |

## Fractal and nonlinear

Long-range correlation and complexity of the activity series ([Costa et
al., 2002](#ref-costa2002); [Kantelhardt et al.,
2002](#ref-kantelhardt2002); [Peng et al., 1994](#ref-peng1994)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| dfa_alpha | Overall detrended-fluctuation scaling exponent | unitless | ~0.5 = noise, ~1 = 1/f, \> 1 = random-walk-like | peng1994 | fractal.dfa()$`alpha        |
|dfa_alpha1  |Short-time-scale DFA exponent (below the crossover)               |unitless        |interpreted as for dfa_alpha                           |peng1994        |fractal.dfa()`$alpha1 |
| dfa_alpha2 | Long-time-scale DFA exponent (above the crossover) | unitless | interpreted as for dfa_alpha | peng1994 | fractal.dfa()$`alpha2       |
|mfdfa_h2    |Generalized Hurst exponent at q = 2 (equals the DFA exponent)     |unitless        |~0.5 to 1.5                                            |kantelhardt2002 |mfdfa()`$alpha_dfa |
| mfdfa_width | Width of the multifractal singularity spectrum | unitless | non-negative; wider = more multifractal | kantelhardt2002 | mfdfa()$`width              |
|mse_area    |Area under the multiscale sample-entropy curve (complexity index) |unitless        |non-negative; higher = more complex                    |costa2002       |multiscale.entropy()`$area |
| mse_slope | Slope of sample entropy across coarse-graining scales | unitless | sign indicates rising or falling complexity with scale | costa2002 | multiscale.entropy()\$slope |

## Rest-activity transitions

State-transition rates and probabilities that capture fragmentation a
single amplitude cannot ([Lim et al., 2011](#ref-lim2011)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| kRA | Rest-to-active state transition rate | per epoch | non-negative | lim2011 | state.transitions()$`kRA |
|kAR    |Active-to-rest state transition rate        |per epoch       |non-negative           |lim2011   |state.transitions()`$kAR |
| pRA | Rest-to-active state transition probability | unitless | 0 to 1 | lim2011 | state.transitions()$`pRA |
|pAR    |Active-to-rest state transition probability |unitless        |0 to 1                 |lim2011   |state.transitions()`$pAR |

## Sleep timing and regularity

Metrics that need a sleep label or sleep periods: the Sleep Regularity
Index, social jet lag, and LIDS ([Phillips et al.,
2017](#ref-phillips2017); [Winnebeck et al., 2018](#ref-winnebeck2018);
[Wittmann et al., 2006](#ref-wittmann2006)).

| Metric | Definition | Formula / units | Range / interpretation | Reference | Output object |
|:---|:---|:---|:---|:---|:---|
| SRI | Sleep Regularity Index: probability two epochs 24 h apart share the same state (requires sleep_state) | percent | approx 0 to 100; higher = more regular (100 = perfectly regular) | phillips2017 | sleep.regularity.index() |
| social_jet_lag_hours | Social jet lag: mid-sleep difference between free and work days (requires sleep_periods) | hours | typically 0 to 3; larger = greater misalignment | wittmann2006 | social.jet.lag()$`social_jet_lag_hours |
|MSW                  |Mid-sleep on work days                                                                                |clock hours (0-24) |0 to 24                                                          |wittmann2006  |social.jet.lag()`$MSW |
| MSF | Mid-sleep on free days | clock hours (0-24) | 0 to 24 | wittmann2006 | social.jet.lag()$`MSF                  |
|lids_period_min      |Mean Locomotor Inactivity During Sleep ultradian period (requires sleep_periods)                      |minutes            |typically 60 to 120 min                                          |winnebeck2018 |lids()`$mean_period_min |
| lids_MRI | Mean Munich Rhythmicity Index of the LIDS cycles (requires sleep_periods) | unitless | 0 to 1; higher = stronger ultradian rhythm | winnebeck2018 | lids()\$mean_MRI |

## References

Full citations for the method keys in the table above: ([Baluev,
2008](#ref-baluev2008); [Cornelissen, 2014](#ref-cornelissen2014);
[Costa et al., 2002](#ref-costa2002); [Kantelhardt et al.,
2002](#ref-kantelhardt2002); [Lim et al., 2011](#ref-lim2011); [Lomb,
1976](#ref-lomb1976); [Nelson et al., 1979](#ref-nelson1979); [Peng et
al., 1994](#ref-peng1994); [Phillips et al., 2017](#ref-phillips2017);
[Politis & Romano, 1992](#ref-politis1992); [Sokolove & Bushell,
1978](#ref-sokolove1978); [Van Someren et al.,
1999](#ref-vansomeren1999); [Winnebeck et al.,
2018](#ref-winnebeck2018); [Witting et al., 1990](#ref-witting1990);
[Wittmann et al., 2006](#ref-wittmann2006)).

Baluev, R. V. (2008). Assessing the statistical significance of
periodogram peaks. *Monthly Notices of the Royal Astronomical Society*,
*385*(3), 1279–1285. <https://doi.org/10.1111/j.1365-2966.2008.12689.x>

Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical Biology
and Medical Modelling*, *11*, 16.
<https://doi.org/10.1186/1742-4682-11-16>

Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Multiscale entropy
analysis of complex physiologic time series. *Physical Review Letters*,
*89*(6), 068102. <https://doi.org/10.1103/PhysRevLett.89.068102>

Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation
analysis of nonstationary time series. *Physica A: Statistical Mechanics
and Its Applications*, *316*(1-4), 87–114.
<https://doi.org/10.1016/S0378-4371(02)01383-3>

Lim, A. S. P., Yu, L., Costa, M. D., Leurgans, S. E., Buchman, A. S.,
Bennett, D. A., & Saper, C. B. (2011). Quantification of the
fragmentation of rest-activity patterns in elderly individuals using a
state transition analysis. *Sleep*, *34*(11), 1569–1581.
<https://doi.org/10.5665/sleep.1400>

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
