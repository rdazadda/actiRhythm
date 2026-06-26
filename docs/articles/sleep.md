# Sleep and rest detection

``` r

library(actiRhythm)
library(ggplot2)
```

## The idea, and a rule to remember

Sleep and rest detection turns a stream of activity counts into an
answer to two different questions. The **epoch scorers** – Cole-Kripke
and Sadeh – label every single minute as sleep or wake, giving a
per-epoch state you can total into sleep time. The **bout detectors** –
change-point, Roenneberg, Crespo, and the state-space model – step back
and mark the *spans* of rest, the onset and offset clock times of each
consolidated period. One tells you *how much* sleep, the other *when* it
happened.

The rule to carry through this article: **decide first whether you want
one main sleep window per night or every rest bout the recording
contains, because that choice – not the algorithm’s accuracy – is what
separates these functions.**
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
finds the single dominant rest bout of each circadian cycle;
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
and
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
return *every* consolidated bout, naps included. Asking a one-bout
method for naps, or a many-bout method for “the night”, is the most
common way to misread the output.

## The math

The two epoch scorers are weighted windows over the counts.
**Cole-Kripke** ([Cole et al., 1992](#ref-cole1992)) forms a sleep index
from a seven-epoch window – four epochs before the current one
($`P_4 \dots P_1`$), the current epoch $`C`$, and two after
($`N_1, N_2`$) – with counts scaled by 100 and capped:
``` math
D = 0.001\,(106\,P_4 + 54\,P_3 + 58\,P_2 + 76\,P_1 + 230\,C + 74\,N_1 + 67\,N_2),
```
and scores the epoch **sleep when $`D < 1`$**. Webster’s rescoring rules
([Webster et al., 1982](#ref-webster1982)) then re-label short sleep
bouts bracketed by sustained wake.

**Sadeh** ([Sadeh et al., 1994](#ref-sadeh1994)) uses an eleven-epoch
window (five each side) with counts capped at 300:
``` math
SI = 7.601 - 0.065\,\overline{A} - 1.08\,\mathrm{NAT} - 0.056\,s - 0.703\,\log(C+1),
```
where $`\overline{A}`$ is the window mean, $`\mathrm{NAT}`$ the count of
epochs in $`[50,100)`$, $`s`$ the standard deviation over the current
and five preceding epochs, and $`C`$ the current count. The epoch is
scored **sleep when $`SI > -4`$**.

The change-point detector ([Chen & Sun, 2024](#ref-chensun2024)) is
geometric rather than weighted: it fits a 24-hour cosinor to bound each
rest and active span roughly, then places the exact transition inside
each bound at the split $`k`$ minimising the within-segment residual sum
of squares, \$\$ k^\\ = \arg\min_k \\\sum\_{i \le k}(x_i - \bar x\_{\le
k})^2 + \sum\_{i \> k}(x_i - \bar x\_{\>k})^2 . \$\$

## Assumptions, and when they break

- **One-minute epochs.** Cole-Kripke and Sadeh were validated on
  60-second epochs and their coefficients are tied to that resolution;
  both warn if `epoch_length` differs. Re-bin to one minute before
  scoring.
- **Counts, not raw acceleration.** The scorers expect ActiGraph-style
  activity counts on the vertical axis; feed them counts (`agd$axis1`),
  not raw g.
- **Non-wear is not sleep.** A taken-off device reads as deep rest. The
  scorers’ `na_action` argument governs this – the default `"na"` emits
  `NA` for missing-count epochs so a device-off gap is not silently
  scored as sleep.
- **A clear day-night contrast for the bout detectors.**
  [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md),
  [`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md),
  and
  [`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
  all need enough rhythm to separate rest from activity; on a flat or
  very short record they return an empty result rather than erroring.

## Recovering known truth

Before trusting any of this on real data, plant an answer and check it
comes back. We build a seven-day recording with a **known nightly sleep
window from 23:00 to 07:00** – near-zero activity at night, a busy day –
and ask whether the epoch scorers label the window sleep and whether
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
recovers the 23:00 onset and 07:00 wake.

``` r

ts  <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 7 * 1440)
hod <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
set.seed(1)
asleep <- hod >= 23 | hod < 7                       # the planted sleep window
counts <- ifelse(asleep, 0, 250) + pmax(0, rnorm(length(ts), 0, 8))

ck <- sleep.cole.kripke(counts)
sd <- sleep.sadeh(counts)
knitr::kable(
  data.frame(
    scorer        = c("Cole-Kripke", "Sadeh"),
    sleep_recall  = c(mean(ck[asleep] == "S"),  mean(sd[asleep] == "S")),
    wake_recall   = c(mean(ck[!asleep] == "W"), mean(sd[!asleep] == "W"))),
  digits = 3,
  caption = "Both scorers label the planted 23:00-07:00 window sleep and the day wake."
)
```

| scorer      | sleep_recall | wake_recall |
|:------------|-------------:|------------:|
| Cole-Kripke |        0.992 |       0.999 |
| Sadeh       |        0.993 |       1.000 |

Both scorers label the planted 23:00-07:00 window sleep and the day
wake. {.table}

Both scorers tag essentially all of the night as sleep and essentially
all of the day as wake. Now the change-point detector, which should
place onset near 23:00 and wake near 07:00 every night:

``` r

cp <- sleep.changepoints(counts, ts)
cps_h <- as.numeric(format(cp$changepoints$time, "%H")) +
         as.numeric(format(cp$changepoints$time, "%M")) / 60
knitr::kable(
  data.frame(
    transition = c("sleep onset", "wake onset"),
    planted    = c(23, 7),
    recovered  = c(mean(cps_h[cp$changepoints$type == "sleep onset"]),
                   mean(cps_h[cp$changepoints$type == "wake onset"]))),
  digits = 2,
  caption = "Mean recovered transition time across the seven nights vs the planted values."
)
```

| transition  | planted | recovered |
|:------------|--------:|----------:|
| sleep onset |      23 |     22.98 |
| wake onset  |       7 |      6.98 |

Mean recovered transition time across the seven nights vs the planted
values. {.table}

``` r

c(n_episodes = cp$n_episodes, mean_sleep_duration_h = cp$mean_sleep_duration)
#>            n_episodes mean_sleep_duration_h 
#>                     6                     8
```

The recovered onset and wake land within a minute of the planted 23:00
and 07:00, the detector finds exactly seven nightly episodes, and the
mean duration is the eight hours we built in. The method recovers the
truth we planted.

## On a real recording

The bundled recording runs the same calls. The two scorers return a
per-epoch state vector;
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
returns a per-night episode table.

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
table(cole_kripke = sleep.cole.kripke(agd$axis1))
#> cole_kripke
#>    S    W 
#> 7646 2273
table(sadeh       = sleep.sadeh(agd$axis1))
#> sadeh
#>    S    W 
#> 7190 2729
cp_real <- sleep.changepoints(agd$axis1, agd$timestamp)
cp_real
#> Change-Point Sleep/Wake Detection (CircaCP)
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

The two scorers agree closely on the total sleep fraction, and the
change-point detector resolves one main sleep episode per night with
onset and wake clock times.

## Reading the numbers

- **`sleep_state` (“S”/“W”)** from the scorers is a per-epoch label. Its
  useful summaries are the *count* of sleep epochs (multiply by the
  epoch length for sleep time) and the *runs* of consecutive sleep
  (sleep periods and their lengths).
- **Onset and wake times** from
  [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
  are clock times; read them as “when did the main rest bout begin and
  end”, and read `mean_sleep_duration` as the average nightly rest
  length.
- **Cole-Kripke vs Sadeh.** Cole-Kripke was validated for adults ([Cole
  et al., 1992](#ref-cole1992)), Sadeh for children and adolescents
  ([Sadeh et al., 1994](#ref-sadeh1994)); pick by population. They will
  not agree epoch-for-epoch, and the disagreement is largest at the
  sleep-wake boundary, not in the middle of the night.

## One bout per night, or every bout? – the rest detectors compared

This is the rule in action, and it stands in for the usual “wrong-way”
demo: there is no wrong algorithm here, only a wrong *question*. The
same recording, run through four detectors, gives four different bout
counts – because they answer different questions, not because three of
them are wrong.

``` r

cp <- sleep.changepoints(agd$axis1, agd$timestamp)   # one main bout per night
rp <- rest.periods(agd$axis1, agd$timestamp)         # every bout (Roenneberg/MASDA)
rc <- rest.crespo(agd$axis1, agd$timestamp)          # every bout (Crespo morphology)
hm <- rest.hmm(agd$axis1, agd$timestamp, seed = 1)   # state-space alternative

knitr::kable(
  data.frame(
    detector = c("sleep.changepoints", "rest.periods", "rest.crespo", "rest.hmm"),
    question = c("one main bout / night", "every bout (naps incl.)",
                 "every bout (naps incl.)", "latent rest state"),
    bouts    = c(cp$n_episodes, rp$n_bouts, rc$n_rest_periods, NA),
    family   = c("CircaCP", "Roenneberg", "Crespo", "Gaussian HMM")),
  caption = "Four detectors, four counts: the difference is the question, not the accuracy."
)
```

| detector           | question                | bouts | family       |
|:-------------------|:------------------------|------:|:-------------|
| sleep.changepoints | one main bout / night   |     7 | CircaCP      |
| rest.periods       | every bout (naps incl.) |     9 | Roenneberg   |
| rest.crespo        | every bout (naps incl.) |     4 | Crespo       |
| rest.hmm           | latent rest state       |    NA | Gaussian HMM |

Four detectors, four counts: the difference is the question, not the
accuracy. {.table}

[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
([Chen & Sun, 2024](#ref-chensun2024)) returns one episode per circadian
cycle: the single dominant night-time rest bout.
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
([Loock et al., 2021](#ref-loock2021); [Roenneberg et al.,
2015](#ref-roenneberg2015)) and
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
([Crespo et al., 2012](#ref-crespo2012)) instead consolidate **every**
rest bout, so a daytime nap or a fragmented night becomes its own row –
they reach that result by two independent routes, Roenneberg’s
trend-and-correlation consolidation versus Crespo’s
median-filter-and-morphology pipeline, which is why their bout counts
differ even though both aim to catch all rest.
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
([Huang et al., 2018](#ref-huang2018hmm)) is the state-space
alternative: rather than thresholding each epoch, it fits a Gaussian
hidden Markov model whose latent rest state *persists*, and reports a
24-hour rest-occupancy profile.

``` r

ggplot(hm$tod_profile, aes(hour, p_rest)) +
  geom_col(fill = "#236192") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(x = "Hour of day", y = "P(rest state)") +
  theme_actiRhythm()
```

![The HMM's probability of being in the rest state across the day. The
state-space model infers the rest band -- centred on the night --
without a fixed count
threshold.](sleep_files/figure-html/hmm-profile-1.png)

The HMM’s probability of being in the rest state across the day. The
state-space model infers the rest band – centred on the night – without
a fixed count threshold.

The lesson is the rule: if you report
[`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
you are reporting *the night*; if you report
[`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
or
[`rest.crespo()`](https://rdazadda.github.io/actiRhythm/reference/rest.crespo.md)
you are reporting *all rest*, naps and all. Choose the one whose
question matches yours, rather than treating the differing bout counts
as a contradiction.

## The wider sleep-and-rest family

**The state-space view in detail.** Beyond the occupancy profile above,
[`rest.hmm()`](https://rdazadda.github.io/actiRhythm/reference/rest.hmm.md)
returns the per-state emission means, the transition matrix, and a
decoded `state_path` with its own `sleep_state` vector, plus AIC/BIC for
choosing between a two- and three-state model ([Huang et al.,
2018](#ref-huang2018hmm)).

``` r

hm$emission[, c("state", "label", "mean_transformed")]
#>   state  label mean_transformed
#> 1     1   rest    6.326781e-295
#> 2     2 active     2.812874e+01
c(time_at_rest = mean(hm$state_path == 1L),
  rest_persistence = hm$transition[1, 1], AIC = hm$AIC)
#>     time_at_rest rest_persistence              AIC 
#>     7.373727e-01     9.514564e-01    -1.604721e+05
```

**The ultradian structure within sleep.** Once you have sleep windows,
[`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md)
describes the sleep-cycle oscillation *inside* them. It needs a
`sleep_periods` data frame with `in_bed_time` and `out_bed_time` columns
– which we build directly from the change-point episodes – transforms
activity to locomotor inactivity ($`100/(\text{activity}+1)`$), and fits
an ultradian cosine, reporting the best period and its Munich
Rhythmicity Index ([Winnebeck et al., 2018](#ref-winnebeck2018)).

``` r

sleep_periods <- data.frame(
  in_bed_time  = cp$sleep_episodes$sleep_onset,
  out_bed_time = cp$sleep_episodes$wake_onset)
li <- lids(agd$axis1, agd$timestamp, sleep_periods)
li
#> Locomotor Inactivity During Sleep (LIDS)
#> 
#>   Sleep periods:    7
#>   Mean LIDS period: 167.9 min
#>   Mean MRI:         2.844
```

The mean LIDS period lands in the ninety-minute-to-three-hour band
expected of human sleep cycles, recovered from the inactivity signal
alone.

## Limitations

- **Epoch length is fixed for the scorers.** Cole-Kripke and Sadeh are
  only valid at one-minute epochs; do not port their coefficients to
  other resolutions.
- **Bout counts are not comparable across detectors.** A larger bout
  count from
  [`rest.periods()`](https://rdazadda.github.io/actiRhythm/reference/rest.periods.md)
  than from
  [`sleep.changepoints()`](https://rdazadda.github.io/actiRhythm/reference/sleep.changepoints.md)
  is by design, not a sign one is better; compare a detector only
  against itself or against its own reference implementation.
- **No in-bed truth.** None of these methods see a diary or a true
  lights-off; they infer rest from movement, so a still-but-awake
  subject reads as rest and an active sleeper as wake. Where a diary
  exists, gate on it.
- **[`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md)
  depends on the supplied sleep windows.** Feeding it poorly placed
  windows propagates straight into the period estimate; build the
  windows from a detector you trust on the recording at hand.

## Reference and validation

The epoch scorers follow Cole et al. (1992) ([Cole et al.,
1992](#ref-cole1992)) with the Webster rescoring rules ([Webster et al.,
1982](#ref-webster1982)) and Sadeh et al. (1994) ([Sadeh et al.,
1994](#ref-sadeh1994)). The bout detectors implement the CircaCP
change-point method of Chen and Sun (2024) ([Chen & Sun,
2024](#ref-chensun2024)), the Roenneberg / MASDA consolidation
([Roenneberg et al., 2015](#ref-roenneberg2015)) with the open
re-implementation of Loock et al. (2021) ([Loock et al.,
2021](#ref-loock2021)), the Crespo et al. (2012) morphology pipeline
([Crespo et al., 2012](#ref-crespo2012)), and the hidden-Markov
rest-activity model of Huang et al. (2018) ([Huang et al.,
2018](#ref-huang2018hmm));
[`lids()`](https://rdazadda.github.io/actiRhythm/reference/lids.md)
follows Winnebeck et al. (2018) ([Winnebeck et al.,
2018](#ref-winnebeck2018)). actiRhythm’s scorers and detectors are
cross-checked against their reference implementations – to the printed
precision – in the
[Validation](https://rdazadda.github.io/actiRhythm/articles/validation.md)
article and the package’s test suite.

## References

Chen, S., & Sun, X. (2024). Validating CircaCP: A generic sleep-wake
cycle detection algorithm for unlabelled actigraphy data. *Royal Society
Open Science*, *11*(5), 231468. <https://doi.org/10.1098/rsos.231468>

Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
(1992). Automatic sleep/wake identification from wrist activity.
*Sleep*, *15*(5), 461–469. <https://doi.org/10.1093/sleep/15.5.461>

Crespo, C., Aboy, M., Fernandez, J. R., & Mojon, A. (2012). Automatic
identification of activity-rest periods based on actigraphy. *Medical &
Biological Engineering & Computing*, *50*(4), 329–340.
<https://doi.org/10.1007/s11517-012-0875-y>

Huang, Q., Cohen, D., Komarzynski, S., Li, X.-M., Innominato, P., Levi,
F., & Finkenstadt, B. (2018). Hidden Markov models for monitoring
circadian rhythmicity in telemetric activity data. *Journal of the Royal
Society Interface*, *15*(139), 20170885.
<https://doi.org/10.1098/rsif.2017.0885>

Loock, A.-S., Khan Sullivan, A., Reis, C., Paiva, T., Ghotbi, N., Pilz,
L. K., Biller, A. M., Molenda, C., Vuori-Brodowski, M. T., Roenneberg,
T., & Winnebeck, E. C. (2021). Validation of the Munich actimetry sleep
detection algorithm for estimating sleep-wake patterns from activity
recordings. *Journal of Sleep Research*, *30*(6), e13371.
<https://doi.org/10.1111/jsr.13371>

Roenneberg, T., Keller, L. K., Fischer, D., Matera, J. L., Vetter, C., &
Winnebeck, E. C. (2015). Human activity and rest in situ. *Methods in
Enzymology*, *552*, 257–283.
<https://doi.org/10.1016/bs.mie.2014.11.028>

Sadeh, A., Sharkey, K. M., & Carskadon, M. A. (1994). Activity-based
sleep-wake identification: An empirical test of methodological issues.
*Sleep*, *17*(3), 201–207. <https://doi.org/10.1093/sleep/17.3.201>

Webster, J. B., Kripke, D. F., Messin, S., Mullaney, D. J., & Wyborney,
G. (1982). An activity based sleep monitor system for ambulatory use.
*Sleep*, *5*(4), 389–399. <https://doi.org/10.1093/sleep/5.4.389>

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and ultradian structure of human sleep in real life. *Current
Biology*, *28*(1), 49–59. <https://doi.org/10.1016/j.cub.2017.11.063>
