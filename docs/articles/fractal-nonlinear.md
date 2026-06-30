# Fractal and nonlinear dynamics

``` r

library(actiRhythm)
library(ggplot2)
```

## The idea, and a rule to remember

The fractal and nonlinear metrics ask a different question from the
cosinor or the nonparametric battery. Those describe the **shape and
depth** of the average day; the fractal tools describe the **temporal
texture** of the minute-to-minute fluctuations: how this minute’s
activity remembers the last, and how that memory stretches across scales
from minutes to hours. Healthy activity is neither smooth nor purely
random: it carries a self-similar, “1/f” structure, and that structure
breaks down in aging and disease in ways the amplitude metrics cannot
see ([Hu et al., 2009](#ref-hu2009)).

The rule to carry through this article: **the scaling exponent measures
correlation, not size.** A signal can be large and regular yet have the
same $`\alpha`$ as a tiny one; $`\alpha`$ reports how the fluctuations
are organised across time, independent of their amplitude. Read it as
the rhythm’s memory, not its strength.

## The math

**Detrended fluctuation analysis (DFA).** Start from activity
$`x_1,\dots,x_N`$ with mean $`\bar x`$, and integrate the mean-centred
signal into a profile ([Hu et al., 2001](#ref-hu2001); [Peng et al.,
1994](#ref-peng1994)):
``` math
Y(k) = \sum_{i=1}^{k} \bigl(x_i - \bar x\bigr).
```
Split $`Y`$ into $`M = \lfloor N/n \rfloor`$ non-overlapping windows of
length $`n`$, fit and remove a least-squares polynomial (linear by
default, quadratic for DFA-2) inside each window, and pool the residuals
into the fluctuation
``` math
F(n) = \sqrt{\tfrac{1}{M n}\sum_{k=1}^{M n}\bigl(Y(k) - \widehat{Y}_n(k)\bigr)^2}.
```
If the signal is self-similar, $`F(n)`$ grows as a power law,
$`F(n)\propto
n^{\alpha}`$, and the **scaling exponent** $`\alpha`$ is the slope of
$`\log F(n)`$ on $`\log n`$. The benchmarks are exact:
$`\alpha\approx 0.5`$ for uncorrelated white noise,
$`\alpha\approx 1.0`$ for $`1/f`$ (pink) noise, and
$`\alpha\approx 1.5`$ for a Brownian random walk.

**Multifractal DFA (MF-DFA).** A single $`\alpha`$ assumes one scaling
law governs the whole signal. MF-DFA relaxes that by weighting the
windowed variances $`F^2(v,n)`$ by a moment order $`q`$([Kantelhardt et
al., 2002](#ref-kantelhardt2002)):
``` math
F_q(n) = \Bigl[\tfrac{1}{2 N_s}\sum_{v}\bigl(F^2(v,n)\bigr)^{q/2}\Bigr]^{1/q}
        \propto n^{\,h(q)} ,
```
with windows taken from both ends (hence $`2 N_s`$), the $`q = 0`$ case
using the geometric mean, and zero-variance windows dropped. $`h(2)`$ is
the ordinary DFA exponent. If $`h(q)`$ is flat the signal is
**monofractal**; if $`h(q)`$ decreases with $`q`$ the signal is
**multifractal**, and the **spectrum width** $`\Delta\alpha`$ (the
spread of Holder exponents) measures how much. The width depends on the
$`q`$ range, so compare it only across recordings analysed with the same
grid (the default is $`q \in [-5, 5]`$).

**Multiscale entropy (MSE).** Complexity is not the same as scaling.
Sample entropy counts how unpredictable a series is: the negative log of
the chance that two sub-sequences matching for $`m`$ points still match
for $`m+1`$ points, within a tolerance $`r`$([Richman & Moorman,
2000](#ref-richman2000)):
``` math
\mathrm{SampEn}(m, r) = -\ln \frac{A}{B}.
```
MSE computes $`\mathrm{SampEn}`$ on coarse-grained copies of the series
across scales $`\tau`$([Costa et al., 2002](#ref-costa2002)). Noise
loses entropy as it is averaged; genuinely complex, correlated signals
hold their entropy across scales.

## Assumptions, and when they break

- **A long, gap-free segment.** Scaling laws need many octaves of window
  sizes; all three functions analyse the **longest continuous non-NA
  run** and return NA (never an error) on series too short to fit a
  slope. A few days of minute-level data is the practical minimum.
- **Stationary fluctuations, not trends.** DFA removes a linear trend
  inside each window, so a slow drift is tolerated, but a strong
  embedded rhythm or a device-off step can masquerade as long-range
  correlation. Gate on wear time first.
- **Amplitude-blind, by design.** $`\alpha`$ is invariant to scaling the
  signal, so it adds to the amplitude metrics rather than replacing
  them; report it alongside the cosinor and nonparametric numbers.
- **Heavy zeros distort the negative moments.** Activity counts have
  long runs of exact zeros, which give near-zero windowed variances.
  Raised to a negative power in MF-DFA these blow up, inflating the
  spectrum width; the $`q>0`$ side and $`h(2)`$ stay trustworthy.

## Recovering known truth

Before trusting $`\alpha`$ on a real recording, watch it recover
dynamics whose answer is fixed by construction. We build three series of
the same length and variance but different memory: **white noise** (no
correlation, planted $`\alpha=0.5`$), **$`1/f`$ pink noise** synthesised
so its power spectrum falls as $`1/f`$ (planted $`\alpha=1.0`$), and a
**Brownian random walk**, the running sum of white noise (planted
$`\alpha=1.5`$).

``` r

set.seed(42)
n <- 8192

white <- rnorm(n)                                   # planted alpha = 0.5
brown <- cumsum(rnorm(n))                            # planted alpha = 1.5
pink  <- {                                           # planted alpha = 1.0
  amp  <- 1 / sqrt(1:(n / 2))                        # 1/f amplitude spectrum
  ph   <- runif(n / 2, 0, 2 * pi)
  half <- complex(modulus = c(0, amp), argument = c(0, ph))
  Re(fft(c(half, Conj(rev(half[2:(n / 2)]))), inverse = TRUE))[1:n]
}

planted <- c(white = 0.5, pink = 1.0, brown = 1.5)
recovered <- c(white = fractal.dfa(white)$alpha,
               pink  = fractal.dfa(pink)$alpha,
               brown = fractal.dfa(brown)$alpha)
knitr::kable(
  data.frame(process = names(planted), planted = planted, recovered_alpha = recovered),
  row.names = FALSE, digits = 3,
  caption = "DFA recovers the planted scaling exponent: 0.5 for white noise, 1.0 for 1/f, 1.5 for a random walk."
)
```

| process | planted | recovered_alpha |
|:--------|--------:|----------------:|
| white   |     0.5 |           0.502 |
| pink    |     1.0 |           0.997 |
| brown   |     1.5 |           1.487 |

DFA recovers the planted scaling exponent: 0.5 for white noise, 1.0 for
1/f, 1.5 for a random walk. {.table}

The exponent tracks the planted process almost exactly. The fluctuation
plot for the random-walk (brown) series shows the scaling directly: a
single straight line on log-log axes, whose slope is the exponent.

``` r

plot_dfa(brown)
```

![DFA fluctuation curve for the random-walk (brown) series on log-log
axes. The slope is the scaling exponent (~1.5 here); one straight line
is monofractal
scaling.](fractal-nonlinear_files/figure-html/dfa-plot-1.png)

DFA fluctuation curve for the random-walk (brown) series on log-log
axes. The slope is the scaling exponent (~1.5 here); one straight line
is monofractal scaling.

## On a real recording

The bundled recording runs the same way.
[`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
returns a typed object that prints its own exponents and an interpretive
guide; `alpha1` and `alpha2` split the short and long timescales at the
`breakpoint_min` window (90 minutes by default).

``` r

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
dfa <- fractal.dfa(agd$axis1)
dfa
#> 
#> Detrended Fluctuation Analysis (Peng et al., 1994)
#> 
#>   Samples analyzed:   9919
#>   Window sizes:       29 (range 4-2479)
#>   alpha  (overall):   0.9264
#>   alpha1 (n < 90):     0.9581
#>   alpha2 (n >= 90):    0.8658
#> 
#>   Guide: near 0.5 = uncorrelated, near 1.0 = 1/f, near 1.5 = Brownian
```

The overall $`\alpha`$ sits in the **0.9-1.0 healthy band** (the $`1/f`$
signature of intact activity regulation) rather than near the 0.5 of
noise or the 1.5 of a random walk. That band, and the aging and disease
reductions below, were established by Hu et al. ([2009](#ref-hu2009))
with quadratic detrending (DFA-2);
[`fractal.dfa()`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md)
defaults to linear DFA-1, so set `detrend_order = 2` to compare directly
against those reference values.

## Reading the numbers

State each exponent in human terms:

- **$`\alpha`$** runs from about 0.5 (uncorrelated) through 1.0
  ($`1/f`$, the healthy target) to 1.5 (random-walk smoothness).
  Reductions toward 0.5 are reported with aging and Alzheimer’s disease
  ([Hu et al., 2009](#ref-hu2009)).
- **$`\alpha_1`$ vs $`\alpha_2`$** compare memory below and above the
  breakpoint; a drop from short to long timescales is a crossover,
  common in activity where the circadian band changes the scaling.
- **Spectrum width $`\Delta\alpha`$** (MF-DFA) is near zero for a
  monofractal and grows with multifractality; read it only on the
  $`q>0`$ side for raw counts.
- **Complexity / area** (MSE) is the entropy summed across scales:
  higher means the signal stays unpredictable even when smoothed, the
  hallmark of physiological complexity rather than noise.

A useful pairing is $`\alpha`$ with the MSE slope: $`1/f`$ scaling and
sustained entropy across scales together say “complex and healthy”,
which neither number asserts alone.

## The wider fractal and nonlinear family

**From one exponent to a spectrum.**
[`mfdfa()`](https://rdazadda.github.io/actiRhythm/reference/mfdfa.md)
asks whether a *single* scaling law fits, or whether different-sized
fluctuations scale differently. We contrast a **monofractal** signal
(white noise, one law) against a **multifractal** one (a binomial
multiplicative cascade, the textbook multifractal): the spectrum width
separates them cleanly while $`h(2)`$ stays interpretable as the DFA
exponent.

``` r

set.seed(7)
mono <- rnorm(8192)
cascade <- Reduce(function(v, .) as.numeric(rbind(v * 0.3, v * 0.7)),
                  1:13, accumulate = FALSE)        # binomial multiplicative cascade
multi <- cascade / sd(cascade)

knitr::kable(
  data.frame(
    signal     = c("monofractal (white noise)", "multifractal (cascade)"),
    h2         = c(mfdfa(mono)$alpha_dfa, mfdfa(multi)$alpha_dfa),
    width      = c(mfdfa(mono)$width,     mfdfa(multi)$width)
  ),
  col.names = c("signal", "h(2)", "spectrum width"),
  row.names = FALSE, digits = 3,
  caption = "A narrow spectrum width flags a monofractal; a wide one flags multifractality."
)
```

| signal                    |  h(2) | spectrum width |
|:--------------------------|------:|---------------:|
| monofractal (white noise) | 0.495 |          0.131 |
| multifractal (cascade)    | 0.806 |          1.228 |

A narrow spectrum width flags a monofractal; a wide one flags
multifractality. {.table}

On the real recording, MF-DFA confirms the DFA exponent through
$`h(2)`$; the raw width is inflated by the zero-variance windows noted
above, so we report $`h(2)`$ as the reliable summary.

``` r

mf <- mfdfa(agd$axis1)
c(h2 = unname(mf$alpha_dfa))
#>        h2 
#> 0.9175375
```

``` r

plot_mfdfa(agd$axis1)
```

![MF-DFA on the recording: the generalized Hurst exponent h(q) with h(2)
dashed, and the singularity spectrum. The wide negative-q side reflects
the zero-variance windows noted above, so read h(2) and the q \> 0
side.](fractal-nonlinear_files/figure-html/mfdfa-plot-1.png)

MF-DFA on the recording: the generalized Hurst exponent h(q) with h(2)
dashed, and the singularity spectrum. The wide negative-q side reflects
the zero-variance windows noted above, so read h(2) and the q \> 0 side.

**Complexity across scales.**
[`multiscale.entropy()`](https://rdazadda.github.io/actiRhythm/reference/multiscale.entropy.md)
applies sample entropy to coarse-grained copies of the recording,
separating signals that stay complex when smoothed from noise that loses
entropy with scale.

``` r

mse <- multiscale.entropy(agd$axis1, scales = 1:20)
mse
#> 
#> Multiscale Sample Entropy (Costa et al., 2002)
#> 
#>   Samples analyzed:   9919
#>   Scales:             20 (1-20)
#>   SampEn @ scale 1:   0.0727
#>   Complexity (area):  1.4636
#>   Slope (mse on scale):0.0005
#> 
#>   Negative slope => noise-like; flat/positive => complex
```

``` r

plot_mse(agd$axis1, scales = 1:20)
```

![Sample entropy across coarse-graining scales. A flat or rising curve
is the signature of physiological complexity; a curve that falls toward
zero is noise-like.](fractal-nonlinear_files/figure-html/mse-plot-1.png)

Sample entropy across coarse-graining scales. A flat or rising curve is
the signature of physiological complexity; a curve that falls toward
zero is noise-like.

The complexity index (`area`, the entropy summed across scales)
condenses that curve to a single number for comparison between
recordings.

## Limitations

- **Data hunger.** Scaling exponents are unstable on short or gappy
  series; the functions fall back to the longest clean run and return NA
  when even that is too short. Treat $`\alpha`$ from under a few days
  with caution.
- **Embedded rhythm leaks into the slope.** A strong 24-hour component
  adds curvature to the log-log plot near the daily scale; the
  $`\alpha_1/\alpha_2`$ split exposes it, but a single pooled $`\alpha`$
  can hide it.
- **Negative-moment fragility.** MF-DFA spectrum width on raw activity
  counts is inflated by zero-variance windows; lean on $`h(2)`$ and the
  $`q>0`$ side, or pre-process the zeros, before reporting
  $`\Delta\alpha`$.
- **Complementary, not diagnostic.** A reduced $`\alpha`$ is a marker,
  not a diagnosis; interpret it next to the amplitude, timing, and
  fragmentation metrics, never in isolation.

## Reference and validation

DFA follows Peng et al. ([1994](#ref-peng1994)), with trend and
crossover handling from Hu et al. ([2001](#ref-hu2001)) and the
activity-specific framing of Hu et al. ([2009](#ref-hu2009)); the
multifractal generalisation is Kantelhardt et al.
([2002](#ref-kantelhardt2002)); and multiscale entropy combines the
sample-entropy estimator of Richman & Moorman ([2000](#ref-richman2000))
with the multiscale construction of Costa et al.
([2002](#ref-costa2002)). actiRhythm’s DFA and MSE are implemented in
base R against these reference algorithms and exercised against their
known limits ($`\alpha\to0.5`$ for white noise, $`\to1.5`$ for a random
walk, entropy declining with scale for noise) in the package’s test
suite; see also the
[Validation](https://rdazadda.github.io/actiRhythm/articles/validation.md)
article.

## References

Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Multiscale entropy
analysis of complex physiologic time series. *Physical Review Letters*,
*89*(6), 068102. <https://doi.org/10.1103/PhysRevLett.89.068102>

Hu, K., Ivanov, P. C., Chen, Z., Carpena, P., & Stanley, H. E. (2001).
Effect of trends on detrended fluctuation analysis. *Physical Review E*,
*64*(1), 011114. <https://doi.org/10.1103/PhysRevE.64.011114>

Hu, K., Van Someren, E. J. W., Shea, S. A., & Scheer, F. A. J. L.
(2009). Reduction of scale invariance of activity fluctuations with
aging and Alzheimer’s disease. *Proceedings of the National Academy of
Sciences*, *106*(8), 2490–2494.
<https://doi.org/10.1073/pnas.0806087106>

Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation
analysis of nonstationary time series. *Physica A: Statistical Mechanics
and Its Applications*, *316*(1-4), 87–114.
<https://doi.org/10.1016/S0378-4371(02)01383-3>

Peng, C. K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E., &
Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
*Physical Review E*, *49*(2), 1685–1689.
<https://doi.org/10.1103/PhysRevE.49.1685>

Richman, J. S., & Moorman, J. R. (2000). Physiological time-series
analysis using approximate entropy and sample entropy. *American Journal
of Physiology - Heart and Circulatory Physiology*, *278*(6),
H2039–H2049. <https://doi.org/10.1152/ajpheart.2000.278.6.H2039>
