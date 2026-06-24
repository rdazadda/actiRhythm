# Circadian Rhythm Visualizations

ggplot2 visualizations for the circadian / chronobiology metrics
produced by actiRhythm. These functions are thin plotting wrappers
around the existing analytic engines and add no new fitting logic:

- [`plot_periodogram`](https://rdazadda.github.io/actiRhythm/reference/plot_periodogram.md)
  draws the Lomb-Scargle periodogram and highlights the endogenous
  period estimated by
  [`circadian.period`](https://rdazadda.github.io/actiRhythm/reference/circadian.period.md).

- [`plot_extended_cosinor`](https://rdazadda.github.io/actiRhythm/reference/plot_extended_cosinor.md)
  overlays the Marler extended (anti-logistic) cosinor fit from
  [`cosinor.antilogistic`](https://rdazadda.github.io/actiRhythm/reference/cosinor.antilogistic.md)
  and the ordinary cosinor from
  [`cosinor.analysis`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  on the averaged 24-hour activity profile.

- [`plot_dfa`](https://rdazadda.github.io/actiRhythm/reference/plot_dfa.md)
  draws the detrended-fluctuation log-log scaling plot from
  [`fractal.dfa`](https://rdazadda.github.io/actiRhythm/reference/fractal.dfa.md).

Every function returns a `ggplot` object, never errors on degenerate or
insufficient input (it returns an annotated empty plot instead), and
falls back to
[`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
when the package theme helpers are unavailable.

## References

Lomb NR (1976). Least-squares frequency analysis of unequally spaced
data. *Astrophysics and Space Science*, 39(2):447-462.

Scargle JD (1982). Studies in astronomical time series analysis. II.
Statistical aspects of spectral analysis of unevenly spaced data. *The
Astrophysical Journal*, 263:835-853.

Marler MR, Gehrman P, Martin JL, Ancoli-Israel S (2006). The sigmoidally
transformed cosine curve: a mathematical model for circadian rhythms
with symmetric non-sinusoidal shapes. *Statistics in Medicine*,
25(22):3893-3904.

Peng CK, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL
(1994). Mosaic organization of DNA nucleotides. *Physical Review E*,
49(2):1685-1689.
