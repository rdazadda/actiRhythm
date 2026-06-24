# Fractal and Complexity Metrics for Activity Time Series

Two nonlinear complexity and scaling metrics for accelerometer activity
time series: Detrended Fluctuation Analysis (DFA) and Multiscale Sample
Entropy (MSE). Both run in base R (only stats is used) so they add no
new runtime dependency. They measure the temporal structure (long-range
correlations and information-theoretic complexity) of minute-level
activity counts, alongside the amplitude- and timing-based circadian
metrics in
[`circadian.rhythm`](https://rdazadda.github.io/actiRhythm/reference/circadian.rhythm.md).

## References

**Detrended Fluctuation Analysis:**

- Peng CK, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL
  (1994). Mosaic organization of DNA nucleotides. Physical Review E,
  49(2):1685-1689.

- Hu K, Ivanov PC, Chen Z, Carpena P, Stanley HE (2001). Effect of
  trends on detrended fluctuation analysis. Physical Review E,
  64(1):011114.

- Hu K, Van Someren EJW, Shea SA, Scheer FAJL (2009). Reduction of scale
  invariance of activity fluctuations with aging and Alzheimer's
  disease. PNAS, 106(8):2490-2494.

**Multiscale Sample Entropy:**

- Richman JS, Moorman JR (2000). Physiological time-series analysis
  using approximate entropy and sample entropy. American Journal of
  Physiology - Heart and Circulatory Physiology, 278(6):H2039-H2049.

- Costa M, Goldberger AL, Peng CK (2002). Multiscale entropy analysis of
  complex physiologic time series. Physical Review Letters,
  89(6):068102.

- Costa M, Goldberger AL, Peng CK (2005). Multiscale entropy analysis of
  biological signals. Physical Review E, 71(2):021906.
