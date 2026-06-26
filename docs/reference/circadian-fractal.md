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

Peng CK, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL
(1994). “Mosaic organization of DNA nucleotides.” *Physical Review E*,
**49**(2), 1685–1689.
[doi:10.1103/PhysRevE.49.1685](https://doi.org/10.1103/PhysRevE.49.1685)
.

Hu K, Ivanov PC, Chen Z, Carpena P, Stanley HE (2001). “Effect of trends
on detrended fluctuation analysis.” *Physical Review E*, **64**(1),
011114.
[doi:10.1103/PhysRevE.64.011114](https://doi.org/10.1103/PhysRevE.64.011114)
.

Hu K, Van Someren EJW, Shea SA, Scheer FAJL (2009). “Reduction of scale
invariance of activity fluctuations with aging and Alzheimer's disease.”
*Proceedings of the National Academy of Sciences*, **106**(8),
2490–2494.
[doi:10.1073/pnas.0806087106](https://doi.org/10.1073/pnas.0806087106) .

Richman JS, Moorman JR (2000). “Physiological time-series analysis using
approximate entropy and sample entropy.” *American Journal of
Physiology - Heart and Circulatory Physiology*, **278**(6), H2039–H2049.
[doi:10.1152/ajpheart.2000.278.6.H2039](https://doi.org/10.1152/ajpheart.2000.278.6.H2039)
.

Costa M, Goldberger AL, Peng C (2002). “Multiscale entropy analysis of
complex physiologic time series.” *Physical Review Letters*, **89**(6),
068102.
[doi:10.1103/PhysRevLett.89.068102](https://doi.org/10.1103/PhysRevLett.89.068102)
.

Costa M, Goldberger AL, Peng C (2005). “Multiscale entropy analysis of
biological signals.” *Physical Review E*, **71**(2), 021906.
[doi:10.1103/PhysRevE.71.021906](https://doi.org/10.1103/PhysRevE.71.021906)
.
