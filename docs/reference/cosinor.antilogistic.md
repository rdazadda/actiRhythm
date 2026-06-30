# Anti-Logistic (Extended) Cosinor Analysis (Marler et al. 2006)

Fits the sigmoidally transformed cosine ("anti-logistic" extended
cosinor) model of Marler et al. (2006) to the averaged 24-hour activity
profile. The extended model relaxes the symmetric shape of the ordinary
cosinor by adding two shape parameters: `alpha` sets the relative width
of the active versus rest phase, and `beta` sets the steepness of the
rest-to-active transitions. It is the same model as
[`ActCR::ActExtendCosinor`](https://rdrr.io/pkg/ActCR/man/ActExtendCosinor.html)
and reproduces its parameter estimates numerically.

## Usage

``` r
cosinor.antilogistic(counts, timestamps, period = 24)
```

## Arguments

- counts:

  Numeric vector of activity counts (one value per epoch).

- timestamps:

  POSIXct vector of timestamps, the same length as `counts`.

- period:

  Numeric period of the rhythm in hours (default `24`).

## Value

A list with class `"actiRhythm_cosinor_ext"` containing:

- minimum:

  Lower asymptote of the fitted curve (rest-phase level).

- amplitude:

  Vertical span of the sigmoidal transition (`amp`); the asymptotic
  ceiling is `minimum + amplitude`, with the fitted-curve maximum
  slightly below it at a finite steepness.

- alpha:

  Width-asymmetry parameter in \\\[-1, 1\]\\.

- beta:

  Steepness parameter (\\\ge 0\\).

- acrophase:

  Time of peak activity in clock hours (`acrotime`).

- acrotime:

  Alias of `acrophase` (ActCR naming).

- UpMesor:

  Clock time of the rest-to-active (rising) transition,
  \\-\arccos(\alpha) / (2\pi/T) + acrotime\\.

- DownMesor:

  Clock time of the active-to-rest (falling) transition,
  \\\arccos(\alpha) / (2\pi/T) + acrotime\\.

- MESOR:

  Midline statistic `minimum + amplitude / 2`.

- F_pseudo:

  Pseudo-F improvement of the extended fit over the ordinary cosinor:
  \\((RSS\_{cos} - RSS\_{ext})/2)/(RSS\_{ext}/(n-5))\\.

- rss_cosinor:

  Residual sum of squares of the ordinary cosinor.

- rss_extended:

  Residual sum of squares of the extended cosinor.

- converged:

  Logical; whether the nonlinear fit converged.

- period:

  Period used (hours).

- n_profile_hours:

  Number of profile hours used in the fit.

On non-convergence or insufficient data all numeric parameters are
returned as `NA` with `converged = FALSE`.

## Details

The model fitted to the averaged profile is \$\$f(t) = minimum +
amplitude \cdot expit(\beta (cos(2\pi (t - acrotime)/T) - \alpha))\$\$
where \\expit(z) = 1/(1+e^{-z})\\. This is the ActCR/Marler
parameterization, in which `amplitude` is the raw multiplier of the
logistic transform (it is *not* renormalized by
\\expit(\beta(1-\alpha))\\). The closely related "normalized" Marler
form, \\f(t) = min + amp \cdot expit(\beta(cos(\cdot) - \alpha)) /
expit(\beta(1-\alpha))\\, rescales `amplitude` so that the peak equals
exactly `minimum + amplitude`; the two forms describe the same curve
shape and the (`alpha`, `beta`, `acrotime`) parameters are identical.
The normalized peak level is reported as `peak`.

Starting values are taken from an ordinary cosinor fit
(`minimum = max(mesor - amp, 0)`, `amplitude = 2 * amp`, `alpha = 0`,
`beta = 2`, `acrotime = ordinary acrophase`). The nonlinear
least-squares problem is solved with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) using the
box-constrained L-BFGS-B method (bounds `lower = c(0, 0, -1, 0, -3)`,
`upper = c(Inf, Inf, 1, Inf, 27)`, matching ActCR), so no extra package
dependency is required.

## References

Marler MR, Gehrman P, Martin JL, Ancoli-Israel S (2006). “The
sigmoidally transformed cosine curve: a mathematical model for circadian
rhythms with symmetric non-sinusoidal shapes.” *Statistics in Medicine*,
**25**(22), 3893–3904.
[doi:10.1002/sim.2466](https://doi.org/10.1002/sim.2466) .

Wang J, Xian H, Di J, Zipunnikov V (2021). *ActCR: Extract Circadian
Rhythms Metrics from Actigraphy Data*. R package.

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01 00:00:00"), by = 60, length.out = 1440 * 3)
hour <- as.numeric(format(ts, "%H")) + as.numeric(format(ts, "%M")) / 60
counts <- 50 + 300 * plogis(4 * (cos((hour - 14) * 2 * pi / 24) - 0.2))
fit <- cosinor.antilogistic(counts, ts)
print(fit)
#> $minimum
#> [1] 49.62602
#> 
#> $amplitude
#> [1] 301.3645
#> 
#> $alpha
#> [1] 0.2014903
#> 
#> $beta
#> [1] 3.903463
#> 
#> $acrophase
#> [1] 14.00833
#> 
#> $acrotime
#> [1] 14.00833
#> 
#> $peak
#> [1] 338.2087
#> 
#> $UpMesor
#> [1] 8.783275
#> 
#> $DownMesor
#> [1] 19.23339
#> 
#> $MESOR
#> [1] 200.3083
#> 
#> $F_pseudo
#> [1] 632333.8
#> 
#> $rss_cosinor
#> [1] 13776.46
#> 
#> $rss_extended
#> [1] 0.2069704
#> 
#> $converged
#> [1] TRUE
#> 
#> $period
#> [1] 24
#> 
#> $n_profile_hours
#> [1] 24
#> 
#> attr(,"class")
#> [1] "actiRhythm_cosinor_ext"
# }
```
