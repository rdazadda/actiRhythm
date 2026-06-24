# Joint Amplitude-Acrophase Confidence Ellipse (Bingham et al. 1982)

Computes the joint confidence region for the cosinor amplitude and
acrophase, following Bingham et al. (1982). The region is an ellipse in
the \\(\beta_1, \beta_2)\\ (cosine/sine coefficient) plane. When the
ellipse excludes the origin, the population amplitude differs
significantly from zero and a rhythm is detected.

## Usage

``` r
cosinor.confidence.ellipse(cosinor_result, level = 0.95, n_points = 200)
```

## Arguments

- cosinor_result:

  A list returned by
  [`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
  (class `"actiRhythm_cosinor"`). It must contain `amplitude`,
  `acrophase` (clock hours), `se_amplitude`, `period` and
  `n_profile_hours`. If the object additionally carries an explicit
  `vcov_beta` (2x2 covariance of \\(\beta_1, \beta_2)\\) and/or
  `beta1`/`beta2`, those are used directly.

- level:

  Confidence level for the region (default `0.95`).

- n_points:

  Number of points used to trace the ellipse boundary (default `200`).

## Value

A list with class `"actiRhythm_cosinor_ellipse"` containing:

- ellipse:

  Data frame with columns `x` (\\\beta_1\\) and `y` (\\\beta_2\\) giving
  the ellipse boundary vertices.

- center:

  Numeric vector `c(beta1, beta2)` of the point estimate.

- excludes_origin:

  Logical; `TRUE` when the ellipse excludes `(0, 0)` (i.e. a rhythm is
  detected at the requested level).

- rhythm_detected:

  Alias of `excludes_origin`.

- distance_stat:

  Mahalanobis-type statistic of the origin relative to the fitted
  ellipse.

- critical_value:

  Threshold \\2 F\_{2, df, level}\\ the statistic is compared against.

- level:

  Confidence level used.

When the input is missing required fields or is degenerate, the boundary
is returned as `NA` and `excludes_origin = NA`.

## Details

Under the cosinor model \\Y(t) = M + \beta_1 cos(\omega t) + \beta_2
sin(\omega t)\\, the \\100(1-\alpha)\\\\ joint confidence region for
\\(\beta_1, \beta_2)\\ is the set of points \\b\\ satisfying \$\$(b -
\hat b)^\top \Sigma^{-1} (b - \hat b) \le 2 F\_{2, df, 1-\alpha}\$\$
where \\\Sigma\\ is the covariance of the estimated coefficients and
\\df\\ the residual degrees of freedom. The averaged-profile cosinor
design used by
[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
has orthogonal cosine and sine columns, so \\\Sigma\\ is (very nearly)
diagonal with equal variances \\\sigma^2\\; in that case `se_amplitude`
(delta method) equals \\\sigma\\, which is how \\\Sigma\\ is
reconstructed when it is not supplied explicitly. The point estimate is
recovered from amplitude and acrophase as \\\beta_1 = A cos(\phi)\\,
\\\beta_2 = A sin(\phi)\\ with \\\phi = acrophase \cdot 2\pi / T\\.

## References

Bingham C, Arbogast B, Cornelissen Guillaume G, Lee JK, Halberg F
(1982). Inferential statistical methods for estimating and comparing
cosinor parameters. *Chronobiologia*, 9(4):397-439.

## Examples

``` r
# \donttest{
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 1440 * 3)
counts <- 100 + 80 * cos(2 * pi * (as.numeric(format(ts, "%H")) - 14) / 24)
cos <- cosinor.analysis(counts, ts)
cosinor.confidence.ellipse(cos)$rhythm_detected
#> [1] NA
# }
```
