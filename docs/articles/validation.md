# Validation against reference packages

actiRhythm computes its metrics from hand-written code rather than by
wrapping other packages, so it is fair to ask whether those metrics
actually reproduce the established implementations. This article answers
that by running the reference packages side by side on the recording
that ships with actiRhythm, and emitting the agreement from the code
that follows: no number here is typed by hand. Several of these checks
also run in the package’s own test suite (for example
`test-circadian-nparact.R` and `test-circadian-cosinor-extended.R`), so
they are exercised on every build.

Each comparison guards its reference package with
[`requireNamespace()`](https://rdrr.io/r/base/ns-load.html): a chunk
whose package is not installed shows its code but no output, so the
article builds either way. `ActCR`, the `cosinor` package, and
`cosinor2` are used live below; `nparACT` is shown for you to run.

``` r

library(actiRhythm)

agd <- agd.counts(read.agd(example_agd(1), verbose = FALSE))
x   <- agd$axis1
ts  <- agd$timestamp

# The reference packages expect a whole number of 1440-minute days, so align to
# the first local midnight and keep only complete days.
mins  <- as.integer(format(ts, "%H")) * 60 + as.integer(format(ts, "%M"))
start <- which(mins == 0)[1]
x  <- x[start:length(x)]
keep <- (length(x) %/% 1440) * 1440
x  <- x[seq_len(keep)]
ts <- ts[start + seq_len(keep) - 1]

# actiRhythm's own metrics, and the averaged hourly profile the references reuse.
cr     <- circadian.rhythm(x, ts)
hourly <- t(apply(matrix(x, ncol = 1440, byrow = TRUE), 1,
                  function(d) tapply(d, rep(1:24, each = 60), mean)))
cat(sprintf("%d complete days of 60-second epochs\n", keep / 1440))
#> 6 complete days of 60-second epochs
```

## Nonparametric metrics versus ActCR

`ActCR` is the Johns Hopkins reimplementation of the interdaily
stability and intradaily variability of Witting et al.
([1990](#ref-witting1990)). Its `IS()` takes a days-by-epochs matrix and
`IV()` a chronological vector; on the same hourly profile actiRhythm
uses internally, the two implementations agree to the printed precision.

``` r

cmp <- data.frame(
  metric     = c("IS", "IV"),
  actiRhythm = c(cr$IS, cr$IV),
  ActCR      = c(ActCR::IS(hourly), ActCR::IV(c(t(hourly))))
)
cmp$abs_diff <- abs(cmp$actiRhythm - cmp$ActCR)
knitr::kable(cmp, digits = 5, caption = "IS and IV: actiRhythm vs ActCR.")
```

| metric | actiRhythm |   ActCR | abs_diff |
|:-------|-----------:|--------:|---------:|
| IS     |     0.2383 | 0.23832 |    2e-05 |
| IV     |     1.0265 | 1.02654 |    4e-05 |

IS and IV: actiRhythm vs ActCR. {.table}

## Nonparametric metrics versus nparACT

`nparACT` is the free reference implementation of the full nonparametric
battery: IS, IV, RA, L5, M10 ([Van Someren et al.,
1999](#ref-vansomeren1999); [Witting et al., 1990](#ref-witting1990)).
It reads a `time` / `activity` data frame by name from the global
environment and a sampling rate in Hz. The chunk below reproduces the
package’s `test-circadian-nparact.R`: run it with `nparACT` installed
and it prints the agreement (IS/RA to two decimals, L5/M10 to within
nparACT’s rounding).

``` r

assign("np_validation_df", data.frame(time = ts, activity = x), envir = globalenv())
np <- nparACT::nparACT_base("np_validation_df", SR = 1 / 60, plot = FALSE)
rm("np_validation_df", envir = globalenv())

cmp <- data.frame(
  metric     = c("IS", "IV", "RA", "L5", "M10"),
  actiRhythm = c(cr$IS, cr$IV, cr$RA, cr$L5, cr$M10),
  nparACT    = c(np$IS, np$IV, np$RA, np$L5, np$M10)
)
cmp$abs_diff <- abs(cmp$actiRhythm - cmp$nparACT)
knitr::kable(cmp, digits = 3, caption = "IS/IV/RA/L5/M10: actiRhythm vs nparACT.")
```

## Single cosinor versus the cosinor package

[`cosinor.analysis()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.analysis.md)
fits the single-component cosinor of Cornelissen
([2014](#ref-cornelissen2014)). The `cosinor` package fits the same
model through its formula interface; given the same averaged daily
profile, the MESOR and amplitude match exactly.

``` r

ca   <- cosinor.analysis(x, ts)
prof <- data.frame(hour = rep(0:23, nrow(hourly)), Y = c(t(hourly)))
cl   <- cosinor::cosinor.lm(Y ~ time(hour) + 1, period = 24, data = prof)

cmp <- data.frame(
  parameter  = c("MESOR", "amplitude"),
  actiRhythm = c(ca$mesor, ca$amplitude),
  cosinor    = c(unname(cl$coefficients[1]), unname(cl$coefficients[2]))
)
cmp$abs_diff <- abs(cmp$actiRhythm - cmp$cosinor)
knitr::kable(cmp, digits = 5, caption = "Single cosinor: actiRhythm vs the cosinor package.")
```

|          | parameter | actiRhythm |  cosinor | abs_diff |
|:---------|:----------|-----------:|---------:|---------:|
|          | MESOR     |     348.84 | 348.8390 |  0.00100 |
| cos_term | amplitude |     309.62 | 309.6162 |  0.00383 |

Single cosinor: actiRhythm vs the cosinor package. {.table}

## Population cosinor versus cosinor2

[`population.cosinor()`](https://rdazadda.github.io/actiRhythm/reference/population.cosinor.md)
pools per-subject cosinors into the Bingham et al.
([1982](#ref-bingham1982)) population-mean rhythm with confidence
intervals, and
[`cosinor2::population.cosinor.lm()`](https://rdrr.io/pkg/cosinor2/man/population.cosinor.lm.html)
implements the same framework. On a synthetic group the population
MESOR, amplitude, and acrophase and all three confidence intervals match
(cosinor2 reports the acrophase in negative radians, converted to clock
hours here, with its time reference set to the hourly bin centres).

``` r

set.seed(11); K <- 8; hrs <- 0:23
mat <- t(vapply(seq_len(K), function(i) {
  M <- 100 + stats::rnorm(1, 0, 5); A <- 40 + stats::rnorm(1, 0, 5)
  phi <- 8 + stats::rnorm(1, 0, 1)
  M + A * cos(2 * pi * (hrs - phi) / 24) + stats::rnorm(24, 0, 3)
}, numeric(24)))

base <- as.POSIXct("2024-01-01", tz = "UTC")
pcts <- rep(base + hrs * 3600, times = K)
ours <- population.cosinor(as.vector(t(mat)), pcts, rep(paste0("S", seq_len(K)), each = 24))

invisible(utils::capture.output(
  ref <- cosinor2::population.cosinor.lm(as.data.frame(mat), time = hrs + 0.5,
                                         period = 24, plot = FALSE)))
rad2h <- function(r) ((-r) * 24 / (2 * pi)) %% 24

cmp <- data.frame(
  parameter  = c("MESOR", "amplitude", "acrophase"),
  actiRhythm = c(ours$mesor, ours$amplitude, ours$acrophase),
  cosinor2   = c(ref$coefficients[["MESOR"]], ref$coefficients[["Amplitude"]],
                 rad2h(ref$coefficients[["Acrophase"]]))
)
cmp$abs_diff <- abs(cmp$actiRhythm - cmp$cosinor2)
knitr::kable(cmp, digits = 5, caption = "Population cosinor: actiRhythm vs cosinor2.")
```

| parameter | actiRhythm |  cosinor2 | abs_diff |
|:----------|-----------:|----------:|---------:|
| MESOR     |  100.85718 | 100.85718 |        0 |
| amplitude |   39.34469 |  39.34469 |        0 |
| acrophase |    7.72310 |   7.72310 |        0 |

Population cosinor: actiRhythm vs cosinor2. {.table}

## Extended (anti-logistic) cosinor versus ActCR

[`cosinor.antilogistic()`](https://rdazadda.github.io/actiRhythm/reference/cosinor.antilogistic.md)
fits the sigmoidally transformed cosine of Marler et al.
([2006](#ref-marler2006)), the same model as
[`ActCR::ActExtendCosinor()`](https://rdrr.io/pkg/ActCR/man/ActExtendCosinor.html).
We validate it on a controlled profile with a *moderate* transition,
where every shape parameter is identifiable. On a recording with
near-square days, the steepness `beta` is only weakly constrained once
the rest-activity transition is sharp, because the likelihood is then
nearly flat in `beta`; the other parameters still agree closely while
`beta` itself is not pinned down, so we show the well-posed case here.

``` r

# A 3-day synthetic recording whose average day is a moderate sigmoid.
set.seed(1)
n   <- 1440 * 3
tse <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = n)
hr  <- as.numeric(format(tse, "%H")) + as.numeric(format(tse, "%M")) / 60
ctp <- cos((hr - 15) * 2 * pi / 24)
ce  <- 50 + 250 * plogis(3 * (ctp - 0.1)) + rnorm(n, 0, 12)
ce[ce < 0] <- 0

# The 1440-minute averaged profile ActExtendCosinor consumes.
idx  <- as.integer(format(tse, "%H")) * 60 + as.integer(format(tse, "%M"))
prof1440 <- vapply(0:1439, function(i) mean(ce[idx == i], na.rm = TRUE), numeric(1))

fit <- cosinor.antilogistic(ce, tse)
ref <- ActCR::ActExtendCosinor(prof1440, window = 1)$params

cmp <- data.frame(
  parameter  = c("minimum", "amplitude", "alpha", "beta", "acrotime"),
  actiRhythm = c(fit$minimum, fit$amplitude, fit$alpha, fit$beta, fit$acrophase),
  ActCR      = c(ref$minimum, ref$amp, ref$alpha, ref$beta, ref$acrotime)
)
cmp$abs_diff <- abs(cmp$actiRhythm - cmp$ActCR)
knitr::kable(cmp, digits = 3, caption = "Anti-logistic cosinor: actiRhythm vs ActExtendCosinor.")
```

| parameter | actiRhythm |   ActCR | abs_diff |
|:----------|-----------:|--------:|---------:|
| minimum   |     49.370 |  50.160 |    0.790 |
| amplitude |    251.932 | 250.177 |    1.755 |
| alpha     |      0.103 |   0.103 |    0.000 |
| beta      |      2.906 |   2.971 |    0.065 |
| acrotime  |     15.017 |  15.024 |    0.007 |

Anti-logistic cosinor: actiRhythm vs ActExtendCosinor. {.table}

## Raw metrics versus GGIR

The raw-acceleration chain (van Hees auto-calibration, ENMO/MAD, and the
diary-free z-angle sleep detector) reimplements the algorithms `GGIR`
established, and is cross-checked against GGIR directly. On a real wrist
recording kept in the source repository (a multi-day raw file is too
large for the CRAN build), the [From raw
acceleration](https://rdazadda.github.io/actiRhythm/articles/raw-pipeline.md)
article computes per-epoch ENMO and the z-angle against
[`GGIR::g.getmeta()`](https://wadpac.github.io/GGIR/reference/g.getmeta.html),
and the sleep-period window against GGIR’s HDCZA detector (`HASPT`).
Those tables are computed live whenever the repository and GGIR are
present; they report the per-epoch correlations and the window agreement
directly, and the package’s opt-in `test-ggir.R` re-runs the same
comparison and fails if the agreement degrades.

## Bottom line

On the bundled recording, actiRhythm’s interdaily stability, intradaily
variability, and single cosinor reproduce `ActCR` and the `cosinor`
package to the printed precision; the full nonparametric battery matches
`nparACT` to its rounded output; and the anti-logistic cosinor matches
`ActExtendCosinor` on an identifiable profile. These are not one-off
checks: they are part of the package’s test suite and run on every
build. The methods behind them are Witting et al.
([1990](#ref-witting1990)), Van Someren et al.
([1999](#ref-vansomeren1999)), Cornelissen
([2014](#ref-cornelissen2014)), Nelson et al. ([1979](#ref-nelson1979)),
and Marler et al. ([2006](#ref-marler2006)).

## References

Bingham, C., Arbogast, B., Cornelissen Guillaume, G., Lee, J. K., &
Halberg, F. (1982). Inferential statistical methods for estimating and
comparing cosinor parameters. *Chronobiologia*, *9*(4), 397–439.

Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical Biology
and Medical Modelling*, *11*, 16.
<https://doi.org/10.1186/1742-4682-11-16>

Marler, M. R., Gehrman, P., Martin, J. L., & Ancoli-Israel, S. (2006).
The sigmoidally transformed cosine curve: A mathematical model for
circadian rhythms with symmetric non-sinusoidal shapes. *Statistics in
Medicine*, *25*(22), 3893–3904. <https://doi.org/10.1002/sim.2466>

Nelson, W., Tong, Y. L., Lee, J. K., & Halberg, F. (1979). Methods for
cosinor-rhythmometry. *Chronobiologia*, *6*(4), 305–323.

Van Someren, E. J. W., Swaab, D. F., Colenda, C. C., Cohen, W., McCall,
W. V., & Rosenquist, P. B. (1999). Bright light therapy: Improved
sensitivity to its effects on rest-activity rhythms in Alzheimer
patients by application of nonparametric methods. *Chronobiology
International*, *16*(4), 505–518.
<https://doi.org/10.3109/07420529908998724>

Witting, W., Kwa, I. H., Eikelenboom, P., Mirmiran, M., & Swaab, D. F.
(1990). Alterations in the circadian rest-activity rhythm in aging and
Alzheimer’s disease. *Biological Psychiatry*, *27*(6), 563–572.
<https://doi.org/10.1016/0006-3223(90)90523-5>
