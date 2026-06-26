# Raw-Acceleration Epoch Metrics (ENMO, MAD, z-Angle)

Reads a raw accelerometer file and returns per-epoch raw activity and
posture metrics, the gravity-preserving signals that counts cannot
represent: ENMO (Euclidean Norm Minus One, a raw activity metric), MAD
(Mean Amplitude Deviation), and the z-angle (arm/posture angle).
Auto-calibration (van Hees 2014) is applied first by default. Requires
the relevant raw reader (read.gt3x for `.gt3x`, GGIRread for
`.cwa`/`.bin`).

## Usage

``` r
raw.metrics(
  x,
  device = "auto",
  epoch = 60,
  metrics = c("ENMO", "MAD", "anglez"),
  calibrate = TRUE,
  tz = "UTC"
)
```

## Arguments

- x:

  A path to a raw file (`.gt3x`, `.cwa`, `.bin`) or a raw data frame
  with `x`/`y`/`z` columns in g and an `fs` attribute (e.g. from
  [`example_raw`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)
  or your own device).

- device:

  One of `"auto"`, `"gt3x"`, `"axivity"`, `"geneactiv"` (default
  `"auto"`, inferred from the extension; used only when `x` is a file
  path).

- epoch:

  Epoch length in seconds (default 60).

- metrics:

  Which metrics to return; any of `"ENMO"`, `"MAD"`, `"anglez"` (default
  all three).

- calibrate:

  Apply van Hees auto-calibration first (default `TRUE`).

- tz:

  Time zone for the timestamps (default `"UTC"`).

## Value

A data frame with `time` and the requested metrics (ENMO and MAD in mg,
anglez in degrees), one row per epoch. The calibration result is
attached as the `"calibration"` attribute.

## References

van Hees VT, Gorzelniak L, Dean Leon EC, Eder M, Pias M, Taherian S,
Ekelund U, Renstrom F, Franks PW, Horsch A, Brage S (2013). “Separating
movement and gravity components in an acceleration signal and
implications for the assessment of human daily physical activity.” *PLoS
ONE*, **8**(4), e61691.
[doi:10.1371/journal.pone.0061691](https://doi.org/10.1371/journal.pone.0061691)
.

Vaha-Ypya H, Vasankari T, Husu P, Suni J, Sievanen H (2015). “A
universal, accurate intensity-based classification of different physical
activities using raw data of accelerometer.” *Clinical Physiology and
Functional Imaging*, **35**(1), 64–70.
[doi:10.1111/cpf.12127](https://doi.org/10.1111/cpf.12127) .

## See also

[`auto.calibrate`](https://rdazadda.github.io/actiRhythm/reference/auto.calibrate.md),
[`circadian.raw`](https://rdazadda.github.io/actiRhythm/reference/circadian.raw.md),
[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md),
[`example_raw`](https://rdazadda.github.io/actiRhythm/reference/example_raw.md)

## Examples

``` r
# On a synthetic raw recording (no file needed); pass a path for a real file
# \donttest{
m <- raw.metrics(example_raw(days = 1), epoch = 60)
head(m)
#>                  time     ENMO      MAD    anglez
#> 1 2024-01-01 12:00:00 43.95883 24.78791  8.280445
#> 2 2024-01-01 12:01:00 43.13860 24.89003 14.666095
#> 3 2024-01-01 12:02:00 44.56626 25.55921  8.287697
#> 4 2024-01-01 12:03:00 44.25995 24.90778 14.705513
#> 5 2024-01-01 12:04:00 44.79032 25.36485  8.292428
#> 6 2024-01-01 12:05:00 44.17071 24.72920 14.755796
# }
```
