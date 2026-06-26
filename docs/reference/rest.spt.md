# Sleep-Period-Time Window from the z-Angle (HDCZA)

Detects the main sleep-period-time (SPT) window per day directly from
the z-angle, with no sleep diary, using the van Hees et al. (2018)
Heuristic algorithm based on the Distribution of Change in Z-Angle
(HDCZA): the absolute 5-second change in arm angle is smoothed over 5
minutes, thresholded at a fraction of its own distribution, and the
longest sustained low-change block per noon-to-noon day becomes the SPT.
This is raw-native: there is no count-based equivalent.

## Usage

``` r
rest.spt(
  anglez,
  timestamps,
  epoch_length = 5,
  pct = 10,
  mult = 15,
  clamp = c(0.13, 0.5),
  min_block = 30,
  max_gap = 60,
  algo = c("HDCZA", "HorAngle"),
  wear = NULL
)
```

## Arguments

- anglez:

  Numeric z-angle (degrees) per epoch, e.g. from
  `raw.metrics(..., metrics = "anglez")` at a short epoch.

- timestamps:

  POSIXct timestamps, one per epoch.

- epoch_length:

  Epoch length in seconds (default 5, the validated value).

- pct:

  Percentile of the change distribution for the threshold (default 10).

- mult:

  Multiplier on that percentile (default 15).

- clamp:

  Lower/upper clamp (degrees) on the threshold (default
  `c(0.13, 0.50)`).

- min_block:

  Minimum SPT block length in minutes (default 30).

- max_gap:

  Maximum movement gap to bridge within the SPT, minutes (default 60).

- algo:

  `"HDCZA"` (wrist, change-in-angle) or `"HorAngle"` (a hip variant
  thresholding the absolute angle at 60 degrees).

- wear:

  Optional logical wear mask, one per epoch (e.g. from
  [`detect.nonwear.raw`](https://rdazadda.github.io/actiRhythm/reference/detect.nonwear.raw.md));
  non-wear epochs are excluded from the SPT so a stationary, taken-off
  device is not scored as sleep.

## Value

An object of class `actiRhythm_spt`: a data frame with one row per day
(`date`, `onset`, `offset`, `duration` in hours, and the `threshold`
used). Never errors; returns no rows if no day has a detectable window.

## References

van Hees VT, Sabia S, Jones SE, Wood AR, Anderson KN, Kivimaki M,
Frayling TM, Pack AI, Bucan M, Trenell MI, Mazzotti DR, Gehrman PR,
Singh-Manoux BA, Weedon MN (2018). “Estimating sleep parameters using an
accelerometer without sleep diary.” *Scientific Reports*, **8**, 12975.
[doi:10.1038/s41598-018-31266-z](https://doi.org/10.1038/s41598-018-31266-z)
.

## See also

[`sib.vanhees`](https://rdazadda.github.io/actiRhythm/reference/sib.vanhees.md),
[`sleep.from.spt`](https://rdazadda.github.io/actiRhythm/reference/sleep.from.spt.md),
[`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)

## Examples

``` r
# One still night inside an active day yields one SPT window
ts <- seq(as.POSIXct("2024-01-01 12:00", tz = "UTC"), by = 5, length.out = 17280)
h <- as.numeric(format(ts, "%H")); night <- h >= 23 | h < 7
set.seed(1)
anglez <- ifelse(night, -60, -30) + rnorm(17280, 0, ifelse(night, 0.02, 20))
rest.spt(anglez, ts, epoch_length = 5)
#>         date               onset              offset duration threshold
#> 1 2024-01-01 2024-01-01 23:00:05 2024-01-02 06:59:55 7.997222 0.2670371
```
