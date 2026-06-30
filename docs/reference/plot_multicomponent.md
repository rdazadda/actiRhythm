# Single Cosine vs Multicomponent Cosinor Fit

Overlays the averaged daily profile with the single 24-hour cosine and
the selected multi-harmonic cosinor (Cornelissen 2014), showing the
structure a single symmetric cosine cannot follow. The number of
harmonics and the fit are taken from
[`cosinor.multicomponent`](https://rdazadda.github.io/actiRhythm/reference/cosinor.multicomponent.md);
the curves are refit on the hour-of-day profile so they align with the
plotted points. Returns a `ggplot` object and never errors.

## Usage

``` r
plot_multicomponent(counts, timestamps, period = 24, max_harmonics = 3)
```

## Arguments

- counts:

  Numeric activity vector.

- timestamps:

  POSIXct timestamps, one per value.

- period:

  Fundamental period in hours (default 24).

- max_harmonics:

  Largest number of harmonics to consider (default 3).

## Value

A `ggplot` object.

## References

Cornelissen G (2014). “Cosinor-based rhythmometry.” *Theoretical Biology
and Medical Modelling*, **11**, 16.
[doi:10.1186/1742-4682-11-16](https://doi.org/10.1186/1742-4682-11-16) .

## Examples

``` r
set.seed(2)
ts <- seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 60, length.out = 4 * 1440)
h  <- as.numeric(format(ts, "%H"))
siesta <- pmax(0, 120 + 70 * cos(2 * pi * (h - 14) / 24) +
  55 * cos(2 * pi * 2 * (h - 14) / 24) + rnorm(length(ts), 0, 12))
plot_multicomponent(siesta, ts)

```
