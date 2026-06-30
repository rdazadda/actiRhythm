# z-Angle Sleep-Period-Time Detection

Plots the z-angle of a raw recording with the HDCZA sleep-period-time
window(s) shaded and non-wear marked, the companion figure to
[`rest.spt`](https://rdazadda.github.io/actiRhythm/reference/rest.spt.md).
Computes the angle, non-wear, and SPT window from the raw input. Returns
a `ggplot` object and never errors.

## Usage

``` r
plot_spt(x, epoch_length = 5)
```

## Arguments

- x:

  A path to a raw file or a raw data frame (as in
  [`raw.metrics`](https://rdazadda.github.io/actiRhythm/reference/raw.metrics.md)).

- epoch_length:

  Epoch length in seconds (default 5, matching `rest.spt`).

## Value

A `ggplot` object.

## References

van Hees VT, Sabia S, Jones SE, Wood AR, Anderson KN, Kivimaki M,
Frayling TM, Pack AI, Bucan M, Trenell MI, Mazzotti DR, Gehrman PR,
Singh-Manoux BA, Weedon MN (2018). “Estimating sleep parameters using an
accelerometer without sleep diary.” *Scientific Reports*, **8**, 12975.
[doi:10.1038/s41598-018-31266-z](https://doi.org/10.1038/s41598-018-31266-z)
.

## Examples

``` r
# \donttest{
plot_spt(example_raw(days = 2))

# }
```
