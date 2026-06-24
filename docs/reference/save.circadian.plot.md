# Save a actiRhythm Plot to a File

Convenience wrapper around
[`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)
for the `plot_*` circadian figures; the output format is taken from the
file extension (`.png`, `.pdf`, `.svg`).

## Usage

``` r
save.circadian.plot(plot, file, width = 8, height = 5, dpi = 300)
```

## Arguments

- plot:

  A ggplot object (e.g. from
  [`plot_actogram`](https://rdazadda.github.io/actiRhythm/reference/plot_actogram.md)).

- file:

  Output path; its extension sets the format.

- width, height:

  Size in inches (defaults 8 x 5).

- dpi:

  Resolution for raster formats (default 300).

## Value

The `file` path, invisibly.
