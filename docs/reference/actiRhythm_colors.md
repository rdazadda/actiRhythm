# actiRhythm Color Generator

Generates colors from the actiRhythm palettes for visualizations. These
colors match the CSS design system.

## Usage

``` r
actiRhythm_colors(n = NULL, type = "categorical")
```

## Arguments

- n:

  Integer. Number of colors to return. If NULL, returns all.

- type:

  Character. One of "categorical", "sequential", "diverging", or
  "intensity"

## Value

A character vector of hex color codes

## Examples

``` r
actiRhythm_colors(4)
#> [1] "#0072B2" "#E69F00" "#009E73" "#CC79A7"
# \donttest{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point() +
  scale_color_manual(values = actiRhythm_colors(3))

# }
```
