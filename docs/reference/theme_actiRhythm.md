# actiRhythm ggplot2 Theme

A clean, professional theme for actiRhythm visualizations.
Publication-ready with excellent readability and consistent with the
dashboard's 6-size typography system.

## Usage

``` r
theme_actiRhythm(base_size = 14, base_family = "", grid = TRUE, dark = FALSE)
```

## Arguments

- base_size:

  Numeric. Base font size in points (default: 14)

- base_family:

  Character. Base font family (default: "")

- grid:

  Logical. Show major grid lines? (default: TRUE)

- dark:

  Logical. Use dark mode? (default: FALSE)

## Value

A ggplot2 theme object

## Details

Typography sizes (matching CSS):

- Caption/footnotes: 11px (0.79x base)

- Labels/metadata: 12px (0.86x base)

- Body text: 14px (base_size)

- Emphasis/titles: 16px (1.14x base)

- Section headings: 20px (1.43x base)

- Page titles: 24px (1.71x base)

## Examples

``` r
# \donttest{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  theme_actiRhythm()

# }
```
