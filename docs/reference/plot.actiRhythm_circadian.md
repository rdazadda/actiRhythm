# Plot Circadian Rhythm Profile

Plots the circadian rhythm analysis.

## Usage

``` r
# S3 method for class 'actiRhythm_circadian'
plot(x, type = "profile", ...)
```

## Arguments

- x:

  actiRhythm_circadian object from circadian.rhythm()

- type:

  Type of plot: "profile" (default), "daily", or "all"

- ...:

  Additional arguments passed to plotting functions

## Value

ggplot object (or list of ggplot objects if type="all")
