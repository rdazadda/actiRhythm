# Coerce sleep/wake state to binary (1 = sleep, 0 = wake, NA = unscored)

Internal helper for
[`sri.matrix()`](https://rdazadda.github.io/actiRhythm/reference/sri.matrix.md).
Accepts the common encodings used across the package: character
`"S"`/`"W"` (case-insensitive), numeric/integer 1/0, or logical
TRUE/FALSE. Anything not recognised as sleep or wake (and any `NA`)
becomes `NA_real_`.

## Usage

``` r
.sri.to.binary(x)
```

## Arguments

- x:

  Vector of sleep states (character, numeric, integer, or logical).

## Value

Numeric vector the same length as `x` with values 1, 0, or `NA_real_`.
