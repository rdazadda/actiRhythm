# Calculate LX/MX using Standard Average-Profile Method (van Someren 1999)

Follows the van Someren (1999) method:

1.  Builds an average 24-hour activity profile across all days

2.  Slides a circular window over that profile to find L5/M10

## Usage

``` r
.calculate.LX.sliding(
  counts,
  timestamps,
  X,
  find_minimum = TRUE,
  epoch_length = 60
)
```

## Arguments

- counts:

  Numeric vector of activity counts

- timestamps:

  POSIXct timestamps

- X:

  Window size in hours (5 for L5, 10 for M10)

- find_minimum:

  TRUE for LX (least active), FALSE for MX (most active)

- epoch_length:

  Epoch length in seconds

## Value

List with value, start_time, start_hour (decimal)

## Details

Matches GGIR, nparACT, and the original paper.

## References

Van Someren EJW, Swaab DF, Colenda CC, Cohen W, McCall WV, Rosenquist PB
(1999). “Bright light therapy: improved sensitivity to its effects on
rest-activity rhythms in Alzheimer patients by application of
nonparametric methods.” *Chronobiology International*, **16**(4),
505–518.
[doi:10.3109/07420529908998724](https://doi.org/10.3109/07420529908998724)
.
