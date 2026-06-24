# Calculate Interdaily Stability (IS) and Intradaily Variability (IV)

IS measures consistency of activity patterns across days. IV measures
fragmentation of the activity rhythm within days.

## Usage

``` r
.calculate.IS.IV(counts, timestamps, epoch_length = 60)
```

## Arguments

- counts:

  Numeric vector of activity counts

- timestamps:

  POSIXct timestamps

- epoch_length:

  Epoch length in seconds

## Value

List with IS and IV values

## Details

**IS Formula (Witting et al., 1990):** IS = (n \* sum((Xh - Xmean)^2)) /
(p \* sum((Xi - Xmean)^2))

Where:

- n = total number of hourly data points

- p = number of hours per day (24)

- Xh = mean activity for hour h across all days

- Xmean = overall mean activity

**IV Formula (Witting et al., 1990):** IV = (n \* sum((Xi - Xi-1)^2)) /
((n-1) \* sum((Xi - Xmean)^2))

## References

Witting W, Kwa IH, Eikelenboom P, Mirmiran M, Swaab DF (1990).
Alterations in the circadian rest-activity rhythm in aging and
Alzheimer's disease. Biological Psychiatry, 27(6):563-572.
[doi:10.1016/0006-3223(90)90523-5](https://doi.org/10.1016/0006-3223%2890%2990523-5)
