# Calculate Social Jet Lag

Difference between sleep midpoint on work days and free days.

## Usage

``` r
social.jet.lag(sleep_periods, work_days = NULL)
```

## Arguments

- sleep_periods:

  Data frame with in_bed_time, out_bed_time columns

- work_days:

  Optional. Logical vector, date vector, or NULL (uses Mon-Fri default)

## Value

List with social jet lag metrics, including `MSW`/`MSF` (mid-sleep on
work and free days), `social_jet_lag_hours` (MSF - MSW), the
sleep-debt-corrected chronotype `MSFsc` and corrected
`social_jet_lag_sc_hours` (SJLsc), and the number of work/free nights.

## Details

Social jet lag is the discrepancy between social and biological time,
computed as MSF - MSW (midpoint on free days minus work days).

`MSFsc` applies the MCTQ sleep-debt correction (Roenneberg et al. 2012):
when free-day sleep exceeds work-day sleep (sleep catch-up), the
free-day mid-sleep is adjusted down by half the excess over the weekly
average sleep duration, and `SJLsc = MSFsc - MSW`. These are based on
time in bed (the supplied in-bed/out-bed times), not on sleep onset and
offset.

Positive values (most common) indicate later sleep timing on free days.
Values \> 1 hour are associated with increased health risks.

## References

Wittmann M, Dinich J, Merrow M, Roenneberg T (2006). “Social jetlag:
misalignment of biological and social time.” *Chronobiology
International*, **23**(1-2), 497–509.
[doi:10.1080/07420520500545979](https://doi.org/10.1080/07420520500545979)
.

Roenneberg T, Allebrandt KV, Merrow M, Vetter C (2012). “Social jetlag
and obesity.” *Current Biology*, **22**(10), 939–943.
[doi:10.1016/j.cub.2012.03.038](https://doi.org/10.1016/j.cub.2012.03.038)
.
