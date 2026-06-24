# Convert a Raw .gt3x File to an .agd File

Computes counts with
[`gt3x.counts`](https://rdazadda.github.io/actiRhythm/reference/gt3x.counts.md)
and writes them to a small `.agd`, carrying device/subject metadata from
the `.gt3x` header. Requires agcounts.

## Usage

``` r
gt3x.to.agd(gt3x_path, agd_path = NULL, epoch = 60, lfe = FALSE, tz = "UTC")
```

## Arguments

- gt3x_path:

  Path to the `.gt3x` file.

- agd_path:

  Output `.agd` path; defaults next to the input.

- epoch:

  Epoch length in seconds (default 60).

- lfe:

  Use the low-frequency extension filter (default `FALSE`).

- tz:

  Time zone (default `"UTC"`).

## Value

The output `.agd` path.
