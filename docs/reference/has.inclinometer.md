# Check if AGD Data Has Inclinometer Information

Checks whether an AGD file contains inclinometer data, which can be used
to enhance sleep detection algorithms.

## Usage

``` r
has.inclinometer(agd_data)
```

## Arguments

- agd_data:

  List returned from read.agd()

## Value

Logical. TRUE if inclinometer data is available
