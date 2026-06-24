# Get Path to Example AGD Files

Get Path to Example AGD Files

## Usage

``` r
example_agd(file = 1)
```

## Arguments

- file:

  Character. Name of the example file or "list" to see available files.

## Value

Character. Full path to the example AGD file.

## Examples

``` r
example_agd()
#> [1] "D:/Rtemp/RtmpmIrHgl/temp_libpath838477b86d0e/actiRhythm/extdata/MOS2E39230594_60sec.agd"
example_agd("list")
#> [1] "MOS2E39230594_60sec.agd" "MOS2E3923063660sec.agd" 
agd_path <- example_agd(1)
```
