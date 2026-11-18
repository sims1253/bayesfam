# Check if a vector is numeric, has no na entries and has length len

Check if a vector is numeric, has no na entries and has length len

## Usage

``` r
isNum_len(x, len = 1)
```

## Arguments

- x:

  Numeric vector to be checked

- len:

  Length of vector, default argument is 1

## Value

Boolean, whether x was numeric and of correct size

## Examples

``` r
bayesfam:::isNum_len(c(1.1, 2.2), 2) # should be TRUE
#> [1] TRUE
bayesfam:::isNum_len(0.2) # should be TRUE
#> [1] TRUE
```
