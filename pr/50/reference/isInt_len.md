# Integer vector check

Integer vector check

## Usage

``` r
isInt_len(int, len = 1)
```

## Arguments

- int:

  Integer vector to be checked

- len:

  Length of vector, default argument is 1

## Value

Boolean, whether int was Integer and of correct size

## Examples

``` r
bayesfam:::isInt_len(c(1, 2), 2) # should be TRUE
#> [1] TRUE
bayesfam:::isInt_len(1, 2) # should be FALSE, wrong length
#> [1] FALSE
```
