# Boolean vector check

Boolean vector check

## Usage

``` r
isLogic_len(logic, len = 1)
```

## Arguments

- logic:

  Logic vector to be checked

- len:

  Length of vector, default argument is 1

## Value

Boolean, whether logic was Boolean and of correct size

## Examples

``` r
bayesfam:::isLogic_len(c(TRUE, FALSE), 2) # should be TRUE
#> [1] TRUE
bayesfam:::isLogic_len(0, len = 1) # should be FALSE, 0 and 1 are numeric
#> [1] FALSE
```
