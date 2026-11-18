# Natural number vector check n \>= 0.

Natural number vector check n \>= 0.

## Usage

``` r
isNat_len(int, len = 1)
```

## Arguments

- int:

  Integer vector to be checked

- len:

  Length of vector, default argument is 1

## Value

Boolean, whether int was Integer \>= 0 and of correct size

## Examples

``` r
bayesfam:::isNat_len(c(1, 2), 2) # should be TRUE
#> [1] TRUE
bayesfam:::isNat_len(-1) # should be FALSE
#> [1] FALSE
```
